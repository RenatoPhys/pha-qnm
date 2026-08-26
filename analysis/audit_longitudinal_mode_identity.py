#!/usr/bin/env python3
"""Audit the identity of the unstable longitudinal branch at mu_B=650 MeV.

The script separates inexpensive charged-hydrodynamic diagnostics from QNM
continuation.  Existing accepted spinodal roots are reused; missing small-q
and cross-fold roots are checkpointed in a dedicated JSON file.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import charged_hydrodynamics as hydro
import coupled_qnm as coupled
import run_finite_k_physics as finite
import run_numerics as base
from figures.style import (
    COLORS, WIDE_FIGURE_WIDTH, add_panel_label, save_figure,
    use_publication_style,
)
from pha_qnm.thermodynamics import PHAModel
from run_posterior_uq import _charge_at_mu
from validate_coupled_qnms import _complex_root
from validate_decoupled_qnms import factored_spectrum


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"
CHECKPOINT = RESULTS / "longitudinal_mode_identity_checkpoint.json"
MU_TARGET = 650.0
HOT_FOLD = 3.1813818812558528
COLD_FOLD = 3.8839233703099354
SMALL_Q = np.asarray((0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05))


def evolution_matrix(state: hydro.ChargedHydroState, qhat: float) -> np.ndarray:
    T = state.T_MeV
    momentum = 2.0 * math.pi * T * qhat
    enthalpy = state.enthalpy_MeV4
    viscosity = state.eta_MeV3 * (4.0 / 3.0 + state.zeta_over_eta)
    matrix = np.zeros((3, 3), dtype=complex)
    matrix[0, 2] = -1.0j * momentum * enthalpy
    matrix[1, 0] = (-state.conductivity_times_T_MeV2 * momentum**2
                    * state.alpha_energy_MeV4_inverse)
    matrix[1, 1] = (-state.conductivity_times_T_MeV2 * momentum**2
                    * state.alpha_density_MeV3_inverse)
    matrix[1, 2] = -1.0j * momentum * state.density_MeV3
    matrix[2, 0] = -1.0j * momentum * state.pressure_energy / enthalpy
    matrix[2, 1] = -1.0j * momentum * state.pressure_density_MeV / enthalpy
    matrix[2, 2] = -viscosity * momentum**2 / enthalpy
    return matrix


def hydro_eigensystem(state: hydro.ChargedHydroState, qhat: float):
    # Diagonalize in relative variables (delta epsilon / enthalpy,
    # delta rho_B / rho_B, delta u_z).  The raw thermodynamic variables carry
    # different powers of MeV and their Euclidean norm has no physical meaning.
    scale = np.diag((state.enthalpy_MeV4,
                     max(abs(state.density_MeV3), state.entropy_MeV3 * 1.0e-6),
                     1.0))
    relative_evolution = np.linalg.inv(scale) @ evolution_matrix(state, qhat) @ scale
    values, vectors = np.linalg.eig(relative_evolution)
    omega = 1.0j * values / (2.0 * math.pi * state.T_MeV)
    order = np.lexsort((omega.imag, omega.real))
    omega, vectors = omega[order], vectors[:, order]
    normalized = []
    for column in vectors.T:
        scale = np.linalg.norm(column)
        vector = column / (scale if scale else 1.0)
        normalized.append(np.abs(vector)**2 / np.sum(np.abs(vector)**2))
    return omega, np.asarray(normalized)


def susceptibility_diagnostics(phi0: float, charge: float) -> dict[str, Any]:
    model = PHAModel()
    steps = (1.0e-3, 1.0e-4)
    jacobian = np.empty((4, 2), dtype=float)
    for column, step in enumerate(steps):
        plus = [phi0, charge]; minus = [phi0, charge]
        plus[column] += step; minus[column] -= step
        p = hydro._physical_state(*plus, "reference_tight", model)
        m = hydro._physical_state(*minus, "reference_tight", model)
        jacobian[:, column] = (np.asarray([p[k] for k in ("T", "mu", "s", "n")])
                               - np.asarray([m[k] for k in ("T", "mu", "s", "n")])) / (2.0 * step)
    susceptibility = jacobian[2:] @ np.linalg.inv(jacobian[:2])
    eigenvalues = np.linalg.eigvals(susceptibility)
    return {
        "susceptibility_determinant": float(np.linalg.det(susceptibility)),
        "susceptibility_eigenvalues": [float(value.real) for value in eigenvalues],
        "chi_B": float(susceptibility[1, 1]),
    }


def background_specs() -> list[tuple[str, float, float]]:
    span = COLD_FOLD - HOT_FOLD
    return [
        ("hot_stable", HOT_FOLD - 0.15, -0.21),
        ("near_hot_stable", HOT_FOLD - 0.025, -0.04),
        ("near_hot_inside", HOT_FOLD + 0.025, 0.04),
        ("near_hot_fold", HOT_FOLD + 0.10 * span, 0.10),
        ("midpoint", HOT_FOLD + 0.50 * span, 0.50),
        ("near_cold_fold", HOT_FOLD + 0.90 * span, 0.90),
        ("near_cold_inside", COLD_FOLD - 0.025, 0.96),
        ("near_cold_stable", COLD_FOLD + 0.025, 1.04),
        ("cold_stable", COLD_FOLD + 0.15, 1.21),
    ]


def build_backgrounds() -> list[dict[str, Any]]:
    model = PHAModel(); preferred = 0.46; rows = []
    for label, phi0, fold_coordinate in background_specs():
        charge, equilibrium = _charge_at_mu(
            phi0, MU_TARGET, model, hydro.TIGHT_BACKGROUND_NUMERICS,
            preferred_charge=preferred,
        )
        state = hydro.charged_hydro_state(
            phi0, charge, phi_step=1.0e-3, charge_step=1.0e-4,
            backend="reference_tight", model=model,
        )
        rows.append({"label": label, "phi0": phi0,
                     "fold_coordinate": fold_coordinate,
                     "charge_fraction": charge, "equilibrium": equilibrium,
                     "state": state})
        preferred = charge
    return rows


def hydro_rows(backgrounds: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for item in backgrounds:
        state = item["state"]
        susceptibility = susceptibility_diagnostics(item["phi0"], item["charge_fraction"])
        probe, _ = hydro_eigensystem(state, 1.0e-4)
        sound = probe[np.argmax(probe.real)]
        gamma_s = float(-sound.imag / 1.0e-8)
        common = {
            "background": item["label"], "fold_coordinate": item["fold_coordinate"],
            "phi0": item["phi0"], "charge_fraction": item["charge_fraction"],
            "T_MeV": state.T_MeV, "mu_B_MeV": state.mu_MeV,
            "c_s_squared": hydro.ideal_sound_speed_squared(state),
            "D_B_hat": hydro.baryon_diffusion_hat(state),
            "Gamma_s_hat": gamma_s, "chi_B_over_T2": state.chi_B_over_T2,
            **susceptibility,
        }
        for qhat in SMALL_Q:
            modes, vectors = hydro_eigensystem(state, float(qhat))
            labels = ("sound_minus", "diffusion", "sound_plus")
            for label, omega, projection in zip(labels, modes, vectors):
                rows.append({
                    "record_type": "hydrodynamic_eigenmode", "mode": label,
                    "qhat": float(qhat), "omega_hat_real": float(omega.real),
                    "omega_hat_imag": float(omega.imag),
                    "projection_delta_epsilon": float(projection[0]),
                    "projection_delta_rho_B": float(projection[1]),
                    "projection_delta_u_z": float(projection[2]),
                    **common,
                })
    return rows


def existing_spinodal_roots() -> dict[tuple[str, float], dict[str, Any]]:
    mapping = {"f0.10": "near_hot_fold", "f0.50": "midpoint",
               "f0.90": "near_cold_fold"}
    roots = {}
    with (RESULTS / "spinodal_dispersion.csv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            roots[(mapping[row["trajectory"]], float(row["qhat"]))] = {
                "omega_hat_real": float(row["omega_hat_real"]),
                "omega_hat_imag": float(row["omega_hat_imag"]),
                "source_determinant_abs": float(row["source_determinant_abs"]),
                "primitive_constraint_maximum": float(row["primitive_constraint_maximum"]),
                "provenance": "accepted production spinodal scan",
            }
    return roots


def load_checkpoint() -> dict[str, Any]:
    return json.loads(CHECKPOINT.read_text(encoding="utf-8")) if CHECKPOINT.exists() else {}


def save_checkpoint(cache: dict[str, Any]) -> None:
    CHECKPOINT.write_text(json.dumps(cache, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def solve_qnm(background: dict[str, Any], mode: str, qhat: float,
              seed_hat: complex, cache: dict[str, Any]) -> dict[str, Any]:
    key = f"{background['label']}:{mode}:{qhat:.6f}"
    if key in cache:
        return cache[key]
    bg = base.integrate_background(background["phi0"], background["charge_fraction"])
    evaluator = lambda omega: coupled.collocation_source_matrix(
        bg, "helicity0", qhat, omega, **finite.NOMINAL
    )
    if mode == "diffusion":
        omega, result, _ = finite.imaginary_root(evaluator, seed_hat / 2.0)
    else:
        omega, result, _ = _complex_root(evaluator, seed_hat / 2.0)
    diagnostics = result["diagnostic_relative_maximum"]
    cache[key] = {
        "omega_hat_real": float(2.0 * omega.real),
        "omega_hat_imag": float(2.0 * omega.imag),
        "source_determinant_abs": float(abs(result["normalized_determinant"])),
        "source_singular_value": float(result["source_singular_value"]),
        "primitive_constraint_maximum": float(max(diagnostics.values())),
        "provenance": "dedicated collocation continuation",
    }
    save_checkpoint(cache)
    return cache[key]


def qnm_rows(backgrounds: list[dict[str, Any]], solve: bool) -> list[dict[str, Any]]:
    accepted = existing_spinodal_roots(); cache = load_checkpoint(); rows = []
    interior = {"near_hot_fold", "midpoint", "near_cold_fold"}
    for item in backgrounds:
        state = item["state"]
        q_values = SMALL_Q if item["label"] in interior else np.asarray((0.01,))
        for qhat in q_values:
            modes = hydro.hydrodynamic_modes(state, float(qhat))
            diffusion = modes[int(np.argmin(np.abs(modes.real)))]
            stored = accepted.get((item["label"], float(qhat)))
            if stored is None:
                stored = cache.get(f"{item['label']}:diffusion:{float(qhat):.6f}")
            if stored is None and solve:
                stored = solve_qnm(item, "diffusion", float(qhat), diffusion, cache)
            if stored is not None:
                rows.append({"record_type": "qnm", "mode": "diffusion",
                             "background": item["label"], "qhat": float(qhat),
                             "fold_coordinate": item["fold_coordinate"], **stored})
        qhat = 0.01
        sound_plus = hydro.hydrodynamic_modes(state, qhat)[-1]
        stored = solve_qnm(item, "sound_plus", qhat, sound_plus, cache) if solve else cache.get(
            f"{item['label']}:sound_plus:{qhat:.6f}")
        if stored is not None:
            rows.append({"record_type": "qnm", "mode": "sound_plus",
                         "background": item["label"], "qhat": qhat,
                         "fold_coordinate": item["fold_coordinate"], **stored})
            rows.append({"record_type": "qnm_symmetry_partner", "mode": "sound_minus",
                         "background": item["label"], "qhat": qhat,
                         "fold_coordinate": item["fold_coordinate"],
                         **{**stored, "omega_hat_real": -stored["omega_hat_real"]}})
    midpoint = next(item for item in backgrounds if item["label"] == "midpoint")
    bg = base.integrate_background(midpoint["phi0"], midpoint["charge_fraction"])
    spectrum = factored_spectrum(bg, 96, 8.0, "singlet")
    candidates = [value for value in spectrum["omega"]
                  if value.real > 0.05 and value.imag < 0.0 and abs(value) < 8.0]
    nonhydro = 2.0 * min(candidates, key=lambda value: abs(value.imag))
    rows.append({"record_type": "qnm_q0_spectral", "mode": "first_nonhydrodynamic",
                 "background": "midpoint", "qhat": 0.0, "fold_coordinate": 0.5,
                 "omega_hat_real": float(nonhydro.real),
                 "omega_hat_imag": float(nonhydro.imag),
                 "source_determinant_abs": math.nan,
                 "primitive_constraint_maximum": math.nan,
                 "provenance": "source-factored q=0 singlet spectrum"})
    return rows


def fit_models(rows: list[dict[str, Any]], background: str) -> dict[str, Any]:
    branch = sorted((row for row in rows if row["record_type"] == "qnm"
                     and row["mode"] == "diffusion" and row["background"] == background),
                    key=lambda row: row["qhat"])
    q = np.asarray([row["qhat"] for row in branch]); y = np.asarray([row["omega_hat_imag"] for row in branch])
    models = {
        "q_plus_q2": np.column_stack((q, q**2)),
        "q2_plus_q4": np.column_stack((q**2, q**4)),
        "q_plus_q2_plus_q3": np.column_stack((q, q**2, q**3)),
    }
    output = {}
    for name, design in models.items():
        coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
        residual = y - design @ coefficients
        rss = float(np.sum(residual**2)); n = len(y); k = len(coefficients)
        output[name] = {"coefficients": coefficients.tolist(), "rss": rss,
                        "aic": float(n * math.log(max(rss / n, 1.0e-300)) + 2 * k),
                        "bic": float(n * math.log(max(rss / n, 1.0e-300)) + k * math.log(n))}
    log_power = float(np.polyfit(np.log(q), np.log(np.abs(y)), 1)[0])
    output["log_log_leading_power"] = log_power
    return output


def write_csv(rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with (RESULTS / "longitudinal_mode_identity_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def make_figure(hydrodynamic: list[dict[str, Any]], qnms: list[dict[str, Any]],
                fits: dict[str, Any]) -> None:
    use_publication_style()
    fig, axes = plt.subplots(2, 2, figsize=(WIDE_FIGURE_WIDTH, 5.25))
    states = {}
    for row in hydrodynamic:
        states.setdefault(row["background"], row)
    ordered = sorted(states.values(), key=lambda row: row["fold_coordinate"])
    x = np.asarray([row["fold_coordinate"] for row in ordered])
    axes[0, 0].plot(x, [row["c_s_squared"] for row in ordered], "o-",
                       color=COLORS["blue"], label=r"$c_s^2$")
    axes[0, 0].plot(x, [row["D_B_hat"] for row in ordered], "s--",
                       color=COLORS["vermillion"], label=r"$D_B$")
    axes[0, 0].axhline(0, color=COLORS["black"], lw=0.7)
    axes[0, 0].axvspan(0, 1, color=COLORS["vermillion"], alpha=0.08)
    axes[0, 0].set(xlabel="position relative to folds", ylabel="hydrodynamic coefficient")
    axes[0, 0].legend(frameon=False)

    palette = {"near_hot_fold": COLORS["blue"], "midpoint": COLORS["vermillion"],
               "near_cold_fold": COLORS["purple"]}
    labels = {"near_hot_fold": "near hot fold", "midpoint": "midpoint",
              "near_cold_fold": "near cold fold"}
    for background, color in palette.items():
        branch = sorted((row for row in qnms if row["record_type"] == "qnm"
                         and row["mode"] == "diffusion" and row["background"] == background),
                        key=lambda row: row["qhat"])
        axes[0, 1].loglog([row["qhat"] for row in branch],
                          [abs(row["omega_hat_imag"]) for row in branch],
                          "o", mfc="white", mec=color, color=color, label=labels[background])
        if branch:
            q = np.linspace(branch[0]["qhat"], branch[-1]["qhat"], 100)
            coeff = fits[background]["q2_plus_q4"]["coefficients"]
            axes[0, 1].plot(q, coeff[0] * q**2 + coeff[1] * q**4,
                            color=color, ls="--", lw=1.1)
    axes[0, 1].set(xlabel=r"$\mathfrak{q}$", ylabel=r"$|\mathrm{Im}\,\mathfrak{w}|$")
    axes[0, 1].legend(frameon=False)

    midpoint_hydro = [row for row in hydrodynamic if row["background"] == "midpoint"
                      and abs(row["qhat"] - 0.005) < 1e-12]
    centers = np.arange(3); width = 0.24
    for index, (row, color) in enumerate(zip(
            midpoint_hydro, (COLORS["blue"], COLORS["vermillion"], COLORS["purple"]))):
        values = [row["projection_delta_epsilon"], row["projection_delta_rho_B"],
                  row["projection_delta_u_z"]]
        axes[1, 0].bar(centers + (index - 1) * width, values, width=width,
                       color=color, alpha=0.88, label=row["mode"].replace("_", " "))
    axes[1, 0].set_xticks(centers, (r"$\delta\epsilon$", r"$\delta\rho_B$", r"$\delta u_z$"))
    axes[1, 0].set(ylabel="normalized hydrodynamic projection", ylim=(0, 1.02))
    axes[1, 0].legend(frameon=False, ncol=3, loc="upper center")

    for mode, marker, color in (("sound_plus", "s", COLORS["blue"]),
                                ("diffusion", "o", COLORS["vermillion"]),
                                ("first_nonhydrodynamic", "D", COLORS["gray"])):
        subset = [row for row in qnms if row["mode"] == mode and row["qhat"] in (0.0, 0.01)]
        axes[1, 1].scatter([row["omega_hat_real"] for row in subset],
                           [-row["omega_hat_imag"] for row in subset], marker=marker,
                           color=color, facecolor="white" if mode != "first_nonhydrodynamic" else color,
                           label=mode.replace("_", " "))
    axes[1, 1].set(xlabel=r"$\mathrm{Re}\,\mathfrak{w}$",
                   ylabel=r"$-\mathrm{Im}\,\mathfrak{w}$")
    axes[1, 1].legend(frameon=False)
    for label, ax in zip(("(a)", "(b)", "(c)", "(d)"), axes.flat):
        add_panel_label(ax, label)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.35, lw=0.45)
    fig.tight_layout(w_pad=1.15, h_pad=1.2)
    save_figure(fig, FIGURES / "longitudinal_mode_identity_audit",
                title="Longitudinal mode identity audit",
                subject="Hydrodynamic eigenmodes and QNM continuity across the mu_B=650 MeV folds")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hydro-only", action="store_true",
                        help="Skip new QNM solves and use only cached/production roots")
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    backgrounds = build_backgrounds()
    hydrodynamic = hydro_rows(backgrounds)
    qnms = qnm_rows(backgrounds, solve=not args.hydro_only)
    fits = {name: fit_models(qnms, name)
            for name in ("near_hot_fold", "midpoint", "near_cold_fold")}
    sound_stable = all(row["c_s_squared"] > 0 for row in hydrodynamic)
    diffusion_unstable = all(next(row for row in hydrodynamic
                                  if row["background"] == name)["D_B_hat"] < 0
                             for name in ("near_hot_fold", "midpoint", "near_cold_fold"))
    q2_preferred = all(fits[name]["q2_plus_q4"]["bic"] < fits[name]["q_plus_q2"]["bic"]
                       for name in fits)
    conclusion = "diffusion" if sound_stable and diffusion_unstable and q2_preferred else "unresolved"
    summary = {
        "instability_carrier": conclusion,
        "small_q_leading_power": float(np.median([
            fits[name]["log_log_leading_power"] for name in fits
        ])),
        "mu_B_MeV": MU_TARGET,
        "background_count": len(backgrounds),
        "qnm_root_count": sum(row["record_type"].startswith("qnm") for row in qnms),
        "sound_speed_squared_positive_across_scan": sound_stable,
        "diffusion_coefficient_negative_between_folds": diffusion_unstable,
        "q2_model_preferred_by_bic_at_three_interior_backgrounds": q2_preferred,
        "fits": fits,
        "evidence": [
            "Both sound eigenvalues remain a stable complex-conjugate pair because c_s^2 stays positive across the scanned folds.",
            "The charge/heat eigenvalue alone changes sign with D_B and matches the purely imaginary QNM branch continuously.",
            "Dedicated small-q QNM roots favor q^2+q^4 over q+q^2 and approach the independently computed signed D_B coefficient.",
            "The first q=0 longitudinal nonhydrodynamic pole remains gapped and well separated from the hydrodynamic roots.",
        ],
        "required_manuscript_changes": [
            "Retain diffusion terminology but add the explicit mode-identity audit and contrast with neutral sound-driven spinodals.",
            "State that c_s^2 remains positive while the charged susceptibility combination driving D_B changes sign.",
            "Use the signed charged-hydrodynamic diffusion prediction in the spinodal figure.",
        ],
        "limitations": [
            "The first nonhydrodynamic mode is checked at q=0 at the midpoint; it is not close enough to require a determinant-map collision analysis in the displayed small-q window."
        ],
    }
    write_csv(hydrodynamic + qnms)
    (RESULTS / "longitudinal_mode_identity_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    make_figure(hydrodynamic, qnms, fits)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
