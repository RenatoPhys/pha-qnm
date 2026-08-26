#!/usr/bin/env python3
"""Production finite-momentum physics scans for the Route-A paper."""

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
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

import charged_hydrodynamics as hydro
import coupled_qnm as coupled
import run_numerics as base
from pha_qnm.thermodynamics import PHAModel
from run_posterior_uq import _charge_at_mu
from validate_coupled_qnms import _complex_root, ROBUST_SUBDOMAINS


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"
CHECKPOINT = RESULTS / "finite_k_checkpoint.json"
NOMINAL = {
    "intervals": 24, "r_uv": 5.0, "r0": 1.0e-6,
    "subdomains": ROBUST_SUBDOMAINS,
}
MAP_CUSP = {
    "phi0": 3.528280966207555,
    "charge_fraction": 0.34637356467181973,
    "T_c_MeV": 103.75528111978664,
    "mu_B_c_MeV": 593.3226560539789,
}

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def load_checkpoint() -> dict[str, dict[str, Any]]:
    if not CHECKPOINT.exists():
        return {}
    return json.loads(CHECKPOINT.read_text(encoding="utf-8"))


def save_checkpoint(cache: dict[str, dict[str, Any]]) -> None:
    CHECKPOINT.write_text(json.dumps(cache, indent=2, sort_keys=True) + "\n",
                          encoding="utf-8")


def root_row(cache: dict[str, dict[str, Any]], key: str,
             bg: dict[str, object], sector: str, qhat: float,
             seed: complex, extra: dict[str, Any],
             imaginary_axis: bool = False) -> dict[str, Any]:
    if key in cache:
        return cache[key]
    evaluator = lambda omega: coupled.collocation_source_matrix(
        bg, sector, qhat, omega, **NOMINAL
    )
    if imaginary_axis:
        omega, result, optimizer = imaginary_root(evaluator, seed)
    else:
        omega, result, optimizer = _complex_root(evaluator, seed)
    diagnostics = result["diagnostic_relative_maximum"]
    row = {
        **extra,
        "sector": sector,
        "qhat": qhat,
        "omega_numeric_real": omega.real,
        "omega_numeric_imag": omega.imag,
        "omega_hat_real": 2.0 * omega.real,
        "omega_hat_imag": 2.0 * omega.imag,
        "source_determinant_abs": abs(result["normalized_determinant"]),
        "source_singular_value": result["source_singular_value"],
        "horizon_gauge_defect": result["horizon_gauge_defect"],
        "primitive_constraint_maximum": max(diagnostics.values()),
        "linear_residual": result["linear_residual"],
        "optimizer_success": optimizer["optimizer_success"],
        "accepted_by_residual": optimizer["accepted_by_residual"],
    }
    cache[key] = row
    save_checkpoint(cache)
    return row


def imaginary_root(evaluator, seed: complex) -> tuple[complex, dict[str, object], dict[str, object]]:
    """Find a symmetry-protected purely imaginary pole without a rank-1 2D solve."""
    if seed.imag == 0.0:
        raise ValueError("an imaginary-axis root requires a nonzero imaginary seed")
    sign = math.copysign(1.0, seed.imag)
    scale = abs(seed.imag)
    evaluations: dict[float, dict[str, object]] = {}

    def result_at(value: float) -> dict[str, object]:
        key = float(value)
        if key not in evaluations:
            evaluations[key] = evaluator(1.0j * key)
        return evaluations[key]

    def scan(ordinates: np.ndarray):
        determinants = np.asarray([
            result_at(value)["normalized_determinant"] for value in ordinates
        ])
        component = (np.real if np.ptp(determinants.real) >= np.ptp(determinants.imag)
                     else np.imag)
        values = component(determinants)
        brackets = [
            (ordinates[index], ordinates[index + 1])
            for index in range(len(ordinates) - 1)
            if values[index] * values[index + 1] <= 0.0
        ]
        return component, brackets

    # Usually the new root remains on the seed's side of the axis.
    ordinates = sign * scale * np.geomspace(0.2, 5.0, 23)
    component, brackets = scan(ordinates)
    if not brackets:
        # At a spinodal-band edge the continuously tracked pole changes sign.
        # Only then expand to a two-sided scan.
        magnitudes = scale * np.geomspace(0.03, 30.0, 31)
        ordinates = np.sort(np.concatenate((-magnitudes, magnitudes)))
        component, brackets = scan(ordinates)
        if not brackets:
            raise RuntimeError(f"could not bracket the imaginary pole near {seed}")
    brackets.sort(key=lambda pair: abs(abs(0.5 * sum(pair)) - abs(seed.imag)))
    left, right = sorted(brackets[0])
    zero = brentq(
        lambda value: float(component(result_at(value)["normalized_determinant"])),
        left, right, xtol=2.0e-13, rtol=2.0e-12,
    )
    # The two determinant components acquire slightly different roundoff zeros
    # close to a critical pole.  Refine the bracket by minimizing the smallest
    # source singular value, which is phase-independent and is the actual QNM
    # boundary condition.
    refinement = minimize_scalar(
        lambda value: math.log(max(
            float(result_at(value)["source_singular_value"]), 1.0e-300
        )),
        bounds=(left, right), method="bounded",
        options={"xatol": 2.0e-13, "maxiter": 80},
    )
    candidates = (zero, float(refinement.x))
    ordinate = min(candidates, key=lambda value: result_at(value)["source_singular_value"])
    omega = 1.0j * ordinate
    result = result_at(ordinate)
    determinant = abs(result["normalized_determinant"])
    if determinant >= 2.0e-8:
        raise RuntimeError(
            f"imaginary root at {omega} has |det S|={determinant:.3e}"
        )
    return omega, result, {
        "optimizer_success": True,
        "optimizer_message": "bracketed imaginary-axis zero",
        "function_evaluations": len(evaluations),
        "accepted_by_residual": True,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def stable_scans(cache: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    backgrounds = (
        ("neutral", 1.0, 0.0),
        ("charged", 2.0, 0.1869733791),
    )
    q_values = (0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12)
    rows = []
    summaries = {}
    for label, phi0, charge in backgrounds:
        bg = base.integrate_background(phi0, charge)
        state = hydro.charged_hydro_state(phi0, charge, backend="reference_tight")
        c_s = math.sqrt(hydro.ideal_sound_speed_squared(state))
        D_eta = hydro.shear_diffusion_hat(state)
        D_B = hydro.baryon_diffusion_hat(state)
        hydro_info = {
            "T_MeV": state.T_MeV, "mu_MeV": state.mu_MeV,
            "c_s": c_s, "D_eta_hat": D_eta, "D_B_hat": D_B,
            "zeta_over_eta": state.zeta_over_eta,
            "sigma_Q_over_T": state.sigma_Q_over_T,
        }
        for qhat in q_values:
            hydro_modes = hydro.hydrodynamic_modes(state, qhat)
            sound_hydro = hydro_modes[np.argmax(hydro_modes.real)]
            diffusion_hydro = hydro_modes[np.argmin(np.abs(hydro_modes.real))]
            mode_specs = (
                ("shear", "helicity1", -0.5j * D_eta * qhat**2),
                ("sound", "helicity0", sound_hydro / 2.0),
                ("diffusion", "helicity0", diffusion_hydro / 2.0),
            )
            for mode, sector, seed in mode_specs:
                row = root_row(
                    cache, f"stable:{label}:{mode}:{qhat:.6f}",
                    bg, sector, qhat, complex(seed),
                    {"trajectory": label, "mode": mode, "phi0": phi0,
                     "charge_fraction": charge, **hydro_info},
                )
                rows.append(row)
        subset = lambda mode: [row for row in rows
                               if row["trajectory"] == label and row["mode"] == mode]
        shear = subset("shear"); sound = subset("sound"); diffusion = subset("diffusion")
        q2 = np.square([row["qhat"] for row in shear])
        shear_coeff = np.polyfit(q2, [-row["omega_hat_imag"] for row in shear], 2)
        diffusion_coeff = np.polyfit(
            q2, [-row["omega_hat_imag"] for row in diffusion], 2
        )
        sound_real = np.polyfit(
            q2, [row["omega_hat_real"] / row["qhat"] for row in sound], 1
        )
        sound_attenuation = np.polyfit(
            q2, [-row["omega_hat_imag"] / row["qhat"]**2 for row in sound], 1
        )
        summaries[label] = {
            **hydro_info,
            "qnm_D_eta_hat": float(shear_coeff[-2]),
            "qnm_D_B_hat": float(diffusion_coeff[-2]),
            "qnm_c_s": float(sound_real[-1]),
            "qnm_sound_attenuation_hat": float(sound_attenuation[-1]),
            "maximum_source_determinant": max(
                row["source_determinant_abs"] for row in rows
                if row["trajectory"] == label
            ),
            "maximum_primitive_constraint": max(
                row["primitive_constraint_maximum"] for row in rows
                if row["trajectory"] == label
            ),
        }
    return rows, summaries


def map_critical_scaling() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    model = PHAModel()
    phi_c = MAP_CUSP["phi0"]
    charge_c = MAP_CUSP["charge_fraction"]
    T_c, mu_c = MAP_CUSP["T_c_MeV"], MAP_CUSP["mu_B_c_MeV"]
    offsets = (0.50, 0.35, 0.25, 0.18, 0.12, 0.08, 0.05, 0.03,
               0.02, 0.012, 0.008, 0.005)
    rows = []
    preferred = charge_c
    for offset in offsets:
        phi0 = phi_c - offset
        charge, equilibrium = _charge_at_mu(
            phi0, mu_c, model, hydro.TIGHT_BACKGROUND_NUMERICS,
            preferred_charge=preferred,
        )
        state = hydro.charged_hydro_state(
            phi0, charge, phi_step=1.0e-3, charge_step=1.0e-4,
            backend="reference_tight", model=model,
        )
        rows.append({
            "phi_offset": offset, "phi0": phi0,
            "charge_fraction": charge,
            "T_MeV": equilibrium["T_MeV"], "mu_MeV": equilibrium["mu_MeV"],
            "reduced_distance": math.hypot(
                equilibrium["T_MeV"] / T_c - 1.0,
                equilibrium["mu_MeV"] / mu_c - 1.0,
            ),
            "D_B_hat": hydro.baryon_diffusion_hat(state),
            "chi_B_over_T2": state.chi_B_over_T2,
            "thermo_condition": state.thermo_jacobian_condition,
            "maxwell_error": state.maxwell_relation_relative_error,
        })
        preferred = charge
    window_rows = []
    for count in range(4, 9):
        fit = rows[-count:]
        log_r = np.log([row["reduced_distance"] for row in fit])
        log_D = np.log([row["D_B_hat"] for row in fit])
        log_chi = np.log([row["chi_B_over_T2"] for row in fit])
        slope_D_r = np.polyfit(log_r, log_D, 1)[0]
        slope_chi_r = np.polyfit(log_r, log_chi, 1)[0]
        slope_D_chi = np.polyfit(log_chi, log_D, 1)[0]
        window_rows.append({
            "point_count": count,
            "D_vs_distance_exponent": float(slope_D_r),
            "chi_vs_distance_exponent": float(slope_chi_r),
            "D_vs_chi_exponent": float(slope_D_chi),
            "z_eta0": float(2.0 - 2.0 * slope_D_chi),
        })
    return rows, {
        "fit_windows": window_rows,
        "z_eta0_midpoint": float(np.median([row["z_eta0"] for row in window_rows])),
        "z_eta0_half_range": float(
            0.5 * np.ptp([row["z_eta0"] for row in window_rows])
        ),
    }


def cep_q4_scan(cache: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    bg = base.integrate_background(MAP_CUSP["phi0"], MAP_CUSP["charge_fraction"])
    q_values = (0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.15)
    rows = []
    for qhat in q_values:
        seed = -0.12j * qhat**4
        rows.append(root_row(
            cache, f"cepv2:q4:{qhat:.6f}", bg, "helicity0", qhat, seed,
            {"trajectory": "MAP_CEP", "mode": "diffusion",
             "phi0": MAP_CUSP["phi0"],
             "charge_fraction": MAP_CUSP["charge_fraction"]},
            imaginary_axis=True,
        ))
    fit_windows = []
    for count in range(4, len(rows) + 1):
        subset = rows[:count]
        exponent, intercept = np.polyfit(
            np.log([row["qhat"] for row in subset]),
            np.log([-row["omega_hat_imag"] for row in subset]), 1,
        )
        fit_windows.append({"point_count": count, "z": float(exponent),
                            "log_amplitude": float(intercept)})
    local = []
    for left, right in zip(rows[:-1], rows[1:]):
        local.append({
            "qhat_geometric_mean": math.sqrt(left["qhat"] * right["qhat"]),
            "z_local": float(math.log(
                (-right["omega_hat_imag"]) / (-left["omega_hat_imag"])
            ) / math.log(right["qhat"] / left["qhat"])),
        })
    return rows, {
        "fit_windows": fit_windows,
        "local_exponents": local,
        "z_small_q": fit_windows[0]["z"],
        "maximum_source_determinant": max(row["source_determinant_abs"] for row in rows),
        "maximum_primitive_constraint": max(
            row["primitive_constraint_maximum"] for row in rows
        ),
    }


def _spinodal_horizon_points() -> list[dict[str, float]]:
    # Fold positions on the independently continued mu_B=650 MeV MAP line.
    hot_fold = 3.1813818812558528
    cold_fold = 3.8839233703099354
    model = PHAModel()
    points = []
    preferred = 0.43
    for fraction in (0.10, 0.50, 0.90):
        phi0 = hot_fold + fraction * (cold_fold - hot_fold)
        charge, state = _charge_at_mu(
            phi0, 650.0, model, hydro.TIGHT_BACKGROUND_NUMERICS,
            preferred_charge=preferred,
        )
        points.append({
            "spinodal_fraction": fraction, "phi0": phi0,
            "charge_fraction": charge, "T_MeV": state["T_MeV"],
            "mu_MeV": state["mu_MeV"],
        })
        preferred = charge
    return points


def spinodal_scans(cache: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    q_values = (0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.16, 0.20,
                0.24, 0.26, 0.28, 0.30, 0.32, 0.36, 0.40, 0.42,
                0.45, 0.50, 0.60)
    rows = []
    summaries = {}
    for point in _spinodal_horizon_points():
        phi0, charge = point["phi0"], point["charge_fraction"]
        label = f"f{point['spinodal_fraction']:.2f}"
        bg = base.integrate_background(phi0, charge)
        state = hydro.charged_hydro_state(
            phi0, charge, phi_step=1.0e-3, charge_step=1.0e-4,
            backend="reference_tight",
        )
        D_B = hydro.baryon_diffusion_hat(state)
        previous = None
        branch = []
        for qhat in q_values:
            seed = previous if previous is not None else -0.5j * D_B * qhat**2
            row = root_row(
                cache, f"spinodalv2:{label}:{qhat:.6f}", bg, "helicity0", qhat,
                complex(seed),
                {**point, "trajectory": label, "mode": "diffusion",
                 "hydro_D_B_hat": D_B,
                 "hydro_omega_hat_imag": -D_B * qhat**2},
                imaginary_axis=True,
            )
            branch.append(row); rows.append(row)
            previous = complex(row["omega_numeric_real"], row["omega_numeric_imag"])
        q = np.asarray([row["qhat"] for row in branch])
        growth = np.asarray([row["omega_hat_imag"] for row in branch])
        positive = np.flatnonzero(growth > 0.0)
        if positive.size < 4 or positive[-1] == len(growth) - 1:
            raise RuntimeError(f"unstable band did not close for {label}")
        last = int(positive[-1])
        edge = brentq(CubicSpline(q, growth), q[last], q[last + 1])
        spline = CubicSpline(q[:last + 2], growth[:last + 2])
        optimum = minimize_scalar(lambda value: -float(spline(value)),
                                  bounds=(q[0], edge), method="bounded")
        q_star = float(optimum.x)
        gamma_star = float(spline(q_star))
        relative = np.abs(
            growth - np.asarray([row["hydro_omega_hat_imag"] for row in branch])
        ) / np.maximum(np.abs(growth), 1.0e-300)
        hydro_valid = np.flatnonzero((growth > 0.0) & (relative < 0.10))
        summaries[label] = {
            **point, "hydro_D_B_hat": D_B,
            "q_edge": float(edge), "q_star": q_star,
            "Gamma_star_hat": gamma_star,
            "hydrodynamic_10percent_qmax": (
                float(q[hydro_valid[-1]]) if hydro_valid.size else float("nan")
            ),
            "maximum_source_determinant": max(
                row["source_determinant_abs"] for row in branch
            ),
            "maximum_primitive_constraint": max(
                row["primitive_constraint_maximum"] for row in branch
            ),
        }
    return rows, summaries


def validation_crosschecks(cache: dict[str, dict[str, Any]],
                           spinodal_summary: dict[str, Any]) -> dict[str, Any]:
    checks = {}
    cases = [
        ("cep_q0p06", MAP_CUSP["phi0"], MAP_CUSP["charge_fraction"], 0.06,
         complex(cache["cepv2:q4:0.060000"]["omega_numeric_real"],
                 cache["cepv2:q4:0.060000"]["omega_numeric_imag"])),
    ]
    middle = spinodal_summary["f0.50"]
    q_mid = min((0.26, 0.28, 0.30), key=lambda value: abs(value - middle["q_star"]))
    spin_row = cache[f"spinodalv2:f0.50:{q_mid:.6f}"]
    cases.append((
        "spinodal_fastest", middle["phi0"], middle["charge_fraction"], q_mid,
        complex(spin_row["omega_numeric_real"], spin_row["omega_numeric_imag"]),
    ))
    for name, phi0, charge, qhat, seed in cases:
        bg = base.integrate_background(phi0, charge)
        shoot_eval = lambda omega: coupled.shooting_source_matrix(
            bg, "helicity0", qhat, omega, r0=1.0e-6, r_uv=5.0,
            rtol=2.0e-9, atol=2.0e-11, max_step=0.025,
        )
        shot, shot_result, _ = _complex_root(shoot_eval, seed)
        fine_eval = lambda omega: coupled.collocation_source_matrix(
            bg, "helicity0", qhat, omega, intervals=28, r0=5.0e-7,
            r_uv=5.5, subdomains=ROBUST_SUBDOMAINS,
        )
        fine, fine_result, _ = _complex_root(fine_eval, seed)
        checks[name] = {
            "qhat": qhat,
            "nominal_omega_hat": [2.0 * seed.real, 2.0 * seed.imag],
            "shooting_omega_hat": [2.0 * shot.real, 2.0 * shot.imag],
            "fine_omega_hat": [2.0 * fine.real, 2.0 * fine.imag],
            "shooting_distance": abs(2.0 * (shot - seed)),
            "fine_distance": abs(2.0 * (fine - seed)),
            "shooting_source_determinant": abs(shot_result["normalized_determinant"]),
            "fine_source_determinant": abs(fine_result["normalized_determinant"]),
            "fine_primitive_constraint": max(
                fine_result["diagnostic_relative_maximum"].values()
            ),
        }
    return checks


def make_figures(stable: list[dict[str, Any]], critical: list[dict[str, Any]],
                 cep: list[dict[str, Any]], spinodal: list[dict[str, Any]]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.2))
    colors = {"neutral": "#0072B2", "charged": "#D55E00"}
    for trajectory in colors:
        for mode, marker in (("shear", "o"), ("sound", "s"), ("diffusion", "^")):
            rows = [row for row in stable
                    if row["trajectory"] == trajectory and row["mode"] == mode]
            if mode == "sound":
                axes[0].plot([row["qhat"] for row in rows],
                             [row["omega_hat_real"] for row in rows], marker,
                             color=colors[trajectory], ms=3,
                             label=f"{trajectory} {mode}")
            axes[1].plot([row["qhat"] for row in rows],
                         [row["omega_hat_imag"] for row in rows], marker,
                         color=colors[trajectory], ms=3,
                         label=f"{trajectory} {mode}")
    axes[0].set(xlabel=r"$\hat q$", ylabel=r"$\mathrm{Re}\,\hat\omega$")
    axes[1].set(xlabel=r"$\hat q$", ylabel=r"$\mathrm{Im}\,\hat\omega$")
    axes[0].legend(frameon=False, fontsize=7)
    axes[2].loglog([row["reduced_distance"] for row in critical],
                   [row["D_B_hat"] for row in critical], "o-", ms=3,
                   label=r"$D_B$")
    axes[2].loglog([row["reduced_distance"] for row in critical],
                   np.reciprocal([row["chi_B_over_T2"] for row in critical]),
                   "s--", ms=3, label=r"$1/(\chi_B/T^2)$")
    axes[2].set(xlabel="reduced distance", ylabel="critical observable")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    for extension in ("pdf", "png"):
        fig.savefig(FIGURES / f"hydrodynamic_critical_dispersion.{extension}", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))
    axes[0].loglog([row["qhat"] for row in cep],
                   [-row["omega_hat_imag"] for row in cep], "o-", color="#0072B2")
    reference_q = np.asarray([cep[0]["qhat"], cep[-1]["qhat"]])
    amplitude = -cep[0]["omega_hat_imag"] / cep[0]["qhat"]**4
    axes[0].loglog(reference_q, amplitude * reference_q**4, "k--", lw=1,
                   label=r"$\hat q^4$")
    axes[0].set(xlabel=r"$\hat q$", ylabel=r"$-\mathrm{Im}\,\hat\omega_D$")
    axes[0].legend(frameon=False)
    palette = ("#0072B2", "#CC79A7", "#D55E00")
    for color, trajectory in zip(palette, ("f0.10", "f0.50", "f0.90")):
        rows = [row for row in spinodal if row["trajectory"] == trajectory]
        axes[1].plot([row["qhat"] for row in rows],
                     [row["omega_hat_imag"] for row in rows], "o-", ms=3,
                     color=color, label=trajectory.replace("f", "fraction "))
        axes[1].plot([row["qhat"] for row in rows],
                     [row["hydro_omega_hat_imag"] for row in rows], "--",
                     color=color, alpha=0.65)
    axes[1].axhline(0.0, color="0.25", lw=0.8)
    axes[1].set(xlabel=r"$\hat q$", ylabel=r"$\mathrm{Im}\,\hat\omega_D$")
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    for extension in ("pdf", "png"):
        fig.savefig(FIGURES / f"cep_spinodal_dynamics.{extension}", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("stable", "critical", "spinodal", "all"),
                        default="all")
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    if args.fresh and CHECKPOINT.exists():
        CHECKPOINT.unlink()
    cache = load_checkpoint()
    summary: dict[str, Any] = {}
    stable_rows = critical_rows = cep_rows = spinodal_rows = []
    if args.stage in ("stable", "all"):
        stable_rows, summary["stable"] = stable_scans(cache)
        write_csv(RESULTS / "finite_k_hydrodynamic_dispersion.csv", stable_rows)
    if args.stage in ("critical", "all"):
        critical_rows, summary["critical_scaling"] = map_critical_scaling()
        cep_rows, summary["cep_q4"] = cep_q4_scan(cache)
        write_csv(RESULTS / "cep_critical_scaling.csv", critical_rows)
        write_csv(RESULTS / "cep_q4_dispersion.csv", cep_rows)
    if args.stage in ("spinodal", "all"):
        spinodal_rows, summary["spinodal"] = spinodal_scans(cache)
        write_csv(RESULTS / "spinodal_dispersion.csv", spinodal_rows)
        summary["crosschecks"] = validation_crosschecks(cache, summary["spinodal"])
    (RESULTS / "finite_k_physics_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if args.stage == "all":
        make_figures(stable_rows, critical_rows, cep_rows, spinodal_rows)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
