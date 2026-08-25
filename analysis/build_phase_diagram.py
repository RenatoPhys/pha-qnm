#!/usr/bin/env python3
"""Build a branch-resolved MAP phase diagram from a horizon-data surface."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator
from scipy.optimize import brentq

from pha_qnm.thermodynamics import PHAModel, solve_background


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


@dataclass
class ConstantMuLine:
    mu: float
    phi: np.ndarray
    charge: np.ndarray
    temperature: np.ndarray
    entropy: np.ndarray
    pressure_relative: np.ndarray


def _surface_point(arguments: tuple[int, int, float, float]) -> tuple[int, int, dict | None]:
    i, j, phi0, charge = arguments
    try:
        return i, j, solve_background(phi0, charge, model=PHAModel())
    except (RuntimeError, ValueError, OverflowError):
        return i, j, None


def compute_surface(phi_values: np.ndarray, charge_values: np.ndarray,
                    workers: int = 1) -> dict[str, np.ndarray]:
    shape = (phi_values.size, charge_values.size)
    arrays = {name: np.full(shape, np.nan) for name in (
        "T", "mu", "s", "rho", "constraint", "gauss", "r_uv")}
    total = phi_values.size * charge_values.size
    completed = 0
    jobs = [(i, j, float(phi0), float(charge))
            for i, phi0 in enumerate(phi_values)
            for j, charge in enumerate(charge_values)]

    def accept(result: tuple[int, int, dict | None]) -> None:
        nonlocal completed
        i, j, state = result
        completed += 1
        if state is not None:
            arrays["T"][i, j] = state["T_MeV"]
            arrays["mu"][i, j] = state["mu_MeV"]
            arrays["s"][i, j] = state["s_MeV3"]
            arrays["rho"][i, j] = state["rho_MeV3"]
            arrays["constraint"][i, j] = state["constraint"]
            arrays["gauss"][i, j] = state["gauss_relative_drift"]
            arrays["r_uv"][i, j] = state["r_uv"]
        if completed % max(1, total // 20) == 0:
            print(f"surface {completed}/{total}", flush=True)

    if workers == 1:
        for job in jobs:
            accept(_surface_point(job))
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            for result in executor.map(_surface_point, jobs, chunksize=16):
                accept(result)
    return arrays


def surface_derivatives(phi: np.ndarray, charge: np.ndarray,
                        temperature: np.ndarray, chemical_potential: np.ndarray) -> np.ndarray:
    logT = np.log(temperature)
    logmu = np.log(np.where(chemical_potential > 0.0, chemical_potential, np.nan))
    T_phi, T_charge = np.gradient(logT, phi, charge, edge_order=2)
    mu_phi, mu_charge = np.gradient(logmu, phi, charge, edge_order=2)
    return T_phi * mu_charge - T_charge * mu_phi


def constant_mu_line(mu_target: float, phi: np.ndarray, charge: np.ndarray,
                     interpolators: dict[str, RegularGridInterpolator]) -> ConstantMuLine | None:
    q_values = []
    valid_phi = []
    for value in phi:
        points = np.column_stack((np.full(charge.size, value), charge))
        mu_row = interpolators["mu"](points)
        finite = np.isfinite(mu_row)
        if np.count_nonzero(finite) < 2:
            continue
        row_mu = mu_row[finite]
        row_q = charge[finite]
        order = np.argsort(row_mu)
        row_mu = row_mu[order]
        row_q = row_q[order]
        if not row_mu[0] <= mu_target <= row_mu[-1]:
            continue
        q_values.append(float(np.interp(mu_target, row_mu, row_q)))
        valid_phi.append(float(value))
    if len(valid_phi) < 12:
        return None
    valid_phi = np.asarray(valid_phi)
    q_values = np.asarray(q_values)
    points = np.column_stack((valid_phi, q_values))
    temperature = interpolators["T"](points)
    entropy = interpolators["s"](points)
    finite = np.isfinite(temperature) & np.isfinite(entropy)
    valid_phi, q_values = valid_phi[finite], q_values[finite]
    temperature, entropy = temperature[finite], entropy[finite]
    if valid_phi.size < 12:
        return None
    dense_phi = np.linspace(valid_phi[0], valid_phi[-1], max(1001, 4 * valid_phi.size + 1))
    q_values = CubicSpline(valid_phi, q_values)(dense_phi)
    temperature = CubicSpline(valid_phi, temperature)(dense_phi)
    entropy = CubicSpline(valid_phi, entropy)(dense_phi)
    valid_phi = dense_phi
    # Along dmu=0, the first law gives dp=s dT.  The integration constant is
    # irrelevant for pressure differences between branches at the same mu.
    pressure = np.concatenate(([0.0], np.cumsum(
        0.5 * (entropy[1:] + entropy[:-1]) * (temperature[1:] - temperature[:-1])
    )))
    return ConstantMuLine(mu_target, valid_phi, q_values, temperature, entropy, pressure)


def phase_points_on_constant_mu(line: ConstantMuLine) -> dict | None:
    phi = line.phi
    temperature = line.temperature
    derivative = np.gradient(temperature, phi, edge_order=2)
    sign_changes = np.flatnonzero(derivative[:-1] * derivative[1:] < 0.0)
    if sign_changes.size < 2:
        return None
    spinodal_phi = []
    for index in sign_changes[:2]:
        spline = CubicSpline(phi[max(0, index - 2):min(phi.size, index + 4)],
                             derivative[max(0, index - 2):min(phi.size, index + 4)])
        try:
            root = brentq(spline, phi[index], phi[index + 1])
        except ValueError:
            root = 0.5 * (phi[index] + phi[index + 1])
        spinodal_phi.append(root)
    spinodal_phi.sort()
    T_spline = CubicSpline(phi, temperature)
    q_spline = CubicSpline(phi, line.charge)
    p_spline = CubicSpline(phi, line.pressure_relative)
    T_spinodal = [float(T_spline(value)) for value in spinodal_phi]

    left_interval = (phi[0], spinodal_phi[0])
    right_interval = (spinodal_phi[1], phi[-1])
    left_end_temperatures = [float(T_spline(left_interval[0])), float(T_spline(left_interval[1]))]
    right_end_temperatures = [float(T_spline(right_interval[0])), float(T_spline(right_interval[1]))]
    lower_T = max(min(left_end_temperatures), min(right_end_temperatures))
    upper_T = min(max(left_end_temperatures), max(right_end_temperatures))
    if not lower_T < upper_T:
        return None

    def outer_root(target_T: float, interval: tuple[float, float]) -> float:
        return brentq(lambda value: float(T_spline(value) - target_T), *interval)

    def pressure_difference(target_T: float) -> float:
        left = outer_root(target_T, left_interval)
        right = outer_root(target_T, right_interval)
        return float(p_spline(right) - p_spline(left))

    sample_T = np.linspace(lower_T + 1.0e-7, upper_T - 1.0e-7, 80)
    try:
        differences = np.array([pressure_difference(value) for value in sample_T])
    except ValueError:
        return {
            "mu_MeV": line.mu,
            "spinodal_low_T_MeV": min(T_spinodal),
            "spinodal_high_T_MeV": max(T_spinodal),
            "coexistence": None,
        }
    brackets = np.flatnonzero(differences[:-1] * differences[1:] <= 0.0)
    if brackets.size == 0:
        return {
            "mu_MeV": line.mu,
            "spinodal_low_T_MeV": min(T_spinodal),
            "spinodal_high_T_MeV": max(T_spinodal),
            "coexistence": None,
        }
    bracket = brackets[0]
    coexistence_T = brentq(pressure_difference, sample_T[bracket], sample_T[bracket + 1])
    phi_left = outer_root(coexistence_T, left_interval)
    phi_right = outer_root(coexistence_T, right_interval)
    return {
        "mu_MeV": line.mu,
        "spinodal_low_T_MeV": min(T_spinodal),
        "spinodal_high_T_MeV": max(T_spinodal),
        "coexistence": {
            "T_MeV": float(coexistence_T),
            "phi0_hot": float(phi_left),
            "phi0_cold": float(phi_right),
            "charge_hot": float(q_spline(phi_left)),
            "charge_cold": float(q_spline(phi_right)),
            "latent_entropy_MeV3": float(CubicSpline(phi, line.entropy)(phi_left)
                                         - CubicSpline(phi, line.entropy)(phi_right)),
            "pressure_residual_MeV4": float(pressure_difference(coexistence_T)),
        },
    }


def save_surface(phi: np.ndarray, charge: np.ndarray, arrays: dict[str, np.ndarray],
                 determinant: np.ndarray) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(RESULTS / "phase_surface.npz", phi0=phi, charge_fraction=charge,
                        determinant=determinant, **arrays)
    fields = ["phi0", "charge_fraction", "T_MeV", "mu_MeV", "s_MeV3", "rho_MeV3",
              "det_log_map", "constraint", "gauss_relative_drift", "r_uv"]
    with (RESULTS / "phase_surface.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i, phi0 in enumerate(phi):
            for j, q in enumerate(charge):
                writer.writerow({
                    "phi0": phi0, "charge_fraction": q,
                    "T_MeV": arrays["T"][i, j], "mu_MeV": arrays["mu"][i, j],
                    "s_MeV3": arrays["s"][i, j], "rho_MeV3": arrays["rho"][i, j],
                    "det_log_map": determinant[i, j],
                    "constraint": arrays["constraint"][i, j],
                    "gauss_relative_drift": arrays["gauss"][i, j],
                    "r_uv": arrays["r_uv"][i, j],
                })


def plot_phase_diagram(phase_rows: list[dict]) -> None:
    usable = [row for row in phase_rows if row is not None]
    spinodal_mu = np.array([row["mu_MeV"] for row in usable])
    spinodal_low = np.array([row["spinodal_low_T_MeV"] for row in usable])
    spinodal_high = np.array([row["spinodal_high_T_MeV"] for row in usable])
    coexistence = [row for row in usable if row["coexistence"] is not None]
    fig, ax = plt.subplots(figsize=(5.2, 3.5))
    ax.plot(spinodal_mu, spinodal_low, color="#D55E00", ls="--", label="spinodals")
    ax.plot(spinodal_mu, spinodal_high, color="#D55E00", ls="--")
    if coexistence:
        ax.plot([row["mu_MeV"] for row in coexistence],
                [row["coexistence"]["T_MeV"] for row in coexistence],
                color="#0072B2", lw=2.0, label="coexistence")
    ax.scatter([593.323], [103.755], marker="*", s=90, color="#009E73",
               edgecolor="white", linewidth=0.6, zorder=5, label="independent cusp")
    ax.scatter([602.475], [103.898], marker="x", s=42, color="#6B7280",
               zorder=5, label="HDF5 metadata")
    ax.set(xlabel=r"$\mu_B$ [MeV]", ylabel=r"$T$ [MeV]")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / "phase_diagram.pdf")
    fig.savefig(FIGURES / "phase_diagram.png", dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-phi", type=int, default=121)
    parser.add_argument("--n-charge", type=int, default=101)
    parser.add_argument("--phi-min", type=float, default=2.6)
    parser.add_argument("--phi-max", type=float, default=5.2)
    parser.add_argument("--charge-min", type=float, default=0.20)
    parser.add_argument("--charge-max", type=float, default=0.62)
    parser.add_argument("--mu-min", type=float, default=595.0)
    parser.add_argument("--mu-max", type=float, default=850.0)
    parser.add_argument("--mu-count", type=int, default=52)
    parser.add_argument("--reuse-surface", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    phi = np.linspace(args.phi_min, args.phi_max, args.n_phi)
    charge = np.linspace(args.charge_min, args.charge_max, args.n_charge)
    surface_path = RESULTS / "phase_surface.npz"
    arrays = None
    if args.reuse_surface and surface_path.exists():
        stored = np.load(surface_path)
        if np.array_equal(stored["phi0"], phi) and np.array_equal(stored["charge_fraction"], charge):
            arrays = {name: stored[name] for name in ("T", "mu", "s", "rho", "constraint", "gauss", "r_uv")}
            print("reused stored phase surface", flush=True)
    if arrays is None:
        arrays = compute_surface(phi, charge, workers=args.workers)
    determinant = surface_derivatives(phi, charge, arrays["T"], arrays["mu"])
    save_surface(phi, charge, arrays, determinant)
    interpolators = {name: RegularGridInterpolator((phi, charge), values,
                                                    bounds_error=False, fill_value=np.nan)
                     for name, values in arrays.items() if name in {"T", "mu", "s", "rho"}}
    mu_values = np.linspace(args.mu_min, args.mu_max, args.mu_count)
    phase_rows = []
    for index, mu_target in enumerate(mu_values):
        line = constant_mu_line(float(mu_target), phi, charge, interpolators)
        phase_rows.append(phase_points_on_constant_mu(line) if line else None)
        print(f"constant-mu continuation {index + 1}/{mu_values.size}", flush=True)
    serializable = [row for row in phase_rows if row is not None]
    (RESULTS / "phase_lines.json").write_text(json.dumps(serializable, indent=2), encoding="utf-8")
    with (RESULTS / "phase_lines.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["mu_MeV", "spinodal_low_T_MeV", "spinodal_high_T_MeV",
                  "coexistence_T_MeV", "phi0_hot", "phi0_cold", "charge_hot",
                  "charge_cold", "latent_entropy_MeV3", "pressure_residual_MeV4"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in serializable:
            coexistence = row["coexistence"] or {}
            writer.writerow({
                "mu_MeV": row["mu_MeV"],
                "spinodal_low_T_MeV": row["spinodal_low_T_MeV"],
                "spinodal_high_T_MeV": row["spinodal_high_T_MeV"],
                "coexistence_T_MeV": coexistence.get("T_MeV", ""),
                "phi0_hot": coexistence.get("phi0_hot", ""),
                "phi0_cold": coexistence.get("phi0_cold", ""),
                "charge_hot": coexistence.get("charge_hot", ""),
                "charge_cold": coexistence.get("charge_cold", ""),
                "latent_entropy_MeV3": coexistence.get("latent_entropy_MeV3", ""),
                "pressure_residual_MeV4": coexistence.get("pressure_residual_MeV4", ""),
            })
    plot_phase_diagram(serializable)
    summary = {
        "surface_points": int(phi.size * charge.size),
        "surface_grid": {
            "n_phi": int(phi.size), "n_charge": int(charge.size),
            "phi_min": float(phi[0]), "phi_max": float(phi[-1]),
            "charge_min": float(charge[0]), "charge_max": float(charge[-1]),
            "mu_min_MeV": float(mu_values[0]), "mu_max_MeV": float(mu_values[-1]),
            "mu_count": int(mu_values.size),
        },
        "failed_surface_points": int(np.count_nonzero(~np.isfinite(arrays["T"]))),
        "max_abs_constraint": float(np.nanmax(np.abs(arrays["constraint"]))),
        "max_abs_gauss_drift": float(np.nanmax(np.abs(arrays["gauss"]))),
        "constant_mu_lines": len(serializable),
        "coexistence_points": sum(row["coexistence"] is not None for row in serializable),
    }
    (RESULTS / "phase_diagram_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
