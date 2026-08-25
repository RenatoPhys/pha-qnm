#!/usr/bin/env python3
"""Validated homogeneous QNM trajectories on physical thermodynamic paths."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

import run_numerics as base
from validate_decoupled_qnms import (factored_spectrum, nearest_mode,
                                     normalized_overlap, shooting_root)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
SECTORS = ("tensor", "vector", "singlet")
SEEDS = {
    "tensor": 0.9365852064372638 - 0.43344931655843155j,
    "vector": 0.8442710399991116 - 0.05158393566646418j,
    "singlet": 0.7462047610067616 - 0.12804252173623262j,
}


def fixed_ratio_background(phi0: float, target: float) -> tuple[dict[str, object], dict[str, float]]:
    cache: dict[float, tuple[dict[str, object], dict[str, float]]] = {}

    def evaluate(charge_fraction: float) -> float:
        key = round(float(charge_fraction), 13)
        if key not in cache:
            bg = base.integrate_background(phi0, float(charge_fraction))
            cache[key] = (bg, base.extract_uv(bg))
        return cache[key][1]["mu_over_T"] - target

    samples = np.linspace(0.0, 0.92, 47)
    lower = samples[0]
    f_lower = evaluate(lower)
    bracket = None
    for upper in samples[1:]:
        try:
            f_upper = evaluate(float(upper))
        except (RuntimeError, ValueError, OverflowError):
            continue
        if f_lower * f_upper <= 0.0:
            bracket = (float(lower), float(upper))
            break
        lower, f_lower = float(upper), f_upper
    if bracket is None:
        raise RuntimeError(f"No mu/T={target} bracket at phi_H={phi0}")
    fraction = float(brentq(evaluate, *bracket, xtol=2.0e-8, rtol=2.0e-10))
    bg = base.integrate_background(phi0, fraction)
    return bg, base.extract_uv(bg)


def trajectory_backgrounds(kind: str) -> list[tuple[dict[str, object], dict[str, float]]]:
    # Start at phi_H=3.5 for continuation, then move independently toward
    # high and low temperatures.  Output ordering is restored below.
    downward = np.arange(3.5, 0.74, -0.25)
    upward = np.arange(3.75, 5.01, 0.25)
    phis = np.concatenate((downward, upward))
    values = []
    for index, phi0 in enumerate(phis):
        if kind == "mu0":
            bg = base.integrate_background(float(phi0), 0.0)
            uv = base.extract_uv(bg)
        elif kind == "mu_over_T_2":
            bg, uv = fixed_ratio_background(float(phi0), 2.0)
        else:
            raise ValueError(kind)
        print(f"BACKGROUND {kind} {index + 1}/{len(phis)} phi_H={phi0:.2f} "
              f"T={uv['T_MeV']:.3f} mu/T={uv['mu_over_T']:.6f}", flush=True)
        values.append((bg, uv))
    return values


def main() -> None:
    critical = json.loads((RESULTS / "reference_cusp_uv1e5_h10.json").read_text(encoding="utf-8"))
    tc = float(critical["T_c_MeV"])
    muc = float(critical["mu_B_c_MeV"])
    rows: list[dict[str, object]] = []
    summary: dict[str, object] = {"trajectories": {}}
    for kind in ("mu0", "mu_over_T_2"):
        backgrounds = trajectory_backgrounds(kind)
        trajectory_rows: list[dict[str, object]] = []
        for sector in SECTORS:
            previous_root = SEEDS[sector]
            previous_spectrum = None
            previous_index = None
            anchor_root = None
            anchor_spectrum = None
            anchor_index = None
            for path_index, (bg, uv) in enumerate(backgrounds):
                if path_index == 12:
                    previous_root = anchor_root
                    previous_spectrum = anchor_spectrum
                    previous_index = anchor_index
                print(f"QNM {kind} {sector} {path_index + 1}/{len(backgrounds)}", flush=True)
                root = shooting_root(bg, sector, previous_root,
                                     acceptance_reference=previous_root,
                                     acceptance_radius=0.45, search_radius=0.45)
                if not root["success"]:
                    raise RuntimeError(
                        f"Shooting continuation failed: {kind} {sector} phi_H={bg['phi0']} "
                        f"root={root['omega']} residual={root['source_residual']}")
                spectrum = factored_spectrum(bg, 112, 8.0, sector)
                spectral_index = nearest_mode(spectrum, root["omega"])
                overlap = math.nan
                if previous_spectrum is not None and previous_index is not None:
                    overlap = normalized_overlap(
                        previous_spectrum["nodes"], previous_spectrum["right"][:, previous_index],
                        spectrum["nodes"], spectrum["right"][:, spectral_index])
                rho = float(uv["rho_over_T3"])
                row = {
                    "trajectory": kind, "path_index": path_index,
                    "parent_index": (0 if path_index == 12 else
                                     path_index - 1 if path_index else -1),
                    "sector": sector, "mode": 0,
                    "phi_H": bg["phi0"], "charge_fraction": bg["charge_fraction"],
                    "Phi1": bg["Phi1"], "T_MeV": uv["T_MeV"],
                    "mu_B_MeV": uv["mu_MeV"], "mu_B_over_T": uv["mu_over_T"],
                    "s_over_T3": uv["s_over_T3"], "n_B_over_T3": rho,
                    "s_over_n_B": uv["s_over_T3"] / rho if abs(rho) > 1.0e-14 else math.nan,
                    "branch_id": "stable_crossover", "locally_stable": True,
                    "inside_spinodal": False,
                    "reduced_distance_to_CEP": math.hypot(
                        (uv["T_MeV"] - tc) / tc, (uv["mu_MeV"] - muc) / muc),
                    "Re_omega_hat": 2.0 * root["omega"].real,
                    "Im_omega_hat": 2.0 * root["omega"].imag,
                    "shooting_source_residual": root["source_residual"],
                    "spectral_distance_hat": 2.0 * abs(
                        spectrum["omega"][spectral_index] - root["omega"]),
                    "spectral_overlap_parent": overlap,
                    "spectral_pencil_residual": spectrum["residual"][spectral_index],
                }
                trajectory_rows.append(row)
                previous_root = root["omega"]
                previous_spectrum = spectrum
                previous_index = spectral_index
                if path_index == 0:
                    anchor_root = previous_root
                    anchor_spectrum = previous_spectrum
                    anchor_index = previous_index
        rows.extend(trajectory_rows)
        summary["trajectories"][kind] = {
            "background_points": len(backgrounds),
            "qnm_points": len(trajectory_rows),
            "maximum_source_residual": max(row["shooting_source_residual"]
                                           for row in trajectory_rows),
            "maximum_spectral_distance_hat": max(row["spectral_distance_hat"]
                                                  for row in trajectory_rows),
            "minimum_parent_overlap": min(row["spectral_overlap_parent"]
                                          for row in trajectory_rows
                                          if math.isfinite(row["spectral_overlap_parent"])),
        }
    rows.sort(key=lambda row: (row["trajectory"], row["sector"], row["T_MeV"]))
    with (RESULTS / "homogeneous_physical_trajectories.csv").open(
            "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (RESULTS / "homogeneous_physical_trajectories.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
