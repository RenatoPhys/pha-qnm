#!/usr/bin/env python3
"""Grid convergence study for spinodal and coexistence lines."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from build_phase_diagram import constant_mu_line, phase_points_on_constant_mu


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"


def main() -> None:
    surface = np.load(RESULTS / "phase_surface.npz")
    mu_targets = np.arange(598.0, 850.1, 4.0)
    records = []
    by_stride = {}
    for stride in (4, 2, 1):
        phi = surface["phi0"][::stride]
        charge = surface["charge_fraction"][::stride]
        arrays = {name: surface[name][::stride, ::stride] for name in ("T", "mu", "s", "rho")}
        interpolators = {
            name: RegularGridInterpolator((phi, charge), values,
                                          bounds_error=False, fill_value=np.nan)
            for name, values in arrays.items()
        }
        by_stride[stride] = {}
        for mu_target in mu_targets:
            line = constant_mu_line(float(mu_target), phi, charge, interpolators)
            result = phase_points_on_constant_mu(line) if line else None
            if result is None:
                continue
            coexistence = result["coexistence"] or {}
            record = {
                "stride": stride,
                "n_phi": phi.size,
                "n_charge": charge.size,
                "mu_MeV": mu_target,
                "spinodal_low_T_MeV": result["spinodal_low_T_MeV"],
                "spinodal_high_T_MeV": result["spinodal_high_T_MeV"],
                "coexistence_T_MeV": coexistence.get("T_MeV", np.nan),
            }
            records.append(record)
            by_stride[stride][mu_target] = record

    reference = by_stride[1]
    comparisons = []
    for stride in (4, 2):
        for mu_target, record in by_stride[stride].items():
            if mu_target not in reference:
                continue
            ref = reference[mu_target]
            comparisons.append({
                "stride": stride,
                "mu_MeV": mu_target,
                "delta_spinodal_low_T_MeV": record["spinodal_low_T_MeV"] - ref["spinodal_low_T_MeV"],
                "delta_spinodal_high_T_MeV": record["spinodal_high_T_MeV"] - ref["spinodal_high_T_MeV"],
                "delta_coexistence_T_MeV": record["coexistence_T_MeV"] - ref["coexistence_T_MeV"],
            })

    with (RESULTS / "phase_line_grid_convergence.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    with (RESULTS / "phase_line_grid_differences.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]))
        writer.writeheader(); writer.writerows(comparisons)

    finite_coexistence = [abs(row["delta_coexistence_T_MeV"]) for row in comparisons
                          if np.isfinite(row["delta_coexistence_T_MeV"])]
    summary = {
        "reference_grid": {"n_phi": int(surface["phi0"].size),
                           "n_charge": int(surface["charge_fraction"].size)},
        "comparison_rows": len(comparisons),
        "max_abs_spinodal_low_delta_T_MeV": max(abs(row["delta_spinodal_low_T_MeV"])
                                                  for row in comparisons),
        "max_abs_spinodal_high_delta_T_MeV": max(abs(row["delta_spinodal_high_T_MeV"])
                                                   for row in comparisons),
        "max_abs_coexistence_delta_T_MeV": max(finite_coexistence),
    }
    (RESULTS / "phase_line_grid_convergence.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
