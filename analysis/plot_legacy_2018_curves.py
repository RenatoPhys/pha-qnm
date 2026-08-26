#!/usr/bin/env python3
"""Independent curve reproduction using the 2018 EMD parametrization."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import charged_hydrodynamics as hydro
import compare_2018_emd as legacy
import run_numerics as base
from figures.style import (COLORS, WIDE_FIGURE_WIDTH, add_panel_label,
                           save_figure, use_publication_style)
from pha_qnm.thermodynamics import PHAModel
from run_posterior_uq import _charge_at_mu
from validate_decoupled_qnms import factored_spectrum


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"
CSV_PATH = RESULTS / "legacy_2018_qnm_curves.csv"
SUMMARY_PATH = RESULTS / "legacy_2018_qnm_curves_summary.json"
MU_VALUES = (0.0, 300.0, 500.0, 600.0, 700.0)
SECTORS = ("tensor", "vector", "singlet")


def select_mode(spectrum: dict, previous: complex | None) -> tuple[complex, int]:
    candidates = [(index, value) for index, value in enumerate(spectrum["omega"])
                  if value.real > 0.05 and value.imag < 0 and abs(value) < 8]
    if not candidates:
        raise RuntimeError("no physical candidate in factored spectrum")
    if previous is None:
        index, value = min(candidates, key=lambda pair: abs(pair[1].imag))
    else:
        index, value = min(candidates, key=lambda pair: abs(pair[1] - previous))
    return value, index


def compute_rows() -> tuple[list[dict], int]:
    legacy._activate(legacy.OLD_PARAMETERS)
    model = PHAModel(legacy.OLD_PARAMETERS)
    rows = []; rejected = 0
    phi_values = np.linspace(0.65, 4.10, 13)
    for mu_target in MU_VALUES:
        preferred = min(0.44, 0.00055 * mu_target)
        states = []
        for phi0 in phi_values:
            try:
                if mu_target == 0:
                    charge = 0.0
                    background = base.integrate_background(float(phi0), charge)
                    uv = base.extract_uv(background)
                else:
                    charge, reference = _charge_at_mu(
                        float(phi0), mu_target, model, hydro.TIGHT_BACKGROUND_NUMERICS,
                        preferred_charge=preferred,
                    )
                    background = base.integrate_background(float(phi0), charge)
                    uv = base.extract_uv(background)
                if not 75.0 <= float(uv["T_MeV"]) <= 450.0:
                    continue
                states.append((float(phi0), float(charge), background, uv))
                preferred = float(charge)
            except (RuntimeError, ValueError, OverflowError):
                rejected += 1
        if len(states) < 5:
            raise RuntimeError(f"insufficient 2018-model backgrounds at mu_B={mu_target}")
        for sector in SECTORS:
            previous = None
            for sequence, (phi0, charge, background, uv) in enumerate(states):
                spectrum = factored_spectrum(background, 112, 8.0, sector)
                omega, index = select_mode(spectrum, previous)
                overlap = math.nan
                if previous is not None:
                    overlap = 1.0 / (1.0 + abs(omega - previous))
                rows.append({
                    "mu_B_target_MeV": mu_target,
                    "sequence": sequence,
                    "sector": sector,
                    "phi0": phi0,
                    "charge_fraction": charge,
                    "T_MeV": float(uv["T_MeV"]),
                    "mu_B_MeV": float(uv["mu_MeV"]),
                    "Re_omega_hat": float(2.0 * omega.real),
                    "Im_omega_hat": float(2.0 * omega.imag),
                    "pencil_residual": float(spectrum["residual"][index]),
                    "tracking_continuity_proxy": float(overlap),
                    "calculation": "source-factored collocation point",
                    "line_role": "piecewise-linear guide to the eye",
                })
                previous = omega
    return rows, rejected


def write_products(rows: list[dict], rejected: int) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    benchmark = json.loads((RESULTS / "rougemont2018_benchmark.json").read_text(encoding="utf-8"))
    summary = {
        "description": "Independent numerical reproduction using the 2018 EMD parametrization",
        "chemical_potentials_MeV": list(MU_VALUES),
        "sectors": list(SECTORS),
        "background_count": len({(row["mu_B_target_MeV"], row["sequence"]) for row in rows}),
        "curve_point_count": len(rows),
        "rejected_background_attempts": rejected,
        "maximum_curve_pencil_residual": max(row["pencil_residual"] for row in rows),
        "benchmark_maximum_shooting_source_residual": benchmark["maximum_shooting_source_residual"],
        "benchmark_maximum_spectral_distance_hat": benchmark["maximum_spectral_distance_hat"],
        "nine_landmark_values_reproduced": benchmark["all_published_rounded_times_reproduced"],
        "interpolation": "piecewise linear, used only to guide the eye",
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_rows() -> list[dict]:
    with CSV_PATH.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key in ("mu_B_target_MeV", "sequence", "T_MeV", "Re_omega_hat", "Im_omega_hat"):
            row[key] = float(row[key])
    return rows


def plot(rows: list[dict]) -> None:
    use_publication_style()
    colors = {0.0: COLORS["black"], 300.0: COLORS["blue"],
              500.0: COLORS["green"], 600.0: COLORS["orange"],
              700.0: COLORS["vermillion"]}
    labels = {"tensor": "quintuplet", "vector": "triplet", "singlet": "singlet"}
    fig, axes = plt.subplots(3, 2, figsize=(WIDE_FIGURE_WIDTH, 6.75), sharex="col")
    for row_index, sector in enumerate(SECTORS):
        for mu in MU_VALUES:
            group = sorted((row for row in rows if row["sector"] == sector
                            and row["mu_B_target_MeV"] == mu), key=lambda row: row["sequence"])
            T = [row["T_MeV"] for row in group]
            real = [row["Re_omega_hat"] for row in group]
            damping = [-row["Im_omega_hat"] for row in group]
            axes[row_index, 0].plot(T, real, color=colors[mu], lw=0.9, alpha=0.8)
            axes[row_index, 0].plot(T, real, ls="none", marker="o", ms=2.6,
                                    mfc="white", mec=colors[mu])
            axes[row_index, 1].plot(T, damping, color=colors[mu], lw=0.9, alpha=0.8)
            axes[row_index, 1].plot(T, damping, ls="none", marker="o", ms=2.6,
                                    mfc="white", mec=colors[mu])
        axes[row_index, 0].set_ylabel(rf"{labels[sector]}  $\mathrm{{Re}}\,\mathfrak{{w}}$")
        axes[row_index, 1].set_ylabel(rf"{labels[sector]}  $-\mathrm{{Im}}\,\mathfrak{{w}}$")
    axes[-1, 0].set_xlabel(r"$T$ [MeV]")
    axes[-1, 1].set_xlabel(r"$T$ [MeV]")
    handles = [plt.Line2D([], [], color=colors[mu], marker="o", mfc="white",
                          label=rf"$\mu_B={mu:.0f}$ MeV") for mu in MU_VALUES]
    axes[0, 0].legend(handles=handles, frameon=False, ncol=2, loc="best")
    for tag, ax in zip(("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), axes.flat):
        add_panel_label(ax, tag)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.3, lw=0.4)
    fig.tight_layout(w_pad=1.0, h_pad=0.85)
    save_figure(fig, FIGURES / "legacy_2018_qnm_curves",
                title="Independent reproduction of 2018 EMD QNM curves",
                subject="Quintuplet, triplet, and singlet trajectories at five chemical potentials")
    plt.close(fig)


def main() -> None:
    if CSV_PATH.exists():
        rows = load_rows()
    else:
        rows, rejected = compute_rows()
        write_products(rows, rejected)
    plot(rows)


if __name__ == "__main__":
    main()
