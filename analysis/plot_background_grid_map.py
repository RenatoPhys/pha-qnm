#!/usr/bin/env python3
"""Map the computed horizon-data grid into thermodynamic coordinates."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from figures.style import (COLORS, WIDE_FIGURE_WIDTH, add_panel_label,
                           save_figure, use_publication_style)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"


def audit_points() -> list[dict[str, float | str]]:
    path = RESULTS / "longitudinal_mode_identity_audit.csv"
    if not path.exists():
        return []
    points = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("record_type") != "hydrodynamic_eigenmode":
                continue
            points.setdefault(row["background"], {
                "background": row["background"], "phi0": float(row["phi0"]),
                "charge_fraction": float(row["charge_fraction"]),
                "T_MeV": float(row["T_MeV"]), "mu_B_MeV": float(row["mu_B_MeV"]),
            })
    return list(points.values())


def main() -> None:
    use_publication_style()
    surface = np.load(RESULTS / "phase_surface.npz")
    phi, charge = surface["phi0"], surface["charge_fraction"]
    temperature, chemical = surface["T"], surface["mu"]
    determinant = surface["determinant"]
    PHI, CHARGE = np.meshgrid(phi, charge, indexing="ij")
    valid = np.isfinite(temperature) & np.isfinite(chemical)
    phase_lines = json.loads((RESULTS / "phase_lines.json").read_text(encoding="utf-8"))
    dynamic = audit_points()

    fig, axes = plt.subplots(1, 2, figsize=(WIDE_FIGURE_WIDTH, 3.45))
    ax = axes[0]
    ax.scatter(PHI[valid][::5], CHARGE[valid][::5], s=2.0, color=COLORS["sky"],
               alpha=0.24, linewidths=0, rasterized=True, label="successful grid")
    failed = ~valid
    ax.scatter(PHI[failed][::8], CHARGE[failed][::8], s=4.0, marker="x",
               color=COLORS["light_gray"], alpha=0.7, linewidths=0.45,
               rasterized=True, label="failed")
    fold = ax.contour(phi, charge, determinant.T, levels=[0.0],
                      colors=[COLORS["vermillion"]], linewidths=1.25)
    mu_contours = ax.contour(phi, charge, chemical.T,
                             levels=[600, 650, 750], colors=[COLORS["gray"]],
                             linewidths=0.8, linestyles=[":", "-.", ":"])
    ax.clabel(mu_contours, fmt=lambda value: rf"$\mu_B={value:.0f}$", fontsize=6.8)
    ax.scatter(3.528280966, 0.346373565, marker="*", s=70,
               color=COLORS["green"], edgecolor="white", linewidth=0.5,
               zorder=8, label="CEP")
    for point in dynamic:
        if point["background"] in {"near_hot_fold", "midpoint", "near_cold_fold"}:
            ax.scatter(point["phi0"], point["charge_fraction"], s=24,
                       facecolor="white", edgecolor=COLORS["vermillion"], zorder=7)
    ax.plot([], [], color=COLORS["vermillion"], label=r"$\det J=0$")
    ax.set(xlabel=r"$\phi_0$", ylabel=r"$\Phi_1/\Phi_1^{\max}$",
           xlim=(phi.min(), phi.max()), ylim=(charge.min(), charge.max()))
    ax.legend(frameon=False, ncol=2, loc="upper right")
    add_panel_label(ax, "(a)")

    ax = axes[1]
    for i in range(0, len(phi), 14):
        mask = valid[i]
        ax.plot(chemical[i, mask], temperature[i, mask], color=COLORS["sky"],
                alpha=0.40, lw=0.65)
    for j in range(0, len(charge), 15):
        mask = valid[:, j]
        ax.plot(chemical[mask, j], temperature[mask, j], color=COLORS["orange"],
                alpha=0.45, lw=0.65, ls="--")
    mu = np.asarray([row["mu_MeV"] for row in phase_lines])
    low = np.asarray([row["spinodal_low_T_MeV"] for row in phase_lines])
    high = np.asarray([row["spinodal_high_T_MeV"] for row in phase_lines])
    coex = [row for row in phase_lines if row["coexistence"]]
    ax.plot(mu, low, color=COLORS["vermillion"], lw=1.2)
    ax.plot(mu, high, color=COLORS["vermillion"], lw=1.2, label="fold image")
    ax.plot([row["mu_MeV"] for row in coex], [row["coexistence"]["T_MeV"] for row in coex],
            color=COLORS["blue"], lw=1.55, label="coexistence")
    ax.scatter(593.322656, 103.755281, marker="*", s=70, color=COLORS["green"],
               edgecolor="white", linewidth=0.5, zorder=8)
    for point in dynamic:
        if point["background"] in {"near_hot_fold", "midpoint", "near_cold_fold"}:
            ax.scatter(point["mu_B_MeV"], point["T_MeV"], s=24, facecolor="white",
                       edgecolor=COLORS["vermillion"], zorder=7)
    ax.plot([], [], color=COLORS["sky"], lw=0.8, label=r"$\phi_0={\rm const.}$")
    ax.plot([], [], color=COLORS["orange"], lw=0.8, ls="--",
            label=r"$\Phi_1/\Phi_1^{\max}={\rm const.}$")
    ax.set(xlabel=r"$\mu_B$ [MeV]", ylabel=r"$T$ [MeV]", xlim=(0, 900), ylim=(65, 190))
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.text(665, 95, "multivalued image", color=COLORS["vermillion"], fontsize=7.2)
    add_panel_label(ax, "(b)")
    for axis in axes:
        axis.grid(axis="y", color=COLORS["light_gray"], alpha=0.33, lw=0.45)
    fig.tight_layout(w_pad=1.15)
    save_figure(fig, FIGURES / "background_grid_map",
                title="Horizon-data grid and thermodynamic map",
                subject="Successful backgrounds, map folds, and dynamic trajectories")
    plt.close(fig)


if __name__ == "__main__":
    main()
