#!/usr/bin/env python3
"""Tabulate and plot the independent reproduction of the MAP critical point."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from figures.style import (COLORS, DOUBLE_COLUMN_WIDTH, add_panel_label,
                           save_figure, use_publication_style)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"
FILE_T = 103.89777513571676
FILE_MU = 602.4749076542938

def main() -> None:
    use_publication_style()
    rows = []
    for line_count in (100, 200, 400, 800):
        data = json.loads((RESULTS / f"reference_cep_neighbors_N{line_count}.json").read_text())
        rows.append({
            "method": "constant_phi0_crossing",
            "resolution": line_count,
            "T_c_MeV": data["T_c_MeV"],
            "mu_B_c_MeV": data["mu_B_c_MeV"],
            "delta_T_from_HDF5_MeV": data["T_c_MeV"] - FILE_T,
            "delta_mu_from_HDF5_MeV": data["mu_B_c_MeV"] - FILE_MU,
        })
    cusp = json.loads((RESULTS / "reference_cusp_uv1e5_h10.json").read_text())
    rows.append({
        "method": "local_cusp",
        "resolution": "finite_difference",
        "T_c_MeV": cusp["T_c_MeV"],
        "mu_B_c_MeV": cusp["mu_B_c_MeV"],
        "delta_T_from_HDF5_MeV": cusp["T_c_MeV"] - FILE_T,
        "delta_mu_from_HDF5_MeV": cusp["mu_B_c_MeV"] - FILE_MU,
    })
    with (RESULTS / "cep_reproduction.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    crossing = rows[:-1]
    n = np.array([row["resolution"] for row in crossing], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COLUMN_WIDTH, 2.75))
    axes[0].plot(n, [row["T_c_MeV"] - cusp["T_c_MeV"] for row in crossing], "o-", color=COLORS["blue"],
                 label=r"constant-$\phi_0$ crossings")
    axes[1].plot(n, [row["mu_B_c_MeV"] - cusp["mu_B_c_MeV"] for row in crossing], "s-", color=COLORS["blue"])
    axes[0].axhline(0, color=COLORS["green"], ls="--", label="local cusp")
    axes[1].axhline(0, color=COLORS["green"], ls="--")
    axes[0].axhline(FILE_T - cusp["T_c_MeV"], color=COLORS["gray"], ls=":", label="HDF5 offset")
    axes[1].axhline(FILE_MU - cusp["mu_B_c_MeV"], color=COLORS["gray"], ls=":")
    for tag, ax in zip(("(a)", "(b)"), axes):
        ax.set_xscale("log", base=2)
        ax.set_xticks(n, [str(int(value)) for value in n])
        ax.set_xlabel(r"Number of constant-$\phi_0$ lines")
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.35, lw=0.45)
        add_panel_label(ax, tag)
    axes[0].set_ylabel(r"$T_c(N)-T_c^{\rm local}$ [MeV]")
    axes[1].set_ylabel(r"$\mu_{B,c}(N)-\mu_{B,c}^{\rm local}$ [MeV]")
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    save_figure(fig, FIGURES / "cep_reproduction",
                title="Critical-point reconstruction convergence",
                subject="Crossing-method differences relative to the local cusp")
    plt.close(fig)


if __name__ == "__main__":
    main()
