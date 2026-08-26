#!/usr/bin/env python3
"""Publication figure for the validated homogeneous physical trajectories."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures.style import (COLORS as PALETTE, WIDE_FIGURE_WIDTH, add_panel_label,
                           save_figure, use_publication_style)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "python" / "homogeneous_physical_trajectories.csv"
FIGURES = ROOT / "paper" / "figures"

COLORS = {"tensor": "#0072B2", "vector": "#D55E00", "singlet": "#009E73"}
LABELS = {"tensor": "quintuplet", "vector": "triplet", "singlet": "singlet"}
STYLES = {"mu0": "-", "mu_over_T_2": (0, (4, 2))}


def main() -> None:
    use_publication_style()
    data = pd.read_csv(DATA)
    fig, axes = plt.subplots(2, 2, figsize=(WIDE_FIGURE_WIDTH, 5.15))
    for trajectory in ("mu0", "mu_over_T_2"):
        for sector in ("tensor", "vector", "singlet"):
            group = data[(data.trajectory == trajectory) & (data.sector == sector)].sort_values("T_MeV")
            label = LABELS[sector] if trajectory == "mu0" else None
            axes[0, 0].plot(group.T_MeV, group.Re_omega_hat,
                            color=COLORS[sector], ls=STYLES[trajectory], marker="o",
                            markevery=3, ms=2.7, lw=1.35, label=label)
            axes[0, 1].semilogy(group.T_MeV, -group.Im_omega_hat,
                               color=COLORS[sector], ls=STYLES[trajectory], marker="o",
                               markevery=3, ms=2.7, lw=1.35)
            axes[1, 0].plot(group.Re_omega_hat, -group.Im_omega_hat,
                            color=COLORS[sector], ls=STYLES[trajectory], marker="o",
                            markevery=3, ms=2.7, lw=1.35)
            tau = 197.3269804 / (2.0 * np.pi * group.T_MeV * (-group.Im_omega_hat))
            axes[1, 1].semilogy(group.T_MeV, tau, color=COLORS[sector],
                               ls=STYLES[trajectory], marker="o", markevery=3,
                               ms=2.7, lw=1.35)
            # Temperature increases away from the low-temperature endpoint.
            index = max(1, len(group) // 2)
            first = group.iloc[index - 1]
            second = group.iloc[index]
            axes[1, 0].annotate("", xy=(second.Re_omega_hat, -second.Im_omega_hat),
                                xytext=(first.Re_omega_hat, -first.Im_omega_hat),
                                arrowprops={"arrowstyle": "->", "color": COLORS[sector],
                                            "lw": 0.8, "mutation_scale": 7})
    axes[0, 0].set(xlabel=r"$T$ [MeV]", ylabel=r"$\mathrm{Re}\,\mathfrak{w}$")
    axes[0, 1].set(xlabel=r"$T$ [MeV]", ylabel=r"$-\mathrm{Im}\,\mathfrak{w}$")
    axes[1, 0].set(xlabel=r"$\mathrm{Re}\,\mathfrak{w}$",
                   ylabel=r"$-\mathrm{Im}\,\mathfrak{w}$")
    axes[1, 1].set(xlabel=r"$T$ [MeV]", ylabel=r"$\tau$ [fm/$c$]")
    axes[1, 1].text(0.04, 0.94, "long lived, still gapped", transform=axes[1, 1].transAxes,
                    ha="left", va="top", color=COLORS["vector"], fontsize=7.2)
    # Explicit line-style key for the two controlled thermodynamic paths.
    sector_handles = [plt.Line2D([], [], color=COLORS[key], label=LABELS[key])
                      for key in ("tensor", "vector", "singlet")]
    path_handles = [plt.Line2D([], [], color=PALETTE["black"], ls=STYLES["mu0"],
                              label=r"$\mu_B=0$"),
                    plt.Line2D([], [], color=PALETTE["black"], ls=STYLES["mu_over_T_2"],
                              label=r"$\mu_B/T=2$")]
    axes[0, 0].legend(handles=sector_handles, frameon=False, ncol=3)
    axes[0, 1].legend(handles=path_handles, frameon=False, ncol=2)
    for tag, ax in zip(("(a)", "(b)", "(c)", "(d)"), axes.flat):
        add_panel_label(ax, tag)
        ax.grid(axis="y", color=PALETTE["light_gray"], alpha=0.35, lw=0.45)
    fig.tight_layout(w_pad=1.15, h_pad=1.1)
    save_figure(fig, FIGURES / "homogeneous_qnm_trajectories",
                title="Validated homogeneous QNM trajectories",
                subject="Real parts, damping, complex-plane paths, and relaxation times")
    plt.close(fig)


if __name__ == "__main__":
    main()
