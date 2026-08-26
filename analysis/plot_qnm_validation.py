#!/usr/bin/env python3
"""Convergence figure for independently validated homogeneous QNMs."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from figures.style import (COLORS, DOUBLE_COLUMN_WIDTH, add_panel_label,
                           save_figure, use_publication_style)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"

def main() -> None:
    use_publication_style()
    convergence = pd.read_csv(RESULTS / "qnm_factored_convergence.csv")
    accepted = pd.read_csv(RESULTS / "qnm_validated_modes.csv")[["sector", "mode"]]
    data = convergence.merge(accepted, on=["sector", "mode"])
    labels = {("tensor", 0): "quintuplet 0", ("vector", 0): "triplet 0",
              ("vector", 1): "triplet 1", ("singlet", 0): "singlet 0"}
    colors = {("tensor", 0): "#0072B2", ("vector", 0): "#D55E00",
              ("vector", 1): "#CC79A7", ("singlet", 0): "#009E73"}
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COLUMN_WIDTH, 2.82))
    for key, group in data.groupby(["sector", "mode"]):
        resolution = group[group.study == "radial_resolution"].sort_values("N")
        domain = group[group.study == "radial_domain"].sort_values("r_cut")
        axes[0].semilogy(resolution.N, 2 * resolution.distance_to_shooting,
                        color=colors[key], marker="o", ms=3, ls="none",
                        label=labels[key])
        axes[1].semilogy(domain.r_cut, 2 * domain.distance_to_shooting,
                        color=colors[key], marker="s", ms=3, ls="none")
    axes[0].set(xlabel=r"Chebyshev intervals $N$",
                ylabel=r"$|\Delta\hat\omega|$ from shooting",
                )
    axes[1].set(xlabel=r"UV cutoff $r_{\rm cut}$",
                ylabel=r"$|\Delta\hat\omega|$ from shooting",
                )
    axes[0].legend(frameon=False, ncol=2)
    for tag, ax in zip(("(a)", "(b)"), axes):
        ax.axhspan(1.0e-7, 2.0e-2, color=COLORS["green"], alpha=0.06,
                   label="accepted envelope" if tag == "(a)" else None)
        ax.axhline(1.0e-7, color=COLORS["gray"], ls="--", lw=0.75,
                   label="shooting root floor" if tag == "(a)" else None)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.38, lw=0.45)
        add_panel_label(ax, tag)
    handles, legend_labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, legend_labels, frameon=False, ncol=2)
    fig.tight_layout(w_pad=1.35)
    save_figure(fig, FIGURES / "qnm_convergence",
                title="Homogeneous QNM convergence and shooting cross-check",
                subject="Separated radial-resolution and UV-domain studies")
    plt.close(fig)


if __name__ == "__main__":
    main()
