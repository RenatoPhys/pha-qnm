#!/usr/bin/env python3
"""Convergence figure for independently validated homogeneous QNMs."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8.2, "axes.labelsize": 8.8,
    "axes.titlesize": 8.8, "legend.fontsize": 7, "pdf.fonttype": 42,
    "ps.fonttype": 42, "figure.dpi": 180, "savefig.dpi": 300,
})


def main() -> None:
    convergence = pd.read_csv(RESULTS / "qnm_factored_convergence.csv")
    accepted = pd.read_csv(RESULTS / "qnm_validated_modes.csv")[["sector", "mode"]]
    data = convergence.merge(accepted, on=["sector", "mode"])
    labels = {("tensor", 0): "quintuplet 0", ("vector", 0): "triplet 0",
              ("vector", 1): "triplet 1", ("singlet", 0): "singlet 0"}
    colors = {("tensor", 0): "#0072B2", ("vector", 0): "#D55E00",
              ("vector", 1): "#CC79A7", ("singlet", 0): "#009E73"}
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 2.82))
    for key, group in data.groupby(["sector", "mode"]):
        resolution = group[group.study == "radial_resolution"].sort_values("N")
        domain = group[group.study == "radial_domain"].sort_values("r_cut")
        axes[0].semilogy(resolution.N, 2 * resolution.distance_to_shooting,
                        color=colors[key], marker="o", ms=3, lw=1.3,
                        label=labels[key])
        axes[1].semilogy(domain.r_cut, 2 * domain.distance_to_shooting,
                        color=colors[key], marker="s", ms=3, lw=1.3)
    axes[0].set(xlabel=r"Chebyshev intervals $N$",
                ylabel=r"$|\Delta\hat\omega|$ from shooting",
                title="(a)  Fixed UV domain")
    axes[1].set(xlabel=r"UV cutoff $r_{\rm cut}$",
                ylabel=r"$|\Delta\hat\omega|$ from shooting",
                title=r"(b)  Fixed $N=160$")
    axes[0].legend(frameon=False, ncol=2)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#D7DCE0", alpha=0.38, lw=0.45)
    fig.tight_layout(w_pad=1.35)
    fig.savefig(FIGURES / "qnm_convergence.pdf", bbox_inches="tight", pad_inches=0.04,
                metadata={"Title": "Homogeneous QNM convergence and shooting cross-check",
                          "Author": "PHA QNM collaboration"})
    fig.savefig(FIGURES / "qnm_convergence.png", dpi=300, bbox_inches="tight",
                pad_inches=0.04)
    plt.close(fig)


if __name__ == "__main__":
    main()
