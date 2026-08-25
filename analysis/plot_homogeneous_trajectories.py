#!/usr/bin/env python3
"""Publication figure for the validated homogeneous physical trajectories."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "python" / "homogeneous_physical_trajectories.csv"
FIGURES = ROOT / "paper" / "figures"

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 9,
    "legend.fontsize": 7.3, "figure.dpi": 180, "savefig.dpi": 300,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.7,
    "xtick.direction": "out", "ytick.direction": "out",
})

COLORS = {"tensor": "#0072B2", "vector": "#D55E00", "singlet": "#009E73"}
LABELS = {"tensor": "quintuplet", "vector": "triplet", "singlet": "singlet"}
STYLES = {"mu0": "-", "mu_over_T_2": (0, (4, 2))}


def main() -> None:
    data = pd.read_csv(DATA)
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.02))
    for trajectory in ("mu0", "mu_over_T_2"):
        for sector in ("tensor", "vector", "singlet"):
            group = data[(data.trajectory == trajectory) & (data.sector == sector)].sort_values("T_MeV")
            label = LABELS[sector] if trajectory == "mu0" else None
            axes[0].plot(group.Re_omega_hat, -group.Im_omega_hat,
                         color=COLORS[sector], ls=STYLES[trajectory], marker="o",
                         ms=2.7, lw=1.35, label=label)
            axes[1].semilogy(group.T_MeV, -group.Im_omega_hat,
                            color=COLORS[sector], ls=STYLES[trajectory], marker="o",
                            ms=2.7, lw=1.35)
            # Temperature increases away from the low-temperature endpoint.
            index = max(1, len(group) // 2)
            first = group.iloc[index - 1]
            second = group.iloc[index]
            axes[0].annotate("", xy=(second.Re_omega_hat, -second.Im_omega_hat),
                             xytext=(first.Re_omega_hat, -first.Im_omega_hat),
                             arrowprops={"arrowstyle": "->", "color": COLORS[sector],
                                         "lw": 0.8, "mutation_scale": 7})
    axes[0].set(xlabel=r"$\mathrm{Re}\,\hat\omega$",
                ylabel=r"$-\mathrm{Im}\,\hat\omega$",
                title="(a)  Pole trajectories")
    axes[0].legend(frameon=False, loc="best")
    axes[1].set(xlabel=r"$T\;[\mathrm{MeV}]$",
                ylabel=r"$-\mathrm{Im}\,\hat\omega$",
                title="(b)  Dimensionless damping rate")
    # Explicit line-style key for the two controlled thermodynamic paths.
    axes[1].plot([], [], color="#39424E", ls=STYLES["mu0"], label=r"$\mu_B=0$")
    axes[1].plot([], [], color="#39424E", ls=STYLES["mu_over_T_2"],
                 label=r"$\mu_B/T=2$")
    axes[1].legend(frameon=False, loc="best")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#D7DCE0", alpha=0.35, lw=0.45)
    fig.tight_layout(w_pad=1.25)
    metadata = {"Title": "Validated homogeneous QNM trajectories",
                "Author": "PHA QNM collaboration"}
    fig.savefig(FIGURES / "homogeneous_qnm_trajectories.pdf", metadata=metadata,
                bbox_inches="tight", pad_inches=0.04)
    fig.savefig(FIGURES / "homogeneous_qnm_trajectories.png", dpi=300,
                bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


if __name__ == "__main__":
    main()
