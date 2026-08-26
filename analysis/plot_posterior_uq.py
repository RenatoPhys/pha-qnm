#!/usr/bin/env python3
"""Plot the controlled posterior propagation used by the Route-A paper."""

from __future__ import annotations

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

def main() -> None:
    use_publication_style()
    rows = json.loads((RESULTS / "posterior_uq_samples.json").read_text(encoding="utf-8"))
    summary = json.loads((RESULTS / "posterior_uq_summary.json").read_text(encoding="utf-8"))
    rows = [row for row in rows if row["status"] == "success"]
    weights = np.asarray([row["posterior_weight"] for row in rows])
    sizes = 22.0 + 420.0 * weights
    likelihood = np.asarray([row["log_likelihood"] for row in rows])
    fig, axes = plt.subplots(1, 3, figsize=(WIDE_FIGURE_WIDTH, 3.05))
    norm = plt.Normalize(likelihood.min(), likelihood.max())
    cmap = plt.get_cmap("viridis")

    def weighted_median(values):
        values = np.asarray(values); order = np.argsort(values)
        return float(values[order][np.searchsorted(np.cumsum(weights[order]), 0.5)])

    scatter = axes[0].scatter(
        [row["cusp"]["mu_B_c_MeV"] for row in rows],
        [row["cusp"]["T_c_MeV"] for row in rows],
        c=likelihood, s=sizes, cmap=cmap, norm=norm, edgecolor="white", linewidth=0.45,
    )
    map_row = next(row for row in rows if row["sample_index"] == 74)
    axes[0].scatter(map_row["cusp"]["mu_B_c_MeV"], map_row["cusp"]["T_c_MeV"],
                    marker="*", s=105, color=COLORS["vermillion"], edgecolor="white",
                    linewidth=0.6, label="MAP")
    axes[0].scatter(weighted_median([row["cusp"]["mu_B_c_MeV"] for row in rows]),
                    weighted_median([row["cusp"]["T_c_MeV"] for row in rows]),
                    marker="D", s=30, facecolor="white", edgecolor=COLORS["black"],
                    linewidth=0.8, label="weighted median")
    axes[0].set(xlabel=r"$\mu_B^c$ [MeV]", ylabel=r"$T_c$ [MeV]")
    axes[0].legend(frameon=False, fontsize=8)
    order = np.argsort([row["critical_fit"]["z_eta0"] for row in rows])
    z = np.asarray([rows[index]["critical_fit"]["z_eta0"] for index in order])
    zerr = np.asarray([
        rows[index]["critical_fit"]["z_window_half_range"] for index in order
    ])
    rank = np.arange(1, len(rows) + 1)
    for x, value, error in zip(rank, z, zerr):
        axes[1].plot([x, x], [value-error, value+error], color=COLORS["gray"], lw=0.65)
    axes[1].scatter(rank, z, c=likelihood[order], cmap=cmap, norm=norm,
                    s=sizes[order], edgecolor="white", linewidth=0.4, zorder=4)
    interval = summary["credible_intervals"]["z_eta0"]
    axes[1].axhspan(interval["q2p5"], interval["q97p5"], color=COLORS["green"],
                    alpha=0.15, label="95% posterior")
    axes[1].axhline(interval["median"], color=COLORS["green"], lw=1.1,
                    label="weighted median")
    axes[1].axhline(interval["median"] - 0.05, color=COLORS["gray"], ls="--", lw=0.8)
    axes[1].axhline(interval["median"] + 0.05, color=COLORS["gray"], ls="--", lw=0.8,
                    label=r"common window $\pm0.05$")
    axes[1].set(xlabel="medoid rank in $z$", ylabel=r"$z$")
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].scatter(
        [row["spinodal"]["spinodal_temperature_width_MeV"] for row in rows],
        [row["spinodal"]["midpoint_D_B_hat"] for row in rows],
        c=likelihood, cmap=cmap, norm=norm, s=sizes, edgecolor="white", linewidth=0.45,
    )
    xmin, xmax = axes[2].get_xlim()
    ymin = min(row["spinodal"]["midpoint_D_B_hat"] for row in rows) * 1.08
    axes[2].axhspan(ymin, 0.0, color=COLORS["vermillion"], alpha=0.07)
    axes[2].axhline(0.0, color=COLORS["black"], lw=0.8)
    axes[2].scatter(weighted_median([row["spinodal"]["spinodal_temperature_width_MeV"] for row in rows]),
                    weighted_median([row["spinodal"]["midpoint_D_B_hat"] for row in rows]),
                    marker="D", s=30, facecolor="white", edgecolor=COLORS["black"], zorder=5)
    axes[2].text(0.97, 0.05, "25/25 medoids unstable", transform=axes[2].transAxes,
                 ha="right", va="bottom", color=COLORS["vermillion"], fontsize=7.1)
    axes[2].set(xlabel=r"$\Delta T_{\rm sp}$ [MeV] at $\mu_c+50$ MeV",
                ylabel=r"$D_B$ at spinodal midpoint")

    for tag, ax in zip(("(a)", "(b)", "(c)"), axes):
        add_panel_label(ax, tag)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.30, lw=0.45)
    colorbar = fig.colorbar(scatter, ax=axes, pad=0.012, fraction=0.025)
    colorbar.set_label(r"$\log\mathcal{L}$", fontsize=8)
    fig.subplots_adjust(left=0.08, right=0.91, bottom=0.18, top=0.97, wspace=0.36)
    save_figure(fig, FIGURES / "posterior_uq", title="Posterior robustness",
                subject="Weighted CEPs, dynamic exponent, and spinodal instability")
    plt.close(fig)


if __name__ == "__main__":
    main()
