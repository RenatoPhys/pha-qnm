#!/usr/bin/env python3
"""Plot the controlled posterior propagation used by the Route-A paper."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def main() -> None:
    rows = json.loads((RESULTS / "posterior_uq_samples.json").read_text(encoding="utf-8"))
    summary = json.loads((RESULTS / "posterior_uq_summary.json").read_text(encoding="utf-8"))
    rows = [row for row in rows if row["status"] == "success"]
    weights = np.asarray([row["posterior_weight"] for row in rows])
    sizes = 22.0 + 420.0 * weights
    likelihood = np.asarray([row["log_likelihood"] for row in rows])
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.25))

    scatter = axes[0].scatter(
        [row["cusp"]["mu_B_c_MeV"] for row in rows],
        [row["cusp"]["T_c_MeV"] for row in rows],
        c=likelihood, s=sizes, cmap="viridis", edgecolor="white", linewidth=0.45,
    )
    map_row = next(row for row in rows if row["sample_index"] == 74)
    axes[0].scatter(map_row["cusp"]["mu_B_c_MeV"], map_row["cusp"]["T_c_MeV"],
                    marker="*", s=115, color="#D55E00", edgecolor="white",
                    linewidth=0.6, label="MAP")
    axes[0].set(xlabel=r"$\mu_B^c$ [MeV]", ylabel=r"$T_c$ [MeV]")
    axes[0].legend(frameon=False, fontsize=8)
    colorbar = fig.colorbar(scatter, ax=axes[0], pad=0.02, fraction=0.05)
    colorbar.set_label(r"$\log\mathcal{L}$", fontsize=8)

    order = np.argsort([row["critical_fit"]["z_eta0"] for row in rows])
    z = np.asarray([rows[index]["critical_fit"]["z_eta0"] for index in order])
    zerr = np.asarray([
        rows[index]["critical_fit"]["z_window_half_range"] for index in order
    ])
    cumulative = np.cumsum(weights[order]) - 0.5 * weights[order]
    axes[1].errorbar(cumulative, z, yerr=zerr, fmt="o", ms=3.2,
                     color="#0072B2", ecolor="#56B4E9", elinewidth=0.8,
                     capsize=1.5)
    interval = summary["credible_intervals"]["z_eta0"]
    axes[1].axhspan(interval["q2p5"], interval["q97p5"], color="#009E73",
                    alpha=0.15, label="95% posterior")
    axes[1].axhline(interval["median"], color="#009E73", lw=1.1)
    axes[1].set(xlabel="cumulative posterior weight", ylabel=r"$z$")
    axes[1].legend(frameon=False, fontsize=8)

    axes[2].scatter(
        [row["spinodal"]["spinodal_temperature_width_MeV"] for row in rows],
        [row["spinodal"]["midpoint_D_B_hat"] for row in rows],
        color="#0072B2", s=sizes, edgecolor="white", linewidth=0.45,
    )
    axes[2].axhline(0.0, color="0.3", lw=0.8)
    axes[2].set(xlabel=r"$\Delta T_{\rm sp}$ [MeV] at $\mu_c+50$ MeV",
                ylabel=r"$D_B$ at spinodal midpoint")

    for tag, ax in zip(("a", "b", "c"), axes):
        ax.set_title(f"({tag})", loc="left", fontweight="normal")
    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "png"):
        fig.savefig(FIGURES / f"posterior_uq.{extension}", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
