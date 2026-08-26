#!/usr/bin/env python3
"""Render model, background, and thermodynamic appendix diagnostics."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import run_numerics as base
from figures.style import (COLORS, WIDE_FIGURE_WIDTH, add_panel_label,
                           save_figure, use_publication_style)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"


def model_kernels() -> None:
    phi = np.linspace(0.0, 6.6, 1400)
    fig, axes = plt.subplots(1, 2, figsize=(WIDE_FIGURE_WIDTH, 2.85))
    axes[0].semilogy(phi, base.V(phi) / -12.0, color=COLORS["blue"])
    axes[0].axvspan(2.2, 6.5, color=COLORS["blue"], alpha=0.07,
                    label="production horizon range")
    axes[0].set(xlabel=r"$\phi$", ylabel=r"$-V(\phi)/12$")
    axes[0].legend(frameon=False)
    axes[1].plot(phi, base.f(phi), color=COLORS["vermillion"])
    axes[1].axvspan(2.2, 6.5, color=COLORS["blue"], alpha=0.07)
    axes[1].set(xlabel=r"$\phi$", ylabel=r"$f(\phi)$", ylim=(-0.03, 1.04))
    inset = axes[1].inset_axes([0.54, 0.42, 0.42, 0.48])
    phi_uv = np.geomspace(1.0e-6, 3.0e-2, 900)
    inset.semilogx(phi_uv, base.f(phi_uv), color=COLORS["vermillion"], lw=1.0)
    inset.axvline(1.0 / base.P.d2, color=COLORS["gray"], ls="--", lw=0.75)
    inset.set(xlabel=r"near-boundary $\phi$", ylabel=r"$f(\phi)$", ylim=(0.34, 1.03))
    inset.tick_params(labelsize=6.5)
    for tag, ax in zip(("(a)", "(b)"), axes):
        add_panel_label(ax, tag)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.33, lw=0.45)
    fig.tight_layout(w_pad=1.15)
    save_figure(fig, FIGURES / "pha_model_kernels",
                title="PHA model kernels and sampled field range",
                subject="Dilaton potential, Maxwell coupling, and resolved UV layer")
    plt.close(fig)


def thermodynamic_validation() -> None:
    data = pd.read_csv(RESULTS / "thermodynamic_validation.csv")
    fine = data[data.charge_step == data.charge_step.min()].sort_values("T0_MeV")
    fig, axes = plt.subplots(1, 3, figsize=(WIDE_FIGURE_WIDTH, 2.75))
    axes[0].plot(fine.T0_MeV, fine.s_over_T3, color=COLORS["blue"], marker="o",
                 markevery=3, mfc="white", ms=3.0)
    axes[0].set(xlabel=r"$T$ [MeV]", ylabel=r"$s/T^3$")
    step_styles = ((0.005, COLORS["blue"], "o", "-"),
                   (0.010, COLORS["vermillion"], "s", "--"),
                   (0.020, COLORS["green"], "^", ":"))
    for step, color, marker, style in step_styles:
        group = data[np.isclose(data.charge_step, step)].sort_values("T0_MeV")
        label = rf"$\delta={step:.3f}$"
        axes[1].plot(group.T0_MeV, group.chi2_B, color=color, marker=marker,
                     markevery=3, mfc="white", ms=2.8, ls=style, label=label)
        if step > 0.005:
            axes[2].semilogy(group.T0_MeV,
                             np.maximum(100.0 * group.relative_step_shift, 1e-12),
                             color=color, marker=marker, markevery=3, mfc="white",
                             ms=2.8, ls=style)
    axes[1].set(xlabel=r"$T$ [MeV]", ylabel=r"$\chi_2^B$")
    axes[2].set(xlabel=r"$T$ [MeV]", ylabel="relative displacement [%]")
    axes[1].legend(frameon=False, ncol=1)
    for tag, ax in zip(("(a)", "(b)", "(c)"), axes):
        ax.set_xscale("log")
        add_panel_label(ax, tag)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.33, lw=0.45)
    fig.tight_layout(w_pad=1.1)
    save_figure(fig, FIGURES / "thermodynamic_validation",
                title="Neutral thermodynamics and derivative stability",
                subject="Entropy, baryon susceptibility, and finite-difference displacement")
    plt.close(fig)


def representative_background() -> None:
    background, _ = base.load_independent_cusp_background()
    r = np.linspace(0.0, 10.0, 1100)
    values = base.profile_on_grid(background, r)
    A, h, phi, potential = values[:4]
    flux = np.exp(2.0 * A) * base.f(phi) * values[7]
    relative_flux = np.abs(flux / flux[0] - 1.0)
    fig, axes = plt.subplots(1, 3, figsize=(WIDE_FIGURE_WIDTH, 2.75))
    axes[0].plot(r, A / A[-1], color=COLORS["blue"], label=r"$A/A(r_{\rm cut})$")
    axes[0].plot(r, h / h[-1], color=COLORS["vermillion"], ls="--", label=r"$h/h_0$")
    axes[0].set(xlabel=r"radial coordinate $r$", ylabel="normalized metric", ylim=(-0.03, 1.05))
    axes[0].legend(frameon=False)
    axes[1].semilogy(r, np.maximum(phi, 1e-16), color=COLORS["blue"], label=r"$\phi$")
    axes[1].semilogy(r, np.maximum(1.0 - base.f(phi), 1e-16),
                    color=COLORS["vermillion"], ls="--", label=r"$1-f(\phi)$")
    axes[1].set(xlabel=r"radial coordinate $r$", ylabel="scalar / coupling layer")
    axes[1].legend(frameon=False)
    inset = axes[1].inset_axes([0.48, 0.49, 0.47, 0.42])
    inset.semilogy(r, np.maximum(1.0 - base.f(phi), 1e-16),
                   color=COLORS["vermillion"], lw=0.9)
    inset.set(xlim=(5.5, 10), ylim=(1e-16, 1e-2), xlabel=r"$r$", ylabel=r"$1-f$")
    inset.tick_params(labelsize=6.3)
    sampled = slice(None, None, 8)
    axes[2].semilogy(r[sampled], np.maximum(relative_flux[sampled], 1e-17),
                    ls="none", marker="o", ms=1.7, color=COLORS["green"],
                    alpha=0.48, label="pointwise error")
    block = 25
    centers = np.asarray([np.mean(r[index:index+block]) for index in range(0, len(r), block)])
    envelope = np.asarray([np.max(relative_flux[index:index+block])
                           for index in range(0, len(r), block)])
    axes[2].semilogy(centers, np.maximum(envelope, 1e-17),
                    color=COLORS["green"], lw=1.0, label="local maximum envelope")
    axes[2].axhline(max(background["qdrift"], 1e-17), color=COLORS["gray"], ls="--",
                    lw=0.8, label="stored Gauss drift")
    axes[2].set(xlabel=r"radial coordinate $r$", ylabel=r"$|\mathcal{Q}/\mathcal{Q}_H-1|$")
    axes[2].legend(frameon=False, fontsize=6.5)
    for tag, ax in zip(("(a)", "(b)", "(c)"), axes):
        add_panel_label(ax, tag)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.3, lw=0.45)
    fig.tight_layout(w_pad=1.1)
    save_figure(fig, FIGURES / "representative_background",
                title="Representative MAP critical background",
                subject="Metric regularity, UV coupling layer, and relative flux conservation")
    plt.close(fig)


def main() -> None:
    use_publication_style()
    model_kernels()
    thermodynamic_validation()
    representative_background()


if __name__ == "__main__":
    main()
