#!/usr/bin/env python3
"""Render the stable, critical, and spinodal finite-momentum figures from data."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

import charged_hydrodynamics as hydro
from figures.style import (COLORS, WIDE_FIGURE_WIDTH, add_panel_label,
                           save_figure, use_publication_style)


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"


def hydrodynamic_critical() -> None:
    stable = pd.read_csv(RESULTS / "finite_k_hydrodynamic_dispersion.csv")
    critical = pd.read_csv(RESULTS / "cep_critical_scaling.csv")
    summary = json.loads((RESULTS / "finite_k_physics_summary.json").read_text(encoding="utf-8"))
    colors = {"neutral": COLORS["blue"], "charged": COLORS["vermillion"]}
    markers = {"shear": "o", "sound": "s", "diffusion": "^"}
    fig, axes = plt.subplots(2, 2, figsize=(WIDE_FIGURE_WIDTH, 5.0),
                             gridspec_kw={"height_ratios": [1.0, 0.82]})
    residual_records = []
    for trajectory, color in colors.items():
        state = summary["stable"][trajectory]
        reference_row = stable[stable.trajectory == trajectory].iloc[0]
        hydro_state = hydro.charged_hydro_state(
            float(reference_row.phi0), float(reference_row.charge_fraction),
            backend="reference_tight")
        probe_modes = hydro.hydrodynamic_modes(hydro_state, 1.0e-4)
        sound_probe = probe_modes[np.argmax(probe_modes.real)]
        independent_gamma = float(-sound_probe.imag / 1.0e-8)
        sound = stable[(stable.trajectory == trajectory) & (stable["mode"] == "sound")]
        qline = np.linspace(sound.qhat.min(), sound.qhat.max(), 160)
        axes[0, 0].plot(qline, state["c_s"] * qline, color=color, lw=1.5)
        axes[0, 0].plot(sound.qhat, sound.omega_hat_real, ls="none", marker="s",
                        mfc="white", mec=color, mew=0.9)
        for mode, marker in markers.items():
            group = stable[(stable.trajectory == trajectory) & (stable["mode"] == mode)]
            coefficient = {"shear": state["D_eta_hat"], "diffusion": state["D_B_hat"],
                           "sound": independent_gamma}[mode]
            axes[0, 1].plot(qline, coefficient * qline**2, color=color,
                            ls={"shear": "-", "sound": "--", "diffusion": ":"}[mode], lw=1.35)
            damping = -group.omega_hat_imag.to_numpy()
            axes[0, 1].plot(group.qhat, damping, ls="none", marker=marker,
                            mfc="white", mec=color, mew=0.9)
            prediction = coefficient * group.qhat.to_numpy()**2
            residual_records.append((trajectory, mode, group.qhat.to_numpy(),
                                     100.0 * (damping - prediction) / prediction))
    axes[0, 0].set(xlabel=r"$\mathfrak{q}$", ylabel=r"$\mathrm{Re}\,\mathfrak{w}_+$")
    axes[0, 1].set(xlabel=r"$\mathfrak{q}$", ylabel=r"$-\mathrm{Im}\,\mathfrak{w}$")
    state_handles = [Line2D([], [], color=color, lw=1.5, label=name)
                     for name, color in colors.items()]
    source_handles = [Line2D([], [], marker="o", mfc="white", mec=COLORS["black"],
                             ls="none", label="QNM"),
                      Line2D([], [], color=COLORS["black"], lw=1.4,
                             label="independent hydro")]
    axes[0, 0].legend(handles=state_handles + source_handles, frameon=False, ncol=2)
    mode_handles = [Line2D([], [], marker=marker, mfc="white", mec=COLORS["black"],
                           ls={"shear": "-", "sound": "--", "diffusion": ":"}[mode],
                           color=COLORS["black"], label=mode)
                    for mode, marker in markers.items()]
    axes[0, 1].legend(handles=mode_handles, frameon=False, ncol=3)

    for trajectory, mode, q, residual in residual_records:
        axes[1, 0].plot(q, residual, ls="none", marker=markers[mode], ms=3.0,
                        mfc="white", mec=colors[trajectory], alpha=0.9)
    axes[1, 0].axhline(0, color=COLORS["black"], lw=0.7)
    axes[1, 0].set(xlabel=r"$\mathfrak{q}$", ylabel="QNM-hydro difference [%]")
    axes[1, 0].set_ylim(-5, 5)

    x = critical.chi_B_over_T2.to_numpy(); y = critical.D_B_hat.to_numpy()
    order = np.argsort(x)
    axes[1, 1].loglog(x, y, "o", mfc="white", mec=COLORS["green"], mew=0.9,
                      label="thermodynamic path")
    fit_rows = critical.iloc[-6:]
    slope, intercept = np.polyfit(np.log(fit_rows.chi_B_over_T2),
                                  np.log(fit_rows.D_B_hat), 1)
    xfit = np.geomspace(fit_rows.chi_B_over_T2.min(), fit_rows.chi_B_over_T2.max(), 100)
    axes[1, 1].loglog(xfit, np.exp(intercept) * xfit**slope,
                      color=COLORS["green"], lw=1.4,
                      label=rf"fit: $z={2-2*slope:.2f}$")
    axes[1, 1].set(xlabel=r"$\chi_B/T^2$", ylabel=r"$D_B$")
    axes[1, 1].legend(frameon=False)
    for label, ax in zip(("(a)", "(b)", "(c)", "(d)"), axes.flat):
        add_panel_label(ax, label)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.33, lw=0.45)
    fig.tight_layout(w_pad=1.1, h_pad=1.15)
    save_figure(fig, FIGURES / "hydrodynamic_critical_dispersion",
                title="Stable hydrodynamics and critical scaling",
                subject="Independent hydrodynamic curves, QNM points, residuals, and D_B-chi_B fit")
    plt.close(fig)


def cep_spinodal() -> None:
    cep = pd.read_csv(RESULTS / "cep_q4_dispersion.csv").sort_values("qhat")
    spinodal = pd.read_csv(RESULTS / "spinodal_dispersion.csv")
    summary = json.loads((RESULTS / "finite_k_physics_summary.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(1, 2, figsize=(WIDE_FIGURE_WIDTH, 3.25))
    q = cep.qhat.to_numpy(); damping = -cep.omega_hat_imag.to_numpy()
    fit_q, fit_y = q[:4], damping[:4]
    exponent, intercept = np.polyfit(np.log(fit_q), np.log(fit_y), 1)
    qline = np.geomspace(q.min(), q.max(), 200)
    axes[0].loglog(q, damping, "o", mfc="white", mec=COLORS["green"], mew=0.9,
                   label="full QNM")
    axes[0].loglog(qline, np.exp(intercept) * qline**exponent,
                   color=COLORS["green"], lw=1.5,
                   label=rf"fit window: $z={exponent:.2f}$")
    amplitude = fit_y[0] / fit_q[0]**4
    axes[0].loglog(qline, amplitude * qline**4, color=COLORS["gray"], ls="--",
                   lw=1.0, label=r"$\mathfrak{q}^4$ reference")
    axes[0].axvspan(fit_q.min(), fit_q.max(), color=COLORS["green"], alpha=0.08)
    axes[0].set(xlabel=r"$\mathfrak{q}$", ylabel=r"$-\mathrm{Im}\,\mathfrak{w}_D$")
    axes[0].legend(frameon=False)
    inset = axes[0].inset_axes([0.56, 0.16, 0.40, 0.32])
    local = summary["cep_q4"]["local_exponents"]
    inset.plot([row["qhat_geometric_mean"] for row in local],
               [row["z_local"] for row in local], "o-", ms=2.3,
               color=COLORS["green"], lw=0.9)
    inset.axhline(4, color=COLORS["gray"], ls="--", lw=0.7)
    inset.set(xlabel=r"$\mathfrak{q}$", ylabel=r"$z_{\rm eff}$", ylim=(3.88, 4.02))
    inset.tick_params(labelsize=6.5)

    trajectories = (("f0.10", "near hot fold", COLORS["blue"], "o"),
                    ("f0.50", "midpoint", COLORS["vermillion"], "s"),
                    ("f0.90", "near cold fold", COLORS["purple"], "^"))
    ymax = spinodal.omega_hat_imag.max() * 1.18
    axes[1].axhspan(0, ymax, color=COLORS["vermillion"], alpha=0.07)
    axes[1].text(0.49, ymax * 0.90, "unstable", color=COLORS["vermillion"], ha="right")
    for key, label, color, marker in trajectories:
        group = spinodal[spinodal.trajectory == key].sort_values("qhat")
        qg, growth = group.qhat.to_numpy(), group.omega_hat_imag.to_numpy()
        axes[1].plot(qg, growth, ls="none", marker=marker, mfc="white", mec=color,
                     mew=0.9, label=label)
        dense = np.linspace(qg.min(), qg.max(), 250)
        axes[1].plot(dense, CubicSpline(qg, growth)(dense), color=color, lw=1.05,
                     alpha=0.8)
        axes[1].plot(qg, group.hydro_omega_hat_imag, color=color, ls="--", lw=1.1)
    middle = summary["spinodal"]["f0.50"]
    axes[1].plot(middle["q_star"], middle["Gamma_star_hat"], marker="*", ms=7,
                 color=COLORS["vermillion"], zorder=8)
    axes[1].axvline(middle["q_edge"], color=COLORS["vermillion"], ls=":", lw=0.9)
    axes[1].annotate(r"$(\mathfrak{q}_*,\Gamma_*)$",
                     (middle["q_star"], middle["Gamma_star_hat"]),
                     xytext=(8, 7), textcoords="offset points", fontsize=7.2)
    axes[1].text(middle["q_edge"] + 0.008, -0.00023, r"$\mathfrak{q}_{\rm edge}$",
                 color=COLORS["vermillion"], fontsize=7.1, rotation=90, va="bottom")
    axes[1].axhline(0, color=COLORS["black"], lw=0.7)
    axes[1].set(xlabel=r"$\mathfrak{q}$", ylabel=r"$\mathrm{Im}\,\mathfrak{w}_D$",
                ylim=(-0.00035, ymax))
    state_handles = [Line2D([], [], marker=marker, ls="none", mfc="white", mec=color,
                            color=color, label=label)
                     for _, label, color, marker in trajectories]
    axes[1].legend(handles=state_handles, frameon=False, fontsize=6.8, ncol=3,
                   loc="upper center", bbox_to_anchor=(0.5, 1.16))
    axes[1].text(0.02, 0.04, "markers: full QNM\nsolid: interpolation\ndashed: charged hydro",
                 transform=axes[1].transAxes, ha="left", va="bottom", fontsize=6.7,
                 color=COLORS["gray"])
    for label, ax in zip(("(a)", "(b)"), axes):
        add_panel_label(ax, label)
        ax.grid(axis="y", color=COLORS["light_gray"], alpha=0.33, lw=0.45)
    fig.tight_layout(w_pad=1.25)
    save_figure(fig, FIGURES / "cep_spinodal_dynamics",
                title="Critical and spinodal longitudinal dynamics",
                subject="Critical q4 dispersion and the closed diffusive unstable band")
    plt.close(fig)


def main() -> None:
    use_publication_style()
    hydrodynamic_critical()
    cep_spinodal()


if __name__ == "__main__":
    main()
