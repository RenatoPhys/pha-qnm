#!/usr/bin/env python3
"""Reproducible Python analysis for the PHA background and derived QNM sectors."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import eig
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"
RESULTS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans",
    "font.size": 8.5,
    "axes.labelsize": 9.0,
    "axes.titlesize": 9.0,
    "axes.titleweight": "semibold",
    "axes.titlelocation": "left",
    "legend.fontsize": 7.5,
    "figure.dpi": 180,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.axisbelow": True,
    "axes.grid": True,
    "grid.color": "#DCE3E9",
    "grid.alpha": 0.8,
    "grid.linewidth": 0.5,
    "axes.edgecolor": "#9AA7B2",
    "axes.linewidth": 0.7,
    "xtick.color": "#536273",
    "ytick.color": "#536273",
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "lines.linewidth": 1.8,
    "lines.markersize": 4.0,
})

BLUE = "#0072B2"
VERMILLION = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
INK = "#1F2933"
MID = "#64707D"


def style_axes(axes):
    """Apply the common publication treatment without touching plot data."""
    for ax in np.atleast_1d(axes).flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(width=0.7)
        ax.title.set_color(INK)
        ax.xaxis.label.set_color(INK)
        ax.yaxis.label.set_color(INK)


def panel_title(ax, tag, title):
    ax.set_title(f"{tag}  {title}", loc="left", pad=5.5, color=INK,
                 fontweight="semibold")


@dataclass(frozen=True)
class Parameters:
    Lambda_MeV: float = 1123.8704743960677
    kappa2: float = 11.37047249736345
    gamma: float = 0.5904238801815029
    b2: float = 0.32722698958661095
    b4: float = -0.04801342482759139
    b6: float = 0.0009385148172188441
    c1: float = 0.008800483500598342
    c2: float = 0.15579934656052294
    c3: float = 0.05281851198008199
    d1: float = 1.7177450982418923
    d2: float = 1679.1777562624359


P = Parameters()
M2 = -12.0 * P.gamma**2 + 2.0 * P.b2
DELTA = 2.0 + math.sqrt(4.0 + M2)
NU = 4.0 - DELTA


def sech(x):
    x = np.asarray(x, dtype=float)
    a = np.abs(x)
    e = np.exp(-a)
    return 2.0 * e / (1.0 + e * e)


def V(phi):
    return (-12.0 * np.cosh(P.gamma * phi) + P.b2 * phi**2
            + P.b4 * phi**4 + P.b6 * phi**6)


def Vp(phi):
    return (-12.0 * P.gamma * np.sinh(P.gamma * phi) + 2.0 * P.b2 * phi
            + 4.0 * P.b4 * phi**3 + 6.0 * P.b6 * phi**5)


def Vpp(phi):
    return (-12.0 * P.gamma**2 * np.cosh(P.gamma * phi) + 2.0 * P.b2
            + 12.0 * P.b4 * phi**2 + 30.0 * P.b6 * phi**4)


def g(phi):
    return P.c1 * phi + P.c2 * phi**2 + P.c3 * phi**3


def gp(phi):
    return P.c1 + 2.0 * P.c2 * phi + 3.0 * P.c3 * phi**2


def gpp(phi):
    return 2.0 * P.c2 + 6.0 * P.c3 * phi


def f(phi):
    return (sech(g(phi)) + P.d1 * sech(P.d2 * phi)) / (1.0 + P.d1)


def fp(phi):
    sg = sech(g(phi))
    sd = sech(P.d2 * phi)
    return (-sg * np.tanh(g(phi)) * gp(phi)
            - P.d1 * P.d2 * sd * np.tanh(P.d2 * phi)) / (1.0 + P.d1)


def fpp(phi):
    gg = g(phi)
    sg = sech(gg)
    tg = np.tanh(gg)
    sd = sech(P.d2 * phi)
    td = np.tanh(P.d2 * phi)
    term_g = sg * ((tg * tg - sg * sg) * gp(phi)**2 - tg * gpp(phi))
    term_d = P.d1 * P.d2**2 * sd * (td * td - sd * sd)
    return (term_g + term_d) / (1.0 + P.d1)


def horizon_data(phi0: float, phi1_electric: float):
    f0, fp0 = float(f(phi0)), float(fp(phi0))
    vmax = math.sqrt(-2.0 * float(V(phi0)) / f0)
    if abs(phi1_electric) >= vmax:
        raise ValueError("extremal or super-extremal horizon datum")
    A1 = -(2.0 * float(V(phi0)) + f0 * phi1_electric**2) / 6.0
    scalar1 = float(Vp(phi0)) - 0.5 * fp0 * phi1_electric**2
    A2 = -scalar1**2 / 12.0
    h2 = 0.5 * (f0 * phi1_electric**2 - 4.0 * A1)
    Phi2 = -0.5 * phi1_electric * (2.0 * A1 + fp0 * scalar1 / f0)
    electro1 = ((-2.0 * A1 * fp0 + float(fpp(phi0)) * scalar1)
                * phi1_electric**2 + 4.0 * fp0 * phi1_electric * Phi2)
    scalar2 = (float(Vpp(phi0)) * scalar1 - (4.0 * A1 + 2.0 * h2) * scalar1
               - 0.5 * electro1) / 4.0
    return {
        "vmax": vmax, "A1": A1, "scalar1": scalar1, "A2": A2,
        "h2": h2, "Phi2": Phi2, "scalar2": scalar2,
    }


def integrate_background(phi0: float, charge_fraction: float, r_max: float = 12.0):
    vmax = math.sqrt(-2.0 * float(V(phi0)) / float(f(phi0)))
    Phi1 = charge_fraction * vmax
    hd = horizon_data(phi0, Phi1)
    eps = 1.0e-7
    conserved_charge = float(f(phi0)) * Phi1
    y0 = np.array([
        hd["A1"] * eps + hd["A2"] * eps**2,
        eps + hd["h2"] * eps**2,
        phi0 + hd["scalar1"] * eps + hd["scalar2"] * eps**2,
        Phi1 * eps + hd["Phi2"] * eps**2,
        hd["A1"] + 2.0 * hd["A2"] * eps,
        1.0 + 2.0 * hd["h2"] * eps,
        hd["scalar1"] + 2.0 * hd["scalar2"] * eps,
    ])

    def rhs(_r, y):
        A, h, phi, _Phi, Ap, hp, phip = y
        fv = float(f(phi))
        em2A = math.exp(-2.0 * A)
        Phip = conserved_charge * em2A / fv
        return np.array([
            Ap, hp, phip, Phip,
            -phip * phip / 6.0,
            -4.0 * Ap * hp + em2A * fv * Phip * Phip,
            -(4.0 * Ap + hp / h) * phip
            + (float(Vp(phi)) - 0.5 * em2A * float(fp(phi)) * Phip * Phip) / h,
        ])

    sol = solve_ivp(rhs, (eps, r_max), y0, method="DOP853", rtol=3.0e-12,
                    atol=3.0e-14, dense_output=True, max_step=0.025)
    if not sol.success:
        raise RuntimeError(sol.message)
    r = sol.t
    raw = sol.y
    Phip = conserved_charge * np.exp(-2.0 * raw[0]) / f(raw[2])
    y = np.vstack((raw, Phip))
    constraint = (y[1] * (24.0 * y[4]**2 - y[6]**2) + 6.0 * y[4] * y[5]
                  + 2.0 * V(y[2]) + np.exp(-2.0 * y[0]) * f(y[2]) * y[7]**2)
    scale = 1.0 + np.abs(2.0 * V(y[2])) + np.abs(6.0 * y[4] * y[5])
    q = np.exp(2.0 * y[0]) * f(y[2]) * y[7]
    qdrift = np.max(np.abs(q / conserved_charge - 1.0)) if abs(conserved_charge) > 1e-14 else 0.0
    return {
        "phi0": phi0, "charge_fraction": charge_fraction, "Phi1": Phi1,
        "r_max": r_max, "solution": sol, "r": r, "y": y,
        "charge": conserved_charge,
        "constraint": float(np.max(np.abs(constraint) / scale)),
        "qdrift": float(qdrift), "horizon": hd,
    }


def background_state(bg, r):
    raw = bg["solution"].sol(r)
    Phip = bg["charge"] * np.exp(-2.0 * raw[0]) / f(raw[2])
    if raw.ndim == 1:
        return np.append(raw, Phip)
    return np.vstack((raw, Phip))


def extract_uv(bg):
    # UV fits must not inherit the nonuniform point density selected by the
    # adaptive ODE driver.  A fixed sampling grid makes the inverse map from
    # (phi_0, Phi_1/Phi_1^max) to (T, mu_B) smooth enough for root finding.
    r = np.linspace(float(bg["r"][0]), float(bg["r_max"]), 2400)
    y = background_state(bg, r)
    A, h, phi, Phi = y[0], y[1], y[2], y[3]
    mask = (phi < 1.0e-4) & (phi > 1.0e-10)
    if np.count_nonzero(mask) < 30:
        mask = (A > A[-1] - 8.0) & (phi > 2.0e-13)
    rr, AA = r[mask], A[mask]
    slope, intercept = np.polyfit(rr, AA, 1)
    alpha = slope * r + intercept
    h0 = float(np.median(h[-max(20, r.size // 10):]))
    phiA_samples = phi[mask] * np.exp(NU * alpha[mask])
    phiA = float(np.median(phiA_samples))
    charge = np.exp(2.0 * A) * f(phi) * y[7]
    charge_uv = float(np.median(charge[mask]))
    Phi2 = -charge_uv / (2.0 * slope)
    Phi0_samples = Phi[mask] - Phi2 * np.exp(-2.0 * alpha[mask])
    Phi0 = float(np.median(Phi0_samples))
    T = P.Lambda_MeV / (4.0 * math.pi * phiA**(1.0 / NU) * math.sqrt(h0))
    mu = Phi0 * P.Lambda_MeV / (phiA**(1.0 / NU) * math.sqrt(h0))
    s_over_T3 = 128.0 * math.pi**4 * h0**1.5 / P.kappa2
    rho_over_T3 = -64.0 * math.pi**3 * Phi2 * h0 / P.kappa2
    return {
        "A_slope": float(slope), "A_intercept": float(intercept), "h0": h0,
        "phiA": phiA, "Phi0": float(Phi0), "Phi2": float(Phi2),
        "T_MeV": float(T), "mu_MeV": float(mu), "mu_over_T": float(mu / T),
        "s_over_T3": float(s_over_T3), "rho_over_T3": float(rho_over_T3),
        "slope_identity_error": float(abs(slope * math.sqrt(h0) - 1.0)),
        "phiA_relative_spread": float(np.std(phiA_samples) / abs(np.mean(phiA_samples))),
    }


def chebyshev_lobatto(intervals: int, lower: float, upper: float):
    j = np.arange(intervals + 1)
    x = -np.cos(np.pi * j / intervals)
    nodes = 0.5 * ((upper - lower) * x + upper + lower)
    bary = np.where((j == 0) | (j == intervals), 0.5, 1.0) * np.where(j % 2 == 0, 1.0, -1.0)
    D = np.zeros((intervals + 1, intervals + 1))
    for i in range(intervals + 1):
        for k in range(intervals + 1):
            if i != k:
                D[i, k] = bary[k] / (bary[i] * (nodes[i] - nodes[k]))
        D[i, i] = -np.sum(D[i])
    return nodes, D, D @ D


def profile_on_grid(bg, nodes):
    phi0, Phi1 = bg["phi0"], bg["Phi1"]
    hd = bg["horizon"]
    values = np.empty((8, nodes.size))
    positive = nodes > 0.0
    values[:, positive] = background_state(bg, nodes[positive])
    values[:, 0] = [0.0, 0.0, phi0, 0.0, hd["A1"], 1.0, hd["scalar1"], Phi1]
    return values


def assemble_pencil(bg, intervals: int, r_cut: float, sector: str, q_value: float = 0.0):
    nodes, D1, D2 = chebyshev_lobatto(intervals, 0.0, r_cut)
    y = profile_on_grid(bg, nodes)
    A, h, phi, _Phi, Ap, hp, phip, Phip = y
    fv, fpv = f(phi), fp(phi)
    n = nodes.size
    M0 = np.zeros((n, n), dtype=complex)
    M1 = np.zeros((n, n), dtype=complex)
    emA = np.exp(-A)
    if sector == "tensor":
        h0 = extract_uv(bg)["h0"]
        k_numeric = q_value / (2.0 * math.sqrt(h0))
        for i in range(n):
            M0[i] = h[i] * D2[i] + (4.0 * h[i] * Ap[i] + hp[i]) * D1[i]
            M0[i, i] -= emA[i]**2 * k_numeric**2
            M1[i] = -2.0j * emA[i] * D1[i]
            M1[i, i] += -3.0j * emA[i] * Ap[i]
    elif sector == "vector":
        ratio = fpv * phip / fv
        for i in range(n):
            M0[i] = h[i] * D2[i] + (2.0 * h[i] * Ap[i] + hp[i] + h[i] * ratio[i]) * D1[i]
            M0[i, i] -= emA[i]**2 * fv[i] * Phip[i]**2
            M1[i] = -2.0j * emA[i] * D1[i]
            M1[i, i] += -1.0j * emA[i] * (Ap[i] + ratio[i])
    else:
        raise ValueError(sector)
    M0[-1] = 0.0
    M1[-1] = 0.0
    M0[-1, -1] = 1.0
    return M0, M1


def qnm_spectrum(bg, intervals: int, r_cut: float, sector: str, q_value: float = 0.0):
    M0, M1 = assemble_pencil(bg, intervals, r_cut, sector, q_value)
    omega, vr = eig(M0, -M1, right=True, check_finite=False)
    finite = np.isfinite(omega) & (np.abs(omega) < 80.0) & (omega.imag < 1.0e-7)
    omega, vr = omega[finite], vr[:, finite]
    residuals = []
    for z, v in zip(omega, vr.T):
        matrix = M0 + z * M1
        residuals.append(np.linalg.norm(matrix @ v) /
                         (max(np.linalg.norm(matrix, ord=np.inf) * np.linalg.norm(v), 1e-300)))
    residuals = np.asarray(residuals)
    order = np.lexsort((np.abs(omega.real), np.abs(omega.imag)))
    return omega[order], residuals[order]


def converged_modes(bg, sector: str, q_value: float = 0.0, count: int = 6):
    settings = [(104, 9.0), (112, 10.0), (128, 10.0), (144, 11.0)]
    spectra = [qnm_spectrum(bg, n, rc, sector, q_value) for n, rc in settings]
    reference, reference_res = spectra[-1]
    accepted = []
    for idx, z in enumerate(reference):
        if z.real < -1.0e-7 or z.imag >= -1.0e-5 or reference_res[idx] > 2.0e-10:
            continue
        distances = []
        for omega, _res in spectra[:-1]:
            distances.append(float(np.min(np.abs(omega - z))))
        if max(distances) < 8.0e-3:
            accepted.append({
                "omega": z, "residual": float(reference_res[idx]),
                "resolution_spread": max(distances),
            })
        if len(accepted) >= count:
            break
    return accepted, spectra


def save_figure(fig, stem):
    metadata = {
        "Title": stem.replace("_", " ").title(),
        "Author": "PHA QNM collaboration",
        "Subject": "Reproducible figure generated from analysis/run_numerics.py",
    }
    fig.savefig(FIGURES / f"{stem}.pdf", metadata=metadata)
    fig.savefig(FIGURES / f"{stem}.png", metadata={"Software": "Matplotlib"})
    plt.close(fig)


def plot_potentials():
    phi = np.linspace(0.0, 6.0, 1200)
    phi_uv = np.geomspace(1.0e-6, 3.0e-2, 900)
    fig, axes = plt.subplots(1, 3, figsize=(7.35, 2.68))
    axes[0].plot(phi, V(phi) / -12.0, color=BLUE)
    axes[0].set(xlabel=r"Scalar field $\phi$", ylabel=r"$-V(\phi)/12$")
    axes[0].set_yscale("log")
    panel_title(axes[0], "(a)", "Dilaton potential")
    axes[1].plot(phi, f(phi), color=VERMILLION)
    axes[1].set(xlabel=r"Scalar field $\phi$", ylabel=r"$f(\phi)$", ylim=(-0.03, 1.04))
    panel_title(axes[1], "(b)", "Maxwell coupling")
    axes[2].semilogx(phi_uv, f(phi_uv), color=VERMILLION)
    axes[2].axvline(1.0 / P.d2, color=MID, ls=(0, (4, 2)), lw=1.0)
    axes[2].annotate(r"$d_2^{-1}$", xy=(1.0 / P.d2, 0.78), xytext=(7, 5),
                     textcoords="offset points", color=MID, fontsize=7.5)
    axes[2].set(xlabel=r"Near-boundary $\phi$", ylabel=r"$f(\phi)$", ylim=(0.34, 1.03))
    panel_title(axes[2], "(c)", "UV boundary layer")
    style_axes(axes)
    fig.tight_layout(w_pad=1.25)
    save_figure(fig, "pha_model_kernels")


def scan_backgrounds():
    phi_values = np.geomspace(0.18, 5.5, 30)
    charge_fractions = [0.0, 0.18, 0.32]
    rows, backgrounds = [], {}
    for cf in charge_fractions:
        for phi0 in phi_values:
            bg = integrate_background(float(phi0), cf)
            uv = extract_uv(bg)
            row = {"phi0": float(phi0), "charge_fraction": cf, "Phi1": bg["Phi1"],
                   "constraint": bg["constraint"], "gauss_drift": bg["qdrift"], **uv}
            rows.append(row)
            backgrounds[(round(float(phi0), 8), cf)] = bg
    with (RESULTS / "background_scan.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    return rows, backgrounds


def plot_background_scan(rows, critical):
    fig, axes = plt.subplots(1, 3, figsize=(7.35, 2.92))
    colors = [BLUE, VERMILLION, GREEN]
    for color, cf in zip(colors, [0.0, 0.18, 0.32]):
        group = [r for r in rows if r["charge_fraction"] == cf]
        T = np.array([r["T_MeV"] for r in group])
        axes[0].plot([r["phi0"] for r in group], T, color=color,
                     label=fr"$\Phi_1/\Phi_1^{{\rm max}}={cf:.2f}$")
        axes[1].plot(T, [r["s_over_T3"] for r in group], color=color)
        if cf > 0:
            axes[2].plot(T, [r["mu_over_T"] for r in group], color=color)
    axes[0].set(xscale="log", yscale="log", xlabel=r"Horizon field $\phi_0$", ylabel=r"$T$ [MeV]")
    axes[1].set(xscale="log", xlabel=r"$T$ [MeV]", ylabel=r"$s/T^3$")
    axes[2].set(xscale="log", yscale="log", xlabel=r"$T$ [MeV]", ylabel=r"$\mu_B/T$")
    panel_title(axes[0], "(a)", "Horizon-data scan")
    panel_title(axes[1], "(b)", "Dimensionless entropy")
    panel_title(axes[2], "(c)", "Charged slices")
    axes[0].scatter([critical["phi0"]], [critical["T_MeV"]], marker="*", s=60,
                    facecolor="#F0E442", edgecolor=INK, linewidth=0.7, zorder=5)
    axes[1].scatter([critical["T_MeV"]], [critical["s_over_T3"]], marker="*", s=60,
                    facecolor="#F0E442", edgecolor=INK, linewidth=0.7, zorder=5)
    axes[2].scatter([critical["T_MeV"]], [critical["mu_over_T"]], marker="*", s=60,
                    facecolor="#F0E442", edgecolor=INK, linewidth=0.7, zorder=5)
    axes[0].set_xticks([0.2, 0.5, 1.0, 2.0, 5.0], ["0.2", "0.5", "1", "2", "5"])
    for ax in axes[1:]:
        ax.set_xticks([50.0, 100.0, 200.0, 500.0, 1000.0],
                      ["50", "100", "200", "500", "1000"])
    for ax in axes:
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.015), handlelength=2.2, columnspacing=1.3)
    style_axes(axes)
    fig.tight_layout(rect=(0, 0, 1, 0.90), w_pad=1.25)
    save_figure(fig, "background_thermodynamics")


def thermodynamic_validation(rows):
    """Reproduce neutral thermodynamics and test the mu_B derivative defining chi_2."""
    neutral = sorted((row for row in rows if row["charge_fraction"] == 0.0),
                     key=lambda row: row["T_MeV"])
    charge_steps = (0.005, 0.010, 0.020)
    records = []
    for base in neutral:
        for step in charge_steps:
            plus_bg = integrate_background(base["phi0"], step)
            minus_bg = integrate_background(base["phi0"], -step)
            plus, minus = extract_uv(plus_bg), extract_uv(minus_bg)
            delta_x = plus["mu_over_T"] - minus["mu_over_T"]
            chi2 = (plus["rho_over_T3"] - minus["rho_over_T3"]) / delta_x
            records.append({
                "phi0": base["phi0"], "T0_MeV": base["T_MeV"],
                "s_over_T3": base["s_over_T3"], "charge_step": step,
                "mu_over_T_plus": plus["mu_over_T"],
                "mu_over_T_minus": minus["mu_over_T"],
                "chi2_B": chi2,
                "max_constraint": max(plus_bg["constraint"], minus_bg["constraint"]),
                "max_gauss_drift": max(plus_bg["qdrift"], minus_bg["qdrift"]),
            })

    finest = {row["phi0"]: row["chi2_B"] for row in records
              if row["charge_step"] == charge_steps[0]}
    for row in records:
        reference = finest[row["phi0"]]
        row["relative_step_shift"] = abs(row["chi2_B"] - reference) / abs(reference)

    with (RESULTS / "thermodynamic_validation.csv").open(
            "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)

    fig, axes = plt.subplots(1, 3, figsize=(7.35, 2.75))
    temperature = np.array([row["T_MeV"] for row in neutral])
    axes[0].plot(temperature, [row["s_over_T3"] for row in neutral],
                 color=BLUE, marker="o", ms=3.0)
    axes[0].set(xlabel=r"$T$ [MeV]", ylabel=r"$s/T^3$")
    panel_title(axes[0], "(a)", "Neutral entropy")

    step_colors = (BLUE, VERMILLION, GREEN)
    for step, color in zip(charge_steps, step_colors):
        group = sorted((row for row in records if row["charge_step"] == step),
                       key=lambda row: row["T0_MeV"])
        label = fr"$\delta={step:.3f}$"
        axes[1].plot([row["T0_MeV"] for row in group],
                     [row["chi2_B"] for row in group],
                     color=color, marker="o", ms=2.7, label=label)
        if step > charge_steps[0]:
            axes[2].semilogy([row["T0_MeV"] for row in group],
                             np.maximum([row["relative_step_shift"] for row in group], 1e-16),
                             color=color, marker="o", ms=2.7, label=label)
    axes[1].set(xlabel=r"$T$ [MeV]", ylabel=r"$\chi_2^B$")
    axes[2].set(xlabel=r"$T$ [MeV]", ylabel="Relative step shift")
    panel_title(axes[1], "(b)", "Baryon susceptibility")
    panel_title(axes[2], "(c)", "Finite-difference stability")
    axes[1].legend(frameon=False, loc="best")
    axes[2].legend(frameon=False, loc="best")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks([50.0, 100.0, 200.0, 500.0, 1000.0],
                      ["50", "100", "200", "500", "1000"])
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    style_axes(axes)
    fig.tight_layout(w_pad=1.25)
    save_figure(fig, "thermodynamic_validation")

    fine_rows = [row for row in records if row["charge_step"] == charge_steps[0]]
    return records, {
        "temperature_range_MeV": [min(temperature), max(temperature)],
        "entropy_range": [min(row["s_over_T3"] for row in neutral),
                          max(row["s_over_T3"] for row in neutral)],
        "chi2_B_range": [min(row["chi2_B"] for row in fine_rows),
                         max(row["chi2_B"] for row in fine_rows)],
        "max_relative_step_shift": max(row["relative_step_shift"] for row in records),
        "max_constraint": max(row["max_constraint"] for row in records),
        "max_gauss_drift": max(row["max_gauss_drift"] for row in records),
        "charge_steps": list(charge_steps),
    }


def representative_background_plot(bg):
    r = np.linspace(0.0, 10.0, 900)
    y = profile_on_grid(bg, r)
    A, h, phi, Phi = y[:4]
    q = np.exp(2.0 * A) * f(phi) * y[7]
    fig, axes = plt.subplots(1, 3, figsize=(7.35, 2.70))
    axes[0].plot(r, A / A[-1], color=BLUE, label=r"$A/A(r_{\rm cut})$")
    axes[0].plot(r, h / h[-1], color=VERMILLION, label=r"$h/h_0$")
    axes[0].set(xlabel=r"Radial coordinate $r$", ylim=(-0.03, 1.05))
    panel_title(axes[0], "(a)", "Metric background")
    axes[0].legend(frameon=False, loc="lower right")
    axes[1].semilogy(r, np.maximum(phi, 1e-16), color=BLUE, label=r"$\phi$")
    axes[1].semilogy(r, np.maximum(1.0 - f(phi), 1e-16), color=VERMILLION,
                    label=r"$1-f(\phi)$")
    axes[1].set(xlabel=r"Radial coordinate $r$", ylim=(5e-17, 12))
    panel_title(axes[1], "(b)", "Scalar and UV layer")
    axes[1].legend(frameon=False, loc="upper right")
    axes[2].plot(r, Phi, color=BLUE, label=r"$\Phi$")
    axes[2].plot(r, q / q[0], color=VERMILLION, label=r"$\mathcal{Q}/\mathcal{Q}_H$")
    axes[2].set(xlabel=r"Radial coordinate $r$", ylim=(-0.03, 1.06))
    panel_title(axes[2], "(c)", "Maxwell sector")
    axes[2].legend(frameon=False, loc="lower right")
    style_axes(axes)
    fig.tight_layout(w_pad=1.25)
    save_figure(fig, "representative_background")


def locate_critical_background():
    target_T = 103.89777513571676
    target_mu = 602.4749076542938
    cache = {}

    def evaluate(x):
        key = tuple(np.round(x, 11))
        if key not in cache:
            bg = integrate_background(float(x[0]), float(x[1]))
            cache[key] = (bg, extract_uv(bg))
        return cache[key]

    def residual(x):
        _bg, uv = evaluate(x)
        return [math.log(uv["T_MeV"] / target_T), math.log(uv["mu_MeV"] / target_mu)]

    fit = least_squares(residual, x0=[3.706, 0.3146],
                        bounds=([2.5, 0.05], [6.5, 0.45]),
                        x_scale=[4.0, 0.25], diff_step=2.0e-3,
                        xtol=2e-11, ftol=2e-11, gtol=2e-11, max_nfev=60)
    bg, uv = evaluate(fit.x)
    # Condition the local inverse map using logarithmic horizon coordinates,
    # so both columns of the Jacobian are dimensionless and comparable.
    jac = np.empty((2, 2))
    rel_step = 5.0e-4
    for column in range(2):
        xp, xm = fit.x.copy(), fit.x.copy()
        xp[column] *= math.exp(rel_step)
        xm[column] *= math.exp(-rel_step)
        up, um = evaluate(xp)[1], evaluate(xm)[1]
        jac[:, column] = [math.log(up["T_MeV"] / um["T_MeV"]),
                          math.log(up["mu_MeV"] / um["mu_MeV"])]
        jac[:, column] /= 2.0 * rel_step
    singular_values = np.linalg.svd(jac, compute_uv=False)
    return bg, {"phi0": float(fit.x[0]), "charge_fraction": float(fit.x[1]),
                "target_T_MeV": target_T, "target_mu_MeV": target_mu,
                "delta_T_MeV": float(uv["T_MeV"] - target_T),
                "delta_mu_MeV": float(uv["mu_MeV"] - target_mu),
                "relative_delta_T": float(uv["T_MeV"] / target_T - 1.0),
                "relative_delta_mu": float(uv["mu_MeV"] / target_mu - 1.0),
                "fit_residual_norm": float(np.linalg.norm(fit.fun)),
                "inverse_map_singular_values": singular_values.tolist(),
                "inverse_map_condition_number": float(singular_values[0] / singular_values[-1]),
                **uv}


def run_qnms(bg):
    records = []
    accepted_by_sector = {}
    spectra_by_sector = {}
    for sector in ("tensor", "vector"):
        accepted, spectra = converged_modes(bg, sector)
        accepted_by_sector[sector] = accepted
        spectra_by_sector[sector] = spectra
        for mode, item in enumerate(accepted):
            z = item["omega"]
            records.append({"sector": sector, "mode": mode,
                            "Re_omega_num": z.real, "Im_omega_num": z.imag,
                            "Re_w": 2.0 * z.real, "Im_w": 2.0 * z.imag,
                            "residual": item["residual"],
                            "resolution_spread": item["resolution_spread"]})

    q_grid = np.linspace(0.0, 3.0, 16)
    # Only the two best-converged tensor poles are continued.  The third
    # homogeneous pole is retained in the table but is too close to the
    # acceptance threshold for a finite-momentum claim.
    base = accepted_by_sector["tensor"][:2]
    previous = [x["omega"] for x in base]
    dispersion = []
    for qv in q_grid:
        omega, residual = qnm_spectrum(bg, 128, 10.0, "tensor", float(qv))
        candidates = [(z, residual[i]) for i, z in enumerate(omega)
                      if z.real >= -1e-7 and z.imag < -1e-5 and residual[i] < 2e-10]
        used = set()
        tracked = []
        for mode, target in enumerate(previous):
            choices = [(abs(z - target), i, z, res) for i, (z, res) in enumerate(candidates)
                       if i not in used]
            if not choices:
                continue
            _dist, idx, z, res = min(choices)
            used.add(idx); tracked.append(z)
            dispersion.append({"q": float(qv), "mode": mode, "Re_w": 2.0*z.real,
                               "Im_w": 2.0*z.imag, "residual": float(res)})
        previous = tracked

    with (RESULTS / "qnm_modes.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    with (RESULTS / "tensor_dispersion.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(dispersion[0]))
        writer.writeheader(); writer.writerows(dispersion)

    fig, axes = plt.subplots(1, 3, figsize=(7.35, 2.75))
    markers = {"tensor": "o", "vector": "s"}
    colors = {"tensor": BLUE, "vector": VERMILLION}
    for sector in ("tensor", "vector"):
        items = accepted_by_sector[sector]
        z = 2.0 * np.array([x["omega"] for x in items])
        axes[0].scatter(z.real, z.imag, label=sector, marker=markers[sector],
                        color=colors[sector], s=28)
    axes[0].axhline(0.0, color="0.2", lw=0.7)
    axes[0].set(xlabel=r"$\mathrm{Re}\,\mathfrak{w}$",
                ylabel=r"$\mathrm{Im}\,\mathfrak{w}$")
    panel_title(axes[0], "(a)", r"Converged $k=0$ poles")
    axes[0].legend(frameon=False, loc="lower left")

    mode_colors = [BLUE, VERMILLION, GREEN]
    for mode in range(3):
        group = [x for x in dispersion if x["mode"] == mode]
        if not group:
            continue
        axes[1].plot([x["q"] for x in group], [x["Re_w"] for x in group],
                     marker="o", ms=3.0, color=mode_colors[mode], label=fr"$n={mode}$")
        axes[2].plot([x["q"] for x in group], [-x["Im_w"] for x in group],
                     marker="o", ms=3.0, color=mode_colors[mode])
    axes[1].set(xlabel=r"Momentum $\mathfrak{q}$", ylabel=r"$\mathrm{Re}\,\mathfrak{w}$")
    axes[2].set(xlabel=r"Momentum $\mathfrak{q}$", ylabel=r"$-\mathrm{Im}\,\mathfrak{w}$")
    panel_title(axes[1], "(b)", "Tensor dispersion")
    panel_title(axes[2], "(c)", "Tensor damping")
    axes[1].legend(frameon=False, loc="upper left")
    style_axes(axes)
    fig.tight_layout(w_pad=1.25)
    save_figure(fig, "qnm_spectra_dispersion")

    fig, axes = plt.subplots(1, 2, figsize=(6.45, 2.72))
    for ax, sector in zip(axes, ("tensor", "vector")):
        reference = accepted_by_sector[sector][0]["omega"]
        Ns, errors, residuals = [], [], []
        for (n, _rc), (omega, res) in zip([(104,9.0),(112,10.0),(128,10.0),(144,11.0)],
                                          spectra_by_sector[sector]):
            idx = int(np.argmin(np.abs(omega - reference)))
            Ns.append(n); errors.append(abs(omega[idx] - reference)); residuals.append(res[idx])
        errors[-1] = np.nan
        x = np.arange(len(Ns))
        ax.semilogy(x, np.maximum(errors, 1e-17), "o-", color=BLUE,
                    markerfacecolor="white", markeredgewidth=1.1, label="mode shift")
        ax.semilogy(x, np.maximum(residuals, 1e-17), "s--", color=VERMILLION,
                    markerfacecolor="white", markeredgewidth=1.1, label="pencil residual")
        ax.set_xticks(x, ["104/9", "112/10", "128/10", "144/11"])
        ax.set(xlabel=r"Resolution $N/r_{\rm cut}$", ylim=(3e-17, 3e-4))
    axes[0].set_ylabel("Absolute diagnostic")
    panel_title(axes[0], "(a)", "Tensor leading pole")
    panel_title(axes[1], "(b)", "Vector leading pole")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.015), handlelength=2.4)
    style_axes(axes)
    fig.tight_layout(rect=(0, 0, 1, 0.90), w_pad=1.25)
    save_figure(fig, "qnm_convergence")
    return records, dispersion


def qnm_temperature_scan(charge_fraction=0.18):
    phi_values = np.geomspace(0.55, 5.5, 11)
    rows = []
    for phi0 in phi_values:
        bg = integrate_background(float(phi0), charge_fraction)
        uv = extract_uv(bg)
        for sector in ("tensor", "vector"):
            accepted, _spectra = converged_modes(bg, sector, count=1)
            if not accepted:
                continue
            z = accepted[0]["omega"]
            rows.append({"phi0": float(phi0), "charge_fraction": charge_fraction,
                         "sector": sector, "T_MeV": uv["T_MeV"],
                         "mu_over_T": uv["mu_over_T"], "Re_w": 2.0*z.real,
                         "Im_w": 2.0*z.imag,
                         "resolution_spread": accepted[0]["resolution_spread"]})
    with (RESULTS / "qnm_temperature_scan.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    fig, axes = plt.subplots(1, 2, figsize=(6.45, 2.72))
    for sector, color, marker in (("tensor", BLUE, "o"), ("vector", VERMILLION, "s")):
        group = sorted((x for x in rows if x["sector"] == sector), key=lambda x: x["T_MeV"])
        axes[0].plot([x["T_MeV"] for x in group], [x["Re_w"] for x in group],
                     marker=marker, ms=3.5, color=color, label=sector)
        axes[1].plot([x["T_MeV"] for x in group], [-x["Im_w"] for x in group],
                     marker=marker, ms=3.5, color=color, label=sector)
    for ax in axes:
        ax.set_xscale("log"); ax.set_xlabel(r"$T$ [MeV]")
        ax.set_xticks([50.0, 100.0, 200.0, 500.0, 1000.0],
                      ["50", "100", "200", "500", "1000"])
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    axes[0].set_ylabel(r"$\mathrm{Re}\,\mathfrak{w}_0$")
    axes[1].set_ylabel(r"$-\mathrm{Im}\,\mathfrak{w}_0$")
    panel_title(axes[0], "(a)", "Leading oscillation")
    panel_title(axes[1], "(b)", "Leading damping")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.015), handlelength=2.4)
    style_axes(axes)
    fig.tight_layout(rect=(0, 0, 1, 0.90), w_pad=1.25)
    save_figure(fig, "qnm_temperature_scan")
    return rows


def main():
    plot_potentials()
    rows, _backgrounds = scan_backgrounds()
    representative, critical = locate_critical_background()
    plot_background_scan(rows, critical)
    thermodynamic_rows, thermodynamic_diagnostics = thermodynamic_validation(rows)
    representative_uv = extract_uv(representative)
    representative_background_plot(representative)
    modes, dispersion = run_qnms(representative)
    temperature_modes = qnm_temperature_scan()
    summary = {
        "parameters": asdict(P), "m2": M2, "Delta": DELTA, "nu": NU,
        "background_grid_points": len(rows),
        "representative": {"phi0": critical["phi0"],
                           "charge_fraction": critical["charge_fraction"],
                           "Phi1": representative["Phi1"],
                           "constraint": representative["constraint"],
                           "gauss_drift": representative["qdrift"], **representative_uv},
        "nearest_cep_background": critical,
        "pseudospectral_qnm_candidates": modes,
        "tensor_dispersion_points": len(dispersion),
        "temperature_scan_points": len(temperature_modes),
        "thermodynamic_validation_points": len(thermodynamic_rows),
        "thermodynamic_validation": thermodynamic_diagnostics,
        "scope": "tensor finite-k and homogeneous vector only; no scalar/helicity-1/helicity-0 claims",
    }
    with (RESULTS / "summary.json").open("w", encoding="utf-8") as out:
        json.dump(summary, out, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
