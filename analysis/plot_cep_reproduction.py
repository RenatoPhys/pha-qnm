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


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
FIGURES = ROOT / "paper" / "figures"
FILE_T = 103.89777513571676
FILE_MU = 602.4749076542938

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIX Two Text", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def main() -> None:
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
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.7))
    axes[0].plot(n, [row["T_c_MeV"] for row in crossing], "o-", color="#0072B2",
                 label=r"constant-$\phi_0$ crossings")
    axes[1].plot(n, [row["mu_B_c_MeV"] for row in crossing], "o-", color="#0072B2")
    axes[0].axhline(cusp["T_c_MeV"], color="#009E73", ls="--", label="local cusp")
    axes[1].axhline(cusp["mu_B_c_MeV"], color="#009E73", ls="--")
    axes[0].axhline(FILE_T, color="#6B7280", ls=":", label="HDF5 metadata")
    axes[1].axhline(FILE_MU, color="#6B7280", ls=":")
    for ax in axes:
        ax.set_xscale("log", base=2)
        ax.set_xticks(n, [str(int(value)) for value in n])
        ax.set_xlabel(r"Number of constant-$\phi_0$ lines")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"$T_c$ [MeV]")
    axes[1].set_ylabel(r"$\mu_{B,c}$ [MeV]")
    axes[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / "cep_reproduction.pdf")
    fig.savefig(FIGURES / "cep_reproduction.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
