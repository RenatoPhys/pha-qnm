#!/usr/bin/env python3
"""Quantitative homogeneous-QNM comparison with Phys. Rev. D 98 (2018).

The earlier EMD potential, Maxwell coupling, scale, and nine rounded
equilibration times are taken from arXiv:1804.00189, eq. (3) and the text below
figs. 3, 5, and 6.  The old model is solved anew with the present background,
source-factored collocation, and complex-shooting code.  This separates a
solver benchmark from the physical comparison with the PHA MAP realization.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

import run_numerics as base
from pha_qnm.thermodynamics import (
    BackgroundNumerics,
    PHAParameters,
    PHAModel,
    locate_cusp_critical_point,
)
from validate_decoupled_qnms import factored_spectrum, nearest_mode, shooting_root


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"
HBARC_MEV_FM = 197.3269804
SECTORS = ("tensor", "vector", "singlet")

SOURCE_URL = "https://arxiv.org/abs/1804.00189"
SOURCE_PDF_SHA256 = "f3a69353443faf13e02a9eed9ae43aabf6e7eb79e6aac3f6bc3a06980cffc55e"

OLD_PARAMETERS = PHAParameters(
    Lambda_MeV=1058.83,
    kappa2=8.0 * math.pi * 0.46,
    gamma=0.63,
    b2=0.65,
    b4=-0.05,
    b6=0.003,
    c1=-0.27,
    c2=0.4,
    c3=0.0,
    d1=1.7,
    d2=100.0,
)

# Values explicitly described as examples read from the plots in the source.
# Their displayed precision, not the present solver, sets their uncertainty.
PUBLISHED_TAU_FM_C = {
    "T400_mu0": {"tensor": 0.06, "vector": 0.08, "singlet": 0.11},
    "T145_mu0": {"tensor": 0.26, "vector": 1.00, "singlet": 0.70},
    "model_CEP": {"tensor": 0.33, "vector": 2.15, "singlet": 1.26},
}


def _activate(parameters: PHAParameters) -> None:
    """Select model parameters used by the shared background/QNM routines."""
    base.P = base.Parameters(**parameters.__dict__)
    base.M2 = -12.0 * base.P.gamma**2 + 2.0 * base.P.b2
    base.DELTA = 2.0 + math.sqrt(4.0 + base.M2)
    base.NU = 4.0 - base.DELTA


def _mu0_background(target_temperature: float) -> tuple[dict[str, object], dict[str, float]]:
    cache: dict[float, tuple[dict[str, object], dict[str, float]]] = {}

    def evaluate(phi0: float) -> float:
        key = round(float(phi0), 13)
        if key not in cache:
            background = base.integrate_background(float(phi0), 0.0)
            cache[key] = (background, base.extract_uv(background))
        return float(cache[key][1]["T_MeV"]) - target_temperature

    samples = np.linspace(0.15, 6.0, 80)
    bracket = None
    previous = float(samples[0])
    previous_value = evaluate(previous)
    for current in samples[1:]:
        current = float(current)
        current_value = evaluate(current)
        if previous_value * current_value <= 0.0:
            bracket = (previous, current)
            break
        previous, previous_value = current, current_value
    if bracket is None:
        raise RuntimeError(f"No neutral background brackets T={target_temperature} MeV")
    phi0 = float(brentq(evaluate, *bracket, xtol=1.0e-12, rtol=1.0e-12))
    background = base.integrate_background(phi0, 0.0)
    return background, base.extract_uv(background)


def _leading_mode(background: dict[str, object], sector: str) -> dict[str, float]:
    spectrum = factored_spectrum(background, 160, 8.0, sector)
    candidates = sorted(
        (value for value in spectrum["omega"]
         if value.real > 0.05 and value.imag < 0.0 and abs(value) < 8.0),
        key=lambda value: abs(value.imag),
    )
    for candidate in candidates[:8]:
        shot = shooting_root(
            background,
            sector,
            candidate,
            acceptance_reference=candidate,
            acceptance_radius=0.20,
            search_radius=0.20,
        )
        if not shot["success"]:
            continue
        spectral_index = nearest_mode(spectrum, shot["omega"])
        omega_hat = 2.0 * shot["omega"]
        return {
            "Re_omega_hat": float(omega_hat.real),
            "Im_omega_hat": float(omega_hat.imag),
            "shooting_source_residual": float(shot["source_residual"]),
            "spectral_distance_hat": float(
                2.0 * abs(spectrum["omega"][spectral_index] - shot["omega"])
            ),
            "spectral_pencil_residual": float(spectrum["residual"][spectral_index]),
        }
    raise RuntimeError(f"No independently accepted leading {sector} mode")


def _row(model: str, landmark: str, sector: str, background: dict[str, object],
         uv: dict[str, float], published_tau: float | None) -> dict[str, object]:
    mode = _leading_mode(background, sector)
    damping = -float(mode["Im_omega_hat"])
    tau = HBARC_MEV_FM / (2.0 * math.pi * float(uv["T_MeV"]) * damping)
    return {
        "model": model,
        "landmark": landmark,
        "sector": sector,
        "T_MeV": float(uv["T_MeV"]),
        "mu_B_MeV": float(uv["mu_MeV"]),
        "phi_H": float(background["phi0"]),
        "charge_fraction": float(background["charge_fraction"]),
        **mode,
        "tau_linear_fm_c": float(tau),
        "published_2018_tau_fm_c": published_tau,
        "delta_from_published_fm_c": (
            float(tau - published_tau) if published_tau is not None else math.nan
        ),
    }


def _pha_rows() -> list[dict[str, object]]:
    _activate(PHAParameters())
    landmarks = []
    for name, temperature in (("T400_mu0", 400.0), ("T145_mu0", 145.0)):
        background, uv = _mu0_background(temperature)
        landmarks.append((name, background, uv))
    critical_background, critical_uv = base.load_independent_cusp_background()
    landmarks.append(("model_CEP", critical_background, critical_uv))
    return [
        _row("PHA_MAP", landmark, sector, background, uv, None)
        for landmark, background, uv in landmarks
        for sector in SECTORS
    ]


def _old_cusp() -> dict[str, object]:
    # The 2018 paper reports only rounded critical coordinates.  A local cusp
    # calculation is therefore used as a reconstruction, not fitted to 89/724.
    result = locate_cusp_critical_point(
        initial_guess=(3.3, 0.44),
        model=PHAModel(OLD_PARAMETERS),
        background_numerics=BackgroundNumerics(
            phiA_tolerance=1.0e-5,
            ricci_tolerance=1.0e-3,
        ),
        phi_step=0.003,
        charge_step=0.0003,
    )
    residual_norm = float(np.linalg.norm(result["cusp_equation_residual"]))
    benchmark_accepted = bool(
        abs(float(result["T_c_MeV"]) - 89.0) < 0.25
        and abs(float(result["mu_B_c_MeV"]) - 724.0) < 1.0
        and residual_norm < 1.0e-3
    )
    if not benchmark_accepted:
        raise RuntimeError(f"The independent 2018 CEP reconstruction failed: {result}")
    result["benchmark_acceptance"] = True
    result["benchmark_residual_norm"] = residual_norm
    result["published_rounded_T_MeV"] = 89.0
    result["published_rounded_mu_B_MeV"] = 724.0
    return result


def _old_rows(cusp: dict[str, object]) -> list[dict[str, object]]:
    _activate(OLD_PARAMETERS)
    landmarks = []
    for name, temperature in (("T400_mu0", 400.0), ("T145_mu0", 145.0)):
        background, uv = _mu0_background(temperature)
        landmarks.append((name, background, uv))
    background = base.integrate_background(
        float(cusp["phi0"]), float(cusp["charge_fraction"])
    )
    landmarks.append(("model_CEP", background, base.extract_uv(background)))
    return [
        _row(
            "Rougemont2018",
            landmark,
            sector,
            background,
            uv,
            PUBLISHED_TAU_FM_C[landmark][sector],
        )
        for landmark, background, uv in landmarks
        for sector in SECTORS
    ]


def main() -> None:
    source_pdf = ROOT / "external" / "rougemont2018" / "paper.pdf"
    local_source_verified = False
    if source_pdf.exists():
        local_hash = hashlib.sha256(source_pdf.read_bytes()).hexdigest()
        if local_hash != SOURCE_PDF_SHA256:
            raise RuntimeError("The local arXiv:1804.00189 PDF checksum changed")
        local_source_verified = True

    pha_rows = _pha_rows()
    cusp = _old_cusp()
    old_rows = _old_rows(cusp)
    rows = old_rows + pha_rows

    maximum_published_difference = max(
        abs(float(row["delta_from_published_fm_c"])) for row in old_rows
    )
    # These numbers were explicitly read from plotted curves ("sim" in the
    # source text), rather than supplied as a numerical table.  A 0.015 fm/c
    # envelope covers the last displayed digit plus plot-reading uncertainty.
    published_tau_acceptance = 0.015
    all_rounded_values_reproduced = all(
        abs(float(row["delta_from_published_fm_c"])) <= published_tau_acceptance
        for row in old_rows
    )
    if not all_rounded_values_reproduced:
        raise RuntimeError("At least one published 2018 time was not reproduced within rounding")
    if max(float(row["shooting_source_residual"]) for row in rows) >= 1.0e-7:
        raise RuntimeError("A benchmark shooting source residual failed")
    if max(float(row["spectral_distance_hat"]) for row in rows) >= 0.02:
        raise RuntimeError("A benchmark spectral/shooting cross-check failed")

    RESULTS.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS / "rougemont2018_benchmark.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "source": {
            "title": "Nonhydrodynamic quasinormal modes and equilibration of a baryon dense holographic QGP with a critical point",
            "citation": "Phys. Rev. D 98 (2018) 034028",
            "doi": "10.1103/PhysRevD.98.034028",
            "url": SOURCE_URL,
            "official_pdf_sha256": SOURCE_PDF_SHA256,
            "local_source_pdf_verified": local_source_verified,
            "published_values_are_plot_readings": True,
        },
        "old_model_parameters": OLD_PARAMETERS.__dict__,
        "independent_old_model_cusp": cusp,
        "row_count": len(rows),
        "old_model_benchmark_count": len(old_rows),
        "all_published_rounded_times_reproduced": all_rounded_values_reproduced,
        "published_tau_acceptance_fm_c": published_tau_acceptance,
        "maximum_abs_tau_difference_from_published_fm_c": maximum_published_difference,
        "maximum_shooting_source_residual": max(
            float(row["shooting_source_residual"]) for row in rows
        ),
        "maximum_spectral_distance_hat": max(
            float(row["spectral_distance_hat"]) for row in rows
        ),
        "table": rows,
    }
    (RESULTS / "rougemont2018_benchmark.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: summary[key] for key in (
        "row_count",
        "old_model_benchmark_count",
        "all_published_rounded_times_reproduced",
        "maximum_abs_tau_difference_from_published_fm_c",
        "maximum_shooting_source_residual",
        "maximum_spectral_distance_hat",
    )}, indent=2))


if __name__ == "__main__":
    main()
