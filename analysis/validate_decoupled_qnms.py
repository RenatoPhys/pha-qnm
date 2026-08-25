#!/usr/bin/env python3
"""Publication-grade checks for the homogeneous decoupled QNM benchmarks."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eig
from scipy.optimize import least_squares

import run_numerics as base


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "python"


def normalizable_power(sector: str) -> float:
    return {"tensor": 4.0, "vector": 2.0,
            "singlet": 2.0 * base.DELTA - 4.0}[sector]


def factored_spectrum(bg: dict[str, object], intervals: int, r_cut: float,
                      sector: str) -> dict[str, object]:
    nodes, D1, _D2 = base.chebyshev_lobatto(intervals, 0.0, r_cut)
    state = base.profile_on_grid(bg, nodes)
    A, h, phi, _Phi, Ap, hp, phip, Phip = state
    emA = np.exp(-A)
    fv, fpv = base.f(phi), base.fp(phi)
    power = normalizable_power(sector)
    App = -phip**2 / 6.0
    if sector == "tensor":
        B0 = 4.0 * h * Ap + hp
        B1 = -2.0j * emA
        C0 = np.zeros_like(h)
        C1 = -3.0j * emA * Ap
    elif sector == "vector":
        ratio = fpv * phip / fv
        B0 = 2.0 * h * Ap + hp + h * ratio
        B1 = -2.0j * emA
        C0 = -emA**2 * fv * Phip**2
        C1 = -1.0j * emA * (Ap + ratio)
    elif sector == "singlet":
        hd = bg["horizon"]
        phipp = np.empty_like(phip)
        phipp[0] = 2.0 * hd["scalar2"]
        phipp[1:] = ((base.Vp(phi[1:]) - 0.5 * emA[1:]**2 * fpv[1:] * Phip[1:]**2
                      - (4.0 * h[1:] * Ap[1:] + hp[1:]) * phip[1:]) / h[1:])
        ratio = phipp / phip - App / Ap
        B0 = 4.0 * h * Ap + hp + 2.0 * h * ratio
        B1 = -2.0j * emA
        C0 = -hp * ratio + emA**2 * Phip**2 / phip * (3.0 * Ap * fpv - fv * phip)
        C1 = -1.0j * emA * (3.0 * Ap + 2.0 * ratio)
    else:
        raise ValueError(sector)
    n = nodes.size
    M0 = np.zeros((n, n), dtype=complex)
    M1 = np.zeros((n, n), dtype=complex)
    for index in range(n):
        M0[index] = (h[index] * _D2[index]
                     + (B0[index] - 2.0 * power * h[index] * Ap[index]) * D1[index])
        M0[index, index] += (h[index] * (power**2 * Ap[index]**2 - power * App[index])
                             - power * Ap[index] * B0[index] + C0[index])
        M1[index] = B1[index] * D1[index]
        M1[index, index] += C1[index] - power * Ap[index] * B1[index]
    # After F=exp(-p A)y, a source-free solution has y'=0 at the exact UV
    # boundary. At finite cutoff this is the leading asymptotic condition.
    M0[-1] = D1[-1]
    M1[-1] = 0.0

    row_norm = np.maximum(np.sum(np.abs(M0), axis=1) + np.sum(np.abs(M1), axis=1), 1e-300)
    M0 = M0 / row_norm[:, np.newaxis]
    M1 = M1 / row_norm[:, np.newaxis]
    omega, right = eig(M0, -M1, left=False, right=True, check_finite=False)
    finite = np.isfinite(omega) & (np.abs(omega) < 80.0) & (omega.imag < 1.0e-7)
    omega, right = omega[finite], right[:, finite]
    residual = []
    for value, rv in zip(omega, right.T):
        matrix = M0 + value * M1
        residual.append(np.linalg.norm(matrix @ rv) /
                        max(np.linalg.norm(matrix, ord=np.inf) * np.linalg.norm(rv), 1e-300))
    order = np.lexsort((np.abs(omega.real), np.abs(omega.imag)))
    return {"nodes": nodes, "omega": omega[order], "right": right[:, order],
            "residual": np.asarray(residual)[order], "M0": M0, "M1": M1}


def eigenvalue_condition(spectrum: dict[str, object], index: int) -> float:
    value = spectrum["omega"][index]
    right = spectrum["right"][:, index]
    matrix = spectrum["M0"] + value * spectrum["M1"]
    left = np.linalg.svd(matrix, full_matrices=False)[0][:, -1]
    denominator = abs(np.vdot(left, spectrum["M1"] @ right))
    return float(math.sqrt(1.0 + abs(value)**2) * np.linalg.norm(left) * np.linalg.norm(right)
                 / max(denominator, 1e-300))


def scalar_ratio(bg: dict[str, object], r: float, state: np.ndarray) -> float:
    A, h, phi, _Phi, Ap, hp, phip, Phip = state
    if r <= 1.01e-7:
        hd = bg["horizon"]
        return 2.0 * hd["scalar2"] / hd["scalar1"] - 2.0 * hd["A2"] / hd["A1"]
    phipp = ((float(base.Vp(phi)) - 0.5 * math.exp(-2.0 * A) * float(base.fp(phi)) * Phip**2
              - (4.0 * h * Ap + hp) * phip) / h)
    return phipp / phip + phip**2 / (6.0 * Ap)


def ode_coefficients(bg: dict[str, object], r: float, sector: str) -> tuple[complex, complex, complex, complex]:
    state = base.background_state(bg, r)
    A, h, phi, _Phi, Ap, hp, phip, Phip = state
    emA = math.exp(-A)
    fv, fpv = float(base.f(phi)), float(base.fp(phi))
    if sector == "tensor":
        return 4.0 * h * Ap + hp, -2.0j * emA, 0.0, -3.0j * emA * Ap
    if sector == "vector":
        ratio = fpv * phip / fv
        return (2.0 * h * Ap + hp + h * ratio, -2.0j * emA,
                -emA**2 * fv * Phip**2, -1.0j * emA * (Ap + ratio))
    if sector == "singlet":
        ratio = scalar_ratio(bg, r, state)
        return (4.0 * h * Ap + hp + 2.0 * h * ratio, -2.0j * emA,
                -hp * ratio + emA**2 * Phip**2 / phip * (3.0 * Ap * fpv - fv * phip),
                -1.0j * emA * (3.0 * Ap + 2.0 * ratio))
    raise ValueError(sector)


def shooting_source(bg: dict[str, object], sector: str, omega: complex,
                    r_uv: float = 10.0, r0: float = 1.0e-5) -> complex:
    B0, B1, C0, C1 = ode_coefficients(bg, r0, sector)
    derivative = -(C0 + omega * C1) / (B0 + omega * B1)
    initial = np.asarray([1.0 + derivative * r0, derivative], dtype=complex)

    def rhs(r: float, value: np.ndarray) -> np.ndarray:
        state = base.background_state(bg, r)
        h = state[1]
        B0r, B1r, C0r, C1r = ode_coefficients(bg, r, sector)
        return np.asarray([value[1],
                           -((B0r + omega * B1r) * value[1]
                             + (C0r + omega * C1r) * value[0]) / h], dtype=complex)

    solution = solve_ivp(rhs, (r0, r_uv), initial, method="DOP853",
                         rtol=3.0e-9, atol=3.0e-11, max_step=0.05,
                         t_eval=np.asarray([r_uv]))
    if not solution.success:
        raise RuntimeError(solution.message)
    state = base.background_state(bg, r_uv)
    return solution.y[0, -1] + solution.y[1, -1] / (normalizable_power(sector) * state[4])


def shooting_root(bg: dict[str, object], sector: str, initial: complex,
                  r_uv: float = 10.0, r0: float = 1.0e-5,
                  acceptance_reference: complex | None = None,
                  acceptance_radius: float = 0.05,
                  search_radius: float = 0.25) -> dict[str, object]:
    def objective(value: np.ndarray) -> np.ndarray:
        source = shooting_source(bg, sector, complex(value[0], value[1]), r_uv=r_uv, r0=r0)
        return np.asarray([source.real, source.imag])

    start = np.asarray([initial.real, initial.imag])
    result = least_squares(objective, start,
                           bounds=(start - search_radius, start + search_radius),
                           xtol=2.0e-7, ftol=2.0e-9, gtol=2.0e-9,
                           max_nfev=36, x_scale="jac")
    value = complex(result.x[0], result.x[1])
    source = shooting_source(bg, sector, value, r_uv=r_uv, r0=r0)
    reference = initial if acceptance_reference is None else acceptance_reference
    accepted = bool(result.success and abs(source) < 1.0e-7
                    and abs(value - reference) < acceptance_radius)
    return {"omega": value, "success": accepted, "optimizer_success": bool(result.success),
            "distance_from_initial": float(abs(value - initial)), "evaluations": int(result.nfev),
            "source_residual": float(abs(source)), "message": str(result.message)}


def nearest_mode(spectrum: dict[str, object], target: complex) -> int:
    return int(np.argmin(np.abs(spectrum["omega"] - target)))


def normalized_overlap(first_nodes: np.ndarray, first_vector: np.ndarray,
                       second_nodes: np.ndarray, second_vector: np.ndarray) -> float:
    common = np.linspace(0.0, min(first_nodes[-1], second_nodes[-1]), 600)
    a = np.interp(common, first_nodes, first_vector.real) + 1.0j * np.interp(
        common, first_nodes, first_vector.imag)
    b = np.interp(common, second_nodes, second_vector.real) + 1.0j * np.interp(
        common, second_nodes, second_vector.imag)
    return float(abs(np.vdot(a, b)) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-300))


def main() -> None:
    bg, _critical = base.load_independent_cusp_background()
    baseline = json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))
    candidates = baseline["pseudospectral_qnm_candidates"]
    rows: list[dict[str, object]] = []
    shooting: list[dict[str, object]] = []
    shooting_convergence: list[dict[str, object]] = []
    spectral_cache: dict[tuple[str, int, float], dict[str, object]] = {}
    for item in candidates:
        sector = str(item["sector"])
        target = complex(float(item["Re_omega_num"]), float(item["Im_omega_num"]))
        print(f"SHOOT {sector} mode={item['mode']} start={target}", flush=True)
        shot = shooting_root(bg, sector, target)
        print(f"SHOT  {sector} mode={item['mode']} root={shot['omega']} "
              f"residual={shot['source_residual']:.3e} success={shot['success']}", flush=True)
        shooting.append({"sector": sector, "mode": int(item["mode"]),
                         "Re_omega": shot["omega"].real, "Im_omega": shot["omega"].imag,
                         "source_residual": shot["source_residual"],
                         "success": shot["success"], "optimizer_success": shot["optimizer_success"],
                         "distance_from_initial": shot["distance_from_initial"],
                         "evaluations": shot["evaluations"]})
        if shot["success"]:
            for study, settings in (
                ("horizon_start", [(value, 10.0) for value in (5.0e-6, 1.0e-5, 2.0e-5)]),
                ("uv_extraction", [(1.0e-5, value) for value in (9.0, 10.0, 11.0)]),
            ):
                for r0_value, ruv_value in settings:
                    variant = shooting_root(bg, sector, shot["omega"], r_uv=ruv_value,
                                            r0=r0_value,
                                            acceptance_reference=shot["omega"])
                    shooting_convergence.append({
                        "sector": sector, "mode": int(item["mode"]), "study": study,
                        "r0": r0_value, "r_uv": ruv_value,
                        "Re_omega": variant["omega"].real,
                        "Im_omega": variant["omega"].imag,
                        "distance_to_reference": abs(variant["omega"] - shot["omega"]),
                        "source_residual": variant["source_residual"],
                        "success": variant["success"],
                    })
        previous = None
        for study, settings in (
            ("radial_resolution", [(n, 10.0) for n in (64, 80, 96, 112, 128, 144, 160)]),
            ("radial_domain", [(160, rc) for rc in (8.0, 9.0, 10.0, 11.0, 12.0)]),
        ):
            for intervals, cutoff in settings:
                key = (sector, intervals, cutoff)
                if key not in spectral_cache:
                    print(f"QZ    {sector} N={intervals} r_cut={cutoff}", flush=True)
                    spectral_cache[key] = factored_spectrum(bg, intervals, cutoff, sector)
                spectrum = spectral_cache[key]
                index = nearest_mode(spectrum, shot["omega"])
                overlap = float("nan")
                if previous is not None:
                    overlap = normalized_overlap(previous["nodes"], previous["right"][:, previous["index"]],
                                                 spectrum["nodes"], spectrum["right"][:, index])
                condition = float("nan")
                if study == "radial_resolution" and intervals == 160 and cutoff == 10.0:
                    condition = eigenvalue_condition(spectrum, index)
                rows.append({"sector": sector, "mode": int(item["mode"]), "study": study,
                             "N": intervals, "r_cut": cutoff,
                             "Re_omega": float(spectrum["omega"][index].real),
                             "Im_omega": float(spectrum["omega"][index].imag),
                             "distance_to_shooting": float(abs(spectrum["omega"][index] - shot["omega"])),
                             "pencil_residual": float(spectrum["residual"][index]),
                             "eigenvalue_condition": condition,
                             "right_overlap_previous": overlap})
                previous = {**spectrum, "index": index}
            previous = None

    with (RESULTS / "qnm_factored_convergence.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (RESULTS / "qnm_shooting_crosscheck.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(shooting[0]))
        writer.writeheader()
        writer.writerows(shooting)
    with (RESULTS / "qnm_shooting_convergence.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(shooting_convergence[0]))
        writer.writeheader()
        writer.writerows(shooting_convergence)

    validated: list[dict[str, object]] = []
    for shot in shooting:
        if not shot["success"]:
            continue
        matching_spectral = [row for row in rows
                             if row["sector"] == shot["sector"] and row["mode"] == shot["mode"]
                             and (row["study"] == "radial_domain" or row["N"] >= 128)]
        matching_shooting = [row for row in shooting_convergence
                             if row["sector"] == shot["sector"] and row["mode"] == shot["mode"]]
        re_error = max(
            max(abs(row["Re_omega"] - shot["Re_omega"]) for row in matching_spectral),
            max(abs(row["Re_omega"] - shot["Re_omega"]) for row in matching_shooting),
        )
        im_error = max(
            max(abs(row["Im_omega"] - shot["Im_omega"]) for row in matching_spectral),
            max(abs(row["Im_omega"] - shot["Im_omega"]) for row in matching_shooting),
        )
        validated.append({
            "sector": shot["sector"], "mode": shot["mode"],
            "Re_omega_num": shot["Re_omega"], "Im_omega_num": shot["Im_omega"],
            "Re_omega_hat": 2.0 * shot["Re_omega"],
            "Im_omega_hat": 2.0 * shot["Im_omega"],
            "Re_error_hat": 2.0 * re_error, "Im_error_hat": 2.0 * im_error,
            "shooting_source_residual": shot["source_residual"],
            "minimum_overlap": min(row["right_overlap_previous"] for row in matching_spectral
                                   if math.isfinite(row["right_overlap_previous"])),
            "maximum_condition": max(row["eigenvalue_condition"] for row in matching_spectral
                                     if math.isfinite(row["eigenvalue_condition"])),
        })
    with (RESULTS / "qnm_validated_modes.csv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=list(validated[0]))
        writer.writeheader()
        writer.writerows(validated)

    summary = {
        "factoring": "F=exp(-p A)y with y'(r_cut)=0",
        "normalizable_powers": {sector: normalizable_power(sector)
                                for sector in ("tensor", "vector", "singlet")},
        "shooting": shooting,
        "maximum_shooting_source_residual": max(row["source_residual"] for row in shooting),
        "maximum_factored_distance_to_shooting": max(
            row["distance_to_shooting"] for row in rows
            if any(item["sector"] == row["sector"] and item["mode"] == row["mode"]
                   and item["success"] for item in shooting)),
        "minimum_overlap": min(row["right_overlap_previous"] for row in rows
                               if math.isfinite(row["right_overlap_previous"])),
        "maximum_eigenvalue_condition": max(row["eigenvalue_condition"] for row in rows
                                            if math.isfinite(row["eigenvalue_condition"])),
        "all_shooting_roots_successful": all(row["success"] for row in shooting),
        "accepted_shooting_roots": sum(row["success"] for row in shooting),
        "maximum_accepted_shooting_parameter_shift": max(
            row["distance_to_reference"] for row in shooting_convergence),
        "all_accepted_shooting_variants_successful": all(
            row["success"] for row in shooting_convergence),
        "validated_modes": validated,
    }
    (RESULTS / "qnm_validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
