#!/usr/bin/env python3
"""Coupled helicity-one/zero collocation pencils from the generated EMD equations.

The primitive fields are kept in radial gauge.  At the UV cutoff the leading
normalizable Robin condition F' + p_F A' F = 0 removes the source coefficient.
The frequency convention is exp(-i omega v + i k z).
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.linalg import eig

import run_numerics as base


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "derivations" / "generated"

SECTOR_SPEC = {
    "helicity1": {
        "fields": ("Hvx", "Hzx", "ax"),
        # E_vx, E_zx and M_x have characteristic second-derivative
        # coefficients proportional to h and therefore impose ingoing
        # regularity directly at the EF horizon.  E_rx is retained as the
        # radial constraint, not used as a bulk evolution equation.
        "equations": ("E_vx", "evolution_E_zx", "evolution_M_x"),
        "powers": (4.0, 4.0, 2.0),
        "diagnostics": ("evolution_E_rx", "constraint_Er_x"),
    },
    "helicity0": {
        "fields": ("Hvv", "Hvz", "Hzz", "Haa", "av", "az", "varphi"),
        # The four Einstein equations with h-weighted radial principal part
        # are characteristic at the future horizon.  M_v is a radial Gauss
        # equation; its horizon row is replaced by M_r below.
        "equations": ("diagnostic_E_vv", "diagnostic_E_vz", "diagnostic_E_zz",
                      "evolution_E_aa", "evolution_M_v", "evolution_M_z",
                      "evolution_scalar"),
        "powers": (4.0, 4.0, 4.0, 4.0, 2.0, 2.0, base.DELTA),
        "diagnostics": ("evolution_E_rr", "evolution_E_rv", "evolution_E_rz",
                        "constraint_Er_r", "constraint_Er_v", "constraint_Er_z",
                        "constraint_M_r"),
        "horizon_constraints": {4: "constraint_M_r"},
        # Four UV Dirichlet rows select a representative of the four
        # residual radial-gauge orbits.  The remaining three Robin rows are
        # precisely the scalar, charge and sound gauge-invariant sources in
        # this gauge.
        "uv_gauge_fields": (0, 2, 3, 5),
    },
}


@lru_cache(maxsize=2)
def generated_data(sector: str) -> dict[str, object]:
    return json.loads((GENERATED / f"{sector}_equations.json").read_text(encoding="utf-8"))


def background_coefficients(bg: dict[str, object], nodes: np.ndarray) -> dict[str, np.ndarray]:
    state = base.profile_on_grid(bg, nodes)
    A, h, phi, Phi, Ar, hr, phir, Phir = state
    fv, fphi, fphiphi = base.f(phi), base.fp(phi), base.fpp(phi)
    Arr = -phir**2 / 6.0
    hrr = -4.0 * Ar * hr + np.exp(-2.0 * A) * fv * Phir**2
    phirr = np.empty_like(phir)
    Phirr = np.empty_like(Phir)
    hd = bg["horizon"]
    phirr[0] = 2.0 * hd["scalar2"]
    Phirr[0] = 2.0 * hd["Phi2"]
    positive = nodes > 0.0
    phirr[positive] = ((base.Vp(phi[positive])
                        - 0.5 * np.exp(-2.0 * A[positive]) * fphi[positive]
                          * Phir[positive]**2
                        - (4.0 * h[positive] * Ar[positive] + hr[positive])
                          * phir[positive]) / h[positive])
    Phirr[positive] = Phir[positive] * (-2.0 * Ar[positive]
                                        - fphi[positive] * phir[positive] / fv[positive])
    return {
        "A": A, "Ar": Ar, "Arr": Arr,
        "h": h, "hr": hr, "hrr": hrr,
        "phi": phi, "phir": phir, "phirr": phirr,
        "Phi": Phi, "Phir": Phir, "Phirr": Phirr,
        "f": fv, "fr": fphi * phir, "fphi": fphi,
        "fphir": fphiphi * phir, "fphiphi": fphiphi,
        "V": base.V(phi), "Vphi": base.Vp(phi), "Vphiphi": base.Vpp(phi),
    }


def coefficient_array(expression: str, env: dict[str, object]) -> np.ndarray:
    # Expressions are generated locally by SymPy and contain only arithmetic
    # plus exp.  Disabling builtins keeps evaluation deliberately narrow.
    value = eval(expression, {"__builtins__": {}}, {**env, "exp": np.exp})
    return np.asarray(value, dtype=complex)


def local_coefficients(sector: str, bg_values: dict[str, np.ndarray],
                       sigma: complex, momentum: complex,
                       equation_names: tuple[str, ...] | None = None) -> np.ndarray:
    data = generated_data(sector)
    fields = tuple(data["fields"])
    equations = data["equations"]
    names = equation_names or tuple(equations)
    n = len(next(iter(bg_values.values())))
    out = np.zeros((len(names), len(fields), 3, n), dtype=complex)
    env: dict[str, object] = {**bg_values, "sigma": sigma, "p": momentum}
    for ieq, name in enumerate(names):
        for coefficient, expression in equations[name]["coefficients"].items():
            for field_index, field in enumerate(fields):
                if coefficient == field:
                    derivative = 0
                elif coefficient == field + "r":
                    derivative = 1
                elif coefficient == field + "rr":
                    derivative = 2
                else:
                    continue
                value = coefficient_array(expression, env)
                out[ieq, field_index, derivative] = value
                break
            else:
                raise RuntimeError(f"Cannot map generated coefficient {coefficient}")
    return out


def assemble_pencil(bg: dict[str, object], sector: str, qhat: float,
                    intervals: int = 64, r_cut: float = 10.0) -> dict[str, object]:
    spec = SECTOR_SPEC[sector]
    fields = tuple(spec["fields"])
    equations = tuple(spec["equations"])
    nodes, D1, D2 = base.chebyshev_lobatto(intervals, 0.0, r_cut)
    bg_values = background_coefficients(bg, nodes)
    h0 = base.extract_uv(bg)["h0"]
    k_numeric = qhat / (2.0 * math.sqrt(h0))
    momentum = 1.0j * k_numeric
    c0 = local_coefficients(sector, bg_values, 0.0j, momentum, equations)
    c1 = local_coefficients(sector, bg_values, 1.0 + 0.0j, momentum, equations) - c0
    field_count, n = len(fields), nodes.size
    M0 = np.zeros((field_count * n, field_count * n), dtype=complex)
    M1 = np.zeros_like(M0)
    for ieq in range(field_count):
        rows = slice(ieq * n, (ieq + 1) * n)
        for jfield in range(field_count):
            cols = slice(jfield * n, (jfield + 1) * n)
            M0[rows, cols] = (c0[ieq, jfield, 2, :, None] * D2
                              + c0[ieq, jfield, 1, :, None] * D1
                              + np.diag(c0[ieq, jfield, 0]))
            M1[rows, cols] = (c1[ieq, jfield, 2, :, None] * D2
                              + c1[ieq, jfield, 1, :, None] * D1
                              + np.diag(c1[ieq, jfield, 0]))
    # Non-characteristic radial equations require one independent regularity
    # relation at the future horizon.  The mapping names the generated
    # constraint that replaces that equation's first collocation row.
    for equation_index, constraint_name in spec.get("horizon_constraints", {}).items():
        d0 = local_coefficients(sector, bg_values, 0.0j, momentum,
                                (constraint_name,))[0]
        d1 = (local_coefficients(sector, bg_values, 1.0 + 0.0j, momentum,
                                 (constraint_name,))[0] - d0)
        row = equation_index * n
        M0[row] = 0.0
        M1[row] = 0.0
        for jfield in range(field_count):
            cols = slice(jfield * n, (jfield + 1) * n)
            M0[row, cols] = (d0[jfield, 2, 0] * D2[0]
                             + d0[jfield, 1, 0] * D1[0])
            M0[row, jfield * n] += d0[jfield, 0, 0]
            M1[row, cols] = (d1[jfield, 2, 0] * D2[0]
                             + d1[jfield, 1, 0] * D1[0])
            M1[row, jfield * n] += d1[jfield, 0, 0]
    # Replace one UV equation per primitive field by the source-free Robin row.
    # Metric, Maxwell and scalar normalizable powers are respectively 4, 2,
    # and Delta in u=exp(-A).
    uv_gauge_fields = set(spec.get("uv_gauge_fields", ()))
    for field_index, power in enumerate(spec["powers"]):
        row = field_index * n + (n - 1)
        M0[row] = 0.0
        M1[row] = 0.0
        if field_index in uv_gauge_fields:
            M0[row, field_index * n + n - 1] = 1.0
        else:
            cols = slice(field_index * n, (field_index + 1) * n)
            M0[row, cols] = D1[-1]
            M0[row, field_index * n + n - 1] += power * bg_values["Ar"][-1]
    row_norm = np.maximum(np.sum(np.abs(M0), axis=1) + np.sum(np.abs(M1), axis=1), 1.0e-300)
    M0 /= row_norm[:, None]
    M1 /= row_norm[:, None]
    return {"M0": M0, "M1": M1, "nodes": nodes, "background": bg_values,
            "k_numeric": k_numeric, "qhat": qhat, "sector": sector}


def spectrum(bg: dict[str, object], sector: str, qhat: float,
             intervals: int = 64, r_cut: float = 10.0,
             max_abs_omega: float = 20.0,
             imaginary_ceiling: float = 1.0e-7) -> dict[str, object]:
    pencil = assemble_pencil(bg, sector, qhat, intervals, r_cut)
    # L(sigma)=M0+sigma M1 and sigma=-i omega, hence omega=i sigma.
    sigma, right = eig(pencil["M0"], -pencil["M1"], right=True, check_finite=False)
    omega = 1.0j * sigma
    finite = (np.isfinite(omega) & (np.abs(omega) < max_abs_omega)
              & (omega.imag < imaginary_ceiling))
    omega, right = omega[finite], right[:, finite]
    residual = np.empty(omega.size)
    for index, (value, vector) in enumerate(zip(omega, right.T)):
        sigma_value = -1.0j * value
        matrix = pencil["M0"] + sigma_value * pencil["M1"]
        residual[index] = np.linalg.norm(matrix @ vector) / max(
            np.linalg.norm(matrix, ord=np.inf) * np.linalg.norm(vector), 1.0e-300)
    order = np.lexsort((np.abs(omega.real), np.abs(omega.imag)))
    return {**pencil, "omega": omega[order], "right": right[:, order],
            "residual": residual[order]}


def field_content(result: dict[str, object], mode_index: int) -> dict[str, float]:
    fields = SECTOR_SPEC[str(result["sector"])]["fields"]
    n = len(result["nodes"])
    vector = result["right"][:, mode_index]
    norms = np.asarray([np.linalg.norm(vector[i * n:(i + 1) * n]) for i in range(len(fields))])
    norms /= max(np.linalg.norm(norms), 1.0e-300)
    return {field: float(value) for field, value in zip(fields, norms)}


def diagnostic_residual(result: dict[str, object], mode_index: int) -> dict[str, float]:
    """Residuals of generated equations not used to construct the pencil."""
    sector = str(result["sector"])
    names = tuple(SECTOR_SPEC[sector]["diagnostics"])
    nodes = result["nodes"]
    n = len(nodes)
    vector = result["right"][:, mode_index]
    omega = result["omega"][mode_index]
    sigma = -1.0j * omega
    coefficients = local_coefficients(sector, result["background"], sigma,
                                      1.0j * result["k_numeric"], names)
    _nodes, D1, D2 = base.chebyshev_lobatto(n - 1, float(nodes[0]), float(nodes[-1]))
    values: dict[str, float] = {}
    for ieq, name in enumerate(names):
        operator = np.zeros((n, len(vector)), dtype=complex)
        for jfield in range(len(SECTOR_SPEC[sector]["fields"])):
            cols = slice(jfield * n, (jfield + 1) * n)
            operator[:, cols] = (coefficients[ieq, jfield, 2, :, None] * D2
                                 + coefficients[ieq, jfield, 1, :, None] * D1
                                 + np.diag(coefficients[ieq, jfield, 0]))
        interior = operator[1:-1]
        row_norm = np.maximum(np.sum(np.abs(interior), axis=1), 1.0e-300)
        scaled = interior / row_norm[:, None]
        values[name] = float(np.linalg.norm(scaled @ vector)
                             / max(np.linalg.norm(vector), 1.0e-300))
    return values
