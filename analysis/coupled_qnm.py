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
from scipy.integrate import cumulative_trapezoid
from scipy.integrate import solve_ivp
from scipy.linalg import eig, null_space
from scipy.optimize import minimize

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
        # E_rx is the independent radial Einstein constraint.  The raised
        # combination E^r_x=e^-A E_vx+h E_rx is algebraically equivalent but
        # loses precision in the UV through cancellation of large terms.
        "diagnostics": ("evolution_E_rx",),
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
        # Retain the primitive radial Einstein/Maxwell equations as the
        # independent diagnostics.  The raised combinations are deliberately
        # excluded here because their expanded UV coefficients suffer
        # catastrophic cancellation even when the primitive equations vanish.
        "diagnostics": ("evolution_E_rr", "evolution_E_rv", "evolution_E_rz",
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
    value = eval(_compiled_expression(expression), {"__builtins__": {}},
                 {**env, "exp": np.exp})
    return np.asarray(value, dtype=complex)


@lru_cache(maxsize=None)
def _compiled_expression(expression: str):
    return compile(expression, "<generated-emd-coefficient>", "eval")


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


def radial_grid(intervals: int, r_cut: float,
                map_alpha: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Chebyshev grid with an optional exponential UV-clustering map.

    ``map_alpha < 0`` places additional resolution near the AdS boundary while
    keeping a finite nonzero Jacobian at both endpoints.  The second derivative
    includes the complete chain-rule term, so the mapped grid is an independent
    radial-domain convergence control rather than a change of equations.
    """
    if abs(map_alpha) < 1.0e-14:
        return base.chebyshev_lobatto(intervals, 0.0, r_cut)
    t, Dt, D2t = base.chebyshev_lobatto(intervals, 0.0, 1.0)
    denominator = np.expm1(map_alpha)
    exponential = np.exp(map_alpha * t)
    nodes = r_cut * np.expm1(map_alpha * t) / denominator
    dr_dt = r_cut * map_alpha * exponential / denominator
    d2r_dt2 = r_cut * map_alpha**2 * exponential / denominator
    D1 = Dt / dr_dt[:, None]
    D2 = D2t / dr_dt[:, None]**2 - (d2r_dt2 / dr_dt**3)[:, None] * Dt
    return nodes, D1, D2


def assemble_pencil(bg: dict[str, object], sector: str, qhat: float,
                    intervals: int = 64, r_cut: float = 10.0,
                    map_alpha: float = 0.0) -> dict[str, object]:
    spec = SECTOR_SPEC[sector]
    fields = tuple(spec["fields"])
    equations = tuple(spec["equations"])
    nodes, D1, D2 = radial_grid(intervals, r_cut, map_alpha)
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
    return {"M0": M0, "M1": M1, "nodes": nodes, "D1": D1, "D2": D2,
            "background": bg_values, "map_alpha": map_alpha,
            "k_numeric": k_numeric, "qhat": qhat, "sector": sector}


def spectrum(bg: dict[str, object], sector: str, qhat: float,
             intervals: int = 64, r_cut: float = 10.0,
             map_alpha: float = 0.0,
             max_abs_omega: float = 20.0,
             imaginary_ceiling: float = 1.0e-7) -> dict[str, object]:
    pencil = assemble_pencil(bg, sector, qhat, intervals, r_cut, map_alpha)
    # L(sigma)=M0+sigma M1 and sigma=-i omega, hence omega=i sigma.
    sigma, left, right = eig(pencil["M0"], -pencil["M1"], left=True,
                             right=True, check_finite=False)
    with np.errstate(invalid="ignore"):
        omega = 1.0j * sigma
    finite = (np.isfinite(omega) & (np.abs(omega) < max_abs_omega)
              & (omega.imag < imaginary_ceiling))
    omega, sigma = omega[finite], sigma[finite]
    left, right = left[:, finite], right[:, finite]
    residual = np.empty(omega.size)
    for index, (value, vector) in enumerate(zip(omega, right.T)):
        sigma_value = -1.0j * value
        matrix = pencil["M0"] + sigma_value * pencil["M1"]
        residual[index] = np.linalg.norm(matrix @ vector) / max(
            np.linalg.norm(matrix, ord=np.inf) * np.linalg.norm(vector), 1.0e-300)
    condition = np.empty(omega.size)
    for index, (value, lv, rv) in enumerate(zip(sigma, left.T, right.T)):
        denominator = abs(np.vdot(lv, pencil["M1"] @ rv))
        condition[index] = (math.sqrt(1.0 + abs(value)**2)
                            * np.linalg.norm(lv) * np.linalg.norm(rv)
                            / max(denominator, 1.0e-300))
    order = np.lexsort((np.abs(omega.real), np.abs(omega.imag)))
    return {**pencil, "omega": omega[order], "sigma": sigma[order],
            "left": left[:, order], "right": right[:, order],
            "condition": condition[order],
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
    D1, D2 = result["D1"], result["D2"]
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


def _local_operator(coefficients: np.ndarray, D1: np.ndarray,
                    D2: np.ndarray) -> np.ndarray:
    """Assemble a radial operator from local field coefficients."""
    equation_count, field_count, _derivatives, n = coefficients.shape
    operator = np.zeros((equation_count * n, field_count * n), dtype=complex)
    for ieq in range(equation_count):
        rows = slice(ieq * n, (ieq + 1) * n)
        for jfield in range(field_count):
            cols = slice(jfield * n, (jfield + 1) * n)
            operator[rows, cols] = (
                coefficients[ieq, jfield, 2, :, None] * D2
                + coefficients[ieq, jfield, 1, :, None] * D1
                + np.diag(coefficients[ieq, jfield, 0])
            )
    return operator


def _invariant_maps(sector: str, background: dict[str, np.ndarray],
                    sigma: complex, momentum: complex) -> tuple[list[np.ndarray], list[float]]:
    """Pointwise maps from primitive fields to physical gauge invariants."""
    n = len(background["A"])
    field_count = len(SECTOR_SPEC[sector]["fields"])

    def empty_map() -> np.ndarray:
        return np.zeros((n, field_count * n), dtype=complex)

    diagonal = np.arange(n)
    if sector == "helicity1":
        shear = empty_map()
        shear[diagonal, diagonal] = momentum
        shear[diagonal, n + diagonal] = -sigma
        current = empty_map()
        current[diagonal, 2 * n + diagonal] = 1.0
        return [shear, current], [4.0, 2.0]

    Ar = background["Ar"]
    safe_Ar = np.where(np.abs(Ar) > 1.0e-14, Ar, 1.0e-14)
    sound = empty_map()
    sound[diagonal, diagonal] = momentum**2
    sound[diagonal, n + diagonal] = -2.0 * momentum * sigma
    sound[diagonal, 2 * n + diagonal] = sigma**2
    sound[diagonal, 3 * n + diagonal] = (
        momentum**2 * (background["h"] + background["hr"] / (2.0 * safe_Ar))
        - sigma**2
    )
    charge = empty_map()
    charge[diagonal, 3 * n + diagonal] = (
        -momentum * background["Phir"] / (2.0 * safe_Ar)
    )
    charge[diagonal, 4 * n + diagonal] = momentum
    charge[diagonal, 5 * n + diagonal] = -sigma
    scalar = empty_map()
    scalar[diagonal, 3 * n + diagonal] = -background["phir"] / (2.0 * safe_Ar)
    scalar[diagonal, 6 * n + diagonal] = 1.0
    return [sound, charge, scalar], [4.0, 2.0, base.DELTA]


def _residual_gauge_vectors(sector: str, nodes: np.ndarray,
                            background: dict[str, np.ndarray], sigma: complex,
                            momentum: complex) -> np.ndarray:
    """Sample the exact residual radial-gauge orbits on the collocation grid."""
    n = len(nodes)
    if sector == "helicity1":
        vector = np.zeros(3 * n, dtype=complex)
        vector[:n] = -sigma
        vector[n:2 * n] = -momentum
        return vector[:, None]

    A = background["A"]
    exp_minus_A = np.exp(-A)
    integral_plus = cumulative_trapezoid(np.exp(A), nodes, initial=0.0)
    integral_minus = cumulative_trapezoid(exp_minus_A, nodes, initial=0.0)
    bases = (
        (exp_minus_A, np.zeros(n), np.zeros(n), np.zeros(n)),
        (-sigma * exp_minus_A * integral_plus, np.ones(n),
         -momentum * integral_minus, np.zeros(n)),
        (np.zeros(n), np.zeros(n), np.ones(n), np.zeros(n)),
        (np.zeros(n), np.zeros(n), np.zeros(n), np.ones(n)),
    )
    vectors = []
    for R, T, Z, gauge_lambda in bases:
        fields = (
            R * (2.0 * background["Ar"] * background["h"] + background["hr"])
            + 2.0 * background["h"] * sigma * T
            - 2.0 * exp_minus_A * sigma * R,
            -sigma * Z + background["h"] * momentum * T
            - exp_minus_A * momentum * R,
            -2.0 * background["Ar"] * R - 2.0 * momentum * Z,
            -2.0 * background["Ar"] * R,
            -R * background["Phir"] - background["Phi"] * sigma * T
            - sigma * gauge_lambda,
            -background["Phi"] * momentum * T - momentum * gauge_lambda,
            -R * background["phir"],
        )
        vectors.append(np.concatenate(fields))
    return np.column_stack(vectors)


def _residual_gauge_horizon_data(bg: dict[str, object], sector: str,
                                 sigma: complex,
                                 momentum: complex) -> np.ndarray:
    """Exact value/first-derivative data for residual gauges at the horizon."""
    if sector == "helicity1":
        values = np.asarray([[-sigma], [-momentum], [0.0]], dtype=complex)
        return np.vstack((values, np.zeros_like(values)))

    hd = bg["horizon"]
    Ar = float(hd["A1"])
    Arr = 2.0 * float(hd["A2"])
    hrr = 2.0 * float(hd["h2"])
    phir = float(hd["scalar1"])
    phirr = 2.0 * float(hd["scalar2"])
    Phir = float(bg["Phi1"])
    Phirr = 2.0 * float(hd["Phi2"])

    R = np.asarray([1.0, 0.0, 0.0, 0.0], dtype=complex)
    T = np.asarray([0.0, 1.0, 0.0, 0.0], dtype=complex)
    Z = np.asarray([0.0, 0.0, 1.0, 0.0], dtype=complex)
    gauge_lambda = np.asarray([0.0, 0.0, 0.0, 1.0], dtype=complex)
    Rr = -Ar * R - sigma * T
    Zr = -momentum * T

    values = np.vstack((
        R - 2.0 * sigma * R,
        -sigma * Z - momentum * R,
        -2.0 * Ar * R - 2.0 * momentum * Z,
        -2.0 * Ar * R,
        -R * Phir - sigma * gauge_lambda,
        -momentum * gauge_lambda,
        -R * phir,
    ))
    derivatives = np.vstack((
        Rr + R * (2.0 * Ar + hrr) + 2.0 * sigma * T
        - 2.0 * sigma * (Rr - Ar * R),
        -sigma * Zr + momentum * T - momentum * (Rr - Ar * R),
        -2.0 * Arr * R - 2.0 * Ar * Rr - 2.0 * momentum * Zr,
        -2.0 * Arr * R - 2.0 * Ar * Rr,
        -Rr * Phir - R * Phirr - Phir * sigma * T,
        -Phir * momentum * T,
        -Rr * phir - R * phirr,
    ))
    return np.vstack((values, derivatives))


def _horizon_second_coefficient(bg: dict[str, object], sector: str,
                                sigma: complex, momentum: complex,
                                values: np.ndarray,
                                derivatives: np.ndarray,
                                derivative_step: float = 2.0e-5
                                ) -> tuple[np.ndarray, float]:
    """Solve the next regular horizon recurrence coefficient.

    For characteristic rows the leading equation fixes ``y1`` and its radial
    derivative fixes ``y2``.  In helicity zero the noncharacteristic ``M_v``
    row supplies the seventh equation for ``y2``.  This prevents an O(r0)
    truncation from exciting the constraint-violating DAE solution.
    """
    equations = tuple(SECTOR_SPEC[sector]["equations"])

    def coefficients_at(radius: float) -> np.ndarray:
        background = background_coefficients(bg, np.asarray([radius]))
        return local_coefficients(
            sector, background, sigma, momentum, equations
        )[:, :, :, 0]

    at_zero = coefficients_at(0.0)
    at_step = coefficients_at(derivative_step)
    at_twice = coefficients_at(2.0 * derivative_step)
    radial_derivative = (-3.0 * at_zero + 4.0 * at_step - at_twice) / (
        2.0 * derivative_step
    )
    if sector == "helicity1":
        characteristic = list(range(len(equations)))
        recurrence = 2.0 * (
            radial_derivative[characteristic, :, 2]
            + at_zero[characteristic, :, 1]
        )
        right = -(
            radial_derivative[characteristic, :, 1] @ derivatives
            + at_zero[characteristic, :, 0] @ derivatives
            + radial_derivative[characteristic, :, 0] @ values
        )
    else:
        characteristic = [index for index in range(len(equations)) if index != 4]
        recurrence_characteristic = 2.0 * (
            radial_derivative[characteristic, :, 2]
            + at_zero[characteristic, :, 1]
        )
        right_characteristic = -(
            radial_derivative[characteristic, :, 1] @ derivatives
            + at_zero[characteristic, :, 0] @ derivatives
            + radial_derivative[characteristic, :, 0] @ values
        )
        recurrence_maxwell = 2.0 * at_zero[4, :, 2][None, :]
        right_maxwell = -(
            at_zero[4, :, 1] @ derivatives + at_zero[4, :, 0] @ values
        )[None, :]
        recurrence = np.vstack((recurrence_characteristic, recurrence_maxwell))
        right = np.vstack((right_characteristic, right_maxwell))
    condition = float(np.linalg.cond(recurrence))
    second = np.linalg.solve(recurrence, right)
    return second, condition


def shooting_source_matrix(bg: dict[str, object], sector: str, qhat: float,
                           omega: complex, r_uv: float = 6.0,
                           r0: float = 1.0e-5, rtol: float = 2.0e-9,
                           atol: float = 2.0e-11,
                           max_step: float = 0.025,
                           horizon_einstein_constraint: tuple[int, str] | None = None
                           ) -> dict[str, object]:
    """Integrate an independent infalling basis and extract its UV sources.

    The regular horizon data are obtained as the null space of the exact
    characteristic equations.  Analytic residual-gauge data are quotiented in
    that finite-dimensional space, leaving two helicity-one or three
    helicity-zero physical solutions.  No collocation or generalized
    eigenvalue matrix enters this calculation.
    """
    if not (0.0 < r0 < r_uv < float(bg["r_max"])):
        raise ValueError("shooting radii must satisfy 0 < r0 < r_uv < r_max")
    spec = SECTOR_SPEC[sector]
    equations = tuple(spec["equations"])
    field_count = len(spec["fields"])
    horizon_background = background_coefficients(bg, np.asarray([0.0]))
    h0 = base.extract_uv(bg)["h0"]
    k_numeric = qhat / (2.0 * math.sqrt(h0))
    momentum = 1.0j * k_numeric
    sigma = -1.0j * omega
    coefficients = local_coefficients(
        sector, horizon_background, sigma, momentum, equations
    )[:, :, :, 0]

    if sector == "helicity1":
        horizon_zero = coefficients[:, :, 0]
        horizon_first = coefficients[:, :, 1]
    else:
        # M_v contains the noncharacteristic a_v'' row.  Its Gauss-law
        # counterpart M_r supplies the missing regular horizon relation.
        characteristic = [index for index in range(len(equations)) if index != 4]
        maxwell_constraint = local_coefficients(
            sector, horizon_background, sigma, momentum, ("constraint_M_r",)
        )[0, :, :, 0]
        extra_constraints = [maxwell_constraint]
        if horizon_einstein_constraint is not None:
            dropped, constraint_name = horizon_einstein_constraint
            if dropped not in characteristic or dropped > 3:
                raise ValueError("the dropped Einstein horizon row must be in 0..3")
            characteristic.remove(dropped)
            extra_constraints.append(local_coefficients(
                sector, horizon_background, sigma, momentum, (constraint_name,)
            )[0, :, :, 0])
        horizon_zero = np.vstack((coefficients[characteristic, :, 0],
                                  *(item[:, 0] for item in extra_constraints)))
        horizon_first = np.vstack((coefficients[characteristic, :, 1],
                                   *(item[:, 1] for item in extra_constraints)))
    horizon_operator = np.hstack((horizon_zero, horizon_first))
    regular_basis = null_space(horizon_operator)
    gauge_data = _residual_gauge_horizon_data(bg, sector, sigma, momentum)
    gauge_defect = float(np.linalg.norm(horizon_operator @ gauge_data)
                         / max(np.linalg.norm(horizon_operator)
                               * np.linalg.norm(gauge_data), 1.0e-300))
    if gauge_defect > 2.0e-9:
        raise RuntimeError(f"analytic horizon gauge defect {gauge_defect:.3e}")
    gauge_coordinates = regular_basis.conj().T @ gauge_data
    physical_coordinates = null_space(gauge_coordinates.conj().T)
    initial_data = regular_basis @ physical_coordinates
    physical_count = physical_coordinates.shape[1]
    expected_count = 2 if sector == "helicity1" else 3
    if physical_count != expected_count:
        raise RuntimeError(
            f"expected {expected_count} physical horizon solutions, got {physical_count}"
        )
    values0 = initial_data[:field_count]
    derivatives0 = initial_data[field_count:]
    second0, horizon_recurrence_condition = _horizon_second_coefficient(
        bg, sector, sigma, momentum, values0, derivatives0
    )
    initial = np.vstack((
        values0 + r0 * derivatives0 + r0**2 * second0,
        derivatives0 + 2.0 * r0 * second0,
    ))

    maximum_c2_condition = 0.0

    def rhs(radius: float, flat: np.ndarray) -> np.ndarray:
        nonlocal maximum_c2_condition
        state = flat.reshape(2 * field_count, physical_count)
        values, derivatives = state[:field_count], state[field_count:]
        background = background_coefficients(bg, np.asarray([radius]))
        local = local_coefficients(
            sector, background, sigma, momentum, equations
        )[:, :, :, 0]
        C0, C1, C2 = local[:, :, 0], local[:, :, 1], local[:, :, 2]
        row_scale = np.maximum(
            np.max(np.abs(C0), axis=1) + np.max(np.abs(C1), axis=1)
            + np.max(np.abs(C2), axis=1), 1.0e-300
        )
        C0, C1, C2 = (item / row_scale[:, None] for item in (C0, C1, C2))
        condition = float(np.linalg.cond(C2))
        maximum_c2_condition = max(maximum_c2_condition, condition)
        second = np.linalg.solve(C2, -(C1 @ derivatives + C0 @ values))
        return np.vstack((derivatives, second)).reshape(-1)

    integration = solve_ivp(
        rhs, (r0, r_uv), initial.reshape(-1), method="DOP853",
        rtol=rtol, atol=atol, max_step=max_step, dense_output=True
    )
    if not integration.success:
        raise RuntimeError(integration.message)

    derivative_step = min(2.0e-4, 0.2 * (r_uv - r0))

    def invariant_values(radius: float) -> tuple[np.ndarray, dict[str, np.ndarray]]:
        state = integration.sol(radius).reshape(2 * field_count, physical_count)
        background = background_coefficients(bg, np.asarray([radius]))
        maps, _powers = _invariant_maps(sector, background, sigma, momentum)
        values = np.vstack([invariant_map @ state[:field_count]
                            for invariant_map in maps])
        return values, background

    invariant, uv_background = invariant_values(r_uv)
    invariant_plus, _ = invariant_values(r_uv + derivative_step)
    invariant_minus, _ = invariant_values(r_uv - derivative_step)
    invariant_derivative = (invariant_plus - invariant_minus) / (2.0 * derivative_step)
    _maps, powers = _invariant_maps(sector, uv_background, sigma, momentum)
    Ar_uv = float(uv_background["Ar"][0])
    source_rows = []
    for index, power in enumerate(powers):
        value = invariant[index]
        derivative = invariant_derivative[index]
        if sector == "helicity0" and index == len(powers) - 1:
            u_uv = math.exp(-float(uv_background["A"][0]))
            denominator = ((base.DELTA - base.NU) * Ar_uv * u_uv**base.NU)
            source_rows.append((derivative + base.DELTA * Ar_uv * value) / denominator)
        else:
            source_rows.append(value + derivative / (power * Ar_uv))
    source_matrix = np.vstack(source_rows)
    # Use analytic kinematic scales, never the instantaneous row norm: the
    # latter would divide out an entire vanishing source row in a decoupled
    # channel and erase the QNM zero itself.
    if sector == "helicity1":
        source_scales = np.asarray([
            max(abs(momentum), abs(sigma), 1.0e-12), 1.0
        ])
    else:
        source_scales = np.asarray([
            max(abs(momentum)**2, abs(momentum * sigma), abs(sigma)**2, 1.0e-12),
            max(abs(momentum), abs(sigma), 1.0e-12),
            1.0,
        ])
    normalized_source = source_matrix / source_scales[:, None]
    _left, singular_values, right_h = np.linalg.svd(normalized_source)
    determinant = complex(np.linalg.det(normalized_source))
    return {
        "omega": omega, "sigma": sigma, "qhat": qhat,
        "k_numeric": k_numeric, "sector": sector,
        "source_matrix": source_matrix,
        "normalized_source_matrix": normalized_source,
        "source_singular_values": singular_values,
        "source_singular_value": float(singular_values[-1]),
        "normalized_determinant": determinant,
        "source_vector": right_h[-1].conj(),
        "horizon_gauge_defect": gauge_defect,
        "horizon_recurrence_condition": horizon_recurrence_condition,
        "maximum_c2_condition": maximum_c2_condition,
        "integration_steps": int(integration.t.size),
        "r0": r0, "r_uv": r_uv,
        "horizon_einstein_constraint": horizon_einstein_constraint,
    }


def collocation_source_matrix(bg: dict[str, object], sector: str, qhat: float,
                              omega: complex, intervals: int = 64,
                              r_uv: float = 6.0,
                              r0: float = 1.0e-4,
                              map_alpha: float = 0.0,
                              subdomains: tuple[float, ...] | None = None,
                              horizon_einstein_constraint: tuple[int, str] | None = None
                              ) -> dict[str, object]:
    """Chebyshev integration of the physical horizon source basis.

    This is an independent discretization of the first-order radial initial
    value system.  It shares only the generated differential kernels and the
    analytic horizon quotient with :func:`shooting_source_matrix`; it does not
    call an ODE integrator or construct the primitive generalized pencil.
    """
    if not (0.0 < r0 < r_uv < float(bg["r_max"])):
        raise ValueError("collocation radii must satisfy 0 < r0 < r_uv < r_max")
    spec = SECTOR_SPEC[sector]
    equations = tuple(spec["equations"])
    field_count = len(spec["fields"])
    h0 = base.extract_uv(bg)["h0"]
    k_numeric = qhat / (2.0 * math.sqrt(h0))
    momentum = 1.0j * k_numeric
    sigma = -1.0j * omega

    horizon_background = background_coefficients(bg, np.asarray([0.0]))
    horizon_coefficients = local_coefficients(
        sector, horizon_background, sigma, momentum, equations
    )[:, :, :, 0]
    if sector == "helicity1":
        horizon_zero = horizon_coefficients[:, :, 0]
        horizon_first = horizon_coefficients[:, :, 1]
    else:
        characteristic = [index for index in range(len(equations)) if index != 4]
        maxwell_constraint = local_coefficients(
            sector, horizon_background, sigma, momentum, ("constraint_M_r",)
        )[0, :, :, 0]
        extra_constraints = [maxwell_constraint]
        if horizon_einstein_constraint is not None:
            dropped, constraint_name = horizon_einstein_constraint
            if dropped not in characteristic or dropped > 3:
                raise ValueError("the dropped Einstein horizon row must be in 0..3")
            characteristic.remove(dropped)
            extra_constraints.append(local_coefficients(
                sector, horizon_background, sigma, momentum, (constraint_name,)
            )[0, :, :, 0])
        horizon_zero = np.vstack((horizon_coefficients[characteristic, :, 0],
                                  *(item[:, 0] for item in extra_constraints)))
        horizon_first = np.vstack((horizon_coefficients[characteristic, :, 1],
                                   *(item[:, 1] for item in extra_constraints)))
    horizon_operator = np.hstack((horizon_zero, horizon_first))
    regular_basis = null_space(horizon_operator)
    gauge_data = _residual_gauge_horizon_data(bg, sector, sigma, momentum)
    gauge_defect = float(np.linalg.norm(horizon_operator @ gauge_data)
                         / max(np.linalg.norm(horizon_operator)
                               * np.linalg.norm(gauge_data), 1.0e-300))
    gauge_coordinates = regular_basis.conj().T @ gauge_data
    physical_coordinates = null_space(gauge_coordinates.conj().T)
    horizon_data = regular_basis @ physical_coordinates
    physical_count = physical_coordinates.shape[1]
    values0 = horizon_data[:field_count]
    derivatives0 = horizon_data[field_count:]
    second0, horizon_recurrence_condition = _horizon_second_coefficient(
        bg, sector, sigma, momentum, values0, derivatives0
    )
    initial = np.vstack((
        values0 + r0 * derivatives0 + r0**2 * second0,
        derivatives0 + 2.0 * r0 * second0,
    ))

    component_count = 2 * field_count
    maximum_c2_condition = 0.0
    identity = np.eye(field_count, dtype=complex)
    zero = np.zeros_like(identity)
    if subdomains is None:
        boundaries = (r0, r_uv)
    else:
        boundaries = (r0, *subdomains, r_uv)
        if any(right <= left for left, right in zip(boundaries, boundaries[1:])):
            raise ValueError("collocation subdomains must be strictly increasing")
    incoming = initial
    residuals = []
    domain_records = []
    for lower, upper in zip(boundaries, boundaries[1:]):
        mapped_nodes, D1, _D2 = radial_grid(intervals, upper - lower, map_alpha)
        nodes = mapped_nodes + lower
        node_count = len(nodes)
        background = background_coefficients(bg, nodes)
        coefficients = local_coefficients(
            sector, background, sigma, momentum, equations
        )
        first_order = np.empty(
            (node_count, component_count, component_count), dtype=complex
        )
        for node_index in range(node_count):
            C0 = coefficients[:, :, 0, node_index]
            C1 = coefficients[:, :, 1, node_index]
            C2 = coefficients[:, :, 2, node_index]
            row_scale = np.maximum(
                np.max(np.abs(C0), axis=1) + np.max(np.abs(C1), axis=1)
                + np.max(np.abs(C2), axis=1), 1.0e-300
            )
            C0, C1, C2 = (item / row_scale[:, None] for item in (C0, C1, C2))
            maximum_c2_condition = max(maximum_c2_condition,
                                       float(np.linalg.cond(C2)))
            first_order[node_index] = np.block([
                [zero, identity],
                [-np.linalg.solve(C2, C0), -np.linalg.solve(C2, C1)],
            ])

        size = component_count * node_count
        operator = np.zeros((size, size), dtype=complex)
        for component in range(component_count):
            rows = slice(component * node_count, (component + 1) * node_count)
            columns = slice(component * node_count, (component + 1) * node_count)
            operator[rows, columns] = D1
            for coupled in range(component_count):
                coupled_columns = slice(coupled * node_count,
                                        (coupled + 1) * node_count)
                operator[rows, coupled_columns] -= np.diag(
                    first_order[:, component, coupled]
                )
        right_hand_side = np.zeros((size, physical_count), dtype=complex)
        for component in range(component_count):
            row = component * node_count
            operator[row] = 0.0
            operator[row, component * node_count] = 1.0
            right_hand_side[row] = incoming[component]
        row_norm = np.maximum(np.sum(np.abs(operator), axis=1), 1.0e-300)
        operator /= row_norm[:, None]
        right_hand_side /= row_norm[:, None]
        solution = np.linalg.solve(operator, right_hand_side)
        residuals.append(float(
            np.linalg.norm(operator @ solution - right_hand_side)
            / max(np.linalg.norm(operator) * np.linalg.norm(solution), 1.0e-300)
        ))
        domain_records.append((nodes, background, solution, first_order))
        incoming = np.vstack([
            solution[component * node_count + node_count - 1]
            for component in range(component_count)
        ])

    values_uv = incoming[:field_count]
    derivatives_uv = incoming[field_count:]
    uv_background = background_coefficients(bg, np.asarray([r_uv]))
    invariant_maps, powers = _invariant_maps(
        sector, uv_background, sigma, momentum
    )
    derivative_step = min(2.0e-5, 0.1 * (r_uv - r0))
    plus_background = background_coefficients(bg, np.asarray([r_uv + derivative_step]))
    minus_background = background_coefficients(bg, np.asarray([r_uv - derivative_step]))
    plus_maps, _ = _invariant_maps(sector, plus_background, sigma, momentum)
    minus_maps, _ = _invariant_maps(sector, minus_background, sigma, momentum)
    Ar_uv = float(uv_background["Ar"][0])
    source_rows = []
    for index, (invariant_map, power) in enumerate(zip(invariant_maps, powers)):
        value = invariant_map @ values_uv
        map_derivative = (plus_maps[index] - minus_maps[index]) / (2.0 * derivative_step)
        derivative = invariant_map @ derivatives_uv + map_derivative @ values_uv
        if sector == "helicity0" and index == len(powers) - 1:
            u_uv = math.exp(-float(uv_background["A"][0]))
            denominator = ((base.DELTA - base.NU) * Ar_uv * u_uv**base.NU)
            source_rows.append((derivative + base.DELTA * Ar_uv * value) / denominator)
        else:
            source_rows.append(value + derivative / (power * Ar_uv))
    source_matrix = np.vstack(source_rows)
    if sector == "helicity1":
        source_scales = np.asarray([
            max(abs(momentum), abs(sigma), 1.0e-12), 1.0
        ])
    else:
        source_scales = np.asarray([
            max(abs(momentum)**2, abs(momentum * sigma), abs(sigma)**2, 1.0e-12),
            max(abs(momentum), abs(sigma), 1.0e-12),
            1.0,
        ])
    normalized_source = source_matrix / source_scales[:, None]
    _left, singular_values, right_h = np.linalg.svd(normalized_source)
    source_vector = right_h[-1].conj()
    diagnostic_names = tuple(spec["diagnostics"])
    diagnostic_maximum = {name: 0.0 for name in diagnostic_names}
    diagnostic_rms_samples = {name: [] for name in diagnostic_names}
    for nodes, domain_background, domain_solution, first_order in domain_records:
        diagnostic_coefficients = local_coefficients(
            sector, domain_background, sigma, momentum, diagnostic_names
        )
        count = len(nodes)
        for node_index in range(1, count - 1):
            state = np.asarray([
                domain_solution[component * count + node_index] @ source_vector
                for component in range(component_count)
            ])
            values, derivatives = state[:field_count], state[field_count:]
            second = first_order[node_index, field_count:] @ state
            for equation_index, name in enumerate(diagnostic_names):
                local = diagnostic_coefficients[equation_index, :, :, node_index]
                terms = (local[:, 0] * values + local[:, 1] * derivatives
                         + local[:, 2] * second)
                state_norm = np.linalg.norm(np.concatenate((values, derivatives, second)))
                relative = float(abs(np.sum(terms))
                                 / max(np.linalg.norm(local) * state_norm, 1.0e-300))
                diagnostic_maximum[name] = max(diagnostic_maximum[name], relative)
                diagnostic_rms_samples[name].append(relative)
    diagnostic_rms = {
        name: float(np.sqrt(np.mean(np.square(samples))))
        for name, samples in diagnostic_rms_samples.items()
    }
    return {
        "omega": omega, "qhat": qhat, "sector": sector,
        "source_matrix": source_matrix,
        "normalized_source_matrix": normalized_source,
        "source_singular_values": singular_values,
        "source_singular_value": float(singular_values[-1]),
        "normalized_determinant": complex(np.linalg.det(normalized_source)),
        "source_vector": source_vector,
        "horizon_gauge_defect": gauge_defect,
        "horizon_recurrence_condition": horizon_recurrence_condition,
        "maximum_c2_condition": maximum_c2_condition,
        "linear_residual": max(residuals),
        "intervals": intervals, "r0": r0, "r_uv": r_uv,
        "map_alpha": map_alpha,
        "subdomain_boundaries": boundaries,
        "horizon_einstein_constraint": horizon_einstein_constraint,
        "diagnostic_relative_maximum": diagnostic_maximum,
        "diagnostic_relative_rms": diagnostic_rms,
    }


def quotient_matrix(bg: dict[str, object], sector: str, qhat: float,
                    omega: complex, intervals: int = 48, r_cut: float = 8.0,
                    map_alpha: float = 0.0,
                    append_diagnostics: bool = True) -> dict[str, object]:
    """Gauge-quotiented primitive operator at a fixed complex frequency.

    The physical UV rows impose source freedom directly on gauge invariants.
    Four (one in helicity one) global orthogonality rows select a representative
    transverse to the analytically known residual radial-gauge orbits.  This
    removes the near-null gauge directions that make the primitive QZ pencil
    ill-conditioned at hydrodynamic frequency.  By default the independent
    primitive radial equations are appended on the open interval.  The result
    is deliberately rectangular: a candidate is null only when it satisfies
    evolution, UV sources, gauge fixing, and every constraint simultaneously,
    without making an arbitrary choice of a Bianchi-redundant Einstein row.
    """
    spec = SECTOR_SPEC[sector]
    equations = tuple(spec["equations"])
    nodes, D1, D2 = radial_grid(intervals, r_cut, map_alpha)
    background = background_coefficients(bg, nodes)
    h0 = base.extract_uv(bg)["h0"]
    k_numeric = qhat / (2.0 * math.sqrt(h0))
    momentum = 1.0j * k_numeric
    sigma = -1.0j * omega
    coefficients = local_coefficients(sector, background, sigma, momentum, equations)
    matrix = _local_operator(coefficients, D1, D2)
    n = len(nodes)

    for equation_index, constraint_name in spec.get("horizon_constraints", {}).items():
        constraint = local_coefficients(sector, background, sigma, momentum,
                                        (constraint_name,))
        constraint_operator = _local_operator(constraint, D1, D2)
        matrix[equation_index * n] = constraint_operator[0]

    invariant_maps, powers = _invariant_maps(sector, background, sigma, momentum)
    boundary_derivative = D1[-1]
    uv_selector = np.zeros(n)
    uv_selector[-1] = 1.0
    physical_rows = [
        (boundary_derivative + power * background["Ar"][-1] * uv_selector) @ invariant_map
        for invariant_map, power in zip(invariant_maps, powers)
    ]
    gauge_vectors = _residual_gauge_vectors(sector, nodes, background, sigma, momentum)
    gauge_basis, _ = np.linalg.qr(gauge_vectors, mode="reduced")
    gauge_rows = [gauge_basis[:, index].conj() for index in range(gauge_basis.shape[1])]
    replacement_rows = physical_rows + gauge_rows
    uv_rows = [equation_index * n + n - 1 for equation_index in range(len(equations))]
    if len(replacement_rows) != len(uv_rows):
        raise RuntimeError("physical/gauge boundary row count does not close the system")
    for row_index, row in zip(uv_rows, replacement_rows):
        matrix[row_index] = row

    square_row_count = matrix.shape[0]
    if append_diagnostics:
        diagnostic_names = tuple(spec["diagnostics"])
        diagnostic_coefficients = local_coefficients(
            sector, background, sigma, momentum, diagnostic_names
        )
        diagnostic_operator = _local_operator(diagnostic_coefficients, D1, D2)
        interior_rows = []
        for equation_index in range(len(diagnostic_names)):
            start = equation_index * n
            interior_rows.append(diagnostic_operator[start + 1:start + n - 1])
        matrix = np.vstack((matrix, *interior_rows))

    row_norm = np.maximum(np.sum(np.abs(matrix), axis=1), 1.0e-300)
    matrix /= row_norm[:, None]
    return {"matrix": matrix, "nodes": nodes, "D1": D1, "D2": D2,
            "background": background, "sigma": sigma, "omega": omega,
            "momentum": momentum, "k_numeric": k_numeric, "qhat": qhat,
            "sector": sector, "gauge_vectors": gauge_vectors,
            "square_row_count": square_row_count,
            "diagnostics_appended": bool(append_diagnostics)}


def quotient_singular_value(bg: dict[str, object], sector: str, qhat: float,
                            omega: complex, intervals: int = 48,
                            r_cut: float = 8.0, map_alpha: float = 0.0,
                            return_vector: bool = False):
    result = quotient_matrix(bg, sector, qhat, omega, intervals, r_cut, map_alpha)
    _left, singular_values, right_h = np.linalg.svd(result["matrix"], full_matrices=False)
    value = float(singular_values[-1])
    if return_vector:
        return value, right_h[-1].conj(), result
    return value


def quotient_mode(bg: dict[str, object], sector: str, qhat: float,
                  initial: complex, intervals: int = 48, r_cut: float = 8.0,
                  map_alpha: float = 0.0, search_radius: float = 0.05,
                  max_iterations: int = 100) -> dict[str, object]:
    """Refine a physical QNM by minimizing the gauge-quotiented source singular value."""
    start = np.asarray([initial.real, initial.imag], dtype=float)

    def objective(value: np.ndarray) -> float:
        if np.linalg.norm(value - start) > search_radius:
            return 10.0 + float(np.linalg.norm(value - start))
        singular = quotient_singular_value(
            bg, sector, qhat, complex(value[0], value[1]), intervals, r_cut, map_alpha
        )
        return math.log10(max(singular, 1.0e-300))

    fit = minimize(objective, start, method="Nelder-Mead",
                   options={"maxiter": max_iterations, "xatol": 2.0e-9,
                            "fatol": 2.0e-6, "adaptive": True})
    omega = complex(fit.x[0], fit.x[1])
    singular, vector, result = quotient_singular_value(
        bg, sector, qhat, omega, intervals, r_cut, map_alpha, return_vector=True
    )
    return {**result, "omega": omega, "right": vector,
            "source_singular_value": singular, "optimizer_success": bool(fit.success),
            "optimizer_message": str(fit.message), "evaluations": int(fit.nfev),
            "distance_from_initial": float(abs(omega - initial))}
