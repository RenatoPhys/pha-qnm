"""Reference-compatible PHA thermodynamics and critical-point search.

This module follows the public MUSES implementation at commit
52b093230a8b5f36a83fc37c6fb5e8b2043c82c4 (2023-08-24), the revision
contemporaneous with the ``sample74`` HDF5 timestamp.  It deliberately keeps
the reference horizon normalization, relaxational UV field, and the
constant-phi_H line-crossing definition of the critical point separate from
the QNM calculation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


@dataclass(frozen=True)
class PHAParameters:
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


@dataclass(frozen=True)
class BackgroundNumerics:
    relaxation_rate: float = 8.0
    ricci_tolerance: float = 1.2e-3
    phiA_tolerance: float = 1.0e-4
    initial_step: float = 2.0e-9
    maximum_radius: float = 80.0
    rtol: float = 2.0e-10
    atol: float = 2.0e-12
    max_step: float = 0.05


@dataclass(frozen=True)
class CriticalPointNumerics:
    phi0_min: float = 0.5
    phi0_max: float = 7.0
    line_count: int = 800
    maximum_iterations: int = 1000
    initial_charge_step: float = 1.0e-4
    maximum_charge_step: float = 1.0
    relative_T_tolerance: float = 1.0e-5
    relative_mu_tolerance: float = 1.0e-5
    absolute_T_tolerance_MeV: float = 1.0e-3
    absolute_mu_tolerance_MeV: float = 1.0e-3


DEFAULT_PARAMETERS = PHAParameters()


def _sech(value: float | np.ndarray) -> float | np.ndarray:
    value = np.asarray(value, dtype=float)
    exponential = np.exp(-np.abs(value))
    result = 2.0 * exponential / (1.0 + exponential * exponential)
    return float(result) if result.ndim == 0 else result


class PHAModel:
    def __init__(self, parameters: PHAParameters = DEFAULT_PARAMETERS):
        self.parameters = parameters
        mass2 = -12.0 * parameters.gamma**2 + 2.0 * parameters.b2
        self.Delta = 2.0 + math.sqrt(4.0 + mass2)
        self.nu = 4.0 - self.Delta

    def V(self, phi: float) -> float:
        p = self.parameters
        return (-12.0 * math.cosh(p.gamma * phi) + p.b2 * phi**2
                + p.b4 * phi**4 + p.b6 * phi**6)

    def Vp(self, phi: float) -> float:
        p = self.parameters
        return (-12.0 * p.gamma * math.sinh(p.gamma * phi) + 2.0 * p.b2 * phi
                + 4.0 * p.b4 * phi**3 + 6.0 * p.b6 * phi**5)

    def f(self, phi: float) -> float:
        p = self.parameters
        polynomial = p.c1 * phi + p.c2 * phi**2 + p.c3 * phi**3
        return (_sech(polynomial) + p.d1 * _sech(p.d2 * phi)) / (1.0 + p.d1)

    def fp(self, phi: float) -> float:
        p = self.parameters
        polynomial = p.c1 * phi + p.c2 * phi**2 + p.c3 * phi**3
        polynomial_prime = p.c1 + 2.0 * p.c2 * phi + 3.0 * p.c3 * phi**2
        return (-_sech(polynomial) * math.tanh(polynomial) * polynomial_prime
                - p.d1 * p.d2 * _sech(p.d2 * phi) * math.tanh(p.d2 * phi)) / (1.0 + p.d1)

    def maximum_electric_slope(self, phi0: float) -> float:
        return math.sqrt(-2.0 * self.V(phi0) / self.f(phi0))


def _ricci_scalar(model: PHAModel, state: np.ndarray) -> float:
    A, Ap, h, hp, phi, phip, _Phi, Phip, _C, _D, _E = state
    App = -phip * phip / 6.0
    hpp = -4.0 * Ap * hp + math.exp(-2.0 * A) * model.f(phi) * Phip * Phip
    return -(20.0 * h * Ap * Ap + 9.0 * Ap * hp + 8.0 * h * App + hpp)


def solve_background(
    phi0: float,
    charge_fraction: float,
    model: PHAModel | None = None,
    numerics: BackgroundNumerics = BackgroundNumerics(),
) -> dict:
    """Integrate one black brane using the public MUSES UV prescription."""
    model = model or PHAModel()
    if not 0.0 <= charge_fraction < 1.0:
        raise ValueError("charge_fraction must lie in [0, 1)")
    Phi1 = charge_fraction * model.maximum_electric_slope(phi0)
    f0 = model.f(phi0)
    fp0 = model.fp(phi0)
    V0 = model.V(phi0)
    Vp0 = model.Vp(phi0)

    phi1 = (-Phi1 * Phi1 * fp0 + 2.0 * Vp0) / 2.0
    A1 = (-f0 * Phi1 * Phi1 - 2.0 * V0) / 6.0
    A2 = -(Phi1 * Phi1 * fp0 - 2.0 * Vp0) ** 2 / 48.0
    h2 = (5.0 * f0 * Phi1 * Phi1 + 4.0 * V0) / 6.0
    Phi2 = Phi1 * (2.0 * f0 * Phi1 * Phi1 + 4.0 * V0
                   + 3.0 * Phi1 * Phi1 * fp0 * fp0 / f0
                   - 6.0 * fp0 * Vp0 / f0) / 12.0

    eps = numerics.initial_step
    y0 = np.array([
        A1 * eps + A2 * eps**2,
        A1 + 2.0 * A2 * eps,
        eps + h2 * eps**2,
        1.0 + 2.0 * h2 * eps,
        phi0 + phi1 * eps,
        phi1,
        Phi1 * eps + Phi2 * eps**2,
        Phi1 + 2.0 * Phi2 * eps,
        phi0,
        Phi1,
        0.0,
    ], dtype=float)

    def rhs(_radius: float, state: np.ndarray) -> np.ndarray:
        A, Ap, h, hp, phi, phip, _Phi, Phip, C, D, _E = state
        fv = model.f(phi)
        fpv = model.fp(phi)
        em2A = math.exp(-2.0 * A)
        return np.array([
            Ap,
            -phip * phip / 6.0,
            hp,
            -4.0 * Ap * hp + em2A * fv * Phip * Phip,
            phip,
            -(hp / h + 4.0 * Ap) * phip
            + (model.Vp(phi) - 0.5 * em2A * Phip * Phip * fpv) / h,
            Phip,
            -(2.0 * Ap + fpv * phip / fv) * Phip,
            numerics.relaxation_rate * (phi * math.exp(model.nu * A) - C),
            numerics.relaxation_rate * (Phip * math.exp(2.0 * A) - D),
            em2A / fv,
        ])

    def convergence_event(_radius: float, state: np.ndarray) -> float:
        A, _Ap, _h, _hp, phi, _phip, _Phi, _Phip, C, _D, _E = state
        phi_precision = abs((C - phi * math.exp(model.nu * A)) / max(abs(C), 1.0e-300))
        ricci_precision = abs(-_ricci_scalar(model, state) / 20.0 - 1.0)
        return max(phi_precision / numerics.phiA_tolerance,
                   ricci_precision / numerics.ricci_tolerance) - 1.0

    convergence_event.terminal = True
    convergence_event.direction = -1.0
    solution = solve_ivp(
        rhs, (eps, numerics.maximum_radius), y0, method="DOP853",
        rtol=numerics.rtol, atol=numerics.atol, max_step=numerics.max_step,
        events=convergence_event,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    state = solution.y[:, -1]
    A, Ap, h0, hp, phi, phip, Phi, Phip, phiA, D, E = state
    ricci_precision = abs(-_ricci_scalar(model, state) / 20.0 - 1.0)
    phiA_precision = abs((phiA - phi * math.exp(model.nu * A)) / max(abs(phiA), 1.0e-300))
    converged = bool(ricci_precision <= numerics.ricci_tolerance * (1.0 + 1.0e-8)
                     and phiA_precision <= numerics.phiA_tolerance * (1.0 + 1.0e-8))
    if not converged:
        raise RuntimeError(
            f"UV convergence failed at r={solution.t[-1]:.6g}: "
            f"Ricci={ricci_precision:.3e}, phiA={phiA_precision:.3e}"
        )

    Phi2_far = -math.sqrt(h0) * f0 * Phi1 / (2.0 * model.f(0.0))
    Phi0_far = Phi - Phi2_far * math.exp(-2.0 * A)
    scale = phiA ** (1.0 / model.nu)
    p = model.parameters
    temperature = p.Lambda_MeV / (4.0 * math.pi * scale * math.sqrt(h0))
    chemical_potential = p.Lambda_MeV * Phi0_far / (scale * math.sqrt(h0))
    entropy = p.Lambda_MeV**3 * 2.0 * math.pi / (phiA ** (3.0 / model.nu) * p.kappa2)
    density = -p.Lambda_MeV**3 * Phi2_far / (phiA ** (3.0 / model.nu)
                                               * math.sqrt(h0) * p.kappa2)
    constraint = (h0 * (24.0 * Ap * Ap - phip * phip) + 6.0 * Ap * hp
                  + 2.0 * model.V(phi) + math.exp(-2.0 * A) * model.f(phi) * Phip * Phip) / (24.0 * h0)
    charge_reference = f0 * Phi1
    charge_uv = model.f(phi) * math.exp(2.0 * A) * Phip
    return {
        "phi0": float(phi0),
        "charge_fraction": float(charge_fraction),
        "Phi1": float(Phi1),
        "T_MeV": float(temperature),
        "mu_MeV": float(chemical_potential),
        "s_MeV3": float(entropy),
        "rho_MeV3": float(density),
        "r_uv": float(solution.t[-1]),
        "A_uv": float(A),
        "A_slope_uv": float(Ap),
        "h0": float(h0),
        "phiA": float(phiA),
        "Phi0": float(Phi0_far),
        "Phi2": float(Phi2_far),
        "relaxation_D": float(D),
        "relaxation_E": float(E),
        "ricci_precision": float(ricci_precision),
        "phiA_precision": float(phiA_precision),
        "constraint": float(constraint),
        "gauss_relative_drift": float(charge_uv / charge_reference - 1.0) if charge_reference else 0.0,
        "nu_relative_error": float(-phip / (phi * Ap) / model.nu - 1.0),
        "converged": converged,
        "parameters": asdict(p),
    }


class _QuadraticWindow:
    def __init__(self, x: list[float], y: list[float]):
        self.x = list(x)
        self.y = list(y)
        self._fit()

    def _fit(self) -> None:
        self.coefficients = np.polyfit(self.x, self.y, 2)

    def __call__(self, value: float) -> float:
        return float(np.polyval(self.coefficients, value))

    def update(self, x: float, y: float) -> None:
        self.x = self.x[1:] + [float(x)]
        self.y = self.y[1:] + [float(y)]
        self._fit()


def _closest_positive_crossing(first: _QuadraticWindow, second: _QuadraticWindow) -> float:
    coefficients = first.coefficients - second.coefficients
    roots = np.roots(coefficients)
    positive = [float(root.real) for root in roots
                if abs(root.imag) <= 1.0e-8 * max(1.0, abs(root.real)) and root.real > 0.0]
    return min(positive) if positive else math.nan


def _equally_spaced_phi0_tracks(
    evaluate: Callable[[float, float], dict], options: CriticalPointNumerics,
) -> tuple[np.ndarray, np.ndarray]:
    # This reproduces lines 90--183 of the 2023 MUSES implementation,
    # including its deliberately modest half-spacing acceptance criterion.
    state_at_max_phi = evaluate(options.phi0_max, 0.0)
    state = evaluate(options.phi0_min, 0.0)
    Tmin = state_at_max_phi["T_MeV"]
    Tmax = state["T_MeV"]
    T_step = (Tmax - Tmin) / options.line_count
    exponential_scale = (options.phi0_max - options.phi0_min) / math.log(Tmax / Tmin)
    phi_step = -0.01 * math.log(1.0 - T_step / Tmax) / exponential_scale
    phi_tracks = np.empty(options.line_count)
    temperatures = np.empty(options.line_count)
    track = 0
    attempts = 0
    phi0 = options.phi0_min
    previous_T = Tmax
    while track < options.line_count:
        target_T = Tmax - track * T_step
        if state["T_MeV"] != target_T:
            state = evaluate(phi0, 0.0)
            denominator = state["T_MeV"] - previous_T
            if denominator != 0.0:
                dphi_dT = phi_step / denominator
                phi_step = dphi_dT * (target_T - state["T_MeV"])
        if (abs(state["T_MeV"] - target_T) < T_step / 2.0
                or attempts > options.maximum_iterations):
            attempts = 0
            phi_tracks[track] = state["phi0"]
            temperatures[track] = state["T_MeV"]
            track += 1
        else:
            attempts += 1
        phi0 += phi_step
        previous_T = state["T_MeV"]
    return phi_tracks, temperatures


def locate_critical_point(
    model: PHAModel | None = None,
    background_numerics: BackgroundNumerics = BackgroundNumerics(),
    options: CriticalPointNumerics = CriticalPointNumerics(),
    progress: Callable[[str], None] | None = None,
    include_non_neighboring: bool = True,
) -> dict:
    """Locate the first crossing of constant-phi0 curves as in MUSES 2023."""
    model = model or PHAModel()
    cache: dict[tuple[float, float], dict] = {}

    def evaluate(phi0: float, charge: float) -> dict:
        key = (float(phi0), float(charge))
        if key not in cache:
            cache[key] = solve_background(key[0], key[1], model, background_numerics)
        return cache[key]

    phi_tracks, zero_mu_T = _equally_spaced_phi0_tracks(evaluate, options)
    T_windows: list[_QuadraticWindow] = []
    charge_windows: list[_QuadraticWindow] = []
    last_charges: list[float] = []
    for index, (phi0, T0) in enumerate(zip(phi_tracks, zero_mu_T)):
        charges = [0.0, options.initial_charge_step, 2.0 * options.initial_charge_step]
        states = [evaluate(float(phi0), charge) for charge in charges]
        mus = [state["mu_MeV"] for state in states]
        temperatures = [float(T0), states[1]["T_MeV"], states[2]["T_MeV"]]
        T_windows.append(_QuadraticWindow(mus, temperatures))
        charge_windows.append(_QuadraticWindow(mus, charges))
        last_charges.append(charges[-1])
        if progress and (index + 1) % max(1, options.line_count // 20) == 0:
            progress(f"initialized {index + 1}/{options.line_count} constant-phi0 lines")

    crossings: list[dict] = []
    good_tracks: list[int] = []

    def refine_pair(upper_index: int, lower_index: int) -> tuple[dict, dict] | None:
        T_upper = _QuadraticWindow(T_windows[upper_index].x, T_windows[upper_index].y)
        T_lower = _QuadraticWindow(T_windows[lower_index].x, T_windows[lower_index].y)
        q_upper = _QuadraticWindow(charge_windows[upper_index].x, charge_windows[upper_index].y)
        q_lower = _QuadraticWindow(charge_windows[lower_index].x, charge_windows[lower_index].y)
        previous_upper = last_charges[upper_index]
        previous_lower = last_charges[lower_index]
        for _iteration in range(options.maximum_iterations):
            candidate_mu = _closest_positive_crossing(T_upper, T_lower)
            if not math.isfinite(candidate_mu):
                return None
            charge_upper = q_upper(candidate_mu)
            charge_lower = q_lower(candidate_mu)
            if charge_upper < 0.0 or charge_lower < 0.0:
                return None
            charge_upper = min(charge_upper, 0.99, previous_upper + options.maximum_charge_step)
            charge_lower = min(charge_lower, 0.99, previous_lower + options.maximum_charge_step)
            previous_upper, previous_lower = charge_upper, charge_lower
            try:
                state_upper = evaluate(float(phi_tracks[upper_index]), charge_upper)
                state_lower = evaluate(float(phi_tracks[lower_index]), charge_lower)
            except (RuntimeError, ValueError, OverflowError):
                return None
            T_upper.update(state_upper["mu_MeV"], state_upper["T_MeV"])
            q_upper.update(state_upper["mu_MeV"], charge_upper)
            T_lower.update(state_lower["mu_MeV"], state_lower["T_MeV"])
            q_lower.update(state_lower["mu_MeV"], charge_lower)
            success_T = abs(state_upper["T_MeV"] - state_lower["T_MeV"]) <= (
                options.absolute_T_tolerance_MeV
                + options.relative_T_tolerance * state_upper["T_MeV"])
            success_mu = abs(state_upper["mu_MeV"] - state_lower["mu_MeV"]) <= (
                options.absolute_mu_tolerance_MeV
                + options.relative_mu_tolerance * state_upper["mu_MeV"])
            if success_T and success_mu:
                return state_lower, state_upper
        return None

    for upper in range(1, options.line_count):
        result = refine_pair(upper, upper - 1)
        if result is not None:
            crossings.extend(result)
            good_tracks.append(upper)
        if progress and (upper + 1) % max(1, options.line_count // 20) == 0:
            progress(f"tested {upper + 1}/{options.line_count} neighboring pairs")

    neighboring_crossing_count = len(crossings)
    neighboring_best = min(crossings, key=lambda state: state["mu_MeV"]) if crossings else None
    if include_non_neighboring:
        for upper in good_tracks:
            for lower in good_tracks:
                if lower >= upper or lower == upper - 1:
                    continue
                result = refine_pair(upper, lower)
                if result is not None:
                    crossings.extend(result)

    if not crossings:
        raise RuntimeError("no constant-phi0 line crossing was found")
    crossings.sort(key=lambda state: state["mu_MeV"])
    best = crossings[0]
    coarse = crossings[3::4]
    coarse_best = min(coarse, key=lambda state: state["mu_MeV"]) if coarse else crossings[1]
    return {
        "status": "success",
        "T_c_MeV": best["T_MeV"],
        "mu_B_c_MeV": best["mu_MeV"],
        "error_T_c_MeV": abs(best["T_MeV"] - coarse_best["T_MeV"]),
        "error_mu_B_c_MeV": abs(best["mu_MeV"] - coarse_best["mu_MeV"]),
        "critical_background": best,
        "crossing_count": len(crossings),
        "neighboring_crossing_count": neighboring_crossing_count,
        "neighboring_minimum": ({"T_MeV": neighboring_best["T_MeV"],
                                  "mu_MeV": neighboring_best["mu_MeV"],
                                  "phi0": neighboring_best["phi0"],
                                  "charge_fraction": neighboring_best["charge_fraction"]}
                                 if neighboring_best else None),
        "good_neighboring_tracks": len(good_tracks),
        "included_non_neighboring": include_non_neighboring,
        "background_evaluations": len(cache),
        "options": asdict(options),
        "background_numerics": asdict(background_numerics),
    }


def locate_cusp_critical_point(
    initial_guess: tuple[float, float] = (3.68, 0.32),
    model: PHAModel | None = None,
    background_numerics: BackgroundNumerics = BackgroundNumerics(),
    phi_step: float = 1.0e-3,
    charge_step: float = 1.0e-4,
) -> dict:
    """Independently locate the cusp of the horizon-to-thermodynamics map.

    A spinodal is a fold where ``det d(log T, log mu)/d(phi0,q)`` vanishes.
    At its endpoint the null direction of that Jacobian is tangent to the fold,
    so the directional derivative of the determinant also vanishes.  These two
    local equations define the cusp without using the HDF5 critical-point
    coordinates as a target.
    """
    model = model or PHAModel()
    cache: dict[tuple[float, float], dict] = {}

    def evaluate(point: np.ndarray) -> dict:
        key = (float(point[0]), float(point[1]))
        if key not in cache:
            cache[key] = solve_background(key[0], key[1], model, background_numerics)
        return cache[key]

    steps = np.array([phi_step, charge_step])

    def jacobian(point: np.ndarray) -> np.ndarray:
        matrix = np.empty((2, 2), dtype=float)
        for column in range(2):
            plus = point.copy(); plus[column] += steps[column]
            minus = point.copy(); minus[column] -= steps[column]
            state_plus = evaluate(plus)
            state_minus = evaluate(minus)
            matrix[:, column] = [
                math.log(state_plus["T_MeV"] / state_minus["T_MeV"]),
                math.log(state_plus["mu_MeV"] / state_minus["mu_MeV"]),
            ]
            matrix[:, column] /= 2.0 * steps[column]
        return matrix

    def determinant(point: np.ndarray) -> float:
        return float(np.linalg.det(jacobian(point)))

    def equations(point: np.ndarray) -> np.ndarray:
        J = jacobian(point)
        detJ = float(np.linalg.det(J))
        _u, singular_values, vh = np.linalg.svd(J)
        null_direction = vh[-1]
        gradient = np.empty(2, dtype=float)
        for axis in range(2):
            plus = point.copy(); plus[axis] += steps[axis]
            minus = point.copy(); minus[axis] -= steps[axis]
            gradient[axis] = (determinant(plus) - determinant(minus)) / (2.0 * steps[axis])
        directional = float(gradient @ null_direction)
        # Fixed scales only condition the nonlinear solve; they do not change
        # either zero.  Values are representative of the MAP cusp region.
        return np.array([detJ / 0.2, directional / 2.0])

    fit = least_squares(
        equations, np.asarray(initial_guess, dtype=float),
        bounds=([2.5, 0.15], [5.0, 0.55]),
        xtol=2.0e-11, ftol=2.0e-11, gtol=2.0e-11,
        max_nfev=120, verbose=0,
    )
    point = fit.x
    state = evaluate(point)
    J = jacobian(point)
    singular_values = np.linalg.svd(J, compute_uv=False)
    residual = equations(point)
    return {
        "status": "success" if fit.success and np.linalg.norm(residual) < 2.0e-5 else "failed",
        "phi0": float(point[0]),
        "charge_fraction": float(point[1]),
        "T_c_MeV": state["T_MeV"],
        "mu_B_c_MeV": state["mu_MeV"],
        "critical_background": state,
        "jacobian_log_thermodynamics": J.tolist(),
        "jacobian_singular_values": singular_values.tolist(),
        "determinant": float(np.linalg.det(J)),
        "cusp_equation_residual": residual.tolist(),
        "optimizer_success": bool(fit.success),
        "optimizer_status": int(fit.status),
        "optimizer_message": str(fit.message),
        "optimizer_optimality": float(fit.optimality),
        "active_bounds": fit.active_mask.tolist(),
        "function_evaluations": int(fit.nfev),
        "background_evaluations": len(cache),
        "finite_difference_steps": {"phi0": phi_step, "charge_fraction": charge_step},
        "background_numerics": asdict(background_numerics),
    }
