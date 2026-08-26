#!/usr/bin/env python3
"""Independent charged-hydrodynamic coefficients for the PHA EMD model."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

import run_numerics as base
from pha_qnm.thermodynamics import (
    BackgroundNumerics, PHAModel, solve_background as reference_background,
)


TIGHT_BACKGROUND_NUMERICS = BackgroundNumerics(
    ricci_tolerance=1.0e-5,
    phiA_tolerance=1.0e-5,
)


@dataclass(frozen=True)
class ChargedHydroState:
    phi0: float
    charge_fraction: float
    T_MeV: float
    mu_MeV: float
    entropy_MeV3: float
    density_MeV3: float
    enthalpy_MeV4: float
    pressure_energy: float
    pressure_density_MeV: float
    alpha_energy_MeV4_inverse: float
    alpha_density_MeV3_inverse: float
    eta_MeV3: float
    zeta_over_eta: float
    sigma_Q_over_T: float
    conductivity_times_T_MeV2: float
    chi_B_over_T2: float
    thermo_jacobian_condition: float
    maxwell_relation_relative_error: float
    finite_difference_phi0: float
    finite_difference_charge: float
    backend: str

    def as_dict(self) -> dict[str, float | str]:
        return asdict(self)


def _physical_state(phi0: float, charge_fraction: float,
                    backend: str,
                    model: PHAModel) -> dict[str, float]:
    if backend == "reference":
        state = reference_background(phi0, charge_fraction, model=model)
        return {
            "T": float(state["T_MeV"]),
            "mu": float(state["mu_MeV"]),
            "s": float(state["s_MeV3"]),
            "n": float(state["rho_MeV3"]),
            "h0": float(state["h0"]),
        }
    if backend == "reference_tight":
        state = reference_background(
            phi0, charge_fraction, model=model,
            numerics=TIGHT_BACKGROUND_NUMERICS,
        )
        return {
            "T": float(state["T_MeV"]),
            "mu": float(state["mu_MeV"]),
            "s": float(state["s_MeV3"]),
            "n": float(state["rho_MeV3"]),
            "h0": float(state["h0"]),
        }
    if backend == "fast":
        if model.parameters != PHAModel().parameters:
            raise ValueError("the fast backend is available only for the MAP model")
        uv = base.extract_uv(base.integrate_background(phi0, charge_fraction))
        temperature = float(uv["T_MeV"])
        return {
            "T": temperature,
            "mu": float(uv["mu_MeV"]),
            "s": float(uv["s_over_T3"] * temperature**3),
            "n": float(uv["rho_over_T3"] * temperature**3),
            "h0": float(uv["h0"]),
        }
    raise ValueError("backend must be 'reference', 'reference_tight', or 'fast'")


def charged_hydro_state(phi0: float, charge_fraction: float,
                        phi_step: float = 1.0e-3,
                        charge_step: float = 1.0e-4,
                        backend: str = "reference",
                        model: PHAModel | None = None) -> ChargedHydroState:
    """Construct the local EOS and horizon transport matrix.

    The EOS derivatives are transformed from horizon coordinates to
    ``(T, mu_B)``.  Bulk viscosity uses the charged Eling--Oz horizon formula,
    and ``sigma_Q`` is the finite incoherent conductivity of the translationally
    invariant EMD plasma.
    """
    model = model or PHAModel()
    center = _physical_state(phi0, charge_fraction, backend, model)
    jacobian = np.empty((4, 2), dtype=float)
    for column, step in enumerate((phi_step, charge_step)):
        if column == 1 and charge_fraction - step < 0.0:
            state_step = _physical_state(phi0, charge_fraction + step, backend, model)
            state_twice = _physical_state(
                phi0, charge_fraction + 2.0 * step, backend, model
            )
            vector_step = np.asarray([
                state_step["T"], state_step["mu"], state_step["s"], state_step["n"]
            ])
            vector_twice = np.asarray([
                state_twice["T"], state_twice["mu"],
                state_twice["s"], state_twice["n"]
            ])
            # T and s are even, while mu and n are odd under charge reversal.
            jacobian[:, column] = 0.0
            jacobian[[1, 3], column] = (
                8.0 * vector_step[[1, 3]] - vector_twice[[1, 3]]
            ) / (6.0 * step)
            continue
        plus = [phi0, charge_fraction]
        minus = [phi0, charge_fraction]
        plus[column] += step
        minus[column] -= step
        state_plus = _physical_state(*plus, backend, model)
        state_minus = _physical_state(*minus, backend, model)
        vector_plus = np.asarray([
            state_plus["T"], state_plus["mu"], state_plus["s"], state_plus["n"]
        ])
        vector_minus = np.asarray([
            state_minus["T"], state_minus["mu"], state_minus["s"], state_minus["n"]
        ])
        jacobian[:, column] = (vector_plus - vector_minus) / (2.0 * step)

    T, mu, entropy, density = (
        center["T"], center["mu"], center["s"], center["n"]
    )
    horizon_to_thermo = jacobian[:2]
    entropy_density_Tmu = jacobian[2:] @ np.linalg.inv(horizon_to_thermo)
    maxwell_error = abs(entropy_density_Tmu[0, 1] - entropy_density_Tmu[1, 0]) / max(
        abs(entropy_density_Tmu[0, 1]), abs(entropy_density_Tmu[1, 0]), 1.0
    )
    energy_density_Tmu = np.vstack((
        np.asarray([T, mu]) @ entropy_density_Tmu,
        entropy_density_Tmu[1],
    ))
    thermo_to_energy_density = np.linalg.inv(energy_density_Tmu)
    pressure_gradient = np.asarray([entropy, density]) @ thermo_to_energy_density
    alpha_gradient = np.asarray([-mu / T**2, 1.0 / T]) @ thermo_to_energy_density

    entropy_density_to_horizon = np.linalg.inv(jacobian[2:])
    phi_gradient = np.asarray([1.0, 0.0]) @ entropy_density_to_horizon
    zeta_over_eta = float(
        (entropy * phi_gradient[0] + density * phi_gradient[1])**2
    )
    enthalpy = T * entropy + mu * density
    eta = entropy / (4.0 * math.pi)
    sigma_Q_over_T = (
        2.0 * math.pi * math.sqrt(center["h0"]) * float(model.f(phi0))
        / model.parameters.kappa2 * (T * entropy / enthalpy)**2
    )
    return ChargedHydroState(
        phi0=float(phi0), charge_fraction=float(charge_fraction),
        T_MeV=T, mu_MeV=mu, entropy_MeV3=entropy, density_MeV3=density,
        enthalpy_MeV4=enthalpy,
        pressure_energy=float(pressure_gradient[0]),
        pressure_density_MeV=float(pressure_gradient[1]),
        alpha_energy_MeV4_inverse=float(alpha_gradient[0]),
        alpha_density_MeV3_inverse=float(alpha_gradient[1]),
        eta_MeV3=float(eta), zeta_over_eta=zeta_over_eta,
        sigma_Q_over_T=float(sigma_Q_over_T),
        conductivity_times_T_MeV2=float(sigma_Q_over_T * T**2),
        chi_B_over_T2=float(entropy_density_Tmu[1, 1] / T**2),
        thermo_jacobian_condition=float(np.linalg.cond(horizon_to_thermo)),
        maxwell_relation_relative_error=float(maxwell_error),
        finite_difference_phi0=float(phi_step),
        finite_difference_charge=float(charge_step), backend=backend,
    )


def hydrodynamic_modes(state: ChargedHydroState,
                       qhat: float) -> np.ndarray:
    """Return the three charged-hydrodynamic roots in ``omega/(2 pi T)``."""
    T = state.T_MeV
    momentum = 2.0 * math.pi * T * qhat
    enthalpy = state.enthalpy_MeV4
    longitudinal_viscosity = state.eta_MeV3 * (4.0 / 3.0 + state.zeta_over_eta)
    evolution = np.zeros((3, 3), dtype=complex)
    evolution[0, 2] = -1.0j * momentum * enthalpy
    evolution[1, 0] = (
        -state.conductivity_times_T_MeV2 * momentum**2
        * state.alpha_energy_MeV4_inverse
    )
    evolution[1, 1] = (
        -state.conductivity_times_T_MeV2 * momentum**2
        * state.alpha_density_MeV3_inverse
    )
    evolution[1, 2] = -1.0j * momentum * state.density_MeV3
    evolution[2, 0] = -1.0j * momentum * state.pressure_energy / enthalpy
    evolution[2, 1] = -1.0j * momentum * state.pressure_density_MeV / enthalpy
    evolution[2, 2] = -longitudinal_viscosity * momentum**2 / enthalpy
    omega_hat = 1.0j * np.linalg.eigvals(evolution) / (2.0 * math.pi * T)
    return omega_hat[np.lexsort((omega_hat.imag, omega_hat.real))]


def ideal_sound_speed_squared(state: ChargedHydroState) -> float:
    return float(
        state.pressure_energy
        + state.density_MeV3 / state.enthalpy_MeV4 * state.pressure_density_MeV
    )


def shear_diffusion_hat(state: ChargedHydroState) -> float:
    """Return ``(2 pi T) eta/(epsilon+P)``."""
    return float(2.0 * math.pi * state.T_MeV * state.eta_MeV3
                 / state.enthalpy_MeV4)


def baryon_diffusion_hat(state: ChargedHydroState,
                         probe_qhat: float = 1.0e-4) -> float:
    """Return the signed coefficient in ``omega_D=-i D_hat qhat^2``."""
    modes = hydrodynamic_modes(state, probe_qhat)
    diffusion = modes[int(np.argmin(np.abs(modes.real)))]
    return float(-diffusion.imag / probe_qhat**2)
