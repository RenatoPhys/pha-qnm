# Numerical method

Backgrounds are initialized from a regular horizon series and integrated in
horizon-data coordinates `(phi0, Phi1)`. The high-accuracy formulation evolves
the conserved Maxwell charge `Q=exp(2A) f Phi'`; UV observables are accepted
only after scalar and Ricci asymptotics and the Einstein/Gauss diagnostics
pass.

The reference critical-point route follows neighboring constant-`phi0` lines.
An independent route solves the fold equation
`det d(log T,log mu)/d(phi0,charge_fraction)=0` together with the vanishing
derivative of that determinant along the Jacobian null direction.
Fixed-chemical-potential lines are continued through both folds. Pressure is
integrated from `dp=s dT`, with coexistence at equal pressure. The production
`141 x 121` surface and two coarsened grids set the quoted phase-line errors.

QNMs use ingoing Eddington--Finkelstein coordinates and the convention
`exp(-i omega v+i k z)`. Residual radial-gauge orbits are quotiented at the
horizon and the evolved variables are physical gauge invariants. A
second-order regular-horizon recurrence initializes a basis of ingoing
solutions. The primary solver integrates that basis with complex DOP853 and
extracts the UV source matrix; roots minimize its smallest singular value and
are checked through the normalized determinant.

The independent radial method uses multidomain Chebyshev--Lobatto
differentiation of the equivalent first-order system, enforces interface
continuity, and reads the same physical UV sources without sharing the DOP853
integration. Resolution, radial domain, domain split, UV extraction, and
horizon start are varied independently. Unused primitive Einstein/Maxwell
equations are evaluated on reconstructed fields as constraint diagnostics.

Charged hydrodynamics is computed separately from local EOS derivatives,
`eta/s=1/(4 pi)`, the Eling--Oz bulk-viscosity formula, and the incoherent
horizon conductivity. These coefficients predict the small-momentum shear,
sound, and diffusion poles before any QNM fit. Critical and spinodal drivers
are checkpointed after each accepted root. Posterior medoids reconstruct their
own cusp and spinodal pair; failed adaptive steps are retried and never dropped
silently.
