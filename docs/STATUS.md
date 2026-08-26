# Status and validation gates

Last updated: 2026-08-25.

## Implemented and validated

- C++23/CMake reference kernels, exact checksum-verified PHA HDF5 selection,
  fourth-order horizon recurrence, adaptive background integration, and
  decoupled source-factored spectral operators.
- Independent MAP cusp reconstruction, a `141 x 121` phase surface,
  fixed-chemical-potential continuation through both folds, equal-pressure
  coexistence, stability labels, and grid-coarsening tests.
- Complete symbolic EMD linearization in ingoing EF coordinates with generated
  helicity-one/helicity-zero JSON and C++ kernels, residual gauge modes,
  constraints, source/vev maps, and ten passing identity/limit tests.
- Gauge-invariant coupled QNM solver with second-order horizon recurrence,
  DOP853 source-matrix shooting, independent multidomain Chebyshev integration,
  conditioning, source singular values, gauge defects, and unused primitive
  Einstein/Maxwell constraint diagnostics.
- Six neutral/charged shear, sound, and baryon/heat-diffusion benchmarks with
  independent variations of resolution, domain, partition, UV extraction, and
  horizon start.
- Charged-hydrodynamic coefficients from equilibrium derivatives, the
  Eling--Oz bulk-viscosity relation, and incoherent horizon conductivity.
- Checkpointed stable, critical, and spinodal finite-momentum scans; critical
  fit-window systematics; a direct `q^4` critical QNM test; and closed spinodal
  unstable bands at three distances between the folds.
- Posterior uncertainty propagation through 25 deterministic weighted medoids
  representing all 1589 successful draws and 100% of posterior weight, with
  25/25 successful independent cusp/spinodal reconstructions.
- Homogeneous three-sector cusp modes, physical trajectories, and the
  independently reproduced 2018 EMD relaxation-time benchmark retained as
  nonhydrodynamic and cross-model validation.
- A 25-page Route-A JHEP manuscript with four central physics figures, seven
  restored validation/benchmark figures, and machine-readable claim evidence.

## Current numerical acceptance summary

- Coupled shooting/collocation distance: at most `3.63e-12`.
- Continuum displacement across the benchmark suite: at most `3.64e-11`.
- Gauge defect: at most `4.05e-16`.
- Source determinant: at most `1.60e-12`.
- Primitive constraint: at most `2.41e-5` in the doubled-horizon-start stress
  test, below the stored `3e-5` threshold.
- Stable hydrodynamic/QNM coefficient difference at `q=0.05`: at most `0.166%`.
- Direct critical-QNM source determinant: at most `1.67e-11`; primitive
  constraint: at most `9.09e-7`.
- Spinodal midpoint shooting and fine-grid cross-checks agree within
  `1.49e-10` and displayed precision, respectively.

## Scope boundaries

- The independently located cusp differs from the public HDF5 metadata by
  about `0.143 MeV` in temperature and `9.152 MeV` in chemical potential. This
  provenance discrepancy is reported and is not folded into numerical error.
- The weighted posterior study propagates equilibrium and horizon-transport
  dynamics (`D_B`, `z`, and the sign of spinodal diffusion). Full coupled-QNM
  band scales are claim-bearing MAP results and are explicitly labeled as
  such; no posterior full-QNM credible band is invented.
- Linear QNMs characterize relaxation of small perturbations, not nonlinear
  far-from-equilibrium thermalization.
- The results are predictions of a classical large-`N_c`, strongly coupled,
  bottom-up holographic model calibrated to QCD thermodynamics, not QCD itself.

## External submission fields

Scientific, numerical, and editorial Route-A gates are closed. Authors must
still supply names, affiliations, email, ORCID, funding/acknowledgments, arXiv
identifier, authorship approval, and persistent-release metadata. The guarded
submission packager fails until those fields are present.
