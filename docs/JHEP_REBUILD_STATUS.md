# JHEP reconstruction status

The manuscript now follows Route A: **critical relaxation and spinodal
instabilities in Bayesian-calibrated holographic QCD**. The former homogeneous
Route-B draft is retained only through its validated nonhydrodynamic and 2018
cross-model benchmarks; it no longer defines the paper's scope.

## Scientific completion

- The exact public MAP realization is checksum tracked. An independent
  fold/null-direction solution locates the cusp at
  `(T_c, mu_B^c)=(103.7553, 593.3227) MeV`; coexistence and both spinodals are
  continued from 594 to 850 MeV with grid-coarsening errors stored explicitly.
- The complete EMD equations were linearized in ingoing EF coordinates.
  Generated helicity-one and helicity-zero artifacts pass ten symbolic tests
  covering Bianchi/Maxwell identities, residual gauge modes, parity, neutral
  and zero-momentum limits, source/vev powers, and degree counting.
- Physical gauge-invariant horizon bases are evolved by DOP853 source-matrix
  shooting. A distinct multidomain Chebyshev first-order integrator provides
  the independent radial method. Unused primitive equations provide
  constraint diagnostics.
- Six neutral/charged shear, sound, and diffusion benchmarks pass separate
  resolution, UV-domain, horizon-start, source, gauge, constraint, and method
  comparisons. Their small-momentum coefficients agree with independently
  evaluated charged hydrodynamics.
- Along the MAP critical approach, `D_B` vanishes while `chi_B` diverges. Five
  nested windows give `z=4.08+/-0.05`; the direct critical QNM dispersion gives
  `z=3.99(2)`.
- Three full longitudinal scans inside the `mu_B=650 MeV` spinodal branch give
  closed unstable bands. At the midpoint,
  `q*=0.2826`, `Gamma*=1.414e-3`, and `q_edge=0.4175`.
- Twenty-five deterministic weighted posterior medoids represent all 1589
  successful HDF5 draws and 100% of their posterior weight. All 25 local
  cusps/spinodal pairs reconstruct successfully; the 95% posterior interval is
  `4.0625 < z < 4.0848`, before the common fit-window systematic, and every
  medoid has negative midpoint spinodal diffusion.

## Reproducible artifacts

The claim-bearing products are:

- `results/python/coupled_qnm_validation.csv` and summary JSON;
- `finite_k_hydrodynamic_dispersion.csv`;
- `cep_critical_scaling.csv` and `cep_q4_dispersion.csv`;
- `spinodal_dispersion.csv` and `finite_k_physics_summary.json`;
- `posterior_selection.csv`, `posterior_uq_samples.json`, and summary JSON;
- `paper/figures/phase_diagram.pdf`,
  `hydrodynamic_critical_dispersion.pdf`, `cep_spinodal_dynamics.pdf`, and
  `posterior_uq.pdf`.

The finite-momentum driver is checkpointed after every root. The posterior
driver includes an adaptive cusp fallback, so no failed medoid is silently
discarded or manually retuned.

## Manuscript and submission status

`paper/main.pdf` is a 16-page JHEP-style working manuscript with a first-page
abstract, controlled keywords, four vector figures, an explicit numerical
error appendix, reproducibility commands, and AI-assisted-technology
disclosure. It has no undefined references/citations, duplicate labels,
overfull boxes, or invented author information.

Scientific and numerical Route-A gates are closed. Literal submission remains
blocked only on author-controlled names, affiliations, email, ORCID, funding,
arXiv, authorship approval, and persistent-release metadata. The guarded
archive builder deliberately refuses to package until those fields exist.
