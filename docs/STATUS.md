# Status and validation gates

Last updated: 2026-08-25.

## Implemented

- C++23/CMake repository with portable, sanitizer, native, production-data, MKL, and HPC presets.
- Zenodo data manifest and exact MD5-verified PHA HDF5 file (raw file ignored by Git).
- Generic native HDF5 maximum-log-likelihood selector and YAML freezer (compiled when HDF5 is enabled).
- Exact MAP model kernels `V`, `V_phi`, `V_phiphi`, `f`, `f_phi`, and `f_phiphi` with stable hyperbolic evaluation.
- Automatic fourth-order horizon Taylor recurrence in compiled C++.
- Reference background equations, constraint/Gauss diagnostics, and adaptive integration.
- Chebyshev--Lobatto matrices and decoupled tensor/vector pencil assembly without forming a matrix inverse.
- Charged-hydrodynamic cubic polynomial helper.
- JHEP-format manuscript with numerical tables/figures and verified bibliography metadata.
- User-authorized Python 3.12 analysis layer with a conserved-charge DOP853 background solver, fixed-grid UV extraction, SciPy generalized QZ solve, convergence scans, and Matplotlib output.
- First JHEP reconstruction pass: MAP-only framing, honest horizon-trajectory labeling, condensed critical-scope discussion, expanded literature, and serif figure style.
- Reference-compatible constant-`phi0` critical-point reproduction, an independent fold/null-direction cusp locator, and a documented HDF5 metadata discrepancy.
- A `141 x 121` thermodynamic surface with 65 fixed-`mu_B` continuations, both spinodals, equal-pressure coexistence, and two-level grid-coarsening checks.
- Complete symbolic EMD linearization in ingoing EF coordinates, generated C++23 helicity-one/helicity-zero kernels, residual gauge modes, constraints, and source/vev maps.
- Homogeneous singlet master-field equivalence to arXiv:1804.00189.
- Analytically source-factored homogeneous operators, separate resolution/domain convergence, left/right conditioning, and eigenfunction overlaps.
- Independent complex shooting for all claim-bearing cusp modes, including horizon-start and UV-extraction variations.
- Physical homogeneous trajectories on `mu_B=0` and `mu_B/T=2`, with full thermodynamic state, branch/stability labels, parent indices, and source residuals.
- Independent reconstruction of the 2018 EMD model and a three-sector comparison
  at two common neutral states and the two model-specific critical endpoints.

## Validated builds

MSVC 19.44.35222, debug, release-portable, and release-data builds, 2026-08-23:

- Potential analytic second derivatives agree with complex-step checks.
- Fourth-order horizon recurrence reproduces all supplied first/second coefficients and remains finite across the tested regular state.
- Representative background (`phi0=1`, `Phi1=0.05`) reached normalized Einstein-constraint residual `8.31151802019e-13` and relative Gauss drift `1.64001264346e-12` over the CLI validation interval.
- Chebyshev differentiation error for the CLI polynomial check was `2.22044604925e-15`.
- The tensor pencil has the exact `k -> 0` reduction and a source-free UV tau row in the unit test.
- All six test groups passed under CTest in every enabled build.
- The dependency-enabled `release-data` build verified the 60,044,400-byte HDF5 checksum and selected `posterior_samples/sample74` natively in C++.
- All 16 floating-point MAP/critical-point fields regenerated from HDF5 are bit-identical to the frozen configuration after parsing as IEEE-754 binary64 values.
- A 90-background Python scan covers 50.34--1003.70 MeV; Gauss drift is below `4.45e-16` and the maximum normalized constraint residual is `3.50e-9` over the full exploratory grid.
- At the independently located cusp, the reduced-background constraint residual is `1.88e-11`; four homogeneous poles pass both factored collocation and independent shooting. Three former discrete-pencil candidates are rejected.
- The two physical trajectories contain 108 pointwise shooting roots. Source residuals stay below `2.2e-8`, parent overlaps above `0.9985`, and shooting/spectral differences below `0.012` in `omega_hat`.
- All nine approximate relaxation times quoted in Phys. Rev. D 98 (2018)
  034028 are reproduced within `0.0111 fm/c`; across the 18-row comparison,
  source residuals stay below `7.42e-9` and spectral/shooting distances below
  `0.00731` in `omega_hat`.
- The reconstruction build passes the release-portable C++ test target and `pha-qnm validate all`; regenerated scientific values are identical to the frozen baseline for the compared background, CEP diagnostic, QNM, dispersion, trajectory-count, and thermodynamic fields.
- The revised 19-page manuscript compiles with Tectonic 0.17.0 without undefined references, undefined citations, duplicate labels, or overfull boxes; all pages were rendered with Poppler and visually inspected on 2026-08-25.
- The local `jheppub.sty` and `JHEP.bst` are byte-identical to the official
  JHEP downloads checked on 2026-08-25. Controlled JHEP keywords and the
  required AI-assisted-technology disclosure are present, and the guarded
  submission packager requires `main.bbl` and complete author/arXiv metadata.

## Scope boundaries and external submission metadata

- The reference-compatible and independent cusp routes agree near `(103.755, 593.323) MeV`, while the posterior HDF5 metadata remain higher by about `0.143 MeV` in temperature and `9.152 MeV` in chemical potential. This is retained as a provenance discrepancy, not folded into numerical uncertainty.
- Dense LAPACK QZ and SLEPc backends require dependency-enabled builds and numerical smoke tests. The pinned vcpkg `lapack-reference` port failed its Windows Fortran/C compatibility check on this host; the MKL/HPC presets remain the production routes.
- The coupled finite-`k` equations are generated and pass symbolic identity/limit tests. A primitive-field prototype reproduces the helicity-one shear limit, but helicity-zero constraint residuals and continuum convergence fail the acceptance gate; no coupled spectrum is claimed.
- No coupled helicity-1/helicity-0 frequency, critical exponent, diffusion constant, pole collision, or spinodal growth rate is claimed. The manuscript follows the homogeneous Route B.
- Posterior uncertainty is not propagated; every title/abstract/result claim is explicitly restricted to the MAP realization.
- Author, affiliation, email, ORCID, funding, acknowledgment, and arXiv metadata must be supplied by the authors before submission.
- The public GitHub repository still contains the pre-reconstruction commit;
  this branch and its machine-readable outputs must be published there no
  later than manuscript submission, then archived with a persistent DOI for
  the JHEP proofreading stage.

## Next commands

```text
cmake --preset release-portable
cmake --build --preset release-portable --parallel
ctest --preset release-portable --output-on-failure
build/release-portable/pha-qnm validate all
```

For the native data gate, configure `release-data` with a vcpkg toolchain, then run `data verify` and `model extract-map`. Configure `release-mkl` or `hpc-slepc` before any QNM scan.
