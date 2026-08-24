# Status and validation gates

Last updated: 2026-08-23.

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
- At the computed background nearest the posterior CEP metadata, the constraint residual is `3.39e-11`; six homogeneous tensor/vector candidates pass the four-grid pseudospectral gate.
- Two tensor candidates were continued over 16 momentum points from `q=0` to `q=3`; all remain in the stable half plane.

## Not yet scientifically complete

- Full thermodynamic surface/branch continuation and UV extraction still require validation against arXiv:2603.20482. The nearest independently computed background differs from the CEP metadata by `-0.620 MeV` in temperature and `-1.156 MeV` in chemical potential; the local inverse-map condition number is about `1.79e3`.
- Dense LAPACK QZ and SLEPc backends require dependency-enabled builds and numerical smoke tests. The pinned vcpkg `lapack-reference` port failed its Windows Fortran/C compatibility check on this host; the MKL/HPC presets remain the production routes.
- The homogeneous scalar equation still requires an independent symbolic derivation.
- Finite-k helicity-1 and helicity-0 systems, Bianchi/Maxwell identity checks, and source analysis are not implemented.
- The reported tensor/vector values are pseudospectral candidates, not publication-grade accepted poles: independent shooting, split-grid source factoring, and overlap-based tracking remain outstanding.
- No scalar, coupled helicity-1/helicity-0 frequency, critical exponent, diffusion constant, pole collision, or spinodal growth rate is claimed.

## Next commands

```text
cmake --preset release-portable
cmake --build --preset release-portable --parallel
ctest --preset release-portable --output-on-failure
build/release-portable/pha-qnm validate all
```

For the native data gate, configure `release-data` with a vcpkg toolchain, then run `data verify` and `model extract-map`. Configure `release-mkl` or `hpc-slepc` before any QNM scan.
- Two tensor candidates were continued over 16 momentum points from `q=0` to `q=3`; all remain in the stable half plane.
- The expanded JHEP manuscript was compiled with Tectonic 0.17.0 into a 13-page PDF; all 13 pages were rendered at 130 dpi with Poppler and visually inspected.
