# PHA QNM

Native C++23/Python research code and a JHEP-format manuscript for critical relaxation and spinodal instabilities in the Bayesian-calibrated polynomial--hyperbolic Ansatz (PHA) Einstein--Maxwell--dilaton model.

The repository contains the checksum-tracked posterior and MAP realization,
continued equilibrium branches, complete finite-momentum helicity-one and
helicity-zero systems, coupled source-matrix shooting, an independent
multidomain Chebyshev integrator, charged-hydrodynamic cross-checks, critical
and spinodal scans, and posterior propagation through 25 weighted medoids.
Every claim-bearing result is stored in a machine-readable CSV or JSON file;
the exact acceptance boundary is recorded in `docs/STATUS.md` and
`docs/claim_evidence_audit.md`.

## Quick start

On a machine with CMake 3.28+ and a C++23 compiler:

```text
cmake --preset release-portable
cmake --build --preset release-portable --parallel
ctest --preset release-portable --output-on-failure
build/release-portable/pha-qnm validate all
```

The production data configuration requires HDF5, CURL, and OpenSSL:

```text
cmake --preset release-data -DCMAKE_TOOLCHAIN_FILE=%VCPKG_ROOT%/scripts/buildsystems/vcpkg.cmake
cmake --build --preset release-data --parallel
build/release-data/pha-qnm data verify data/raw/Bayesian_polyhyper_muses.hdf5
build/release-data/pha-qnm model extract-map data/raw/Bayesian_polyhyper_muses.hdf5 configs/pha_map.generated.yaml
```

The dense QZ path additionally requires LAPACK/MKL (`release-mkl`); the distributed path uses PETSc/SLEPc (`hpc-slepc`).

Run the numerical analysis and regenerate the figures with:

```text
python -m pip install -r analysis/requirements.txt
python analysis/reproduce_reference_cep.py --lines 800 --neighbors-only --output results/python/reference_cep_neighbors_N800.json
python analysis/locate_reference_cusp.py --phi-step 0.01 --charge-step 0.001 --uv-tolerance 1e-5 --output results/python/reference_cusp_uv1e5_h10.json
python analysis/build_phase_diagram.py --n-phi 141 --n-charge 121 --phi-min 2.2 --phi-max 6.5 --charge-min 0.02 --charge-max 0.62 --mu-min 594 --mu-max 850 --mu-count 65
python analysis/validate_phase_diagram.py
python analysis/plot_cep_reproduction.py
python analysis/run_numerics.py
python derivations/generate_linearized_emd.py
python derivations/validate_linearized_emd.py
python analysis/validate_decoupled_qnms.py
python analysis/build_homogeneous_trajectories.py
python analysis/compare_2018_emd.py
python analysis/plot_qnm_validation.py
python analysis/plot_homogeneous_trajectories.py
python analysis/validate_coupled_qnms.py
python analysis/run_finite_k_physics.py --stage all
python analysis/run_posterior_uq.py --samples 25 --workers 4
python analysis/plot_posterior_uq.py
```

The production phase-surface run evaluates 17,061 backgrounds and supports
parallel workers through `--workers`.  After the surface exists, pass
`--reuse-surface` to rebuild the fixed-chemical-potential branches and plots
without reintegrating it.

The C++ code remains the independent reference implementation; the numerical
paper results are written to `results/python/` and `paper/figures/`.

The repository includes the generated manuscript at [`paper/main.pdf`](paper/main.pdf).
Appendix B reproduces the neutral entropy and baryon susceptibility with three
centered finite-difference steps. Its 90 machine-readable evaluations and
diagnostics are stored in
[`results/python/thermodynamic_validation.csv`](results/python/thermodynamic_validation.csv).

Compile the JHEP manuscript on Windows with the pinned, checksum-verified
Tectonic bootstrapper:

```text
powershell -ExecutionPolicy Bypass -File paper/build_pdf.ps1
```

## Scientific scope

The Fourier convention is `exp(-i omega v + i k z)`; stable modes therefore
have negative imaginary frequency.  The Route-A manuscript identifies shear,
sound, and baryon/heat diffusion on neutral and charged stable backgrounds,
extracts critical slowing at the independently located MAP cusp, resolves the
closed longitudinal spinodal band at three locations between its folds, and
propagates the thermodynamic/horizon-transport observables through a weighted
posterior ensemble.  Full finite-momentum QNM band scales are reported for the
MAP realization; the posterior study establishes robustness of the critical
exponent and the sign of spinodal diffusion without presenting uncomputed
posterior full-QNM intervals.

As an external model benchmark, the same numerical pipeline reproduces all
nine approximate relaxation times quoted for the EMD realization of
Phys. Rev. D 98 (2018) 034028. The versioned comparison is stored in
`results/python/rougemont2018_benchmark.{csv,json}`.

See `docs/STATUS.md` for the exact implemented/validated boundary and `docs/DECISIONS.md` for recorded choices.

## Citation and license

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). The source code
is distributed under the BSD 3-Clause License; see [`LICENSE`](LICENSE).
