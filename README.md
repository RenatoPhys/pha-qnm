# PHA QNM

Native C++23 research code and a JHEP-format working manuscript for quasinormal modes of the maximum-likelihood realization of the Bayesian-calibrated polynomial--hyperbolic Ansatz (PHA) Einstein--Maxwell--dilaton model.

The repository contains an exact, checksum-tracked MAP realization, analytic model kernels, fourth-order horizon data, a reference background integrator, source-factored Chebyshev operators, independent complex shooting, and controlled homogeneous QNM trajectories. Coupled helicity-1/helicity-0 kernels are generated and symbolically tested but remain outside the homogeneous paper's claims; the precise boundary is recorded in `docs/STATUS.md`.

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

The Fourier convention is `exp(-i omega v + i k z)`; stable modes therefore have negative imaginary frequency. Four homogeneous poles at the independently located MAP cusp pass source-factored continuum studies and independent shooting. The leading quintuplet, triplet, and singlet poles are continued on `mu_B=0` and `mu_B/T=2`. The paper does not claim diffusion, sound, spinodal growth, or posterior uncertainty.

As an external model benchmark, the same numerical pipeline reproduces all
nine approximate relaxation times quoted for the EMD realization of
Phys. Rev. D 98 (2018) 034028. The versioned comparison is stored in
`results/python/rougemont2018_benchmark.{csv,json}`.

See `docs/STATUS.md` for the exact implemented/validated boundary and `docs/DECISIONS.md` for recorded choices.

## Citation and license

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). The source code
is distributed under the BSD 3-Clause License; see [`LICENSE`](LICENSE).
