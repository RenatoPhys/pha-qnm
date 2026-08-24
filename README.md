# PHA QNM

Native C++23 research code and a JHEP-format manuscript for quasinormal modes of the Bayesian polynomial--hyperbolic Ansatz (PHA) Einstein--Maxwell--dilaton model.

The repository is deliberately honest about validation. It contains an exact, checksum-tracked MAP realization, analytic model kernels, fourth-order horizon data, a reference background integrator, Chebyshev operators, and decoupled QNM pencil assembly. A user-authorized Python layer now runs the background scans, convergence studies, QNM candidate extraction, and manuscript plots. Coupled helicity-1/helicity-0 kernels and publication-grade acceptance remain gated in `docs/STATUS.md`.

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
python analysis/run_numerics.py
```

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

The Fourier convention is `exp(-i omega v + i k z)`; stable modes therefore have negative imaginary frequency. The current tensor and homogeneous-vector values pass a four-grid pseudospectral candidate gate, but remain pending an independent shooting check. Bulk viscosity is not a novelty target and is used only as optional external hydrodynamic input.

See `docs/STATUS.md` for the exact implemented/validated boundary and `docs/DECISIONS.md` for recorded choices.

## Citation and license

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). The source code
is distributed under the BSD 3-Clause License; see [`LICENSE`](LICENSE).
