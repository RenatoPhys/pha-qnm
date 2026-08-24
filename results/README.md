# Result schema

Authoritative spectra are versioned HDF5, written by C++ and chunked by scan. Every accepted record includes:

```text
background_id, phi0, Phi1, T_MeV, muB_MeV, branch_label,
stability_label, channel, helicity, q_dimensionless, k_numerical,
mode_id, omega_real_dimensionless, omega_imag_dimensionless,
resolution, residual, constraint_residual, condition_number, converged,
tracking_parent, solver_backend, blas_lapack_vendor, mpi_size,
omp_threads, compiler_id, compiler_version, code_commit, config_hash
```

CSV exports are inspection-only. No result currently in this directory is publication-authoritative.

