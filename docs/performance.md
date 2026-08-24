# Performance

No optimization claim is made before profiling. The reference path uses contiguous column-major matrices and preallocated state vectors. Planned measurements separately cover background RHS evaluation, horizon recurrence, UV extraction, operator assembly, dense QZ, shift-invert, and HDF5 I/O. Production manifests record compiler flags, BLAS vendor, CPU, MPI ranks, OpenMP threads, and floating-point mode. Unsafe `-ffast-math`/`-Ofast` flags are forbidden.

