# Decisions

## 2026-08-23

- The user's explicit venue request overrides the briefing's provisional RevTeX default: the manuscript uses the JHEP class/style layout.
- The project name is `pha-qnm` and its intended destination is `C:\\Users\\User\\OneDrive\\Documentos\\rnt\\Physics\\pha-qnm`.
- All production computations are compiled C++23. A temporary, untracked read-only HDF5 inspection was used only to cross-check the exact output expected from the native loader; it is not part of any result path.
- The exact maximum-likelihood realization is `posterior_samples/sample74` from Zenodo v1, not the rounded parameter table and not the separate default currently shown in MUSES documentation.
- No novelty or numerical physics claim enters the manuscript before its validation gate is recorded as passing.
- Ingoing Eddington--Finkelstein coordinates and `exp(-i omega v + i k z)` are canonical; stable poles have `Im(omega)<0`.
- Dense generalized eigenproblems retain the pencil `(M0 + omega M1)v=0`; the code must never construct `M1^{-1}M0`.
- The first portable build has no unsafe floating-point flags. Production HDF5/network/LAPACK capabilities are explicit CMake options so a dependency-free reference core remains testable.

