# Symbolic derivation gate

Finite-momentum helicity-1 and helicity-0 kernels may not be handwritten into the production solver. This directory will contain the auditable tensor expansion, unsimplified equations, background substitutions, Bianchi/Maxwell identity reductions, optimized generated C++ kernels, and expression hashes.

No file under `generated/` is authoritative until its debug form agrees pointwise with the unsimplified equations and its `k -> 0`, neutral, and conformal tests pass.

