# Validation protocol

Accepted backgrounds target relative Einstein-constraint residual below `1e-9` and relative Gauss drift below `1e-10` away from singular/critical points. Accepted decoupled QNMs require at least three radial resolutions, stable domain interfaces, a normalized pencil residual, boundary residual, and a shooting cross-check. Coupled modes additionally require unused Einstein/Maxwell constraints below `1e-7`, gauge independence, and correct `k -> 0` representation content.

No returned eigenvalue is physical merely because a numerical eigensolver produced it. Result records retain convergence, residual, constraint, conditioning, solver, compiler, and provenance fields.

