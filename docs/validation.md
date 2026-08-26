# Validation protocol

Accepted backgrounds target relative Einstein-constraint residual below
`1e-9` and Gauss drift below `1e-10`, with stored pointwise diagnostics near
singular/critical regions.

An accepted coupled pole requires all of the following:

1. a small normalized source determinant and source singular value;
2. residual-gauge defect below the stored threshold;
3. convergence of unused primitive Einstein/Maxwell constraints;
4. separate stability under radial resolution, UV domain, domain partition,
   UV extraction, and horizon-start changes;
5. agreement between DOP853 source-matrix shooting and the independent
   multidomain Chebyshev integrator;
6. continuity with the expected neutral, zero-momentum, or hydrodynamic mode;
7. a finite, recorded horizon-recurrence and principal-matrix condition.

The automated six-case benchmark thresholds are versioned in
`coupled_qnm_validation_summary.json`; the deliberately doubled horizon start
uses a primitive-constraint ceiling of `3e-5`. Critical and spinodal headline
roots receive additional fine-grid and direct-shooting checks.

Hydrodynamic coefficients are derived independently and compared
coefficient-by-coefficient. Critical fits use nested windows and quote their
half-range as a systematic. A spinodal claim requires a full-QNM band that
returns to the stable half-plane at finite momentum. Posterior claims require
explicit represented weights and successful reconstruction counts.

No eigenvalue is accepted merely because it solves a discrete pencil or makes
an optimizer report success.
