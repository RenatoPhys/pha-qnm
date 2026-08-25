# Thermodynamic reference audit

Last updated: 2026-08-24.

## Provenance

The public reference implementation was inspected from
`https://gitlab.com/nsf-muses/module-holographic-eos/muses-numrelholo.git`.
The contemporaneous source state used for the comparison is commit
`52b093230a8b5f36a83fc37c6fb5e8b2043c82c4` (2023-08-24); release tag `v1.3.0`
was also retained locally for comparison. The external worktrees live under
the ignored `external/` directory and are not repository inputs.

The raw posterior file is `data/raw/Bayesian_polyhyper_muses.hdf5`. Its maximum
log-likelihood PHA realization is `posterior_samples/sample74`, with
`logL = 31.075137487467885`. All model parameters used in the audit are the
exact binary64 values from that group.

## Reconstructed reference route

The 2023 algorithm uses:

- a second-order regular horizon series;
- the relaxation equation `C' = 8[phi exp(nu A) - C]` for the UV scalar
  coefficient;
- joint scalar-coefficient and Ricci convergence conditions in the UV;
- the conserved Gauss charge for `Phi2`;
- constant-`phi0` trajectories over `0.5 <= phi0 <= 7`, initialized with a
  charge-fraction step of `1e-4`;
- critical-point selection from crossings of neighboring trajectories.

The reproduced neighboring-crossing estimates are:

| number of trajectories | `T_c` [MeV] | `mu_B^c` [MeV] |
|---:|---:|---:|
| 100 | 103.665610 | 594.761441 |
| 200 | 103.746458 | 593.605539 |
| 400 | 103.749470 | 593.562252 |
| 800 | 103.761142 | 593.384499 |

## Independent cusp and phase construction

The independent locator solves the fold condition on the logarithmic
thermodynamic map together with the vanishing derivative of the fold
determinant along the Jacobian null direction. With finite-difference steps
`dphi0 = 0.01`, `d(charge_fraction) = 0.001`, and UV tolerance `1e-5`, it gives:

- `phi0 = 3.528280966207555`;
- `charge_fraction = 0.34637356467181973`;
- `T_c = 103.75528111978664 MeV`;
- `mu_B^c = 593.3226560539789 MeV`.

The result is `0.0059 MeV` and `0.0618 MeV` from the 800-trajectory estimate.
A production `141 x 121` horizon-data surface was then interpolated into 65
constant-chemical-potential continuations over `594--850 MeV`. Both folds and
an equal-pressure coexistence point were found on every line. Grid coarsening
by factors two and four changes the spinodal temperatures by at most
`0.0178 MeV` and the coexistence temperatures by at most `0.0085 MeV`.

## Metadata finding

The HDF5 group stores:

- `T_c = 103.89777513571676 MeV`;
- `mu_B^c = 602.4749076542938 MeV`;
- numerical error estimates `0.0011068 MeV` and `0.0028784 MeV`.

Relative to the independent cusp, these coordinates are higher by
`0.142494 MeV` and `9.152252 MeV`. Because the two locally implemented routes
agree and the exact sample parameters were used, the file metadata are not
used to select the QNM background. The discrepancy is preserved explicitly;
no causal explanation is claimed without further provenance information.

## Machine-readable evidence

- `results/python/reference_cep_neighbors_N800.json`
- `results/python/reference_cusp_uv1e5_h10.json`
- `results/python/phase_surface.npz`
- `results/python/phase_lines.csv`
- `results/python/phase_line_grid_convergence.json`
- `results/python/phase_diagram_summary.json`
