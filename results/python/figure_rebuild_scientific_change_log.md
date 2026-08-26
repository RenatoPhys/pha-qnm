# Figure-rebuild scientific change log

Date: 2026-08-26

## Existing validated results

No tracked production result was changed. The pre/post Git blob hashes are
identical for:

| product | Git blob hash |
|---|---|
| `finite_k_physics_summary.json` | `9730f1e10883b43b66d1991c2ce296818a04aa65` |
| `posterior_uq_summary.json` | `3cc713e65a71d41d3b585b5958ed1bbb1ccef33d` |
| `coupled_qnm_validation.csv` | `0457cf534021bef8b406573b53f79ebe5640494f` |
| `thermodynamic_validation.csv` | `2d877e2ff20b09dd0f05330f017d32fa9cf8e436` |
| `phase_diagram_summary.json` | `6d9618d7dc4b6ae68c8e3de2464494c2b78c1719` |
| `rougemont2018_benchmark.json` | `2860577b781da1346edfc79a0af64529f9774400` |

The MAP critical exponent, critical QNM exponent, spinodal maximum, band edge,
posterior intervals, phase lines, and validation thresholds are therefore
numerically invariant under the editorial rebuild.

## New scientific audit

`longitudinal_mode_identity_audit.csv/json` adds a dedicated, non-destructive
mode-identity study along `mu_B=650 MeV`:

- nine ordered backgrounds cross both folds and include stable points on both
  sides;
- the complete charged-hydrodynamic eigensystem is stored at eight small
  momenta, with eigenvectors normalized in relative thermodynamic variables;
- accepted production roots are reused and missing QNM roots are continued in
  a dedicated checkpoint;
- `c_s^2` stays positive and the sound pair stays stable across the scan;
- signed `D_B` alone changes sign between the folds;
- QNM leading powers are 1.990 (near hot fold), 1.994 (midpoint), and 1.984
  (near cold fold);
- `q^2+q^4` is preferred over models with a leading linear term at all three
  interior backgrounds;
- the first longitudinal nonhydrodynamic pole remains gapped at the midpoint.

Decision: `instability_carrier = diffusion`. The main text may retain
`baryon/heat-diffusion branch`, provided it states the positive sound speed and
contrasts this charged instability with sound-driven neutral spinodals.

The largest finite audit source determinant is `1.08e-9` and the largest
primitive constraint is `3.26e-6`, both inside the established acceptance
envelopes.

## New 2018-model curves

`legacy_2018_qnm_curves.csv/json` adds 195 source-factored collocation points on
65 independently solved backgrounds at `mu_B = 0, 300, 500, 600, 700 MeV` in
the quintuplet, triplet, and singlet sectors. Piecewise-linear segments are
declared guides to the eye. The maximum pencil residual is `1.12e-15`; the
existing nine-point shooting benchmark remains unchanged and supplies the
independent radial cross-check. No digitized published curve or universality
claim is introduced.

