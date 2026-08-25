# JHEP reconstruction status

The submission manuscript follows Route B: a complete homogeneous-response
paper for the maximum-likelihood PHA realization. Route A remains a separate
future project because the primitive helicity-zero finite-momentum pencil did
not pass constraint and continuum-convergence tests. No failed coupled result
is used in the manuscript.

## Completed in the first reconstruction pass

- Froze commit `c3c71a7` as annotated tag `pre-jhep-rebuild` and recorded the
  baseline blobs in `docs/baseline_manifest.txt`.
- Created branch `codex/jhep-rebuild`.
- Reframed the title and abstract as a maximum-likelihood, lattice-calibrated
  study rather than a Bayesian posterior prediction.
- Removed the main section that purported to discuss critical/spinodal dynamics
  without the coupled helicity-zero calculation.
- Recast the fixed-horizon scan as a horizon-data trajectory survey and stated
  explicitly that it is not temperature dependence at fixed thermodynamic
  control variables.
- Moved dataset checks and software metadata to the reproducibility appendix.
- Reduced parameter/table precision in the paper while preserving exact values
  in YAML/JSON.
- Expanded the literature basis for dynamic universality, critical QNMs,
  spinodals, charged-fluid pole collisions, and pseudospectral sensitivity.
- Changed plot typography to a serif/STIX system with restrained grid use.
- Rebuilt the C++ release-portable target, passed its test suite and CLI
  validation, regenerated the Python outputs without changing the compared
  scientific values, and visually inspected all 15 pages of the final PDF.

## Completed thermodynamic gate (P1)

- Reconstructed the public 2023 reference-compatible critical-point route and
  followed neighboring constant-`phi0` crossings through 800 trajectories.
- Independently solved the fold/null-direction cusp equations, obtaining
  `(T_c, mu_B^c) = (103.7553, 593.3227) MeV` at
  `(phi0, charge_fraction) = (3.528281, 0.346374)`.
- Built a `141 x 121` horizon-data surface and continued 65 fixed-chemical-
  potential lines through their stable, metastable, and unstable branches.
- Located both spinodals and the equal-pressure coexistence point on every
  line from `mu_B = 594` to `850 MeV`.
- Coarsened-grid comparisons bound spinodal-temperature changes by `0.018 MeV`
  and coexistence-temperature changes by `0.0085 MeV`.
- Preserved the provenance discrepancy: the HDF5 metadata chemical potential
  is `9.152 MeV` above the independently located cusp. No parameters were
  rounded or retuned to hide it.

## Completed equation gate (P2)

- Independently linearized the complete EMD equations in ingoing EF
  coordinates with exact first-order dual expressions.
- Generated the helicity-one and helicity-zero radial coefficient maps and
  C++23 kernels with SHA-256 provenance.
- Retained three evolution equations plus two diagnostics in helicity one and
  seven frequency-linear characteristic equations plus seven diagnostics in
  helicity zero.
- Verified Bianchi and Maxwell identities, residual gauge modes, transverse
  parity, neutral and `k -> 0` limits, UV source/vev powers, and the physical
  degree-of-freedom count.
- Proved that the homogeneous singlet pencil is the published gauge-invariant
  master equation after `psi=S/(phi'/A')`, and obtained its first
  four-grid-stable candidate at the independent cusp.

## Completed homogeneous solver gate (P3) and physical trajectories

- Assembled the UV-factored quintuplet, triplet, and singlet operators
  analytically and varied radial resolution and UV domain independently.
- Computed left/right condition estimates and eigenfunction overlaps.
- Implemented complex shooting with independent horizon integration and UV
  source extraction; varied both the horizon start and UV extraction radius.
- Accepted four cusp modes and rejected three old discrete-pencil roots that
  fail the independent test.
- Continued the leading pole in all three sectors through 18 backgrounds on
  each of the physical paths `mu_B=0` and `mu_B/T=2`.
- Achieved source residuals below `2.2e-8`, minimum parent overlap `0.9985`,
  and pointwise shooting/spectral agreement within `0.012` in `omega_hat`.
- Rewrote the manuscript around the validated homogeneous result and removed
  the old fixed-horizon QNM survey and finite-momentum claims.
- Audited the 2024--2026 literature and updated published journal metadata.
- Recomputed the earlier 2018 EMD model with its published potential, coupling,
  and scale. All nine quoted relaxation times are reproduced within
  `0.0111 fm/c`; the maximum source residual is `7.42e-9` and the maximum
  spectral/shooting distance is `0.00731` in `omega_hat`.
- Added a direct comparison at `(T,mu_B)=(400,0)` and `(145,0) MeV` and at the
  two model-specific critical endpoints.
- Compiled the JHEP PDF with Tectonic; no undefined references, undefined
  citations, duplicate labels, or overfull boxes remain. All pages were
  rendered with Poppler and visually inspected on 2026-08-25.
- Compared `jheppub.sty` and `JHEP.bst` byte-for-byte with the official JHEP
  downloads on 2026-08-25; both match, including bibliography style 2.18.
- Replaced free-form keywords by four entries from the current controlled JHEP
  keyword list and added the journal-required disclosure of AI-assisted work.
- Added a minimal source-archive builder that includes `main.bbl`, the BibTeX
  inputs, the exact official styles, and only the eight referenced PDF figures;
  it refuses to package a manuscript with incomplete author/arXiv metadata.
- Re-emitted the two remaining Matplotlib Type-3 figures with embedded
  TrueType fonts. All manuscript and figure PDFs now pass the font audit, and a
  13-file draft source archive compiles in isolation.
- Added `docs/claim_evidence_audit.md`, mapping every claim class to its
  machine-readable evidence and recording the explicit exclusion boundary for
  failed or unimplemented coupled physics.

## Remaining acceptance gates

- [x] Separate radial-resolution and radial-domain convergence studies.
- [x] Source factoring for every homogeneous master field.
- [x] Left/right eigenvectors and eigenvalue condition estimates.
- [x] Eigenfunction-overlap trajectory tracking.
- [x] Independent complex-shooting confirmation.
- [x] Homogeneous singlet derivation and symbolic validation.
- [x] Coupled helicity-one and helicity-zero equations with symbolic constraint tests
  (reproducibility artifact; numerical spectra excluded from Route B).
- [x] Physical thermodynamic trajectories and branch labels.
- [x] Explicit MAP-only framing in title, abstract, results, and conclusions.
- [x] Quantitative solver benchmark and cross-model comparison with the 2018
  homogeneous EMD analysis.
- [x] Official JHEP styles, controlled keywords, AI disclosure, and source
  archive contents verified against the current author instructions.
- [ ] Author, affiliation, email, ORCID, funding, and arXiv metadata supplied by
  the authors.

The scientific and numerical Route-B gates are closed. Submission still needs
author-supplied names, affiliations, emails, ORCIDs, acknowledgments/funding,
and arXiv metadata, followed by publication of this reconstructed branch and
its machine-readable outputs at the repository named in the manuscript.
