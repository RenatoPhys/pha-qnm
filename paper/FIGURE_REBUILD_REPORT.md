# PHA-QNM figure and editorial rebuild report

Date: 2026-08-26  
Branch: `codex/figure-editorial-rebuild` from `origin/codex/jhep-rebuild`

## Outcome

The figure/editorial rebuild is complete. The longitudinal carrier is audited
as the charged baryon/heat-diffusion branch, all reader-facing figures are
script-generated PDF/PNG pairs, the exploratory background figure is replaced
by a horizon-grid/thermodynamic-map figure, and the 2018 appendix now contains
curves rather than only nine tabulated landmarks. The manuscript builds to a
29-page PDF with no undefined references, undefined citations, overfull boxes,
duplicate labels, or missing figures.

## Scientific decision

The `mu_B=650 MeV` audit uses nine backgrounds, the complete charged-fluid
matrix, relative-variable eigenvectors, QNM continuation down to
`mathfrak q=0.0025`, and an explicit first nonhydrodynamic check. It finds:

- `c_s^2 > 0` throughout and a stable sound pair;
- `D_B < 0` only between the folds;
- median QNM small-momentum power `1.990`;
- decisive BIC preference for `q^2+q^4` over `q+q^2` at all three interior
  backgrounds;
- a separated, gapped first nonhydrodynamic pole.

Terminology is therefore `baryon/heat-diffusion branch`, with the charged
mechanism and contrast to neutral sound-driven spinodals stated explicitly.

## Figure changes

| Figure | Reader-facing change | Inputs |
|---|---|---|
| `phase_diagram` | spinodal envelope, stability semantics, both dynamical paths, three scan backgrounds, CEP inset | phase lines; finite-k summary |
| `background_grid_map` | new horizon grid, failures, `det J=0`, constant-mu contours, mapped grid families | phase surface; phase lines; audit |
| `hydrodynamic_critical_dispersion` | independent hydro curves, QNM markers, residual panel, direct `D_B`-`chi_B` fit | stable dispersion; critical scaling |
| `cep_spinodal_dynamics` | visible fit window, `z_eff` inset, unstable half-plane, full-QNM/interpolation/hydro separation, MAP band landmarks | CEP and spinodal CSVs |
| `posterior_uq` | consistent weight/color encoding, weighted medians, rank forest, separate common systematic | posterior samples and summary |
| `longitudinal_mode_identity_audit` | new coefficient, power, eigenvector, and pole-separation evidence | dedicated audit CSV/JSON |
| `homogeneous_qnm_trajectories` | 2x2 real/damping/complex/tau layout | validated homogeneous trajectories |
| `legacy_2018_qnm_curves` | new 3x2 reproduction at five chemical potentials | new 2018-model curve CSV; unchanged benchmark |
| `cep_reproduction` | differences relative to local cusp | stored crossing/cusp results |
| `qnm_convergence` | separated studies, accepted envelope, root floor | convergence and validated-mode CSVs |
| `pha_model_kernels` | sampled horizon-field range and resolved UV inset | analytic MAP kernels; phase-surface range |
| `representative_background` | normalized profiles, UV-layer inset, pointwise and envelope flux error | MAP critical background |
| `thermodynamic_validation` | redundant step encoding and percent displacement panel | unchanged thermodynamic validation CSV |

The exact semantic encodings, axes, fit windows, interpolation roles, captions,
and dependencies are machine-readable in `analysis/figure_manifest.json`.

## Editorial changes

- Abstract reduced to the requested scale and limited to three numerical
  groups.
- A three-claim novelty paragraph was added at the end of the introduction.
- The mode audit and its limitations are included before the thermodynamic
  reproduction appendix.
- All principal captions now distinguish points, predictions, fits,
  interpolations, and bands.
- The model-specific nature of `z approximately 4` is explicit; the manuscript
  does not present it as a direct QCD determination.
- `background_thermodynamics.pdf` remains in the repository as historical
  provenance but is no longer included in the manuscript.

## New data

- `results/python/longitudinal_mode_identity_audit.csv`
- `results/python/longitudinal_mode_identity_audit.json`
- `results/python/longitudinal_mode_identity_checkpoint.json`
- `results/python/legacy_2018_qnm_curves.csv`
- `results/python/legacy_2018_qnm_curves_summary.json`

The legacy reproduction contains 195 QNM points on 65 backgrounds; its maximum
pencil residual is `1.12e-15`. The audit's largest source determinant is
`1.08e-9` and its largest primitive constraint is `3.26e-6`.

## Numerical regression

No pre-existing production CSV/JSON changed. Baseline and final Git blob hashes
are identical for the finite-k, posterior, coupled-validation,
thermodynamic-validation, phase-summary, and 2018-benchmark products. Full
details are in `results/python/figure_rebuild_scientific_change_log.md`.

## Commands executed

```text
python analysis/audit_longitudinal_mode_identity.py
python analysis/plot_legacy_2018_curves.py
python analysis/plot_appendix_diagnostics.py
python analysis/render_all_figures.py --gallery --check
powershell -ExecutionPolicy Bypass -File paper/build_pdf.ps1
pdftoppm -png -r 110 paper/main.pdf tmp/pdfs/manuscript-rebuild/page
pdfinfo paper/main.pdf
python -m pytest -q tests/test_figure_contracts.py
powershell -ExecutionPolicy Bypass -File paper/build_submission_archive.ps1
```

The submission-archive dry run correctly stops before writing an archive
because author, affiliation, email, arXiv, and `CITATION.cff` author metadata
are intentionally absent.

## QA

- Figure contracts: `12 passed`.
- Manuscript: 29 A4 pages, 459,675 bytes, PDF 1.5.
- All manuscript fonts are embedded (recursive Type-0/descendant-font check).
- Every page was rendered with Poppler and inspected in three contact sheets.
- All manifested figures were inspected at normal size, thumbnail size, and
  grayscale in `results/figure_gallery/index.html`.
- No clipped labels, overlapping legends, corrupt pages, or unreadable
  grayscale encodings remain.

## Before/after gallery

Tracked before snapshots are in `results/figure_gallery/before/`. Final normal,
thumbnail, and grayscale renders plus `contact_sheet.png` are in
`results/figure_gallery/`. Open `results/figure_gallery/index.html` for direct
comparison.

## Remaining warnings

- Tectonic reports a missing host Fontconfig default file, but uses embedded
  document fonts and completes successfully.
- TeX reports two underfull horizontal boxes; visual inspection shows no
  defect. There are no underfull vertical or overfull boxes.
- `jheppub` warns that author e-mail metadata are missing, as intentionally
  documented in the manuscript.
- `hyperref` reports that the legacy `pagecolor` option is unavailable; links
  and PDF metadata render correctly.
- The source archive remains gated on author-controlled identity and release
  metadata. No push or merge was performed.
