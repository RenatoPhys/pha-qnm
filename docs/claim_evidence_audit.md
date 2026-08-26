# Claim-to-evidence audit

This is the internal hostile-referee audit for the Route-A manuscript.

| Claim class | Primary evidence | Independent/checking evidence | Acceptance boundary |
|---|---|---|---|
| Exact PHA posterior and MAP provenance | `configs/pha_map*.yaml`, `data/manifest.json`, `posterior_selection.csv` | Native HDF5 selector and binary64 comparison in `STATUS.md` | Upstream CEP metadata difference is disclosed, not absorbed into error bars. |
| MAP cusp, coexistence, and spinodals | `reference_cusp_*.json`, `phase_lines.csv`, `phase_diagram_summary.json` | Constant-`phi0` route, independent fold/null equations, coarsened grids, `validate_phase_diagram.py` | Thermodynamic and numerical errors remain distinct. |
| Neutral thermodynamic reproduction | `thermodynamic_validation.csv` | Three centered charge steps, charge-conjugation parity, Einstein/Gauss diagnostics, appendix figure | This is a neutral EOS/susceptibility validation, not a finite-density critical claim. |
| Complete coupled field content | Generated helicity-zero/one JSON and C++ kernels | `validate_linearized_emd.py` (10/10) | Gauge or constraint failures exclude a root. |
| Physical coupled QNMs | `coupled_qnm.py`, `coupled_qnm_validation.csv` | DOP853 source-matrix shooting versus multidomain Chebyshev integration; independent `N`, domain, partition, and `r0` variations | Stored thresholds cover source determinant, gauge defect, primitive constraints, continuum, and cross-method distance. |
| Stable shear, sound, and diffusion | `finite_k_hydrodynamic_dispersion.csv` | Independent thermodynamic/horizon formulae in `charged_hydrodynamics.py` and six benchmark cases | Hydrodynamic coefficients are never fitted back from the QNMs. |
| Critical slowing and direct dynamic exponent | `cep_critical_scaling.csv`, `cep_q4_dispersion.csv` | Five thermodynamic windows, four direct-QNM windows, homogeneous gapped modes, fine-grid/shooting cross-check | The paper quotes the dominant `+/-0.05` window systematic. |
| Closed spinodal unstable band | `spinodal_dispersion.csv` | Three positions between folds, direct shooting and fine-grid checks near the maximum | MAP band scales are not promoted to posterior full-QNM intervals. |
| Bayesian robustness | `posterior_uq_samples.json`, `posterior_uq_summary.json` | 25/25 independently reconstructed medoids, full represented weight, adaptive-step audit | Claims cover `T_c`, `mu_c`, `z`, `D_B`, and spinodal width; posterior full-QNM bands are not claimed. |
| Nonhydrodynamic benchmark | `qnm_validated_modes.csv`, homogeneous trajectories, 2018 benchmark CSV/JSON | Source-factored collocation and independent shooting | Used as validation and context, not to infer spinodal growth. |
| JHEP artifact integrity | `paper/main.tex`, `main.bbl`, `main.pdf` | compilation log, font/structure/render checks, guarded packager | Author/arXiv/funding/release metadata remain external submission gates. |

## Numerical headline values

- Maximum shooting/collocation distance in the coupled benchmark:
  `3.63e-12` in `omega/(2 pi T)`.
- Maximum continuum displacement: `3.64e-11`.
- Maximum source determinant: `1.60e-12`.
- Maximum primitive-constraint diagnostic: `2.41e-5` under the deliberately
  doubled horizon-start stress test; the production roots are lower.
- Largest stable hydrodynamic comparison difference at `q=0.05`: `0.166%`.
- Critical thermodynamic exponent: `z=4.08+/-0.05`; direct QNM: `z=3.99(2)`.
- Posterior medoids: `25/25` successful, representing all successful draws.

No scientific claim in the manuscript depends only on an eigensolver residual,
a single discretization, or an unversioned terminal observation.
