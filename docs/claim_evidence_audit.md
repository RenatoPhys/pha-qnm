# Claim-to-evidence audit

This is the internal hostile-referee audit for the Route-B manuscript. It maps
each claim class to an authoritative, versioned artifact and states what was
actually accepted. A claim is excluded from the paper when the corresponding
gate is not closed.

| Claim in the manuscript | Primary evidence | Independent/checking evidence | Acceptance boundary |
|---|---|---|---|
| The calculation uses the maximum-likelihood PHA realization (`sample74`). | `configs/pha_map.yaml`, `configs/pha_map.generated.yaml`, `data/manifest.json` | Native HDF5 selection and binary64 comparison recorded in `docs/STATUS.md` | MAP point only; no posterior credible interval is claimed. |
| The spectral background is an independently located cusp near `(103.755, 593.323) MeV`. | `results/python/reference_cusp_uv1e5_h10.json` | Constant-`phi0` crossings through 800 lines in `results/python/reference_cep_neighbors_N800.json` | The distinct HDF5 metadata value is retained as provenance, not treated as the spectral background. |
| The MAP model has continued coexistence and two spinodal branches above the cusp. | `results/python/phase_lines.csv`, `results/python/phase_diagram_summary.json` | `analysis/validate_phase_diagram.py` and the two coarsened-grid comparison tables | Only equilibrium branch structure is claimed; no spinodal growth rate or unstable QNM band is inferred. |
| All three homogeneous `SO(3)` sectors have derived source/vev maps and tested identities. | `derivations/generated/helicity0_equations.json`, `derivations/generated/helicity1_equations.json`, generated C++ headers | Ten tests in `derivations/validate_linearized_emd.py` cover constraints, Bianchi/Maxwell identities, gauge modes, neutral and `k=0` limits, and UV powers | Generated finite-momentum systems are reproducibility artifacts; their numerical spectra are excluded. |
| Four cusp poles pass both continuum and source tests; three former candidates do not. | `results/python/qnm_validated_modes.csv`, `results/python/qnm_validation_summary.json` | `results/python/qnm_shooting_crosscheck.csv`, separate resolution/domain tables, and horizon/UV shooting variants | Only the four independently confirmed poles enter table 2 or scientific claims. |
| The leading homogeneous poles are continuously followed on two physical paths. | `results/python/homogeneous_physical_trajectories.csv` and `.json` | 108 pointwise shooting roots, parent eigenfunction overlaps, and spectral distances stored row by row | Paths are `mu_B=0` and `mu_B/T=2`; the discarded fixed-horizon scan is not presented as temperature dependence. |
| The earlier 2018 EMD relaxation times are reproduced before comparing models. | `results/python/rougemont2018_benchmark.csv` and `.json` | Local cusp reconstruction, independent shooting, and the checksum-verified official source/PDF recorded by `analysis/compare_2018_emd.py` | Nine rounded published values are reproduced within `0.0111 fm/c`; the comparison is not a universality or critical-exponent claim. |
| Neutral entropy and baryon susceptibility are numerically stable. | `results/python/thermodynamic_validation.csv` | Three centered finite-difference steps and background constraint/Gauss diagnostics | This is a neutral-equilibrium cross-check, not a finite-density CEP susceptibility claim. |
| The manuscript does not claim shear, sound, diffusion, critical slowing, or spinodal dynamics. | Scope statements in `paper/main.tex` and `docs/derivation_finite_k.md` | Failed primitive finite-momentum constraint/continuum gate documented in `docs/STATUS.md` | Route A remains future work; no failed coupled result appears in a table or figure. |
| The submitted artifact follows current JHEP formatting requirements. | `paper/main.tex`, `paper/main.bbl`, `paper/main.pdf` | Official-style SHA-256 comparison, PDF structural/font/link checks, Poppler rendering, grayscale inspection, and guarded source-archive builder | Author/arXiv metadata and publication of the reconstructed code release remain external gates. |

## Numerical assertions rerun on 2026-08-25

- accepted independently shot cusp roots: `4`;
- physical-trajectory QNM rows: `108`;
- maximum trajectory source residual: `2.1872e-8`;
- minimum parent eigenfunction overlap: `0.998545`;
- 2018/PHA comparison rows: `18`;
- maximum difference from the nine rounded 2018 relaxation times:
  `0.011018 fm/c`;
- maximum comparison source residual: `7.4160e-9`;
- symbolic perturbation tests: `10/10` passed;
- portable C++ test target: `1/1` passed.

## Submission-only gates

No scientific claim remains unsupported by the stored evidence. Literal
submission is nevertheless prohibited until the author-controlled checklist in
`paper/SUBMISSION_CHECKLIST.md` is closed. In particular, an unpushed local
branch is not treated as public data availability.
