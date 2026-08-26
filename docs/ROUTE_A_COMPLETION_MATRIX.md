# Route-A completion matrix

This matrix implements the Route-A definition of done in
`JHEP_UPGRADE_PLAN_PHA_QNM.md`. A closed item has a versioned,
machine-readable artifact and a validation gate; terminal observations alone
do not count.

| Priority | Required outcome | Evidence | Status |
|---|---|---|---|
| P0 | Preserve the pre-rebuild state and provenance. | `pre-jhep-rebuild` tag, `docs/baseline_manifest.txt`, data manifest. | **Closed.** |
| P1 | Reconstruct MAP thermodynamics, independently locate the cusp, and continue coexistence and both spinodals. | `reference_cusp_*.json`, `phase_surface.*`, `phase_lines.*`, `phase_diagram_summary.json`, `validate_phase_diagram.py`. | **Closed.** The 9.152 MeV difference from the HDF5 CEP metadata is retained as an upstream provenance discrepancy rather than hidden by retuning. |
| P2 | Derive complete finite-momentum helicity-one/helicity-zero systems, gauge transformations, constraints, and source/vev maps. | Generated JSON/C++ kernels, `validate_linearized_emd.py` (10/10), `derivation_finite_k.md`, gauge invariants in `coupled_qnm.py`. | **Closed.** |
| P3 | Obtain continuum coupled QNMs with source factoring, constraints, conditioning, tracking, and an independent method. | `coupled_qnm.py`, `validate_coupled_qnms.py`, `coupled_qnm_validation.csv` and summary. | **Closed.** Six neutral/charged shear, sound, and diffusion roots pass DOP853 source-matrix shooting and independent multidomain Chebyshev integration, with independent resolution, domain, and horizon-start variations. |
| P4 | Identify stable shear, sound, and baryon/heat diffusion and match charged hydrodynamics. | `charged_hydrodynamics.py`, `finite_k_hydrodynamic_dispersion.csv`, `finite_k_physics_summary.json`, figure 2. | **Closed.** The neutral and charged coefficients agree at the reported small momenta; the largest benchmark relative difference is below 0.17%. |
| P5a | Quantify critical slowing with controlled windows and nonhydrodynamic behavior. | `cep_critical_scaling.csv`, `cep_q4_dispersion.csv`, JSON fit windows, homogeneous cusp modes, figures 2--3. | **Closed.** Thermodynamic scaling gives `z=4.08+/-0.05`; direct low-momentum QNMs give `z=3.99(2)`, while homogeneous gapped modes remain finite. |
| P5b | Resolve the spinodal unstable band, fastest scale, and hydrodynamic breakdown. | `spinodal_dispersion.csv`, three branch locations in `finite_k_physics_summary.json`, shooting/fine-grid cross-checks, figure 3. | **Closed.** All three scans form closed bands; the midpoint gives `q*=0.283`, `Gamma*=1.41e-3`, and `q_edge=0.417`. |
| P6 | Propagate posterior uncertainty or remove Bayesian claims. | `posterior_selection.csv`, `posterior_uq_samples.json`, `posterior_uq_summary.json`, figure 4. | **Closed.** Twenty-five deterministic weighted medoids represent all 1589 successful released draws and full posterior weight; 25/25 reconstruct their local cusp and spinodal pair. |
| P7 | Rewrite around Route A and pass the JHEP/internal-referee audit. | `paper/main.tex`, `paper/main.pdf`, four claim-bearing figures, `claim_evidence_audit.md`, submission checklist and guarded packager. | **Scientifically and editorially closed.** Author identity, affiliation, funding, ORCID, arXiv, and release metadata remain external author-controlled submission fields. |

## Acceptance rules enforced

1. A coupled pole enters the paper only after its source determinant, gauge
   defect, primitive constraints, continuum envelope, and independent radial
   method pass their stored thresholds.
2. Hydrodynamic coefficients come independently from equilibrium derivatives
   and horizon formulae; they are not fitted back from QNM data.
3. Critical exponents retain their window systematic and are not quoted with
   spurious precision.
4. Spinodal growth is claimed only from a closed full-QNM unstable band.
5. Posterior and numerical/window uncertainties are reported separately.

The remaining author-controlled fields do not alter the scientific completion
of Route A and are intentionally not replaced by fictitious metadata.
