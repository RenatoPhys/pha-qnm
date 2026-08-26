# JHEP submission checklist

Checked against the official JHEP author instructions and tool downloads on
2026-08-25.

## Closed locally

- `jheppub.sty` and `JHEP.bst` are byte-identical to the official downloads;
  bibliography style 2.18 is used.
- The Route-A abstract fits on page 1 and contains no citations or displayed
  formulae.
- The four keywords are controlled JHEP entries.
- The acknowledgments contain the required AI-assisted-technology disclosure.
- The complete coupled finite-momentum, critical, spinodal, and posterior
  claims map to versioned CSV/JSON evidence and independent checks.
- The source build generates `main.bbl`; the guarded packager requires it.
- The master file is `main.tex`, with exactly four referenced PDF figures under
  `figures/`: `phase_diagram.pdf`, `hydrodynamic_critical_dispersion.pdf`,
  `cep_spinodal_dynamics.pdf`, and `posterior_uq.pdf`.
- The manuscript compiles without undefined citations/references, duplicate
  labels, overfull boxes, or TeX errors.
- All manuscript/figure fonts are embedded and no Type 3 fonts are present.
- All 16 pages have been rendered and visually inspected for clipping,
  overlap, figure legibility, and bibliography layout.

## Author-controlled items still required

- Insert every author's name, affiliation, email, and corresponding-author
  designation in `main.tex`; update `CITATION.cff` and add ORCIDs.
- Add funding identifiers and additional acknowledgments, or confirm none.
- Obtain and insert the arXiv identifier; keep the JHEP and arXiv versions in
  sync.
- Approve authorship, originality, institutional authorization, and exclusive
  submission declarations.
- Select the submission-form keywords and data/code-availability statement.
- Archive the release persistently and provide its DOI when required.

## Final packaging command

```powershell
powershell -ExecutionPolicy Bypass -File paper/build_submission_archive.ps1
```

The command deliberately fails while author, affiliation, email, arXiv, or
`CITATION.cff` authorship metadata are missing.

Official sources:

- <https://jhep.sissa.it/jhep/help/JHEP/JHEP_author.jsp>
- <https://jhep.sissa.it/jhep/help/JHEP_TeXclass.jsp>
- <https://jhep.sissa.it/jhep/help/keywordsList.jsp>
