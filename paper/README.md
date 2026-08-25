# Manuscript

`main.tex` uses the standard JHEP article layout and `JHEP.bst`. On Windows,
the reproducible build command is:

```text
powershell -ExecutionPolicy Bypass -File build_pdf.ps1
```

The manuscript contains the reconstructed phase structure, independently
validated quintuplet/triplet/singlet poles, and controlled homogeneous
trajectories at `mu_B=0` and `mu_B/T=2`. It also includes a reproducible
three-sector comparison with Phys. Rev. D 98 (2018) 034028 at two common
neutral states and at the two model-specific endpoints. It is framed as a
maximum-likelihood (MAP) analysis, not as posterior uncertainty quantification.
Scientific and numerical Route-B gates are closed; author, affiliation, ORCID,
funding, acknowledgment, and arXiv metadata remain for the authors to supply.

The current JHEP compliance audit and author-controlled submission items are
listed in `SUBMISSION_CHECKLIST.md`. After those fields are filled, create the
minimal source archive with:

```text
powershell -ExecutionPolicy Bypass -File build_submission_archive.ps1
```

The packager refuses incomplete author/arXiv metadata and includes the required
generated `main.bbl` plus only files needed to compile the article.
