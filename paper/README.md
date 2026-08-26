# Manuscript

`main.tex` is the Route-A JHEP manuscript **Critical relaxation and spinodal
instabilities in Bayesian holographic QCD**. It presents the complete coupled
helicity-one/helicity-zero calculation, stable charged-hydrodynamic matches,
critical scaling, closed spinodal unstable bands, and a weighted 25-medoid
posterior analysis. Four referenced vector figures carry the scientific
narrative; detailed diagnostics remain machine readable.

Build on Windows with:

```text
powershell -ExecutionPolicy Bypass -File build_pdf.ps1
```

The scientific, numerical, layout, and local JHEP gates are closed. Author,
affiliation, email, ORCID, funding, acknowledgment, arXiv, and release metadata
must be supplied by the authors and are intentionally not fabricated.

After those fields are filled, create the minimal source archive with:

```text
powershell -ExecutionPolicy Bypass -File build_submission_archive.ps1
```

The packager reads figure dependencies directly from `main.tex`, requires the
generated `main.bbl`, refuses incomplete author/arXiv metadata, and includes
only compilation inputs.
