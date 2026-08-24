# Manuscript

`main.tex` uses the standard JHEP article layout and `JHEP.bst`. Build from this directory with:

```text
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

The manuscript contains numerical tensor and homogeneous-vector candidates, background scans, and validation diagnostics. It is not submission-ready until the remaining shooting, scalar, and coupled-helicity gates in `docs/STATUS.md` are complete. Author metadata is intentionally left for the user to supply.
