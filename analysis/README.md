# Python numerical analysis

This directory is an explicit user-authorized Python analysis layer. It does
not replace the C++ reference kernels; it independently integrates the same
background equations, assembles the currently derived decoupled pencils, runs
resolution/cutoff tests, and produces the manuscript figures and tables.

```powershell
python -m pip install -r analysis/requirements.txt
python analysis/run_numerics.py
```

Outputs are written to `results/python/` and `paper/figures/`. The script does
not report scalar, helicity-one, or helicity-zero spectra because those coupled
equations have not passed the derivation gates.
