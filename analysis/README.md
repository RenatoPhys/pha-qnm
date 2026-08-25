# Python numerical analysis

This directory is an explicit user-authorized Python analysis layer. It does
not replace the C++ reference kernels; it independently integrates the same
background equations, assembles the currently derived decoupled pencils, runs
resolution/cutoff and shooting tests, and produces the manuscript figures and
tables.

```powershell
python -m pip install -r analysis/requirements.txt
python analysis/run_numerics.py
python analysis/validate_decoupled_qnms.py
python analysis/build_homogeneous_trajectories.py
python analysis/compare_2018_emd.py
```

Outputs are written to `results/python/` and `paper/figures/`. The accepted
homogeneous quintuplet, triplet, and singlet modes are source factored and
confirmed by complex shooting. Generated finite-momentum helicity-one and
helicity-zero equations are retained as symbolic artifacts, but their coupled
spectra are not reported because the numerical gates have not passed.
