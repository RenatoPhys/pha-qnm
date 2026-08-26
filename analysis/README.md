# Python numerical analysis

This directory is the reproducible Python analysis layer.  It independently
integrates the same backgrounds as the C++ reference kernels, solves the
decoupled and coupled fluctuation systems, validates the poles with two radial
methods, and produces the manuscript tables and figures.

```powershell
python -m pip install -r analysis/requirements.txt
python analysis/run_numerics.py
python analysis/validate_decoupled_qnms.py
python analysis/build_homogeneous_trajectories.py
python analysis/compare_2018_emd.py
python analysis/validate_coupled_qnms.py
python analysis/run_finite_k_physics.py --stage all
python analysis/run_posterior_uq.py --samples 25 --workers 4
python analysis/plot_posterior_uq.py
```

Outputs are written to `results/python/` and `paper/figures/`.  Accepted
finite-momentum roots pass source-determinant, gauge, primitive-constraint,
continuum, horizon-start, UV-domain, and independent-method gates.  The main
drivers are checkpointed and store their thresholds and error budgets in JSON
summaries.  The older homogeneous pipeline remains as a nonhydrodynamic and
cross-model benchmark rather than as the scientific scope boundary.
