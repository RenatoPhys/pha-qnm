# Result products

The publication-authoritative Route-B outputs are the versioned CSV and JSON
files under `results/python/`. Exact binary64 values, numerical settings,
background coordinates, branch labels, QNM residuals, condition estimates,
tracking parents, overlaps, and cross-check envelopes remain machine readable;
the manuscript rounds them according to the observed numerical uncertainty.

Principal products are:

- `phase_surface.csv` and `phase_surface.npz`: the 141 by 121 horizon-data
  thermodynamic surface;
- `phase_lines.{csv,json}` and `phase_line_grid_convergence.{csv,json}`:
  fixed-chemical-potential branches, spinodals, coexistence, and grid errors;
- `qnm_validated_modes.csv`, `qnm_factored_convergence.csv`, and
  `qnm_shooting_*.csv`: accepted homogeneous modes and independent checks;
- `homogeneous_physical_trajectories.{csv,json}`: the two controlled PHA paths;
- `rougemont2018_benchmark.{csv,json}`: the independently recomputed 2018 EMD
  benchmark and its direct comparison with the PHA MAP realization.

The C++/HDF5 schema documented in the original prototype remains a future
finite-momentum database design; it is not used to support claims in the
homogeneous manuscript.
