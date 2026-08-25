# Generated kernels

`helicity1_equations.json` and `helicity0_equations.json` are the auditable raw
radial coefficient maps emitted by `../generate_linearized_emd.py`.
`manifest.json` records SHA-256 hashes for these files, the generator, and the
generated C++23 headers under `include/pha_qnm/perturbations/generated/`.

Regenerate and validate with:

```text
python derivations/generate_linearized_emd.py
python derivations/validate_linearized_emd.py
```

Generated kernels are accepted only together with the symbolic validation and
compiled `generated_kernels` unit test. Coupled spectra remain a separate
solver/constraint acceptance gate.
