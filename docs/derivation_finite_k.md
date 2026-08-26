# Finite-momentum perturbations

Momentum points along `z`, leaving `SO(2)` helicities. Helicity two reduces
exactly to the tensor equation with the additional `-exp(-2A) k^2/h` term.

The coupled systems are derived in ingoing EF coordinates with radial gauge:

- helicity one: `{H_vx,H_zx,a_x}`;
- helicity zero: `{H_vv,H_vz,H_zz,H_aa,a_v,a_z,varphi}`, with
  `H_aa=H_xx=H_yy`.

The generator `derivations/generate_linearized_emd.py` expands the complete
trace-reversed Einstein equations, Maxwell equations, and scalar equation
using exact first-order dual expressions. It emits raw coefficient JSON,
content hashes, and C++23 local radial kernels. The helicity-one output contains
three evolution equations plus `E_vx` and the raised radial Einstein
constraint. The helicity-zero output contains seven frequency-linear
characteristic equations and seven unused Einstein/Maxwell diagnostics.

Residual radial-gauge transformations are retained and tested. In the
transverse sector, `X'=0` gives the pure-gauge mode
`(H_vx,H_zx,a_x)=(-sigma X,-p X,0)`. In helicity zero the residual parameters
satisfy

```text
T' = 0,
R' = -A' R - sigma T,
Z' = -p exp(-A) T,
lambda' = 0.
```

All generated equations annihilate these modes on a background satisfying the
EMD equations. The same validation suite checks the geometric Bianchi identity,
the Maxwell identity, transverse parity, the neutral limit, the `k -> 0`
triplet and singlet reductions, and the correct UV source/vev powers.

The primitive UV fields have source/vev powers `u^0/u^4` for metric
perturbations, `u^0/u^2` for Maxwell perturbations, and `u^nu/u^Delta` for the
dilaton. Radial-gauge residuals reduce the three helicity-one primitives to two
physical degrees of freedom and the seven helicity-zero primitives to three.
The numerical implementation constructs physical gauge invariants from these
primitives, quotients the residual-gauge horizon data, and evaluates the unused
primitive equations as constraints. Source-matrix shooting and the independent
multidomain Chebyshev integrator are implemented in `analysis/coupled_qnm.py`;
their automated continuum and constraint gates are recorded in
`results/python/coupled_qnm_validation_summary.json`.
