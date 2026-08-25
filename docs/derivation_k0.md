# Homogeneous perturbations

At zero momentum, rotational invariance gives the `SO(3)` quintuplet, triplet,
and singlet. The tensor and vector equations are stored as operator pencils
only after multiplication by `h`, preserving a regular horizon row and
linearity in frequency.

The singlet invariant is

```text
S = varphi - [phi'/(2 A')] (H_xx + H_yy + H_zz)/3 .
```

The implemented master variable is `psi=S/Y`, where `Y=phi'/A'`. Substituting
this field redefinition into eq. (16) of arXiv:1804.00189 reproduces the
compiled scalar pencil exactly after the background equations are used. That
comparison is automated in `derivations/validate_linearized_emd.py`; it is no
longer an unverified equation from the research brief.

In the ultraviolet, `S` has powers `u^nu` and `u^Delta`. Consequently `psi`
has a constant source and a normalizable power `u^(2 Delta - 4)`. The QNM tau
row sets the constant source to zero. At nonzero frequency the homogeneous
energy and charge perturbations are constrained, leaving this single physical
singlet degree of freedom.
