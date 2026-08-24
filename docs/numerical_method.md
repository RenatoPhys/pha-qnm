# Numerical method

Backgrounds are initialized from a compiled fourth-order regular horizon series and integrated in horizon-data coordinates `(phi0, Phi1)`. The reference formulation evolves the Maxwell field and monitors `Q=exp(2A) f Phi'`; the reduced formulation replaces `Phi'` by the conserved charge. UV observables are not accepted from arbitrary unconstrained polynomial fits.

QNMs use ingoing Eddington--Finkelstein coordinates. After multiplying by the blackening factor, decoupled equations are linear pencils `M0 + omega M1`; source and regular horizon tau rows are imposed directly. Multidomain Chebyshev interfaces are placed around the `d2 phi ~ 1` boundary layer. The production dense path calls LAPACK QZ and the large-scan path targets the same pencil with SLEPc shift-invert.

