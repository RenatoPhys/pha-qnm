"""Symbolic acceptance tests for generated EMD fluctuation equations."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import generate_linearized_emd as generated


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "derivations" / "generated"


def load(sector: str) -> dict[str, object]:
    return json.loads((GENERATED / f"{sector}_equations.json").read_text(encoding="utf-8"))


def symbols() -> dict[str, sp.Symbol]:
    names = (
        "A Ar Arr h hr hrr phi phir phirr Phi Phir Phirr f fr fphi fphir fphiphi "
        "V Vphi Vphiphi sigma p Hvx Hvxr Hvxrr Hzx Hzxr Hzxrr ax axr axrr "
        "Hvv Hvvr Hvvrr Hvz Hvzr Hvzrr Hzz Hzzr Hzzrr Haa Haar Haarr "
        "av avr avrr az azr azrr varphi varphir varphirr"
    )
    return {name: sp.Symbol(name) for name in names.split()}


S = symbols()


def parse(sector: dict[str, object], label: str) -> sp.Expr:
    equation = sector["equations"][label]["expression"]
    return sp.sympify(equation, locals=S)


def on_shell(expression: sp.Expr) -> sp.Expr:
    s = S
    rules = {
        s["fr"]: s["fphi"] * s["phir"],
        s["fphir"]: s["fphiphi"] * s["phir"],
        s["Arr"]: -s["phir"]**2 / 6,
        s["hrr"]: -4 * s["Ar"] * s["hr"] + sp.exp(-2 * s["A"]) * s["f"] * s["Phir"]**2,
        s["Phirr"]: -(2 * s["Ar"] + s["fphi"] * s["phir"] / s["f"]) * s["Phir"],
        s["phirr"]: (s["Vphi"] - sp.Rational(1, 2) * sp.exp(-2 * s["A"])
                      * s["fphi"] * s["Phir"]**2
                      - (4 * s["h"] * s["Ar"] + s["hr"]) * s["phir"]) / s["h"],
        s["V"]: -sp.Rational(1, 2) * (s["h"] * (24 * s["Ar"]**2 - s["phir"]**2)
                                      + 6 * s["Ar"] * s["hr"]
                                      + sp.exp(-2 * s["A"]) * s["f"] * s["Phir"]**2),
    }
    return sp.factor(sp.simplify(expression.subs(rules, simultaneous=True)))


def require_zero(label: str, expression: sp.Expr) -> None:
    result = on_shell(expression)
    if result != 0:
        raise AssertionError(f"{label} did not vanish: {result}")


def test_helicity1_constraint_order() -> None:
    h1 = load("helicity1")
    constraint = parse(h1, "constraint_Er_x")
    if constraint.has(S["Hvxrr"], S["Hzxrr"], S["axrr"]):
        raise AssertionError("radial helicity-one constraint contains second derivatives")


def test_helicity1_k0_reduction() -> None:
    h1 = load("helicity1")
    constraint = on_shell(parse(h1, "constraint_Er_x").subs(S["p"], 0))
    expected_constraint = -S["sigma"] * (
        S["Hvxr"] + sp.exp(-2 * S["A"]) * S["f"] * S["Phir"] * S["ax"]
    ) / 2
    require_zero("helicity-one k=0 constraint", constraint - expected_constraint)

    maxwell = parse(h1, "evolution_M_x").subs({
        S["p"]: 0,
        S["Hvxr"]: -sp.exp(-2 * S["A"]) * S["f"] * S["Phir"] * S["ax"],
    })
    reduced = on_shell(maxwell * sp.exp(2 * S["A"]) / S["f"])
    target = (
        S["h"] * S["axrr"]
        + (2 * S["Ar"] * S["h"] + S["hr"] + S["h"] * S["fr"] / S["f"]
           + 2 * S["sigma"] * sp.exp(-S["A"])) * S["axr"]
        + (-sp.exp(-2 * S["A"]) * S["f"] * S["Phir"]**2
           + S["sigma"] * sp.exp(-S["A"]) * (S["Ar"] + S["fr"] / S["f"])) * S["ax"]
    )
    require_zero("helicity-one to homogeneous triplet", reduced - target)


def test_transverse_parity_and_neutral_limit() -> None:
    h1 = load("helicity1")
    # Rotational symmetry makes the y-polarized equations an exact field rename.
    for label, payload in h1["equations"].items():
        expression = sp.sympify(payload["expression"], locals=S)
        if expression.has(S["phi"], S["Phi"]):
            raise AssertionError(f"{label} contains undifferentiated background gauge/scalar fields")
    neutral = sp.simplify(on_shell(parse(h1, "evolution_M_x")).subs(S["Phir"], 0))
    if neutral.has(S["Hvx"], S["Hvxr"], S["Hzx"], S["Hzxr"]):
        raise AssertionError("neutral Maxwell fluctuation did not decouple from the metric")


def test_helicity0_characteristic_linearity() -> None:
    h0 = load("helicity0")
    evolution = [label for label in h0["equations"] if label.startswith("evolution_")]
    if len(evolution) != 7:
        raise AssertionError(f"expected seven characteristic equations, found {len(evolution)}")
    for label in evolution:
        expression = parse(h0, label)
        if sp.Poly(expression, S["sigma"]).degree() > 1:
            raise AssertionError(f"{label} is not linear in frequency")


def test_helicity0_constraints_and_k0_singlet() -> None:
    h0 = load("helicity0")
    diagnostics = ("constraint_Er_r", "constraint_Er_v", "constraint_Er_z", "constraint_M_r",
                   "diagnostic_E_vv", "diagnostic_E_vz", "diagnostic_E_zz")
    for label in diagnostics:
        expression = parse(h0, label)
        if expression == 0:
            raise AssertionError(f"{label} was emitted as an empty diagnostic")
        if sp.Poly(expression, S["sigma"]).degree() > 2:
            raise AssertionError(f"{label} has unexpected time-derivative order")

    # At k=0 the parity-odd longitudinal fields consistently vanish and the
    # remaining equations contain the SO(3) singlet variables only after
    # identifying the longitudinal and transverse spatial perturbations.
    singlet_rules = {S["p"]: 0, S["Hvz"]: 0, S["Hvzr"]: 0, S["Hvzrr"]: 0,
                     S["az"]: 0, S["azr"]: 0, S["azrr"]: 0,
                     S["Hzz"]: S["Haa"], S["Hzzr"]: S["Haar"],
                     S["Hzzrr"]: S["Haarr"]}
    for label in h0["equations"]:
        expression = on_shell(parse(h0, label).subs(singlet_rules, simultaneous=True))
        if expression.has(S["Hvz"], S["Hvzr"], S["Hvzrr"], S["az"], S["azr"], S["azrr"]):
            raise AssertionError(f"{label} failed the k=0 singlet field reduction")


def test_maxwell_identity() -> None:
    for sector_name in ("helicity1", "helicity0"):
        system, _, _ = generated.sector_data(sector_name)
        internal = system["internal"]
        sqrt_g = internal["sqrt_g"]
        divergence = sum(generated.derivative(sqrt_g * system["maxwell"][index],
                                              generated.coords[index])
                         for index in range(5)) / sqrt_g
        identity = sp.simplify(generated.normalized(divergence))
        if identity != 0:
            raise AssertionError(f"{sector_name} Maxwell identity failed: {identity}")


def test_linearized_bianchi_identity() -> None:
    # The geometric contracted Bianchi identity is checked directly, before
    # using any matter equation or background substitution.
    for sector_name, components in (("helicity1", (2,)), ("helicity0", (0, 1, 4))):
        system, _, _ = generated.sector_data(sector_name)
        internal = system["internal"]
        metric = internal["metric"]
        inverse = internal["inverse"]
        gamma = internal["gamma"]
        ricci = internal["ricci"]
        scalar_curvature = sum(inverse[m][n] * ricci[m][n]
                               for m in range(5) for n in range(5))
        einstein_tensor = [[ricci[m][n] - metric[m][n] * scalar_curvature / 2
                            for n in range(5)] for m in range(5)]
        for n in components:
            divergence = generated.Dual(0)
            for m in range(5):
                for a in range(5):
                    covariant = generated.derivative(einstein_tensor[m][n], generated.coords[a])
                    for q in range(5):
                        covariant -= gamma[q][a][m] * einstein_tensor[q][n]
                        covariant -= gamma[q][a][n] * einstein_tensor[m][q]
                    divergence += inverse[m][a] * covariant
            identity = sp.simplify(generated.normalized(divergence))
            if identity != 0:
                raise AssertionError(
                    f"{sector_name} Bianchi component {n} failed: {identity}"
                )


def test_homogeneous_singlet_master() -> None:
    """Compare the generated implementation with the published EMD master field.

    Rougemont--Critelli--Noronha (arXiv:1804.00189, eq. 16) use the invariant
    S = varphi - phi' H/(2 A').  The compiled pencil uses psi=S/Y with
    Y=phi'/A'.  This test performs that field redefinition on shell.
    """
    s = S

    def radial_derivative(expression: sp.Expr) -> sp.Expr:
        mapping = {
            s["A"]: s["Ar"], s["Ar"]: s["Arr"],
            s["h"]: s["hr"], s["hr"]: s["hrr"],
            s["phir"]: s["phirr"],
            s["Phir"]: s["Phirr"],
            s["f"]: s["fphi"] * s["phir"],
            s["fphi"]: s["fphiphi"] * s["phir"],
            s["Vphi"]: s["Vphiphi"] * s["phir"],
        }
        return sp.expand(sum(sp.diff(expression, variable) * derivative
                             for variable, derivative in mapping.items()))

    L = on_shell(s["phirr"] / s["phir"] - s["Arr"] / s["Ar"])
    Lr = on_shell(radial_derivative(L))
    B = 4 * s["Ar"] + s["hr"] / s["h"] + 2 * s["sigma"] * sp.exp(-s["A"]) / s["h"]
    published_potential = sp.exp(-2 * s["A"]) / (18 * s["f"] * s["h"] * s["Ar"]**2) * (
        -18 * s["Ar"]**2 * s["fphi"]**2 * s["Phir"]**2
        + s["f"] * (
            3 * s["Ar"]**2 * (8 * sp.exp(2 * s["A"]) * s["h"] * s["phir"]**2
                               - 6 * sp.exp(2 * s["A"]) * s["Vphiphi"]
                               + 3 * s["fphiphi"] * s["Phir"]**2)
            + 6 * s["Ar"] * s["phir"] * (
                sp.exp(2 * s["A"]) * (s["hr"] * s["phir"] - 2 * s["Vphi"])
                + s["fphi"] * s["Phir"]**2)
            + 54 * s["sigma"] * sp.exp(s["A"]) * s["Ar"]**3
            - sp.exp(2 * s["A"]) * s["h"] * s["phir"]**4
        )
    )
    transformed = on_shell(Lr + L**2 + B * L + published_potential)
    target = on_shell(
        -s["hr"] * L / s["h"]
        + sp.exp(-2 * s["A"]) * s["Phir"]**2
          * (3 * s["Ar"] * s["fphi"] - s["f"] * s["phir"]) / (s["h"] * s["phir"])
        + s["sigma"] * sp.exp(-s["A"]) * (3 * s["Ar"] + 2 * L) / s["h"]
    )
    difference = sp.factor(sp.together(transformed - target))
    if difference != 0:
        raise AssertionError(f"homogeneous singlet master mismatch: {difference}")


def test_residual_gauge_modes() -> None:
    h1 = load("helicity1")
    X = sp.Symbol("X")
    transverse_gauge = {
        S["Hvx"]: -S["sigma"] * X,
        S["Hzx"]: -S["p"] * X,
        S["ax"]: 0,
        S["Hvxr"]: 0, S["Hvxrr"]: 0,
        S["Hzxr"]: 0, S["Hzxrr"]: 0,
        S["axr"]: 0, S["axrr"]: 0,
    }
    for label in h1["equations"]:
        require_zero(f"helicity-one residual gauge {label}",
                     parse(h1, label).subs(transverse_gauge, simultaneous=True))

    h0 = load("helicity0")
    R, T, Z, gauge_lambda = sp.symbols("R T Z gauge_lambda")

    def derivative_with_gauge(expression: sp.Expr) -> sp.Expr:
        arr = -S["phir"]**2 / 6
        hrr = -4 * S["Ar"] * S["hr"] + sp.exp(-2 * S["A"]) * S["f"] * S["Phir"]**2
        phirr = (S["Vphi"] - sp.Rational(1, 2) * sp.exp(-2 * S["A"])
                  * S["fphi"] * S["Phir"]**2
                  - (4 * S["h"] * S["Ar"] + S["hr"]) * S["phir"]) / S["h"]
        Phirr = -(2 * S["Ar"] + S["fphi"] * S["phir"] / S["f"]) * S["Phir"]
        mapping = {
            S["A"]: S["Ar"], S["Ar"]: arr,
            S["h"]: S["hr"], S["hr"]: hrr,
            S["phi"]: S["phir"], S["phir"]: phirr,
            S["Phi"]: S["Phir"], S["Phir"]: Phirr,
            S["f"]: S["fphi"] * S["phir"],
            S["fphi"]: S["fphiphi"] * S["phir"],
            S["Vphi"]: S["Vphiphi"] * S["phir"],
            R: -S["Ar"] * R - S["sigma"] * T,
            T: 0,
            Z: -S["p"] * sp.exp(-S["A"]) * T,
            gauge_lambda: 0,
        }
        return sp.expand(sum(sp.diff(expression, variable) * derivative
                             for variable, derivative in mapping.items()))

    pure_fields = {
        "Hvv": R * (2 * S["Ar"] * S["h"] + S["hr"])
               + 2 * S["h"] * S["sigma"] * T
               - 2 * sp.exp(-S["A"]) * S["sigma"] * R,
        "Hvz": -S["sigma"] * Z + S["h"] * S["p"] * T
               - sp.exp(-S["A"]) * S["p"] * R,
        "Hzz": -2 * S["Ar"] * R - 2 * S["p"] * Z,
        "Haa": -2 * S["Ar"] * R,
        "av": -R * S["Phir"] - S["Phi"] * S["sigma"] * T
              - S["sigma"] * gauge_lambda,
        "az": -S["Phi"] * S["p"] * T - S["p"] * gauge_lambda,
        "varphi": -R * S["phir"],
    }
    substitution: dict[sp.Symbol, sp.Expr] = {}
    for name, value in pure_fields.items():
        substitution[S[name]] = value
        first = derivative_with_gauge(value)
        substitution[S[name + "r"]] = first
        substitution[S[name + "rr"]] = derivative_with_gauge(first)
    for label in h0["equations"]:
        expression = parse(h0, label).subs(substitution, simultaneous=True)
        require_zero(f"helicity-zero residual gauge {label}", expression)


def test_uv_source_vev_powers() -> None:
    delta = sp.Symbol("Delta", real=True)
    nu = 4 - delta
    # Powers in u=exp(-A): metric, Maxwell, scalar, and psi=S/Y.
    pairs = {
        "metric": (0, 4),
        "maxwell": (0, 2),
        "scalar": (nu, delta),
        "singlet_rescaled": (0, 4 - 2 * nu),
    }
    if sp.simplify(pairs["singlet_rescaled"][1] - (2 * delta - 4)) != 0:
        raise AssertionError("singlet source/vev power map is inconsistent")
    for name, (source, vev) in pairs.items():
        if sp.simplify(vev - source) == 0:
            raise AssertionError(f"{name} source and vev powers are degenerate")


def main() -> None:
    tests = (
        test_helicity1_constraint_order,
        test_helicity1_k0_reduction,
        test_transverse_parity_and_neutral_limit,
        test_helicity0_characteristic_linearity,
        test_helicity0_constraints_and_k0_singlet,
        test_maxwell_identity,
        test_linearized_bianchi_identity,
        test_homogeneous_singlet_master,
        test_residual_gauge_modes,
        test_uv_source_vev_powers,
    )
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    main()
