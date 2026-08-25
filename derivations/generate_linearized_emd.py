"""Generate auditable linearized EMD equations in ingoing EF coordinates.

The generator uses first-order dual expressions, so every tensor operation is
truncated exactly at linear order in the perturbations.  It emits the raw
equations, local radial-operator coefficients, and content hashes.  Generated
files are evidence; production kernels are promoted only after the identity
and limiting-case tests pass.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "derivations" / "generated"
CPP_OUTPUT = ROOT / "include" / "pha_qnm" / "perturbations" / "generated"

v, r, x, y, z = sp.symbols("v r x y z", real=True)
coords = (v, r, x, y, z)
sv, pz = sp.symbols("sigma p", complex=True)
wave = sp.exp(sv * v + pz * z)

A = sp.Function("A")(r)
h = sp.Function("h")(r)
phi = sp.Function("phi")(r)
Phi = sp.Function("Phi")(r)
f0 = sp.Function("f")(r)
fp = sp.Function("f_phi")(r)
fpp = sp.Function("f_phiphi")(r)
V0 = sp.Function("V")(r)
Vp = sp.Function("V_phi")(r)
Vpp = sp.Function("V_phiphi")(r)


@dataclass(frozen=True)
class Dual:
    background: sp.Expr
    linear: sp.Expr = sp.S.Zero

    @staticmethod
    def lift(value: object) -> "Dual":
        return value if isinstance(value, Dual) else Dual(sp.sympify(value))

    def __add__(self, other: object) -> "Dual":
        other = Dual.lift(other)
        return Dual(self.background + other.background, self.linear + other.linear)

    __radd__ = __add__

    def __neg__(self) -> "Dual":
        return Dual(-self.background, -self.linear)

    def __sub__(self, other: object) -> "Dual":
        return self + (-Dual.lift(other))

    def __rsub__(self, other: object) -> "Dual":
        return Dual.lift(other) - self

    def __mul__(self, other: object) -> "Dual":
        other = Dual.lift(other)
        return Dual(self.background * other.background,
                    self.linear * other.background + self.background * other.linear)

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "Dual":
        other = Dual.lift(other)
        return Dual(self.background / other.background,
                    self.linear / other.background
                    - self.background * other.linear / other.background**2)

    def __rtruediv__(self, other: object) -> "Dual":
        return Dual.lift(other) / self


def derivative(value: Dual, coordinate: sp.Symbol) -> Dual:
    return Dual(sp.diff(value.background, coordinate), sp.diff(value.linear, coordinate))


def background_metric() -> tuple[list[list[sp.Expr]], list[list[sp.Expr]]]:
    eA = sp.exp(A)
    e2A = eA**2
    g = [[sp.S.Zero for _ in range(5)] for _ in range(5)]
    g[0][0] = -e2A * h
    g[0][1] = g[1][0] = eA
    for index in (2, 3, 4):
        g[index][index] = e2A
    inverse = sp.Matrix(g).inv().tolist()
    return g, inverse


G0, G0INV = background_metric()


def make_metric(perturbations: dict[tuple[int, int], sp.Expr]) -> tuple[list[list[Dual]], list[list[Dual]]]:
    e2A = sp.exp(2 * A)
    dg = [[sp.S.Zero for _ in range(5)] for _ in range(5)]
    for (i, j), amplitude in perturbations.items():
        dg[i][j] = e2A * amplitude * wave
        dg[j][i] = dg[i][j]
    metric = [[Dual(G0[i][j], dg[i][j]) for j in range(5)] for i in range(5)]
    inverse = []
    for i in range(5):
        row = []
        for j in range(5):
            correction = -sum(G0INV[i][a] * dg[a][b] * G0INV[b][j]
                              for a in range(5) for b in range(5))
            row.append(Dual(G0INV[i][j], correction))
        inverse.append(row)
    return metric, inverse


def christoffel(metric: list[list[Dual]], inverse: list[list[Dual]]) -> list[list[list[Dual]]]:
    return [[[sum(inverse[a][d] * (derivative(metric[d][c], coords[b])
                                    + derivative(metric[d][b], coords[c])
                                    - derivative(metric[b][c], coords[d]))
                 for d in range(5)) / 2
              for c in range(5)] for b in range(5)] for a in range(5)]


def ricci(gamma: list[list[list[Dual]]]) -> list[list[Dual]]:
    result = [[Dual(0) for _ in range(5)] for _ in range(5)]
    for m in range(5):
        for n in range(5):
            value = Dual(0)
            for a in range(5):
                value += derivative(gamma[a][m][n], coords[a])
                value -= derivative(gamma[a][m][a], coords[n])
                for b in range(5):
                    value += gamma[a][a][b] * gamma[b][m][n]
                    value -= gamma[a][n][b] * gamma[b][m][a]
            result[m][n] = value
    return result


def field_strength(gauge_linear: dict[int, sp.Expr]) -> list[list[Dual]]:
    gauge = [Dual(Phi if i == 0 else 0,
                  gauge_linear.get(i, sp.S.Zero) * wave) for i in range(5)]
    return [[derivative(gauge[n], coords[m]) - derivative(gauge[m], coords[n])
             for n in range(5)] for m in range(5)]


def raised_field_strength(F: list[list[Dual]], inverse: list[list[Dual]]) -> list[list[Dual]]:
    return [[sum(inverse[m][a] * inverse[n][b] * F[a][b]
                 for a in range(5) for b in range(5))
             for n in range(5)] for m in range(5)]


def equations(metric_perturbations: dict[tuple[int, int], sp.Expr],
              gauge_perturbations: dict[int, sp.Expr], scalar: sp.Expr) -> dict[str, object]:
    metric, inverse = make_metric(metric_perturbations)
    gamma = christoffel(metric, inverse)
    R = ricci(gamma)
    F = field_strength(gauge_perturbations)
    Fup = raised_field_strength(F, inverse)
    scalar_field = Dual(phi, scalar * wave)
    coupling = Dual(f0, fp * scalar * wave)
    potential = Dual(V0, Vp * scalar * wave)
    potential_phi = Dual(Vp, Vpp * scalar * wave)
    coupling_phi = Dual(fp, fpp * scalar * wave)

    F2 = sum(F[m][n] * Fup[m][n] for m in range(5) for n in range(5))
    einstein = [[Dual(0) for _ in range(5)] for _ in range(5)]
    for m in range(5):
        for n in range(5):
            scalar_term = derivative(scalar_field, coords[m]) * derivative(scalar_field, coords[n]) / 2
            maxwell_mn = sum(F[m][a] * inverse[a][b] * F[n][b]
                             for a in range(5) for b in range(5))
            einstein[m][n] = (R[m][n] - scalar_term - potential * metric[m][n] / 3
                              - coupling * (maxwell_mn - metric[m][n] * F2 / 6) / 2)

    trace_metric_perturbation = sum(G0INV[m][n] * metric[m][n].linear
                                    for m in range(5) for n in range(5))
    sqrt_g = Dual(sp.exp(4 * A), sp.exp(4 * A) * trace_metric_perturbation / 2)
    maxwell = []
    for n in range(5):
        maxwell.append(sum(derivative(sqrt_g * coupling * Fup[m][n], coords[m])
                           for m in range(5)) / sqrt_g)

    box_phi = Dual(0)
    for m in range(5):
        for n in range(5):
            second = derivative(derivative(scalar_field, coords[n]), coords[m])
            connection = sum(gamma[a][m][n] * derivative(scalar_field, coords[a])
                             for a in range(5))
            box_phi += inverse[m][n] * (second - connection)
    scalar_eq = box_phi - potential_phi - coupling_phi * F2 / 4
    return {"einstein": einstein, "maxwell": maxwell, "scalar": scalar_eq,
            "internal": {"metric": metric, "inverse": inverse, "gamma": gamma,
                         "ricci": R, "field_strength": F, "sqrt_g": sqrt_g}}


def normalized(expression: Dual) -> sp.Expr:
    result = sp.expand(expression.linear / wave).xreplace({wave: sp.S.One})
    return sp.factor_terms(result)


def symbol_substitutions(fields: Iterable[sp.Expr]) -> tuple[dict[sp.Expr, sp.Symbol], list[sp.Symbol]]:
    substitutions: dict[sp.Expr, sp.Symbol] = {}
    for name, function in (("A", A), ("h", h), ("phi", phi), ("Phi", Phi),
                           ("f", f0), ("fphi", fp), ("fphiphi", fpp),
                           ("V", V0), ("Vphi", Vp), ("Vphiphi", Vpp)):
        substitutions[function] = sp.Symbol(name, real=True)
        for order in (1, 2, 3):
            substitutions[sp.diff(function, r, order)] = sp.Symbol(name + "r" * order, real=True)
    basis: list[sp.Symbol] = []
    for field in fields:
        name = str(field.func)
        for order in (0, 1, 2):
            source = field if order == 0 else sp.diff(field, r, order)
            target = sp.Symbol(name + "r" * order, complex=True)
            substitutions[source] = target
            basis.append(target)
    # Replace derivatives before their underlying functions.
    ordered = dict(sorted(substitutions.items(), key=lambda item: len(str(item[0])), reverse=True))
    return ordered, basis


def encode_equation(expression: Dual, substitutions: dict[sp.Expr, sp.Symbol],
                    basis: list[sp.Symbol]) -> dict[str, str]:
    raw = normalized(expression).subs(substitutions, simultaneous=True)
    expanded = sp.expand(raw)
    coefficient_pairs = [(item, sp.factor_terms(expanded.coeff(item))) for item in basis
                         if expanded.coeff(item) != 0]
    remainder = sp.simplify(expanded - sum(value * item for item, value in coefficient_pairs))
    if remainder != 0:
        raise RuntimeError(f"nonlinear or uncollected remainder: {remainder}")
    return {"expression": sp.sstr(raw),
            "coefficients": {str(item): sp.sstr(value) for item, value in coefficient_pairs}}


def write_cpp_types() -> Path:
    path = CPP_OUTPUT / "emd_kernel_types.generated.hpp"
    path.write_text("""// Generated by derivations/generate_linearized_emd.py; do not edit.
#pragma once

#include "pha_qnm/core/types.hpp"

#include <array>
#include <string_view>

namespace pha_qnm::generated {

struct EmdPoint {
  double A{}, Ar{}, Arr{};
  double h{}, hr{}, hrr{};
  double phi{}, phir{}, phirr{};
  double Phi{}, Phir{}, Phirr{};
  double f{}, fr{}, fphi{}, fphir{}, fphiphi{};
  double V{}, Vphi{}, Vphiphi{};
};

template <std::size_t EquationCount, std::size_t FieldCount>
struct LocalRadialCoefficients {
  std::array<Complex, EquationCount * FieldCount * 3> values{};

  [[nodiscard]] constexpr Complex operator()(std::size_t equation,
                                              std::size_t field,
                                              std::size_t derivative) const {
    return values[(equation * FieldCount + field) * 3 + derivative];
  }
};

} // namespace pha_qnm::generated
""", encoding="utf-8")
    return path


def cpp_expression(text: str) -> str:
    expression = sp.sympify(text)
    background_names = ("A", "Ar", "Arr", "h", "hr", "hrr", "phi", "phir", "phirr",
                        "Phi", "Phir", "Phirr", "f", "fr", "fphi", "fphir", "fphiphi",
                        "V", "Vphi", "Vphiphi")
    replacements = {sp.Symbol(name): sp.Symbol(f"b_{name}")
                    for name in background_names}
    expression = expression.subs(replacements, simultaneous=True)
    code = sp.cxxcode(expression, standard="c++17")
    # MSVC cannot deduce std::complex<double> multiplication from an integer
    # literal. Promote direct integer factors of the complex Fourier symbols.
    code = re.sub(r"(?<![\w.])(\d+)\*(sigma|p)\b", r"\1.0*\2", code)
    code = re.sub(r"\b(sigma|p)\*(\d+)(?![\w.])", r"\1*\2.0", code)
    for name in sorted(background_names, key=len, reverse=True):
        code = code.replace(f"b_{name}", f"b.{name}")
    return code


def derivative_index(field: str, coefficient: str) -> int | None:
    for index, suffix in ((0, ""), (1, "r"), (2, "rr")):
        if coefficient == field + suffix:
            return index
    return None


def write_cpp_sector(name: str, fields: list[sp.Expr],
                     equations: dict[str, dict[str, object]]) -> Path:
    field_names = [str(field.func) for field in fields]
    equation_names = list(equations)
    path = CPP_OUTPUT / f"{name}_kernels.generated.hpp"
    function = f"{name}_coefficients"
    lines = [
        "// Generated by derivations/generate_linearized_emd.py; do not edit.",
        "#pragma once", "", '#include "emd_kernel_types.generated.hpp"', "",
        "#include <array>", "#include <cmath>", "#include <string_view>", "",
        "namespace pha_qnm::generated {", "",
        f"inline constexpr std::array<std::string_view, {len(equation_names)}> {name}_equations = {{",
    ]
    lines.extend(f'    "{item}",' for item in equation_names)
    lines.extend(["};", "",
                  f"inline constexpr std::array<std::string_view, {len(field_names)}> {name}_fields = {{"])
    lines.extend(f'    "{item}",' for item in field_names)
    lines.extend(["};", "",
                  f"inline LocalRadialCoefficients<{len(equation_names)}, {len(field_names)}> {function}(",
                  "    const EmdPoint& b, Complex sigma, Complex p) {",
                  f"  LocalRadialCoefficients<{len(equation_names)}, {len(field_names)}> out;",
                  "  using std::exp;", "  using std::pow;"])
    for equation_index, label in enumerate(equation_names):
        coefficients = equations[label]["coefficients"]
        for coefficient_name, value in coefficients.items():
            for field_index, field in enumerate(field_names):
                order = derivative_index(field, coefficient_name)
                if order is not None:
                    flat = (equation_index * len(field_names) + field_index) * 3 + order
                    lines.append(f"  out.values[{flat}] = {cpp_expression(value)};")
                    break
            else:
                raise RuntimeError(f"unmapped coefficient {coefficient_name}")
    lines.extend(["  return out;", "}", "", "} // namespace pha_qnm::generated", ""])
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def sector_data(name: str) -> tuple[dict[str, object], list[sp.Expr], dict[str, Dual]]:
    if name == "helicity1":
        Hvx, Hzx, ax = (sp.Function(item)(r) for item in ("Hvx", "Hzx", "ax"))
        system = equations({(0, 2): Hvx, (4, 2): Hzx}, {2: ax}, sp.S.Zero)
        selected = {
            "evolution_E_rx": system["einstein"][1][2],
            "E_vx": system["einstein"][0][2],
            "evolution_E_zx": system["einstein"][4][2],
            "evolution_M_x": system["maxwell"][2],
            "constraint_Er_x": (sp.exp(-A) * system["einstein"][0][2]
                                 + h * system["einstein"][1][2]),
        }
        return system, [Hvx, Hzx, ax], selected
    if name == "helicity0":
        names = ("Hvv", "Hvz", "Hzz", "Haa", "av", "az", "varphi")
        Hvv, Hvz, Hzz, Haa, av, az, varphi = (sp.Function(item)(r) for item in names)
        system = equations({(0, 0): Hvv, (0, 4): Hvz, (4, 4): Hzz,
                            (2, 2): Haa, (3, 3): Haa}, {0: av, 4: az}, varphi)
        selected = {
            "evolution_E_rr": system["einstein"][1][1],
            "evolution_E_rv": system["einstein"][1][0],
            "evolution_E_rz": system["einstein"][1][4],
            "evolution_E_aa": (system["einstein"][2][2] + system["einstein"][3][3]) / 2,
            "evolution_M_v": system["maxwell"][0],
            "evolution_M_z": system["maxwell"][4],
            "evolution_scalar": system["scalar"],
            "diagnostic_E_vv": system["einstein"][0][0],
            "diagnostic_E_vz": system["einstein"][0][4],
            "diagnostic_E_zz": system["einstein"][4][4],
            "constraint_Er_r": (sp.exp(-A) * system["einstein"][0][1]
                                 + h * system["einstein"][1][1]),
            "constraint_Er_v": (sp.exp(-A) * system["einstein"][0][0]
                                 + h * system["einstein"][1][0]),
            "constraint_Er_z": (sp.exp(-A) * system["einstein"][0][4]
                                 + h * system["einstein"][1][4]),
            "constraint_M_r": system["maxwell"][1],
        }
        return system, [Hvv, Hvz, Hzz, Haa, av, az, varphi], selected
    raise ValueError(name)


def write_sector(name: str) -> dict[str, str]:
    _, fields, selected = sector_data(name)
    substitutions, basis = symbol_substitutions(fields)
    encoded = {label: encode_equation(value, substitutions, basis)
               for label, value in selected.items()}
    data = {
        "coordinate_order": ["v", "r", "x", "y", "z"],
        "fourier_replacements": {"sigma": "-i omega", "p": "i k"},
        "metric_normalization": "delta g_MN = exp(2 A) H_MN exp(sigma v + p z)",
        "fields": [str(field.func) for field in fields],
        "equations": encoded,
    }
    path = OUTPUT / f"{name}_equations.json"
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    cpp_path = write_cpp_sector(name, fields, encoded)
    return {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest(),
            cpp_path.relative_to(ROOT).as_posix(): hashlib.sha256(cpp_path.read_bytes()).hexdigest()}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    CPP_OUTPUT.mkdir(parents=True, exist_ok=True)
    hashes: dict[str, str] = {}
    type_path = write_cpp_types()
    hashes[type_path.relative_to(ROOT).as_posix()] = hashlib.sha256(type_path.read_bytes()).hexdigest()
    for sector in ("helicity1", "helicity0"):
        hashes.update(write_sector(sector))
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    validator = Path(__file__).with_name("validate_linearized_emd.py")
    hashes[validator.name] = hashlib.sha256(validator.read_bytes()).hexdigest()
    manifest = {"generator": Path(__file__).name, "sympy": sp.__version__, "sha256": hashes}
    (OUTPUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
