#include "pha_qnm/background/horizon_series.hpp"

#include "pha_qnm/core/series.hpp"
#include "pha_qnm/model/potentials.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pha_qnm {
namespace {

using Jet = Series<6>;

struct Jets {
  Jet A;
  Jet h;
  Jet phi;
  Jet Phi;
};

Jet residual_A(const Jets& j) {
  const auto Ap = derivative(j.A);
  (void)Ap;
  const auto App = derivative(derivative(j.A));
  const auto phip = derivative(j.phi);
  return App + (phip * phip) / Jet(6.0);
}

Jet residual_h(const Jets& j, const PhaParameters& p) {
  const auto Ap = derivative(j.A);
  const auto hp = derivative(j.h);
  const auto Phip = derivative(j.Phi);
  return derivative(hp) + 4.0 * Ap * hp -
         math::exp(-2.0 * j.A) * coupling(j.phi, p) * Phip * Phip;
}

Jet residual_Phi(const Jets& j, const PhaParameters& p) {
  const auto Ap = derivative(j.A);
  const auto phip = derivative(j.phi);
  const auto Phip = derivative(j.Phi);
  return derivative(Phip) +
         (2.0 * Ap + coupling_prime(j.phi, p) * phip / coupling(j.phi, p)) * Phip;
}

Jet residual_phi_regular(const Jets& j, const PhaParameters& p) {
  const auto Ap = derivative(j.A);
  const auto hp = derivative(j.h);
  const auto phip = derivative(j.phi);
  const auto Phip = derivative(j.Phi);
  return j.h * derivative(phip) + (4.0 * j.h * Ap + hp) * phip -
         potential_prime(j.phi, p) +
         0.5 * math::exp(-2.0 * j.A) * coupling_prime(j.phi, p) * Phip * Phip;
}

template <class Residual>
double solve_linear_coefficient(Jet& target, std::size_t power, std::size_t equation_power,
                                Residual&& residual) {
  target[power] = 0.0;
  const double r0 = residual()[equation_power];
  target[power] = 1.0;
  const double slope = residual()[equation_power] - r0;
  if (!std::isfinite(slope) || std::abs(slope) < 1.0e-14)
    throw std::runtime_error("singular horizon recurrence");
  target[power] = -r0 / slope;
  return target[power];
}

} // namespace

double maximum_horizon_electric_field(double phi0, const PhaParameters& p) {
  const double ratio = -2.0 * potential(phi0, p) / coupling(phi0, p);
  if (!(ratio > 0.0)) throw std::domain_error("no regular charged horizon at phi0");
  return std::sqrt(ratio);
}

HorizonSeries make_horizon_series(double phi0, double Phi1, const PhaParameters& p) {
  if (!(phi0 >= 0.0) || !std::isfinite(phi0) || !std::isfinite(Phi1))
    throw std::invalid_argument("non-finite horizon datum");
  if (std::abs(Phi1) >= maximum_horizon_electric_field(phi0, p))
    throw std::domain_error("horizon electric field is extremal or super-extremal");

  Jets j;
  j.A[0] = 0.0;
  j.h[0] = 0.0;
  j.h[1] = 1.0;
  j.phi[0] = phi0;
  j.Phi[0] = 0.0;
  j.Phi[1] = Phi1;
  j.A[1] = -(2.0 * potential(phi0, p) + coupling(phi0, p) * Phi1 * Phi1) / 6.0;
  j.phi[1] = potential_prime(phi0, p) - 0.5 * coupling_prime(phi0, p) * Phi1 * Phi1;

  for (std::size_t n = 2; n <= 4; ++n) {
    solve_linear_coefficient(j.A, n, n - 2, [&] { return residual_A(j); });
    solve_linear_coefficient(j.h, n, n - 2, [&] { return residual_h(j, p); });
    solve_linear_coefficient(j.Phi, n, n - 2, [&] { return residual_Phi(j, p); });
    solve_linear_coefficient(j.phi, n, n - 1, [&] { return residual_phi_regular(j, p); });
  }

  HorizonSeries out;
  std::copy_n(j.A.c.begin(), 5, out.A.begin());
  std::copy_n(j.h.c.begin(), 5, out.h.begin());
  std::copy_n(j.phi.c.begin(), 5, out.phi.begin());
  std::copy_n(j.Phi.c.begin(), 5, out.Phi.begin());
  return out;
}

std::array<double, 8> HorizonSeries::state(double r) const {
  auto value = [r](const std::array<double, 5>& c) {
    double y = c.back();
    for (std::size_t i = c.size() - 1; i-- > 0;) y = y * r + c[i];
    return y;
  };
  auto slope = [r](const std::array<double, 5>& c) {
    double y = 4.0 * c[4];
    for (int i = 3; i >= 1; --i) y = y * r + static_cast<double>(i) * c[static_cast<std::size_t>(i)];
    return y;
  };
  return {value(A), value(h), value(phi), value(Phi),
          slope(A), slope(h), slope(phi), slope(Phi)};
}

} // namespace pha_qnm
