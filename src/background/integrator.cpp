#include "pha_qnm/background/integrator.hpp"

#include "pha_qnm/background/horizon_series.hpp"
#include "pha_qnm/model/potentials.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace pha_qnm {

BackgroundState background_rhs(double, const BackgroundState& y, const PhaParameters& p) {
  const auto [A, h, phi, Phi, Ap, hp, phip, Phip] = y;
  (void)Phi;
  if (!(h > 0.0)) throw std::domain_error("background RHS requires h>0");
  const double f = coupling(phi, p);
  const double fp = coupling_prime(phi, p);
  if (!(f > 0.0)) throw std::domain_error("non-positive Maxwell coupling");
  const double em2A = std::exp(-2.0 * A);
  return {Ap,
          hp,
          phip,
          Phip,
          -phip * phip / 6.0,
          -4.0 * Ap * hp + em2A * f * Phip * Phip,
          -(4.0 * Ap + hp / h) * phip +
              (potential_prime(phi, p) - 0.5 * em2A * fp * Phip * Phip) / h,
          -(2.0 * Ap + fp * phip / f) * Phip};
}

double einstein_constraint(const BackgroundState& y, const PhaParameters& p) {
  const auto [A, h, phi, Phi, Ap, hp, phip, Phip] = y;
  (void)Phi;
  return h * (24.0 * Ap * Ap - phip * phip) + 6.0 * Ap * hp +
         2.0 * potential(phi, p) + std::exp(-2.0 * A) * coupling(phi, p) * Phip * Phip;
}

double gauss_charge(const BackgroundState& y, const PhaParameters& p) {
  return std::exp(2.0 * y[0]) * coupling(y[2], p) * y[7];
}

namespace {

BackgroundState add_scaled(const BackgroundState& y, const BackgroundState& k, double scale) {
  BackgroundState out{};
  for (std::size_t i = 0; i < out.size(); ++i) out[i] = y[i] + scale * k[i];
  return out;
}

BackgroundState rk4_step(double r, const BackgroundState& y, double step,
                         const PhaParameters& p) {
  const auto k1 = background_rhs(r, y, p);
  const auto k2 = background_rhs(r + 0.5 * step, add_scaled(y, k1, 0.5 * step), p);
  const auto k3 = background_rhs(r + 0.5 * step, add_scaled(y, k2, 0.5 * step), p);
  const auto k4 = background_rhs(r + step, add_scaled(y, k3, step), p);
  BackgroundState out{};
  for (std::size_t i = 0; i < out.size(); ++i)
    out[i] = y[i] + step * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
  return out;
}

bool finite(const BackgroundState& y) {
  return std::all_of(y.begin(), y.end(), [](double x) { return std::isfinite(x); });
}

} // namespace

BackgroundSolution integrate_background(double phi0, double Phi1, const PhaParameters& p,
                                        const IntegratorOptions& opt) {
  if (!(opt.epsilon_horizon > 0.0) || !(opt.r_max > opt.epsilon_horizon))
    throw std::invalid_argument("invalid radial interval");
  const auto series = make_horizon_series(phi0, Phi1, p);
  double r = opt.epsilon_horizon;
  BackgroundState y = series.state(r);
  double step = std::min(opt.initial_step, opt.r_max - r);
  const double q0 = gauss_charge(y, p);

  BackgroundSolution out;
  out.points.push_back({r, y, einstein_constraint(y, p), q0});

  while (r < opt.r_max && out.accepted_steps + out.rejected_steps < opt.max_steps) {
    step = std::min({step, opt.max_step, opt.r_max - r});
    const auto full = rk4_step(r, y, step, p);
    const auto half1 = rk4_step(r, y, 0.5 * step, p);
    const auto half2 = rk4_step(r + 0.5 * step, half1, 0.5 * step, p);
    double error = 0.0;
    for (std::size_t i = 0; i < y.size(); ++i) {
      const double scale = opt.absolute_tolerance +
                           opt.relative_tolerance * std::max(std::abs(y[i]), std::abs(half2[i]));
      error = std::max(error, std::abs(half2[i] - full[i]) / (15.0 * scale));
    }
    if (!finite(half2)) throw std::runtime_error("non-finite background state");
    if (error <= 1.0) {
      for (std::size_t i = 0; i < y.size(); ++i)
        y[i] = half2[i] + (half2[i] - full[i]) / 15.0;
      r += step;
      if (!(y[1] > 0.0)) throw std::runtime_error("blackening function left physical branch");
      const double constraint = einstein_constraint(y, p);
      const double q = gauss_charge(y, p);
      const double constraint_scale = 1.0 + std::abs(2.0 * potential(y[2], p)) +
                                      std::abs(6.0 * y[4] * y[5]);
      out.maximum_normalized_constraint =
          std::max(out.maximum_normalized_constraint, std::abs(constraint) / constraint_scale);
      out.maximum_relative_charge_drift =
          std::max(out.maximum_relative_charge_drift,
                   std::abs(q - q0) / std::max(1.0e-300, std::abs(q0)));
      out.points.push_back({r, y, constraint, q});
      ++out.accepted_steps;
    } else {
      ++out.rejected_steps;
    }
    const double factor = error == 0.0 ? 2.0 : std::clamp(0.9 * std::pow(error, -0.2), 0.2, 2.0);
    step *= factor;
    if (step < 10.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, r))
      throw std::runtime_error("background step size underflow");
  }
  if (r < opt.r_max) throw std::runtime_error("background exceeded maximum step count");
  return out;
}

} // namespace pha_qnm

