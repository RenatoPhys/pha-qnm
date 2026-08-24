#include "test_support.hpp"

#include "pha_qnm/background/horizon_series.hpp"
#include "pha_qnm/model/potentials.hpp"

#include <cmath>

void test_horizon_series() {
  const pha_qnm::PhaParameters p;
  const double phi0 = 1.0;
  const double Phi1 = 0.05;
  const auto s = pha_qnm::make_horizon_series(phi0, Phi1, p);
  const double A1 = -(2.0 * pha_qnm::potential(phi0, p) +
                      pha_qnm::coupling(phi0, p) * Phi1 * Phi1) / 6.0;
  const double phi1 = pha_qnm::potential_prime(phi0, p) -
                      0.5 * pha_qnm::coupling_prime(phi0, p) * Phi1 * Phi1;
  expect_near(s.A[1], A1, 1.0e-14, "A1");
  expect_near(s.phi[1], phi1, 1.0e-14, "phi1");
  expect_near(s.A[2], -phi1 * phi1 / 12.0, 1.0e-12, "A2");
  expect_near(s.h[2], 0.5 * (pha_qnm::coupling(phi0, p) * Phi1 * Phi1 - 4.0 * A1),
              1.0e-12, "h2");
  for (double value : s.A) expect_true(std::isfinite(value), "non-finite A horizon coefficient");
  for (double value : s.phi) expect_true(std::isfinite(value), "non-finite phi horizon coefficient");
}

