#include "test_support.hpp"

#include "pha_qnm/model/potentials.hpp"

#include <complex>

void test_potentials() {
  const pha_qnm::PhaParameters p;
  const double x = 1.0e-4;
  const double step = 1.0e-20;
  const std::complex<double> z{x, step};
  const double vpp_cs = std::imag(pha_qnm::potential_prime(z, p)) / step;
  const double fpp_cs = std::imag(pha_qnm::coupling_prime(z, p)) / step;
  expect_near(pha_qnm::potential_second(x, p), vpp_cs, 2.0e-12,
              "V second derivative disagrees with complex step");
  expect_near(pha_qnm::coupling_second(x, p), fpp_cs, 2.0e-8,
              "f second derivative disagrees with complex step");
  expect_near(pha_qnm::coupling(0.0, p), 1.0, 1.0e-15, "f(0) must be one");
  const auto d = pha_qnm::operator_dimension(p);
  expect_true(d.delta > 2.0 && d.delta < 4.0, "unexpected operator dimension");
}

