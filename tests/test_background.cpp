#include "test_support.hpp"

#include "pha_qnm/background/integrator.hpp"

void test_background() {
  pha_qnm::IntegratorOptions options;
  options.r_max = 0.25;
  options.initial_step = 1.0e-5;
  options.max_step = 1.0e-3;
  options.absolute_tolerance = 1.0e-12;
  options.relative_tolerance = 1.0e-11;
  const auto solution = pha_qnm::integrate_background(1.0, 0.05,
                                                       pha_qnm::PhaParameters{}, options);
  expect_true(solution.points.size() > 20, "background stored too few points");
  expect_true(solution.maximum_normalized_constraint < 1.0e-7,
              "background constraint did not converge");
  expect_true(solution.maximum_relative_charge_drift < 1.0e-8,
              "Gauss charge drift did not converge");
}

