#include "test_support.hpp"

#include "pha_qnm/spectral/chebyshev.hpp"

#include <algorithm>
#include <cmath>

void test_chebyshev() {
  const auto grid = pha_qnm::chebyshev_lobatto(18, -0.3, 1.7);
  double e1 = 0.0;
  double e2 = 0.0;
  for (std::size_t i = 0; i < grid.nodes.size(); ++i) {
    double d1{};
    double d2{};
    for (std::size_t j = 0; j < grid.nodes.size(); ++j) {
      const double f = std::pow(grid.nodes[j], 6);
      d1 += grid.D1(i, j) * f;
      d2 += grid.D2(i, j) * f;
    }
    e1 = std::max(e1, std::abs(d1 - 6.0 * std::pow(grid.nodes[i], 5)));
    e2 = std::max(e2, std::abs(d2 - 30.0 * std::pow(grid.nodes[i], 4)));
  }
  expect_true(e1 < 1.0e-9, "Chebyshev D1 polynomial error");
  expect_true(e2 < 1.0e-7, "Chebyshev D2 polynomial error");
}

