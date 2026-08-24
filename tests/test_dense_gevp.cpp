#include "test_support.hpp"

#include "pha_qnm/spectral/dense_gevp.hpp"

#include <algorithm>
#include <complex>
#include <vector>

using namespace pha_qnm;

void test_dense_gevp() {
  ComplexMatrix m0(2, 2);
  ComplexMatrix m1(2, 2);
  m0(0, 0) = Complex{1.0, 0.0};
  m0(1, 1) = Complex{2.0, 0.0};
  m1(0, 0) = Complex{1.0, 0.0};
  m1(1, 1) = Complex{1.0, 0.0};

  const auto result = generalized_eigenvalues(m0, m1, false);
  std::vector<double> eigenvalues;
  eigenvalues.reserve(result.alpha.size());
  for (std::size_t i = 0; i < result.alpha.size(); ++i) {
    expect_true(std::abs(result.beta[i]) > 0.0, "finite test eigenvalue");
    const auto lambda = result.alpha[i] / result.beta[i];
    expect_near(lambda.imag(), 0.0, 1e-12, "real diagonal eigenvalue");
    eigenvalues.push_back(lambda.real());
  }
  std::sort(eigenvalues.begin(), eigenvalues.end());
  expect_near(eigenvalues.at(0), 1.0, 1e-12, "first generalized eigenvalue");
  expect_near(eigenvalues.at(1), 2.0, 1e-12, "second generalized eigenvalue");
}
