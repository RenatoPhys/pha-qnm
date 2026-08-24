#include "test_support.hpp"

#include "pha_qnm/perturbations/decoupled.hpp"
#include "pha_qnm/spectral/chebyshev.hpp"

#include <cmath>

void test_operators() {
  const auto grid = pha_qnm::chebyshev_lobatto(10, 0.0, 1.0);
  const std::size_t n = grid.nodes.size();
  pha_qnm::RadialProfile b;
  b.A.resize(n); b.h.resize(n); b.phi.resize(n); b.Phi_prime.resize(n);
  b.A_prime.resize(n); b.h_prime.resize(n); b.phi_prime.resize(n);
  b.f.resize(n); b.f_prime.resize(n); b.Y_prime_over_Y.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double r = grid.nodes[i];
    b.A[i] = r; b.h[i] = r; b.phi[i] = 0.1 * r; b.Phi_prime[i] = 0.02;
    b.A_prime[i] = 1.0; b.h_prime[i] = 1.0; b.phi_prime[i] = 0.1;
    b.f[i] = 1.0; b.f_prime[i] = 0.0; b.Y_prime_over_Y[i] = 0.2;
  }
  b.horizon_index = 0;
  b.boundary_index = n - 1;
  const auto k0 = pha_qnm::assemble_tensor_pencil(b, grid.D1, grid.D2, 0.0);
  const auto tiny = pha_qnm::assemble_tensor_pencil(b, grid.D1, grid.D2, 1.0e-8);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      expect_true(std::abs(k0.M1(i, j) - tiny.M1(i, j)) < 1.0e-15,
                  "tensor k->0 changed M1");
      expect_true(std::abs(k0.M0(i, j) - tiny.M0(i, j)) < 1.0e-14,
                  "tensor k->0 limit failed");
    }
  expect_near(k0.M0(b.boundary_index, b.boundary_index).real(), 1.0, 0.0,
              "UV source row diagonal");
  const auto scalar = pha_qnm::assemble_scalar_pencil(b, grid.D1, grid.D2);
  expect_near(scalar.M0(b.boundary_index, b.boundary_index).real(), 1.0, 0.0,
              "scalar UV source row diagonal");
}
