#include "pha_qnm/perturbations/decoupled.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pha_qnm {
namespace {

void validate(const RadialProfile& b, const RealMatrix& D1, const RealMatrix& D2) {
  const std::size_t n = b.A.size();
  if (n == 0 || D1.rows() != n || D1.cols() != n || D2.rows() != n || D2.cols() != n)
    throw std::invalid_argument("operator/profile dimension mismatch");
  const std::vector<std::size_t> sizes = {b.h.size(), b.phi.size(), b.Phi_prime.size(),
      b.A_prime.size(), b.h_prime.size(), b.phi_prime.size(), b.f.size(), b.f_prime.size()};
  if (!std::all_of(sizes.begin(), sizes.end(), [n](std::size_t x) { return x == n; }))
    throw std::invalid_argument("incomplete radial profile");
  if (b.horizon_index >= n || b.boundary_index >= n || b.horizon_index == b.boundary_index)
    throw std::invalid_argument("invalid boundary indices");
}

void impose_source_free_row(OperatorPencil& p, std::size_t row) {
  for (std::size_t j = 0; j < p.M0.cols(); ++j) {
    p.M0(row, j) = Complex{};
    p.M1(row, j) = Complex{};
  }
  p.M0(row, row) = Complex{1.0, 0.0};
}

} // namespace

OperatorPencil assemble_tensor_pencil(const RadialProfile& b, const RealMatrix& D1,
                                      const RealMatrix& D2, double k) {
  validate(b, D1, D2);
  const std::size_t n = b.A.size();
  OperatorPencil p{ComplexMatrix(n, n), ComplexMatrix(n, n)};
  const Complex minus_i{0.0, -1.0};
  for (std::size_t i = 0; i < n; ++i) {
    const double emA = std::exp(-b.A[i]);
    const double em2A = emA * emA;
    for (std::size_t j = 0; j < n; ++j) {
      p.M0(i, j) = b.h[i] * D2(i, j) +
                   (4.0 * b.h[i] * b.A_prime[i] + b.h_prime[i]) * D1(i, j);
      p.M1(i, j) = 2.0 * minus_i * emA * D1(i, j);
    }
    p.M0(i, i) -= em2A * k * k;
    p.M1(i, i) += 3.0 * minus_i * emA * b.A_prime[i];
  }
  impose_source_free_row(p, b.boundary_index);
  return p;
}

OperatorPencil assemble_vector_pencil(const RadialProfile& b, const RealMatrix& D1,
                                      const RealMatrix& D2) {
  validate(b, D1, D2);
  const std::size_t n = b.A.size();
  OperatorPencil p{ComplexMatrix(n, n), ComplexMatrix(n, n)};
  const Complex minus_i{0.0, -1.0};
  for (std::size_t i = 0; i < n; ++i) {
    if (!(b.f[i] > 0.0)) throw std::domain_error("non-positive coupling in vector pencil");
    const double emA = std::exp(-b.A[i]);
    const double ratio = b.f_prime[i] * b.phi_prime[i] / b.f[i];
    for (std::size_t j = 0; j < n; ++j) {
      p.M0(i, j) = b.h[i] * D2(i, j) +
                   (2.0 * b.h[i] * b.A_prime[i] + b.h_prime[i] + b.h[i] * ratio) * D1(i, j);
      p.M1(i, j) = 2.0 * minus_i * emA * D1(i, j);
    }
    p.M0(i, i) -= emA * emA * b.f[i] * b.Phi_prime[i] * b.Phi_prime[i];
    p.M1(i, i) += minus_i * emA * (b.A_prime[i] + ratio);
  }
  impose_source_free_row(p, b.boundary_index);
  return p;
}

OperatorPencil assemble_scalar_pencil(const RadialProfile& b, const RealMatrix& D1,
                                      const RealMatrix& D2) {
  validate(b, D1, D2);
  const std::size_t n = b.A.size();
  if (b.Y_prime_over_Y.size() != n)
    throw std::invalid_argument("scalar pencil requires Y'/Y profile");
  OperatorPencil p{ComplexMatrix(n, n), ComplexMatrix(n, n)};
  const Complex minus_i{0.0, -1.0};
  for (std::size_t i = 0; i < n; ++i) {
    if (std::abs(b.phi_prime[i]) < 1.0e-300)
      throw std::domain_error("scalar pencil encountered phi'=0");
    const double emA = std::exp(-b.A[i]);
    const double yratio = b.Y_prime_over_Y[i];
    for (std::size_t j = 0; j < n; ++j) {
      p.M0(i, j) = b.h[i] * D2(i, j) +
                   (4.0 * b.h[i] * b.A_prime[i] + b.h_prime[i] +
                    2.0 * b.h[i] * yratio) * D1(i, j);
      p.M1(i, j) = 2.0 * minus_i * emA * D1(i, j);
    }
    p.M0(i, i) += -b.h_prime[i] * yratio +
        emA * emA * b.Phi_prime[i] * b.Phi_prime[i] / b.phi_prime[i] *
        (3.0 * b.A_prime[i] * b.f_prime[i] - b.f[i] * b.phi_prime[i]);
    p.M1(i, i) += minus_i * emA * (3.0 * b.A_prime[i] + 2.0 * yratio);
  }
  impose_source_free_row(p, b.boundary_index);
  return p;
}

double normalized_pencil_residual(const OperatorPencil& p, Complex omega,
                                  const std::vector<Complex>& v) {
  if (p.M0.rows() != v.size() || p.M0.cols() != v.size())
    throw std::invalid_argument("pencil residual dimension mismatch");
  double residual2 = 0.0;
  double v2 = 0.0;
  double matrix_scale = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    Complex row{};
    double row_scale = 0.0;
    for (std::size_t j = 0; j < v.size(); ++j) {
      const Complex a = p.M0(i, j) + omega * p.M1(i, j);
      row += a * v[j];
      row_scale += std::abs(a);
    }
    residual2 += std::norm(row);
    v2 += std::norm(v[i]);
    matrix_scale = std::max(matrix_scale, row_scale);
  }
  return std::sqrt(residual2) /
         std::max(1.0e-300, matrix_scale * std::sqrt(v2));
}

} // namespace pha_qnm
