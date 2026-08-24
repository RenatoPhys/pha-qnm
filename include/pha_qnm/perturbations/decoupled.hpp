#pragma once

#include "pha_qnm/core/types.hpp"

#include <vector>

namespace pha_qnm {

struct RadialProfile {
  std::vector<double> A;
  std::vector<double> h;
  std::vector<double> phi;
  std::vector<double> Phi_prime;
  std::vector<double> A_prime;
  std::vector<double> h_prime;
  std::vector<double> phi_prime;
  std::vector<double> f;
  std::vector<double> f_prime;
  std::vector<double> Y_prime_over_Y;
  std::size_t horizon_index{};
  std::size_t boundary_index{};
};

struct OperatorPencil {
  ComplexMatrix M0;
  ComplexMatrix M1;
};

OperatorPencil assemble_tensor_pencil(const RadialProfile& background,
                                      const RealMatrix& D1, const RealMatrix& D2,
                                      double momentum);
OperatorPencil assemble_vector_pencil(const RadialProfile& background,
                                      const RealMatrix& D1, const RealMatrix& D2);
OperatorPencil assemble_scalar_pencil(const RadialProfile& background,
                                      const RealMatrix& D1, const RealMatrix& D2);
double normalized_pencil_residual(const OperatorPencil& pencil, Complex omega,
                                  const std::vector<Complex>& eigenvector);

} // namespace pha_qnm
