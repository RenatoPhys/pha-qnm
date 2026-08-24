#pragma once

#include "pha_qnm/core/types.hpp"

#include <vector>

namespace pha_qnm {

struct GeneralizedEigenResult {
  std::vector<Complex> alpha;
  std::vector<Complex> beta;
  ComplexMatrix left_eigenvectors;
  ComplexMatrix right_eigenvectors;
};

GeneralizedEigenResult generalized_eigenvalues(ComplexMatrix M0, ComplexMatrix M1,
                                               bool compute_eigenvectors = true);

} // namespace pha_qnm

