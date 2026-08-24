#pragma once

#include "pha_qnm/core/types.hpp"

#include <vector>

namespace pha_qnm {

struct ChebyshevGrid {
  std::vector<double> nodes;
  RealMatrix D1;
  RealMatrix D2;
};

ChebyshevGrid chebyshev_lobatto(std::size_t intervals, double lower, double upper);

} // namespace pha_qnm

