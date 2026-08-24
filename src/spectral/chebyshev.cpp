#include "pha_qnm/spectral/chebyshev.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace pha_qnm {

ChebyshevGrid chebyshev_lobatto(std::size_t intervals, double lower, double upper) {
  if (intervals < 2) throw std::invalid_argument("Chebyshev grid needs at least two intervals");
  if (!(upper > lower)) throw std::invalid_argument("invalid Chebyshev interval");
  const std::size_t count = intervals + 1;
  ChebyshevGrid out{{}, RealMatrix(count, count), RealMatrix(count, count)};
  out.nodes.resize(count);
  std::vector<double> barycentric(count);
  for (std::size_t j = 0; j < count; ++j) {
    const double x = -std::cos(std::numbers::pi * static_cast<double>(j) /
                               static_cast<double>(intervals));
    out.nodes[j] = 0.5 * ((upper - lower) * x + upper + lower);
    const double endpoint = (j == 0 || j == intervals) ? 0.5 : 1.0;
    barycentric[j] = (j % 2 == 0 ? 1.0 : -1.0) * endpoint;
  }
  for (std::size_t i = 0; i < count; ++i) {
    double row_sum = 0.0;
    for (std::size_t j = 0; j < count; ++j) {
      if (i == j) continue;
      const double value = barycentric[j] /
                           (barycentric[i] * (out.nodes[i] - out.nodes[j]));
      out.D1(i, j) = value;
      row_sum += value;
    }
    out.D1(i, i) = -row_sum;
  }
  for (std::size_t i = 0; i < count; ++i)
    for (std::size_t j = 0; j < count; ++j)
      for (std::size_t k = 0; k < count; ++k)
        out.D2(i, j) += out.D1(i, k) * out.D1(k, j);
  return out;
}

} // namespace pha_qnm

