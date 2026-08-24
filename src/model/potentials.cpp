#include "pha_qnm/model/potentials.hpp"

#include <cmath>
#include <stdexcept>

namespace pha_qnm {

OperatorDimension operator_dimension(const PhaParameters& p) {
  const double mass2 = -12.0 * p.gamma * p.gamma + 2.0 * p.b2;
  const double radicand = 4.0 + mass2;
  if (!(radicand >= 0.0)) throw std::domain_error("PHA scalar violates the real-dimension bound");
  const double delta = 2.0 + std::sqrt(radicand);
  return {mass2, delta, 4.0 - delta};
}

} // namespace pha_qnm

