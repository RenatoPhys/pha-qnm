#pragma once

#include "pha_qnm/core/series.hpp"
#include "pha_qnm/model/pha_parameters.hpp"

#include <cmath>

namespace pha_qnm {

template <class T>
T potential(const T& phi, const PhaParameters& p) {
  const T phi2 = phi * phi;
  return -12.0 * math::cosh(p.gamma * phi) + p.b2 * phi2 +
         p.b4 * phi2 * phi2 + p.b6 * phi2 * phi2 * phi2;
}

template <class T>
T potential_prime(const T& phi, const PhaParameters& p) {
  const T phi2 = phi * phi;
  return -12.0 * p.gamma * math::sinh(p.gamma * phi) + 2.0 * p.b2 * phi +
         4.0 * p.b4 * phi2 * phi + 6.0 * p.b6 * phi2 * phi2 * phi;
}

template <class T>
T potential_second(const T& phi, const PhaParameters& p) {
  const T phi2 = phi * phi;
  return -12.0 * p.gamma * p.gamma * math::cosh(p.gamma * phi) + 2.0 * p.b2 +
         12.0 * p.b4 * phi2 + 30.0 * p.b6 * phi2 * phi2;
}

template <class T>
T coupling_argument(const T& phi, const PhaParameters& p) {
  return p.c1 * phi + p.c2 * phi * phi + p.c3 * phi * phi * phi;
}

template <class T>
T coupling(const T& phi, const PhaParameters& p) {
  return (math::sech(coupling_argument(phi, p)) + p.d1 * math::sech(p.d2 * phi)) /
         (1.0 + p.d1);
}

template <class T>
T coupling_prime(const T& phi, const PhaParameters& p) {
  const T g = coupling_argument(phi, p);
  const T gp = p.c1 + 2.0 * p.c2 * phi + 3.0 * p.c3 * phi * phi;
  return (-math::sech(g) * math::tanh(g) * gp -
          p.d1 * p.d2 * math::sech(p.d2 * phi) * math::tanh(p.d2 * phi)) /
         (1.0 + p.d1);
}

template <class T>
T coupling_second(const T& phi, const PhaParameters& p) {
  const T g = coupling_argument(phi, p);
  const T gp = p.c1 + 2.0 * p.c2 * phi + 3.0 * p.c3 * phi * phi;
  const T gpp = 2.0 * p.c2 + 6.0 * p.c3 * phi;
  const T sg = math::sech(g);
  const T tg = math::tanh(g);
  const T sd = math::sech(p.d2 * phi);
  const T td = math::tanh(p.d2 * phi);
  const T first = sg * (tg * tg - sg * sg) * gp * gp - sg * tg * gpp;
  const T second = p.d1 * p.d2 * p.d2 * sd * (td * td - sd * sd);
  return (first + second) / (1.0 + p.d1);
}

OperatorDimension operator_dimension(const PhaParameters& parameters);

} // namespace pha_qnm

