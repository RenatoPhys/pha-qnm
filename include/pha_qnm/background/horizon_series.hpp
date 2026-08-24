#pragma once

#include "pha_qnm/model/pha_parameters.hpp"

#include <array>

namespace pha_qnm {

struct HorizonSeries {
  std::array<double, 5> A{};
  std::array<double, 5> h{};
  std::array<double, 5> phi{};
  std::array<double, 5> Phi{};

  [[nodiscard]] std::array<double, 8> state(double r) const;
};

HorizonSeries make_horizon_series(double phi0, double Phi1,
                                  const PhaParameters& parameters);
double maximum_horizon_electric_field(double phi0, const PhaParameters& parameters);

} // namespace pha_qnm

