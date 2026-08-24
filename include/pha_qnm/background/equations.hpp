#pragma once

#include "pha_qnm/model/pha_parameters.hpp"

#include <array>

namespace pha_qnm {

using BackgroundState = std::array<double, 8>;

BackgroundState background_rhs(double r, const BackgroundState& y,
                               const PhaParameters& parameters);
double einstein_constraint(const BackgroundState& y, const PhaParameters& parameters);
double gauss_charge(const BackgroundState& y, const PhaParameters& parameters);

} // namespace pha_qnm

