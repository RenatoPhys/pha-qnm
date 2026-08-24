#pragma once

#include "pha_qnm/background/equations.hpp"

#include <vector>

namespace pha_qnm {

struct IntegratorOptions {
  double epsilon_horizon = 1.0e-6;
  double r_max = 12.0;
  double initial_step = 1.0e-4;
  double max_step = 5.0e-2;
  double absolute_tolerance = 1.0e-11;
  double relative_tolerance = 1.0e-10;
  std::size_t max_steps = 2'000'000;
};

struct BackgroundPoint {
  double r{};
  BackgroundState state{};
  double constraint{};
  double charge{};
};

struct BackgroundSolution {
  std::vector<BackgroundPoint> points;
  double maximum_normalized_constraint{};
  double maximum_relative_charge_drift{};
  std::size_t accepted_steps{};
  std::size_t rejected_steps{};
};

BackgroundSolution integrate_background(double phi0, double Phi1,
                                        const PhaParameters& parameters,
                                        const IntegratorOptions& options = {});

} // namespace pha_qnm

