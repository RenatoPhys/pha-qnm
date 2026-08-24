#pragma once

#include <string>

namespace pha_qnm {

struct PhaParameters {
  double Lambda_MeV = 1123.8704743960677;
  double kappa2 = 11.37047249736345;
  double gamma = 0.5904238801815029;
  double b2 = 0.32722698958661095;
  double b4 = -0.04801342482759139;
  double b6 = 0.0009385148172188441;
  double c1 = 0.008800483500598342;
  double c2 = 0.15579934656052294;
  double c3 = 0.05281851198008199;
  double d1 = 1.7177450982418923;
  double d2 = 1679.1777562624359;
  std::string sample_id = "sample74";
  double log_likelihood = 31.075137487467885;
};

struct OperatorDimension {
  double mass_squared{};
  double delta{};
  double nu{};
};

} // namespace pha_qnm

