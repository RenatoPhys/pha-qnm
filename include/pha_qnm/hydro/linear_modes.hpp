#pragma once

#include "pha_qnm/core/types.hpp"

namespace pha_qnm {

struct ChargedHydroParameters {
  double enthalpy{};
  double density{};
  double pressure_energy{};
  double pressure_density{};
  double alpha_energy{};
  double alpha_density{};
  double conductivity_times_temperature{};
  double longitudinal_viscosity{};
};

ComplexMatrix charged_hydrodynamic_matrix(Complex omega, double momentum,
                                          const ChargedHydroParameters& parameters);
Complex determinant3(const ComplexMatrix& matrix);

} // namespace pha_qnm

