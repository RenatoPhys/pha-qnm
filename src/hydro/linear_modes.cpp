#include "pha_qnm/hydro/linear_modes.hpp"

#include <stdexcept>

namespace pha_qnm {

ComplexMatrix charged_hydrodynamic_matrix(Complex omega, double k,
                                          const ChargedHydroParameters& p) {
  ComplexMatrix m(3, 3);
  const Complex i{0.0, 1.0};
  m(0, 0) = -i * omega;
  m(0, 2) = i * k * p.enthalpy;
  m(1, 0) = i * k * p.pressure_energy;
  m(1, 1) = i * k * p.pressure_density;
  m(1, 2) = -i * omega * p.enthalpy + p.longitudinal_viscosity * k * k;
  m(2, 0) = p.conductivity_times_temperature * k * k * p.alpha_energy;
  m(2, 1) = -i * omega + p.conductivity_times_temperature * k * k * p.alpha_density;
  m(2, 2) = i * k * p.density;
  return m;
}

Complex determinant3(const ComplexMatrix& m) {
  if (m.rows() != 3 || m.cols() != 3) throw std::invalid_argument("determinant3 requires 3x3");
  return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) -
         m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
         m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
}

} // namespace pha_qnm

