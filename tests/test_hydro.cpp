#include "test_support.hpp"

#include "pha_qnm/hydro/linear_modes.hpp"

void test_hydro() {
  const pha_qnm::ChargedHydroParameters p{10.0, 2.0, 0.3, 0.2, -0.01, 0.05, 0.4, 0.6};
  const auto zero = pha_qnm::charged_hydrodynamic_matrix({}, 0.0, p);
  expect_true(std::abs(pha_qnm::determinant3(zero)) == 0.0,
              "conservation roots must meet at omega=k=0");
  pha_qnm::ComplexMatrix identity(3, 3);
  identity(0, 0) = identity(1, 1) = identity(2, 2) = {1.0, 0.0};
  expect_near(pha_qnm::determinant3(identity).real(), 1.0, 0.0, "3x3 determinant");
}

