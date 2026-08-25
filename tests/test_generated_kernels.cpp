#include "test_support.hpp"

#include "pha_qnm/perturbations/generated/helicity0_kernels.generated.hpp"
#include "pha_qnm/perturbations/generated/helicity1_kernels.generated.hpp"

#include <algorithm>
#include <cmath>

namespace {

pha_qnm::generated::EmdPoint on_shell_point() {
  pha_qnm::generated::EmdPoint b;
  b.A = 0.3;
  b.Ar = 0.7;
  b.h = 0.2;
  b.hr = 0.4;
  b.phir = 0.5;
  b.Phir = 0.1;
  b.f = 1.2;
  b.fphi = 0.2;
  b.fphiphi = -0.1;
  b.fr = b.fphi * b.phir;
  b.fphir = b.fphiphi * b.phir;
  b.Arr = -b.phir * b.phir / 6.0;
  b.hrr = -4.0 * b.Ar * b.hr + std::exp(-2.0 * b.A) * b.f * b.Phir * b.Phir;
  b.Phirr = -(2.0 * b.Ar + b.fphi * b.phir / b.f) * b.Phir;
  b.Vphi = -0.8;
  b.phirr = (b.Vphi - 0.5 * std::exp(-2.0 * b.A) * b.fphi * b.Phir * b.Phir
             - (4.0 * b.h * b.Ar + b.hr) * b.phir) / b.h;
  b.V = -0.5 * (b.h * (24.0 * b.Ar * b.Ar - b.phir * b.phir)
                + 6.0 * b.Ar * b.hr
                + std::exp(-2.0 * b.A) * b.f * b.Phir * b.Phir);
  b.Vphiphi = -2.0;
  return b;
}

template <class Coefficients>
void expect_finite(const Coefficients& coefficients, const std::string& label) {
  const bool finite = std::ranges::all_of(coefficients.values, [](const auto& value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  });
  expect_true(finite, label + " emitted a non-finite coefficient");
}

} // namespace

void test_generated_kernels() {
  using pha_qnm::Complex;
  using namespace pha_qnm::generated;
  const auto b = on_shell_point();
  const Complex sigma{-0.17, 0.31};
  const Complex momentum{0.0, 0.4};
  const auto h1 = helicity1_coefficients(b, sigma, momentum);
  const auto h0 = helicity0_coefficients(b, sigma, momentum);
  expect_finite(h1, "helicity one");
  expect_finite(h0, "helicity zero");
  expect_true(helicity1_fields.size() == 3 && helicity1_equations.size() == 5,
              "helicity-one generated dimensions changed");
  expect_true(helicity0_fields.size() == 7 && helicity0_equations.size() == 14,
              "helicity-zero generated dimensions changed");

  const auto k0 = helicity1_coefficients(b, sigma, Complex{});
  const std::size_t constraint = 4;
  const auto hvx_value = k0(constraint, 0, 0);
  expect_true(std::abs(hvx_value) < 1.0e-13,
              "on-shell k=0 transverse constraint retained an algebraic metric term");
  expect_true(std::abs(k0(constraint, 0, 1) + sigma / 2.0) < 1.0e-13,
              "k=0 transverse metric constraint coefficient mismatch");
  const auto expected_ax = -0.5 * b.Phir * b.f * sigma * std::exp(-2.0 * b.A);
  expect_true(std::abs(k0(constraint, 2, 0) - expected_ax) < 1.0e-13,
              "k=0 transverse Maxwell mixing coefficient mismatch");

  auto neutral = b;
  neutral.Phir = 0.0;
  neutral.Phirr = 0.0;
  const auto neutral_h1 = helicity1_coefficients(neutral, sigma, momentum);
  const std::size_t maxwell_x = 3;
  expect_true(std::abs(neutral_h1(maxwell_x, 0, 0)) < 1.0e-14
              && std::abs(neutral_h1(maxwell_x, 0, 1)) < 1.0e-14,
              "neutral Maxwell equation did not decouple from metric fields");

  const auto h0_zero = helicity0_coefficients(b, Complex{}, momentum);
  const auto h0_one = helicity0_coefficients(b, Complex{1.0, 0.0}, momentum);
  const auto h0_two = helicity0_coefficients(b, Complex{2.0, 0.0}, momentum);
  for (std::size_t equation = 0; equation < 7; ++equation) {
    for (std::size_t field = 0; field < 7; ++field) {
      for (std::size_t derivative = 0; derivative < 3; ++derivative) {
        const auto second_difference = h0_two(equation, field, derivative)
                                     - 2.0 * h0_one(equation, field, derivative)
                                     + h0_zero(equation, field, derivative);
        expect_true(std::abs(second_difference) < 2.0e-13,
                    "helicity-zero characteristic kernel is nonlinear in frequency");
      }
    }
  }
}
