#pragma once

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

inline void expect_true(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

inline void expect_near(double actual, double expected, double tolerance,
                        const std::string& message) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
    throw std::runtime_error(message + ": actual=" + std::to_string(actual) +
                             " expected=" + std::to_string(expected));
}

void test_potentials();
void test_horizon_series();
void test_background();
void test_chebyshev();
void test_operators();
void test_hydro();
#ifdef PHA_QNM_HAS_LAPACK
void test_dense_gevp();
#endif
