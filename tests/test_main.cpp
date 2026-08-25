#include "test_support.hpp"

#include <exception>
#include <iostream>
#include <utility>
#include <vector>

int main() {
  const std::vector<std::pair<const char*, void (*)()>> tests = {
      {"potentials", test_potentials},
      {"horizon_series", test_horizon_series},
      {"background", test_background},
      {"chebyshev", test_chebyshev},
      {"operators", test_operators},
      {"generated_kernels", test_generated_kernels},
      {"hydro", test_hydro}
#ifdef PHA_QNM_HAS_LAPACK
      , {"dense_gevp", test_dense_gevp}
#endif
  };
  int failures = 0;
  for (const auto& [name, test] : tests) {
    try { test(); std::cout << "PASS " << name << '\n'; }
    catch (const std::exception& e) { ++failures; std::cerr << "FAIL " << name << ": " << e.what() << '\n'; }
  }
  return failures == 0 ? 0 : 1;
}
