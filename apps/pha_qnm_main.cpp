#include "pha_qnm/background/integrator.hpp"
#include "pha_qnm/io/posterior.hpp"
#ifdef PHA_QNM_HAS_NETWORK
#include "pha_qnm/io/download.hpp"
#endif
#include "pha_qnm/model/potentials.hpp"
#include "pha_qnm/spectral/chebyshev.hpp"

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

void usage() {
  std::cout << "pha-qnm 0.1.0\n"
               "Usage:\n"
               "  pha-qnm model info\n"
               "  pha-qnm model extract-map <posterior.h5> <output.yaml>\n"
               "  pha-qnm data download <output.h5>\n"
               "  pha-qnm data verify <posterior.h5>\n"
               "  pha-qnm background point <phi0> <Phi1> [r_max]\n"
               "  pha-qnm validate all\n"
               "Global workflow flags reserved for scans: --dry-run --resume --log-level\n";
}

double number(const char* text, std::string_view label) {
  std::size_t used{};
  const double value = std::stod(text, &used);
  if (used != std::string(text).size() || !std::isfinite(value))
    throw std::invalid_argument("invalid " + std::string(label));
  return value;
}

int validate_all() {
  const pha_qnm::PhaParameters p;
  const auto dimension = pha_qnm::operator_dimension(p);
  if (!(dimension.delta > 2.0 && dimension.delta < 4.0))
    throw std::runtime_error("operator dimension validation failed");
  const auto grid = pha_qnm::chebyshev_lobatto(16, 0.0, 1.0);
  double polynomial_error = 0.0;
  for (std::size_t i = 0; i < grid.nodes.size(); ++i) {
    double derivative{};
    for (std::size_t j = 0; j < grid.nodes.size(); ++j)
      derivative += grid.D1(i, j) * std::pow(grid.nodes[j], 4);
    polynomial_error = std::max(polynomial_error,
                                std::abs(derivative - 4.0 * std::pow(grid.nodes[i], 3)));
  }
  if (polynomial_error > 1.0e-10) throw std::runtime_error("Chebyshev validation failed");
  const auto background = pha_qnm::integrate_background(1.0, 0.05, p,
      {.epsilon_horizon = 1.0e-6, .r_max = 0.5, .initial_step = 1.0e-5,
       .max_step = 2.0e-3, .absolute_tolerance = 1.0e-11,
       .relative_tolerance = 1.0e-10, .max_steps = 1'000'000});
  std::cout << std::setprecision(12)
            << "validation.operator_dimension=" << dimension.delta << '\n'
            << "validation.chebyshev_error=" << polynomial_error << '\n'
            << "validation.background_constraint=" << background.maximum_normalized_constraint << '\n'
            << "validation.gauss_drift=" << background.maximum_relative_charge_drift << '\n';
  return 0;
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) { usage(); return 0; }
    const std::string command = argv[1];
    if (command == "model" && argc >= 3 && std::string(argv[2]) == "info") {
      const pha_qnm::PhaParameters p;
      const auto d = pha_qnm::operator_dimension(p);
      std::cout << std::setprecision(17)
                << "sample_id=" << p.sample_id << '\n'
                << "log_likelihood=" << p.log_likelihood << '\n'
                << "m2=" << d.mass_squared << '\n'
                << "Delta=" << d.delta << '\n'
                << "nu=" << d.nu << '\n';
      return 0;
    }
    if (command == "model" && argc == 5 && std::string(argv[2]) == "extract-map") {
      const auto selection = pha_qnm::select_maximum_likelihood(argv[3]);
      pha_qnm::write_frozen_yaml(selection, argv[4], "7fd567dccfaea48095ca5df53a8c17d6");
      std::cout << "selected=" << selection.parameters.sample_id << '\n';
      return 0;
    }
#ifdef PHA_QNM_HAS_NETWORK
    if (command == "data" && argc == 4 && std::string(argv[2]) == "download") {
      constexpr const char* url = "https://zenodo.org/records/13830379/files/Bayesian_polyhyper_muses.hdf5?download=1";
      pha_qnm::download_file(url, argv[3]);
      pha_qnm::require_md5(argv[3], "7fd567dccfaea48095ca5df53a8c17d6");
      std::cout << "md5=7fd567dccfaea48095ca5df53a8c17d6\n";
      return 0;
    }
    if (command == "data" && argc == 4 && std::string(argv[2]) == "verify") {
      pha_qnm::require_md5(argv[3], "7fd567dccfaea48095ca5df53a8c17d6");
      std::cout << "md5=7fd567dccfaea48095ca5df53a8c17d6\n";
      return 0;
    }
#endif
    if (command == "background" && argc >= 5 && std::string(argv[2]) == "point") {
      pha_qnm::IntegratorOptions options;
      if (argc >= 6) options.r_max = number(argv[5], "r_max");
      const auto solution = pha_qnm::integrate_background(number(argv[3], "phi0"),
          number(argv[4], "Phi1"), pha_qnm::PhaParameters{}, options);
      const auto& last = solution.points.back();
      std::cout << std::setprecision(17) << "r=" << last.r << "\nA=" << last.state[0]
                << "\nh=" << last.state[1] << "\nphi=" << last.state[2]
                << "\nPhi=" << last.state[3]
                << "\nmax_constraint=" << solution.maximum_normalized_constraint
                << "\nmax_gauss_drift=" << solution.maximum_relative_charge_drift << '\n';
      return 0;
    }
    if (command == "validate" && argc >= 3 && std::string(argv[2]) == "all")
      return validate_all();
    usage();
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
