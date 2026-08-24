#pragma once

#include "pha_qnm/model/pha_parameters.hpp"

#include <filesystem>
#include <string>

namespace pha_qnm {

struct CriticalPointMetadata {
  std::string status;
  double temperature_MeV{};
  double chemical_potential_MeV{};
  double temperature_error_MeV{};
  double chemical_potential_error_MeV{};
};

struct PosteriorSelection {
  PhaParameters parameters;
  CriticalPointMetadata critical_point;
};

PosteriorSelection select_maximum_likelihood(const std::filesystem::path& hdf5_path);
void write_frozen_yaml(const PosteriorSelection& selection,
                       const std::filesystem::path& output_path,
                       const std::string& source_md5);

} // namespace pha_qnm

