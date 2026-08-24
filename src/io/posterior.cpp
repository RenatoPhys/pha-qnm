#include "pha_qnm/io/posterior.hpp"

#ifndef PHA_QNM_HAS_HDF5
#include <stdexcept>

namespace pha_qnm {

PosteriorSelection select_maximum_likelihood(const std::filesystem::path&) {
  throw std::runtime_error("native HDF5 support is disabled; configure PHA_QNM_ENABLE_HDF5=ON");
}

void write_frozen_yaml(const PosteriorSelection&, const std::filesystem::path&,
                       const std::string&) {
  throw std::runtime_error("native HDF5 support is disabled; configure PHA_QNM_ENABLE_HDF5=ON");
}

} // namespace pha_qnm
#endif

