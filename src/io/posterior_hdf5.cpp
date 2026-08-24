#include "pha_qnm/io/posterior.hpp"

#include <hdf5.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace pha_qnm {
namespace {

class H5Handle {
public:
  enum class Kind { file, group, attribute, datatype };
  H5Handle(hid_t id, Kind kind) : id_(id), kind_(kind) {
    if (id_ < 0) throw std::runtime_error("HDF5 operation failed");
  }
  ~H5Handle() {
    if (id_ < 0) return;
    if (kind_ == Kind::file) H5Fclose(id_);
    else if (kind_ == Kind::group) H5Gclose(id_);
    else if (kind_ == Kind::attribute) H5Aclose(id_);
    else H5Tclose(id_);
  }
  H5Handle(const H5Handle&) = delete;
  H5Handle& operator=(const H5Handle&) = delete;
  operator hid_t() const { return id_; }
private:
  hid_t id_;
  Kind kind_;
};

double number_attribute(hid_t object, const char* name) {
  H5Handle attr(H5Aopen(object, name, H5P_DEFAULT), H5Handle::Kind::attribute);
  double value{};
  if (H5Aread(attr, H5T_NATIVE_DOUBLE, &value) < 0) throw std::runtime_error("cannot read numeric HDF5 attribute");
  return value;
}

std::string string_attribute(hid_t object, const char* name) {
  H5Handle attr(H5Aopen(object, name, H5P_DEFAULT), H5Handle::Kind::attribute);
  H5Handle type(H5Aget_type(attr), H5Handle::Kind::datatype);
  if (H5Tis_variable_str(type) > 0) {
    char* raw = nullptr;
    if (H5Aread(attr, type, &raw) < 0) throw std::runtime_error("cannot read string HDF5 attribute");
    std::string value = raw ? raw : "";
    if (raw) H5free_memory(raw);
    return value;
  }
  const std::size_t size = H5Tget_size(type);
  std::vector<char> buffer(size + 1, '\0');
  if (H5Aread(attr, type, buffer.data()) < 0) throw std::runtime_error("cannot read string HDF5 attribute");
  return buffer.data();
}

std::vector<std::string> child_names(hid_t group) {
  hsize_t count{};
  if (H5Gget_num_objs(group, &count) < 0) throw std::runtime_error("cannot enumerate HDF5 group");
  std::vector<std::string> names;
  for (hsize_t i = 0; i < count; ++i) {
    const ssize_t size = H5Gget_objname_by_idx(group, i, nullptr, 0);
    if (size < 0) throw std::runtime_error("cannot read HDF5 child name");
    std::string name(static_cast<std::size_t>(size) + 1, '\0');
    H5Gget_objname_by_idx(group, i, name.data(), name.size());
    name.resize(static_cast<std::size_t>(size));
    names.push_back(std::move(name));
  }
  return names;
}

} // namespace

PosteriorSelection select_maximum_likelihood(const std::filesystem::path& path) {
  H5Handle file(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Handle::Kind::file);
  H5Handle posterior(H5Gopen2(file, "/posterior_samples", H5P_DEFAULT), H5Handle::Kind::group);
  double best = -std::numeric_limits<double>::infinity();
  std::string best_name;
  for (const auto& name : child_names(posterior)) {
    H5Handle sample(H5Gopen2(posterior, name.c_str(), H5P_DEFAULT), H5Handle::Kind::group);
    const double value = number_attribute(sample, "log_likelihood");
    if (value > best) { best = value; best_name = name; }
  }
  if (best_name.empty()) throw std::runtime_error("posterior contains no samples");
  H5Handle sample(H5Gopen2(posterior, best_name.c_str(), H5P_DEFAULT), H5Handle::Kind::group);
  H5Handle params(H5Gopen2(sample, "parameters", H5P_DEFAULT), H5Handle::Kind::group);
  H5Handle critical(H5Gopen2(sample, "critical_point", H5P_DEFAULT), H5Handle::Kind::group);

  PosteriorSelection out;
  out.parameters.sample_id = best_name;
  out.parameters.log_likelihood = best;
  out.parameters.Lambda_MeV = number_attribute(params, "Lambda");
  out.parameters.kappa2 = number_attribute(params, "kappa2");
  out.parameters.gamma = number_attribute(params, "gamma");
  out.parameters.b2 = number_attribute(params, "b2");
  out.parameters.b4 = number_attribute(params, "b4");
  out.parameters.b6 = number_attribute(params, "b6");
  out.parameters.c1 = number_attribute(params, "c1");
  out.parameters.c2 = number_attribute(params, "c2");
  out.parameters.c3 = number_attribute(params, "c3");
  out.parameters.d1 = number_attribute(params, "d1");
  out.parameters.d2 = number_attribute(params, "d2");
  out.critical_point.status = string_attribute(critical, "status");
  out.critical_point.temperature_MeV = number_attribute(critical, "Tc_in_MeV");
  out.critical_point.chemical_potential_MeV = number_attribute(critical, "muc_in_MeV");
  out.critical_point.temperature_error_MeV = number_attribute(critical, "err_Tc_in_MeV");
  out.critical_point.chemical_potential_error_MeV = number_attribute(critical, "err_muc_in_MeV");
  return out;
}

void write_frozen_yaml(const PosteriorSelection& s, const std::filesystem::path& path,
                       const std::string& md5) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open frozen YAML output");
  out << std::setprecision(17)
      << "schema_version: 1\nmodel: polynomial_hyperbolic\nselection:\n"
      << "  rule: maximum_log_likelihood\n  group: posterior_samples\n  sample_id: "
      << s.parameters.sample_id << "\n  log_likelihood: " << s.parameters.log_likelihood
      << "\nsource:\n  doi: 10.5281/zenodo.13830379\n  file: Bayesian_polyhyper_muses.hdf5\n"
      << "  md5: " << md5 << "\nparameters:\n"
      << "  Lambda_MeV: " << s.parameters.Lambda_MeV << "\n  kappa2: " << s.parameters.kappa2
      << "\n  gamma: " << s.parameters.gamma << "\n  b2: " << s.parameters.b2
      << "\n  b4: " << s.parameters.b4 << "\n  b6: " << s.parameters.b6
      << "\n  c1: " << s.parameters.c1 << "\n  c2: " << s.parameters.c2
      << "\n  c3: " << s.parameters.c3 << "\n  d1: " << s.parameters.d1
      << "\n  d2: " << s.parameters.d2 << "\ncritical_point_metadata:\n"
      << "  status: " << s.critical_point.status << "\n  T_c_MeV: " << s.critical_point.temperature_MeV
      << "\n  mu_B_c_MeV: " << s.critical_point.chemical_potential_MeV
      << "\n  error_T_c_MeV: " << s.critical_point.temperature_error_MeV
      << "\n  error_mu_B_c_MeV: " << s.critical_point.chemical_potential_error_MeV << '\n';
}

} // namespace pha_qnm

