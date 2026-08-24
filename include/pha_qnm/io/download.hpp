#pragma once

#include <filesystem>
#include <string>

namespace pha_qnm {

void download_file(const std::string& url, const std::filesystem::path& output);
std::string md5_file(const std::filesystem::path& path);
void require_md5(const std::filesystem::path& path, const std::string& expected_lowercase);

} // namespace pha_qnm

