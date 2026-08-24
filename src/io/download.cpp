#include "pha_qnm/io/download.hpp"

#include <curl/curl.h>
#include <openssl/evp.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace pha_qnm {
namespace {

size_t write_chunk(char* data, size_t size, size_t count, void* user) {
  auto& stream = *static_cast<std::ofstream*>(user);
  const size_t bytes = size * count;
  stream.write(data, static_cast<std::streamsize>(bytes));
  return stream ? bytes : 0;
}

} // namespace

void download_file(const std::string& url, const std::filesystem::path& output) {
  std::filesystem::create_directories(output.parent_path());
  const auto temporary = output.string() + ".part";
  std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
  if (!stream) throw std::runtime_error("cannot open download output");
  std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> curl(curl_easy_init(), curl_easy_cleanup);
  if (!curl) throw std::runtime_error("curl initialization failed");
  curl_easy_setopt(curl.get(), CURLOPT_URL, url.c_str());
  curl_easy_setopt(curl.get(), CURLOPT_FOLLOWLOCATION, 1L);
  curl_easy_setopt(curl.get(), CURLOPT_FAILONERROR, 1L);
  curl_easy_setopt(curl.get(), CURLOPT_WRITEFUNCTION, write_chunk);
  curl_easy_setopt(curl.get(), CURLOPT_WRITEDATA, &stream);
  curl_easy_setopt(curl.get(), CURLOPT_USERAGENT, "pha-qnm/0.1.0");
  const CURLcode status = curl_easy_perform(curl.get());
  stream.close();
  if (status != CURLE_OK) {
    std::filesystem::remove(temporary);
    throw std::runtime_error(std::string("download failed: ") + curl_easy_strerror(status));
  }
  std::filesystem::rename(temporary, output);
}

std::string md5_file(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open file for checksum");
  std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> context(EVP_MD_CTX_new(), EVP_MD_CTX_free);
  if (!context || EVP_DigestInit_ex(context.get(), EVP_md5(), nullptr) != 1)
    throw std::runtime_error("OpenSSL MD5 initialization failed");
  std::array<char, 1 << 16> buffer{};
  while (input) {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const auto count = input.gcount();
    if (count > 0 && EVP_DigestUpdate(context.get(), buffer.data(), static_cast<std::size_t>(count)) != 1)
      throw std::runtime_error("OpenSSL MD5 update failed");
  }
  std::array<unsigned char, EVP_MAX_MD_SIZE> digest{};
  unsigned int length{};
  if (EVP_DigestFinal_ex(context.get(), digest.data(), &length) != 1)
    throw std::runtime_error("OpenSSL MD5 finalization failed");
  std::ostringstream out;
  out << std::hex << std::setfill('0');
  for (unsigned int i = 0; i < length; ++i) out << std::setw(2) << static_cast<unsigned>(digest[i]);
  return out.str();
}

void require_md5(const std::filesystem::path& path, const std::string& expected) {
  std::string normalized = expected;
  std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  const std::string actual = md5_file(path);
  if (actual != normalized)
    throw std::runtime_error("MD5 mismatch: expected " + normalized + ", obtained " + actual);
}

} // namespace pha_qnm

