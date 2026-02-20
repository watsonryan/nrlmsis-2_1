/**
 * @file fortran_unformatted_reader.cpp
 * @brief Fortran unformatted sequential reader implementation.
 * @author Watosn
 */

#include "msis21/detail/fortran_unformatted_reader.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <stdexcept>

namespace msis21::detail {

FortranUnformattedReader::FortranUnformattedReader(std::filesystem::path path,
                                                   MarkerSize marker_size,
                                                   ByteOrder byte_order)
    : stream_(std::move(path), std::ios::binary), marker_size_(marker_size), needs_swap_(false) {
  if (!stream_) {
    throw std::runtime_error("failed to open file");
  }
  const bool file_little = (byte_order == ByteOrder::LittleEndian);
  needs_swap_ = (host_is_little_endian() != file_little);
}

std::vector<std::byte> FortranUnformattedReader::read_record() {
  const auto record_size = read_marker();
  std::vector<std::byte> bytes(record_size);
  stream_.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
  if (stream_.gcount() != static_cast<std::streamsize>(bytes.size())) {
    throw std::runtime_error("failed to read record payload");
  }
  const auto trailing_size = read_marker();
  if (record_size != trailing_size) {
    throw std::runtime_error("record marker mismatch");
  }
  return bytes;
}

bool FortranUnformattedReader::eof() { return stream_.peek() == std::ifstream::traits_type::eof(); }

bool FortranUnformattedReader::host_is_little_endian() {
#if defined(__cpp_lib_endian)
  return std::endian::native == std::endian::little;
#else
  const std::uint16_t marker = 1;
  return *reinterpret_cast<const std::uint8_t*>(&marker) == 1;
#endif
}

void FortranUnformattedReader::byteswap_inplace(std::byte* data, std::size_t size) {
  std::reverse(data, data + static_cast<std::ptrdiff_t>(size));
}

std::uint64_t FortranUnformattedReader::read_marker() {
  if (marker_size_ == MarkerSize::Bits32) {
    std::array<std::byte, 4> raw{};
    stream_.read(reinterpret_cast<char*>(raw.data()), static_cast<std::streamsize>(raw.size()));
    if (stream_.gcount() != static_cast<std::streamsize>(raw.size())) {
      throw std::runtime_error("failed to read 32-bit marker");
    }
    if (needs_swap_) {
      byteswap_inplace(raw.data(), raw.size());
    }
    std::uint32_t value{};
    std::memcpy(&value, raw.data(), sizeof(value));
    return value;
  }

  std::array<std::byte, 8> raw{};
  stream_.read(reinterpret_cast<char*>(raw.data()), static_cast<std::streamsize>(raw.size()));
  if (stream_.gcount() != static_cast<std::streamsize>(raw.size())) {
    throw std::runtime_error("failed to read 64-bit marker");
  }
  if (needs_swap_) {
    byteswap_inplace(raw.data(), raw.size());
  }
  std::uint64_t value{};
  std::memcpy(&value, raw.data(), sizeof(value));
  return value;
}

}  // namespace msis21::detail
