/**
 * @file fortran_unformatted_reader.hpp
 * @brief Reader for Fortran unformatted sequential records.
 * @author Watosn
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace msis21::detail {

class FortranUnformattedReader {
 public:
  enum class MarkerSize : std::uint8_t { Bits32, Bits64 };
  enum class ByteOrder : std::uint8_t { LittleEndian, BigEndian };

  explicit FortranUnformattedReader(std::filesystem::path path,
                                    MarkerSize marker_size = MarkerSize::Bits32,
                                    ByteOrder byte_order = ByteOrder::LittleEndian);

  [[nodiscard]] std::vector<std::byte> read_record();
  [[nodiscard]] bool eof();

  template <class T>
  [[nodiscard]] T read_scalar() {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto record = read_record();
    if (record.size() != sizeof(T)) {
      throw std::runtime_error("record size mismatch");
    }
    T value{};
    std::memcpy(&value, record.data(), sizeof(T));
    if constexpr (sizeof(T) == 4 || sizeof(T) == 8) {
      if (needs_swap_) {
        byteswap_inplace(reinterpret_cast<std::byte*>(&value), sizeof(T));
      }
    }
    return value;
  }

  template <class T>
  [[nodiscard]] std::vector<T> read_array(std::size_t count) {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto record = read_record();
    if (record.size() != sizeof(T) * count) {
      throw std::runtime_error("record size mismatch");
    }
    std::vector<T> values(count);
    std::memcpy(values.data(), record.data(), record.size());
    if constexpr (sizeof(T) == 4 || sizeof(T) == 8) {
      if (needs_swap_) {
        for (auto& value : values) {
          byteswap_inplace(reinterpret_cast<std::byte*>(&value), sizeof(T));
        }
      }
    }
    return values;
  }

 private:
  static bool host_is_little_endian();
  static void byteswap_inplace(std::byte* data, std::size_t size);
  std::uint64_t read_marker();

  std::ifstream stream_;
  MarkerSize marker_size_;
  bool needs_swap_;
};

}  // namespace msis21::detail
