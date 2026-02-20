/**
 * @file test_fortran_reader.cpp
 * @brief Unit tests for Fortran unformatted reader primitives.
 * @author Watosn
 */

#include "msis21/detail/fortran_unformatted_reader.hpp"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <vector>

#include <gtest/gtest.h>

namespace {

void write_u32(std::ofstream& stream, std::uint32_t value) {
  stream.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void write_record(std::ofstream& stream, const std::vector<std::uint8_t>& payload) {
  const auto n = static_cast<std::uint32_t>(payload.size());
  write_u32(stream, n);
  stream.write(reinterpret_cast<const char*>(payload.data()), static_cast<std::streamsize>(payload.size()));
  write_u32(stream, n);
}

}  // namespace

TEST(FortranUnformattedReader, ReadsBackSingleRecord) {
  const auto path = std::filesystem::temp_directory_path() / "fortran_reader_single.bin";
  {
    std::ofstream out(path, std::ios::binary);
    ASSERT_TRUE(static_cast<bool>(out));
    write_record(out, {1, 2, 3, 4, 5});
  }

  msis21::detail::FortranUnformattedReader reader(path);
  const auto rec = reader.read_record();

  ASSERT_EQ(rec.size(), 5U);
  EXPECT_EQ(std::to_integer<unsigned>(rec[0]), 1U);
  EXPECT_EQ(std::to_integer<unsigned>(rec[4]), 5U);
}

TEST(FortranUnformattedReader, ReadsArrayRecord) {
  const auto path = std::filesystem::temp_directory_path() / "fortran_reader_array.bin";
  {
    std::ofstream out(path, std::ios::binary);
    ASSERT_TRUE(static_cast<bool>(out));
    std::vector<std::uint8_t> payload(4 * sizeof(std::uint32_t));
    const std::uint32_t values[4] = {10U, 20U, 30U, 40U};
    std::memcpy(payload.data(), values, payload.size());
    write_record(out, payload);
  }

  msis21::detail::FortranUnformattedReader reader(path);
  const auto values = reader.read_array<std::uint32_t>(4);
  ASSERT_EQ(values.size(), 4U);
  EXPECT_EQ(values[0], 10U);
  EXPECT_EQ(values[3], 40U);
}
