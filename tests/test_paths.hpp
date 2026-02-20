/**
 * @file test_paths.hpp
 * @brief Path helpers for test data access.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <string>

#ifndef MSIS21_SOURCE_DIR
#define MSIS21_SOURCE_DIR "."
#endif

inline std::filesystem::path msis21_source_dir() {
  return std::filesystem::path{MSIS21_SOURCE_DIR};
}

inline std::filesystem::path msis21_data_path(const std::string& name) {
  return msis21_source_dir() / "data" / name;
}
