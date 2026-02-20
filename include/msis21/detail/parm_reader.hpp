/**
 * @file parm_reader.hpp
 * @brief Parameter-file reader interfaces for model initialization.
 * @author Watosn
 */
#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

#include "msis21/status.hpp"

namespace msis21::detail {

struct Parameters {
  std::size_t rows{};
  std::size_t cols{};
  std::vector<double> beta;
};

Status load_parameters(const std::filesystem::path& parm_path, Parameters& out);

}  // namespace msis21::detail
