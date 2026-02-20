/**
 * @file parm_reader.hpp
 * @brief Parameter-file reader interfaces for model initialization.
 * @author Watosn
 */
#pragma once

#include <filesystem>

#include "msis21/status.hpp"

namespace msis21::detail {

Status validate_parm_file(const std::filesystem::path& parm_path);

}  // namespace msis21::detail
