/**
 * @file parm_reader.cpp
 * @brief Parameter file reader scaffolding.
 * @author Watosn
 */

#include "msis21/detail/parm_reader.hpp"

#include <filesystem>

namespace msis21::detail {

Status validate_parm_file(const std::filesystem::path& parm_path) {
  std::error_code ec;
  if (!std::filesystem::exists(parm_path, ec)) {
    return Status::ParmFileNotFound;
  }
  if (ec) {
    return Status::ParmFileFormatError;
  }
  return Status::ParmFileFormatError;
}

}  // namespace msis21::detail
