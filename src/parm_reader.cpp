/**
 * @file parm_reader.cpp
 * @brief Parameter file reader scaffolding.
 * @author Watosn
 */

#include "msis21/detail/parm_reader.hpp"

#include <cstring>
#include <fstream>

#include "msis21/detail/constants.hpp"

namespace msis21::detail {

Status load_parameters(const std::filesystem::path& parm_path, Parameters& out) {
  std::error_code ec;
  if (!std::filesystem::exists(parm_path, ec)) {
    return Status::ParmFileNotFound;
  }
  if (ec) {
    return Status::ParmFileFormatError;
  }

  std::ifstream stream(parm_path, std::ios::binary | std::ios::ate);
  if (!stream) {
    return Status::ParmFileFormatError;
  }

  const auto file_size = stream.tellg();
  if (file_size <= 0) {
    return Status::ParmFileFormatError;
  }
  constexpr std::size_t element_size = sizeof(double);
  const auto size = static_cast<std::size_t>(file_size);
  if (size % element_size != 0) {
    return Status::ParmFileFormatError;
  }

  const std::size_t element_count = size / element_size;
  if (element_count % kMaxBasisFunctions != 0) {
    return Status::ParmFileFormatError;
  }
  const std::size_t column_count = element_count / kMaxBasisFunctions;
  if (column_count != kParmColumnCount) {
    return Status::ParmFileFormatError;
  }

  out.rows = kMaxBasisFunctions;
  out.cols = column_count;
  out.beta.resize(element_count);

  stream.seekg(0, std::ios::beg);
  stream.read(reinterpret_cast<char*>(out.beta.data()), static_cast<std::streamsize>(size));
  if (stream.gcount() != static_cast<std::streamsize>(size)) {
    return Status::ParmFileFormatError;
  }

  return Status::Ok;
}

}  // namespace msis21::detail
