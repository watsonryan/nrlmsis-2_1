/**
 * @file fortran_backend.hpp
 * @brief Bridge to upstream Fortran reference implementation.
 * @author Watosn
 */
#pragma once

#include <filesystem>

#include "msis21/detail/calc.hpp"

namespace msis21::detail {

Status initialize_fortran_backend(const std::filesystem::path& parm_file_path);
[[nodiscard]] bool fortran_backend_available();
[[nodiscard]] CalcResult evaluate_fortran_backend(const Input& in);

}  // namespace msis21::detail
