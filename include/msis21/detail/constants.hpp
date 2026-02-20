/**
 * @file constants.hpp
 * @brief Core constants for NRLMSIS 2.1 parameter handling.
 * @author Watosn
 */
#pragma once

#include <cstddef>

namespace msis21::detail {

inline constexpr std::size_t kMaxBasisFunctions = 512;
inline constexpr std::size_t kParmColumnCount = 131;

}  // namespace msis21::detail
