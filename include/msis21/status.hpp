/**
 * @file status.hpp
 * @brief Status codes returned by model operations.
 * @author Watosn
 */
#pragma once

namespace msis21 {

enum class Status {
  Ok,
  NotInitialized,
  InvalidInput,
  ParmFileNotFound,
  ParmFileFormatError,
  NumericalError,
};

}  // namespace msis21
