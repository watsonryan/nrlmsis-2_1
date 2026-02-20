/**
 * @file status.hpp
 * @brief Status codes returned by model operations.
 * @author Watosn
 */
#pragma once

#include <string_view>

namespace msis21 {

enum class Status {
  Ok,
  NotInitialized,
  InvalidInput,
  ParmFileNotFound,
  ParmFileFormatError,
  NumericalError,
};

[[nodiscard]] constexpr std::string_view status_to_string(Status status) noexcept {
  switch (status) {
    case Status::Ok:
      return "Ok";
    case Status::NotInitialized:
      return "NotInitialized";
    case Status::InvalidInput:
      return "InvalidInput";
    case Status::ParmFileNotFound:
      return "ParmFileNotFound";
    case Status::ParmFileFormatError:
      return "ParmFileFormatError";
    case Status::NumericalError:
      return "NumericalError";
  }
  return "UnknownStatus";
}

}  // namespace msis21
