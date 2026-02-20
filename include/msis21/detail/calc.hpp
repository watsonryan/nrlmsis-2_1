/**
 * @file calc.hpp
 * @brief MSISCALC evaluation plumbing.
 * @author Watosn
 */
#pragma once

#include "msis21/options.hpp"
#include "msis21/status.hpp"
#include "msis21/types.hpp"

namespace msis21::detail {

struct Parameters;

struct CalcResult {
  Status status{Status::Ok};
  Output out{};
};

[[nodiscard]] CalcResult evaluate_msiscalc(const Input& in,
                                           const Options& options,
                                           const Parameters& parameters);

}  // namespace msis21::detail
