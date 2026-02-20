/**
 * @file calc.cpp
 * @brief MSISCALC evaluation plumbing implementation.
 * @author Watosn
 */

#include "msis21/detail/calc.hpp"

#include <array>

#include "msis21/detail/constants.hpp"
#include "msis21/detail/gfn.hpp"
#include "msis21/detail/parm_reader.hpp"

namespace msis21::detail {

namespace {

bool validate_input(const Input& in) {
  if (in.sec < 0.0 || in.sec > 86400.0) {
    return false;
  }
  if (in.glat_deg < -90.0 || in.glat_deg > 90.0) {
    return false;
  }
  if (in.stl_hr < 0.0 || in.stl_hr > 24.0) {
    return false;
  }
  if (in.f107a < 0.0 || in.f107 < 0.0 || in.ap < 0.0) {
    return false;
  }
  return true;
}

}  // namespace

CalcResult evaluate_msiscalc(const Input& in, const Options& options, const Parameters& /*parameters*/) {
  if (!validate_input(in)) {
    return CalcResult{.status = Status::InvalidInput};
  }

  GlobeCalculator globe;
  GlobeInput globe_input;
  globe_input.doy = static_cast<double>(in.iyd % 1000);
  globe_input.utsec = in.sec;
  globe_input.lat = in.glat_deg;
  globe_input.lon = in.glon_deg;
  globe_input.sfluxavg = in.f107a;
  globe_input.sflux = in.f107;
  globe_input.ap = {in.ap, in.ap, in.ap, in.ap, in.ap, in.ap, in.ap};

  std::array<bool, kMaxBasisFunctions> switches;
  switches.fill(true);
  for (int i = 0; i < 25; ++i) {
    if (options.switches.legacy[static_cast<std::size_t>(i)] == 0) {
      switches[static_cast<std::size_t>(i)] = false;
    }
  }

  const auto basis = globe.globe(globe_input, switches);
  (void)basis;

  return CalcResult{.status = Status::NumericalError};
}

}  // namespace msis21::detail
