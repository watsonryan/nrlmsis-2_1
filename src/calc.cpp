/**
 * @file calc.cpp
 * @brief MSISCALC evaluation plumbing implementation.
 * @author Watosn
 */

#include "msis21/detail/calc.hpp"

#include <array>
#include <cmath>
#include <vector>

#include "msis21/detail/constants.hpp"
#include "msis21/detail/fortran_backend.hpp"
#include "msis21/detail/gfn.hpp"
#include "msis21/detail/log.hpp"
#include "msis21/detail/parm_reader.hpp"
#include "msis21/detail/utils.hpp"
#include <spdlog/spdlog.h>

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
  if (in.alt_km < -1.0 || in.alt_km > 2000.0) {
    return false;
  }
  return true;
}

double beta_at(const Parameters& parameters, int row, int col) {
  return parameters.beta[static_cast<std::size_t>(row) +
                         static_cast<std::size_t>(col) * parameters.rows];
}

double linear_dot(const Parameters& parameters, const std::array<double, kMaxBasisFunctions>& basis, int col) {
  double sum = 0.0;
  for (int row = 0; row <= kMbf; ++row) {
    sum += beta_at(parameters, row, col) * basis[static_cast<std::size_t>(row)];
  }
  return sum;
}

void set_range(std::array<bool, kMaxBasisFunctions>& swg, int begin, int end, bool value) {
  for (int i = begin; i <= end; ++i) {
    swg[static_cast<std::size_t>(i)] = value;
  }
}

void set_indices(std::array<bool, kMaxBasisFunctions>& swg, const std::vector<int>& idx, bool value) {
  for (const int i : idx) {
    swg[static_cast<std::size_t>(i)] = value;
  }
}

std::array<bool, kMaxBasisFunctions> legacy_switches_to_basis(const Options& options) {
  std::array<bool, kMaxBasisFunctions> swg;
  swg.fill(true);

  std::array<double, 25> swleg{};
  std::array<double, 25> swc{};
  for (int i = 0; i < 25; ++i) {
    const double sv = static_cast<double>(options.switches.legacy[static_cast<std::size_t>(i)]);
    swleg[static_cast<std::size_t>(i)] = std::fmod(sv, 2.0);
    swc[static_cast<std::size_t>(i)] =
        (std::abs(sv) == 1.0 || std::abs(sv) == 2.0) ? 1.0 : 0.0;
  }

  swg[0] = true;
  set_range(swg, kCsfx, kCsfx + kNsfx - 1, swleg[0] == 1.0);
  swg[310] = (swleg[0] == 1.0);
  set_range(swg, 1, 6, swleg[1] == 1.0);
  set_range(swg, 304, 305, swleg[1] == 1.0);
  set_range(swg, 311, 312, swleg[1] == 1.0);
  set_range(swg, 313, 314, swleg[1] == 1.0);
  set_indices(swg, {7, 8, 11, 12, 15, 16, 19, 20}, swleg[2] == 1.0);
  set_range(swg, 306, 307, swleg[2] == 1.0);
  set_indices(swg, {21, 22, 25, 26, 29, 30, 33, 34}, swleg[3] == 1.0);
  set_range(swg, 308, 309, swleg[3] == 1.0);
  set_indices(swg, {9, 10, 13, 14, 17, 18}, swleg[4] == 1.0);
  set_indices(swg, {23, 24, 27, 28, 31, 32}, swleg[5] == 1.0);
  set_range(swg, 35, 94, swleg[6] == 1.0);
  set_range(swg, 300, 303, swleg[6] == 1.0);
  set_range(swg, 95, 144, swleg[7] == 1.0);
  set_range(swg, 145, 184, swleg[13] == 1.0);
  set_range(swg, kCmag, kCmag + 1, false);
  if ((swleg[8] > 0.0) || (swleg[12] == 1.0)) {
    swg[static_cast<std::size_t>(kCmag)] = true;
    swg[static_cast<std::size_t>(kCmag + 1)] = true;
  }
  if (swleg[8] < 0.0) {
    swg[static_cast<std::size_t>(kCmag)] = false;
    swg[static_cast<std::size_t>(kCmag + 1)] = true;
  }
  set_range(swg, kCmag + 2, kCmag + 12, swleg[8] == 1.0);
  set_range(swg, kCmag + 28, kCmag + 40, swleg[8] == -1.0);
  set_range(swg, kCspw, kCsfx - 1, (swleg[10] == 1.0) && (swleg[9] == 1.0));
  set_range(swg, kCut, kCut + kNut - 1, (swleg[11] == 1.0) && (swleg[9] == 1.0));
  set_range(swg, kCmag + 13, kCmag + 25, (swleg[12] == 1.0) && (swleg[9] == 1.0));
  set_range(swg, kCmag + 41, kCmag + 53, (swleg[12] == 1.0) && (swleg[9] == 1.0));

  set_range(swg, kCsfxMod, kCsfxMod + kNsfxMod - 1, swc[0] == 1.0);
  if (swc[0] == 0.0) {
    set_range(swg, 302, 303, false);
    set_range(swg, 304, 305, false);
    set_range(swg, 306, 307, false);
    set_range(swg, 308, 309, false);
    set_range(swg, 311, 314, false);
    swg[447] = false;
    swg[454] = false;
  }
  if (swc[1] == 0.0) {
    set_range(swg, 9, 20, false);
    set_range(swg, 23, 34, false);
    set_range(swg, 35, 184, false);
    set_range(swg, 185, 294, false);
    set_range(swg, 392, 414, false);
    set_range(swg, 420, 442, false);
    set_range(swg, 449, 453, false);
  }
  if (swc[2] == 0.0) {
    set_range(swg, 201, 204, false);
    set_range(swg, 209, 212, false);
    set_range(swg, 217, 220, false);
    set_range(swg, 255, 258, false);
    set_range(swg, 263, 266, false);
    set_range(swg, 271, 274, false);
    set_range(swg, 306, 307, false);
  }
  if (swc[3] == 0.0) {
    set_range(swg, 225, 228, false);
    set_range(swg, 233, 236, false);
    set_range(swg, 241, 244, false);
    set_range(swg, 275, 278, false);
    set_range(swg, 283, 286, false);
    set_range(swg, 291, 294, false);
    set_range(swg, 308, 309, false);
  }
  if (swc[4] == 0.0) {
    set_range(swg, 47, 70, false);
    set_range(swg, 105, 124, false);
    set_range(swg, 153, 168, false);
    set_range(swg, 197, 200, false);
    set_range(swg, 205, 216, false);
    set_range(swg, 259, 270, false);
    set_range(swg, 394, 397, false);
    set_range(swg, 407, 410, false);
    set_range(swg, 422, 425, false);
    set_range(swg, 435, 438, false);
    swg[446] = false;
  }
  if (swc[5] == 0.0) {
    set_range(swg, 221, 224, false);
    set_range(swg, 229, 232, false);
    set_range(swg, 237, 240, false);
    set_range(swg, 279, 282, false);
    set_range(swg, 287, 290, false);
  }
  if (swc[6] == 0.0) {
    set_range(swg, 398, 401, false);
    set_range(swg, 426, 429, false);
  }
  if (swc[10] == 0.0) {
    set_range(swg, 402, 410, false);
    set_range(swg, 430, 438, false);
    set_range(swg, 452, 453, false);
  }
  if (swc[11] == 0.0) {
    set_range(swg, 411, 414, false);
    set_range(swg, 439, 440, false);
  }
  return swg;
}

}  // namespace

CalcResult evaluate_msiscalc(const Input& in, const Options& options, const Parameters& parameters) {
  configure_logging_once();
  if (!validate_input(in)) {
    spdlog::warn("Invalid input for evaluate: iyd={} sec={} alt_km={}", in.iyd, in.sec, in.alt_km);
    return CalcResult{.status = Status::InvalidInput};
  }
  if (fortran_backend_available()) {
    return evaluate_fortran_backend(in);
  }
  if (parameters.rows != kMaxBasisFunctions || parameters.cols != kParmColumnCount) {
    return CalcResult{.status = Status::ParmFileFormatError};
  }

  GlobeCalculator globe;
  GlobeInput globe_input;
  globe_input.doy = static_cast<double>(in.iyd % 1000) + in.sec / 86400.0;
  globe_input.utsec = in.sec;
  globe_input.lat = in.glat_deg;
  globe_input.lon = in.glon_deg;
  globe_input.sfluxavg = in.f107a;
  globe_input.sflux = in.f107;
  globe_input.ap = {in.ap, in.ap, in.ap, in.ap, in.ap, in.ap, in.ap};

  const auto switches = legacy_switches_to_basis(options);

  const auto basis = globe.globe(globe_input, switches);
  // Temperature anchors from the TN subset: [itb0=21, itgb0=22, itex=23].
  double tb0 = linear_dot(parameters, basis, 21);
  double tgb0 = linear_dot(parameters, basis, 22);
  double tex = linear_dot(parameters, basis, 23);

  if (!std::isfinite(tb0) || tb0 < 120.0 || tb0 > 600.0) {
    tb0 = 190.0;
  }
  if (!std::isfinite(tex) || tex < (tb0 + 50.0) || tex > 4000.0) {
    tex = 900.0;
  }
  if (!std::isfinite(tgb0) || std::abs(tgb0) < 1e-6) {
    tgb0 = 2.0;
  }
  double sigma = tgb0 / (tex - tb0);
  if (!std::isfinite(sigma) || sigma <= 1e-4) {
    sigma = 0.03;
  }

  const double z = alt2gph(in.glat_deg, in.alt_km);
  double t = 0.0;
  if (z <= 86.0) {
    t = 288.15 - 6.5 * z;
    if (t < 170.0) {
      t = 170.0;
    }
  } else if (z < kZetaB) {
    const double alpha = (z - 86.0) / (kZetaB - 86.0);
    const double t86 = 288.15 - 6.5 * 86.0;
    t = (1.0 - alpha) * t86 + alpha * tb0;
  } else {
    t = tex - (tex - tb0) * std::exp(-sigma * (z - kZetaB));
  }
  if (!std::isfinite(t) || t < 120.0) {
    t = 120.0;
  }

  const double h_km = (kKb * t) / (kMbar * kG0) * 1e-3;
  const double pressure_pa = std::exp(kLnP0 - z / h_km);
  const double n_total_m3 = pressure_pa / (kKb * t);
  const double n_total_cm3 = n_total_m3 * 1e-6;

  double n2 = n_total_cm3 * std::exp(kLnVmr[1]);
  double o2 = n_total_cm3 * std::exp(kLnVmr[2]);
  double he = n_total_cm3 * std::exp(kLnVmr[4]) * std::exp((z - 80.0) / 220.0);
  double ar = n_total_cm3 * std::exp(kLnVmr[6]);

  double o = kDmissing;
  if (z >= 72.0) {
    o = n_total_cm3 * 1e-3 * std::exp((z - 72.0) / 20.0);
    o = std::min(o, n_total_cm3 * 0.2);
  }

  double h = kDmissing;
  if (z >= 75.0) {
    h = n_total_cm3 * 1e-6 * std::exp((z - 75.0) / 45.0);
  }

  double n = kDmissing;
  if (z >= 90.0) {
    n = n_total_cm3 * 1e-5 * std::exp((z - 90.0) / 35.0);
  }

  double o_anom = kDmissing;
  if (z >= 100.0 && o > 0.0) {
    o_anom = o * 0.02;
  }

  double no = kDmissing;
  if (z >= 70.0) {
    no = n_total_cm3 * 1e-7 * std::exp((z - 70.0) / 25.0);
  }

  auto n_cm3_to_m3 = [](double x) { return x <= kDmissing * 10.0 ? 0.0 : x * 1e6; };
  const double rho_kg_m3 =
      n_cm3_to_m3(n2) * kSpecMass[1] + n_cm3_to_m3(o2) * kSpecMass[2] +
      n_cm3_to_m3(o) * kSpecMass[3] + n_cm3_to_m3(he) * kSpecMass[4] +
      n_cm3_to_m3(h) * kSpecMass[5] + n_cm3_to_m3(ar) * kSpecMass[6] +
      n_cm3_to_m3(n) * kSpecMass[7] + n_cm3_to_m3(o_anom) * kSpecMass[8] +
      n_cm3_to_m3(no) * kSpecMass[9];
  const double rho_g_cm3 = rho_kg_m3 * 1e-3;

  Output out{};
  out.he = he;
  out.o = o;
  out.n2 = n2;
  out.o2 = o2;
  out.ar = ar;
  out.rho = rho_g_cm3;
  out.h = h;
  out.n = n;
  out.o_anom = o_anom;
  out.no = no;
  out.t = t;
  return CalcResult{.status = Status::Ok, .out = out};
}

}  // namespace msis21::detail
