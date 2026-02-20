/**
 * @file calc.cpp
 * @brief MSISCALC evaluation plumbing implementation.
 * @author Watosn
 */

#include "msis21/detail/calc.hpp"

#include <array>
#include <cmath>

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

}  // namespace

CalcResult evaluate_msiscalc(const Input& in, const Options& options, const Parameters& parameters) {
  if (!validate_input(in)) {
    return CalcResult{.status = Status::InvalidInput};
  }
  if (parameters.rows != kMaxBasisFunctions || parameters.cols != kParmColumnCount) {
    return CalcResult{.status = Status::ParmFileFormatError};
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

  const double z = in.alt_km;
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
