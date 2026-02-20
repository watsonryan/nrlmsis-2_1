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
#include "msis21/detail/dfn.hpp"
#include "msis21/detail/gfn.hpp"
#include "msis21/detail/log.hpp"
#include "msis21/detail/parm_reader.hpp"
#include "msis21/detail/tfn.hpp"
#include "msis21/detail/utils.hpp"
#include <spdlog/spdlog.h>

namespace msis21::detail {

namespace {

#ifdef MSIS21_USE_FORTRAN_REFERENCE_BACKEND
extern "C" void msis21_fortran_eval(int iyd,
                                    float sec,
                                    float alt,
                                    float glat,
                                    float glon,
                                    float f107a,
                                    float f107,
                                    float apd,
                                    double* tn,
                                    double dn[10],
                                    int* status);
#endif

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

[[maybe_unused]] EtaTable make_eta_tn() {
  EtaTable eta{};
  for (int k = 2; k <= 6; ++k) {
    for (int j = 0; j <= kNl; ++j) {
      eta[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] =
          1.0 / (kNodesTN[static_cast<std::size_t>(j + k - 1)] - kNodesTN[static_cast<std::size_t>(j)]);
    }
  }
  return eta;
}

[[maybe_unused]] double spline_sum_rel(const std::array<double, kNl + 1>& coeff,
                                       const BsplineResult& b,
                                       int order,
                                       int rel_min) {
  double sum = 0.0;
  for (int rel = rel_min; rel <= 0; ++rel) {
    const int idx = b.i + rel;
    if (idx >= 0 && idx <= kNl) {
      sum += coeff[static_cast<std::size_t>(idx)] * b.value(rel, order);
    }
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

[[maybe_unused]] std::array<bool, kMaxBasisFunctions> legacy_switches_to_basis(const Options& options) {
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
    set_range(swg, 205, 208, false);
    set_range(swg, 213, 216, false);
    set_range(swg, 259, 262, false);
    set_range(swg, 267, 270, false);
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
  if (parameters.rows != kMaxBasisFunctions || parameters.cols != kParmColumnCount) {
    return CalcResult{.status = Status::ParmFileFormatError};
  }

#ifdef MSIS21_USE_FORTRAN_REFERENCE_BACKEND
  const float sec = static_cast<float>(in.sec);
  const float alt = static_cast<float>(in.alt_km);
  const float glat = static_cast<float>(in.glat_deg);
  const float glon = static_cast<float>(in.glon_deg);
  const float f107a = static_cast<float>(in.f107a);
  const float f107 = static_cast<float>(in.f107);
  const float apd = static_cast<float>(in.ap);
  double tn = 0.0;
  double dn_m3[10]{};
  int fstatus = 0;
  msis21_fortran_eval(in.iyd, sec, alt, glat, glon, f107a, f107, apd, &tn, dn_m3, &fstatus);
  if (fstatus != 0) {
    return CalcResult{.status = Status::NumericalError};
  }
  auto m3_to_cm3 = [](double x) { return x == kDmissing ? kDmissing : x * 1.0e-6; };
  Output out{};
  out.he = m3_to_cm3(dn_m3[4]);
  out.o = m3_to_cm3(dn_m3[3]);
  out.n2 = m3_to_cm3(dn_m3[1]);
  out.o2 = m3_to_cm3(dn_m3[2]);
  out.ar = m3_to_cm3(dn_m3[6]);
  out.rho = (dn_m3[0] == kDmissing ? kDmissing : dn_m3[0] * 1.0e-3);
  out.h = m3_to_cm3(dn_m3[5]);
  out.n = m3_to_cm3(dn_m3[7]);
  out.o_anom = m3_to_cm3(dn_m3[8]);
  out.no = m3_to_cm3(dn_m3[9]);
  out.t = tn;
  (void)options;
  return CalcResult{.status = Status::Ok, .out = out};
#else

  GlobeCalculator globe;
  const auto f32 = [](double x) { return static_cast<double>(static_cast<float>(x)); };
  GlobeInput globe_input;
  globe_input.doy = static_cast<double>(in.iyd % 1000);
  globe_input.utsec = f32(in.sec);
  globe_input.lat = f32(in.glat_deg);
  globe_input.lon = f32(in.glon_deg);
  globe_input.sfluxavg = f32(in.f107a);
  globe_input.sflux = f32(in.f107);
  const double apd = f32(in.ap);
  globe_input.ap = {apd, apd, apd, apd, apd, apd, apd};

  const auto switches = legacy_switches_to_basis(options);

  const auto basis = globe.globe(globe_input, switches);
  const auto tpro = tfnparm(basis, switches, parameters);

  const double z = alt2gph(f32(in.glat_deg), f32(in.alt_km));
  static const EtaTable eta_tn = make_eta_tn();
  double t = 0.0;
  BsplineResult bs{};
  bool have_bs = false;
  if (z < kZetaB) {
    const int kmax = (z < kZetaF) ? 5 : 6;
    bs = bspline(z, kNodesTN, kNd + 2, kmax, eta_tn);
    have_bs = true;
    std::array<double, 4> w{};
    w[0] = bs.value(-3, 4);
    w[1] = bs.value(-2, 4);
    w[2] = bs.value(-1, 4);
    w[3] = bs.value(0, 4);
    t = tfnx(z, bs.i, w, tpro);
  } else {
    t = tpro.tex - (tpro.tex - tpro.tb0) * std::exp(-tpro.sigma * (z - kZetaB));
  }

  const double delz = z - kZetaB;
  double vz = 0.0;
  double wz = 0.0;
  double lndtotz = 0.0;
  if (z < kZetaF) {
    if (!have_bs) {
      bs = bspline(z, kNodesTN, kNd + 2, 5, eta_tn);
    }
    vz = spline_sum_rel(tpro.beta, bs, 5, -4) + tpro.cvs;
    lndtotz = kLnP0 - kMbarG0DivKb * (vz - tpro.vzeta0) - std::log(kKb * t);
  } else {
    if (z < kZetaB) {
      if (!have_bs) {
        bs = bspline(z, kNodesTN, kNd + 2, 6, eta_tn);
      }
      vz = spline_sum_rel(tpro.beta, bs, 5, -4) + tpro.cvs;
      wz = spline_sum_rel(tpro.gamma, bs, 6, -5) + tpro.cvs * delz + tpro.cws;
    } else {
      vz = (delz + std::log(t / tpro.tex) / tpro.sigma) / tpro.tex + tpro.cvb;
      wz = (0.5 * delz * delz + dilog(tpro.b * std::exp(-tpro.sigma * delz)) / tpro.sigmasq) / tpro.tex +
           tpro.cvb * delz + tpro.cwb;
    }
    const double lnpz = kLnP0 - kMbarG0DivKb * (vz - tpro.vzeta0);
    lndtotz = lnpz - std::log(kKb * t);
  }
  const double hrfact = 0.5 * (1.0 + std::tanh(kHgamma * (z - kZetagamma)));
  const double n2_m3 = dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(2, basis, switches, parameters, tpro));
  const double o2_m3 = dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(3, basis, switches, parameters, tpro));
  const double o_m3 = dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(4, basis, switches, parameters, tpro));
  const double he_m3 = dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(5, basis, switches, parameters, tpro));
  const double h_m3 = dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(6, basis, switches, parameters, tpro));
  const double ar_m3 = dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(7, basis, switches, parameters, tpro));
  const double n_m3 = dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(8, basis, switches, parameters, tpro));
  const double o_anom_m3 =
      dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(9, basis, switches, parameters, tpro));
  const double no_m3 =
      dfnx(z, t, lndtotz, vz, wz, hrfact, tpro, dfnparm(10, basis, switches, parameters, tpro));

  auto safe_m3 = [](double x) { return x == kDmissing ? 0.0 : x; };
  const double rho_kg_m3 =
      safe_m3(n2_m3) * kSpecMass[1] + safe_m3(o2_m3) * kSpecMass[2] + safe_m3(o_m3) * kSpecMass[3] +
      safe_m3(he_m3) * kSpecMass[4] + safe_m3(h_m3) * kSpecMass[5] + safe_m3(ar_m3) * kSpecMass[6] +
      safe_m3(n_m3) * kSpecMass[7] + safe_m3(o_anom_m3) * kSpecMass[8];
  const double rho_g_cm3 = rho_kg_m3 * 1e-3;
  auto m3_to_cm3 = [](double x) { return x == kDmissing ? kDmissing : x * 1.0e-6; };
  const double he = m3_to_cm3(he_m3);
  const double o = m3_to_cm3(o_m3);
  const double n2 = m3_to_cm3(n2_m3);
  const double o2 = m3_to_cm3(o2_m3);
  const double ar = m3_to_cm3(ar_m3);
  const double h = m3_to_cm3(h_m3);
  const double n = m3_to_cm3(n_m3);
  const double o_anom = m3_to_cm3(o_anom_m3);
  const double no = m3_to_cm3(no_m3);

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
#endif
}

}  // namespace msis21::detail
