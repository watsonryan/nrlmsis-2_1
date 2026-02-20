/**
 * @file dfn.cpp
 * @brief Density profile helper implementations from msis_dfn.F90.
 * @author Watosn
 */

#include "msis21/detail/dfn.hpp"

#include <cmath>

#include "msis21/detail/gfn.hpp"
#include "msis21/detail/utils.hpp"

namespace msis21::detail {

namespace {

inline constexpr int kN2Start = 25;
inline constexpr int kO2Start = 35;
inline constexpr int kO1Start = 45;
inline constexpr int kHeStart = 63;
inline constexpr int kH1Start = 73;
inline constexpr int kArStart = 83;
inline constexpr int kN1Start = 93;
inline constexpr int kOaStart = 103;
inline constexpr int kNoStart = 113;

double subset_beta_at(const Parameters& parameters, int subset_start, int row, int col_local) {
  const int col = subset_start + col_local;
  return parameters.beta[static_cast<std::size_t>(row) +
                         static_cast<std::size_t>(col) * parameters.rows];
}

double subset_linear_dot(const Parameters& parameters,
                         int subset_start,
                         int col_local,
                         const std::array<double, kMaxBasisFunctions>& gf) {
  double sum = 0.0;
  for (int row = 0; row <= kMbf; ++row) {
    sum += subset_beta_at(parameters, subset_start, row, col_local) * gf[static_cast<std::size_t>(row)];
  }
  return sum;
}

std::array<double, kMaxBasisFunctions> subset_column_beta(const Parameters& parameters,
                                                          int subset_start,
                                                          int col_local) {
  std::array<double, kMaxBasisFunctions> out{};
  for (int row = 0; row < static_cast<int>(kMaxBasisFunctions); ++row) {
    out[static_cast<std::size_t>(row)] = subset_beta_at(parameters, subset_start, row, col_local);
  }
  return out;
}

std::array<double, kNmag> subset_geomag_params(const Parameters& parameters, int subset_start, int col_local) {
  std::array<double, kNmag> p{};
  for (int i = 0; i < kNmag; ++i) {
    p[static_cast<std::size_t>(i)] = subset_beta_at(parameters, subset_start, kCmag + i, col_local);
  }
  return p;
}

std::array<double, kNut> subset_ut_params(const Parameters& parameters, int subset_start, int col_local) {
  std::array<double, kNut> p{};
  for (int i = 0; i < kNut; ++i) {
    p[static_cast<std::size_t>(i)] = subset_beta_at(parameters, subset_start, kCut + i, col_local);
  }
  return p;
}

double sflux_dffact_for_subset(const Parameters& parameters, int subset_start, int col_local) {
  const double b0 = subset_beta_at(parameters, subset_start, 0, col_local);
  if (!std::isfinite(b0) || std::abs(b0) < 1e-12) {
    return 0.0;
  }
  return 1.0 / b0;
}

std::array<double, 13> geomag_bf_slice(const std::array<double, kMaxBasisFunctions>& basis) {
  std::array<double, 13> out{};
  for (int i = 0; i < 13; ++i) {
    out[static_cast<std::size_t>(i)] = basis[static_cast<std::size_t>(kCmag + i)];
  }
  return out;
}

std::array<double, 14> geomag_plg_slice(const std::array<double, kMaxBasisFunctions>& basis) {
  std::array<double, 14> out{};
  for (int i = 0; i < 14; ++i) {
    out[static_cast<std::size_t>(i)] = basis[static_cast<std::size_t>(kCmag + 13 + i)];
  }
  return out;
}

std::array<double, 9> ut_bf_slice(const std::array<double, kMaxBasisFunctions>& basis) {
  std::array<double, 9> out{};
  for (int i = 0; i < 9; ++i) {
    out[static_cast<std::size_t>(i)] = basis[static_cast<std::size_t>(kCut + i)];
  }
  return out;
}

double spline_sum_rel(const std::array<double, kNl + 1>& coeff, const BsplineResult& b, int order, int rel_min) {
  double sum = 0.0;
  for (int rel = rel_min; rel <= 0; ++rel) {
    const int idx = b.i + rel;
    if (idx >= 0 && idx <= kNl) {
      sum += coeff[static_cast<std::size_t>(idx)] * b.value(rel, order);
    }
  }
  return sum;
}

double wz_at(double z, const TnParm& tpro, const EtaTable& eta_tn) {
  const double delz = z - kZetaB;
  if (z < kZetaB) {
    const auto b = bspline(z, kNodesTN, kNd + 2, 6, eta_tn);
    return spline_sum_rel(tpro.gamma, b, 6, -5) + tpro.cvs * delz + tpro.cws;
  }
  return (0.5 * delz * delz + dilog(tpro.b * std::exp(-tpro.sigma * delz)) / tpro.sigmasq) / tpro.tex +
         tpro.cvb * delz + tpro.cwb;
}

EtaTable make_eta_tn() {
  EtaTable eta{};
  for (int k = 2; k <= 6; ++k) {
    for (int j = 0; j <= kNl; ++j) {
      eta[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] =
          1.0 / (kNodesTN[static_cast<std::size_t>(j + k - 1)] - kNodesTN[static_cast<std::size_t>(j)]);
    }
  }
  return eta;
}

}  // namespace

double pwmp(double z,
            const std::array<double, 5>& zm,
            const std::array<double, 5>& m,
            const std::array<double, 5>& dmdz) {
  if (z >= zm[4]) {
    return m[4];
  }
  if (z <= zm[0]) {
    return m[0];
  }

  for (int inode = 0; inode <= 3; ++inode) {
    if (z < zm[static_cast<std::size_t>(inode + 1)]) {
      return m[static_cast<std::size_t>(inode)] +
             dmdz[static_cast<std::size_t>(inode)] *
                 (z - zm[static_cast<std::size_t>(inode)]);
    }
  }

  return m[4];
}

double dfnx(double z,
            double tnz,
            double lndtotz,
            double vz,
            double wz,
            double hrfact,
            const TnParm& /*tpro*/,
            const DnParm& dpro) {
  if (z < dpro.zmin) {
    return kDmissing;
  }

  if (dpro.ispec == 9) {
    const double ln_d =
        dpro.lndref - (z - dpro.zref) / kHoa - dpro.c * std::exp(-(z - dpro.zeta_c) / dpro.hc);
    return std::exp(ln_d);
  }

  if (dpro.ispec == 10 && dpro.lndref == 0.0) {
    return kDmissing;
  }

  double ccor = 0.0;
  switch (dpro.ispec) {
    case 2:
    case 3:
    case 5:
    case 7:
      ccor = dpro.r * (1.0 + std::tanh((z - dpro.zeta_r) / (hrfact * dpro.hr)));
      break;
    case 4:
    case 6:
    case 8:
    case 10:
      ccor = -dpro.c * std::exp(-(z - dpro.zeta_c) / dpro.hc) +
             dpro.r * (1.0 + std::tanh((z - dpro.zeta_r) / (hrfact * dpro.hr)));
      break;
    default:
      break;
  }

  if (z < dpro.zhyd) {
    switch (dpro.ispec) {
      case 2:
      case 3:
      case 5:
      case 7:
        return std::exp(lndtotz + dpro.lnphi_f + ccor);
      default:
        return kDmissing;
    }
  }

  const double mz = pwmp(z, dpro.zeta_mi, dpro.mi, dpro.ami);
  double ihyd = mz * vz - dpro.izref;
  if (z > dpro.zeta_mi[0] && z < dpro.zeta_mi[4]) {
    int i = 0;
    for (int i1 = 1; i1 <= 3; ++i1) {
      if (z < dpro.zeta_mi[static_cast<std::size_t>(i1)]) {
        break;
      }
      i = i1;
    }
    ihyd -= (dpro.ami[static_cast<std::size_t>(i)] * wz + dpro.xmi[static_cast<std::size_t>(i)]);
  } else if (z >= dpro.zeta_mi[4]) {
    ihyd -= dpro.xmi[4];
  }

  double ln_d = dpro.lndref - ihyd * kG0DivKb + ccor;
  ln_d = std::exp(ln_d) * dpro.tref / tnz;
  return ln_d;
}

DnParm dfnparm(int ispec,
               const std::array<double, kMaxBasisFunctions>& gf,
               const std::array<bool, kMaxBasisFunctions>& swg,
               const Parameters& parameters,
               const TnParm& tpro) {
  DnParm dpro{};
  dpro.ispec = ispec;
  static const EtaTable eta_tn = make_eta_tn();
  const auto bf_mag = geomag_bf_slice(gf);
  const auto plg_mag = geomag_plg_slice(gf);
  const auto bf_ut = ut_bf_slice(gf);
  std::array<bool, kNmag> swg_mag{};
  std::array<bool, kNut> swg_ut{};
  for (int i = 0; i < kNmag; ++i) {
    swg_mag[static_cast<std::size_t>(i)] = swg[static_cast<std::size_t>(kCmag + i)];
  }
  for (int i = 0; i < kNut; ++i) {
    swg_ut[static_cast<std::size_t>(i)] = swg[static_cast<std::size_t>(kCut + i)];
  }

  auto dyn_terms = [&](int subset_start, int col_local, double sflux_fact) {
    double r = subset_linear_dot(parameters, subset_start, col_local, gf);
    r += sfluxmod(col_local, gf, subset_column_beta(parameters, subset_start, col_local), sflux_fact, swg);
    r += geomag(subset_geomag_params(parameters, subset_start, col_local), bf_mag, plg_mag, swg_mag);
    r += utdep(subset_ut_params(parameters, subset_start, col_local), bf_ut, swg_ut);
    return r;
  };

  switch (ispec) {
    case 2:
      dpro.lnphi_f = kLnVmr[2 - 1];
      dpro.lndref = tpro.lndtot_f + dpro.lnphi_f;
      dpro.zref = kZetaF;
      dpro.zmin = -1.0;
      dpro.zhyd = kZetaF;
      dpro.zeta_m = subset_linear_dot(parameters, kN2Start, 1, gf);
      dpro.hml = subset_beta_at(parameters, kN2Start, 0, 2);
      dpro.hmu = subset_beta_at(parameters, kN2Start, 0, 3);
      dpro.r = 0.0;
      dpro.zeta_r = subset_beta_at(parameters, kN2Start, 0, 8);
      dpro.hr = subset_beta_at(parameters, kN2Start, 0, 9);
      break;
    case 3:
      dpro.lnphi_f = kLnVmr[3 - 1];
      dpro.lndref = tpro.lndtot_f + dpro.lnphi_f;
      dpro.zref = kZetaF;
      dpro.zmin = -1.0;
      dpro.zhyd = kZetaF;
      dpro.zeta_m = subset_beta_at(parameters, kO2Start, 0, 1);
      dpro.hml = subset_beta_at(parameters, kO2Start, 0, 2);
      dpro.hmu = subset_beta_at(parameters, kO2Start, 0, 3);
      dpro.r = subset_linear_dot(parameters, kO2Start, 7, gf);
      dpro.r += geomag(subset_geomag_params(parameters, kO2Start, 7), bf_mag, plg_mag, swg_mag);
      dpro.zeta_r = subset_beta_at(parameters, kO2Start, 0, 8);
      dpro.hr = subset_beta_at(parameters, kO2Start, 0, 9);
      break;
    case 4:
      dpro.lndref = subset_linear_dot(parameters, kO1Start, 0, gf);
      dpro.zref = kZetaA;
      dpro.zmin = 50.0;
      dpro.zhyd = kZetaA;
      dpro.zeta_m = subset_beta_at(parameters, kO1Start, 0, 1);
      dpro.hml = subset_beta_at(parameters, kO1Start, 0, 2);
      dpro.hmu = subset_beta_at(parameters, kO1Start, 0, 3);
      dpro.c = subset_linear_dot(parameters, kO1Start, 4, gf);
      dpro.zeta_c = subset_beta_at(parameters, kO1Start, 0, 5);
      dpro.hc = subset_beta_at(parameters, kO1Start, 0, 6);
      dpro.r = dyn_terms(kO1Start, 7, 0.0);
      dpro.zeta_r = subset_beta_at(parameters, kO1Start, 0, 8);
      dpro.hr = subset_beta_at(parameters, kO1Start, 0, 9);
      break;
    case 5:
      dpro.lnphi_f = kLnVmr[5 - 1];
      dpro.lndref = tpro.lndtot_f + dpro.lnphi_f;
      dpro.zref = kZetaF;
      dpro.zmin = -1.0;
      dpro.zhyd = kZetaF;
      dpro.zeta_m = subset_beta_at(parameters, kHeStart, 0, 1);
      dpro.hml = subset_beta_at(parameters, kHeStart, 0, 2);
      dpro.hmu = subset_beta_at(parameters, kHeStart, 0, 3);
      dpro.r = dyn_terms(kHeStart, 7, sflux_dffact_for_subset(parameters, kHeStart, 7));
      dpro.zeta_r = subset_beta_at(parameters, kHeStart, 0, 8);
      dpro.hr = subset_beta_at(parameters, kHeStart, 0, 9);
      break;
    case 6:
      dpro.lndref = subset_linear_dot(parameters, kH1Start, 0, gf);
      dpro.zref = kZetaA;
      dpro.zmin = 75.0;
      dpro.zhyd = kZetaF;
      dpro.zeta_m = subset_beta_at(parameters, kH1Start, 0, 1);
      dpro.hml = subset_beta_at(parameters, kH1Start, 0, 2);
      dpro.hmu = subset_beta_at(parameters, kH1Start, 0, 3);
      dpro.c = subset_linear_dot(parameters, kH1Start, 4, gf);
      dpro.zeta_c = subset_linear_dot(parameters, kH1Start, 5, gf);
      dpro.hc = subset_beta_at(parameters, kH1Start, 0, 6);
      dpro.r = dyn_terms(kH1Start, 7, 0.0);
      dpro.zeta_r = subset_beta_at(parameters, kH1Start, 0, 8);
      dpro.hr = subset_beta_at(parameters, kH1Start, 0, 9);
      break;
    case 7:
      dpro.lnphi_f = kLnVmr[7 - 1];
      dpro.lndref = tpro.lndtot_f + dpro.lnphi_f;
      dpro.zref = kZetaF;
      dpro.zmin = -1.0;
      dpro.zhyd = kZetaF;
      dpro.zeta_m = subset_beta_at(parameters, kArStart, 0, 1);
      dpro.hml = subset_beta_at(parameters, kArStart, 0, 2);
      dpro.hmu = subset_beta_at(parameters, kArStart, 0, 3);
      dpro.r = subset_linear_dot(parameters, kArStart, 7, gf);
      dpro.r += geomag(subset_geomag_params(parameters, kArStart, 7), bf_mag, plg_mag, swg_mag);
      dpro.r += utdep(subset_ut_params(parameters, kArStart, 7), bf_ut, swg_ut);
      dpro.zeta_r = subset_beta_at(parameters, kArStart, 0, 8);
      dpro.hr = subset_beta_at(parameters, kArStart, 0, 9);
      break;
    case 8:
      dpro.lndref = subset_linear_dot(parameters, kN1Start, 0, gf);
      dpro.lndref += sfluxmod(0, gf, subset_column_beta(parameters, kN1Start, 0), 0.0, swg);
      dpro.lndref += geomag(subset_geomag_params(parameters, kN1Start, 0), bf_mag, plg_mag, swg_mag);
      dpro.lndref += utdep(subset_ut_params(parameters, kN1Start, 0), bf_ut, swg_ut);
      dpro.zref = kZetaB;
      dpro.zmin = 90.0;
      dpro.zhyd = kZetaF;
      dpro.zeta_m = subset_beta_at(parameters, kN1Start, 0, 1);
      dpro.hml = subset_beta_at(parameters, kN1Start, 0, 2);
      dpro.hmu = subset_beta_at(parameters, kN1Start, 0, 3);
      dpro.c = subset_beta_at(parameters, kN1Start, 0, 4);
      dpro.zeta_c = subset_beta_at(parameters, kN1Start, 0, 5);
      dpro.hc = subset_beta_at(parameters, kN1Start, 0, 6);
      dpro.r = subset_linear_dot(parameters, kN1Start, 7, gf);
      dpro.zeta_r = subset_beta_at(parameters, kN1Start, 0, 8);
      dpro.hr = subset_beta_at(parameters, kN1Start, 0, 9);
      break;
    case 9:
      dpro.lndref = subset_linear_dot(parameters, kOaStart, 0, gf);
      dpro.lndref += geomag(subset_geomag_params(parameters, kOaStart, 0), bf_mag, plg_mag, swg_mag);
      dpro.zref = kZetaB;
      dpro.zmin = 120.0;
      dpro.zhyd = 0.0;
      dpro.c = subset_beta_at(parameters, kOaStart, 0, 4);
      dpro.zeta_c = subset_beta_at(parameters, kOaStart, 0, 5);
      dpro.hc = subset_beta_at(parameters, kOaStart, 0, 6);
      return dpro;
    case 10:
      if (subset_beta_at(parameters, kNoStart, 0, 0) == 0.0) {
        dpro.lndref = 0.0;
        return dpro;
      }
      dpro.lndref = subset_linear_dot(parameters, kNoStart, 0, gf);
      dpro.lndref += geomag(subset_geomag_params(parameters, kNoStart, 0), bf_mag, plg_mag, swg_mag);
      dpro.zref = kZetaB;
      dpro.zmin = 72.5;
      dpro.zhyd = kZetaB;
      dpro.zeta_m = subset_linear_dot(parameters, kNoStart, 1, gf);
      dpro.hml = subset_linear_dot(parameters, kNoStart, 2, gf);
      dpro.hmu = subset_linear_dot(parameters, kNoStart, 3, gf);
      dpro.c = subset_linear_dot(parameters, kNoStart, 4, gf);
      dpro.c += geomag(subset_geomag_params(parameters, kNoStart, 4), bf_mag, plg_mag, swg_mag);
      dpro.zeta_c = subset_linear_dot(parameters, kNoStart, 5, gf);
      dpro.hc = subset_linear_dot(parameters, kNoStart, 6, gf);
      dpro.r = subset_linear_dot(parameters, kNoStart, 7, gf);
      dpro.zeta_r = subset_linear_dot(parameters, kNoStart, 8, gf);
      dpro.hr = subset_linear_dot(parameters, kNoStart, 9, gf);
      break;
    default:
      return dpro;
  }

  dpro.zeta_mi[0] = dpro.zeta_m - 2.0 * dpro.hml;
  dpro.zeta_mi[1] = dpro.zeta_m - dpro.hml;
  dpro.zeta_mi[2] = dpro.zeta_m;
  dpro.zeta_mi[3] = dpro.zeta_m + dpro.hmu;
  dpro.zeta_mi[4] = dpro.zeta_m + 2.0 * dpro.hmu;
  dpro.mi[0] = kMbar;
  dpro.mi[4] = kSpecMass[static_cast<std::size_t>(ispec)];
  dpro.mi[2] = (dpro.mi[0] + dpro.mi[4]) / 2.0;
  const double delm = std::tanh(1.0) * (dpro.mi[4] - dpro.mi[0]) / 2.0;
  dpro.mi[1] = dpro.mi[2] - delm;
  dpro.mi[3] = dpro.mi[2] + delm;
  for (int i = 0; i <= 3; ++i) {
    dpro.ami[static_cast<std::size_t>(i)] =
        (dpro.mi[static_cast<std::size_t>(i + 1)] - dpro.mi[static_cast<std::size_t>(i)]) /
        (dpro.zeta_mi[static_cast<std::size_t>(i + 1)] - dpro.zeta_mi[static_cast<std::size_t>(i)]);
  }
  for (int i = 0; i <= 4; ++i) {
    dpro.wmi[static_cast<std::size_t>(i)] = wz_at(dpro.zeta_mi[static_cast<std::size_t>(i)], tpro, eta_tn);
  }
  dpro.xmi[0] = -dpro.ami[0] * dpro.wmi[0];
  for (int i = 1; i <= 3; ++i) {
    dpro.xmi[static_cast<std::size_t>(i)] =
        dpro.xmi[static_cast<std::size_t>(i - 1)] -
        dpro.wmi[static_cast<std::size_t>(i)] * (dpro.ami[static_cast<std::size_t>(i)] - dpro.ami[static_cast<std::size_t>(i - 1)]);
  }
  dpro.xmi[4] = dpro.xmi[3] + dpro.wmi[4] * dpro.ami[3];

  if (dpro.zref == kZetaF) {
    dpro.tref = tpro.tzeta_f;
    dpro.izref = kMbar * tpro.vzeta_f;
  } else if (dpro.zref == kZetaB) {
    dpro.tref = tpro.tb0;
    dpro.izref = 0.0;
    if (kZetaB > dpro.zeta_mi[0] && kZetaB < dpro.zeta_mi[4]) {
      int i = 0;
      for (int i1 = 1; i1 <= 3; ++i1) {
        if (kZetaB < dpro.zeta_mi[static_cast<std::size_t>(i1)]) {
          break;
        }
        i = i1;
      }
      dpro.izref -= dpro.xmi[static_cast<std::size_t>(i)];
    } else {
      dpro.izref -= dpro.xmi[4];
    }
  } else if (dpro.zref == kZetaA) {
    const double mzref = pwmp(dpro.zref, dpro.zeta_mi, dpro.mi, dpro.ami);
    dpro.tref = tpro.tzeta_a;
    dpro.izref = mzref * tpro.vzeta_a;
    if (kZetaA > dpro.zeta_mi[0] && kZetaA < dpro.zeta_mi[4]) {
      int i = 0;
      for (int i1 = 1; i1 <= 3; ++i1) {
        if (kZetaA < dpro.zeta_mi[static_cast<std::size_t>(i1)]) {
          break;
        }
        i = i1;
      }
      dpro.izref -= (dpro.ami[static_cast<std::size_t>(i)] * tpro.wzeta_a + dpro.xmi[static_cast<std::size_t>(i)]);
    } else {
      dpro.izref -= dpro.xmi[4];
    }
  }

  return dpro;
}

}  // namespace msis21::detail
