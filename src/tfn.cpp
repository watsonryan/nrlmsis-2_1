/**
 * @file tfn.cpp
 * @brief Vertical temperature profile utilities from msis_tfn.F90.
 * @author Watosn
 */

#include "msis21/detail/tfn.hpp"

#include <algorithm>
#include <cmath>

#include "msis21/detail/gfn.hpp"
#include "msis21/detail/utils.hpp"

namespace msis21::detail {

namespace {

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

std::array<double, kMaxBasisFunctions> column_beta(const Parameters& parameters, int col) {
  std::array<double, kMaxBasisFunctions> out{};
  for (int row = 0; row < static_cast<int>(kMaxBasisFunctions); ++row) {
    out[static_cast<std::size_t>(row)] = beta_at(parameters, row, col);
  }
  return out;
}

bool smod_enabled(const Parameters& parameters, int col) {
  return beta_at(parameters, kCsfxMod, col) != 0.0 || beta_at(parameters, kCsfxMod + 1, col) != 0.0 ||
         beta_at(parameters, kCsfxMod + 2, col) != 0.0;
}

double dot3(const std::array<double, kNl + 1>& v, int start, const std::array<double, 3>& w) {
  return v[static_cast<std::size_t>(start)] * w[0] + v[static_cast<std::size_t>(start + 1)] * w[1] +
         v[static_cast<std::size_t>(start + 2)] * w[2];
}

double dot4(const std::array<double, kNl + 1>& v, int start, const std::array<double, 4>& w) {
  return v[static_cast<std::size_t>(start)] * w[0] + v[static_cast<std::size_t>(start + 1)] * w[1] +
         v[static_cast<std::size_t>(start + 2)] * w[2] + v[static_cast<std::size_t>(start + 3)] * w[3];
}

double dot5(const std::array<double, kNl + 1>& v, int start, const std::array<double, 5>& w) {
  return v[static_cast<std::size_t>(start)] * w[0] + v[static_cast<std::size_t>(start + 1)] * w[1] +
         v[static_cast<std::size_t>(start + 2)] * w[2] + v[static_cast<std::size_t>(start + 3)] * w[3] +
         v[static_cast<std::size_t>(start + 4)] * w[4];
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

std::array<double, kNmag> geomag_params_for_col(const Parameters& parameters, int col) {
  std::array<double, kNmag> out{};
  for (int i = 0; i < kNmag; ++i) {
    out[static_cast<std::size_t>(i)] = beta_at(parameters, kCmag + i, col);
  }
  return out;
}

std::array<double, kNut> ut_params_for_col(const Parameters& parameters, int col) {
  std::array<double, kNut> out{};
  for (int i = 0; i < kNut; ++i) {
    out[static_cast<std::size_t>(i)] = beta_at(parameters, kCut + i, col);
  }
  return out;
}

void fill_wbeta_wgamma(std::array<double, kNl + 1>& wbeta, std::array<double, kNl + 1>& wgamma) {
  for (int i = 0; i <= kNl; ++i) {
    wbeta[static_cast<std::size_t>(i)] =
        (kNodesTN[static_cast<std::size_t>(i + 4)] - kNodesTN[static_cast<std::size_t>(i)]) / 4.0;
    wgamma[static_cast<std::size_t>(i)] =
        (kNodesTN[static_cast<std::size_t>(i + 5)] - kNodesTN[static_cast<std::size_t>(i)]) / 5.0;
  }
}

}  // namespace

std::array<double, 3> continuity_coefficients(double tex, double tgb0, double tb0) {
  const double sigma = tgb0 / (tex - tb0);
  std::array<double, 3> bc{};
  bc[0] = 1.0 / tb0;
  bc[1] = -tgb0 / (tb0 * tb0);
  bc[2] = -bc[1] * (sigma + 2.0 * tgb0 / tb0);
  return bc;
}

double tfnx(double z, int iz, const std::array<double, 4>& wght, const TnParm& tpro) {
  if (z < kZetaB) {
    const int i = std::max(iz - 3, 0);
    const int j = (iz < 3) ? (-iz) : -3;
    double sum = 0.0;
    for (int k = i; k <= iz; ++k) {
      sum += tpro.cf[static_cast<std::size_t>(k)] * wght[static_cast<std::size_t>(k - iz - j)];
    }
    return 1.0 / sum;
  }

  return tpro.tex - (tpro.tex - tpro.tb0) * std::exp(-tpro.sigma * (z - kZetaB));
}

TnParm tfnparm(const std::array<double, kMaxBasisFunctions>& gf,
               const std::array<bool, kMaxBasisFunctions>& swg,
               const Parameters& parameters) {
  TnParm tpro{};
  std::array<double, kNl + 1> wbeta{};
  std::array<double, kNl + 1> wgamma{};
  fill_wbeta_wgamma(wbeta, wgamma);

  for (int ix = 0; ix <= (kItb0 - 1); ++ix) {
    tpro.cf[static_cast<std::size_t>(ix)] = linear_dot(parameters, gf, ix);
  }
  for (int ix = 0; ix <= (kItb0 - 1); ++ix) {
    if (smod_enabled(parameters, ix)) {
      const auto beta_col = column_beta(parameters, ix);
      tpro.cf[static_cast<std::size_t>(ix)] +=
          sfluxmod(ix, gf, beta_col, 1.0 / beta_at(parameters, 0, ix), swg);
    }
  }

  tpro.tex = linear_dot(parameters, gf, kItex);
  tpro.tgb0 = linear_dot(parameters, gf, kItgb0);
  tpro.tb0 = linear_dot(parameters, gf, kItb0);

  const auto beta_tex = column_beta(parameters, kItex);
  const auto beta_tgb0 = column_beta(parameters, kItgb0);
  const auto beta_tb0 = column_beta(parameters, kItb0);
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

  tpro.tex += sfluxmod(kItex, gf, beta_tex, 1.0 / beta_at(parameters, 0, kItex), swg);
  tpro.tex += geomag(geomag_params_for_col(parameters, kItex), bf_mag, plg_mag, swg_mag);
  tpro.tex += utdep(ut_params_for_col(parameters, kItex), bf_ut, swg_ut);

  if (smod_enabled(parameters, kItgb0)) {
    tpro.tgb0 += sfluxmod(kItgb0, gf, beta_tgb0, 1.0 / beta_at(parameters, 0, kItgb0), swg);
  }
  tpro.tgb0 += geomag(geomag_params_for_col(parameters, kItgb0), bf_mag, plg_mag, swg_mag);

  if (smod_enabled(parameters, kItb0)) {
    tpro.tb0 += sfluxmod(kItb0, gf, beta_tb0, 1.0 / beta_at(parameters, 0, kItb0), swg);
  }
  tpro.tb0 += geomag(geomag_params_for_col(parameters, kItb0), bf_mag, plg_mag, swg_mag);
  tpro.sigma = tpro.tgb0 / (tpro.tex - tpro.tb0);

  const auto bc = continuity_coefficients(tpro.tex, tpro.tgb0, tpro.tb0);
  tpro.cf[static_cast<std::size_t>(kItb0)] = bc[0] - 10.0 * bc[1] + 33.333333333333336 * bc[2];
  tpro.cf[static_cast<std::size_t>(kItgb0)] = bc[0] - 16.666666666666668 * bc[2];
  tpro.cf[static_cast<std::size_t>(kItex)] = bc[0] + 10.0 * bc[1] + 33.333333333333336 * bc[2];

  tpro.tzeta_f = 1.0 / dot3(tpro.cf, kIzfx, kS4zetaF);
  tpro.tzeta_a = 1.0 / dot3(tpro.cf, kIzax, kS4zetaA);
  tpro.dlntdz_a = -dot3(tpro.cf, kIzax, kWghtAxdz) * tpro.tzeta_a;

  tpro.beta[0] = tpro.cf[0] * wbeta[0];
  for (int ix = 1; ix <= kNl; ++ix) {
    tpro.beta[static_cast<std::size_t>(ix)] =
        tpro.beta[static_cast<std::size_t>(ix - 1)] + tpro.cf[static_cast<std::size_t>(ix)] * wbeta[static_cast<std::size_t>(ix)];
  }
  tpro.gamma[0] = tpro.beta[0] * wgamma[0];
  for (int ix = 1; ix <= kNl; ++ix) {
    tpro.gamma[static_cast<std::size_t>(ix)] =
        tpro.gamma[static_cast<std::size_t>(ix - 1)] +
        tpro.beta[static_cast<std::size_t>(ix)] * wgamma[static_cast<std::size_t>(ix)];
  }

  tpro.b = 1.0 - tpro.tb0 / tpro.tex;
  tpro.sigmasq = tpro.sigma * tpro.sigma;
  tpro.cvs = -dot4(tpro.beta, kItb0 - 1, kS5zetaB);
  tpro.cws = -dot5(tpro.gamma, kItb0 - 2, kS6zetaB);
  tpro.cvb = -std::log(1.0 - tpro.b) / (tpro.sigma * tpro.tex);
  tpro.cwb = -dilog(tpro.b) / (tpro.sigmasq * tpro.tex);
  tpro.vzeta_f = dot4(tpro.beta, kIzfx - 1, kS5zetaF) + tpro.cvs;
  tpro.vzeta_a = dot4(tpro.beta, kIzax - 1, kS5zetaA) + tpro.cvs;
  tpro.wzeta_a = dot5(tpro.gamma, kIzax - 2, kS6zetaA) + tpro.cvs * (kZetaA - kZetaB) + tpro.cws;
  tpro.vzeta0 = dot3(tpro.beta, 0, kS5zeta0) + tpro.cvs;
  tpro.lndtot_f = kLnP0 - kMbarG0DivKb * (tpro.vzeta_f - tpro.vzeta0) - std::log(kKb * tpro.tzeta_f);

  return tpro;
}

}  // namespace msis21::detail
