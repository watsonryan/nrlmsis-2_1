/**
 * @file gfn.cpp
 * @brief Horizontal/time-dependent basis implementations from msis_gfn.F90.
 * @author Watosn
 */

#include "msis21/detail/gfn.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace msis21::detail {

GlobeCalculator::GlobeCalculator() {
  plg_.fill({});
  cdoy_.fill(0.0);
  sdoy_.fill(0.0);
  clst_.fill(0.0);
  slst_.fill(0.0);
  clon_.fill(0.0);
  slon_.fill(0.0);
}

std::array<double, kMaxBasisFunctions> GlobeCalculator::globe(
    const GlobeInput& in,
    const std::array<bool, kMaxBasisFunctions>& switches) {
  std::array<double, kMaxBasisFunctions> bf{};

  if (in.lat != lastlat_) {
    const double clat = std::sin(in.lat * kDeg2Rad);
    const double slat = std::cos(in.lat * kDeg2Rad);
    const double clat2 = clat * clat;
    const double clat4 = clat2 * clat2;
    const double slat2 = slat * slat;

    plg_[0][0] = 1.0;
    plg_[1][0] = clat;
    plg_[2][0] = 0.5 * (3.0 * clat2 - 1.0);
    plg_[3][0] = 0.5 * (5.0 * clat * clat2 - 3.0 * clat);
    plg_[4][0] = (35.0 * clat4 - 30.0 * clat2 + 3.0) / 8.0;
    plg_[5][0] = (63.0 * clat2 * clat2 * clat - 70.0 * clat2 * clat + 15.0 * clat) / 8.0;
    plg_[6][0] = (11.0 * clat * plg_[5][0] - 5.0 * plg_[4][0]) / 6.0;

    plg_[1][1] = slat;
    plg_[2][1] = 3.0 * clat * slat;
    plg_[3][1] = 1.5 * (5.0 * clat2 - 1.0) * slat;
    plg_[4][1] = 2.5 * (7.0 * clat2 * clat - 3.0 * clat) * slat;
    plg_[5][1] = 1.875 * (21.0 * clat4 - 14.0 * clat2 + 1.0) * slat;
    plg_[6][1] = (11.0 * clat * plg_[5][1] - 6.0 * plg_[4][1]) / 5.0;

    plg_[2][2] = 3.0 * slat2;
    plg_[3][2] = 15.0 * slat2 * clat;
    plg_[4][2] = 7.5 * (7.0 * clat2 - 1.0) * slat2;
    plg_[5][2] = 3.0 * clat * plg_[4][2] - 2.0 * plg_[3][2];
    plg_[6][2] = (11.0 * clat * plg_[5][2] - 7.0 * plg_[4][2]) / 4.0;

    plg_[3][3] = 15.0 * slat2 * slat;
    plg_[4][3] = 105.0 * slat2 * slat * clat;
    plg_[5][3] = (9.0 * clat * plg_[4][3] - 7.0 * plg_[3][3]) / 2.0;
    plg_[6][3] = (11.0 * clat * plg_[5][3] - 8.0 * plg_[4][3]) / 3.0;

    lastlat_ = in.lat;
  }

  if (in.doy != lastdoy_) {
    cdoy_[1] = std::cos(kDoy2Rad * in.doy);
    sdoy_[1] = std::sin(kDoy2Rad * in.doy);
    cdoy_[2] = std::cos(kDoy2Rad * in.doy * 2.0);
    sdoy_[2] = std::sin(kDoy2Rad * in.doy * 2.0);
    lastdoy_ = in.doy;
  }

  const double lst = std::fmod(in.utsec / 3600.0 + in.lon / 15.0 + 24.0, 24.0);
  if (lst != lastlst_) {
    clst_[1] = std::cos(kLst2Rad * lst);
    slst_[1] = std::sin(kLst2Rad * lst);
    clst_[2] = std::cos(kLst2Rad * lst * 2.0);
    slst_[2] = std::sin(kLst2Rad * lst * 2.0);
    clst_[3] = std::cos(kLst2Rad * lst * 3.0);
    slst_[3] = std::sin(kLst2Rad * lst * 3.0);
    lastlst_ = lst;
  }

  if (in.lon != lastlon_) {
    clon_[1] = std::cos(kDeg2Rad * in.lon);
    slon_[1] = std::sin(kDeg2Rad * in.lon);
    clon_[2] = std::cos(kDeg2Rad * in.lon * 2.0);
    slon_[2] = std::sin(kDeg2Rad * in.lon * 2.0);
    lastlon_ = in.lon;
  }

  int c = kCtimeInd;
  for (int n = 0; n <= kAmaxN; ++n) {
    bf[static_cast<std::size_t>(c)] = plg_[n][0];
    ++c;
  }

  if (c != kCintAnn) {
    throw std::runtime_error("problem with basis definitions");
  }
  for (int s = 1; s <= kAmaxS; ++s) {
    const double cosdoy = cdoy_[s];
    const double sindoy = sdoy_[s];
    for (int n = 0; n <= kAmaxN; ++n) {
      const double pl = plg_[n][0];
      bf[static_cast<std::size_t>(c)] = pl * cosdoy;
      bf[static_cast<std::size_t>(c + 1)] = pl * sindoy;
      c += 2;
    }
  }

  if (c != kCtide) {
    throw std::runtime_error("problem with basis definitions");
  }
  for (int l = 1; l <= kTmaxL; ++l) {
    const double coslst = clst_[l];
    const double sinlst = slst_[l];
    for (int n = l; n <= kTmaxN; ++n) {
      const double pl = plg_[n][l];
      bf[static_cast<std::size_t>(c)] = pl * coslst;
      bf[static_cast<std::size_t>(c + 1)] = pl * sinlst;
      c += 2;
    }
    for (int s = 1; s <= kTmaxS; ++s) {
      const double cosdoy = cdoy_[s];
      const double sindoy = sdoy_[s];
      for (int n = l; n <= kTmaxN; ++n) {
        const double pl = plg_[n][l];
        bf[static_cast<std::size_t>(c)] = pl * coslst * cosdoy;
        bf[static_cast<std::size_t>(c + 1)] = pl * sinlst * cosdoy;
        bf[static_cast<std::size_t>(c + 2)] = pl * coslst * sindoy;
        bf[static_cast<std::size_t>(c + 3)] = pl * sinlst * sindoy;
        c += 4;
      }
    }
  }

  if (c != kCspw) {
    throw std::runtime_error("problem with basis definitions");
  }
  for (int m = 1; m <= kPmaxM; ++m) {
    const double coslon = clon_[m];
    const double sinlon = slon_[m];
    for (int n = m; n <= kPmaxN; ++n) {
      const double pl = plg_[n][m];
      bf[static_cast<std::size_t>(c)] = pl * coslon;
      bf[static_cast<std::size_t>(c + 1)] = pl * sinlon;
      c += 2;
    }
    for (int s = 1; s <= kPmaxS; ++s) {
      const double cosdoy = cdoy_[s];
      const double sindoy = sdoy_[s];
      for (int n = m; n <= kPmaxN; ++n) {
        const double pl = plg_[n][m];
        bf[static_cast<std::size_t>(c)] = pl * coslon * cosdoy;
        bf[static_cast<std::size_t>(c + 1)] = pl * sinlon * cosdoy;
        bf[static_cast<std::size_t>(c + 2)] = pl * coslon * sindoy;
        bf[static_cast<std::size_t>(c + 3)] = pl * sinlon * sindoy;
        c += 4;
      }
    }
  }

  if (c != kCsfx) {
    throw std::runtime_error("problem with basis definitions");
  }
  const double dfa = in.sfluxavg - kSfluxAvgReference;
  const double df = in.sflux - in.sfluxavg;
  bf[static_cast<std::size_t>(c)] = dfa;
  bf[static_cast<std::size_t>(c + 1)] = dfa * dfa;
  bf[static_cast<std::size_t>(c + 2)] = df;
  bf[static_cast<std::size_t>(c + 3)] = df * df;
  bf[static_cast<std::size_t>(c + 4)] = df * dfa;
  c += kNsfx;

  if (c != kCextra) {
    throw std::runtime_error("problem with basis definitions");
  }
  const double sza = solzen(in.doy, lst, in.lat, in.lon);
  bf[static_cast<std::size_t>(c)] = -0.5 * std::tanh((sza - 98.0) / 6.0);
  bf[static_cast<std::size_t>(c + 1)] = -0.5 * std::tanh((sza - 101.5) / 20.0);
  bf[static_cast<std::size_t>(c + 2)] = dfa * bf[static_cast<std::size_t>(c)];
  bf[static_cast<std::size_t>(c + 3)] = dfa * bf[static_cast<std::size_t>(c + 1)];
  bf[static_cast<std::size_t>(c + 4)] = dfa * plg_[2][0];
  bf[static_cast<std::size_t>(c + 5)] = dfa * plg_[4][0];
  bf[static_cast<std::size_t>(c + 6)] = dfa * plg_[0][0] * cdoy_[1];
  bf[static_cast<std::size_t>(c + 7)] = dfa * plg_[0][0] * sdoy_[1];
  bf[static_cast<std::size_t>(c + 8)] = dfa * plg_[0][0] * cdoy_[2];
  bf[static_cast<std::size_t>(c + 9)] = dfa * plg_[0][0] * sdoy_[2];
  if (in.sfluxavg <= kSfluxAvgQuadCutoff) {
    bf[static_cast<std::size_t>(c + 10)] = dfa * dfa;
  } else {
    bf[static_cast<std::size_t>(c + 10)] =
        (kSfluxAvgQuadCutoff - kSfluxAvgReference) *
        (2.0 * dfa - (kSfluxAvgQuadCutoff - kSfluxAvgReference));
  }
  bf[static_cast<std::size_t>(c + 11)] = bf[static_cast<std::size_t>(c + 10)] * plg_[2][0];
  bf[static_cast<std::size_t>(c + 12)] = bf[static_cast<std::size_t>(c + 10)] * plg_[4][0];
  bf[static_cast<std::size_t>(c + 13)] = df * plg_[2][0];
  bf[static_cast<std::size_t>(c + 14)] = df * plg_[4][0];

  c = kCnonLin;
  if (c != kCsfxMod) {
    throw std::runtime_error("problem with basis set");
  }
  bf[static_cast<std::size_t>(c)] = dfa;
  bf[static_cast<std::size_t>(c + 1)] = dfa * dfa;
  bf[static_cast<std::size_t>(c + 2)] = df;
  bf[static_cast<std::size_t>(c + 3)] = df * df;
  bf[static_cast<std::size_t>(c + 4)] = df * dfa;
  c += kNsfxMod;

  if (c != kCmag) {
    throw std::runtime_error("problem with basis set");
  }
  for (int i = 0; i < 7; ++i) {
    bf[static_cast<std::size_t>(c + i)] = in.ap[static_cast<std::size_t>(i)] - 4.0;
  }
  bf[static_cast<std::size_t>(c + 8)] = kDoy2Rad * in.doy;
  bf[static_cast<std::size_t>(c + 9)] = kLst2Rad * lst;
  bf[static_cast<std::size_t>(c + 10)] = kDeg2Rad * in.lon;
  bf[static_cast<std::size_t>(c + 11)] = kLst2Rad * in.utsec / 3600.0;
  bf[static_cast<std::size_t>(c + 12)] = std::abs(in.lat);
  c += 13;
  for (int m = 0; m <= 1; ++m) {
    for (int n = 0; n <= kAmaxN; ++n) {
      bf[static_cast<std::size_t>(c)] = plg_[n][m];
      ++c;
    }
  }

  c = kCut;
  bf[static_cast<std::size_t>(c)] = kLst2Rad * in.utsec / 3600.0;
  bf[static_cast<std::size_t>(c + 1)] = kDoy2Rad * in.doy;
  bf[static_cast<std::size_t>(c + 2)] = dfa;
  bf[static_cast<std::size_t>(c + 3)] = kDeg2Rad * in.lon;
  bf[static_cast<std::size_t>(c + 4)] = plg_[1][0];
  bf[static_cast<std::size_t>(c + 5)] = plg_[3][0];
  bf[static_cast<std::size_t>(c + 6)] = plg_[5][0];
  bf[static_cast<std::size_t>(c + 7)] = plg_[3][2];
  bf[static_cast<std::size_t>(c + 8)] = plg_[5][2];

  for (std::size_t i = 0; i <= static_cast<std::size_t>(kMbf); ++i) {
    if (!switches[i]) {
      bf[i] = 0.0;
    }
  }

  return bf;
}

double GlobeCalculator::solzen(double ddd, double lst, double lat, double lon) {
  static constexpr double humr = kPi / 12.0;
  static constexpr std::array<double, 5> p{0.017203534, 0.034407068, 0.051610602, 0.068814136,
                                           0.103221204};

  const double wlon = 360.0 - lon;
  const double teqnx = ddd + (lst + wlon / 15.0) / 24.0 + 0.9369;

  double dec = 23.256 * std::sin(p[0] * (teqnx - 82.242)) +
               0.381 * std::sin(p[1] * (teqnx - 44.855)) +
               0.167 * std::sin(p[2] * (teqnx - 23.355)) -
               0.013 * std::sin(p[3] * (teqnx + 11.97)) +
               0.011 * std::sin(p[4] * (teqnx - 10.410)) + 0.339137;
  dec *= kDeg2Rad;

  const double tf = teqnx - 0.5;
  const double teqt = -7.38 * std::sin(p[0] * (tf - 4.0)) - 9.87 * std::sin(p[1] * (tf + 9.0)) +
                      0.27 * std::sin(p[2] * (tf - 53.0)) -
                      0.2 * std::cos(p[3] * (tf - 17.0));

  const double phi = humr * (lst - 12.0) + teqt * kDeg2Rad / 4.0;
  const double rlat = lat * kDeg2Rad;

  double cosx = std::sin(rlat) * std::sin(dec) + std::cos(rlat) * std::cos(dec) * std::cos(phi);
  if (std::abs(cosx) > 1.0) {
    cosx = std::copysign(1.0, cosx);
  }

  return std::acos(cosx) / kDeg2Rad;
}

}  // namespace msis21::detail
