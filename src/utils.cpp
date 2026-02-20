/**
 * @file utils.cpp
 * @brief Utility math function implementations ported from msis_utils.F90.
 * @author Watosn
 */

#include "msis21/detail/utils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace msis21::detail {

double BsplineResult::value(int rel_index, int order) const {
  return s[rel_index + 5][order];
}

void BsplineResult::set(int rel_index, int order, double v) { s[rel_index + 5][order] = v; }

double alt2gph(double lat_deg, double alt_km) {
  constexpr double deg2rad = 0.017453292519943295;
  constexpr double a = 6378.1370 * 1e3;
  constexpr double finv = 298.257223563;
  constexpr double w = 7292115e-11;
  const double gm = static_cast<double>(398600.4418f) * 1.0e9;

  constexpr double asq = a * a;
  constexpr double wsq = w * w;
  constexpr double f = 1.0 / finv;
  constexpr double esq = 2.0 * f - f * f;
  const double e = std::sqrt(esq);
  const double elin = a * e;
  const double elinsq = elin * elin;
  const double epr = e / (1 - f);
  const double q0 = ((1.0 + 3.0 / (epr * epr)) * std::atan(epr) - 3.0 / epr) / 2.0;
  const double u0 = -gm * std::atan(epr) / elin - wsq * asq / 3.0;
  constexpr double g0 = 9.80665;
  const double gm_div_elin = gm / elin;

  constexpr double x0sq = 2e7 * 2e7;
  constexpr double hsq = 1.2e7 * 1.2e7;

  const double altm = alt_km * 1000.0;
  const double sinlat = std::sin(lat_deg * deg2rad);
  const double sinsq = sinlat * sinlat;
  const double v = a / std::sqrt(1 - esq * sinsq);
  const double xp = v + altm;
  const double zp = v * (1 - esq) + altm;
  const double xsq = xp * xp * (1 - sinsq);
  const double zsq = zp * zp * sinsq;
  const double rsq_min_elin = xsq + zsq - elinsq;
  const double usq = rsq_min_elin / 2.0 + std::sqrt((rsq_min_elin * rsq_min_elin) / 4.0 + elinsq * zsq);
  const double cossqdelta = zsq / usq;

  const double epru = elin / std::sqrt(usq);
  const double atan_epru = std::atan(epru);
  const double q = ((1 + 3.0 / (epru * epru)) * atan_epru - 3.0 / epru) / 2.0;
  double u = -gm_div_elin * atan_epru - wsq * (asq * q * (cossqdelta - 1 / 3.0) / q0) / 2.0;

  const double vc = (xsq <= x0sq) ? ((wsq / 2.0) * xsq)
                                  : ((wsq / 2.0) * (hsq * std::tanh((xsq - x0sq) / hsq) + x0sq));
  u -= vc;

  return (u - u0) / g0 / 1000.0;
}

double gph2alt(double lat_deg, double gph_km) {
  constexpr int max_iter = 10;
  constexpr double epsilon = 0.0005;

  double x = gph_km;
  int n = 0;
  double dx = epsilon + epsilon;
  while (std::abs(dx) > epsilon && n < max_iter) {
    const double y = alt2gph(lat_deg, x);
    const double dydz = (alt2gph(lat_deg, x + dx) - y) / dx;
    dx = (gph_km - y) / dydz;
    x += dx;
    ++n;
  }
  return x;
}

BsplineResult bspline(double x, std::span<const double> nodes, int nd, int kmax, const EtaTable& eta) {
  if (kmax < 2 || kmax > 6) {
    throw std::invalid_argument("kmax must be in [2,6]");
  }
  if (static_cast<int>(nodes.size()) < nd + 1) {
    throw std::invalid_argument("nodes too short");
  }

  BsplineResult result;

  if (x >= nodes[nd]) {
    result.i = nd;
    return result;
  }
  if (x <= nodes[0]) {
    result.i = -1;
    return result;
  }

  int low = 0;
  int high = nd;
  int i = (low + high) / 2;
  while (x < nodes[i] || x >= nodes[i + 1]) {
    if (x < nodes[i]) {
      high = i;
    } else {
      low = i;
    }
    i = (low + high) / 2;
  }
  result.i = i;

  std::array<double, 5> w{};

  result.set(0, 2, (x - nodes[i]) * eta[i][2]);
  if (i > 0) {
    result.set(-1, 2, 1 - result.value(0, 2));
  }
  if (i >= nd - 1) {
    result.set(0, 2, 0.0);
  }

  w[4] = (x - nodes[i]) * eta[i][3];
  if (i != 0) {
    w[3] = (x - nodes[i - 1]) * eta[i - 1][3];
  }
  if (i < nd - 2) {
    result.set(0, 3, w[4] * result.value(0, 2));
  }
  if ((i - 1) >= 0 && (i - 1) < (nd - 2)) {
    result.set(-1, 3, w[3] * result.value(-1, 2) + (1.0 - w[4]) * result.value(0, 2));
  }
  if ((i - 2) >= 0) {
    result.set(-2, 3, (1.0 - w[3]) * result.value(-1, 2));
  }

  for (int l = 0; l >= -2; --l) {
    const int j = i + l;
    if (j < 0) {
      break;
    }
    w[l + 4] = (x - nodes[j]) * eta[j][4];
  }
  if (i < nd - 3) {
    result.set(0, 4, w[4] * result.value(0, 3));
  }
  for (int l = -1; l >= -2; --l) {
    if ((i + l) >= 0 && (i + l) < (nd - 3)) {
      result.set(l, 4, w[l + 4] * result.value(l, 3) + (1.0 - w[l + 5]) * result.value(l + 1, 3));
    }
  }
  if ((i - 3) >= 0) {
    result.set(-3, 4, (1.0 - w[2]) * result.value(-2, 3));
  }

  for (int l = 0; l >= -3; --l) {
    const int j = i + l;
    if (j < 0) {
      break;
    }
    w[l + 4] = (x - nodes[j]) * eta[j][5];
  }
  if (i < nd - 4) {
    result.set(0, 5, w[4] * result.value(0, 4));
  }
  for (int l = -1; l >= -3; --l) {
    if ((i + l) >= 0 && (i + l) < (nd - 4)) {
      result.set(l, 5, w[l + 4] * result.value(l, 4) + (1.0 - w[l + 5]) * result.value(l + 1, 4));
    }
  }
  if ((i - 4) >= 0) {
    result.set(-4, 5, (1.0 - w[1]) * result.value(-3, 4));
  }
  if (kmax == 5) {
    return result;
  }

  for (int l = 0; l >= -4; --l) {
    const int j = i + l;
    if (j < 0) {
      break;
    }
    w[l + 4] = (x - nodes[j]) * eta[j][6];
  }
  if (i < nd - 5) {
    result.set(0, 6, w[4] * result.value(0, 5));
  }
  for (int l = -1; l >= -4; --l) {
    if ((i + l) >= 0 && (i + l) < (nd - 5)) {
      result.set(l, 6, w[l + 4] * result.value(l, 5) + (1.0 - w[l + 5]) * result.value(l + 1, 5));
    }
  }
  if ((i - 5) >= 0) {
    result.set(-5, 6, (1.0 - w[0]) * result.value(-4, 5));
  }

  return result;
}

double dilog(double x) {
  constexpr double pi = 3.1415926535897932384626433832795;
  constexpr double pi2_6 = pi * pi / 6.0;

  if (x > 0.5) {
    const double lnx = std::log(x);
    x = 1.0 - x;
    const double xx = x * x;
    const double x4 = 4.0 * x;
    return pi2_6 - lnx * std::log(x) -
           (4.0 * xx *
                (23.0 / 16.0 + x / 36.0 + xx / 576.0 + xx * x / 3600.0) +
            x4 + 3.0 * (1.0 - xx) * lnx) /
               (1.0 + x4 + xx);
  }

  const double xx = x * x;
  const double x4 = 4.0 * x;
  return (4.0 * xx *
              (23.0 / 16.0 + x / 36.0 + xx / 576.0 + xx * x / 3600.0) +
          x4 + 3.0 * (1.0 - xx) * std::log(1.0 - x)) /
         (1.0 + x4 + xx);
}

}  // namespace msis21::detail
