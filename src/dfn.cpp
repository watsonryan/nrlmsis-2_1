/**
 * @file dfn.cpp
 * @brief Density profile helper implementations from msis_dfn.F90.
 * @author Watosn
 */

#include "msis21/detail/dfn.hpp"

#include <cmath>

namespace msis21::detail {

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

}  // namespace msis21::detail
