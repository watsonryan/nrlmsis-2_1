/**
 * @file dfn.hpp
 * @brief Species density profile helpers from msis_dfn.F90.
 * @author Watosn
 */
#pragma once

#include <array>

#include "msis21/detail/tfn.hpp"

namespace msis21::detail {

struct DnParm {
  double lnphi_f{};
  double lndref{};
  double zeta_m{};
  double hml{};
  double hmu{};
  double c{};
  double zeta_c{};
  double hc{};
  double r{};
  double zeta_r{};
  double hr{};
  std::array<double, kNsplO1 + 2> cf{};
  double zref{};
  std::array<double, 5> mi{};
  std::array<double, 5> zeta_mi{};
  std::array<double, 5> ami{};
  std::array<double, 5> wmi{};
  std::array<double, 5> xmi{};
  double izref{};
  double tref{};
  double zmin{};
  double zhyd{};
  int ispec{};
};

[[nodiscard]] double pwmp(double z,
                          const std::array<double, 5>& zm,
                          const std::array<double, 5>& m,
                          const std::array<double, 5>& dmdz);

[[nodiscard]] double dfnx(double z,
                          double tnz,
                          double lndtotz,
                          double vz,
                          double wz,
                          double hrfact,
                          const TnParm& tpro,
                          const DnParm& dpro);

}  // namespace msis21::detail
