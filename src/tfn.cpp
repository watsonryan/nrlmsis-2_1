/**
 * @file tfn.cpp
 * @brief Vertical temperature profile utilities from msis_tfn.F90.
 * @author Watosn
 */

#include "msis21/detail/tfn.hpp"

#include <algorithm>
#include <cmath>

namespace msis21::detail {

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

}  // namespace msis21::detail
