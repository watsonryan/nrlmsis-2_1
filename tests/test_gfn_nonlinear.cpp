/**
 * @file test_gfn_nonlinear.cpp
 * @brief Unit tests for legacy nonlinear basis terms.
 * @author Watosn
 */

#include "msis21/detail/gfn.hpp"

#include <array>
#include <cmath>

#include <gtest/gtest.h>

namespace {

double g0fn_ref(double a, double k00r, double k00s) {
  return a + (k00r - 1.0) * (a + (std::exp(-a * k00s) - 1.0) / k00s);
}

}  // namespace

TEST(GfnNonlinear, SfluxModUsesExpectedBuckets) {
  std::array<double, msis21::detail::kMaxBasisFunctions> gf{};
  std::array<double, msis21::detail::kMaxBasisFunctions> beta_col{};
  std::array<bool, msis21::detail::kMaxBasisFunctions> swg{};
  gf.fill(1.0);
  beta_col.fill(1.0);
  swg.fill(true);

  const double val = msis21::detail::sfluxmod(0, gf, beta_col, 2.0, swg);
  EXPECT_DOUBLE_EQ(val, 840.0);
}

TEST(GfnNonlinear, GeomagDailyModeSimpleCase) {
  std::array<double, msis21::detail::kNmag> p0{};
  std::array<double, 13> bf{};
  std::array<double, 14> plg{};
  std::array<bool, msis21::detail::kNmag> swg_mag{};
  swg_mag.fill(false);
  swg_mag[0] = true;
  swg_mag[1] = true;
  swg_mag[2] = true;
  p0[0] = 0.5;
  p0[1] = 0.1;
  p0[2] = 2.0;
  bf[0] = 3.0;
  plg[0] = 1.0;

  const double expected = 2.0 * g0fn_ref(3.0, 0.5, 0.1);
  const double got = msis21::detail::geomag(p0, bf, plg, swg_mag);
  EXPECT_NEAR(got, expected, 1e-12);
}

TEST(GfnNonlinear, UtdepSwitchesDisableTerms) {
  std::array<double, msis21::detail::kNut> p0{};
  std::array<double, 9> bf{};
  std::array<bool, msis21::detail::kNut> swg_ut{};
  p0.fill(1.0);
  bf.fill(1.0);
  swg_ut.fill(false);
  const double got = msis21::detail::utdep(p0, bf, swg_ut);
  EXPECT_DOUBLE_EQ(got, 0.0);
}
