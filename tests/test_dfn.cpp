/**
 * @file test_dfn.cpp
 * @brief Unit tests for density profile helper functions.
 * @author Watosn
 */

#include "msis21/detail/dfn.hpp"

#include <array>

#include <gtest/gtest.h>

TEST(Dfn, PiecewiseMassInterpolation) {
  const std::array<double, 5> z = {70.0, 80.0, 90.0, 100.0, 110.0};
  const std::array<double, 5> m = {28.0, 25.0, 22.0, 20.0, 18.0};
  const std::array<double, 5> dmdz = {-0.3, -0.3, -0.2, -0.2, 0.0};

  EXPECT_DOUBLE_EQ(msis21::detail::pwmp(60.0, z, m, dmdz), 28.0);
  EXPECT_DOUBLE_EQ(msis21::detail::pwmp(120.0, z, m, dmdz), 18.0);
  EXPECT_NEAR(msis21::detail::pwmp(85.0, z, m, dmdz), 23.5, 1e-12);
}

TEST(Dfn, AnomalousOxygenPath) {
  msis21::detail::TnParm tpro;
  msis21::detail::DnParm dpro;
  dpro.ispec = 9;
  dpro.lndref = 30.0;
  dpro.zref = 120.0;
  dpro.c = 0.1;
  dpro.zeta_c = 100.0;
  dpro.hc = 20.0;
  dpro.zmin = -1.0;

  const double dens = msis21::detail::dfnx(140.0, 700.0, 40.0, 0.0, 0.0, 1.0, tpro, dpro);
  EXPECT_GT(dens, 0.0);
}
