/**
 * @file test_tfn.cpp
 * @brief Unit tests for tfn helper functions.
 * @author Watosn
 */

#include "msis21/detail/tfn.hpp"

#include <array>
#include <cmath>

#include <gtest/gtest.h>

TEST(Tfn, ContinuityCoefficientsFinite) {
  const auto bc = msis21::detail::continuity_coefficients(1000.0, 5.0, 200.0);
  EXPECT_TRUE(std::isfinite(bc[0]));
  EXPECT_TRUE(std::isfinite(bc[1]));
  EXPECT_TRUE(std::isfinite(bc[2]));
}

TEST(Tfn, BatesRegionTemperature) {
  msis21::detail::TnParm tpro;
  tpro.tex = 900.0;
  tpro.tb0 = 300.0;
  tpro.sigma = 0.05;

  const double t = msis21::detail::tfnx(130.0, 0, std::array<double, 4>{0.0, 0.0, 0.0, 0.0}, tpro);
  EXPECT_GT(t, tpro.tb0);
  EXPECT_LT(t, tpro.tex);
}
