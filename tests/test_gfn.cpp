/**
 * @file test_gfn.cpp
 * @brief Unit tests for horizontal basis function generation.
 * @author Watosn
 */

#include "msis21/detail/gfn.hpp"

#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>

TEST(Gfn, GeneratesFiniteBasis) {
  msis21::detail::GlobeCalculator calc;
  msis21::detail::GlobeInput in;
  in.doy = 172.0;
  in.utsec = 36000.0;
  in.lat = 45.0;
  in.lon = -120.0;
  in.sfluxavg = 150.0;
  in.sflux = 140.0;
  in.ap = {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0};

  std::array<bool, msis21::detail::kMaxBasisFunctions> switches;
  switches.fill(true);

  const auto bf = calc.globe(in, switches);
  EXPECT_TRUE(std::all_of(bf.begin(), bf.end(), [](double x) { return std::isfinite(x); }));
  EXPECT_NEAR(bf[0], 1.0, 1e-12);
}
