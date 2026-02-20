/**
 * @file test_utils.cpp
 * @brief Unit tests for msis utility functions.
 * @author Watosn
 */

#include "msis21/detail/utils.hpp"

#include <array>
#include <cmath>

#include <gtest/gtest.h>

TEST(Utils, AltitudeGeopotentialRoundTrip) {
  constexpr double lat = 45.0;
  constexpr double alt = 400.0;

  const double gph = msis21::detail::alt2gph(lat, alt);
  const double alt_back = msis21::detail::gph2alt(lat, gph);

  EXPECT_NEAR(alt_back, alt, 0.01);
}

TEST(Utils, DilogSanity) {
  const double d0 = msis21::detail::dilog(0.0);
  const double d1 = msis21::detail::dilog(0.5);
  EXPECT_NEAR(d0, 0.0, 1e-12);
  EXPECT_TRUE(std::isfinite(d1));
  EXPECT_GT(d1, 0.5);
}

TEST(Utils, BsplineProducesNonNegativeBasis) {
  constexpr int nd = 8;
  std::array<double, nd + 1> nodes = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  msis21::detail::EtaTable eta{};

  for (int k = 2; k <= 6; ++k) {
    for (int j = 0; j <= nd - k + 1; ++j) {
      eta[j][k] = 1.0 / (nodes[static_cast<std::size_t>(j + k - 1)] - nodes[static_cast<std::size_t>(j)]);
    }
  }

  const auto res = msis21::detail::bspline(3.4, nodes, nd, 6, eta);
  EXPECT_GE(res.i, 0);

  double sum = 0.0;
  for (int rel = -5; rel <= 0; ++rel) {
    const double v = res.value(rel, 6);
    EXPECT_GE(v, 0.0);
    sum += v;
  }
  EXPECT_GT(sum, 0.8);
  EXPECT_LT(sum, 1.2);
}
