/**
 * @file test_invariants.cpp
 * @brief Physical and reproducibility invariants for msis21 outputs.
 * @author Watosn
 */

#include "msis21/msis21.hpp"

#include <array>
#include <cmath>

#include <gtest/gtest.h>

#include "test_paths.hpp"

TEST(Invariants, DensitiesNonNegativeAndFinite) {
  auto model = msis21::Model::load_from_file(msis21_data_path("msis21.parm"), msis21::Options{});

  msis21::Input in{};
  in.iyd = 70178;
  in.sec = 43200.0;
  in.glat_deg = 25.0;
  in.glon_deg = -45.0;
  in.stl_hr = 12.0;
  in.f107a = 150.0;
  in.f107 = 150.0;
  in.ap = 8.0;

  for (double z = 250.0; z <= 2000.0; z += 125.0) {
    in.alt_km = z;
    const auto out = model.evaluate(in);
    ASSERT_EQ(out.status, msis21::Status::Ok);
    EXPECT_TRUE(std::isfinite(out.out.t));
    EXPECT_GT(out.out.t, 100.0);
    EXPECT_GE(out.out.rho, 0.0);
    EXPECT_GE(out.out.he, 0.0);
    EXPECT_GE(out.out.o, 0.0);
    EXPECT_GE(out.out.n2, 0.0);
    EXPECT_GE(out.out.o2, 0.0);
    EXPECT_GE(out.out.ar, 0.0);
    EXPECT_GE(out.out.h, 0.0);
    EXPECT_GE(out.out.n, 0.0);
    EXPECT_GE(out.out.o_anom, 0.0);
    EXPECT_GE(out.out.no, 0.0);
  }
}

TEST(Invariants, EvaluateDeterministicRepeated) {
  auto model = msis21::Model::load_from_file(msis21_data_path("msis21.parm"), msis21::Options{});

  msis21::Input in{};
  in.iyd = 70200;
  in.sec = 30000.0;
  in.alt_km = 400.0;
  in.glat_deg = 10.0;
  in.glon_deg = 120.0;
  in.stl_hr = 16.0;
  in.f107a = 130.0;
  in.f107 = 120.0;
  in.ap = 6.0;

  const auto a = model.evaluate(in);
  const auto b = model.evaluate(in);
  ASSERT_EQ(a.status, msis21::Status::Ok);
  ASSERT_EQ(b.status, msis21::Status::Ok);
  EXPECT_DOUBLE_EQ(a.out.rho, b.out.rho);
  EXPECT_DOUBLE_EQ(a.out.o, b.out.o);
  EXPECT_DOUBLE_EQ(a.out.t, b.out.t);
}

TEST(Invariants, BatchEvaluationMatchesScalar) {
  auto model = msis21::Model::load_from_file(msis21_data_path("msis21.parm"), msis21::Options{});

  std::array<msis21::Input, 3> in{};
  in[0] = msis21::Input{70178, 40000.0, 300.0, 0.0, 0.0, 12.0, 150.0, 150.0, 4.0};
  in[1] = msis21::Input{70178, 45000.0, 500.0, 20.0, 40.0, 13.0, 160.0, 155.0, 12.0};
  in[2] = msis21::Input{70178, 50000.0, 900.0, -10.0, -70.0, 14.0, 120.0, 110.0, 20.0};

  std::array<msis21::Output, 3> out_batch{};
  const auto status = model.evaluate_batch(in, out_batch);
  ASSERT_EQ(status, msis21::Status::Ok);

  for (std::size_t i = 0; i < in.size(); ++i) {
    const auto scalar = model.evaluate(in[i]);
    ASSERT_EQ(scalar.status, msis21::Status::Ok);
    EXPECT_DOUBLE_EQ(out_batch[i].rho, scalar.out.rho);
    EXPECT_DOUBLE_EQ(out_batch[i].o, scalar.out.o);
    EXPECT_DOUBLE_EQ(out_batch[i].t, scalar.out.t);
  }
}
