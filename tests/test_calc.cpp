/**
 * @file test_calc.cpp
 * @brief Unit tests for msiscalc plumbing.
 * @author Watosn
 */

#include "msis21/detail/calc.hpp"
#include "msis21/detail/constants.hpp"
#include "msis21/detail/parm_reader.hpp"

#include <gtest/gtest.h>

TEST(Calc, InvalidInputRejected) {
  msis21::detail::Parameters params;
  params.rows = msis21::detail::kMaxBasisFunctions;
  params.cols = msis21::detail::kParmColumnCount;
  params.beta.resize(params.rows * params.cols, 0.0);

  msis21::Input in;
  in.iyd = 70178;
  in.sec = -1.0;

  msis21::Options options;
  const auto res = msis21::detail::evaluate_msiscalc(in, options, params);
  EXPECT_EQ(res.status, msis21::Status::InvalidInput);
}
