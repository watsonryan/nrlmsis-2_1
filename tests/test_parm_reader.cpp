/**
 * @file test_parm_reader.cpp
 * @brief Tests for parameter file loading.
 * @author Watosn
 */

#include "msis21/detail/constants.hpp"
#include "msis21/detail/parm_reader.hpp"

#include <filesystem>

#include <gtest/gtest.h>

TEST(ParmReader, LoadsReferenceParmFile) {
  msis21::detail::Parameters params;
  const auto status =
      msis21::detail::load_parameters(std::filesystem::path{"data/msis21.parm"}, params);
  ASSERT_EQ(status, msis21::Status::Ok);
  EXPECT_EQ(params.rows, msis21::detail::kMaxBasisFunctions);
  EXPECT_EQ(params.cols, msis21::detail::kParmColumnCount);
  EXPECT_EQ(params.beta.size(), params.rows * params.cols);
}
