/**
 * @file test_gtd8d.cpp
 * @brief Unit tests for legacy wrappers.
 * @author Watosn
 */

#include "msis21/detail/gtd8d.hpp"

#include <gtest/gtest.h>

TEST(Gtd8d, RejectsWhenNotInitialized) {
  msis21::LegacyState state;
  msis21::Input in{};
  const auto out = msis21::gtd8d(state, in);
  EXPECT_EQ(out.status, msis21::Status::NotInitialized);
}
