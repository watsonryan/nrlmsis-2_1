/**
 * @file test_smoke.cpp
 * @brief Basic smoke tests for project wiring.
 * @author Watosn
 */

#include <gtest/gtest.h>

TEST(Smoke, BasicMath) {
  EXPECT_EQ(2 + 2, 4);
}
