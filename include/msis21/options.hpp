/**
 * @file options.hpp
 * @brief Runtime options for model initialization.
 * @author Watosn
 */
#pragma once

#include <array>

namespace msis21 {

struct Switches {
  std::array<int, 25> legacy{};
};

struct Options {
  Switches switches{};
};

}  // namespace msis21
