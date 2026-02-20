/**
 * @file utils.hpp
 * @brief Utility math functions ported from msis_utils.F90.
 * @author Watosn
 */
#pragma once

#include <array>
#include <span>

namespace msis21::detail {

using EtaTable = std::array<std::array<double, 7>, 31>;
using BsplineTable = std::array<std::array<double, 7>, 6>;

struct BsplineResult {
  BsplineTable s{};
  int i{-1};

  [[nodiscard]] double value(int rel_index, int order) const;
  void set(int rel_index, int order, double v);
};

[[nodiscard]] double alt2gph(double lat_deg, double alt_km);
[[nodiscard]] double gph2alt(double lat_deg, double gph_km);
[[nodiscard]] BsplineResult bspline(double x,
                                    std::span<const double> nodes,
                                    int nd,
                                    int kmax,
                                    const EtaTable& eta);
[[nodiscard]] double dilog(double x);

}  // namespace msis21::detail
