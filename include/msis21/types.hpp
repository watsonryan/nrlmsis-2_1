/**
 * @file types.hpp
 * @brief Core request/response structures for MSIS evaluation.
 * @author Watosn
 */
#pragma once

namespace msis21 {

struct Input {
  int iyd{};
  double sec{};
  double alt_km{};
  double glat_deg{};
  double glon_deg{};
  double stl_hr{};
  double f107a{};
  double f107{};
  double ap{};
};

struct Output {
  double he{};
  double o{};
  double n2{};
  double o2{};
  double ar{};
  double rho{};
  double h{};
  double n{};
  double o_anom{};
  double no{};
  double t{};
};

}  // namespace msis21
