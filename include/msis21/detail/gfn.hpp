/**
 * @file gfn.hpp
 * @brief Horizontal/time-dependent basis functions from msis_gfn.F90.
 * @author Watosn
 */
#pragma once

#include <array>
#include <cstddef>

#include "msis21/detail/constants.hpp"

namespace msis21::detail {

struct GlobeInput {
  double doy{};
  double utsec{};
  double lat{};
  double lon{};
  double sfluxavg{};
  double sflux{};
  std::array<double, 7> ap{};
};

class GlobeCalculator {
 public:
  GlobeCalculator();

  [[nodiscard]] std::array<double, kMaxBasisFunctions> globe(
      const GlobeInput& in,
      const std::array<bool, kMaxBasisFunctions>& switches);

  [[nodiscard]] static double solzen(double ddd, double lst, double lat, double lon);

 private:
  std::array<std::array<double, kMaxLegendreDegree + 1>, kMaxLegendreDegree + 1> plg_{};
  std::array<double, 3> cdoy_{};
  std::array<double, 3> sdoy_{};
  std::array<double, 4> clst_{};
  std::array<double, 4> slst_{};
  std::array<double, 3> clon_{};
  std::array<double, 3> slon_{};

  double lastlat_{-999.9};
  double lastdoy_{-999.9};
  double lastlst_{-999.9};
  double lastlon_{-999.9};
};

[[nodiscard]] double sfluxmod(
    int iz,
    const std::array<double, kMaxBasisFunctions>& gf,
    const std::array<double, kMaxBasisFunctions>& beta_col,
    double dffact,
    const std::array<bool, kMaxBasisFunctions>& swg);

[[nodiscard]] double geomag(
    const std::array<double, kNmag>& p0,
    const std::array<double, 13>& bf,
    const std::array<double, 14>& plg,
    const std::array<bool, kNmag>& swg_mag);

[[nodiscard]] double utdep(
    const std::array<double, kNut>& p0,
    const std::array<double, 9>& bf,
    const std::array<bool, kNut>& swg_ut);

}  // namespace msis21::detail
