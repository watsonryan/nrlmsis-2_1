/**
 * @file tfn.hpp
 * @brief Vertical temperature profile utilities from msis_tfn.F90.
 * @author Watosn
 */
#pragma once

#include <array>

#include "msis21/detail/constants.hpp"
#include "msis21/detail/parm_reader.hpp"

namespace msis21::detail {

struct TnParm {
  std::array<double, kNl + 1> cf{};
  double tzeta_f{};
  double tzeta_a{};
  double dlntdz_a{};
  double lndtot_f{};
  double tex{};
  double tgb0{};
  double tb0{};
  double sigma{};
  double sigmasq{};
  double b{};
  std::array<double, kNl + 1> beta{};
  std::array<double, kNl + 1> gamma{};
  double cvs{};
  double cvb{};
  double cws{};
  double cwb{};
  double vzeta_f{};
  double vzeta_a{};
  double wzeta_a{};
  double vzeta0{};
};

[[nodiscard]] std::array<double, 3> continuity_coefficients(double tex, double tgb0, double tb0);
[[nodiscard]] double tfnx(double z, int iz, const std::array<double, 4>& wght, const TnParm& tpro);
[[nodiscard]] TnParm tfnparm(const std::array<double, kMaxBasisFunctions>& gf,
                             const std::array<bool, kMaxBasisFunctions>& swg,
                             const Parameters& parameters);

}  // namespace msis21::detail
