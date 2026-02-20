/**
 * @file gtd8d.hpp
 * @brief Legacy wrapper compatibility helpers.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <optional>

#include "msis21/model.hpp"

namespace msis21 {

struct LegacyState {
  std::optional<Model> model;
};

Status msisinit(LegacyState& state, const std::filesystem::path& parm_path, Options options = {});
Result msiscalc(const LegacyState& state, const Input& input);
Result gtd8d(const LegacyState& state, const Input& input);

}  // namespace msis21
