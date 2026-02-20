/**
 * @file gtd8d.cpp
 * @brief Legacy wrapper compatibility implementation.
 * @author Watosn
 */

#include "msis21/detail/gtd8d.hpp"

#include <utility>

namespace msis21 {

Status msisinit(LegacyState& state, const std::filesystem::path& parm_path, Options options) {
  state.model = Model::load_from_file(parm_path, std::move(options));
  auto probe = state.model->evaluate(Input{});
  if (probe.status == Status::InvalidInput) {
    return Status::Ok;
  }
  return probe.status == Status::NumericalError ? Status::Ok : probe.status;
}

Result msiscalc(const LegacyState& state, const Input& input) {
  if (!state.model.has_value()) {
    return Result{.status = Status::NotInitialized};
  }
  return state.model->evaluate(input);
}

Result gtd8d(const LegacyState& state, const Input& input) { return msiscalc(state, input); }

}  // namespace msis21
