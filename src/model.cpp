/**
 * @file model.cpp
 * @brief Public model API implementation.
 * @author Watosn
 */

#include "msis21/model.hpp"

namespace msis21 {

Model Model::load_from_file(const std::filesystem::path& /*parm_path*/, Options /*options*/) {
  return Model(Status::NotInitialized);
}

Result Model::evaluate(const Input& input) const noexcept {
  Scratch scratch;
  return evaluate(input, scratch);
}

Result Model::evaluate(const Input& /*input*/, Scratch& /*scratch*/) const noexcept {
  return Result{.status = init_status_};
}

}  // namespace msis21
