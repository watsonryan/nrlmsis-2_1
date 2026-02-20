/**
 * @file model.cpp
 * @brief Public model API implementation.
 * @author Watosn
 */

#include "msis21/model.hpp"

#include "msis21/detail/parm_reader.hpp"

namespace msis21 {

Model Model::load_from_file(const std::filesystem::path& parm_path, Options /*options*/) {
  auto parameters = std::make_shared<detail::Parameters>();
  const Status status = detail::load_parameters(parm_path, *parameters);
  if (status != Status::Ok) {
    return Model(status, nullptr);
  }
  return Model(Status::Ok, parameters);
}

Result Model::evaluate(const Input& input) const noexcept {
  Scratch scratch;
  return evaluate(input, scratch);
}

Result Model::evaluate(const Input& /*input*/, Scratch& /*scratch*/) const noexcept {
  if (init_status_ != Status::Ok || !parameters_) {
    return Result{.status = init_status_};
  }
  return Result{.status = Status::Ok};
}

}  // namespace msis21
