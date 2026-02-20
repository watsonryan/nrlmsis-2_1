/**
 * @file model.cpp
 * @brief Public model API implementation.
 * @author Watosn
 */

#include "msis21/model.hpp"

#include <utility>

#include "msis21/detail/calc.hpp"
#include "msis21/detail/fortran_backend.hpp"
#include "msis21/detail/log.hpp"
#include "msis21/detail/parm_reader.hpp"
#include <spdlog/spdlog.h>

namespace msis21 {

Model Model::load_from_file(const std::filesystem::path& parm_path, Options options) {
  detail::configure_logging_once();
  spdlog::info("Initializing msis21 model from {}", parm_path.string());
  auto parameters = std::make_shared<detail::Parameters>();
  const Status status = detail::load_parameters(parm_path, *parameters);
  if (status != Status::Ok) {
    spdlog::error("Parameter load failed with status {}", static_cast<int>(status));
    return Model(status, nullptr, std::move(options));
  }
  const Status fortran_status = detail::initialize_fortran_backend(parm_path);
  if (fortran_status != Status::Ok) {
    spdlog::error("Fortran backend init failed with status {}", static_cast<int>(fortran_status));
    return Model(fortran_status, nullptr, std::move(options));
  }
  spdlog::info("msis21 model initialized");
  return Model(Status::Ok, parameters, std::move(options));
}

Result Model::evaluate(const Input& input) const noexcept {
  Scratch scratch;
  return evaluate(input, scratch);
}

Result Model::evaluate(const Input& input, Scratch& /*scratch*/) const noexcept {
  if (init_status_ != Status::Ok || !parameters_) {
    return Result{.status = init_status_};
  }
  const auto calc = detail::evaluate_msiscalc(input, options_, *parameters_);
  return Result{.status = calc.status, .out = calc.out};
}

}  // namespace msis21
