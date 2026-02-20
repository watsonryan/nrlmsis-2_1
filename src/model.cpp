/**
 * @file model.cpp
 * @brief Public model API implementation.
 * @author Watosn
 */

#include "msis21/model.hpp"

#include <utility>

#include "msis21/detail/calc.hpp"
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
    spdlog::error("Parameter load failed with status {} ({})", status_to_string(status), static_cast<int>(status));
    return Model(status, nullptr, std::move(options));
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

Status Model::evaluate_batch(std::span<const Input> in, std::span<Output> out) const noexcept {
  if (init_status_ != Status::Ok || !parameters_) {
    return init_status_;
  }
  if (in.size() != out.size()) {
    return Status::InvalidInput;
  }
  for (std::size_t i = 0; i < in.size(); ++i) {
    const auto res = evaluate(in[i]);
    if (res.status != Status::Ok) {
      return res.status;
    }
    out[i] = res.out;
  }
  return Status::Ok;
}

}  // namespace msis21
