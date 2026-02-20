/**
 * @file model.hpp
 * @brief Main model object lifecycle and evaluation API.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "msis21/options.hpp"
#include "msis21/status.hpp"
#include "msis21/types.hpp"

namespace msis21::detail {
struct Parameters;
}

namespace msis21 {

struct Result {
  Status status{Status::Ok};
  Output out{};
};

class Model {
 public:
  static Model load_from_file(const std::filesystem::path& parm_path, Options options);

  [[nodiscard]] Result evaluate(const Input& input) const noexcept;

  class Scratch {};

 [[nodiscard]] Result evaluate(const Input& input, Scratch& scratch) const noexcept;

 private:
  explicit Model(Status status, std::shared_ptr<const detail::Parameters> parameters)
      : init_status_(status), parameters_(std::move(parameters)) {}

  Status init_status_{Status::NotInitialized};
  std::shared_ptr<const detail::Parameters> parameters_{};
};

}  // namespace msis21
