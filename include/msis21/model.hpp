/**
 * @file model.hpp
 * @brief Main model object lifecycle and evaluation API.
 * @author Watosn
 */
#pragma once

#include <filesystem>

#include "msis21/options.hpp"
#include "msis21/status.hpp"
#include "msis21/types.hpp"

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
  explicit Model(Status status) : init_status_(status) {}

  Status init_status_{Status::NotInitialized};
};

}  // namespace msis21
