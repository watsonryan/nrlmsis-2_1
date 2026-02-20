/**
 * @file log.cpp
 * @brief Logging bootstrap for msis21 runtime diagnostics.
 * @author Watosn
 */

#include "msis21/detail/log.hpp"

#include <cstdlib>
#include <mutex>

#include <spdlog/spdlog.h>

namespace msis21::detail {

void configure_logging_once() {
  static std::once_flag once;
  std::call_once(once, [] {
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    const char* env_level = std::getenv("MSIS21_LOG_LEVEL");
    if (env_level != nullptr) {
      spdlog::set_level(spdlog::level::from_str(env_level));
      return;
    }
    spdlog::set_level(spdlog::level::info);
  });
}

}  // namespace msis21::detail
