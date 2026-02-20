#!/usr/bin/env bash
# Author: Watosn
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${1:-${ROOT_DIR}/build/macos-clang-debug}"

if [[ ! -d "${BUILD_DIR}" ]]; then
  exit 0
fi

rm -f "${BUILD_DIR}/.ninja_log" "${BUILD_DIR}/.ninja_deps"
