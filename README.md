# NRLMSIS 2.1 C++20 Port

C++20 port and integration layer for NRLMSIS 2.1 with deterministic build presets and reference-parity validation.

## Architecture
```mermaid
flowchart LR
  A[Input: iyd/sec/alt/lat/lon/f107/ap] --> B[msis21::Model]
  B --> C[detail::calc]
  C --> D[detail::gfn/tfn/dfn]
  D --> E[Output: He O N2 O2 Ar rho H N O* NO T]
  C --> F[Golden tests]
  F --> G[msis2.1_test_ref_dp.txt]
```

## Build and Test
```bash
cmake --preset macos-clang-debug
cmake --build --preset macos-clang-debug
ctest --preset macos-clang-debug --output-on-failure
```

## Logging and Printing
- Runtime logging uses `spdlog`.
- Output formatting uses `fmt` where applicable.
