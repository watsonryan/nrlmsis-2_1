# Numeric Reproducibility

## Compiler Flags
- `-ffp-contract=off` is enabled by default for Clang/GCC builds.
- No fast-math flags are enabled.

## Precision
- Core calculations use `double`.
- Parameter file is parsed as little-endian IEEE-754 64-bit doubles.

## Remaining Work
- Lock down exact operation ordering in fully coupled `msiscalc` path.
- Add platform-by-platform tolerance characterization once golden pass is complete.
