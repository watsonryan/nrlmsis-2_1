# NRLMSIS 2.1 C++ Design

## Layers
- Public API: `msis21::Model`, typed `Input`/`Output`, `Options`.
- Legacy compatibility: `msisinit`, `msiscalc`, `gtd8d` wrappers.
- Detail modules: parameter load, basis functions, temperature/density helpers.
- Runtime is C++-only; no Fortran backend is linked.

## Thread Safety
- `Model` is immutable after `load_from_file`.
- `GlobeCalculator` currently caches intermediate values and is instantiated per evaluation path.

## Current State
- Parameter file loading, utility math, basis generation, and profile helper kernels are implemented.
- Full physics parity work is in progress; current solver is deterministic and test-covered.
