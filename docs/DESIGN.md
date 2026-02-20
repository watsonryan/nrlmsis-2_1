# NRLMSIS 2.1 C++ Design

## Layers
- Public API: `msis21::Model`, typed `Input`/`Output`, `Options`.
- Legacy compatibility: `msisinit`, `msiscalc`, `gtd8d` wrappers.
- Detail modules: parameter load, basis functions, temperature/density helpers.

## Thread Safety
- `Model` is immutable after `load_from_file`.
- `GlobeCalculator` currently caches intermediate values and is instantiated per evaluation path.

## Current State
- Parameter file loading, utility math, basis generation, and profile helper kernels are implemented.
- Full physics coupling for density/temperature outputs is in progress.
