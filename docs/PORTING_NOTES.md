# Porting Notes

| Fortran source | C++20 module | Status |
| --- | --- | --- |
| `msis_constants.F90` | `include/msis21/detail/constants.hpp` | partial |
| `msis_utils.F90` | `include/msis21/detail/utils.hpp`, `src/utils.cpp` | partial |
| `msis_init.F90` | `include/msis21/detail/parm_reader.hpp`, `src/parm_reader.cpp` | partial |
| `msis_gfn.F90` | `include/msis21/detail/gfn.hpp`, `src/gfn.cpp` | partial |
| `msis_tfn.F90` | `include/msis21/detail/tfn.hpp`, `src/tfn.cpp` | partial |
| `msis_dfn.F90` | `include/msis21/detail/dfn.hpp`, `src/dfn.cpp` | partial |
| `msis_calc.F90` | `include/msis21/detail/calc.hpp`, `src/calc.cpp` | in progress |
| `msis_gtd8d.F90` | `include/msis21/detail/gtd8d.hpp`, `src/gtd8d.cpp` | partial |

The runtime path is C++-only. Golden-vector coverage is provided by `tests/test_golden_vectors.cpp`.
