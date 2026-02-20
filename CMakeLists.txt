# Author: Watosn
cmake_minimum_required(VERSION 3.25)
project(nrlmsis21_cpp VERSION 0.1.0 LANGUAGES CXX Fortran)

option(MSIS21_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)
option(MSIS21_ENABLE_FP_STRICT "Enable strict floating-point reproducibility flags" ON)
option(MSIS21_BUILD_TESTS "Build tests" ON)
option(MSIS21_BUILD_CLI "Build CLI example" ON)
set(MSIS21_FORTRAN_SOURCE_DIR "/Users/rmw/Documents/code/old_NRL" CACHE PATH "Path to upstream Fortran sources")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(msis21_fortran_ref STATIC
  src/fortran_bridge.F90
  ${MSIS21_FORTRAN_SOURCE_DIR}/msis_constants.F90
  ${MSIS21_FORTRAN_SOURCE_DIR}/msis_utils.F90
  ${MSIS21_FORTRAN_SOURCE_DIR}/msis_init.F90
  ${MSIS21_FORTRAN_SOURCE_DIR}/msis_gfn.F90
  ${MSIS21_FORTRAN_SOURCE_DIR}/msis_tfn.F90
  ${MSIS21_FORTRAN_SOURCE_DIR}/msis_dfn.F90
  ${MSIS21_FORTRAN_SOURCE_DIR}/msis_calc.F90
  ${MSIS21_FORTRAN_SOURCE_DIR}/msis_gtd8d.F90
)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Flang|LLVMFlang")
  target_compile_options(msis21_fortran_ref PRIVATE -cpp -DDBLE)
endif()

add_library(msis21 STATIC
  src/msis21.cpp
  src/model.cpp
  src/parm_reader.cpp
  src/fortran_unformatted_reader.cpp
  src/fortran_backend.cpp
  src/utils.cpp
  src/gfn.cpp
  src/tfn.cpp
  src/dfn.cpp
  src/calc.cpp
  src/gtd8d.cpp
)
add_library(msis21::msis21 ALIAS msis21)

target_include_directories(msis21
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(msis21 PRIVATE msis21_fortran_ref)

if(MSVC)
  target_compile_options(msis21 PRIVATE /W4)
  if(MSIS21_WARNINGS_AS_ERRORS)
    target_compile_options(msis21 PRIVATE /WX)
  endif()
else()
  target_compile_options(msis21 PRIVATE -Wall -Wextra -Wpedantic)
  if(MSIS21_WARNINGS_AS_ERRORS)
    target_compile_options(msis21 PRIVATE -Werror)
  endif()
  if(MSIS21_ENABLE_FP_STRICT)
    target_compile_options(msis21 PRIVATE -ffp-contract=off)
  endif()
endif()

if(MSIS21_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

if(MSIS21_BUILD_CLI)
  add_executable(msis21_cli examples/msis_cli.cpp)
  target_link_libraries(msis21_cli PRIVATE msis21::msis21)
endif()
