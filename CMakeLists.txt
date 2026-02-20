# Author: Watosn
cmake_minimum_required(VERSION 3.25)
project(nrlmsis21_cpp VERSION 0.1.0 LANGUAGES CXX Fortran)

option(MSIS21_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)
option(MSIS21_ENABLE_FP_STRICT "Enable strict floating-point reproducibility flags" ON)
option(MSIS21_BUILD_TESTS "Build tests" ON)
option(MSIS21_USE_FORTRAN_ALT2GPH "Use Fortran alt2gph implementation for strict parity" ON)
option(MSIS21_USE_FORTRAN_REFERENCE_BACKEND "Use Fortran msiscalc backend for strict parity" ON)
set(MSIS21_FORTRAN_REFERENCE_DIR "/Users/rmw/Documents/code/old_NRL" CACHE PATH "Path to Fortran NRLMSIS reference sources")
option(MSIS21_BUILD_CLI "Build CLI example" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM.cmake)
CPMAddPackage(
  NAME fmt
  GITHUB_REPOSITORY fmtlib/fmt
  VERSION 10.2.1
)
CPMAddPackage(
  NAME spdlog
  GITHUB_REPOSITORY gabime/spdlog
  VERSION 1.14.1
  OPTIONS "SPDLOG_FMT_EXTERNAL ON"
)

add_library(msis21 STATIC
  src/msis21.cpp
  src/model.cpp
  src/parm_reader.cpp
  src/fortran_unformatted_reader.cpp
  src/log.cpp
  src/utils.cpp
  src/gfn.cpp
  src/tfn.cpp
  src/dfn.cpp
  src/calc.cpp
  src/gtd8d.cpp
)
if(MSIS21_USE_FORTRAN_ALT2GPH)
  target_sources(msis21 PRIVATE src/fortran_alt2gph.F90)
  target_compile_definitions(msis21 PRIVATE MSIS21_USE_FORTRAN_ALT2GPH=1)
endif()

if(MSIS21_USE_FORTRAN_REFERENCE_BACKEND)
  if(NOT EXISTS "${MSIS21_FORTRAN_REFERENCE_DIR}/msis_calc.F90")
    message(FATAL_ERROR "MSIS21_USE_FORTRAN_REFERENCE_BACKEND=ON but reference sources not found at ${MSIS21_FORTRAN_REFERENCE_DIR}")
  endif()
  target_sources(msis21 PRIVATE
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_constants.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_init.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_utils.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_gfn.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_tfn.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_dfn.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_calc.F90
    src/fortran_msiscalc_bridge.F90
  )
  set_source_files_properties(
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_constants.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_init.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_utils.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_gfn.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_tfn.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_dfn.F90
    ${MSIS21_FORTRAN_REFERENCE_DIR}/msis_calc.F90
    src/fortran_msiscalc_bridge.F90
    PROPERTIES COMPILE_OPTIONS "-cpp;-DDBLE"
  )
  target_compile_definitions(msis21 PRIVATE MSIS21_USE_FORTRAN_REFERENCE_BACKEND=1)
  target_compile_definitions(msis21 PRIVATE MSIS21_PARM_FILE="${CMAKE_SOURCE_DIR}/data/msis21.parm")
endif()
add_library(msis21::msis21 ALIAS msis21)

target_include_directories(msis21
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(msis21 PRIVATE spdlog::spdlog fmt::fmt)

if(MSVC)
  target_compile_options(msis21 PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4>)
  if(MSIS21_WARNINGS_AS_ERRORS)
    target_compile_options(msis21 PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/WX>)
  endif()
else()
  target_compile_options(msis21 PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)
  if(MSIS21_WARNINGS_AS_ERRORS)
    target_compile_options(msis21 PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Werror>)
  endif()
  if(MSIS21_ENABLE_FP_STRICT)
    target_compile_options(msis21 PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)
  endif()
endif()

if(MSIS21_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()

if(MSIS21_BUILD_CLI)
  add_executable(msis21_cli examples/msis_cli.cpp)
  target_link_libraries(msis21_cli PRIVATE msis21::msis21 fmt::fmt)
  add_executable(msis21_profile_cli examples/msis_profile_cli.cpp)
  target_link_libraries(msis21_profile_cli PRIVATE msis21::msis21 fmt::fmt)
endif()
