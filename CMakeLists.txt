# Author: Watosn
cmake_minimum_required(VERSION 3.25)
project(nrlmsis21_cpp VERSION 0.1.0 LANGUAGES CXX)

option(MSIS21_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)
option(MSIS21_ENABLE_FP_STRICT "Enable strict floating-point reproducibility flags" ON)
option(MSIS21_BUILD_TESTS "Build tests" ON)

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
add_library(msis21::msis21 ALIAS msis21)

target_include_directories(msis21
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(msis21 PRIVATE spdlog::spdlog fmt::fmt)

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
