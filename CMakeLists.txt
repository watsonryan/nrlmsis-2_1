# Author: Watosn
cmake_minimum_required(VERSION 3.25)
project(nrlmsis21_cpp VERSION 0.1.0 LANGUAGES CXX)

option(MSIS21_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)
option(MSIS21_ENABLE_FP_STRICT "Enable strict floating-point reproducibility flags" ON)
option(MSIS21_BUILD_TESTS "Build tests" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(msis21 STATIC
  src/msis21.cpp
  src/model.cpp
  src/parm_reader.cpp
  src/fortran_unformatted_reader.cpp
  src/utils.cpp
  src/gfn.cpp
)
add_library(msis21::msis21 ALIAS msis21)

target_include_directories(msis21
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

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
