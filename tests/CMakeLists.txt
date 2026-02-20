# Author: Watosn
include(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/CPM.cmake)

CPMAddPackage(
  NAME googletest
  GITHUB_REPOSITORY google/googletest
  VERSION 1.17.0
  OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt ON"
)

add_executable(msis21_tests
  test_smoke.cpp
  test_fortran_reader.cpp
  test_parm_reader.cpp
  test_utils.cpp
  test_gfn.cpp
  test_tfn.cpp
  test_dfn.cpp
  test_calc.cpp
  test_gtd8d.cpp
  test_golden_vectors.cpp
  test_invariants.cpp
)

target_link_libraries(msis21_tests
  PRIVATE
    msis21::msis21
    GTest::gtest_main
)
target_compile_definitions(msis21_tests PRIVATE MSIS21_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

include(GoogleTest)
gtest_discover_tests(msis21_tests)
