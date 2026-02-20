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
)

target_link_libraries(msis21_tests
  PRIVATE
    msis21::msis21
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(msis21_tests)
