# Author: Watosn
include_guard(GLOBAL)

include(${CMAKE_CURRENT_LIST_DIR}/CPM.cmake)

# fmt formatting library
CPMAddPackage(
  NAME fmt
  GITHUB_REPOSITORY fmtlib/fmt
  GIT_TAG 10.2.1
)

# spdlog logger
CPMAddPackage(
  NAME spdlog
  GITHUB_REPOSITORY gabime/spdlog
  VERSION 1.14.1
  OPTIONS "SPDLOG_FMT_EXTERNAL ON"
)

# Eigen 5 linear algebra headers
set(MSIS21_EIGEN_REPOSITORY "https://gitlab.com/libeigen/eigen.git" CACHE STRING "Eigen repository URL")
set(MSIS21_EIGEN_TAG "5.0.0" CACHE STRING "Eigen git tag or branch")
CPMAddPackage(
  NAME eigen
  GIT_REPOSITORY ${MSIS21_EIGEN_REPOSITORY}
  GIT_TAG ${MSIS21_EIGEN_TAG}
)

set(MSIS21_EIGEN_TARGET "")
if(TARGET Eigen3::Eigen)
  set(MSIS21_EIGEN_TARGET Eigen3::Eigen)
elseif(TARGET eigen)
  set(MSIS21_EIGEN_TARGET eigen)
endif()

if(MSIS21_BUILD_TESTS)
  CPMAddPackage(
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    VERSION 1.17.0
    OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt ON"
  )
endif()
