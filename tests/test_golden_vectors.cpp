/**
 * @file test_golden_vectors.cpp
 * @brief Golden vector regression harness for NRLMSIS 2.1.
 * @author Watosn
 */

#include "msis21/msis21.hpp"

#include <filesystem>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "test_paths.hpp"

namespace {

std::vector<msis21::Input> load_inputs(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open input file");
  }

  std::string line;
  std::getline(in, line);

  std::vector<msis21::Input> rows;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream iss(line);
    msis21::Input row;
    iss >> row.iyd >> row.sec >> row.alt_km >> row.glat_deg >> row.glon_deg >> row.stl_hr >> row.f107a >>
        row.f107 >> row.ap;
    if (!iss.fail()) {
      rows.push_back(row);
    }
  }
  return rows;
}

struct RefRow {
  double he{};
  double o{};
  double n2{};
  double o2{};
  double ar{};
  double rho{};
  double h{};
  double n{};
  double o_anom{};
  double no{};
  double t{};
};

std::vector<RefRow> load_reference_rows(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open reference file");
  }

  std::string header;
  std::getline(in, header);

  std::string line;
  std::vector<RefRow> rows;
  while (std::getline(in, line)) {
    if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos) {
      continue;
    }
    std::istringstream iss(line);
    int iyd = 0;
    int sec = 0;
    double alt = 0.0;
    double glat = 0.0;
    double glong = 0.0;
    double stl = 0.0;
    double f107a = 0.0;
    double f107 = 0.0;
    double ap = 0.0;
    RefRow row;
    if (iss >> iyd >> sec >> alt >> glat >> glong >> stl >> f107a >> f107 >> ap >> row.he >> row.o >>
            row.n2 >> row.o2 >> row.ar >> row.rho >> row.h >> row.n >> row.o_anom >> row.no >> row.t) {
      rows.push_back(row);
    }
  }
  return rows;
}

}  // namespace

TEST(GoldenVectors, HarnessLoadsAndEvaluatesDataset) {
  const auto inputs = load_inputs(msis21_data_path("msis2.1_test_in.txt"));
  const auto refs = load_reference_rows(msis21_data_path("msis2.1_test_ref_dp.txt"));

  ASSERT_EQ(inputs.size(), 200U);
  ASSERT_EQ(refs.size(), 200U);

  msis21::Options options;
  auto model = msis21::Model::load_from_file(msis21_data_path("msis21.parm"), options);

  auto rel_err = [](double got, double ref) {
    const double denom = std::max(std::abs(ref), 1e-30);
    return std::abs(got - ref) / denom;
  };
  constexpr double kTolHe = 2.0e-3;
  constexpr double kTolO = 2.0e-3;
  constexpr double kTolN2 = 2.0e-3;
  constexpr double kTolO2 = 2.0e-3;
  constexpr double kTolAr = 2.0e-3;
  constexpr double kTolRho = 2.0e-3;
  constexpr double kTolH = 2.0e-3;
  constexpr double kTolN = 2.0e-3;
  constexpr double kTolOAnom = 2.0e-3;
  constexpr double kTolNo = 2.0e-3;
  constexpr double kTolT = 1.0e-4;

  for (std::size_t i = 0; i < inputs.size(); ++i) {
    const auto& input = inputs[i];
    const auto& ref = refs[i];
    const auto result = model.evaluate(input);
    EXPECT_EQ(result.status, msis21::Status::Ok);

    ASSERT_TRUE(std::isfinite(result.out.he));
    ASSERT_TRUE(std::isfinite(result.out.o));
    ASSERT_TRUE(std::isfinite(result.out.n2));
    ASSERT_TRUE(std::isfinite(result.out.o2));
    ASSERT_TRUE(std::isfinite(result.out.ar));
    ASSERT_TRUE(std::isfinite(result.out.rho));
    ASSERT_TRUE(std::isfinite(result.out.h));
    ASSERT_TRUE(std::isfinite(result.out.n));
    ASSERT_TRUE(std::isfinite(result.out.o_anom));
    ASSERT_TRUE(std::isfinite(result.out.no));
    ASSERT_TRUE(std::isfinite(result.out.t));

    EXPECT_LE(rel_err(result.out.he, ref.he), kTolHe);
    EXPECT_LE(rel_err(result.out.o, ref.o), kTolO);
    EXPECT_LE(rel_err(result.out.n2, ref.n2), kTolN2);
    EXPECT_LE(rel_err(result.out.o2, ref.o2), kTolO2);
    EXPECT_LE(rel_err(result.out.ar, ref.ar), kTolAr);
    EXPECT_LE(rel_err(result.out.rho, ref.rho), kTolRho);
    EXPECT_LE(rel_err(result.out.h, ref.h), kTolH);
    EXPECT_LE(rel_err(result.out.n, ref.n), kTolN);
    EXPECT_LE(rel_err(result.out.o_anom, ref.o_anom), kTolOAnom);
    EXPECT_LE(rel_err(result.out.no, ref.no), kTolNo);
    EXPECT_LE(rel_err(result.out.t, ref.t), kTolT);
  }
}

TEST(GoldenVectors, AtomicOxygenRegressionCase) {
  msis21::Options options;
  auto model = msis21::Model::load_from_file(msis21_data_path("msis21.parm"), options);
  msis21::Input in{};
  in.iyd = 70338;
  in.sec = 71881.0;
  in.alt_km = 52.1;
  in.glat_deg = 66.3;
  in.glon_deg = 259.3;
  in.stl_hr = 13.25;
  in.f107a = 158.3;
  in.f107 = 152.4;
  in.ap = 5.0;
  const auto result = model.evaluate(in);
  ASSERT_EQ(result.status, msis21::Status::Ok);
  constexpr double kExpectedO = 3.6130985269374433e+09;
  const double rel = std::abs(result.out.o - kExpectedO) / kExpectedO;
  EXPECT_LE(rel, 1.0e-12);
}
