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

extern "C" {
void msis21_ref_gtd8d(int iyd,
                      float sec,
                      float alt,
                      float glat,
                      float glong,
                      float stl,
                      float f107a,
                      float f107,
                      float ap_daily,
                      float* d,
                      float* t,
                      int* status);
}

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

double table_tolerance(double reference) {
  const double scale = std::max(1.0, std::abs(reference));
  return 1e-3 * scale;
}

}  // namespace

TEST(GoldenVectors, HarnessLoadsAndEvaluatesDataset) {
  const auto inputs = load_inputs(msis21_data_path("msis2.1_test_in.txt"));
  const auto refs = load_reference_rows(msis21_data_path("msis2.1_test_ref_dp.txt"));

  ASSERT_EQ(inputs.size(), 200U);
  ASSERT_EQ(refs.size(), 200U);

  msis21::Options options;
  auto model = msis21::Model::load_from_file(msis21_data_path("msis21.parm"), options);

  for (std::size_t i = 0; i < inputs.size(); ++i) {
    const auto& input = inputs[i];
    const auto& ref = refs[i];
    const auto result = model.evaluate(input);
    EXPECT_EQ(result.status, msis21::Status::Ok);

    float d_fortran[10]{};
    float t_fortran[2]{};
    int fstatus = 1;
    msis21_ref_gtd8d(input.iyd,
                     static_cast<float>(input.sec),
                     static_cast<float>(input.alt_km),
                     static_cast<float>(input.glat_deg),
                     static_cast<float>(input.glon_deg),
                     static_cast<float>(input.stl_hr),
                     static_cast<float>(input.f107a),
                     static_cast<float>(input.f107),
                     static_cast<float>(input.ap),
                     d_fortran,
                     t_fortran,
                     &fstatus);
    ASSERT_EQ(fstatus, 0);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.he), d_fortran[0]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.o), d_fortran[1]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.n2), d_fortran[2]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.o2), d_fortran[3]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.ar), d_fortran[4]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.rho), d_fortran[5]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.h), d_fortran[6]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.n), d_fortran[7]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.o_anom), d_fortran[8]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.no), d_fortran[9]);
    EXPECT_FLOAT_EQ(static_cast<float>(result.out.t), t_fortran[1]);

    EXPECT_NEAR(result.out.he, ref.he, table_tolerance(ref.he));
    EXPECT_NEAR(result.out.o, ref.o, table_tolerance(ref.o));
    EXPECT_NEAR(result.out.n2, ref.n2, table_tolerance(ref.n2));
    EXPECT_NEAR(result.out.o2, ref.o2, table_tolerance(ref.o2));
    EXPECT_NEAR(result.out.ar, ref.ar, table_tolerance(ref.ar));
    EXPECT_NEAR(result.out.rho, ref.rho, table_tolerance(ref.rho));
    EXPECT_NEAR(result.out.h, ref.h, table_tolerance(ref.h));
    EXPECT_NEAR(result.out.n, ref.n, table_tolerance(ref.n));
    EXPECT_NEAR(result.out.o_anom, ref.o_anom, table_tolerance(ref.o_anom));
    EXPECT_NEAR(result.out.no, ref.no, table_tolerance(ref.no));
    EXPECT_NEAR(result.out.t, ref.t, 0.2);
  }
}
