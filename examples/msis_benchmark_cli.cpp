/**
 * @file msis_benchmark_cli.cpp
 * @brief Microbenchmark for batch model evaluation throughput.
 * @author Watosn
 */

#include "msis21/msis21.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <fmt/format.h>

namespace {

std::vector<msis21::Input> load_inputs(const std::string& path) {
  std::ifstream in(path);
  std::vector<msis21::Input> rows;
  if (!in) {
    return rows;
  }
  std::string line;
  std::getline(in, line);  // header
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    msis21::Input row{};
    if (std::sscanf(line.c_str(),
                    "%d %lf %lf %lf %lf %lf %lf %lf %lf",
                    &row.iyd,
                    &row.sec,
                    &row.alt_km,
                    &row.glat_deg,
                    &row.glon_deg,
                    &row.stl_hr,
                    &row.f107a,
                    &row.f107,
                    &row.ap) == 9) {
      rows.push_back(row);
    }
  }
  return rows;
}

}  // namespace

int main(int argc, char** argv) {
  int repeats = 2000;
  if (argc >= 2) {
    repeats = std::max(1, std::atoi(argv[1]));
  }

  const auto inputs = load_inputs("data/msis2.1_test_in.txt");
  if (inputs.empty()) {
    fmt::print(stderr, "failed to load benchmark input set\n");
    return 1;
  }

  auto model = msis21::Model::load_from_file("data/msis21.parm", msis21::Options{});

  // Warmup to stabilize caches and branch predictors.
  double sink = 0.0;
  for (const auto& in : inputs) {
    const auto out = model.evaluate(in);
    sink += out.out.rho + out.out.o;
  }

  const auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < repeats; ++r) {
    for (const auto& in : inputs) {
      const auto out = model.evaluate(in);
      sink += out.out.rho + out.out.o;
    }
  }
  const auto t1 = std::chrono::steady_clock::now();

  const auto total_calls = static_cast<double>(repeats) * static_cast<double>(inputs.size());
  const auto elapsed_ns =
      static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
  const double ns_per_call = elapsed_ns / total_calls;
  const double calls_per_sec = 1.0e9 / ns_per_call;

  fmt::print("repeats={} rows={} total_calls={:.0f}\n", repeats, inputs.size(), total_calls);
  fmt::print("elapsed_s={:.6f} ns_per_call={:.3f} calls_per_sec={:.3f}\n",
             elapsed_ns * 1.0e-9,
             ns_per_call,
             calls_per_sec);
  fmt::print("sink={:.16e}\n", sink);
  return std::isfinite(sink) ? 0 : 2;
}
