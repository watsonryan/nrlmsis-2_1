/**
 * @file fortran_backend.cpp
 * @brief Bridge to upstream Fortran reference implementation.
 * @author Watosn
 */

#include "msis21/detail/fortran_backend.hpp"

#include <array>
#include <cmath>
#include <string>

#include "msis21/detail/constants.hpp"

namespace msis21::detail {

namespace {

extern "C" {
void msis21_ref_init(const char* parm_dir, int* status);
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

bool g_fortran_ready = false;

}  // namespace

Status initialize_fortran_backend(const std::filesystem::path& parm_file_path) {
  if (!std::filesystem::exists(parm_file_path)) {
    return Status::ParmFileNotFound;
  }

  auto parm_dir = parm_file_path.parent_path().string();
  if (!parm_dir.empty() && parm_dir.back() != std::filesystem::path::preferred_separator) {
    parm_dir.push_back(std::filesystem::path::preferred_separator);
  }

  int status = 1;
  msis21_ref_init(parm_dir.c_str(), &status);
  g_fortran_ready = (status == 0);
  return g_fortran_ready ? Status::Ok : Status::NumericalError;
}

bool fortran_backend_available() { return g_fortran_ready; }

CalcResult evaluate_fortran_backend(const Input& in) {
  if (!g_fortran_ready) {
    return CalcResult{.status = Status::NotInitialized};
  }

  std::array<float, 10> d{};
  std::array<float, 2> t{};
  int status = 1;
  msis21_ref_gtd8d(in.iyd,
                   static_cast<float>(in.sec),
                   static_cast<float>(in.alt_km),
                   static_cast<float>(in.glat_deg),
                   static_cast<float>(in.glon_deg),
                   static_cast<float>(in.stl_hr),
                   static_cast<float>(in.f107a),
                   static_cast<float>(in.f107),
                   static_cast<float>(in.ap),
                   d.data(),
                   t.data(),
                   &status);

  if (status != 0) {
    return CalcResult{.status = Status::NumericalError};
  }

  for (const auto value : d) {
    if (!std::isfinite(value)) {
      return CalcResult{.status = Status::NumericalError};
    }
  }

  Output out{};
  out.he = static_cast<double>(d[0]);
  out.o = static_cast<double>(d[1]);
  out.n2 = static_cast<double>(d[2]);
  out.o2 = static_cast<double>(d[3]);
  out.ar = static_cast<double>(d[4]);
  out.rho = static_cast<double>(d[5]);
  out.h = static_cast<double>(d[6]);
  out.n = static_cast<double>(d[7]);
  out.o_anom = static_cast<double>(d[8]);
  out.no = static_cast<double>(d[9]);
  out.t = static_cast<double>(t[1]);

  return CalcResult{.status = Status::Ok, .out = out};
}

}  // namespace msis21::detail
