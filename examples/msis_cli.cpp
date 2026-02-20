/**
 * @file msis_cli.cpp
 * @brief Simple CLI for manual model invocation.
 * @author Watosn
 */

#include "msis21/msis21.hpp"

#include <cstdlib>
#include <fmt/format.h>

int main(int argc, char** argv) {
  if (argc != 10) {
    fmt::print(stderr,
               "usage: msis_cli <iyd> <sec> <alt_km> <glat_deg> <glon_deg> <stl_hr> <f107a> "
               "<f107> <ap>\n");
    return 1;
  }

  msis21::Input in;
  in.iyd = std::atoi(argv[1]);
  in.sec = std::atof(argv[2]);
  in.alt_km = std::atof(argv[3]);
  in.glat_deg = std::atof(argv[4]);
  in.glon_deg = std::atof(argv[5]);
  in.stl_hr = std::atof(argv[6]);
  in.f107a = std::atof(argv[7]);
  in.f107 = std::atof(argv[8]);
  in.ap = std::atof(argv[9]);

  auto model = msis21::Model::load_from_file("data/msis21.parm", msis21::Options{});
  const auto out = model.evaluate(in);
  fmt::print("status={}\n", static_cast<int>(out.status));
  fmt::print(
             "he={:.16e} o={:.16e} n2={:.16e} o2={:.16e} ar={:.16e} rho={:.16e} "
             "h={:.16e} n={:.16e} oa={:.16e} no={:.16e} t={:.16e}\n",
             out.out.he,
             out.out.o,
             out.out.n2,
             out.out.o2,
             out.out.ar,
             out.out.rho,
             out.out.h,
             out.out.n,
             out.out.o_anom,
             out.out.no,
             out.out.t);
  return 0;
}
