/**
 * @file msis_profile_cli.cpp
 * @brief Emits total and anomalous-oxygen mass density profiles in g/cm^3.
 * @author Watosn
 */

#include "msis21/msis21.hpp"

#include <cstdlib>
#include <fmt/format.h>

namespace {

double oxygen_mass_g_cm3(double o_cm3) {
  constexpr double avogadro = 6.02214076e23;
  constexpr double oxygen_mass_kg = (31.9988 / 2.0) / (1.0e3 * avogadro);
  return o_cm3 * oxygen_mass_kg * 1.0e3;
}

}  // namespace

int main(int argc, char** argv) {
  double z_start = 2000.0;
  double z_end = 250.0;
  double z_step = 25.0;

  if (argc == 4) {
    z_start = std::atof(argv[1]);
    z_end = std::atof(argv[2]);
    z_step = std::atof(argv[3]);
  } else if (argc != 1) {
    fmt::print(stderr, "usage: msis_profile_cli [start_km end_km step_km]\n");
    return 1;
  }

  if (z_step <= 0.0 || z_start < z_end) {
    fmt::print(stderr, "require start_km >= end_km and step_km > 0\n");
    return 1;
  }

  auto model = msis21::Model::load_from_file("data/msis21.parm", msis21::Options{});

  fmt::print("# alt_km rho_g_cm3 o_g_cm3 ao_g_cm3 status\n");
  for (double z = z_start; z >= z_end - 1e-9; z -= z_step) {
    msis21::Input in{};
    in.iyd = 70178;
    in.sec = 43200.0;
    in.alt_km = z;
    in.glat_deg = 0.0;
    in.glon_deg = 0.0;
    in.stl_hr = 12.0;
    in.f107a = 150.0;
    in.f107 = 150.0;
    in.ap = 4.0;

    const auto out = model.evaluate(in);
    const double o_g_cm3 = oxygen_mass_g_cm3(out.out.o);
    const double ao_g_cm3 = oxygen_mass_g_cm3(out.out.o_anom);
    fmt::print("{:.1f} {:.8e} {:.8e} {:.8e} {}\n",
               z,
               out.out.rho,
               o_g_cm3,
               ao_g_cm3,
               static_cast<int>(out.status));
  }

  return 0;
}
