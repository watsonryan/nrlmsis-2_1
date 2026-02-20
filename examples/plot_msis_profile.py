#!/usr/bin/env python3
"""
@file plot_msis_profile.py
@brief Generate density and anomalous oxygen profile plots from msis21_profile_cli output.
@author Watosn
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt


def run_profile(cli_path: Path, start_km: float, end_km: float, step_km: float) -> tuple[list[float], list[float], list[float]]:
    cmd = [str(cli_path), f"{start_km}", f"{end_km}", f"{step_km}"]
    out = subprocess.check_output(cmd, text=True)

    alt_km: list[float] = []
    rho_g_cm3: list[float] = []
    ao_g_cm3: list[float] = []

    for line in out.splitlines():
        if not line or line.startswith("#"):
            continue
        z_s, rho_s, ao_s, status_s = line.split()
        if status_s != "0":
            continue
        alt_km.append(float(z_s))
        rho_g_cm3.append(float(rho_s))
        ao_g_cm3.append(float(ao_s))

    return alt_km, rho_g_cm3, ao_g_cm3


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot MSIS total density and anomalous oxygen profiles.")
    parser.add_argument("--cli", default="build/macos-clang-debug/msis21_profile_cli", help="Path to msis21_profile_cli binary")
    parser.add_argument("--start-km", type=float, default=2000.0)
    parser.add_argument("--end-km", type=float, default=250.0)
    parser.add_argument("--step-km", type=float, default=25.0)
    parser.add_argument("--output", default="docs/figures/msis_profile_2000_250.png")
    args = parser.parse_args()

    cli_path = Path(args.cli)
    if not cli_path.exists():
        raise FileNotFoundError(f"CLI not found: {cli_path}")

    alt_km, rho_g_cm3, ao_g_cm3 = run_profile(cli_path, args.start_km, args.end_km, args.step_km)
    if not alt_km:
        raise RuntimeError("No valid profile samples returned by CLI")

    fig, axes = plt.subplots(1, 2, figsize=(11, 7), sharey=True)

    axes[0].semilogx(rho_g_cm3, alt_km, color="#0f766e", lw=2.2)
    axes[0].set_xlabel("Total Density [g/cm^3]")
    axes[0].set_ylabel("Altitude [km]")
    axes[0].set_title("NRLMSIS Total Mass Density")
    axes[0].grid(True, which="both", ls="--", alpha=0.35)

    axes[1].semilogx(ao_g_cm3, alt_km, color="#b45309", lw=2.2)
    axes[1].set_xlabel("Anomalous O Mass Density [g/cm^3]")
    axes[1].set_title("NRLMSIS Anomalous Oxygen")
    axes[1].grid(True, which="both", ls="--", alpha=0.35)

    ymin = min(args.start_km, args.end_km)
    ymax = max(args.start_km, args.end_km)
    axes[0].set_ylim(ymin, ymax)

    fig.suptitle("NRLMSIS 2.1 Profile (2000 km to 250 km)", fontsize=13)
    fig.tight_layout()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
