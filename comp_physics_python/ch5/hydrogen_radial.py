"""Hydrogen 1s radial solution via Numerov (Thijssen §5.3.1)."""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(str(pathlib.Path(__file__).resolve().parents[2]))
    from comp_physics_python.ch5.common import RadialGrid, shoot_bound_state  # type: ignore
else:  # pragma: no cover
    from .common import RadialGrid, shoot_bound_state


def solve_hydrogen(h: float, rmax: float) -> dict:
    grid = RadialGrid(step=h, max_radius=rmax)
    potential = -1.0 / grid.safe_r
    energy, u = shoot_bound_state(potential, grid, tail_z=1.0, energy_bounds=(-0.8, -0.3))
    _, density = np.zeros_like(u), np.zeros_like(u)
    density_single = np.zeros_like(u)
    density_single[1:] = (u[1:] / grid.safe_r[1:]) ** 2 / (4.0 * np.pi)
    density_single[0] = density_single[1]
    psi = u / grid.safe_r
    return {
        "grid": grid.r,
        "energy": energy,
        "u": u,
        "psi": psi,
        "density": density_single,
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=float, default=0.02, help="radial step h")
    parser.add_argument("--rmax", type=float, default=40.0, help="maximum radius")
    parser.add_argument("--out", type=str, default="hydrogen_wavefunction.dat")
    args = parser.parse_args(argv)

    result = solve_hydrogen(args.step, args.rmax)
    print(f"1s eigenvalue ≈ {result['energy']:.8f} Ha (exact -0.5)")
    data = np.column_stack([result["grid"], result["psi"], result["density"]])
    np.savetxt(args.out, data, fmt="%12.6f", header="r  psi(r)  density(r)")
    print(f"Saved wavefunction samples to {args.out}")


if __name__ == "__main__":
    main()
