"""Self-consistent Hartree/DFT models for helium (Thijssen §5.3)."""

from __future__ import annotations

import argparse
import pathlib
import sys
from dataclasses import dataclass

import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(str(pathlib.Path(__file__).resolve().parents[2]))
    from comp_physics_python.ch5.common import (  # type: ignore
        RadialGrid,
        ceperley_alder,
        densities_from_u,
        exact_exchange_from_density,
        hartree_potential,
        lda_exchange,
        shoot_bound_state,
    )
else:  # pragma: no cover
    from .common import (
        RadialGrid,
        ceperley_alder,
        densities_from_u,
        exact_exchange_from_density,
        hartree_potential,
        lda_exchange,
        shoot_bound_state,
    )


@dataclass
class ModelConfig:
    mode: str
    step: float
    rmax: float
    max_iter: int = 200
    tol: float = 1e-9

    def __post_init__(self) -> None:
        modes = {"hartree", "slater", "exactx", "lda"}
        if self.mode not in modes:
            raise ValueError(f"Mode must be one of {modes}")


def run_scf(cfg: ModelConfig) -> dict:
    grid = RadialGrid(cfg.step, cfg.rmax)
    Z = 2.0
    electrons = 2.0
    hart = np.zeros_like(grid.r)
    vx = np.zeros_like(grid.r)
    vc = np.zeros_like(grid.r)
    exc_energy = 0.0
    corr_energy = 0.0

    total_energy = 0.0
    for iteration in range(cfg.max_iter):
        potential = -Z / grid.safe_r + hart - vx - vc
        energy, u = shoot_bound_state(potential, grid, tail_z=Z, energy_bounds=(-2.0, -0.3))
        dens_single, dens_total = densities_from_u(u, grid, electrons=1.0)
        dens_total *= electrons

        new_hart = hartree_potential(dens_single, grid, target_charge=1.0)

        if cfg.mode == "hartree":
            new_vx = np.zeros_like(vx)
            exc_energy = 0.0
        elif cfg.mode == "slater":
            new_vx, exc_energy = lda_exchange(dens_total, grid)
        elif cfg.mode == "exactx":
            new_vx, exc_energy = exact_exchange_from_density(dens_total, grid)
        else:  # lda
            new_vx, exc_energy = lda_exchange(dens_total, grid)

        if cfg.mode == "lda":
            eps_c, vc_local = ceperley_alder(dens_total)
            new_vc = vc_local
            corr_energy = grid.integrate(eps_c * dens_total * grid.surface)
        else:
            new_vc = np.zeros_like(vc)
            corr_energy = 0.0

        hart_energy = grid.integrate(dens_single * new_hart * grid.surface)
        new_total = electrons * energy - hart_energy + exc_energy + corr_energy

        if abs(new_total - total_energy) < cfg.tol:
            total_energy = new_total
            hart = new_hart
            vx = new_vx
            vc = new_vc
            break

        total_energy = new_total
        hart = new_hart
        vx = new_vx
        vc = new_vc
    else:
        raise RuntimeError("SCF did not converge")

    return {
        "grid": grid,
        "eigenvalue": energy,
        "total_energy": total_energy,
        "hartree": hart,
        "vx": vx,
        "vc": vc,
        "mode": cfg.mode,
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["hartree", "slater", "exactx", "lda"], help="approximation")
    parser.add_argument("--step", type=float, default=0.02, help="radial step")
    parser.add_argument("--rmax", type=float, default=30.0, help="maximum radius")
    parser.add_argument("--max-iter", type=int, default=200)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--out", type=str, default="helium_potentials.npz")
    args = parser.parse_args(argv)

    cfg = ModelConfig(mode=args.mode, step=args.step, rmax=args.rmax, max_iter=args.max_iter, tol=args.tol)
    result = run_scf(cfg)
    print(f"Converged {args.mode} model: ε0 = {result['eigenvalue']:.8f} Ha")
    print(f"Total energy E = {result['total_energy']:.8f} Ha")
    np.savez(
        args.out,
        r=result["grid"].r,
        hartree=result["hartree"],
        v_exchange=result["vx"],
        v_corr=result["vc"],
    )
    print(f"Saved Hartree/DFT potentials to {args.out}")


if __name__ == "__main__":
    main()
