"""Metropolis Monte Carlo for the Lennard-Jones fluid (Thijssen §10.3.2)."""

from __future__ import annotations

import argparse
import pathlib
import sys
from dataclasses import dataclass

import numpy as np


@dataclass
class LJMCConfig:
    n_atoms: int
    density: float
    temperature: float
    cutoff: float = 2.5
    max_displacement: float = 0.3
    target_acceptance: float = 0.4

    def __post_init__(self) -> None:
        self.beta = 1.0 / self.temperature
        self.box_length = (self.n_atoms / self.density) ** (1.0 / 3.0)
        cells = round((self.n_atoms / 4) ** (1.0 / 3.0))
        if 4 * cells**3 != self.n_atoms:
            raise ValueError("n_atoms must be 4*cells^3 for FCC initialisation")
        self.cells = cells


def make_fcc(cfg: LJMCConfig) -> np.ndarray:
    a = cfg.box_length / cfg.cells
    basis = np.array(
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]], dtype=float
    )
    positions = []
    for ix in range(cfg.cells):
        for iy in range(cfg.cells):
            for iz in range(cfg.cells):
                origin = np.array([ix, iy, iz], dtype=float)
                positions.extend(a * (origin + basis))
    return np.array(positions)


def minimum_image(vec: np.ndarray, box: float) -> np.ndarray:
    return vec - box * np.rint(vec / box)


def pair_energy(displacements: np.ndarray, cutoff2: float) -> float:
    r2 = np.sum(displacements * displacements, axis=1)
    mask = r2 < cutoff2
    if not np.any(mask):
        return 0.0
    r2 = r2[mask]
    inv_r2 = 1.0 / r2
    inv_r6 = inv_r2**3
    inv_r12 = inv_r6**2
    return np.sum(4.0 * (inv_r12 - inv_r6))


def total_energy(positions: np.ndarray, cfg: LJMCConfig) -> float:
    energy = 0.0
    cutoff2 = cfg.cutoff**2
    for i in range(cfg.n_atoms - 1):
        disp = positions[i + 1 :] - positions[i]
        disp = minimum_image(disp, cfg.box_length)
        energy += pair_energy(disp, cutoff2)
    return energy


def metropolis(cfg: LJMCConfig, steps: int, burn_in: int = 1000, seed: int | None = None) -> dict:
    rng = np.random.default_rng(seed)
    positions = make_fcc(cfg)
    energy = total_energy(positions, cfg)
    cutoff2 = cfg.cutoff**2
    energies = []
    acceptance = 0

    for step in range(1, steps + burn_in + 1):
        idx = rng.integers(cfg.n_atoms)
        trial = positions[idx] + rng.uniform(-cfg.max_displacement, cfg.max_displacement, size=3)
        trial %= cfg.box_length
        disp_old = positions - positions[idx]
        disp_new = positions - trial
        disp_old = minimum_image(disp_old, cfg.box_length)
        disp_new = minimum_image(disp_new, cfg.box_length)
        old_e = pair_energy(np.delete(disp_old, idx, axis=0), cutoff2)
        new_e = pair_energy(np.delete(disp_new, idx, axis=0), cutoff2)
        delta = new_e - old_e
        if delta < 0.0 or rng.random() < np.exp(-cfg.beta * delta):
            positions[idx] = trial
            energy += delta
            acceptance += 1
        if step % 100 == 0:
            acc_ratio = acceptance / 100
            if acc_ratio < cfg.target_acceptance:
                cfg.max_displacement *= 0.9
            else:
                cfg.max_displacement *= 1.1
            acceptance = 0
        if step > burn_in:
            energies.append(energy)
    energies_np = np.array(energies)
    return {
        "energy": energies_np.mean(),
        "heat_capacity": cfg.beta**2 * energies_np.var(),
        "positions": positions,
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=256)
    parser.add_argument("--density", type=float, default=0.9)
    parser.add_argument("--temp", type=float, default=1.2)
    parser.add_argument("--cutoff", type=float, default=2.5)
    parser.add_argument("--steps", type=int, default=20000)
    parser.add_argument("--burn", type=int, default=2000)
    args = parser.parse_args(argv)

    cfg = LJMCConfig(
        n_atoms=args.n,
        density=args.density,
        temperature=args.temp,
        cutoff=args.cutoff,
    )
    result = metropolis(cfg, steps=args.steps, burn_in=args.burn)
    print(f"<E>/N = {result['energy']/cfg.n_atoms:.4f}")
    print(f"C_v/N = {result['heat_capacity']/cfg.n_atoms:.4f}")
    np.save(pathlib.Path("lj_mc_positions.npy"), result["positions"])


if __name__ == "__main__":
    main()
