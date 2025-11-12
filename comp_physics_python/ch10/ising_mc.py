"""Metropolis Monte Carlo for the 2D Ising model."""

from __future__ import annotations

import argparse
import numpy as np


def metropolis(size: int, j: float, h: float, steps: int, burn_in: int = 2000) -> dict:
    spins = np.ones((size, size), dtype=int)
    energy = -2 * j * size * size
    magnetisation = size * size
    rng = np.random.default_rng()
    def local_field(x: int, y: int) -> int:
        return spins[(x + 1) % size, y] + spins[(x - 1) % size, y] + spins[x, (y + 1) % size] + spins[x, (y - 1) % size]
    energies = []
    mags = []
    for step in range(steps + burn_in):
        for _ in range(size * size):
            x = rng.integers(size)
            y = rng.integers(size)
            delta = 2 * spins[x, y] * (j * local_field(x, y) + h)
            if delta <= 0 or rng.random() < np.exp(-delta):
                spins[x, y] *= -1
                energy += delta
                magnetisation += 2 * spins[x, y]
        if step >= burn_in:
            energies.append(energy)
            mags.append(magnetisation)
    energies = np.array(energies)
    mags = np.array(mags)
    beta = 1.0
    return {
        "energy": energies.mean() / (size * size),
        "specific_heat": beta**2 * energies.var() / (size * size),
        "magnetisation": mags.mean() / (size * size),
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=32)
    parser.add_argument("--J", type=float, default=0.4406868, help="Reduced coupling (K/kT)")
    parser.add_argument("--H", type=float, default=0.0, help="Reduced field (B/kT)")
    parser.add_argument("--steps", type=int, default=5000)
    args = parser.parse_args(argv)
    result = metropolis(args.size, args.J, args.H, args.steps)
    print(f"<E>/N = {result['energy']:.4f}")
    print(f"C_v = {result['specific_heat']:.4f}")
    print(f"<m> = {result['magnetisation']:.4f}")


if __name__ == "__main__":
    main()
