"""Restricted Hartree–Fock for the helium atom (Thijssen §4.3.2)."""

from __future__ import annotations

import argparse
import math
import pathlib
import sys
from dataclasses import dataclass
from typing import Tuple

import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(str(pathlib.Path(__file__).resolve().parents[2]))
    from comp_physics_python.ch3._linalg import generalized_eigh  # type: ignore
else:  # pragma: no cover
    from ..ch3._linalg import generalized_eigh


@dataclass
class HeliumConfig:
    alphas: Tuple[float, ...] = (0.298073, 1.242567, 5.782948, 38.47497)
    convergence: float = 1e-11
    max_iterations: int = 200

    def __post_init__(self) -> None:
        self.alphas = tuple(float(a) for a in self.alphas)
        self.n = len(self.alphas)
        self.pi = math.pi


def overlap_matrix(cfg: HeliumConfig) -> np.ndarray:
    n = cfg.n
    S = np.zeros((n, n), dtype=float)
    for r in range(n):
        for s in range(n):
            denom = cfg.alphas[r] + cfg.alphas[s]
            factor = cfg.pi / denom
            S[r, s] = factor * math.sqrt(factor)
    return S


def kinetic_element(a: float, b: float, pi: float) -> float:
    alph = a + b
    factor = pi / alph
    return 3.0 * factor * math.sqrt(factor) * a * b / alph


def nuclear_coulomb_element(a: float, b: float, pi: float, Z: float) -> float:
    return -2.0 * Z * pi / (a + b)


def hamilton_matrix(cfg: HeliumConfig, Z: float = 2.0) -> np.ndarray:
    n = cfg.n
    H = np.zeros((n, n), dtype=float)
    for r in range(n):
        for s in range(n):
            H[r, s] = kinetic_element(cfg.alphas[r], cfg.alphas[s], cfg.pi) + nuclear_coulomb_element(
                cfg.alphas[r], cfg.alphas[s], cfg.pi, Z
            )
    return H


def two_electron_tensor(cfg: HeliumConfig) -> np.ndarray:
    n = cfg.n
    q = np.zeros((n, n, n, n), dtype=float)
    pref = 2.0 * (cfg.pi**2) * math.sqrt(cfg.pi)
    for r in range(n):
        for s in range(n):
            for t in range(n):
                for u in range(n):
                    a = cfg.alphas[r] + cfg.alphas[s]
                    b = cfg.alphas[t] + cfg.alphas[u]
                    q[r, s, t, u] = pref / math.sqrt(a + b) / (a * b)
    return q


def density_matrix(c: np.ndarray, S: np.ndarray) -> np.ndarray:
    norm = float(c @ (S @ c))
    dens = 2.0 * np.outer(c, c) / norm
    return dens


def build_g_matrix(dens: np.ndarray, q: np.ndarray) -> np.ndarray:
    n = dens.shape[0]
    G = np.zeros((n, n), dtype=float)
    for r in range(n):
        for s in range(n):
            val = 0.0
            for t in range(n):
                for u in range(t):
                    val += dens[t, u] * q[r, s, t, u]
                val += 0.5 * dens[t, t] * q[r, s, t, t]
            G[r, s] = val
    return G


def total_energy(eps0: float, dens: np.ndarray, H: np.ndarray) -> float:
    n = H.shape[0]
    energy = eps0
    for i in range(n):
        for j in range(i):
            energy += dens[i, j] * H[i, j]
        energy += 0.5 * dens[i, i] * H[i, i]
    return energy


def hartree_fock_helium(cfg: HeliumConfig) -> dict:
    S = overlap_matrix(cfg)
    H = hamilton_matrix(cfg)
    q = two_electron_tensor(cfg)

    c = np.ones(cfg.n, dtype=float)
    dens = density_matrix(c, S)

    energy = 0.0
    for _ in range(cfg.max_iterations):
        G = build_g_matrix(dens, q)
        F = H + G
        evals, evecs = generalized_eigh(F, S)
        c = evecs[:, 0]
        dens = density_matrix(c, S)
        new_energy = total_energy(evals[0], dens, H)
        if abs(new_energy - energy) < cfg.convergence:
            energy = new_energy
            break
        energy = new_energy
    else:
        raise RuntimeError("SCF did not converge for helium")

    grid = np.linspace(0.0, 4.0, 400)
    psi = np.zeros_like(grid)
    for coeff, alpha in zip(c, cfg.alphas):
        psi += coeff * np.exp(-alpha * grid * grid)

    return {
        "energy": energy,
        "orbital_energy": evals[0],
        "coefficients": c,
        "grid_r": grid,
        "radial_wavefunction": psi,
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-iter", type=int, default=200)
    parser.add_argument("--tol", type=float, default=1e-11)
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--out", type=str, default="WaveFunc_he_py")
    args = parser.parse_args(argv)

    cfg = HeliumConfig(convergence=args.tol, max_iterations=args.max_iter)
    result = hartree_fock_helium(cfg)
    print(f"Variational HF energy = {result['energy']:.8f} Ha")
    print("Coefficients:", result["coefficients"])

    grid = np.linspace(0.0, 4.0, args.samples)
    psi = np.interp(grid, result["grid_r"], result["radial_wavefunction"])
    out = pathlib.Path(args.out)
    np.savetxt(out, np.column_stack([grid, psi]), fmt="%12.6f")
    print(f"Radial wave function written to {out}")


if __name__ == "__main__":
    main()
