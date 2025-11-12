"""Restricted Hartree–Fock for the H₂ molecule (Thijssen §4.3.2)."""

from __future__ import annotations

import argparse
import math
import pathlib
import sys
from dataclasses import dataclass
from typing import Iterable, Sequence, Tuple

import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(str(pathlib.Path(__file__).resolve().parents[2]))
    from comp_physics_python.ch3._linalg import generalized_eigh  # type: ignore
else:  # pragma: no cover
    from ..ch3._linalg import generalized_eigh


@dataclass
class H2Config:
    base_alphas: Tuple[float, ...] = (
        13.00773,
        1.962079,
        0.444529,
        0.1219492,
    )
    convergence: float = 1e-10
    max_iterations: int = 200

    def __post_init__(self) -> None:
        self.base_alphas = tuple(float(a) for a in self.base_alphas)
        self.n = len(self.base_alphas)
        self.total_basis = 2 * self.n
        self.full_alphas: Tuple[float, ...] = self.base_alphas + self.base_alphas
        self.pi = math.pi


def center_sign(index: int, half: int) -> int:
    return -1 if index < half else 1


def gaussian_overlap(a: float, b: float, same_center: bool, dist: float) -> float:
    factor = (math.pi / (a + b)) ** 1.5
    if same_center:
        return factor
    gamma = a * b / (a + b)
    return factor * math.exp(-gamma * dist * dist)


def build_overlap(cfg: H2Config, dist: float) -> np.ndarray:
    n = cfg.total_basis
    half = cfg.n
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            same = (i < half) == (j < half)
            S[i, j] = gaussian_overlap(cfg.full_alphas[i], cfg.full_alphas[j], same, dist)
    return S


def kinetic_same(a: float, b: float) -> float:
    K0 = a * b / (a + b)
    return 3.0 * K0 * (math.pi / (a + b)) ** 1.5


def kinetic_cross(a: float, b: float, overlap: float, dist: float) -> float:
    K0 = a * b / (a + b)
    return K0 * (3.0 - 2.0 * K0 * dist * dist) * overlap


def boys_function(t: float) -> float:
    if t <= 0.0:
        return 1.0
    return math.sqrt(math.pi) * math.erf(math.sqrt(t)) / (2.0 * math.sqrt(t))


def coulomb_same(a: float, b: float, dist: float) -> float:
    K0 = math.pi / (a + b)
    return -2.0 * K0 * (1.0 + boys_function((a + b) * dist * dist))


def coulomb_cross(a: float, b: float, overlap: float, dist: float) -> float:
    K0 = math.pi / (a + b)
    K1 = overlap / math.sqrt(K0)
    t1 = a * a * dist * dist / (a + b)
    t2 = b * b * dist * dist / (a + b)
    return -2.0 * K1 * (boys_function(t1) + boys_function(t2))


def build_hamilton(cfg: H2Config, dist: float, S: np.ndarray) -> np.ndarray:
    n = cfg.total_basis
    half = cfg.n
    H = np.zeros((n, n), dtype=float)
    for r in range(n):
        for s in range(n):
            a = cfg.full_alphas[r]
            b = cfg.full_alphas[s]
            same = (r < half) == (s < half)
            if same:
                val = kinetic_same(a, b) + coulomb_same(a, b, dist)
            else:
                overlap = S[r, s]
                val = kinetic_cross(a, b, overlap, dist) + coulomb_cross(a, b, overlap, dist)
            H[r, s] = val
    return H


def boys_argument(alpha_r: float, alpha_s: float, sign_r: int, sign_s: int, dist: float) -> float:
    return (alpha_r * sign_r + alpha_s * sign_s) * dist / (alpha_r + alpha_s) / 2.0


def two_electron_tensor(cfg: H2Config, dist: float, S: np.ndarray) -> np.ndarray:
    n = cfg.total_basis
    half = cfg.n
    q = np.zeros((n, n, n, n), dtype=float)
    for r in range(n):
        for s in range(n):
            for t in range(n):
                for u in range(n):
                    a = cfg.full_alphas[r]
                    b = cfg.full_alphas[s]
                    c = cfg.full_alphas[t]
                    d = cfg.full_alphas[u]
                    p1 = a + b
                    p2 = c + d
                    p3 = p1 + p2
                    P = (center_sign(r, half) * a + center_sign(s, half) * b) / p1 * dist / 2.0
                    Q = (center_sign(t, half) * c + center_sign(u, half) * d) / p2 * dist / 2.0
                    tt = p1 * p2 / p3 * (P - Q) ** 2
                    lam = 2.0 * math.sqrt(p1 * p2 / (math.pi * p3))
                    g = boys_function(tt)
                    q[r, s, t, u] = S[r, s] * S[t, u] * lam * g
    return q


def density_matrix(c: np.ndarray, S: np.ndarray) -> np.ndarray:
    norm = float(c @ (S @ c))
    return 2.0 * np.outer(c, c) / norm


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


def electronic_energy(dens: np.ndarray, H: np.ndarray, F: np.ndarray) -> float:
    n = H.shape[0]
    energy = 0.0
    for i in range(n):
        for j in range(i):
            energy += dens[i, j] * (H[i, j] + F[i, j])
        energy += 0.5 * dens[i, i] * (H[i, i] + F[i, i])
    return energy


def scf_h2(cfg: H2Config, dist: float) -> dict:
    S = build_overlap(cfg, dist)
    H = build_hamilton(cfg, dist, S)
    q = two_electron_tensor(cfg, dist, S)

    c = np.ones(cfg.total_basis, dtype=float)
    dens = density_matrix(c, S)

    energy = 0.0
    for _ in range(cfg.max_iterations):
        G = build_g_matrix(dens, q)
        F = H + G
        evals, evecs = generalized_eigh(F, S)
        c = evecs[:, 0]
        dens = density_matrix(c, S)
        new_energy = electronic_energy(dens, H, F)
        if abs(new_energy - energy) < cfg.convergence:
            energy = new_energy
            break
        energy = new_energy
    else:
        raise RuntimeError("SCF did not converge for H2")

    total_energy = energy + 1.0 / dist
    return {"distance": dist, "electronic_energy": energy, "total_energy": total_energy}


def sweep_h2(cfg: H2Config, distances: Sequence[float]) -> list[dict]:
    return [scf_h2(cfg, d) for d in distances]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=float, default=1.0)
    parser.add_argument("--stop", type=float, default=2.0)
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument("--out", type=str, default="e_vs_dist_py.dat")
    args = parser.parse_args(argv)

    cfg = H2Config()
    dists = np.linspace(args.start, args.stop, args.points)
    results = sweep_h2(cfg, dists)
    for row in results:
        print(f"R = {row['distance']:.3f}  E_total = {row['total_energy']:.8f}")
    data = np.array([[row["distance"], row["total_energy"]] for row in results])
    out = pathlib.Path(args.out)
    np.savetxt(out, data, fmt="%12.6f")
    print(f"Saved energies to {out}")


if __name__ == "__main__":
    main()
