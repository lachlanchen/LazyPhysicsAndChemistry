"""Rosenbluth polymer-growth algorithm (Thijssen §10.6.1)."""

from __future__ import annotations

import argparse
import numpy as np


def lennard_jones_tail(dist2: float, sigma6: float, sigma12: float, epsilon: float) -> float:
    inv_r2 = 1.0 / dist2
    inv_r6 = inv_r2**3
    return 4.0 * epsilon * (-sigma6 * inv_r6 + sigma12 * inv_r6**2)


def potential(new_pos: np.ndarray, polymer: np.ndarray, sigma6: float, sigma12: float, epsilon: float) -> float:
    rel = polymer - new_pos
    dist2 = np.sum(rel * rel, axis=1)
    mask = dist2 > 0.0
    if not np.any(mask):
        return 0.0
    return np.sum([lennard_jones_tail(d, sigma6, sigma12, epsilon) for d in dist2[mask]])


def grow_polymer(length: int, theta_num: int, sigma: float, epsilon: float, rng: np.random.Generator) -> tuple[np.ndarray, float]:
    polymer = np.zeros((length, 2))
    polymer[1] = np.array([1.0, 0.0])
    sigma6 = sigma**6
    sigma12 = sigma**12
    weight = 1.0
    for j in range(1, length - 1):
        theta0 = rng.random() * 2 * np.pi
        weights = []
        candidates = []
        for k in range(theta_num):
            theta = theta0 + 2 * np.pi * k / theta_num
            new_pos = polymer[j] + np.array([np.cos(theta), np.sin(theta)])
            w = np.exp(-potential(new_pos, polymer[: j + 1], sigma6, sigma12, epsilon))
            weights.append(w)
            candidates.append(new_pos)
        weights = np.array(weights)
        total = weights.sum()
        if total == 0:
            return grow_polymer(length, theta_num, sigma, epsilon, rng)
        weights /= total
        idx = rng.choice(theta_num, p=weights)
        polymer[j + 1] = candidates[idx]
        weight *= total / theta_num
    return polymer, weight


def run_rosenbluth(length: int, steps: int, theta_num: int, sigma: float, epsilon: float, seed: int | None = None) -> float:
    rng = np.random.default_rng(seed)
    total_weight = 0.0
    weighted_r2 = 0.0
    for _ in range(steps):
        polymer, weight = grow_polymer(length, theta_num, sigma, epsilon, rng)
        end_to_end = np.sum((polymer[-1] - polymer[0]) ** 2)
        weighted_r2 += weight * end_to_end
        total_weight += weight
    return weighted_r2 / total_weight


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=73)
    parser.add_argument("--steps", type=int, default=30000)
    parser.add_argument("--theta", type=int, default=6)
    parser.add_argument("--sigma", type=float, default=0.8)
    parser.add_argument("--epsilon", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=1234)
    args = parser.parse_args(argv)
    avg = run_rosenbluth(args.length, args.steps, args.theta, args.sigma, args.epsilon, seed=args.seed)
    print(f"<R^2> = {avg:.4f}")


if __name__ == "__main__":
    main()
