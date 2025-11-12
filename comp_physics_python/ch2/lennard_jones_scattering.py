"""Python translation of the Chapter 2 Lennard-Jones scattering code.

The original FORTRAN sources live under
``comp_physics/comp_physics_textbook_code/4556_ch2/ch2`` and compute the
elastic H–Kr cross section that Thijssen uses as his Section 2 example.

This module keeps the same numerical strategy (Numerov integration of the
radial Schrödinger equation, phase shifts from the logarithmic derivative, and
summation over partial waves) but wraps it in a modern Python API so it is easy
to experiment with alternative parameters or plug results into notebooks.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field
from typing import Iterable, List, Sequence

import numpy as np


@dataclass
class ScatteringConfig:
    """Container for the basic Lennard-Jones scattering parameters.

    Values mirror the prompts in ``scatter.f`` but come with sensible defaults
    so ``compute_cross_sections`` can run headless.
    """

    max_energy: float = 3.5
    min_energy: float = 0.1
    energy_points: int = 100  # EnerNum in the original program
    l_max: int = 6
    start_radius: float = 0.75  # starting point for Numerov integration
    step: float = 0.02  # integration step (HStep)
    max_radius: float = 5.0  # MaxDist for first integration window
    epsilon: float = 5.9  # Lennard-Jones well depth (in Rydberg units)
    rm: float = 3.57  # Lennard-Jones minimum position (Angstrom)
    rydberg_prefactor: float | None = None  # Rm^2 * 0.48 by default
    max_grid_size: int = 60_000

    # Derived quantities (filled in __post_init__)
    pi: float = field(init=False)
    ryd_const: float = field(init=False)
    max_steps_primary: int = field(init=False)

    def __post_init__(self) -> None:
        self.pi = math.pi
        if self.rydberg_prefactor is None:
            self.rydberg_prefactor = self.rm * self.rm * 0.48
        self.ryd_const = self.rydberg_prefactor
        # Snap max_radius to the Numerov grid, matching the FORTRAN behaviour.
        steps = max(1, int(round((self.max_radius - self.start_radius) / self.step)))
        self.max_steps_primary = steps
        self.max_radius = self.start_radius + steps * self.step

    @property
    def energies(self) -> np.ndarray:
        """Array of EnerNum+1 equidistant energies (same as FORTRAN loop)."""

        delta_e = (self.max_energy - self.min_energy) / self.energy_points
        values = self.min_energy + delta_e * np.arange(self.energy_points + 1)
        return values


# ---------------------------------------------------------------------------
# Bessel helpers
# ---------------------------------------------------------------------------

def spherical_bessel_j(l: int, x: float) -> float:
    """Spherical Bessel j_l(x) using the same upward recursion as ``special.f``."""

    if abs(x) < 1e-12:
        return 1.0 if l == 0 else 0.0
    if l == 0:
        return math.sin(x) / x
    if l == 1:
        return math.sin(x) / (x * x) - math.cos(x) / x
    jm2 = math.sin(x) / x
    jm1 = math.sin(x) / (x * x) - math.cos(x) / x
    j_val = jm1
    for ell in range(2, l + 1):
        j_val = (2 * ell - 1) / x * jm1 - jm2
        jm2, jm1 = jm1, j_val
    return j_val


def spherical_bessel_n(l: int, x: float) -> float:
    """Spherical Bessel n_l(x) using upward recursion (FORTRAN's ``SphBesN``)."""

    if abs(x) < 1e-12:
        # The Neumann function diverges at the origin; return a large value with
        # the correct sign to avoid numerical issues downstream.
        return -1e12
    if l == 0:
        return -math.cos(x) / x
    if l == 1:
        return -math.cos(x) / (x * x) - math.sin(x) / x
    nm2 = -math.cos(x) / x
    nm1 = -math.cos(x) / (x * x) - math.sin(x) / x
    n_val = nm1
    for ell in range(2, l + 1):
        n_val = (2 * ell - 1) / x * nm1 - nm2
        nm2, nm1 = nm1, n_val
    return n_val


# ---------------------------------------------------------------------------
# Numerov machinery
# ---------------------------------------------------------------------------

def lennard_jones_f(r: np.ndarray, l: int, energy: float, cfg: ScatteringConfig) -> np.ndarray:
    """Value of F(r) from the radial Schrödinger equation.

    With Thijssen's units the radial equation for u(r) = r R(r) reads

        u''(r) = [l(l+1)/r^2 + RydConst * epsilon (1/r^12 - 2/r^6)
                  - RydConst * energy] u(r)

    The Numerov method integrates this second-order ODE.  ``r`` can be either a
    float or an array; ``numpy`` broadcasting keeps the implementation compact.
    """

    r = np.asarray(r, dtype=float)
    r2 = r * r
    r4 = r2 * r2
    r6 = r4 * r2
    r12 = r6 * r6
    centrifugal = l * (l + 1) / r2
    lj = cfg.ryd_const * cfg.epsilon * (1.0 / r12 - 2.0 / r6)
    return centrifugal + lj - cfg.ryd_const * energy


def numerov(step: float, f_vals: np.ndarray, phi0: float, phi1: float) -> np.ndarray:
    """Numerov propagation matching the signature used in ``scatter.f``.

    ``f_vals`` must contain the value of ``F`` at every grid point
    (including the end point).  ``phi0`` and ``phi1`` act as the first two
    solution values at ``r = start`` and ``r = start + step``.
    """

    n = f_vals.size
    if n < 2:
        raise ValueError("Need at least two grid points for Numerov integration")

    solution = np.zeros(n, dtype=float)
    solution[0] = phi0
    solution[1] = phi1

    fac = (step * step) / 12.0
    w_prev = (1.0 - fac * f_vals[0]) * phi0
    w = (1.0 - fac * f_vals[1]) * phi1
    phi = phi1

    for idx in range(1, n - 1):
        w_next = 2.0 * w - w_prev + (step * step) * phi * f_vals[idx]
        w_prev, w = w, w_next
        denominator = 1.0 - fac * f_vals[idx + 1]
        if abs(denominator) < 1e-14:
            raise ZeroDivisionError("Numerov denominator vanished; adjust the grid")
        phi = w_next / denominator
        solution[idx + 1] = phi

    return solution


def _initial_conditions(cfg: ScatteringConfig, l: int, energy: float) -> tuple[float, float]:
    """Replicates the closed-form estimates for Phi(start) and Phi(start+h)."""

    exp_coeff = math.sqrt(cfg.epsilon * cfg.ryd_const / 25.0)
    phi_start = math.exp(-exp_coeff * cfg.start_radius ** -5)
    phi_deriv = 5.0 * exp_coeff * cfg.start_radius ** -6 * phi_start

    a = cfg.step ** 2 * float(lennard_jones_f(cfg.start_radius + cfg.step, l, energy, cfg)) / 12.0
    b = cfg.step ** 2 * float(lennard_jones_f(cfg.start_radius - cfg.step, l, energy, cfg)) / 12.0
    c = cfg.step ** 2 * float(lennard_jones_f(cfg.start_radius, l, energy, cfg))
    c = (2.0 + 5.0 * c / 6.0) * phi_start

    numerator = 2.0 * cfg.step * phi_deriv * (1.0 - b) + c * (1.0 - 2.0 * b)
    denominator = (1.0 - 2.0 * a) * (1.0 - b) + (1.0 - a) * (1.0 - 2.0 * b)
    phi_next = numerator / denominator
    return phi_start, phi_next


def _g_quotient(cfg: ScatteringConfig, l: int, energy: float, k: float) -> tuple[float, float]:
    """Integrate the radial equation and return (G quotient, second radius)."""

    phi_start, phi_next = _initial_conditions(cfg, l, energy)

    quarter_lambda = 0.5 * cfg.pi / k
    sec_r = cfg.max_radius + quarter_lambda
    sec_max_steps = max(1, int(round((sec_r - cfg.start_radius) / cfg.step)))
    if sec_max_steps + 1 > cfg.max_grid_size:
        raise ValueError(
            f"Grid too large: need {sec_max_steps + 1} points but max_grid_size is {cfg.max_grid_size}."
        )
    sec_r = cfg.start_radius + sec_max_steps * cfg.step

    grid = cfg.start_radius + cfg.step * np.arange(sec_max_steps + 1)
    f_vals = lennard_jones_f(grid, l, energy, cfg)
    solution = numerov(cfg.step, f_vals, phi_start, phi_next)

    phi_end_primary = solution[cfg.max_steps_primary]
    phi_end_secondary = solution[sec_max_steps]

    g_quot = (phi_end_secondary * cfg.max_radius) / (phi_end_primary * sec_r)
    return g_quot, sec_r


def phase_shift(cfg: ScatteringConfig, l: int, energy: float, k: float) -> float:
    """Compute the partial-wave phase shift δ_l(E)."""

    g_quot, sec_r = _g_quotient(cfg, l, energy, k)
    numerator = g_quot * spherical_bessel_j(l, k * cfg.max_radius) - spherical_bessel_j(l, k * sec_r)
    denominator = g_quot * spherical_bessel_n(l, k * cfg.max_radius) - spherical_bessel_n(l, k * sec_r)
    return math.atan2(numerator, denominator)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_cross_sections(cfg: ScatteringConfig, energies: Sequence[float] | None = None) -> list[dict]:
    """Return a list of dicts with energy, σ_tot, and the phase shifts."""

    energies = np.array(energies if energies is not None else cfg.energies, dtype=float)
    results: list[dict] = []

    for idx, energy in enumerate(energies, start=1):
        k = math.sqrt(cfg.ryd_const * energy)
        sigma_tot = 0.0
        shifts: List[float] = []
        for l in range(cfg.l_max + 1):
            delta = phase_shift(cfg, l, energy, k)
            shifts.append(delta)
            sigma_tot += 4.0 * cfg.pi / (k * k) * (2 * l + 1) * math.sin(delta) ** 2
        results.append(
            {
                "index": idx,
                "energy": energy,
                "k": k,
                "sigma_total": sigma_tot,
                "phase_shifts": shifts,
                "tan_delta_lmax": math.tan(shifts[-1]),
            }
        )
    return results


def write_sigmadat(results: Iterable[dict], path: str) -> None:
    """Mimic the FORTRAN output file (two-column Ener vs σ_tot)."""

    data = np.array([[row["energy"], row["sigma_total"]] for row in results], dtype=float)
    np.savetxt(path, data, fmt="%12.6f")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-energy", type=float, default=3.5)
    parser.add_argument("--min-energy", type=float, default=0.1)
    parser.add_argument("--energy-points", type=int, default=100)
    parser.add_argument("--l-max", type=int, default=6)
    parser.add_argument("--start", type=float, default=0.75, help="Starting radius for Numerov integration")
    parser.add_argument("--step", type=float, default=0.02, help="Integration step size")
    parser.add_argument("--max-radius", type=float, default=5.0, help="Upper radius for the first Numerov segment")
    parser.add_argument("--epsilon", type=float, default=5.9, help="Lennard-Jones depth ε")
    parser.add_argument("--rm", type=float, default=3.57, help="Lennard-Jones minimum Rm (Å)")
    parser.add_argument(
        "--out",
        type=str,
        default="sigmadat",
        help="Output file (Ener vs σ_tot). Defaults to the FORTRAN filename.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    cfg = ScatteringConfig(
        max_energy=args.max_energy,
        min_energy=args.min_energy,
        energy_points=args.energy_points,
        l_max=args.l_max,
        start_radius=args.start,
        step=args.step,
        max_radius=args.max_radius,
        epsilon=args.epsilon,
        rm=args.rm,
    )

    results = compute_cross_sections(cfg)
    write_sigmadat(results, args.out)

    for row in results:
        print(
            f"{row['index']:4d}  Ener = {row['energy']:8.4f}  "
            f"SigmaTot = {row['sigma_total']:10.6f}  TanDelta(l_max) = {row['tan_delta_lmax']:9.5f}"
        )
    print(f"Saved Ener vs σ_tot to {args.out}")


if __name__ == "__main__":
    main()
