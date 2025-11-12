"""Shared radial-grid utilities for Chapter 5 Hartree/DFT models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass
class RadialGrid:
    """Uniform radial grid starting at 0."""

    step: float
    max_radius: float

    def __post_init__(self) -> None:
        n = int(round(self.max_radius / self.step))
        self.r = np.linspace(0.0, n * self.step, n + 1)
        self.step = self.r[1] - self.r[0]
        self.surface = 4.0 * np.pi * self.r * self.r
        self._safe_r = np.where(self.r > 1e-12, self.r, 1e-12)

    def integrate(self, values: np.ndarray) -> float:
        return np.trapz(values, self.r)

    @property
    def safe_r(self) -> np.ndarray:
        return self._safe_r


def numerov_forward(f_vals: np.ndarray, grid: RadialGrid) -> np.ndarray:
    """Integrate u'' = f(r) u using the Numerov method."""

    h = grid.step
    n = f_vals.size
    u = np.zeros(n)
    u[1] = h  # u ~ r near the origin for l=0
    h2 = h * h
    for i in range(1, n - 1):
        k_minus = 1.0 + h2 * f_vals[i - 1] / 12.0
        k = 1.0 + h2 * f_vals[i] / 12.0
        k_plus = 1.0 + h2 * f_vals[i + 1] / 12.0
        u[i + 1] = (
            (2.0 * u[i] * (1.0 - 5.0 * h2 * f_vals[i] / 12.0) - u[i - 1] * k_minus)
            / k_plus
        )
    return u


def numerov_inward(f_vals: np.ndarray, grid: RadialGrid, z: float) -> np.ndarray:
    """Integrate inward from r_max using the hydrogenic tail ~ r e^{-Z r}."""

    h = grid.step
    n = f_vals.size
    u = np.zeros(n)
    r = grid.safe_r
    u[-1] = r[-1] * np.exp(-z * r[-1])
    u[-2] = r[-2] * np.exp(-z * r[-2])
    h2 = h * h
    for i in range(n - 2, 0, -1):
        k_plus = 1.0 + h2 * f_vals[i + 1] / 12.0
        k = 1.0 + h2 * f_vals[i] / 12.0
        k_minus = 1.0 + h2 * f_vals[i - 1] / 12.0
        u[i - 1] = (
            2.0 * u[i] * (1.0 - 5.0 * h2 * f_vals[i] / 12.0) - k_plus * u[i + 1]
        ) / k_minus
    u[0] = 2.0 * u[1] - u[2] + h2 * f_vals[1] * u[1]
    return u


def shoot_bound_state(
    potential: np.ndarray,
    grid: RadialGrid,
    energy_bounds: Tuple[float, float] = (-3.0, -0.01),
    tol: float = 1e-10,
    tail_z: float | None = None,
) -> Tuple[float, np.ndarray]:
    """Return eigenvalue and radial solution for u'' = 2(V - E)u."""

    def residual(energy: float) -> Tuple[float, np.ndarray]:
        f = 2.0 * (potential - energy)
        if tail_z is None:
            u = numerov_forward(f, grid)
            slope = (u[-1] - u[-3]) / (2.0 * grid.step)
            kappa = np.sqrt(max(-2.0 * energy, 1e-12))
            return slope / u[-1] + kappa, u
        u = numerov_inward(f, grid, tail_z)
        return u[0], u

    e_low, e_high = energy_bounds
    res_low, _ = residual(e_low)
    res_high, _ = residual(e_high)
    if res_low * res_high > 0:
        scan = np.linspace(e_low, e_high, 600)
        prev_e = scan[0]
        prev_res, _ = residual(prev_e)
        for e in scan[1:]:
            curr_res, _ = residual(e)
            if prev_res * curr_res <= 0:
                e_low, e_high = prev_e, e
                res_low, res_high = prev_res, curr_res
                break
            prev_res, prev_e = curr_res, e
        else:
            raise RuntimeError("Failed to bracket eigenvalue")

    if e_low > e_high:
        e_low, e_high = e_high, e_low
        res_low, res_high = res_high, res_low

    for _ in range(160):
        mid = 0.5 * (e_low + e_high)
        res_mid, u_mid = residual(mid)
        if abs(res_mid) < tol or abs(e_high - e_low) < tol:
            norm = np.sqrt(grid.integrate(u_mid**2))
            return mid, u_mid / norm
        if res_low * res_mid < 0:
            e_high, res_high = mid, res_mid
        else:
            e_low, res_low = mid, res_mid
    raise RuntimeError("Numerov shoot failed to converge")


def densities_from_u(u: np.ndarray, grid: RadialGrid, electrons: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """Return single-electron and total densities."""

    density_single = np.zeros_like(u)
    density_single[1:] = (u[1:] / grid.safe_r[1:]) ** 2 / (4.0 * np.pi)
    density_single[0] = density_single[1]
    norm = grid.integrate(density_single * grid.surface)
    if norm > 0:
        density_single /= norm
    density_total = electrons * density_single
    return density_single, density_total


def hartree_potential(density_single: np.ndarray, grid: RadialGrid, target_charge: float = 1.0) -> np.ndarray:
    """Solve Poisson's equation for the Hartree potential of one electron."""

    rhs = -4.0 * np.pi * density_single * grid.safe_r
    h = grid.step
    phi = np.zeros_like(density_single)
    phi[1] = grid.safe_r[1]
    for i in range(1, len(phi) - 1):
        phi[i + 1] = 2.0 * phi[i] - phi[i - 1] + h * h * rhs[i]
    correction = target_charge * grid.safe_r[-1] - phi[-1]
    phi += correction * (grid.safe_r / grid.safe_r[-1])
    V = np.zeros_like(phi)
    V[1:] = phi[1:] / grid.safe_r[1:]
    V[0] = V[1]
    return V


def lda_exchange(density_total: np.ndarray, grid: RadialGrid) -> Tuple[np.ndarray, float]:
    """Local density approximation (Slater) exchange potential and energy."""

    pref = (3.0 / np.pi) ** (1.0 / 3.0)
    n = np.clip(density_total, 0.0, None)
    eps_x = -0.75 * pref * np.power(n, 1.0 / 3.0)
    vx = (4.0 / 3.0) * eps_x
    energy = grid.integrate(eps_x * n * grid.surface)
    vx[0] = vx[1]
    return vx, energy


def exact_exchange_from_density(density_total: np.ndarray, grid: RadialGrid) -> Tuple[np.ndarray, float]:
    coeff = (3.0 / (2.0 * np.pi * np.pi)) ** (1.0 / 3.0)
    n = np.clip(density_total, 0.0, None)
    vx = -coeff * np.power(n, 1.0 / 3.0)
    vx[0] = vx[1]
    energy = -0.5 * grid.integrate(vx * n * grid.surface)
    return vx, energy


def ceperley_alder(density: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (epsilon_c, V_c) using the Perdew–Zunger fit."""

    gamma = -0.1423
    beta1 = 1.0529
    beta2 = 0.3334
    a = 0.0311
    b = -0.048
    c = 0.0020
    d = -0.0116

    rs = np.where(density > 0, (0.75 / (np.pi * density)) ** (1.0 / 3.0), np.inf)
    eps = np.zeros_like(rs)
    vc = np.zeros_like(rs)

    mask = rs >= 1.0
    sqrt_rs = np.sqrt(rs[mask])
    eps[mask] = gamma / (1.0 + beta1 * sqrt_rs + beta2 * rs[mask])
    vc[mask] = eps[mask] * (
        1.0 + 7.0 * beta1 * sqrt_rs / 6.0 + 4.0 * beta2 * rs[mask] / 3.0
    ) / (1.0 + beta1 * sqrt_rs + beta2 * rs[mask])

    mask_inv = ~mask
    rs_inv = rs[mask_inv]
    eps[mask_inv] = a * np.log(rs_inv) + b + c * rs_inv * np.log(rs_inv) + d * rs_inv
    vc[mask_inv] = (
        a * np.log(rs_inv)
        + (b - a / 3.0)
        + 2.0 * c * rs_inv * np.log(rs_inv) / 3.0
        + (2.0 * d - c) * rs_inv / 3.0
    )

    eps[np.isinf(rs)] = 0.0
    vc[np.isinf(rs)] = 0.0

    return eps, vc
