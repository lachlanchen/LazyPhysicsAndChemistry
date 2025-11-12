# Chapter 10 – Monte Carlo Methods in Python

This folder mirrors the Chapter 10 examples from Thijssen’s *Computational Physics*.
Each script keeps the physics and algorithms of the original Fortran versions
but wraps them in short, runnable Python programs.

## 1. Lennard‑Jones Metropolis MC (`lj_mc.py`)

We place `N = 4L³` atoms on an FCC lattice of box length
\(L = (N/\rho)^{1/3}\) and sample the canonical ensemble with the Metropolis
algorithm.  The energy change for moving particle *i* is computed with the
minimum-image convention and the truncated LJ potential

\[ V(r) = 4ε\left[ (σ/r)^{12} − (σ/r)^6 \right], \quad r < r_c. \]

Acceptance probability: `min(1, exp(-β ΔE))`.  Every 100 trial moves we adapt the
maximal displacement so the acceptance hovers around 40 %, reproducing the
behaviour of `mc.F`.

Example:

```bash
python comp_physics_python/ch10/lj_mc.py --n 256 --density 0.85 --temp 1.2 \
       --steps 40000 --burn 5000
```

Outputs mean energy per particle and heat capacity, plus the final configuration
(`lj_mc_positions.npy`).

## 2. 2D Ising Metropolis (`ising_mc.py`)

Implements the lattice update used in `IsMC.F`.  For coupling
`J = K/kT` and field `H = B/kT`, the flip cost is
`ΔE = 2 s_i (J Σ_j s_j + H)`.  Metropolis sampling accumulates \(E/N\), specific
heat, and magnetisation after a burn-in of 2000 sweeps.

```
python comp_physics_python/ch10/ising_mc.py --size 32 --J 0.44 --H 0.0 --steps 6000
```

## 3. Rosenbluth polymer growth (`rosenbluth_polymer.py`)

Translates `rosenbluth.f90`.  A chain is grown segment by segment; for each
monomer we evaluate `θ` trial directions, weight them by
\(\exp(-β V_{\text{LJ}})\), pick one with probability proportional to the weight, and
multiply the chain weight by the sum of the trial weights.  The script then
estimates the average end-to-end distance squared:

\[
\langle R^2 \rangle = \frac{\sum w_i R_i^2}{\sum w_i}.
\]

```
python comp_physics_python/ch10/rosenbluth_polymer.py --length 73 --steps 30000
```
