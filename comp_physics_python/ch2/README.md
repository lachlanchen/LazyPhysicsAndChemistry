# Chapter 2 – Lennard‑Jones Scattering in Python

This folder mirrors the FORTRAN sources in
`comp_physics/comp_physics_textbook_code/4556_ch2/ch2`.  The original
implementation from Thijssen’s *Computational Physics* (Section 2) evaluates the
H–Kr elastic cross section by shooting the radial Schrödinger equation with the
Numerov method and extracting phase shifts from the asymptotic solution.  The
Python rewrite keeps the exact same physics while adding documentation and a
clean API.

## Mathematical overview

We study a single channel scattering problem with an effective Lennard‑Jones
interaction

\[
V(r) = \varepsilon \left[\left(\frac{R_m}{r}\right)^{12} - 2\left(\frac{R_m}{r}\right)^6\right],
\]

expressed in Thijssen’s scaled units where the reduced mass and \(\hbar\) factors
collapse into a single prefactor `RydConst = R_m^2 · 0.48`.  Writing the radial
wave function as \(u(r) = r R(r)\), the Schrödinger equation becomes

\[
\frac{\mathrm{d}^2 u}{\mathrm{d} r^2} = \Bigg[\frac{\ell (\ell + 1)}{r^2}
+ \texttt{RydConst}·\varepsilon\Big(\frac{1}{r^{12}} - \frac{2}{r^{6}}\Big)
- \texttt{RydConst}·E \Bigg] u(r) 
\equiv F(r;\ell,E) u(r).
\]

The Numerov method integrates this second‑order ODE with local error
\(\mathcal{O}(h^6)\).  Following the book we start the integration at
`start_radius` with the same analytic Taylor expansion that appears in
`scatter.f`, so that the first two points \(u(r_0)\) and \(u(r_1)\) are
consistent with the rapidly varying repulsive core.  Two integrations are
performed:

1. from `start_radius` up to `max_radius` (continuous version of the user’s
   `MaxDist`) so we can read \(u(r_{\text{max}})\);
2. from the same starting point to a second radius
   \(r_\mathrm{sec} = r_{\text{max}} + \lambda/4\) where \(\lambda = 2\pi/k\).

The ratio of the two end‑points eliminates the unknown normalization and gives
Thijssen’s

\[
G = \frac{u(r_\mathrm{sec}) · r_{\text{max}}}{u(r_{\text{max}}) · r_\mathrm{sec}}.
\]

Matching this numerical solution to the known spherical Bessel functions outside
the potential then yields the phase shift

\[
\tan\delta_\ell = \frac{G j_\ell(kr_{\text{max}}) - j_\ell(kr_\mathrm{sec})}
{G n_\ell(kr_{\text{max}}) - n_\ell(kr_\mathrm{sec})}.
\]

Finally the elastic cross section follows from the standard partial‑wave sum

\[
\sigma(E) = \frac{4\pi}{k^2} \sum_{\ell=0}^{\ell_{\max}}
(2\ell+1) \sin^2 \delta_\ell(E), \qquad k = \sqrt{\texttt{RydConst} · E}.
\]

## Files

- `lennard_jones_scattering.py` — fully documented translation that exposes
  `ScatteringConfig`, `compute_cross_sections`, and a CLI compatible with the
  FORTRAN prompts.
- `README.md` (this file) — derivation and usage instructions.

## Usage

```bash
# Run with the default parameters from the book and reproduce sigmadat
python comp_physics_python/ch2/lennard_jones_scattering.py --out sigmadat_py

# Play with the potential depth and angular momentum cutoff
python comp_physics_python/ch2/lennard_jones_scattering.py \
    --epsilon 6.5 --l-max 8 --energy-points 150 --out sigmadat_eps6p5
```

Inside a notebook you can work with the returned dictionaries directly:

```python
from comp_physics_python.ch2.lennard_jones_scattering import (
    ScatteringConfig, compute_cross_sections,
)

cfg = ScatteringConfig(max_energy=4.0, energy_points=200, l_max=8)
results = compute_cross_sections(cfg)
energies = [row["energy"] for row in results]
sigma = [row["sigma_total"] for row in results]
```

Each `row["phase_shifts"]` entry stores the full list `[δ_0, …, δ_{ℓmax}]` so you
can analyse individual partial waves or verify Levinson’s theorem.

## Notes on numerical stability

- The Numerov integration uses the same step size for both radii.  Make sure the
  requested `max_radius` keeps `sec_r` within `max_grid_size` = 60 000 grid
  points.
- If the denominator inside the Numerov update or the phase‑shift expression
  vanishes, decrease `step` or increase `l_max` gradually.  The code raises a
  descriptive exception instead of failing silently.
- The spherical Bessel functions follow the upward recursion from `special.f`
  so the Python and FORTRAN numbers agree to floating‑point precision for the
  energy range used in the book.
