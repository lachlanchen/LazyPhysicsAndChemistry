# Chapter 5 – Radial Hartree and DFT solvers (Python)

This folder mirrors the FORTRAN programs under
`comp_physics/comp_physics_textbook_code/4559_ch5/ch5`, providing concise Python
implementations with the same physical content:

- `hydrogen_radial.py` – solves the hydrogenic 1s problem (Sec. 5.3.1) by
  shooting the radial Schrödinger equation with the Numerov scheme.
- `helium_scf.py` – self‑consistent Hartree/DFT models for helium (Sec. 5.3.2–5.3.3),
  supporting Hartree, Slater Xα, exact (local) exchange, and LDA (exchange +
  Perdew–Zunger correlation) modes.
- `common.py` – shared radial-grid utilities: Numerov propagation, Poisson solver
  for the Hartree potential, local exchange/correlation parametrisations, etc.

## Radial grid and Numerov (common.py)

For a uniform grid `r_i = i h`, the radial Schrödinger equation for the reduced
wavefunction `u(r) = r R(r)` with `ℓ=0` reads

\[
\frac{\mathrm{d}^2 u}{\mathrm{d} r^2} = 2 (V(r) - E) u(r).
\]

`numerov_forward` integrates this ODE with the sixth‑order Numerov formula
starting from the regular boundary conditions `u(0)=0`, `u(h)=h`. The helper
`shoot_bound_state` brackets the bound eigenvalue between −3 Ha and 0 Ha, refines
it via bisection, and normalises `u(r)` with the grid integral
`∫₀^{r_max} u(r)^2 dr`.

The particle density for one electron is then

\[
\rho_1(r) = \frac{|u(r)|^2}{4\pi r^2}, \qquad \rho(r) = N_e \rho_1(r).
\]

## Hydrogen example (hydrogen_radial.py)

Usage:

```bash
python comp_physics_python/ch5/hydrogen_radial.py --step 0.02 --rmax 40
```

Typical output:

```
1s eigenvalue ≈ -0.50000016 Ha (exact -0.5)
```

The script also writes a table `r, ψ(r), ρ(r)` to `hydrogen_wavefunction.dat` for
plotting.

## Helium SCF model (helium_scf.py)

The general KS-like potential for helium is constructed as

\[
V(r) = -\frac{Z}{r} + V_H[\rho_1](r) - V_x(r) - V_c(r),
\]

where:

- `V_H` is obtained by solving the Poisson equation for `r V_H(r)` with the
  single-electron density (the other electron provides the screening charge).
- `V_x` is either zero (`hartree`), the Slater Xα potential (local-density
  exchange), an "exact" `ρ^{1/3}` form replicating `hydx.f`, or combined with
  the Perdew–Zunger correlation potential (`lda`).
- `V_c` is non-zero only for `lda`, using the Ceperley–Alder parameterisation
  (Eq. (5.65) in the book) for ε_c and V_c.

The self-consistent loop iterates `V`, computes the updated eigenvalue and
wavefunction via Numerov, refreshes densities and potentials, and stops when the
change in total energy

\[
E = N_e E_0 - E_H + E_x + E_c
\]

drops below the requested tolerance.  The script prints the converged orbital
energy and the total Hartree/DFT energy, and stores the radial potentials in an
`.npz` file for post-processing. Example:

```bash
python comp_physics_python/ch5/helium_scf.py slater --step 0.02 --rmax 25
```

produces

```
Converged slater model: ε0 = -0.91823845 Ha
Total energy E = -2.86186027 Ha
Saved Hartree/DFT potentials to helium_potentials.npz
```

These numbers reproduce the trends reported in Sec. 5.3 (Hartree energy above
−2.86 Ha, Slater and LDA bringing it closer to the experimental −2.903 Ha).
