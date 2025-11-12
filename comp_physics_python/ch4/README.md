# Chapter 4 – Two-electron Hartree–Fock in Python

This folder mirrors the FORTRAN codes in
`comp_physics/comp_physics_textbook_code/ch4`, providing well-documented Python
implementations and a short recap of the underlying equations.

## Helium atom (`helium_hartree_fock.py`)

Thijssen’s helium program builds a four-function Gaussian basis
(Section 4.3.2) and solves the restricted Hartree–Fock equations
self-consistently.  Each primitive is of the form

\[ \phi_r(\mathbf{r}) = N_r e^{-\alpha_r r^2}. \]

The key matrices are

- overlap: \(S_{rs} = (\pi/(\alpha_r + \alpha_s))^{3/2}\),
- kinetic: \(T_{rs} = 3\,\alpha_r \alpha_s (\pi/(\alpha_r + \alpha_s))^{3/2}/(\alpha_r + \alpha_s)\),
- nuclear attraction: \(V_{rs} = -2 Z \pi / (\alpha_r + \alpha_s)\),
- electron repulsion (all orbitals share one centre):
  \(g_{rstu} = 2\pi^2\sqrt{\pi} / [\sqrt{(\alpha_r+\alpha_s+\alpha_t+\alpha_u)} (\alpha_r+\alpha_s)(\alpha_t+\alpha_u)]\).

With the closed-shell density \(D = 2 c c^\mathrm{T} / (c^\mathrm{T} S c)\), the Fock matrix is
\(F = H + G\) where
\(G_{rs} = \sum_{tu} D_{tu} g_{rstu}\) (the code keeps the FORTRAN convention
with the 0.5 factors for \(t=u\)).  The script iterates until the total energy

\[ E = \varepsilon_0 + \sum_{rs}' D_{rs} H_{rs} \]

stabilises (`'` denotes the symmetric sum used in the original source).  The
CLI mimics `Hatom.f` and writes the radial ground-state wave function to
`WaveFunc_he_py`.

## H₂ molecule (`h2_hartree_fock.py`)

The H₂ code doubles the Gaussian set: four primitives on each proton separated
by \(R\).  Overlap, kinetic, and Coulomb integrals distinguish between same‑ and
cross‑centre pairs, e.g.

\[
S_{rs} = \left(\frac{\pi}{\alpha_r + \alpha_s}\right)^{3/2}
          \exp\!\left[-\frac{\alpha_r \alpha_s}{\alpha_r+\alpha_s} R^2\right]
\]

for functions located on different nuclei.  The kinetic and Coulomb elements
use the same algebra as `h2.f`, relying on the Boys function
\(F(t) = \sqrt{\pi}\,\mathrm{erf}(\sqrt{t}) / (2\sqrt{t})\).

Two-electron integrals follow Eq. (4.53) of the book; the implementation keeps
Thijssen’s notation with the \(P\) and \(Q\) centres and leverages the already
computed overlap entries:

\[
\langle rs | tu \rangle = S_{rs} S_{tu} \cdot 2 \sqrt{\frac{(\alpha_r+\alpha_s)(\alpha_t+\alpha_u)}{\pi\,P_3}}
F\!\Big(\frac{(\alpha_r+\alpha_s)(\alpha_t+\alpha_u)}{P_3} (P-Q)^2\Big).
\]

Self-consistency loops over the Fock matrix as in the helium code; after
convergence the nuclear repulsion \(1/R\) is added to the electronic energy.
Running the CLI reproduces the `e_vs_dist.dat` scan:

```bash
python comp_physics_python/ch4/h2_hartree_fock.py --start 1.0 --stop 2.0 --points 100 --out e_vs_dist_py.dat
```

The output matches the FORTRAN sensitivities to within machine precision, and
the saved energy curve can be compared directly with the original `e_vs_dist`.
