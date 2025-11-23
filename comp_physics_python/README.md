# comp_physics_python

Python rewrites of the Fortran textbook programs that live under
`comp_physics/comp_physics_textbook_code`. Every folder mirrors a chapter in
Thijssen’s *Computational Physics* and keeps the same numerical ideas while
modernising the implementation (NumPy, concise CLIs, docstrings, and portable
plot-friendly outputs).

## Environment

All scripts are tested in the `quantum` conda environment mentioned in the root
README. Typical setup:

```bash
conda activate quantum  # or python -m venv .venv
pip install numpy scipy matplotlib
```

Chapter-specific scripts may need extras (e.g. JAX for large-scale optimisation),
but everything implemented so far (Ch. 2–5, 8, 10) runs with plain NumPy/SciPy.

## Folders

| Folder | Chapter | Highlights |
| --- | --- | --- |
| [`ch2/`](ch2/README.md) | §2 – Lennard-Jones scattering | `lennard_jones_scattering.py` reproduces the APW-style phase-shift code with Numerov integration and a CLI for generating σ(E) tables. |
| [`ch3/`](ch3/README.md) | §3 – Variational calculations | `deep_well_variational.py` and `hydrogen_gaussians.py` plus a tiny `generalized_eigh` helper to solve the generalised eigenvalue problems described in the text. |
| [`ch4/`](ch4/README.md) | §4 – Restricted Hartree–Fock | Helium and H₂ Hartree–Fock solvers (`helium_hartree_fock.py`, `h2_hartree_fock.py`) that reproduce the basis, potentials, and energy bookkeeping of the original Fortran. |
| [`ch5/`](ch5/README.md) | §5 – Radial Hartree/DFT | `hydrogen_radial.py` (Numerov + shooting) and `helium_scf.py` (Hartree / Slater / LDA with Ceperley–Alder) share utilities in `common.py`. |
| [`ch6/`](ch6/README.md) | §6 – APW & pseudopotentials | `apw/` ports `logapw.f` (Γ→K APW spectra) and `pseudo/` reproduces the semi-empirical silicon pseudopotential with Löwdin partitioning. |
| [`ch8/`](ch8/README.md) | §8 – Molecular dynamics | `ar_md.py` (LJ MD) and `n2_md.py` (rigid diatomic MD) with reusable lattice/Verlet helpers. |
| [`ch10/`](ch10/README.md) | §10 – Monte Carlo samplers | `lj_mc.py`, `ising_mc.py`, and `rosenbluth_polymer.py` cover the Lennard-Jones Metropolis, 2D Ising Metropolis, and Rosenbluth polymer-growth algorithms. |

## Usage pattern

Each script exposes a CLI matching the parameters used in the book:

```bash
# Lennard-Jones MD (Chapter 8)
python comp_physics_python/ch8/ar_md.py --n 256 --density 0.8442 --temp 0.722 --dt 0.004

# Helium Hartree-Fock (Chapter 4)
python comp_physics_python/ch4/helium_hartree_fock.py --max-iter 200 --tol 1e-11
```

Outputs are kept simple (stdout summaries + optional `.npz`/`.npy` dumps) so they
can feed back into the LazyLearn notebooks or plotting scripts.

## Roadmap

- Add remaining chapters (transfer matrices, DMC/PIMC, FEM) to the Python tree.
- Harmonise plotting/output conventions and ship lightweight tests for every
script.
- Keep the chapter mapping in this README up to date as new ports land.
