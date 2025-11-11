# LazyPhysicsAndChemistry

LazyPhysicsAndChemistry is the code + notebook half of **LazyLearn**—my intentionally slow learning log for physics and chemistry. The living notes, wins, and TODOs surface at [learn.lazying.art](https://learn.lazying.art) (served directly from `docs/` in this repo), while the runnable artifacts stay here so experiments always have a home.

## LazyLearn

- **Home base:** [learn.lazying.art](https://learn.lazying.art) — the public-facing site with weekly focuses, backlog, and highlights.
- **Source of truth:** everything the site links to lives in `examples/`, `comp_physics/`, `multiwfn/`, or `figures/`.
- **Update flow:** ship code/notebooks first, regenerate plots if needed, and then add an entry to `docs/` so the site reflects the latest work.

## What lives here

- `examples/` – focused Python scripts (QAOA + VQE) that run on commodity laptops with Qiskit or PennyLane.
- `comp_physics/` – computational physics notebooks, helper scripts like `numerov.py`, and supporting data/figures that back the notebook chapters.
- `multiwfn/` – upstream Multiwfn 3.8 developer source drop plus the reference PDFs for quick lookup.
- `figures/` – static PNGs used by the scripts/notebooks for reports or slides.
- `docs/` – the LazyLearn microsite that GitHub Pages (or your own host) serves at `learn.lazying.art`.

Everything else (book scans, Gaussian symlinks, checkpoints, etc.) stays ignored so the repo stays small when cloned.

## Python environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

Jupyter notebooks inside `comp_physics/` expect the same environment; launch with `jupyter lab` or `jupyter notebook` after activating the venv.

## Example workflows

- **QAOA with Qiskit** – `python examples/qaoa_qiskit_maxcut.py` (no Aer dependency, pure statevector backend).
- **QAOA with PennyLane** – `python examples/qaoa_pennylane_maxcut.py` (uses `default.qubit`).
- **VQE for H₂** – `python examples/pennylane_chemistry_h2_vqe.py` to reproduce the `figures/pennylane_h2_vqe_convergence.png` curve.

All scripts log their intermediate metrics so you can reuse the resulting plots or extend them with new molecules/graphs.

## Computational physics notebooks

The `comp_physics/` directory mirrors the structure of the working notes:

- `comp_physics_textbook_code/` – reusable routines extracted from the notebooks.
- Standalone chapters such as `chapter1.ipynb`, `chapter2.ipynb`, `numerov.ipynb`, and `numpy_1ddft.ipynb`.
- Topic folders (`bosonscattering/`, `lensless/`, `lightscattering/`, etc.) that gather raw data and helper scripts per experiment.

Open any notebook after activating the Python environment. If you need extra dependencies, record them in `comp_physics/environments.yaml` before sharing.

## Multiwfn references

`multiwfn/` keeps the upstream `Multiwfn_3.8_dev_src_Linux` tree plus the PDF manual and quick-start guide so Gaussian/Multiwfn post-processing stays a short hop away. No binaries are committed; just point your Multiwfn build to that folder when needed.

## Figures

All generated PNGs live under `figures/`. If you regenerate or add plots from notebooks/scripts, drop them here so they stay versioned alongside the code that produced them.

## Reference stacks & related repos

LazyLearn leans on a couple of large companion repos that stay **outside** of version control here but are symlinked in for quick lookups:

- `leonardsusskind/` → symlink to `../leonardsusskind/`, which mirrors the lecture notes + supplemental packets from [github.com/lachlanchen/leonardsusskind](https://github.com/lachlanchen/leonardsusskind).
- `the_theoretical_minimum/` → symlink to `../the_theoretical_minimum/`, keeping the LaTeX sources, figures, and scripts from [github.com/lachlanchen/the_theoretical_minimum](https://github.com/lachlanchen/the_theoretical_minimum).

These repos are hefty (books, PDFs, figures), so they remain ignored via `.gitignore`. When setting up a fresh machine:

```bash
cd ~/ProjectsLFS
git clone git@github.com:lachlanchen/leonardsusskind.git
git clone git@github.com:lachlanchen/the_theoretical_minimum.git
cd QUANTUM
ln -s ../leonardsusskind/
ln -s ../the_theoretical_minimum/
```

Additional personal sandboxes such as `cellist` and `lazealoptix` live under `~/Projects/` and feed future LazyLearn chapters (optics + creative tooling), but they stay decoupled from this repo until specific experiments are ready.

## Version control notes

- Heavy directories such as `books/`, Gaussian symlinks, checkpoint files, and local scratch artifacts remain ignored via `.gitignore`.
- Keep contributions focused on the tracked folders above to preserve the "lazy" cloning experience.
- When updating the website, edit `docs/` locally, test via `python -m http.server --directory docs`, and push—GitHub Pages (or another host) can point the custom domain `learn.lazying.art` at this repo automatically.
