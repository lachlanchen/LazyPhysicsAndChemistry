# LazyPhysics and Chemistry

LazyPhysics and Chemistry is the code + notebook half of **LazyLearn**—my intentionally slow learning log for physics and chemistry. The living notes, wins, and TODOs surface at [learn.lazying.art](https://learn.lazying.art) (served directly from `docs/` in this repo), while the runnable artifacts stay here so experiments always have a home.

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

## Support LazyLearn

Helping LazyLearn keeps the experiments, documentation, and open tooling flowing:

- Cover hosting/inference/storage for the public demos and notebooks.
- Fund focused hack-weeks on EchoMind, LazyEdit, and the quantum/physics utilities here.
- Prototype optics + wearables (IdeasGlass, LightMind) that feed future chapters.
- Sponsor free deployments for students, community labs, and creators.

### Donate

<div align="center">
<table style="margin:0 auto; text-align:center; border-collapse:collapse;">
  <tr>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;">
      <a href="https://chat.lazying.art/donate">https://chat.lazying.art/donate</a>
    </td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;">
      <a href="https://chat.lazying.art/donate"><img src="figures/donate_button.svg" alt="Donate" height="44"></a>
    </td>
  </tr>
  <tr>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;">
      <a href="https://paypal.me/RongzhouChen">
        <img src="https://img.shields.io/badge/PayPal-Donate-003087?logo=paypal&logoColor=white" alt="Donate with PayPal">
      </a>
    </td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;">
      <a href="https://buy.stripe.com/aFadR8gIaflgfQV6T4fw400">
        <img src="https://img.shields.io/badge/Stripe-Donate-635bff?logo=stripe&logoColor=white" alt="Donate with Stripe">
      </a>
    </td>
  </tr>
  <tr>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><strong>WeChat</strong></td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><strong>Alipay</strong></td>
  </tr>
  <tr>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><img alt="WeChat QR" src="figures/donate_wechat.png" width="240"/></td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><img alt="Alipay QR" src="figures/donate_alipay.png" width="240"/></td>
  </tr>
</table>
</div>

**支援 / Donate**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります。  
- 你的支持将用于研发与运维，帮助我持续公开分享更多项目与改进。  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## Version control notes

- Heavy directories such as `books/`, Gaussian symlinks, checkpoint files, and local scratch artifacts remain ignored via `.gitignore`.
- Keep contributions focused on the tracked folders above to preserve the "lazy" cloning experience.
- When updating the website, edit `docs/` locally, test via `python -m http.server --directory docs`, and push—GitHub Pages (or another host) can point the custom domain `learn.lazying.art` at this repo automatically.
