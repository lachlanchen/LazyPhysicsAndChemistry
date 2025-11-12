<p>
  <b>Langues :</b>
  <a href="../README.md">English</a>
  · <a href="README.zh-Hans.md">中文 (简体)</a>
  · <a href="README.zh-Hant.md">中文（繁體）</a>
  · <a href="README.ja.md">日本語</a>
  · <a href="README.ko.md">한국어</a>
  · <a href="README.vi.md">Tiếng Việt</a>
  · <a href="README.ar.md">العربية</a>
  · <a href="README.fr.md">Français</a>
  · <a href="README.es.md">Español</a>
</p>

# LazyPhysics and Chemistry

LazyPhysics and Chemistry constitue la partie code + notebooks de **LazyLearn**, mon journal d’apprentissage volontairement lent en physique et chimie. Les notes publiques, objectifs et réussites sont visibles sur [learn.lazying.art](https://learn.lazying.art) (généré depuis `docs/`), tandis que les artefacts exécutables vivent ici pour rester reproductibles.

## LazyLearn

- **Base principale :** [learn.lazying.art](https://learn.lazying.art) — focus hebdomadaires, backlog et highlights.
- **Source de vérité :** tout lien du site pointe vers `examples/`, `comp_physics/`, `multiwfn/` ou `figures/`.
- **Flux de mise à jour :** publier d’abord le code/notebook, régénérer les graphiques si nécessaire, puis ajouter une entrée dans `docs/` pour synchroniser le site.

## Contenu du dépôt

- `examples/` – scripts Python QAOA/VQE prêts à tourner sur un laptop avec Qiskit ou PennyLane.
- `comp_physics/` – notebooks de physique numérique, scripts utilitaires (`numerov.py`) et données/figures associées.
- `comp_physics_python/` – portage Python des programmes Fortran du manuel *Computational Physics*, classé par chapitre (voir [comp_physics_python/README.md](../comp_physics_python/README.md)).
- `multiwfn/` – dépôt source Multiwfn 3.8 dev et documentation PDF.
- `figures/` – PNG statiques issus des notebooks/scripts pour rapports ou présentations.
- `docs/` – microsite LazyLearn publié sur `learn.lazying.art` via GitHub Pages ou hébergeur perso.

Le reste (scans de livres, symlinks Gaussian, checkpoints, etc.) est ignoré via `.gitignore` afin de garder le dépôt léger.

## Environnement Python

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

Activez l’environnement puis lancez `jupyter lab` ou `jupyter notebook` pour travailler dans `comp_physics/`.

## Workflows d’exemple

- **QAOA (Qiskit)** – `python examples/qaoa_qiskit_maxcut.py` (statevector uniquement).
- **QAOA (PennyLane)** – `python examples/qaoa_pennylane_maxcut.py` (backend `default.qubit`).
- **VQE pour H₂** – `python examples/pennylane_chemistry_h2_vqe.py` pour reproduire `figures/pennylane_h2_vqe_convergence.png`.

Chaque script journalise ses métriques intermédiaires afin de réutiliser ou d’étendre facilement les expériences.

## Notebooks de physique numérique

`comp_physics/` reflète la structure des notes papier :

- `comp_physics_textbook_code/` – routines réutilisables extraites des chapitres.
- Notebooks indépendants : `chapter1.ipynb`, `chapter2.ipynb`, `numerov.ipynb`, `numpy_1ddft.ipynb`, etc.
- Dossiers thématiques (`bosonscattering/`, `lensless/`, `lightscattering/`…) regroupant données brutes et scripts.

Documentez toute dépendance supplémentaire dans `comp_physics/environments.yaml`.

## Références Multiwfn

Le dossier `multiwfn/` contient `Multiwfn_3.8_dev_src_Linux`, le manuel et le guide rapide. Aucun binaire n’est fourni : compilez directement depuis cette source.

## Figures

Toutes les images PNG sont stockées dans `figures/`. Ajoutez-y vos nouveaux graphiques pour qu’ils soient versionnés.

## Portage des codes du manuel

`comp_physics_python/` regroupe les versions Python des programmes Fortran décrits dans *Computational Physics*. Chaque sous-dossier correspond à un chapitre (ex. `ch4/` = Hartree–Fock, `ch8/` = dynamique moléculaire, `ch10/` = Monte Carlo) et chaque script propose une CLI pour relancer les expériences dans un environnement Python/NumPy moderne. Voir [comp_physics_python/README.md](../comp_physics_python/README.md) pour la cartographie détaillée et les instructions d’exécution.

## Soutenir LazyLearn

Soutenir le projet permet de :

- Couvrir l’hébergement, l’inférence et le stockage des démos/publiques.
- Financer des « sprints » open-source sur EchoMind, LazyEdit et les utilitaires quantique/physique.
- Prototyper optique et wearables (IdeasGlass, LightMind) pour les prochains chapitres.
- Offrir des déploiements gratuits aux étudiants, labs communautaires et créateurs.

### Faire un don

<div align="center">
<table style="margin:0 auto; text-align:center; border-collapse:collapse;">
  <tr>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;">
      <a href="https://chat.lazying.art/donate">https://chat.lazying.art/donate</a>
    </td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;">
      <a href="https://chat.lazying.art/donate"><img src="../figures/donate_button.svg" alt="Donate" height="44"></a>
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
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><img alt="WeChat QR" src="../figures/donate_wechat.png" width="240"/></td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><img alt="Alipay QR" src="../figures/donate_alipay.png" width="240"/></td>
  </tr>
</table>
</div>

**Message de remerciement**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります.  
- 你的支持将用于研发与运维，帮助我持续公开分享更多项目与改进.  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## Notes de version

- Les dossiers lourds (`books/`, symlinks Gaussian, checkpoints, temporaires…) sont ignorés via `.gitignore`.
- Concentrez vos contributions sur les répertoires suivis pour garder l’expérience de clonage « légère ».
- Pour mettre à jour le site : modifiez `docs/`, testez via `python -m http.server --directory docs`, puis poussez. GitHub Pages (ou autre) peut alors lier le domaine `learn.lazying.art` à ce dossier.
