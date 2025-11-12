<p>
  <b>Idiomas:</b>
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

LazyPhysics and Chemistry es la mitad de código y cuadernos de **LazyLearn**, mi diario de aprendizaje pausado en física y química. Los avances públicos viven en [learn.lazying.art](https://learn.lazying.art) (generado desde `docs/`), mientras que los experimentos reproducibles permanecen en este repositorio.

## LazyLearn

- **Base:** [learn.lazying.art](https://learn.lazying.art) publica el enfoque semanal, el backlog y los highlights.
- **Fuente única:** cada enlace apunta a `examples/`, `comp_physics/`, `multiwfn/` o `figures/`.
- **Flujo de actualización:** sube el código/cuadernos primero, regenera las figuras si hace falta y registra el cambio en `docs/` para mantener el sitio sincronizado.

## Qué contiene el repo

- `examples/` – scripts Python de QAOA/VQE que corren en un portátil con Qiskit o PennyLane.
- `comp_physics/` – notebooks de física computacional, scripts auxiliares (`numerov.py`) y datos/figuras.
- `comp_physics_python/` – port de los programas Fortran del libro *Computational Physics* organizado por capítulos (ver [comp_physics_python/README.md](../comp_physics_python/README.md)).
- `multiwfn/` – versión upstream de Multiwfn 3.8 dev con PDFs de referencia.
- `figures/` – PNGs generados para informes o presentaciones.
- `docs/` – micrositio de LazyLearn publicado en `learn.lazying.art` (GitHub Pages u otro host).

El resto (libros escaneados, symlinks Gaussian, checkpoints, etc.) queda excluido vía `.gitignore` para mantener el repo liviano.

## Entorno Python

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

Tras activar el entorno, abre `jupyter lab` o `jupyter notebook` para trabajar en `comp_physics/`.

## Flujos de ejemplo

- **QAOA con Qiskit** – `python examples/qaoa_qiskit_maxcut.py` (solo statevector).
- **QAOA con PennyLane** – `python examples/qaoa_pennylane_maxcut.py` (backend `default.qubit`).
- **VQE para H₂** – `python examples/pennylane_chemistry_h2_vqe.py` para reproducir `figures/pennylane_h2_vqe_convergence.png`.

Cada script registra métricas intermedias para facilitar la reutilización.

## Notebooks de física computacional

`comp_physics/` refleja la estructura de las notas:

- `comp_physics_textbook_code/` – rutinas reutilizables extraídas de los capítulos.
- Notebooks individuales: `chapter1.ipynb`, `chapter2.ipynb`, `numerov.ipynb`, `numpy_1ddft.ipynb`, etc.
- Carpetas temáticas (`bosonscattering/`, `lensless/`, `lightscattering/` …) con datos crudos y scripts.

Documenta cualquier dependencia extra en `comp_physics/environments.yaml`.

## Referencias Multiwfn

`multiwfn/` alberga `Multiwfn_3.8_dev_src_Linux`, el manual y la guía rápida. No se suben binarios; compila desde ahí.

## Figuras

Todas las imágenes PNG se encuentran en `figures/`. Incluye versiones nuevas aquí para tenerlas versionadas.

## Código del libro en Python

`comp_physics_python/` reúne las traducciones a Python de los programas Fortran de *Computational Physics*. Cada subcarpeta corresponde a un capítulo (por ejemplo `ch4/` = Hartree–Fock, `ch8/` = dinámica molecular, `ch10/` = Monte Carlo) y cada script ofrece una CLI para repetir los cálculos en un entorno moderno de Python/NumPy. Consulta [comp_physics_python/README.md](../comp_physics_python/README.md) para conocer el estado actual y los comandos.

## Apoya LazyLearn

Tu apoyo permite:

- Cubrir hosting/inferencia/almacenamiento de las demos y notebooks públicos.
- Financiar semanas de desarrollo abierto para EchoMind, LazyEdit y las utilidades cuánticas/físicas del repo.
- Prototipar hardware óptico y wearable (IdeasGlass, LightMind).
- Patrocinar despliegues gratuitos para estudiantes, laboratorios comunitarios y creadores.

### Donar

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

**Mensaje de agradecimiento**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります.  
- 你的支持将用于研发与运维，帮助我持续公开分享更多项目与改进.  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## Notas de control de versiones

- `books/`, symlinks de Gaussian, checkpoints y artefactos locales están excluidos en `.gitignore`.
- Mantén tus contribuciones dentro de los directorios seguidos para conservar una experiencia de clonación ligera.
- Para actualizar el sitio: edita `docs/`, prueba con `python -m http.server --directory docs` y luego haz push. GitHub Pages (u otro host) puede apuntar el dominio `learn.lazying.art` a esa carpeta.
