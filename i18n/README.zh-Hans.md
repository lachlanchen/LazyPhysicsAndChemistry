<p>
  <b>语言:</b>
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

# LazyPhysics and Chemistry（懒惰物理与化学）

LazyPhysics and Chemistry 是 **LazyLearn** 的代码与笔记本一半——我用来缓慢学习物理与化学的日志。对外展示的进度、重点与待办发布在 [learn.lazying.art](https://learn.lazying.art)（由本仓库 `docs/` 目录生成），而所有可运行的实验与脚本都保存在这里，随时复现。

## LazyLearn

- **基地：** [learn.lazying.art](https://learn.lazying.art) 公布每周主题、积压任务与亮点。
- **唯一事实来源：** 网站上提到的内容都能在 `examples/`、`comp_physics/`、`multiwfn/` 或 `figures/` 找到对应代码/图像。
- **更新流程：** 先提交代码或笔记本，必要时重绘图，再在 `docs/` 中追加条目让网站同步。

## 仓库结构

- `examples/` —— 运行于普通笔记本电脑的 QAOA/VQE Python 脚本（依赖 Qiskit 或 PennyLane）。
- `comp_physics/` —— 计算物理笔记本、辅助脚本（如 `numerov.py`）以及配套数据/插图。
- `comp_physics_python/` —— Jos Thijssen 教科书 Fortran 代码的 Python 改写，按章节分目录（详见 [comp_physics_python/README.md](../comp_physics_python/README.md)）。
- `multiwfn/` —— Multiwfn 3.8 dev 版本源码与参考 PDF，方便配合 Gaussian 做后处理。
- `figures/` —— 由脚本或笔记本生成的 PNG，便于复用。
- `docs/` —— LazyLearn 微型站点，可通过 GitHub Pages（或自托管）发布到 `learn.lazying.art`。

其余大体量文件（图书扫描、Gaussian 符号链接、检查点等）被 `.gitignore` 排除，保证仓库克隆轻量。

## Python 环境

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

激活虚拟环境后运行 `jupyter lab` 或 `jupyter notebook`，`comp_physics/` 下的所有笔记本都会使用相同依赖。

## 示例流程

- **QAOA（Qiskit）** —— `python examples/qaoa_qiskit_maxcut.py`（仅用 statevector）。
- **QAOA（PennyLane）** —— `python examples/qaoa_pennylane_maxcut.py`（`default.qubit` 后端）。
- **VQE 求 H₂** —— `python examples/pennylane_chemistry_h2_vqe.py`，即可复现 `figures/pennylane_h2_vqe_convergence.png`。

脚本都会打印中间指标，便于复用或扩展到新分子/图。

## 计算物理笔记本

`comp_physics/` 目录复刻了纸质笔记结构：

- `comp_physics_textbook_code/` —— 从章节中抽出的可复用程序。
- 独立章节：`chapter1.ipynb`、`chapter2.ipynb`、`numerov.ipynb`、`numpy_1ddft.ipynb` 等。
- 主题文件夹：`bosonscattering/`、`lensless/`、`lightscattering/` 等，用于归档原始数据与脚本。

如需额外依赖，请在 `comp_physics/environments.yaml` 记录。

## Multiwfn 参考

`multiwfn/` 收录上游 `Multiwfn_3.8_dev_src_Linux`、使用手册与快速入门文档。没有提交任何二进制，可直接指向该目录进行编译或查阅。

## 图像

全部 PNG 存放在 `figures/`。若重生成或新增图表，请放入此处方便版本管理。

## 教程代码移植

`comp_physics_python/` 收录了《Computational Physics》中 Fortran 程序的 Python 版本。每个子目录对应书里的一个章节（例如 `ch4/` 是 Hartree–Fock，`ch8/` 是分子动力学，`ch10/` 是 Monte Carlo），并提供命令行入口以便在现代环境中直接复现实验。更多详情与运行示例见 [comp_physics_python/README.md](../comp_physics_python/README.md)。

## 支持 LazyLearn

支援 LazyLearn 可以：

- 承担示例/笔记本的托管、推理与存储成本。
- 资助 EchoMind、LazyEdit 以及本仓库量子/物理工具的集中开发周。
- 试制 IdeasGlass、LightMind 等光学/可穿戴硬件，为后续章节铺路。
- 为学生、社区实验室与创作者提供免费部署。

### 捐助渠道

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
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><strong>微信</strong></td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><strong>支付宝</strong></td>
  </tr>
  <tr>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><img alt="WeChat QR" src="../figures/donate_wechat.png" width="240"/></td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><img alt="Alipay QR" src="../figures/donate_alipay.png" width="240"/></td>
  </tr>
</table>
</div>

**支援说明**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります。  
- 你的支持将用于研发与运维，帮助我持续公开分享更多项目与改进。  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## 版本控制说明

- `books/`、Gaussian 符号链接、检查点文件及本地临时数据都在 `.gitignore` 中保持排除。
- 贡献请聚焦上述受控目录，确保克隆体验轻巧。
- 更新网站时，在本地修改 `docs/`，使用 `python -m http.server --directory docs` 预览，然后推送；GitHub Pages（或其他托管）即可将自定义域名 `learn.lazying.art` 指向该目录。
