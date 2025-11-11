<p>
  <b>言語:</b>
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

# LazyPhysics and Chemistry（レイジーフィジクス＆ケミストリー）

LazyPhysics and Chemistry は **LazyLearn**（物理と化学をゆっくり学ぶ記録）のコード＆ノート部分です。公開中の進捗やハイライトは [learn.lazying.art](https://learn.lazying.art)（`docs/` から生成）で確認でき、再現可能な実験とスクリプトは本リポジトリに置いています。

## LazyLearn

- **ホーム:** [learn.lazying.art](https://learn.lazying.art) に週次テーマ、バックログ、ハイライトを掲載。
- **ソース・オブ・トゥルース:** サイト上のリンクは `examples/`、`comp_physics/`、`multiwfn/`、`figures/` に収まります。
- **更新フロー:** まずコード/ノートをコミットし、必要なら図を再生成してから `docs/` にエントリを追加し、サイトを同期。

## リポジトリ構成

- `examples/` — QAOA/VQE の軽量 Python スクリプト（Qiskit or PennyLane）。
- `comp_physics/` — 計算物理ノート、`numerov.py` などの補助スクリプト、関連データ/図。
- `multiwfn/` — Multiwfn 3.8 dev ソースと PDF ドキュメント。Gaussian 後処理に使用。
- `figures/` — ノート/スクリプトから生成した PNG。
- `docs/` — LazyLearn のミニサイト。GitHub Pages（または自前ホスティング）で `learn.lazying.art` に公開。

その他の巨大ファイル（書籍スキャン、Gaussian シンボリックリンク、チェックポイント等）は `.gitignore` で除外し、軽量なクローンを保ちます。

## Python 環境

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

仮想環境を有効化後、`jupyter lab` もしくは `jupyter notebook` で `comp_physics/` 以下のノートを開きます。

## サンプルワークフロー

- **QAOA (Qiskit)** — `python examples/qaoa_qiskit_maxcut.py`（statevector のみ）。
- **QAOA (PennyLane)** — `python examples/qaoa_pennylane_maxcut.py`（`default.qubit`）。
- **H₂ VQE** — `python examples/pennylane_chemistry_h2_vqe.py` で `figures/pennylane_h2_vqe_convergence.png` を再現。

各スクリプトは途中経過をログ出力するため、別の分子やグラフへ展開しやすくなっています。

## 計算物理ノート

`comp_physics/` ディレクトリは紙のノート構成をそのまま再現：

- `comp_physics_textbook_code/` — 章ごとの再利用コード。
- `chapter1.ipynb`、`chapter2.ipynb`、`numerov.ipynb`、`numpy_1ddft.ipynb` などの独立ノート。
- `bosonscattering/`、`lensless/`、`lightscattering/` など、実験ごとのデータ/スクリプトをまとめたフォルダ。

追加依存が必要な場合は `comp_physics/environments.yaml` に追記してください。

## Multiwfn リファレンス

`multiwfn/` には `Multiwfn_3.8_dev_src_Linux`、マニュアル、クイックスタート PDF を収録。バイナリは含めず、必要に応じてここからビルドします。

## 図表

生成した PNG はすべて `figures/` に保存。再生成や新規の図があればここへ追加します。

## LazyLearn を支援

LazyLearn へのサポートが実現すること：

- 公開デモ/ノートのホスティング・推論・ストレージ費用のカバー。
- EchoMind、LazyEdit、量子/物理ツールの集中開発ウィークを後押し。
- IdeasGlass、LightMind などの光学/ウェアラブル試作。
- 学生、コミュニティラボ、クリエイター向けの無償デプロイ。

### 寄付

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

**サポートのメッセージ**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります。  
- 你的支持将用于研发与运维，帮助我持续公开分享更多项目与改进。  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## バージョン管理メモ

- `books/`、Gaussian シンボリックリンク、チェックポイント、ローカル一時ファイルなどは `.gitignore` で除外。
- 上記ディレクトリに集中してコントリビュートし、「ゆるくクローンできる」状態を維持。
- サイトを更新する際は `docs/` を編集し、`python -m http.server --directory docs` でプレビュー後に push。GitHub Pages 等で `learn.lazying.art` をこのディレクトリに向けます。
