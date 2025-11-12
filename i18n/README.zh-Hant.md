<p>
  <b>語言:</b>
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

# LazyPhysics and Chemistry（懶惰物理與化學）

LazyPhysics and Chemistry 是 **LazyLearn** 的程式碼與筆記本半邊——一份不趕進度的物理、化學學習誌。對外公開的每週焦點、成果與待辦記錄在 [learn.lazying.art](https://learn.lazying.art)（由 `docs/` 產生），而所有可重現的實驗都保留在本倉庫。

## LazyLearn

- **基地：** [learn.lazying.art](https://learn.lazying.art) 分享每週主題、積壓清單與亮點。
- **單一真實來源：** 網站上的連結均對應 `examples/`、`comp_physics/`、`multiwfn/` 或 `figures/`。
- **更新流程：** 先提交程式碼或筆記本、必要時重繪圖，再在 `docs/` 留下紀錄讓網站同步。

## 倉庫內容

- `examples/` —— QAOA 與 VQE Python 腳本，搭配 Qiskit 或 PennyLane 即可在筆電上運行。
- `comp_physics/` —— 計算物理筆記本、輔助腳本（如 `numerov.py`）以及配套資料/插圖。
- `comp_physics_python/` —— Jos Thijssen 教科書 Fortran 程式的 Python 版本，每章一個子資料夾（詳見 [comp_physics_python/README.md](../comp_physics_python/README.md)）。
- `multiwfn/` —— 上游 Multiwfn 3.8 dev 原始碼與 PDF 手冊，支援 Gaussian 後處理。
- `figures/` —— 由腳本/筆記本產出的 PNG，用於簡報或論文。
- `docs/` —— LazyLearn 迷你站點，透過 GitHub Pages（或自架）發佈到 `learn.lazying.art`。

其它大型內容（書籍掃描、Gaussian 符號連結、檢查點等）透過 `.gitignore` 排除，保持倉庫精簡。

## Python 環境

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

啟用虛擬環境後執行 `jupyter lab` 或 `jupyter notebook`，即可打開 `comp_physics/` 內的所有筆記本。

## 範例流程

- **QAOA（Qiskit）** —— `python examples/qaoa_qiskit_maxcut.py`（純 statevector）。
- **QAOA（PennyLane）** —— `python examples/qaoa_pennylane_maxcut.py`（`default.qubit`）。
- **VQE for H₂** —— `python examples/pennylane_chemistry_h2_vqe.py` 可重現 `figures/pennylane_h2_vqe_convergence.png`。

腳本會輸出中間指標，方便延伸至其它分子或圖。

## 計算物理筆記本

`comp_physics/` 重現紙本筆記的架構：

- `comp_physics_textbook_code/` —— 從章節萃取的可重用程式。
- 獨立章節：`chapter1.ipynb`、`chapter2.ipynb`、`numerov.ipynb`、`numpy_1ddft.ipynb` 等。
- 主題資料夾：`bosonscattering/`、`lensless/`、`lightscattering/` ……收集原始數據與腳本。

若需要額外依賴，請更新 `comp_physics/environments.yaml`。

## Multiwfn 參考

`multiwfn/` 內含 `Multiwfn_3.8_dev_src_Linux`、操作手冊與快速指南，沒有提交任何二進位檔，直接指向該資料夾即可。

## 圖像

所有 PNG 皆儲存在 `figures/`，若重新繪製或新增圖表，請放入此處。

## 教程程式移植

`comp_physics_python/` 正逐步收錄《Computational Physics》裡的 Fortran 範例，每個子資料夾對應原書一章（例如 `ch4/` 為 Hartree–Fock、`ch8/` 為分子動力學、`ch10/` 為 Monte Carlo）。所有腳本都提供 CLI，可以在現代 Python/NumPy 環境中直接重現書中的數值結果。詳情請見 [comp_physics_python/README.md](../comp_physics_python/README.md)。

## 支援 LazyLearn

支持 LazyLearn 代表：

- 支付公開示例與筆記本的託管、推理與儲存成本。
- 贊助 EchoMind、LazyEdit 以及量子/物理工具的密集開發時段。
- 打造 IdeasGlass、LightMind 等光學與穿戴原型，為後續章節鋪路。
- 為學生、社群實驗室與創作者提供免費部署。

### 捐助

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
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><strong>支付寶</strong></td>
  </tr>
  <tr>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><img alt="WeChat QR" src="../figures/donate_wechat.png" width="240"/></td>
    <td style="text-align:center; vertical-align:middle; padding:6px 12px;"><img alt="Alipay QR" src="../figures/donate_alipay.png" width="240"/></td>
  </tr>
</table>
</div>

**支援聲明**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります。  
- 你的支持將用於研發與運維，幫助我持續公開分享更多專案與改進。  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## 版本控制

- `books/`、Gaussian 符號連結、checkpoint 檔與本地暫存資料皆已在 `.gitignore` 排除。
- 請集中於上述受控目錄貢獻，維持「輕鬆克隆」體驗。
- 更新網站時，先在本地修改 `docs/`，用 `python -m http.server --directory docs` 預覽後再推送；GitHub Pages（或其他托管）即可把自訂網域 `learn.lazying.art` 指向此目錄。
