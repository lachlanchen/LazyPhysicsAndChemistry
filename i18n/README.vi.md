<p>
  <b>Ngôn ngữ:</b>
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

LazyPhysics and Chemistry là phần mã + notebook của **LazyLearn**—nhật ký học vật lý và hóa học theo nhịp chậm. Những cập nhật công khai xuất hiện trên [learn.lazying.art](https://learn.lazying.art) (được build từ `docs/`), còn các thí nghiệm có thể chạy lại đều nằm trong repo này.

## LazyLearn

- **Trạm chính:** [learn.lazying.art](https://learn.lazying.art) — nơi ghi chú trọng tâm tuần, backlog và highlight.
- **Nguồn sự thật:** mọi liên kết đều trỏ về `examples/`, `comp_physics/`, `multiwfn/` hoặc `figures/`.
- **Quy trình cập nhật:** đẩy code/notebook trước, tái sinh biểu đồ nếu cần rồi thêm mục mới ở `docs/` để website đồng bộ.

## Trong repo có gì

- `examples/` – script Python QAOA/VQE chạy được trên laptop với Qiskit hoặc PennyLane.
- `comp_physics/` – notebook vật lý tính toán, script trợ giúp (`numerov.py`) và dữ liệu/hình đi kèm.
- `multiwfn/` – bản nguồn Multiwfn 3.8 dev và tài liệu tham khảo.
- `figures/` – PNG sinh ra từ notebook/script.
- `docs/` – microsite LazyLearn, có thể deploy lên `learn.lazying.art` bằng GitHub Pages hoặc host riêng.

Những thư mục nặng (sách scan, symlink Gaussian, file checkpoint...) đều bị `.gitignore` loại trừ để repo nhẹ khi clone.

## Môi trường Python

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

Kích hoạt môi trường rồi chạy `jupyter lab`/`jupyter notebook` để mở các file trong `comp_physics/`.

## Quy trình mẫu

- **QAOA với Qiskit** – `python examples/qaoa_qiskit_maxcut.py` (chỉ dùng statevector).
- **QAOA với PennyLane** – `python examples/qaoa_pennylane_maxcut.py` (`default.qubit`).
- **VQE cho H₂** – `python examples/pennylane_chemistry_h2_vqe.py` để tái tạo `figures/pennylane_h2_vqe_convergence.png`.

Mọi script đều log chỉ số trung gian để dễ tái sử dụng.

## Notebook vật lý tính toán

`comp_physics/` phản chiếu cấu trúc sổ tay:

- `comp_physics_textbook_code/` – các đoạn mã tái sử dụng từ chương sách.
- Notebook độc lập như `chapter1.ipynb`, `chapter2.ipynb`, `numerov.ipynb`, `numpy_1ddft.ipynb`...
- Folder chủ đề (`bosonscattering/`, `lensless/`, `lightscattering/`...) gom dữ liệu gốc và script.

Nếu cần thêm dependency, hãy ghi vào `comp_physics/environments.yaml`.

## Tài liệu Multiwfn

`multiwfn/` lưu `Multiwfn_3.8_dev_src_Linux`, manual và quick-start. Không đính kèm binary—hãy build trực tiếp từ đây.

## Hình ảnh

Tất cả PNG nằm trong `figures/`. Khi có biểu đồ mới, hãy thả vào thư mục này để được version hóa.

## Ủng hộ LazyLearn

Ủng hộ giúp:

- Trang trải chi phí hosting/inference/storage cho demo công khai và notebook.
- Tài trợ các tuần làm việc tập trung cho EchoMind, LazyEdit, cùng những tiện ích lượng tử/vật lý khác.
- Thử nghiệm phần cứng quang học & wearable như IdeasGlass, LightMind.
- Tài trợ triển khai miễn phí cho sinh viên, phòng lab cộng đồng và creator.

### Donate

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

**Thông điệp ủng hộ**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります.  
- 你的支持将用于研发与运维，帮助我持续公开分享更多项目与改进.  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## Ghi chú version control

- `books/`, symlink Gaussian, file checkpoint và dữ liệu tạm đều bị `.gitignore` loại bỏ.
- Hãy tập trung đóng góp trong các thư mục nói trên để giữ repo nhẹ.
- Muốn cập nhật website: sửa `docs/`, preview bằng `python -m http.server --directory docs`, rồi push. GitHub Pages (hoặc host khác) có thể trỏ domain `learn.lazying.art` về thư mục này.
