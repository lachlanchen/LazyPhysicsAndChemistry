<p>
  <b>언어:</b>
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

LazyPhysics and Chemistry 는 **LazyLearn**(천천히 진행하는 물리·화학 학습 로그)의 코드/노트북 파트입니다. 외부 공개용 하이라이트와 진행 상황은 [learn.lazying.art](https://learn.lazying.art) (리포지터리 `docs/`에서 제공)에 게시되고, 재현 가능한 실험 코드는 이곳에 보관합니다.

## LazyLearn

- **홈베이스:** [learn.lazying.art](https://learn.lazying.art) — 주간 포커스, 백로그, 하이라이트를 확인.
- **단일 소스:** 사이트에 연결된 모든 요소는 `examples/`, `comp_physics/`, `multiwfn/`, `figures/` 안에 있습니다.
- **업데이트 흐름:** 코드/노트북을 먼저 푸시 → 필요한 경우 도표 재생성 → `docs/`에 항목 추가 → 사이트 자동 반영.

## 이 리포에 있는 것

- `examples/` – Qiskit 또는 PennyLane으로 구동되는 QAOA/VQE Python 스크립트.
- `comp_physics/` – 계산물리 노트북, `numerov.py` 같은 헬퍼, 데이터/그림.
- `comp_physics_python/` – Thijssen 교재 Fortran 코드를 장별로 옮긴 Python 모음 (자세한 내용은 [comp_physics_python/README.md](../comp_physics_python/README.md)).
- `multiwfn/` – Multiwfn 3.8 dev 소스와 PDF 매뉴얼.
- `figures/` – 노트/스크립트가 만들어낸 PNG.
- `docs/` – LazyLearn 사이트. GitHub Pages 등으로 `learn.lazying.art`에 배포.

기타 대용량 항목(도서 스캔, Gaussian 심볼릭 링크, 체크포인트 등)은 `.gitignore`에 포함되어 있어 클론이 가볍습니다.

## Python 환경

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

가상환경을 활성화한 뒤 `jupyter lab` 또는 `jupyter notebook`으로 `comp_physics/` 아래 노트북을 엽니다.

## 예제 워크플로

- **QAOA (Qiskit)** – `python examples/qaoa_qiskit_maxcut.py` (statevector 백엔드만 사용).
- **QAOA (PennyLane)** – `python examples/qaoa_pennylane_maxcut.py` (`default.qubit`).
- **H₂ VQE** – `python examples/pennylane_chemistry_h2_vqe.py`로 `figures/pennylane_h2_vqe_convergence.png` 재현.

모든 스크립트는 중간 지표를 기록해 새로운 분자/그래프 실험에 재활용하기 쉽습니다.

## 계산물리 노트북

`comp_physics/` 디렉터리는 노트 구조를 그대로 반영합니다.

- `comp_physics_textbook_code/` – 각 장에서 추출한 재사용 코드.
- 독립 노트북: `chapter1.ipynb`, `chapter2.ipynb`, `numerov.ipynb`, `numpy_1ddft.ipynb` 등.
- 주제 폴더: `bosonscattering/`, `lensless/`, `lightscattering/` 등. 원본 데이터와 스크립트 묶음.

추가 의존성이 필요하면 `comp_physics/environments.yaml`에 기록하세요.

## Multiwfn 참고

`multiwfn/`에는 `Multiwfn_3.8_dev_src_Linux`와 매뉴얼, 퀵스타트 PDF가 들어 있습니다. 바이너리는 포함하지 않습니다.

## 이미지

생성된 모든 PNG는 `figures/`에 저장합니다. 새 그림도 동일 위치에 보관하세요.

## 교재 코드 포팅

`comp_physics_python/`에는 Thijssen의 *Computational Physics* 예제가 장별로 정리된 Python 버전이 들어 있습니다(예: `ch4/` = Hartree–Fock, `ch8/` = MD, `ch10/` = Monte Carlo). 모든 스크립트가 CLI를 제공하므로 최신 Python/NumPy 환경에서 책 속 수치를 재현할 수 있습니다. 자세한 소개와 사용법은 [comp_physics_python/README.md](../comp_physics_python/README.md)를 참고하세요.

## LazyLearn 지원

후원은 다음을 가능하게 합니다.

- 공개 데모와 노트북을 위한 호스팅/추론/스토리지 비용 충당.
- EchoMind, LazyEdit, 양자·물리 유틸 개발 집중 주간 후원.
- IdeasGlass, LightMind 같은 광학·웨어러블 프로토타입 제작.
- 학생/커뮤니티 랩/크리에이터 대상 무상 배포.

### 기부하기

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

**지원 메시지**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります。  
- 你的支持将用于研发与运维，帮助我持续公开分享更多项目与改进。  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## 버전 관리 메모

- `books/` 디렉터리, Gaussian 심볼릭 링크, 체크포인트, 로컬 임시 산출물 등은 `.gitignore` 로 제외했습니다.
- 위에 언급한 디렉터리에 집중해 기여하면 가벼운 클론 경험을 유지할 수 있습니다.
- 사이트를 갱신할 땐 `docs/` 를 수정하고 `python -m http.server --directory docs` 로 미리보기 후 push 하세요. GitHub Pages 등에서 `learn.lazying.art` 를 해당 디렉터리에 연결하면 됩니다.
