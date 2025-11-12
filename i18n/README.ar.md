<p>
  <b>اللغة:</b>
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

مشروع LazyPhysics and Chemistry هو الجانب البرمجي ودفاتر الملاحظات من **LazyLearn**—مذكّرة تعلم فيزياء وكيمياء بوتيرة هادئة. التحديثات العامة والمنجزات تعرض على [learn.lazying.art](https://learn.lazying.art) (يُبنى مباشرة من مجلد `docs/`)، بينما تبقى الأكواد والتجارب القابلة للإعادة في هذا المستودع.

## ما هو LazyLearn؟

- **المقر:** موقع [learn.lazying.art](https://learn.lazying.art) يضم التركيز الأسبوعي، قوائم المهام، وأبرز النتائج.
  
- **مصدر الحقيقة:** كل رابط على الموقع يعود إلى واحدة من المجلدات `examples/`، `comp_physics/`، `multiwfn/`، أو `figures/`.
  
- **طريقة التحديث:** ادفع الكود/الدفاتر أولاً، أعد توليد الرسوم إذا لزم الأمر، ثم حدّث `docs/` ليعكس الموقع آخر الأعمال.

## ماذا يوجد في المستودع؟

- `examples/` — سكربتات Python مركّزة لتجارب QAOA وVQE باستخدام Qiskit أو PennyLane.
- `comp_physics/` — دفاتر فيزياء حسابية، سكربتات مساعدة مثل `numerov.py`، وبيانات داعمة.
- `comp_physics_python/` — نسخ Python لبرامج Fortran الواردة في كتاب *Computational Physics*، مرتبة حسب الفصول (انظر [comp_physics_python/README.md](../comp_physics_python/README.md)).
- `multiwfn/` — إسقاط لمصدر Multiwfn 3.8 dev مع أدلة PDF مرجعية.
- `figures/` — صور PNG ثابتة تستخدم في العروض أو التقارير.
- `docs/` — موقع LazyLearn المصغّر الذي يُنشر إلى `learn.lazying.art` عبر GitHub Pages أو استضافة خاصة.

تظل الملفات الثقيلة (مسوحات الكتب، روابط Gaussian، ملفات checkpoint) خارج Git بفضل `.gitignore` لضمان خفّة المستودع.

## بيئة Python

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install qiskit pennylane numpy matplotlib
```

بعد تفعيل البيئة شغّل `jupyter lab` أو `jupyter notebook` لفتح دفاتر `comp_physics/`.

## أمثلة لتدفقات العمل

- **QAOA مع Qiskit** — الأمر `python examples/qaoa_qiskit_maxcut.py` (يعتمد على statevector فقط).
- **QAOA مع PennyLane** — الأمر `python examples/qaoa_pennylane_maxcut.py` (باستخدام `default.qubit`).
- **VQE لجزيء H₂** — `python examples/pennylane_chemistry_h2_vqe.py` لإعادة إنتاج الرسم `figures/pennylane_h2_vqe_convergence.png`.

كل سكربت يسجل القياسات الوسيطة لتسهيل إعادة الاستخدام.

## دفاتر الفيزياء الحسابية

مجلد `comp_physics/` يعكس بنية الملاحظات:

- `comp_physics_textbook_code/` — روتينات قابلة لإعادة الاستخدام مأخوذة من الفصول.
- دفاتر مستقلة مثل `chapter1.ipynb`، `chapter2.ipynb`، `numerov.ipynb`، `numpy_1ddft.ipynb`.
- مجلدات موضوعية (`bosonscattering/`, `lensless/`, `lightscattering/` ...) تجمع البيانات الخام والسكربتات لكل تجربة.

سجّل أي اعتمادات إضافية في `comp_physics/environments.yaml` قبل المشاركة.

## مراجع Multiwfn

يحتوي `multiwfn/` على الشجرة `Multiwfn_3.8_dev_src_Linux` مع دليل الاستخدام والبدء السريع. لا تُضمَّن ملفات تنفيذية؛ ابنِ المشروع مباشرةً من هناك.

## الرسوم والصور

جميع ملفات PNG محفوظة داخل `figures/`. أضف أي رسومات جديدة لهذا المجلد ليتم تتبعها بالإصدارات.

## ترجمات كود الكتاب

يحتوي مجلد `comp_physics_python/` على نسخ Python من برامج Fortran الواردة في كتاب *Computational Physics*. كل مجلد فرعي يمثّل فصلاً (مثلًا `ch4/` = Hartree–Fock، و `ch8/` = ديناميكيات جزيئية، و `ch10/` = مونت كارلو)، وجميع السكربتات توفر واجهة أوامر لتشغيل التجارب على بيئة Python/NumPy حديثة. راجع [comp_physics_python/README.md](../comp_physics_python/README.md) للاطلاع على حالة التغطية وتعليمات التشغيل.

## دعم LazyLearn

يساعد دعمك على:

- تغطية تكاليف الاستضافة/الاستدلال/التخزين للتجارب والدفاتر العامة.
- تمويل أسابيع التطوير المركّز لـ EchoMind وLazyEdit وأدوات الفيزياء/الكم.
- بناء نماذج أولية للأجهزة الضوئية والقابلة للارتداء مثل IdeasGlass وLightMind.
- توفير نشر مجاني للطلبة ومختبرات المجتمع وصنّاع المحتوى.

### التبرع

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

**رسالة الدعم**

- ご支援は研究・開発と運用の継続に役立ち、より多くのオープンなプロジェクトを皆さんに届ける力になります.  
- 你的支持将用于研发与运维，帮助我持续公开分享更多项目与改进.  
- Your support sustains my research, development, and ops so I can keep sharing more open projects and improvements.

## ملاحظات التحكم بالإصدارات

- مجلد `books/` وروابط Gaussian وملفات checkpoint والمخرجات المؤقتة كلها مستبعدة عبر `.gitignore`.
- ركّز المساهمات على المجلدات المذكورة للحفاظ على سهولة الاستنساخ.
- لتحديث الموقع: عدّل `docs/` محليًا، جرّب عبر `python -m http.server --directory docs`، ثم ادفع التغييرات. يمكن توجيه النطاق `learn.lazying.art` إلى هذا المجلد بواسطة GitHub Pages أو أي استضافة أخرى.
