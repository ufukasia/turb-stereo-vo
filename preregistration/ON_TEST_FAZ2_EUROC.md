# Faz 2 EuRoC ön-testleri — beklenen değerler ÖNCE yazıldı

Tarih: 2026-11-07 (gece). Platform: `pyeco_turb`. Ön-kayıt: `ON_KAYIT_FAZ2.md` (taslak v0.1).
Bu dosyanın §1–§3 bölümleri **hiçbir ön-test sonucu görülmeden** yazıldı. Sonuçlar §4'e eklenir;
§1–§3'te sonradan değişiklik yapılmaz (gerekirse §4'te "sapma" olarak yazılır).

## 0. Ortam ve sabitler

- Python: `turbulans vio/kod/turbvio/bin/python`; `OMP_NUM_THREADS=4`, `SLAM_DATA_ROOT=<data-disk>/sdata`.
- Yapılandırma: `pyeco_turb/config.yaml` (yerinde değiştirilmedi; hücre ayarları `turb/runner.py::_build_cfg` ile).
  `turbulence.psf_legacy_scale: false` (fiziksel PSF) config'te zaten böyle; runner `Cell.psf_legacy_scale=False`.
- **GPU eşlemesi (SAPMA, ölçüldü):** görev metnindeki `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1`
  bu makinede **RTX 3090**'ı seçiyor (torch ile doğrulandı) ve simülatör (Africa_001, bus 2B) o kartta.
  RTX 4060 Ti = `PCI_BUS_ID` sırasında **0**. Bütün ön-testler `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=0`
  ile 4060 Ti'da koşuldu; her kayıtta `fe_gpu` alanı kartın adını taşır.
- `turb/experiments.py::EUROC` bu makinede olmayan bir yolu (`<home>/Desktop/fable/sdata/`) gösteriyor;
  ön-testler `base_path=$SLAM_DATA_ROOT` ile kuruldu.
- Pencereler: `turb/pencere_sec.py` → `pyeco_turb/evidence/pencereler_faz2.json` (yalnız GT; kural betiğin başında).

## 1. Denetimler ve beklenen değerler

| # | denetim | beklenen (koşumdan önce) | geçme ölçütü |
|---|---|---|---|
| T1 | Pencere seçimi (yalnız GT) | 9 sekansın her birinde 3 örtüşmeyen 400 karelik pencere, hepsi ≥ 5 m; MH01/500 (0.075 m) listede yok | her pencere ≥ 5 m, örtüşme yok, başlangıç ≥ 100, GT kapsamı içinde |
| T2 | Regresyon kapısı (runner değişikliği) | Eski varsayılan (`WINDOW_STARTS=(500,900,1300,1600)`, MH01, D/r₀=0, A0, NEES açık) ile değişiklik öncesi ve sonrası `run_cell` kayıtları **bit düzeyinde aynı** (`m_*`, `ms_*`, `n_*` alanları). `code_fp` farklı (runner.py parmak izine dahil; beklenen). `Cell().windows`, `_arm_argv`, `_build_cfg` aynı. | metrik alanlarında fark 0 |
| T3 | Fizik denetimi (`turb.sanity --physics`, fiziksel PSF) | S1/S1b/S2/S4 ✅; S6 kırınım FWHM oranı 1.00±0.15 ✅; S7 EE50 oranı 1.00±0.20 ✅ (D/r₀ ≤ 12'de). S3/S5/S7b uyarı (⚠️) olabilir (önceki koşumlarda bilinen). D/r₀ 16–24 için ayrıca PSF tablosu: EE50 sim/teori oranı 0.8–1.2 | 🛑 yok |
| T4 | Preflight P-maddeleri (`turb.preflight --report_only`, duman kayıtları üzerinde) | P1 ✅ (çekirdek kollarda başarısızlık < %10), P2 ✅, P4 ✅, P7 ✅ (tek code_fp), P8 ✅ (ATE < 10 m), P9/P10/P12 ✅ (kayıttan bağımsız birim denetimleri). P5: D/r₀=0'da A0~A1 aynı (tasarım gereği), ama bütün kayıtlar üzerinde medyan fark > 1e-3 olduğundan ✅ beklenir. P13 (D/r₀=0 A1 vs A5 bağıl fark < 1e-3): A5 yama yolu faktörleri kare sonunda eklediği için **bit-aynı değil**; ✅ ya da küçük farkla 🛑 olabilir → ölçüt: < %10 ise engelleyici değil. P14 (σ̂/oracle 0.5–2): A3'te D/r₀=16'da ✅ beklenir; D/r₀=0'da oracle=0 olduğu için denetim dışı. E-maddeleri bilgi amaçlı (eski robust rejimine göre yazılmış). | P1/P2/P4/P7/P8 🛑 yok |
| T5 | A0C k kalibrasyonu (V1_01, V1_02, V2_01, V2_02; D/r₀=0; 3 pencere/sekans = 12 pencere) | V-serisi temiz A0 medyan konum NEES'i 50–500 (MH01 pilot 367; Vicon GT daha iyi olduğundan biraz daha düşük olabilir) → seçilecek k ∈ {8, 12, 16}. Ölçekleme eşdeğerliği: gerçek A0C(k) koşumunun NEES'i A0/k²'nin ±%30'u içinde (Huber eşiği beyazlatılmış birimde → tam eşdeğer DEĞİL, yaklaşık). k seçimi **gerçek A0C koşumlarından** yapılır (yedi k'nin hepsi koşulur); ölçekleme yalnız karşılaştırma için. | k ızgarada tanımlı; A0C'nin bütün pencereleri `ok` |
| T6 | Duman: çökme yok | 2 sekans (MH02, V1_01) × 3 D/r₀ (0, 8, 16) × 7 kol (A0, A0C(yeni k), A1, A3, A5, AG, A0C+A3) × 3 pencere = 126 kayıt, hepsi `status=ok` | 126/126 ok |
| T7 | Duman: alanlar dolu | her kayıtta `m_ate_se3`, `ms_ate_se3`, `m_rpe_trans`, `n_nees_pos_median`, `n_nees_rot_median`, `n_n`, `t_backend`, `arm_stats` dolu ve sonlu | eksik/NaN 0 |
| T8 | D/r₀=0: A1 == A0 | `m_ate_se3`, `ms_ate_se3`, `n_nees_pos_median` **bit düzeyinde aynı** (backend.py: A1 sigma_oracle=0 → σ = σ₀, aynı faktör sırası) | fark = 0 |
| T9 | D/r₀=0: A5 ≈ A0 | ATE bağıl farkı < %10, NEES bağıl farkı < %25 (faktör sırası farklı, sayısal olarak küçük ayrışma) | bu bantta |
| T10 | Türbülans gerçekten uygulandı | D/r₀=0'da `quality` = None (enjektör kurulmaz). D/r₀=8 ve 16'da `quality` dolu; PSNR(16) < PSNR(8) (en az 1 dB); PSNR(8) < 40 dB. A1/A5 `arm_stats.sigma_px_oracle` > 0 ve D/r₀=16'da 8'dekinden büyük (σ ∝ (D/r₀)^{5/6} → oran ≈ 1.78). Ön-uç dökümü D/r₀ ile değişir (eşleşme sayısı farklı). | hepsi |
| T11 | Temiz ATE makul (hareketli pencere) | D/r₀=0 A0 pencere ATE'si (SE(3) hizalı, gerçek zamanlı): 0.5–6 cm; PyEco'nun tam sekans değerlerinden (MH02 0.047, V1_01 0.050 m) **küçük ya da benzer**, çünkü pencere 20 s. Hiçbir temiz pencere > 0.10 m değil. | ≤ 0.10 m |
| T12 | D/r₀=8 dayanıklılık (yalnız gözlem; H1'in testi değil) | A0 ATE(8)/ATE(0) pencere medyanı ≈ 1 (0.8–1.5); D/r₀=16'da artış olabilir, ıraksama (max_err > 5 m) beklenmez | bilgi |
| T13 | NEES | Temiz A0 medyan konum NEES ≫ 3 (beklenen 50–1000). NEES örnek sayısı `n_n` (cov_every=10, 400 kare): pencere başına 15–40 (anahtar kare kipinde örnek sayısı < 40). A0C(k) NEES ≈ 2.37 civarı (V1_01'de; kalibrasyon kümesinin parçası). AG ve A3 NEES < A0 NEES (her D/r₀'da; ikisi de σ'yı büyütür). D/r₀ ≥ 8'de A1 NEES < A0 NEES. A0C+A3 NEES ≤ A0C NEES. | yön denetimleri |
| T14 | Süre/disk | Ön-uç (4060 Ti, ~1400 kare): 60–150 s. Arka-uç kol başına (~1400 kare): A0/A0C/A1/AG 1–4 dk; A3/A0C+A3 (DIMM) 1.2–2× A0; A5 (yama) 2–5× A0. Kayıt başına < 20 kB; geçici döküm hücre başına < 0.5 GB. | ölçülür, §5'te tahmin |

## 2. Koşum düzeni

1. T2 öncesi koşumu (kod değişikliğinden önce): `faz2/regresyon.py` → `regresyon_once.json`.
2. Runner/deney değişikliği → T2 sonrası koşumu `regresyon_sonra.json`, karşılaştırma.
3. T3: `python -m turb.sanity --physics` + D/r₀ ∈ {16, 20, 24} için PSF tablosu.
4. T5: V-serisi D/r₀=0 ön-uçları (kalıcı dizin, sonra duman testinde V1_01 için yeniden kullanılır);
   A0 (NEES) + A0C k ∈ {1,2,4,8,12,16,24}.
5. T6–T13: MH02, V1_01 × D/r₀ {0, 8, 16} × 7 kol, tohum 0, NEES açık (`cov_every=10`).
   Pencereler: o sekansın seçilmiş 3 penceresinin hepsi (görevdeki "bir pencere"yi kapsar; hücre maliyeti
   zaten son pencerenin sonuyla belirlendiği için ölçülen süre gerçek ızgaradaki hücre süresidir).
6. T4: `python -m turb.preflight --report_only --records <duman kayıtları>`.

## 3. Engelleyici sayılacak durumlar (uzun koşum başlatılmaz)

- Herhangi bir kolda çökme (T6), boş metrik alanı (T7), T8'in bit-aynılığının bozulması, T2 regresyonunun bozulması.
- Türbülansın uygulanmadığına dair işaret (T10).
- Temiz ATE > 0.10 m (T11) ya da ıraksama.
- k'nin ızgara sınırında (1 ya da 24) çıkması (kalibrasyon ızgarası yetersiz demektir).

---

## 4. Sonuçlar (koşumdan SONRA eklendi; §1–§3 değiştirilmedi)

### 4.1 T1 — Pencere seçimi (`evidence/pencereler_faz2.json`, sha256 `70ae3486…8b01`; betik sha256 `a5108f23…4851`)

| sekans | başlangıçlar | GT yol [m] | kapsam (sınır kutusu köşegeni) [m] | son pencere sonu |
|---|---|---|---|---|
| MH01 | 100, 721, 1121 | 5.44, 5.03, 11.30 | **0.96**, 4.16, 7.89 | 1521 |
| MH02 | 100, 600, 1000 | 7.81, 5.03, 10.01 | **0.91**, 4.47, 7.72 | 1400 |
| MH03 | 162, 562, 962 | 5.05, 28.86, 26.11 | 3.81, 9.58, 10.06 | 1362 |
| MH04 | 132, 532, 932 | 5.01, 24.68, 23.10 | 3.39, 14.63, 16.14 | 1332 |
| MH05 | 102, 502, 902 | 5.03, 21.24, 23.69 | 2.51, 14.44, 20.40 | 1302 |
| V1_01_easy | 100, 500, 900 | 6.42, 6.75, 7.89 | 3.61, 2.43, 4.14 | 1300 |
| V1_02_medium | 100, 500, 900 | 19.91, 20.12, 17.02 | 6.27, 5.19, 5.88 | 1300 |
| V2_01_easy | 100, 500, 900 | 8.35, 5.96, 5.92 | 4.86, 2.88, 2.15 | 1300 |
| V2_02_medium | 100, 500, 900 | 15.29, 20.53, 11.75 | 5.84, 6.39, 4.14 | 1300 |

Sonuç: **GEÇTİ** (27 pencere, hepsi ≥ 5 m, örtüşme yok, MH01/500 yok).
**Karar gerektiren bulgu:** MH01/100 ve MH02/100, ön-kayıttaki "yol ≥ 5 m" ölçütünü IMU uyarım sallamasıyla
yerinde salınarak geçiyor (kapsam < 1 m, 20 s'de 250–350° toplam dönme). B6'daki "durağan pencere" sorunu değil
ama "hareketli pencere" amacına da tam uymuyor. Seçenek: `--min_kapsam 1.0` eklenirse MH01 → 721/1121/1521,
MH02 → 600/1000/1400 olur (öteki sekanslar değişmez; `evidence/pencereler_faz2_kapsam1m.json` hazır).
Ön-kayıt dondurulmadan önce kullanıcı seçmeli. Duman testi her iki listede de bulunan pencerelerle
yorumlanmalı (MH02/600, MH02/1000; V1_01'in üçü).

### 4.2 T2 — Regresyon kapısı

`turb/runner.py` değişikliği: yalnız `windows_for(seq, path=None)` eklendi + belge satırı (yedek `runner.py.yedek_faz2`).
`turb/experiments.py`: `cells_FAZ2` + `FAZ2` kaydı (yedek `experiments.py.yedek_faz2`); eski tanımlar dokunulmadı.

| | önce | sonra |
|---|---|---|
| `WINDOW_STARTS`, `Cell().windows`, `_arm_argv`, `_build_cfg` | (500,900,1300,1600) | aynı (4/4 eşit) |
| `run_cell` MH01 D/r₀=0 A0, 4 pencere, 180 alan | — | **0 fark** (yalnız `frontend_stats.total_execution_time_seconds` ve `t_*` süreleri farklı) |
| `code_fp` | dcead3bb7b42 | 65dbe4761b2b (beklenen: runner.py parmak izinde) |
| ATE (500/900/1300/1600) | 0.0004 / 0.0070 / 0.0135 / 0.0127 m | aynı |

Sonuç: **GEÇTİ.** (Not: eski 500 penceresinde yol 0.08 m, ATE 0.4 mm → B6 yeniden üretildi.)

### 4.3 T3 — Fizik denetimi (`evidence/faz2_on_test/sanity_fizik.md`, `psf_tablo_16_24.json`)

S1 0.950 ✅, S1b ✅, S2 üstel 1.736 ✅, S3 0.957/0.860 ✅, S4 en büyük sapma %11 ✅, S5 0.631 ✅,
S6 FWHM oranı 1.011 ✅ (`psf_legacy_scale=False`; eski kip oranı 4.78), S7 EE50 oranı 0.96–1.01 ✅ (D/r₀ ≤ 12),
S7b ⚠️ (FWHM 8–12'de 0.58; beklenen uyarı). Genel: SARI, 🛑 yok → **GEÇTİ.**

Faz 2'nin üst ucu için PSF tablosu (EE50 simülasyon/Fried kısa-poz):

| D/r₀ | 8 | 12 | 16 | 20 | 24 |
|---|---|---|---|---|---|
| EE50 sim / teori [px] | 6.43/6.37 | 9.65/9.86 | 12.77/13.43 | 15.75/17.06 | 18.17/20.72 |
| oran | 1.01 | 0.98 | 0.95 | 0.92 | **0.88** |

Bant (0.8–1.2) içinde → geçti, ama oran D/r₀ ile düzenli düşüyor. Olası nedenler: çekirdek 31 px'te kesiliyor
(D/r₀=24'te EE50 çapı 18 px; hale enerjisi kırpılıyor) ve ince ekranda r₀ = 32/24 = 1.3 örnek (az örnekleme).
Sonuç: D/r₀ 20–24'te bulanıklık teoriden %8–12 dar, yani bozulma hafif **eksik** üretiliyor. Engelleyici
değil, ama makalede sınırlılık olarak yazılmalı ya da uzun koşumdan önce çekirdek boyu büyütülmeli
(büyütmek ön-uç dökümlerini ve süreyi değiştirir).

### 4.4 T5 — A0C k kalibrasyonu (V-serisi, D/r₀=0, 12 pencere; yedi k'nin hepsi gerçekten koşuldu)

| k | gerçek A0C: pencere medyanlarının medyanı (NEES_pos) | havuzlanmış örnek medyanı (≈169 örnek) | A0/k² (ölçekleme) | gerçek / ölçekli (pencere medyanı; min–maks) | NEES_rot | ATE A0C/A0 |
|---|---|---|---|---|---|---|
| 1 | 118.52 | 107.73 | 118.52 | 1.00 (1.00–1.00) | 119.78 | 1.000 |
| 2 | 34.87 | 32.08 | 29.63 | 1.08 (0.91–1.79) | 31.34 | 0.964 |
| 4 | 11.27 | 11.21 | 7.41 | 1.37 (1.02–3.11) | 8.76 | 1.000 |
| 8 | 3.07 | 3.20 | 1.85 | 1.95 (1.01–4.08) | 2.94 | 0.954 |
| **12** | **2.54** | **2.07** | 0.82 | 2.64 (1.01–6.73) | 1.33 | 1.008 |
| 16 | 1.34 | 1.50 | 0.46 | 2.82 (1.36–9.17) | 0.78 | 1.023 |
| 24 | 0.96 | 0.92 | 0.21 | 3.77 (1.81–36.2) | 0.41 | 1.051 |

- **Seçilen k = 12** (2.37'ye en yakın: 2.54; havuzlanmış örnek medyanı tanımıyla da 12 → 2.07). Eski 16 değil.
- **Ölçekleme eşdeğerliği YOK (beklenti T5 ±%30 → KALDI):** k ≥ 4'te gerçek A0C NEES'i A0/k²'nin 1.4–3.8 katı
  (tek pencerede 36×). Neden: Huber eşiği beyazlatılmış birimde sabit; σ k katına çıkınca Huber fiilen kapanır,
  aykırı gözlemler tam ağırlık alır ve göreli poz kestirimi değişir (ATE değişmiyor, ±%5). Dolayısıyla A0C'yi
  A0 dökümünden ölçekleyerek hesaplamak **geçersiz**; k ancak gerçek koşumla seçilebilir (bu yapıldı).
  Bu, ON_KAYIT §0'daki pilot tablosunun ("k=16 A0C dökümünden A0'a geri ölçeklenerek") da yaklaşık olduğunu gösterir.
- V-serisi temiz A0 NEES'i sekanslar arasında çok heterojen: V1_01 pencerelerinde 316–1514 (dönme NEES 617–1331),
  V2'de 25–128. 12 pencerenin 5–23 NEES örneği var (anahtar kare kipi, `cov_every=10`) → k'nin belirsizliği büyük:
  k=8 (3.07) ile k=12 (2.54) arasındaki fark gürültü mertebesinde.
- k=1 A0C, A0 ile birebir aynı (beklenen plasebo).

### 4.5 T6–T13 — Duman testi (MH02 + V1_01_easy × D/r₀ {0, 8, 16} × 7 kol × 3 pencere, tohum 0, k=12)

Kayıtlar: `pyeco_turb/evidence/faz2_on_test/duman_kayitlari.jsonl` (126 kayıt; `_tum` dosyasında ek `cov_every` denemeleri).
Kod parmak izi bütün kayıtlarda `65dbe4761b2b`. Ön-uç kartı her kayıtta `NVIDIA GeForce RTX 4060 Ti`.

A0 temiz ve türbülanslı pencere ATE'si (SE(3) hizalı, gerçek zamanlı yörünge) [m]:

| sekans / pencere | D/r₀=0 | D/r₀=8 | D/r₀=16 | oran 8/0 | oran 16/0 |
|---|---|---|---|---|---|
| MH02/100 (sallama penceresi) | 0.0056 | 0.1039 | 0.6270 | 18.5 | 112 |
| MH02/600 | 0.0124 | 0.0442 | 0.3847 | 3.6 | 31 |
| MH02/1000 | 0.0123 | 0.0725 | 0.2252 | 5.9 | 18 |
| V1_01/100 | 0.0288 | 0.0645 | 0.2057 | 2.2 | 7.1 |
| V1_01/500 | 0.0271 | 0.1139 | 0.3939 | 4.2 | 15 |
| V1_01/900 | 0.0185 | 0.1201 | 0.3565 | 6.5 | 19 |

Medyan konum NEES (pencere medyanı; aralık 6 pencere üzerinden):

| kol | D/r₀=0 | D/r₀=8 | D/r₀=16 |
|---|---|---|---|
| A0 | 101–1514 | 75–1281 | 539–2874 |
| A0C (k=12) | 1.1–14.6 | 1.2–25.1 | 13–59 |
| A1 | = A0 (bit-aynı) | 7–56 | 21–98 |
| A3 | 35–1049 | 8–147 | 100–321 (MH02) |
| A5 | = A0 (bit-aynı) | 3.6–60 | 25–100 (MH02) |
| AG | 65–1219 | 37–905 | 51–771 (MH02) |
| A0C+A3 (k=12) | 0.9–14.5 | 0.9–24 | 25–54 (MH02) |

Denetimler:

- **T6 çökme yok:** 126/126 `ok` → **GEÇTİ.**
- **T7 alanlar dolu:** eksik/NaN 0 → **GEÇTİ.** (A0'da `arm_stats` = None tasarım gereği: A0'da kol nesnesi yok,
  PyEco yolu birebir. Beklenti tablosundaki "her kayıtta arm_stats" ifadesi A0 için yanlıştı.)
- **T8 D/r₀=0 A1 == A0:** 6/6 pencerede ATE, düzleştirilmiş ATE, RPE, NEES_pos, NEES_rot bit düzeyinde aynı → **GEÇTİ.**
- **T9 D/r₀=0 A5 ≈ A0:** 6/6 pencerede bağıl fark 0 (bit-aynı; beklenen bantın da içinde) → **GEÇTİ.**
  (Sigma=0'da yamalar korelasyonsuz → bireysel faktör yolu; preflight P13 de 8e-14 ile geçti.)
- **T10 türbülans uygulandı:** D/r₀=0'da `quality=None`; PSNR MH02 20.8 → 17.8 dB, V1_01 22.9 → 20.7 dB (her ikisi < 40,
  düşüş ≥ 2.2 dB); oracle σ 2.19 → 3.90 px (oran 1.78 = 2^{5/6}); ön-uç ortalama sağ-görüntü kp sayısı ve stereo tetik
  oranı D/r₀ ile artıyor (MH02 stereo_run_ratio 0.28 → 0.51 → 0.98) → **GEÇTİ.**
- **T11 temiz ATE:** MH02 0.6–1.2 cm, V1_01 1.8–2.9 cm; hepsi ≤ 0.10 m → **GEÇTİ.** (Kalibrasyon kümesinde V2_02/100: 9.8 cm,
  sınıra yakın ama altında.)
- **T12 D/r₀=8 dayanıklılığı (yalnız gözlem):** beklenen 0.8–1.5 → ölçülen **2.2–18.5 (6/6 pencere > 1.5)** → **BEKLENTİ TUTMADI.**
  D/r₀=16'da 7–112×; ıraksama yok (max_err ≤ 1.1 m < 5 m). Bu, H1'in ilk yarısının ("D/r₀ ≤ 8'de ≤ 1.2×; > 1.5 ise çürür")
  tam ızgarada çürüyeceğine işaret ediyor ve önceki ön-testteki "D/r₀ ≤ 8'de bozulma yok" gözlemiyle (eski hat) çelişiyor.
  **Dikkat (ön-kayıt bütünlüğü):** ön-kayıt henüz dondurulmadı; H1'in bu sonuç görüldükten sonra değiştirilmesi sonradan
  uydurma olur. Öneri: H1 olduğu gibi dondurulsun, bu duman gözlemi §0'a "pilot olgu" olarak eklensin.
- **T13 NEES:** temiz A0 NEES 101–1514 ≫ 3 ✅. A3 < A0: 18/18 ✅. A1 < A0 (D/r₀ ≥ 8): 12/12 ✅. AG < A0: 17/18 (MH02/8/100'de
  AG 905 > A0 828) — çoğunlukla ✅. A0C+A3 ≤ A0C: 13/18 — **kısmen** (D/r₀ ≥ 8'de A3 bileşeni A0C+A3'ün kestirimini değiştiriyor;
  yön kuralı mekanik bir gereklilik değilmiş, beklenti fazla iddialıydı). **NEES örnek sayısı:** pencere başına 5–34
  (medyan 11) → beklenen 15–40 **TUTMADI**: döküm anahtar kare başına `cov_every`'de bir yapılıyor (≈ 0.26 anahtar kare/kare);
  `cov_every=10` → pencere başına ≈ 10 örnek. Ek deneme (MH02 D/r₀=0 A0): `cov_every=2` → 39–80, `cov_every=1` → 78–159
  örnek/pencere; NEES medyanları değişmiyor (109/101/140 → 147/102/141). Süre farkı yük altında ölçülemedi
  (776 / 485 / 681 s, makine paylaşımlı; aynı hesaplamanın A0 ve A1 kopyaları 776 s ile 394 s sürdü).
- **σ̂ (A3, DIMM) / oracle:** D/r₀=8'de 0.59–0.68, D/r₀=16'da 0.44–0.46 (P14 toplu 0.53, sınırda geçiyor). D/r₀=0'da
  σ̂ = 0.56–0.87 px (türbülans yokken de şişirme; A3'ün temiz NEES'ini A0'ın altına çekiyor).
  A1/A5 toplam σ'yı hypot(1, 1.75)·σ_oracle = 2.02·σ_oracle alıyor → D/r₀=16'da A1 7.9 px, A3 1.8 px.

### 4.6 T4 — Preflight (`turb.preflight --report_only`, duman kayıtları; `evidence/faz2_on_test/preflight_duman.md`)

| madde | sonuç | yorum |
|---|---|---|
| P1 çökme | ✅ (bütün kollar %0) | |
| P2 NEES verimi | 🛑 medyan 11 / beklenen ≥ 30 | eşik kare tabanlı (400/10=40); PyEco dökümü anahtar kare tabanlı. Mekanik kusur değil, **protokol kararı**: `cov_every` 2'ye (ya da 1'e) indirilmeli |
| P3, P4, P6, P7, P8, P9, P10 | ✅ | P7 tek code_fp; P8 en büyük ATE 0.63 m |
| P5 kol ayırt ediciliği | ✅ (10 ikili farklı) | A0~A1~A5 D/r₀=0'da tasarım gereği aynı, toplamda farklı |
| P11 robust simetrisi | 🛑 "sayaç yok" | **eski denetim**: PyEco kolunda `n_fac_robust` sayacı yok. Yerine her kolun `arm_stats.robust` alanı aynı (mono 2.4477, stereo 2.7955, gözlem başına) ve T8/T9 bit-aynılığı robust yolun ortak olduğunu gösteriyor → dolaylı olarak geçti; preflight'a PyEco sayacı eklenmeli |
| P12 Schur | 🛑 "TypeError" | **eski denetim**: PyEco'da `_pose_marginal_schur(gfg, keys, dims, n_pose_keys)` imzası farklı. Aynı test PyEco imzasıyla koşuldu (`faz2_on_test/p12_pyeco_imzasi.txt`): bireysel ve yama (bağlı Hll) yollarında yoğun tersle bağıl fark 2.0e-16 → **geçti**; preflight güncellenmeli |
| P13 | ✅ (8e-14) | |
| P14 | ✅ (σ̂/oracle 0.53, saçılım 0.13) | sınırda |
| E1–E4 | ✅ | bilgi amaçlı (A0/A3 5.3×, A1/A5 1.30×, A0/A0C 57×, A5/A0 ATE 1.07×; 4 eşleşmiş hücre) |

Kapı kararı: "KALDI" (P2, P11, P12). P11 ve P12 denetimin PyEco'ya taşınmamasından (kod yolu doğrulandı); P2 gerçek bir
protokol eksikliği (pencere başına ≈ 10 NEES örneği).

### 4.7 Uçtan uca hat testi (gerçek uzun-koşum yolu)

`experiments.cells_FAZ2()` → `runner.run_experiment_pipelined` (spool, yığın, iki faz), V2_01_easy × D/r₀ {0, 4} × {A0, A0C},
`evidence/faz2_on_test/hat_testi.jsonl`: 12/12 kayıt `ok`, `backend_argv=--turb_cal_factor 12.0` kayıtta, A0C `arm_stats.cal_factor=12.0`.
D/r₀=0 kayıtları ön-test sürücüsünün ürettiği kayıtlarla **bit düzeyinde aynı** (ATE, NEES_pos/rot) → sürücüyle elde edilen
duman sonuçları uzun-koşum yoluna aktarılabilir. (Gözlem: D/r₀=4'te bile A0 ATE oranı 1.2–3.5×, PSNR 27.1 dB.)

### 4.8 T14 — Süre ve disk tahmini (tam ızgara, ON_KAYIT §1)

Izgara: 9 sekans × 3 pencere; D/r₀=0 tek tohum, D/r₀ ∈ {4..24} × tohum 0–4 → **279 hücre**; kollar A0, A0C, A1, A3, AG,
A0C+A3 (+ A5 yalnız MH02, MH04) → **1736 arka-uç koşumu**, 5208 pencere kaydı. Hücre başına kare = son pencere sonu (1300–1521).

Ölçülen hızlar (4060 Ti; 6 arka-uç işçisi, OMP=4; **makine başka bir ajanın 15+ `turb.uret_derinlik` süreci ve AirSim
simülatörüyle paylaşımlıydı, yük ortalaması 20–88** → süreler gürültülü ve şişkin):

| bileşen | s/kare (medyan) | 1350 karelik hücrede |
|---|---|---|
| ön-uç D/r₀=0 | 0.025–0.042 | 35–58 s |
| ön-uç D/r₀=8 / 16 | 0.075–0.099 / 0.133–0.144 | 100–140 s / 180–200 s (türbülans CPU'da: faz ekranı + `cv2.filter2D` 31 px çekirdek) |
| A0, A0C, A1, AG, A0C+A3 | 0.07–0.40 (hafif yükte 0.07–0.09) | 1.6–9 dk |
| **A3** (yama + DIMM) | **0.96–1.64** | **16–50 dk** (A0'ın 8–20 katı) |
| A5 (D/r₀ > 0) | 1.06–1.21 | 24–31 dk |

Toplam:

| senaryo | ön-uç (tek akış) | arka-uç işçi-saati | 6 işçi | 10 işçi |
|---|---|---|---|---|
| ölçülen (yüklü makine) | 12.8 sa | 286 | **~48 sa** | **~29 sa** |
| hafif yük (V1_01 temiz hızları) | ~10 sa | ~170 | ~28 sa | ~17 sa |

- Darboğaz **arka-uç CPU'su**, özellikle A3 (arka-uç süresinin ≈ %55–60'ı). Ön-uç zaten CPU ağırlıklı ve arka-uçla örtüşüyor.
- **Yalnız GPU1 (4060 Ti) ile iki GPU arasında anlamlı fark yok:** ikinci GPU yalnız ön-ucu hızlandırır, o da darboğaz değil.
  4060 Ti'da 2 eşzamanlı ön-uç (≈ 1 GB VRAM/süreç) yeterli. Bir arka-uç süreci ≈ 1.3–2 çekirdek kullanıyor; 24 çekirdekte
  10–12 işçi mümkün, ama makine başka işlerle paylaşılırsa kazanç kaybolur.
- `cov_every` 10 → 2 yapılırsa ek maliyet bu ön-testte yük gürültüsü yüzünden ölçülemedi (tahmin ≤ %20); temiz bir makinede
  tek hücreyle ölçülmeli.
- **Disk:** kayıt ≈ 3.9 kB → 5208 kayıt ≈ 20 MB; geçici döküm ≈ 0.13–0.19 GB/hücre, spool'da en çok ≈ 2 yığın × 24 hücre ≈ 9 GB
  (157 GB boş); stdout günlüğü hücre başına ≈ 1–2 MB → 0.3–0.6 GB. Disk sorun değil.

### 4.9 Ön-test tablosu (beklenen ↔ ölçülen)

| # | denetim | beklenen | ölçülen | sonuç |
|---|---|---|---|---|
| T1 | pencere seçimi | 27 pencere ≥ 5 m, örtüşmesiz | 27/27; MH01/100 ve MH02/100 kapsam < 1 m (sallama) | ✅ (karar gerekiyor) |
| T2 | regresyon kapısı | metrik alanlarında 0 fark | 180 alanda 0 fark | ✅ |
| T3 | fizik (fiziksel PSF) | 🛑 yok; EE50 16–24'te 0.8–1.2 | 🛑 yok (S7b ⚠️); EE50 0.95/0.92/0.88 | ✅ (D/r₀ 20–24'te bulanıklık %8–12 eksik) |
| T4 | preflight P-maddeleri | P1/P2/P4/P7/P8 ✅ | P2 🛑 (≈ 11 örnek/pencere); P11/P12 🛑 eski denetim, PyEco yolu ayrıca doğrulandı | ❌ P2 |
| T5 | A0C k | k ∈ {8,12,16}; ölçekleme ±%30 | **k = 12** (2.54; havuz 2.07); ölçekleme 1.4–3.8× sapıyor | ✅ k / ❌ ölçekleme eşdeğerliği |
| T6 | çökme | 126/126 ok | 126/126 ok (+12/12 hat testi) | ✅ |
| T7 | alanlar dolu | 0 eksik | 0 eksik (A0 arm_stats=None tasarım gereği) | ✅ |
| T8 | D/r₀=0 A1 == A0 | bit-aynı | 6/6 bit-aynı | ✅ |
| T9 | D/r₀=0 A5 ≈ A0 | ATE < %10 | 6/6 bit-aynı | ✅ |
| T10 | türbülans uygulandı | PSNR(16) < PSNR(8) − 1 dB | 20.8→17.8, 22.9→20.7 dB; σ oranı 1.78 | ✅ |
| T11 | temiz ATE | 0.5–6 cm, ≤ 0.10 m | 0.6–2.9 cm | ✅ |
| T12 | D/r₀=8 oranı (gözlem) | 0.8–1.5 | 2.2–18.5 | ❌ beklenti (bilimsel sonuç; H1 riski) |
| T13 | NEES yönleri ve örnek sayısı | A0 ≫ 3; A3, AG, A1 < A0; n 15–40 | A0 101–1514; A3 18/18, A1 12/12, AG 17/18; A0C+A3 ≤ A0C 13/18; **n 5–34 (medyan 11)** | kısmen ❌ |
| T14 | süre | ön-uç 60–150 s; A3 1.2–2× A0; A5 2–5× A0 | ön-uç 35–200 s; **A3 8–20× A0**; A5 ≈ 10× A0 | ❌ (A3 maliyeti) |

### 4.10 Engelleyici durumlar ve öneri

§3'teki engelleyici ölçütlerin **hiçbiri** tetiklenmedi (çökme yok, boş alan yok, T8 ve T2 bit-aynı, türbülans uygulanıyor,
temiz ATE ≤ 0.10 m, k ızgara sınırında değil). Mekanik olarak hat uzun koşuma hazır. Ancak uzun koşumdan önce
kullanıcının karar vermesi gereken dört nokta var:

1. **NEES örnek sayısı (P2):** pencere başına ≈ 10 örnek; ANEES bandı χ²(3n)/n bu n ile geniş. Öneri: `cov_every=2`
   (≈ 40–80 örnek/pencere). Değişirse T5 kalibrasyonu (≈ 1 sa) yinelenmeli. (`cells_FAZ2` şu an `cov_every=10`.)
2. **Pencere listesi:** MH01/100 ve MH02/100 sallama pencereleri; `--min_kapsam 1.0` ile değiştirmek ya da ön-kayıt metnine
   uymak. Ön-kayıt dondurulmadan önce seçilmeli.
3. **Süre:** tam ızgara ≈ 17–48 sa (A3 baskın). Makinenin başka işlerle paylaşılmadığı bir zamanda, 10 arka-uç işçisi ve
   2 ön-uç ile koşulması önerilir (≈ 1–1.5 gün). İkinci GPU'yu kullanmak süreyi kısaltmaz.
4. **H1 ve ön-kayıt bütünlüğü:** duman testi D/r₀=8'de 2–18× ATE artışı gösteriyor (V2_01'de D/r₀=4'te bile 1.2–3.5×).
   H1 değiştirilecekse bu artık sonuç görüldükten sonradır; öneri H1'i değiştirmeden dondurmak ve bu gözlemi §0'a pilot
   olgu olarak yazmak.

Ayrıca (engelleyici değil): preflight P11/P12'nin PyEco'ya taşınması; D/r₀ 20–24'te PSF çekirdeğinin (31 px) büyütülmesi
ya da sınırlılık olarak yazılması; GPU eşlemesinin (4060 Ti = `PCI_BUS_ID` 0) uzun-koşum betiğinde düzeltilmesi.

**Öneri: KOŞULLU GO.** Hat mekanik olarak hazır. Uzun koşum, (1) `cov_every` ve (2) pencere listesi kararları verildikten,
ON_KAYIT dondurulduktan sonra (gerekirse T5 yinelenerek) ve makine boşken başlatılmalı. Ön-test sırasında uzun koşum
BAŞLATILMADI.

### 4.11 Ek (2026-11-08): GPU belirlenimciliği ve temel çizgi tasarımı (koordinatör notu üzerine)

**Bağlam (AirSim bulgusu):** PyEco aynı GPU'da bit-belirlenimci; ön-uç (ALIKED/LightGlue) RTX 4060 Ti ile RTX 3090 arasında
kayan nokta düzeyinde farklı. AF_F15_h5'te bu, temiz ATE'yi 0.108 → 0.159 m (%47) oynattı. Arka-uç aynı dökümle bit-aynı.

**Bu ön-testlerde kullanılan GPU:** bütün ön-uçlar (regresyon, kalibrasyon, duman, hat testleri) **tek kartta, RTX 4060 Ti**
(`CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=0`) koştu. Sürücü kayıtlarında `fe_gpu` = "NVIDIA GeForce RTX 4060 Ti"
(duman 132/132, kalibrasyon 96/96). Farklı GPU'dan gelen bir temiz/türbülanslı karşılaştırması yok.

**Runner değişikliği:** `turb/runner.py` artık her kayda `fe_gpu` yazıyor (`_gpu_name()`; `_frontend_worker` meta'sı ve
`run_cell`). Önbellekten okunan dökümde (`cache_dump=True`) alan `None` kalır; FAZ2 hücreleri `cache_dump=False`.
Yeni code_fp `311d74ec9e42` (duman kayıtları `65dbe4761b2b`). Değişiklik yalnız bir alan ekliyor; aşağıdaki tekrar testinde
metrik alanları eski kodla bit-aynı.

**EuRoC'ta tekrar-koşum belirlenimciliği (aynı GPU, 4060 Ti) — GEÇTİ:**

| test | ne tekrarlandı | sonuç |
|---|---|---|
| MH01 D/r₀=0, 2000 kare (T2) | ön-uç + arka-uç, `run_cell` | 180 alan bit-aynı |
| V2_01 D/r₀=0, 1300 kare | ön-uç + arka-uç; sürücü ile uzun-koşum yolu karşılaştırması | A0, A0C 6/6 kayıt bit-aynı |
| **V2_01 D/r₀=4 (türbülanslı), 1300 kare** | `run_experiment_pipelined` iki kez (eski ve yeni runner) | A0, A0C 6/6 kayıt bit-aynı; ön-uç istatistikleri aynı; PSNR aynı (27.119) |

(`evidence/faz2_on_test/hat_testi_1.jsonl`, `hat_testi_2.jsonl`.) GPU'lar arası belirlenimcilik EuRoC'ta sınanmadı:
3090 simülatörle meşguldü. AirSim bulgusu nedeniyle yok sayılmamalı.

**Kurallar (uzun koşum için zorunlu):**
1. Bir karşılaştırma ailesinin (aynı sekans × bütün D/r₀ × bütün tohumlar × bütün kollar) bütün hücreleri **tek, sabit bir
   GPU'da** koşulur. Her kayıtta `fe_gpu` bulunur ve analiz betiği aile içinde tek GPU olduğunu denetler.
2. Türbülanslı bir koşum, başka GPU'da üretilmiş bir temiz koşumla **asla** karşılaştırılmaz. İki GPU kullanılacaksa
   bölünme sekans düzeyinde yapılır (ör. MH'ler 4060 Ti'da, V'ler 3090'da), D/r₀ ya da tohum düzeyinde yapılmaz.

**Tek temiz koşum geçerli bir temel çizgi DEĞİL (EuRoC'ta da ATE kaotik):**
- Kanıt (bu ön-testten): aynı dökümde yalnız gürültü ağırlığı değişince (A0C, k=2; temiz V-serisi, 12 pencere) pencere ATE'si
  A0'ın **0.63–1.18 katı** oluyor (k=4: 0.51–1.23; k=8: 0.55–1.58). Medyan ≈ 1. Yani küçük bir pertürbasyon pencere ATE'sini
  ±%40–50 oynatıyor. AirSim'deki GPU kaynaklı %47'lik oynamayla aynı mertebe.
- Bu nedenle ON_KAYIT §1'deki "D/r₀=0 tek tohum (türbülans kapalı → hat belirlenimci)" tasarımı **yetersiz**. Belirlenimci
  olması, tek koşumun dağılımı temsil ettiği anlamına gelmez. T12'deki "D/r₀=8 / D/r₀=0 oranı 2.2–18.5" tek temiz koşuma
  bölündüğü için kısmen bu kaotiklikten şişmiş olabilir. Özellikle temiz ATE'si 5–6 mm olan MH02/100'deki 18.5× böyle.
- **Gereken tasarım:** tohumlu bir temiz dağılım + eşli karşılaştırma.
  - Her (sekans, tohum) için temiz koşum aynı tohumla küçük bir sensör gürültüsü (ör. σ = 1 DN, Gauss, tohumlu) içerir.
    Türbülanslı koşumlar da **aynı gürültü tohumuyla, aynı gürültüyü** taşır. Eşleme (sekans, tohum, pencere) üzerinden yapılır.
  - Böylece D/r₀=0'da da 5 tohum olur (tam ızgaraya +36 hücre, ≈ +%13 süre) ve "kırılma noktası" (§2) tek bir temiz sayıya
    değil, temiz dağılımın medyanına göre tanımlanır.
  - **Durum:** EuRoC ön-uç kancasında (`turb/simulator.py` / `turb/hooks.py`) sensör gürültüsü seçeneği **yok**.
    `uret_derinlik --gurultu` yalnız AirSim yolunda var. Ayrıca runner `enabled = D/r₀ > 0` kuruyor, yani D/r₀=0'da kanca
    hiç çalışmıyor. EuRoC karşılığı yazılmalı:
    - `turbulence.sensor_noise_dn` + tohum; türbülansla bağımsız, D/r₀=0'da da etkin;
    - `_fe_key`'e eklenmeli;
    - gürültü kapalıyken E0/T2 bit-aynılık kapısı yeniden koşulmalı;
    - gürültülü temiz dağılımın yayılımı (5 tohum × 3 pencere) ölçülmeli.
    Bu yapılmadan uzun koşum başlatılmamalı.

**Güncellenmiş öneri (§4.10'un yerine geçer): ŞİMDİLİK NO-GO.** Hat mekanik olarak hazır ve aynı GPU'da belirlenimci.
Ancak temel çizgi tasarımı (tohumlu temiz dağılım + eşli tasarım) uygulanmadan ve ön-teste tabi tutulmadan tam ızgara
anlamlı bir D/r₀=0 referansı üretmez. Sıra:
1. EuRoC sensör-gürültüsü seçeneği + bit-aynılık kapısı;
2. temiz dağılım ön-testi (yayılım ölçümü);
3. ON_KAYIT §1/§2 güncellemesi;
4. `cov_every` ve pencere listesi kararları (§4.10);
5. dondurma;
6. tek GPU'da uzun koşum.

---

## 5. Ek ön-test: tohumlu sensör gürültüsü ve temiz dağılım (2026-11-08)

### 5.1 Uygulama

- `turb/hooks.py::sensor_noise`: `turb/uret_derinlik.py::gurultu_ekle` ile aynı matematik. Gauss σ DN, türbülanstan SONRA,
  restorasyondan ÖNCE uygulanıyor; `rint` + [0, 255] kırpma + uint8; üreteç `default_rng([noise_seed, kare, kamera])`.
  `noise_seed` verilmezse türbülans tohumu kullanılıyor.
- Kanca artık `sensor_noise_dn > 0` iken D/r₀=0'da da kuruluyor (`config.apply_frontend` kapısı ve `hooks.params_from_cfg`).
  Türbülans simülatörü pasif kalıyor (`enabled=False`), yalnız gürültü uygulanıyor.
- `turb/runner.py`: `Cell.sensor_noise_dn` (0) ve `Cell.noise_seed` (None) alanları. Gürültü yalnız açıkken
  `cfg["turbulence"]`'a, `Cell.key()`'e (`/n1_…`) ve döküm önbellek anahtarına giriyor; kapalıyken üçü de eskisiyle aynı.
- Yedekler: `config.py.yedek_faz2`, `turb/hooks.py.yedek_faz2`, `turb/runner.py.yedek_faz2_gpu`.

### 5.2 Beklenen değerler (temiz-dağılım koşumundan ÖNCE yazıldı)

Tasarım: MH02 ve V1_01_easy, seçilmiş 3'er pencere; A0; σ = 1 DN; tohum 0–4 (gürültü tohumu = türbülans tohumu).
D/r₀=0 ve D/r₀=8, aynı 5 tohumla, eşli. Tek GPU (4060 Ti), uzun-koşum yolu (`run_experiment_pipelined`), 6 arka-uç işçisi.

| # | denetim | beklenen | ölçüt |
|---|---|---|---|
| N1 | gürültü = 0 regresyonu | kanca çıktısı baytları (MH02, 30 kare × 2 kamera, D/r₀ 0 ve 8) ve hat testi metrik alanları (V2_01, D/r₀ 0 ve 4, A0/A0C) değişiklik öncesiyle **bit-aynı** | fark 0 |
| N2 | gürültü istatistiği | std ≈ 1.04 DN (σ=1 + yuvarlama √(1+1/12)); ortalama ≈ 0; tohumlar arası, kameralar arası ve kareler arası korelasyon < 0.02; aynı tohum aynı baytlar; `uret_derinlik.gurultu_ekle` ile bayt-aynı | hepsi |
| N3 | temiz yayılım (pencere başına, 5 tohum) | ATE değişim katsayısı CV = std/ortalama **%5–30** (medyan ≈ %15). Gerekçe: k=2 pertürbasyonu pencere ATE'sini 0.63–1.18× oynattı. 1 DN de benzer büyüklükte bir pertürbasyon | bilgi; CV > %50 ise temiz referans tek koşumla hiç kurulamaz |
| N4 | gürültünün yanlılığı | 5 tohumun temiz ATE medyanı, gürültüsüz temiz ATE'nin 0.7–1.4 katı (1 DN görüntü SNR'sine göre çok küçük, sistematik etki beklenmez) | bilgi |
| N5 | eşli D/r₀=8 oranı | ATE(8, s) / ATE(0, s), 2 sekans × 3 pencere × 5 tohum = 30 çift: pencere medyanları ≥ 2 ve çiftlerin **≥ %80'inde oran > 1.5** (duman testindeki 2.2–18.5 kısmen tek temiz koşumdan şişmiş olabilir, ama 2–18× tek başına yayılımla açıklanamaz). Oranın 1.5'in altına inen çiftleri varsa en çok MH02/600'de (en küçük duman oranı 3.6) ve V1_01/100'de (2.2) | bilgi (H1'e ön-bakış; ön-kayıt değişmez) |
| N6 | türbülanslı yayılım | D/r₀=8'de tohumlar arası ATE CV, temizinkinden büyük (%20–60): türbülans gerçeklemesi de değişiyor | bilgi |
| N7 | alanlar ve GPU | 20 hücre × 3 pencere = 60 kayıt `ok`; her kayıtta `fe_gpu` = "NVIDIA GeForce RTX 4060 Ti"; kayıt anahtarında `/n1`; D/r₀=0 gürültülü kayıtlarda `quality` = None, D/r₀=8'de PSNR gürültüsüz değere (MH02 20.78, V1_01 22.92 dB) yakın (±0.3 dB) | hepsi |
| N8 | süre | hücre başına ön-uç D/r₀=0 40–70 s, D/r₀=8 100–150 s; A0 arka-uç 2–8 dk (makine yüklü) | ölçülür |

### 5.3 Sonuçlar (koşumdan SONRA eklendi; §5.2 değiştirilmedi)

Not: N1'in bayt kısmı ve N2 birim testi §5.2 tablosu yazılmadan hemen önce koşuldu. O sırada referans olarak teorik değerler
(σ√(1+1/12), bağımsızlık) ve koordinatörün AirSim değerleri (std 1.02, korelasyon 0.008/0.003) kullanıldı. N3–N8 tablo
yazıldıktan sonra koşuldu.

**N1 — gürültü = 0 regresyonu: GEÇTİ.**
- Kanca çıktı baytları (MH02 kare 600–629, cam0+cam1): değişiklik öncesi ve sonrası SHA-1 aynı. D/r₀=0'da `63911aeb…` (kanca
  kurulmuyor, eskisi gibi), D/r₀=8'de `70a29d19…`.
- Uzun-koşum yolu hat testi (V2_01, D/r₀ 0 ve 4, A0/A0C, 12 kayıt): 528 metrik alanında 0 fark; anahtarlar aynı; ön-uç istatistikleri
  aynı. Kayıtlara yalnız iki yeni hücre alanı (`sensor_noise_dn`, `noise_seed`) eklendi.
- code_fp `311d74ec9e42` → `08e28b83407c`.
- Dosyalar: `evidence/faz2_on_test/bayt_once.json`, `bayt_sonra.json`, `hat_testi_3.jsonl`.

**N2 — gürültü istatistiği: GEÇTİ.** std 1.040 DN (kırpılmamış pikseller); ortalama 0.0002; tohum korelasyonu 0.009;
cam0–cam1 korelasyonu 0.001; ardışık kare korelasyonu 0.013; aynı tohumda baytlar aynı (iki ayrı koşum, aynı SHA-1);
farklı tohumda farklı; `uret_derinlik.gurultu_ekle` ile bayt-aynı; σ=0'da girdi nesnesi aynen dönüyor
(`gurultu_birim.py`, `bayt_gurultu1_*.json`).

**N3–N7 — temiz dağılım ve eşli D/r₀=8** (`temiz_dagilim.jsonl`, 60 kayıt, hepsi `ok`; hepsinde
`fe_gpu` = RTX 4060 Ti; anahtarlarda `/n1`):

| sekans/pencere | temiz ATE, tohum 0–4 [mm] | temiz CV | gürültüsüz temiz [mm] | temiz medyan / gürültüsüz | D/r₀=8 ATE, tohum 0–4 [mm] | D/r₀=8 CV | eşli oran 8/0 (tohum 0–4) |
|---|---|---|---|---|---|---|---|
| MH02/100 (sallama) | 6.2, 6.1, 6.2, 6.8, 7.8 | 11 % | 5.6 | 1.10 | 87, 82, 121, 107, 147 | 24 % | 14.1, 13.4, 19.5, 15.8, 18.7 |
| MH02/600 | 15.9, 11.6, 20.9, 13.8, 20.1 | **24 %** | 12.4 | 1.29 | 57, 73, 78, 50, 49 | 22 % | 3.6, 6.3, 3.7, 3.6, 2.4 |
| MH02/1000 | 11.6, 11.9, 13.3, 13.5, 10.6 | 10 % | 12.3 | 0.97 | 75, 70, 101, 246, 98 | 62 % | 6.5, 5.8, 7.6, 18.3, 9.3 |
| V1_01/100 | 28.4, 29.8, 28.5, 30.3, 29.6 | 3 % | 28.8 | 1.03 | 58, 112, 79, 76, 56 | 29 % | 2.1, 3.7, 2.8, 2.5, 1.9 |
| V1_01/500 | 29.3, 26.6, 27.1, 26.6, 29.4 | 5 % | 27.1 | 1.00 | 112, 169, 129, 201, 122 | 25 % | 3.8, 6.4, 4.8, 7.5, 4.1 |
| V1_01/900 | 18.2, 16.2, 17.1, 19.8, 19.1 | 8 % | 18.5 | 0.98 | 111, 90, 66, 216, 71 | 55 % | 6.1, 5.6, 3.9, 10.9, 3.7 |

- **N3 temiz CV:** pencere başına %3–24, medyan %9 → beklenen bant (%5–30) içinde; V1_01/100 alt sınırın altında (%3).
  Hiçbir pencere %50'yi geçmiyor. 1 DN, ağırlık değişikliğinden (k=2: 0.63–1.18×) daha küçük bir pertürbasyon. Yine de
  MH02/600'de tek temiz koşum dağılımın ±%25'i kadar yanılabilir.
- **N4 yanlılık:** temiz medyan / gürültüsüz ATE 0.97–1.29 → bant (0.7–1.4) içinde; sistematik etki yok
  (MH02/600'deki 1.29, o penceredeki büyük yayılımın içinde).
- **N5 eşli oran:** 30/30 çiftte oran > 1.5 (%100; beklenen ≥ %80). Medyan 5.7, en küçük 1.9, p10 2.5, en büyük 19.5.
  Pencere medyanları: MH02 15.8 / 3.6 / 7.6, V1_01 2.5 / 4.8 / 5.6. Hepsi ≥ 2 → **beklenti tuttu**.
  Temiz medyana göre oran (medyan(8)/medyan(0)): 17.3, 3.6, 8.2, 2.6, 4.8, 5.0.
  En küçük oran beklendiği gibi V1_01/100'de. MH02/100'deki büyük oran sallama penceresinin küçük temiz ATE'sinden geliyor
  (§4.1 kapsam kararıyla ilgili).
  → D/r₀=8'deki bozulma tek temiz koşumun yapaylığı değil. Duman testindeki 2–18× bulgusu dağılımla da korunuyor.
- **N6 türbülanslı yayılım:** D/r₀=8 CV %22–62 (medyan %27), temizden büyük → beklenti tuttu. Tohum 3, iki pencerede
  (MH02/1000, V1_01/900) uç değer (≈ 2× medyan) veriyor → türbülans gerçeklemesinin kendisi ağır kuyruklu.
- **N7:** 60/60 `ok`, eksik alan yok, `fe_gpu` tek kart. D/r₀=0 gürültülü kayıtlarda `quality` = None. D/r₀=8 tohum 0 PSNR:
  MH02 20.77, V1_01 22.91 dB (gürültüsüz 20.78 / 22.92) → ±0.3 dB içinde. Tohumlar arası PSNR 20.8–21.8 / 22.9–23.6
  (türbülans gerçeklemesi değişiyor).
- **N8 süre (hücre başına, 4060 Ti, 6 işçi, makine AirSim işleriyle yüklü):**

  | | ön-uç | A0 arka-uç |
  |---|---|---|
  | V1_01 D/r₀=0 | 51 s (39–56) | 110 s (109–112) |
  | V1_01 D/r₀=8 | 86 s (83–87) | 71 s (64–83) |
  | MH02 D/r₀=0 | 60 s (56–67) | **429 s (366–442)** |
  | MH02 D/r₀=8 | 113 s (107–136) | 126 s (106–133) |

  20 hücre toplam 33 dk sürdü (1 GPU + 6 işçi). Gözlem: temiz MH02 arka-ucu türbülanslıdan 3.4× yavaş (duman testinde de
  394–776 s) — temiz veride izler uzun, pencere büyük. §4.8 tahmininde temiz hücreler hafif yük hızıyla hesaplanmıştı.
  Temiz hücre sayısı 9'dan 45'e çıkınca (5 tohum) bu hücreler için süre ≈ 5 kat artar; toplam ≈ +%10–15.

**Tasarıma sonuçlar:**
1. Seçenek çalışıyor ve gürültü kapalıyken bit-aynı. Uzun koşumda D/r₀=0 dahil **her hücre σ = 1 DN, gürültü tohumu = hücre
   tohumu (0–4)** ile koşulmalı. Eşleme (sekans, tohum, pencere) üzerinden yapılmalı ve temiz referans 5 tohumun dağılımı olmalı.
2. Temiz yayılım küçük ama sıfır değil (medyan %9, en çok %24). "Kırılma noktası" ve H1 oranları temiz **medyana** göre
   tanımlanmalı; ON_KAYIT §2'deki "D/r₀=0 medyanı" ifadesi bu tasarımla tutarlı, §1'deki "D/r₀=0 tek tohum" ifadesi
   değişmeli (ön-kayıt güncellemesi kullanıcıya/koordinatöre bırakıldı, dondurulmadı).
3. Eşli D/r₀=8 oranı 30/30'da > 1.5 → H1'in ilk yarısı tam ızgarada büyük olasılıkla çürür. Bu bir ön-bakıştır; H1 metni
   değiştirilmemeli (§4.11).

Uzun koşum başlatılmadı; ON_KAYIT dondurulmadı.

---

## 6. Uzun koşum hazırlığı — RTX 3090 (2026-11-08; kararlar kullanıcı tarafından onaylandı)

Onaylanan kararlar: (1) her hücre σ=1 DN, gürültü tohumu = hücre tohumu, tohum 0–4 (D/r₀=0 dahil), eşleme (sekans, tohum,
pencere); (2) `cov_every=2`; (3) `evidence/pencereler_faz2_kapsam1m.json`; (4) H1 değişmeden (ON_KAYIT'ı koordinatör düzenliyor);
(5) bütün EuRoC ızgarası tek kartta: RTX 3090 = `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1` (güç sınırı 390 W,
CPU governor `performance`). Ön-testlerde en çok 8 arka-uç işçisi; koordinatörün AirSim ızgarası aynı anda 3090'ı ve CPU'yu
kullanıyor.

### 6.1 Beklenen değerler (koşumlardan ÖNCE yazıldı)

| # | denetim | beklenen | ölçüt |
|---|---|---|---|
| R1 | `cells_FAZ2` güncellemesi | 9 sekans × 7 D/r₀ × 5 tohum = **315 hücre** (eskisi 279); kol listesi değişmez (6 kol + A5 yalnız MH02/MH04) → **1960 arka-uç koşumu**, 5880 pencere kaydı; her hücrede `cov_every=2`, `sensor_noise_dn=1`, anahtarda `/n1`; pencereler MH01 721/1121/1521, MH02 600/1000/1400 (son pencere sonu 1921 / 1800 kare → bu iki sekansın hücreleri ≈ %26–29 daha uzun) | sayılar tutmalı |
| R2 | 3090'da belirlenimlilik (V2_01, D/r₀ 0 ve 4, A0/A0C, σ=1, tohum 0, uzun-koşum yolu, iki kez) | iki koşum **bit-aynı** (12 kayıt, bütün metrik alanları ve ön-uç istatistikleri); her kayıtta `fe_gpu` = "NVIDIA GeForce RTX 3090". Risk: aynı anda başka bir PyEco ön-ucu 3090'ı kullanıyorsa (koordinatörün AirSim ızgarası) README'ye göre bit-aynılık bozulabilir; bu durumda iki koşum sırasındaki GPU süreçleri kaydedilir ve sonuç "paylaşımlı GPU'da belirlenimci değil" olarak raporlanır | bit-aynı |
| R3 | 3090 ↔ 4060 Ti farkı (bilgi) | aynı hücrenin ATE'si kartlar arasında bit-aynı DEĞİL; pencere başına fark mertebesi temiz yayılımla (%3–24) benzer, en çok ±%50 (AirSim'de %47) | bilgi |
| R4 | NEES örnek sayısı (`cov_every=2`) | pencere başına **40–80** (MH02 4060 Ti denemesinde 39–80; V-serisinde anahtar kare yoğunluğu benzer) | ≥ 30 her pencerede, medyan 40–80 |
| R5 | A0C k yeniden kalibrasyonu (V-serisi D/r₀=0, σ=1, tohum 0, cov_every=2, 3090, 12 pencere, gerçek A0C koşumları) | A0 pencere-medyan NEES'lerinin medyanı ≈ 100–150 (eskisi 118.5); seçilen k ∈ {8, 12, 16}, en olası **12**; k ızgara sınırında (1 ya da 24) değil | k iç noktada |
| R6 | `cov_every=2` maliyeti (aynı dökümde A0, cov_every 10 ↔ 2, V-serisi 4 sekans, eşzamanlı) | A0 arka-uç süresi +%0–20 | bilgi (tahmine girer) |
| R7 | tam hücre süresi (3090; MH02 D/r₀=0 ve D/r₀=8, tohum 0, 7 kol, cov_every=2, σ=1) | ön-uç D/r₀=0 40–70 s, D/r₀=8 90–150 s (MH02 artık 1800 kare); A0 türü kollar 2–8 dk (temiz MH02 daha yavaş, ≈ 7–10 dk); A3 ve A5 25–60 dk | ölçülür |
| R8 | tam ızgara süresi | 315 hücre; 8 işçide ≈ 30–50 sa, 12 işçide ≈ 20–35 sa (A3 baskın; makine paylaşımı süreyi uzatır) | ölçümden tahmin |
| R9 | başlatma betiği (`turb/faz2_kos.py`) | `--kuru` kipinde: GPU adı 3090 değilse durur; mevcut kayıtlarda `fe_gpu` ≠ 3090 ya da farklı `code_fp` varsa durur; hücre sayısı 315, yapılacak hücre sayısını yazar; çalıştırılmaz | kuru kip doğru |

### 6.2 Sonuçlar (koşumdan SONRA eklendi; §6.1 değiştirilmedi)

**R1 — `cells_FAZ2`: GEÇTİ.** 315 hücre, 1960 arka-uç koşumu, 5880 pencere kaydı. Kol kümeleri: 245 hücrede 6 kol,
70 hücrede (MH02, MH04) 7 kol. `cov_every=2`, `sensor_noise_dn=1.0`, `noise_seed=None` (= hücre tohumu), anahtar
`…/s{t}/n1/faz2`. Pencereler `pencereler_faz2_kapsam1m.json`'dan geliyor (MH01 721/1121/1521, MH02 600/1000/1400). Kol listesi
değişmedi. Yedek: `turb/experiments.py.yedek_faz2_karar`.

**R2 — 3090'da belirlenimlilik: GEÇTİ.** V2_01, D/r₀ 0 ve 4, A0/A0C, σ=1, tohum 0, `cov_every=2`, uzun-koşum yolu iki kez
(`evidence/faz2_on_test/r3090/hat_a.jsonl`, `hat_b.jsonl`). 12 kayıt, 552 alan, **0 fark**; ön-uç istatistikleri aynı.
Her kayıtta `fe_gpu` = "NVIDIA GeForce RTX 3090".
- Not: iki koşumun ön-uçları sırasında koordinatörün PyEco süreçleri de 3090'daydı (`yuk.log`: 11:23–11:29 ve 11:38'de
  ≈ 2 GB'lık iki ek süreç). Yani bu donanım/sürücüde GPU paylaşımı bit-aynılığı bozmadı.
- Yük ortalaması koşum sırasında 4 → 36 (koşum a) ve 8 → 165 (koşum b).
- Ön-uç süresi yükle çok oynuyor: D/r₀=4 hücresi 357 s / 101 s, D/r₀=0 hücresi 64 s / 47 s.

**R3 — kartlar arası fark:** sınanmadı. 4060 Ti'da aynı yapılandırmayla (σ=1, `cov_every=2`) koşulmuş bir V2_01 kaydı yok.
Kural gereği karşılaştırma da yapılmayacak: bütün ızgara 3090'da koşulacak.

**R4 (ön bakış, hat testi):** V2_01'de pencere başına NEES örneği 25–71 (pencere 100: 25–33; 500: 42–43; 900: 59–71).
Tam sonuç kalibrasyonda (aşağıda).

**R9 — başlatma betiği: GEÇTİ.**
- Dosyalar: `pyeco_turb/turb/faz2_kos.py` ve kabuk sarmalayıcısı `pyeco_turb/faz2_kos.sh` (`kuru` | `baslat [işçi]` | `denetle`).
- `./faz2_kos.sh kuru` çıktısı: tek kart RTX 3090; 315 hücre, 1960 kol koşumu, kalan 315; code_fp `08e28b83407c`;
  config sha256 `92e1dc7f…`; pencere dosyası sha256 `cd788381…`; k 12.0 (kalibrasyondan sonra güncellenecek).
- Kapı denemeleri:
  - `CUDA_VISIBLE_DEVICES=0` (4060 Ti) → G1 "DUR", çıkış kodu 1.
  - 4060 Ti kayıtlarını içeren bir çıktı dosyasıyla → G2 "DUR".
  - `--denetle` eski duman kayıtlarında: "fe_gpu tek kart değil" ve "NEES örnek sayısı düşük" yakalandı;
    D/r₀=0'da A1==A0 6/6.
- Koşum başlatılmadı.

---

## 7. Küçük A0 ızgaraları: H2 faktöriyel ve H6 nedensel müdahale (2026-11-08; hazırlık + ön-test, BAŞLATILMADI)

Kurallar v0.3: σ=1 DN (gürültü tohumu = hücre tohumu), tohum 0–4, `pencereler_faz2_kapsam1m.json`, tek kart RTX 3090,
eşli (sekans, tohum, pencere), yalnız A0, `cov_every=2`. Etiket `faz2`; ana ızgarayla aynı yapılandırmadaki hücrelerin anahtarı
ana ızgaranınkiyle aynı.

### 7.1 Tanımlar (`turb/experiments.py`)

- `cells_FAZ2_faktoriyel(seeds=5, levels=(8,16), seqs=None, kosullar=("tilt","blur"))`:
  - koşullar `yok`=(tilt kapalı, blur none), `tilt`=(açık, none), `blur`=(kapalı, speckle_global), `ikisi`=(açık, speckle_global).
  - Varsayılan yalnız `tilt` ve `blur`: `ikisi` ana ızgaranın A0 hücresiyle aynı anahtar ve yapılandırma; `yok` ise
    (aşağıdaki B1 tutarsa) ana ızgaranın D/r₀=0 hücresiyle özdeş.
  - 9 sekans × 2 D/r₀ × 2 koşul × 5 tohum = **180 hücre** (dört koşulun hepsi: 360).
- `cells_FAZ2_h6(kirilma, seeds=5, seqs=None, varsayilan_dahil=False)`:
  - D/r₀ ∈ {k−4, k, k+4} × replenish_ratio_kp ∈ {0.2, 0.4} × max_keypoints ∈ {256, 512}.
  - (0.2, 256) PyEco varsayılanı = ana ızgara hücresi, varsayılan olarak üretilmez.
  - 3 yeni bütçe × 3 D/r₀ × 9 sekans × 5 tohum = **405 hücre** (varsayılan dahil 540).
  - Etiket `faz2_h6_r{r}_kp{kp}`: `Cell.key` cfg_overrides'ı içermiyor, ayrı etiket çakışmayı önlüyor.
- `butce_overrides(kp, r)`: `config.load()` oranları yükleme anında uyguladığı için yalnız `max_keypoints`'i ezmek T / N_rep /
  stereo_uzak_min'i eski bütçede bırakırdı. Bu yüzden türetilmiş değerler de yazılıyor:
  `target_feature_count=round(0.6875·kp)`, `replenish_threshold=round(r·kp)`, `backend.stereo_uzak_min=round(0.083333·kp)`.
  kp ≤ 128 (kare-pencere eşiği) desteklenmiyor.
- Runner değişmedi (kod parmak izi aynı; ana ızgara sırasında da güvenle eklenebilir).

### 7.2 Beklenen değerler (duman testlerinden ÖNCE yazıldı)

| # | denetim | beklenen | ölçüt |
|---|---|---|---|
| B1 | `yok` koşulu (D/r₀=8 ve 16, kanal yok, σ=1) | kanca çıktısı D/r₀=0 gürültülü çıktıyla **bayt-aynı** (simülatör kanal yokken görüntüye dokunmuyor). Tutarsa `yok` = ana ızgara D/r₀=0 | bayt-aynı |
| B2 | kanallar gerçekten ayrı | `tilt`, `blur`, `ikisi` çıktıları D/r₀=0'dan farklı; ortalama mutlak fark `blur` > `tilt` (blur her pikseli değiştirir, tilt alt-piksel kaydırma); `ikisi` ≥ her biri; D/r₀=16 > D/r₀=8 | hepsi |
| B3 | H2 duman (MH02, tohum 0, D/r₀ 8 ve 16, `tilt` ve `blur`, A0; 4 hücre × 3 pencere) | 12/12 `ok`, `fe_gpu`=3090. `blur` ATE'si `tilt` ATE'sinden büyük (H2'nin öngörüsü yönünde; duman testinde yalnız bilgi). PSNR: `blur` < `tilt` | çökme yok, alanlar dolu |
| B4 | H6 bütçe gerçekten uygulandı (MH02, tohum 0, D/r₀=12, 4 bütçe; A0) | kayıt `cell.cfg_overrides`'ta kp/r/T/N_rep/stereo_uzak_min görünür. Ön-uç istatistiği `avg_kp1_count`: kp=256'da ≤ 256, kp=512'de > 256 (≈ 400–512). `replenish_ratio 0.4` → stereo tetik oranı (`stereo_run_ratio`) 0.2'ye göre artar. (0.2, 256) hücresi `config.load()` varsayılanıyla aynı türetilmiş değerleri verir (T=176, N_rep=51, stereo_uzak_min=21) | hepsi |
| B5 | H6 duman | 4/4 hücre `ok`; kp=512'de ön-uç ve arka-uç daha yavaş (≈ 1.5–2×) | ölçülür |
| B6 | süre | H2: 180 hücre × (ön-uç + A0) ≈ 180 × (2–5 dk ön-uç + 2–8 dk A0) → 8 işçide ≈ 5–10 sa; H6: 405 hücre, kp=512 daha pahalı → ≈ 12–25 sa (ön-uç baskın olabilir; tek GPU akışı) | ölçümden tahmin |

**R5 — A0C k yeniden kalibrasyonu (RTX 3090, `cov_every=2`, σ=1, tohum 0; V-serisi D/r₀=0, 12 pencere; gerçek A0C koşumları): GEÇTİ, k = 12 (değişmedi).**
`evidence/faz2_on_test/r3090/kalibrasyon_3090.txt`.

| k | A0C pencere-medyanlarının medyanı (NEES_pos) | A0/k² | NEES_rot | ATE A0C/A0 (medyan) |
|---|---|---|---|---|
| 1 | 113.65 | 113.65 | 107.5 | 1.000 |
| 2 | 34.37 | 28.41 | 27.8 | 0.994 |
| 4 | 11.00 | 7.10 | 8.2 | 1.003 |
| 8 | 4.40 | 1.78 | 2.4 | 1.041 |
| **12** | **2.38** | 0.79 | 1.1 | 0.977 |
| 16 | 1.80 | 0.44 | 0.8 | 0.972 |
| 24 | 1.02 | 0.20 | 0.5 | 1.119 |

- 2.37'ye en yakın değer k=12'de (2.38). k iç noktada. A0'ın pencere-medyan NEES'lerinin medyanı 113.65, beklenen 100–150
  bandında. `FAZ2_K_A0C = 12.0` olarak kaldı; yorum güncellendi.
- k ≥ 4'te ölçekleme eşdeğerliği yine yok: gerçek A0C NEES'i A0/k²'nin 1.5–5 katı (§4.4 ile aynı bulgu).

**R4 — NEES örnek sayısı (`cov_every=2`): KISMEN.**
- Pencere başına 25–123 örnek, medyan 64 (beklenen medyan 40–80 → tuttu).
- İki erken pencere 30'un altında kaldı: V1_01/100'de 29, V2_01/100'de 25 (beklenen "her pencerede ≥ 30" → tutmadı).
  Bu pencereler başlangıçtan 100 kare sonra başlıyor; pencere-göreli NEES için gereken kovaryans geçmişi o noktada henüz kısa.
- V1_02 ve V2_02'de pencere başına 68–123 örnek (beklenen üst sınır 80'in üstünde; anahtar kare yoğunluğu daha yüksek).
- Aynı dökümde `cov_every=10` ile 5–25 örnek. A0'ın NEES medyanı cov10/cov2 oranı 1.05; ATE bit-aynı
  (kovaryans dökümü kestirimi değiştirmiyor).

**R6 — `cov_every=2`'nin maliyeti: BEKLENTİ TUTMADI (beklenen +%0–20; ölçülen ≈ +%85).**
Aynı dökümde, aynı anda çalışan iki A0 işinin süreç CPU süreleri (`/proc` örneklemesi, `r3090/cpu.log`):

| sekans | A0, cov_every=2 [CPU-s] | A0, cov_every=10 [CPU-s] | oran |
|---|---|---|---|
| V1_01 | 470 | 245 | 1.92 |
| V1_02 | 772 | 682 | 1.13 |
| V2_01 | 634 | 249 | 2.54 |
| V2_02 | 966 | 526 | 1.84 |

Kovaryans hesabı (pencere Schur tümleyeni) A0 türü kolların maliyetini yaklaşık ikiye katlıyor (medyan 1.88×). Mutlak ek yük
kol başına ≈ 100–400 CPU-s. A3/A5'te bu ek yük göreli olarak küçük, çünkü onların maliyeti yama faktörlerinden geliyor.
Duvar süreleri yük yüzünden kullanılamadı: yük ortalaması 25–225 arasında gezdi; aynı hesaplamanın kopyaları 1.5–3 kat
farklı sürelerde bitti.

**R7 — tam hücre süresi (RTX 3090; MH02, 1800 kare, tohum 0, 7 kol, `cov_every=2`, σ=1; `r3090/kalib3090.jsonl`).**
150 kayıt, hepsi `ok`, `fe_gpu` hepsinde 3090. Yük ortalaması ölçüm boyunca 20–225 (koordinatörün AirSim ve rakip işleri).

| | ön-uç [s] | A0 | A0C | A1 | AG | A0C+A3 | A3 | A5 |
|---|---|---|---|---|---|---|---|---|
| MH02 D/r₀=0, duvar [s] | 122 | 1542 | 1503 | 1450 | 1037 | 1502 | **5798** | 1345 |
| MH02 D/r₀=8, duvar [s] | 203 | 403 | 345 | 673 | 626 | 544 | **3698** | **4030** |

Örneklenen çekirdek oranı kol başına 0.95–1.35 iş parçacığı. Beklentiyle karşılaştırma: ön-uç bandın üst ucunda ya da üstünde
(GPU ve CPU paylaşımlı); A0 türü kollar temiz MH02'de 17–26 dk (beklenen 7–10 dk), D/r₀=8'de 6–11 dk; A3 62–97 dk ve
A5 67 dk (beklenen 25–60 dk) → beklentiden **yavaş**. Nedenleri: `cov_every=2` (R6), 1800 karelik MH02 ve aşırı yük.
Temiz MH02'de A0 türü kollar türbülanslıdan 2–4 kat yavaş (4060 Ti'daki gözlemle aynı).

**R8 — tam ızgara süre tahmini (CPU-saati tabanlı; `r3090/sure_cpu.py`, `sure_tahmini.txt`).**
Kol başına CPU-s/kare MH02 ölçümlerinden alındı (D/r₀=0 ve 8). V-serisi için temiz A0C CPU oranı kullanıldı (0.58 × MH).
D/r₀ ≥ 12 için ×1.5 (kaba varsayım). Hücre karesi = son pencere sonu.

- Arka-uç toplamı ≈ **740 CPU-saat**; A3 + A5 payı %61.
- Ön-uç toplamı (tek akış, en hafif yükteki hızlarla) ≈ 18 sa. Ön-uç darboğaz değil; `gpu_jobs=1` yeterli ve
  bit-tekrarlanabilir (R2).

| kullanılabilir donanım iş parçacığı (12 çekirdek / 24 iş parçacığı) | önerilen arka-uç işçisi | tahmini süre |
|---|---|---|
| 24 (makine yalnız bu ızgaraya) | 16 | ≈ 31 sa |
| 16 | 12 | ≈ 46 sa |
| 12 (AirSim ana ızgara + rakip ızgarası paralel; koordinatörün 14:39 notu) | 8 | ≈ 62 sa |

- CPU-saatleri ağır yük altında ölçüldü. SMT çekişmesi CPU-s'yi şişirdiği için boş makinede gerçek süre %20–30 kısa olabilir.
- D/r₀ ≥ 12 çarpanı belirsiz.
- Gerçekçi aralık **≈ 1.5–3 gün**.
- İlk 24 hücrelik yığın bittiğinde gerçek hızla yeniden tahmin edilmeli (`./faz2_kos.sh denetle` + kayıtlardaki `t_backend`).
- Disk: kayıt ≈ 4 kB × 5880 ≈ 25 MB; spool ≈ 48 hücre × 0.2 GB ≈ 10 GB (`/tmp`'de 136 GB boş).

**R10 — `uret_sekans.py` rektifikasyon farkı (koordinatör notu): ETKİLENMİYOR.**
- `cells_FAZ2`, H2/H6 hücreleri ve k kalibrasyonu `turb/uret_sekans.py` ile üretilmiş, önceden bozulmuş görüntü kullanmıyor.
- Bütün hücrelerin `base_path`'i özgün EuRoC: `<data-disk>/sdata`, dosya tarihleri 15/19 Ekim; kalibrasyon ve hat
  testi kayıtlarında da bu değer.
- Bozulma ön-uç kancasında, PyEco'nun kendi rektifikasyonundan SONRA uygulanıyor
  (`stereo_front_end.py:965`, `rectify_images` → `self.turb.apply_stereo`).
- Rektifikasyon matrisi `config.yaml: frontend.stereo_extrinsics: sensor_yaml`'dan geliyor
  (`run_frontend.py:27` matcher'a geçiriyor).
- `turb/runner.py`, `turb/experiments.py`, `turb/hooks.py`, `config.py`, `run_frontend.py` ve `stereo_front_end.py`'de
  `uret_sekans` referansı yok (grep boş).

### 7.3 Sonuçlar (koşumdan SONRA eklendi; §7.2 değiştirilmedi)

Duman: MH02, tohum 0, RTX 3090, sürücü (`_build_cfg` + `_run_arms`, uzun-koşum yoluyla bit-aynı olduğu R2/§4.7'de gösterildi).
8 hücre × 3 pencere = 24 kayıt, hepsi `ok`, hepsinde `fe_gpu` = 3090 (`evidence/faz2_on_test/r3090/h2h6.jsonl`).

**B1 — `yok` = D/r₀=0: GEÇTİ.** MH02 30 kare × 2 kamera, σ=1, tohum 0: `yok`@D/r₀=8 ve `yok`@D/r₀=16 kanca çıktısı D/r₀=0
gürültülü çıktıyla bayt-aynı (SHA-1 `2b5dec3e…`, fark 0.0 DN; `r3090/fakt_bayt.json`). → `yok` koşulu ana ızgaranın
D/r₀=0 hücreleriyle eşleniyor, ayrıca koşulmasına gerek yok. `ikisi`@D/r₀=8 SHA-1'i (`2f241ec4…`) daha önceki σ=1 ölçümüyle aynı.

**B2 — kanallar: KISMEN.** Kanallar ayrı ayrı uygulanıyor ve D/r₀=16 > D/r₀=8 (tutuyor).
D/r₀=0'a göre ortalama mutlak fark (DN):

| koşul | D/r₀=8 | D/r₀=16 |
|---|---|---|
| tilt | 13.8 | 21.1 |
| blur | 9.0 | 15.0 |
| ikisi | 13.3 | 19.4 |

Beklenen "blur > tilt" ve "ikisi ≥ her biri" **tutmadı**: piksel farkında tilt baskın; bulanıklık tilt'in keskin kenarlarda
yarattığı farkı biraz yumuşattığı için `ikisi` < `tilt`. Bu, piksel ölçüsü hakkında bir yanlış beklentiydi; kod kusuru değil.

**B3 — H2 duman: GEÇTİ (mekanik); bilimsel ön-bakış beklentinin TERSİ.** A0 ATE [m], pencere 600 / 1000 / 1400:

| koşul | D/r₀=8 | D/r₀=16 | PSNR (8 / 16) [dB] |
|---|---|---|---|
| yalnız tilt | 0.023 / 0.055 / 0.136 | 0.320 / 0.245 / 0.568 | 20.0 / 16.4 |
| yalnız blur | 0.016 / 0.014 / 0.031 | 0.110 / 0.054 / 0.174 | 25.1 / 21.4 |
| ikisi (ana ızgara yapılandırması, `kalib3090`) | 0.049 / 0.073 / 0.090 | — | — |
| temiz (σ=1, D/r₀=0) | 0.016 / 0.013 / 0.015 | | |

- Tek hücrede yalnız-tilt bozulması yalnız-blur'dan büyük. D/r₀=8'de yalnız-blur ATE'si temize yakın.
- Bu, H2'nin öngörüsünün ("artışın ≥ %75'i yalnız-blur'da; yalnız-tilt ≤ %25") tersi, çürütme ölçütüne ("yalnız-tilt payı
  > %50") yakın. Tek sekans ve tek tohumluk bir ön-bakış; H2 metni değiştirilmemeli (§4.11 ile aynı ilke).
- PSNR beklentisi "blur < tilt" de tutmadı: tilt, temiz referansa göre PSNR'ı daha çok düşürüyor.

**B4 — H6 bütçesi gerçekten uygulandı: GEÇTİ.**
- Kayıtta `cell.cfg_overrides` görünüyor:
  - kp=512: T=352, N_rep=102 (r=0.2) / 205 (r=0.4), stereo_uzak_min=43.
  - kp=256, r=0.4: N_rep=102.
- (0.2, 256)'nın türetilmiş değerleri `config.load()` varsayılanıyla aynı (176 / 51 / 21).
- Ön-uç `avg_kp1_count`: kp=256'da 201, kp=512'de 354. > 256 ölçütü tuttu; beklenen 400–512'nin altında, çünkü D/r₀=12
  bulanıklığı algılanan nokta sayısını düşürüyor.
- Sağ görüntü nokta sayısı 187 → 356–364.
- `stereo_run_ratio`: r=0.4'te 1.00, r=0.2'de 0.92–0.98; artıyor ama zaten tavana yakın.

**B5 — H6 süre: BEKLENTİDEN YAVAŞ.** MH02 D/r₀=12:
- kp=256 → A0 arka-uç 631–776 s duvar (1120–1312 CPU-s), ön-uç 321–390 s.
- kp=512 → arka-uç 3275–3292 s (4066–4108 CPU-s; **≈ 3.3–5×**, beklenen 1.5–2×), ön-uç 569–782 s (≈ 1.8×).

**B6 — süre tahmini (CPU tabanlı, MH02 ölçümlerinden; 9 sekans ortalaması için kare oranı 0.8 ve V ölçeği 0.58):**

| ızgara | hücre | arka-uç | ön-uç (tek akış) | 12 iş parçacığı + 1 ön-uç | 12 iş parçacığı + 2 ön-uç |
|---|---|---|---|---|---|
| H2 faktöriyel (tilt + blur) | 180 | ≈ 40 CPU-sa | ≈ 10 sa | ≈ 10 sa (ön-uç sınırlı) | ≈ 5–6 sa |
| H6 nedensel (3 yeni bütçe) | 405 | ≈ 230 CPU-sa | ≈ 30 sa | ≈ 30 sa | ≈ 19 sa (arka-uç sınırlı) |

- H6, ana ızgaranın kırılma noktasına bağlı; ancak ana ızgara bitince başlatılabilir.
- `--gpu_jobs 2` iki ön-ucu aynı 3090'da paralel koşar; R2'de GPU paylaşımı bit-aynılığı bozmadı.
- Başlatma (hazır, ÇALIŞTIRILMADI):
  - `./faz2_kos.sh kuru --deney FAKT`
  - `./faz2_kos.sh kuru --deney H6 --kirilma <K>`
  - Çıktılar `evidence/FAZ2_FAKT/`, `evidence/FAZ2_H6/`.
- Tracker: `evidence/FAZ2`, `FAZ2_FAKT`, `FAZ2_H6` ve `faz2_on_test` dizinlerine boş `.trackerignore` kondu. Spool dizinleri
  `/tmp` altında (~/Desktop dışında).

### 7.4 Ön-test tablosu (§6 + §7)

| # | denetim | beklenen | ölçülen | sonuç |
|---|---|---|---|---|
| R1 | `cells_FAZ2` | 315 hücre, 1960 kol koşumu | 315 / 1960 / 5880 kayıt | ✅ |
| R2 | 3090 belirlenimlilik | bit-aynı, fe_gpu=3090 | 552 alanda 0 fark (GPU paylaşımlıyken bile) | ✅ |
| R3 | kartlar arası fark | bilgi | sınanmadı (eşdeğer 4060 Ti kaydı yok; kural gereği gereksiz) | — |
| R4 | NEES örnek/pencere | medyan 40–80, her pencere ≥ 30 | medyan 64 (25–123); 2 erken pencere < 30 | ⚠️ |
| R5 | A0C k | k ∈ {8,12,16}, iç nokta | **k = 12** (2.38) | ✅ |
| R6 | `cov_every=2` maliyeti | +%0–20 | **+%85** (A0 CPU 1.13–2.54×) | ❌ beklenti |
| R7 | tam hücre süresi | A3/A5 25–60 dk | A3 62–97 dk, A5 67 dk (yüklü makine) | ❌ beklenti |
| R8 | tam ızgara | 8 işçi 30–50 sa | ≈ 740 CPU-sa → 24 iş parçacığı ≈ 31 sa, 12 iş parçacığı ≈ 62 sa | ölçüldü |
| R9 | başlatma betiği | kapılar çalışır | G1/G2/G3/G4 + denetim doğru; çalıştırılmadı | ✅ |
| R10 | `uret_sekans` rektifikasyonu | — | etkilenmiyor (yalnız kanca + PyEco sensor_yaml rektifikasyonu) | ✅ |
| B1 | `yok` = D/r₀=0 | bayt-aynı | bayt-aynı | ✅ |
| B2 | kanal farkları | blur > tilt | tilt > blur | ❌ beklenti (kusur değil) |
| B3 | H2 duman | 12/12 ok; blur ATE > tilt ATE | 12/12 ok; **tilt ATE > blur ATE** | ✅ mekanik / ön-bakış H2 tersine |
| B4 | H6 bütçe uygulandı | kp, T, N_rep, kp1 sayısı | 352/102/205/43; kp1 201 → 354 | ✅ |
| B5 | H6 süre | kp=512 1.5–2× | arka-uç 3.3–5×, ön-uç 1.8× | ❌ beklenti |

## 8. Koşum planlaması (2026-11-09 00:35; koordinatör)
- **Bekçi:** `pyeco_turb/sefer_bekci.sh` makineyi boş bırakmamak için bir ızgara bitince sıradakini başlatır
  (günlük `pyeco_turb/evidence/bekci/bekci.log`).
- **Kırılma hesabı:** `turb/kirilma_hesap.py`, ON_KAYIT §2 tanımının birebir uygulaması.
  - Kısmi veride elle sağlandı: MH01 D/r₀=0 medyan 0.0131, D/r₀=4 medyan 0.0412.
  - Karar yalnız ızgara tamamken verilir.
- **H6 spekülatif ön başlatma:**
  - Yalnız SIRALAMA amaçlıdır. MH01'de D/r₀=4 oranı 3.15 olduğu için K=4 (düzeyler 0/4/8) varsayılarak parça 0/2,
    makinede yer açılınca başlar.
  - Kesin K, tam ızgarada §2 tanımıyla hesaplanır.
  - K≠4 çıkarsa yalnız doğru düzeylerdeki hücreler analize girer; ortak hücreler devam kaydından atlanır.
- **Parçalama:** `faz2_kos.py --parca I/N`, sha1(key) mod N ile hücreleri böler. Bölme K'dan bağımsızdır.
  - Hücre içeriği değişmez; code_fp 08e28b83407c aynı.
  - Çıktılar `evidence/FAZ2_H6/records_p0.jsonl` ve `records_p1.jsonl`; analizde birleştirilir.
  - Yedek: `turb/faz2_kos.py.yedek_20261109_parca`.

## 9. H2 analizinin işlemsel tanımı (2026-11-09 08:05; FAKT sonuçlarına BAKILMADAN yazıldı)

FAKT ızgarası 07:56'da bitti: 180/180 hücre, 540 kayıt, hepsi `ok`, 3090, code_fp 08e28b83407c. Denetim çıktısında yalnız
mekanik özet görüldü: kayıt/durum sayıları, NEES örnek sayısı, "ıraksama 5". ATE değerlerine bakılmadı. Görülen tek bilimsel
sayı, B3 dumanıdır (MH02, tohum 0; §7.3). ON_KAYIT §3 H2 metni değişmez; aşağıdakiler yalnız onu sayıya çevirir.

**Koşullar** (§7.1 ve B1):
- `tilt`, `blur`: FAKT ızgarası.
- `ikisi`: ana ızgaranın A0 kaydı (aynı anahtar).
- `yok`: ana ızgaranın D/r₀=0 A0 kaydı (B1: bayt-aynı).

**Birimler ve metrik:**
- Birim: (sekans, tohum), §1'deki gibi. Birim değeri: pencere ATE'lerinin (`m_ate_se3`) ortalaması.
- Birim için dört koşulun hepsi gerekir; eksikse birim dışlanır ve raporlanır.
- `m_valid` = false olan pencere dışlanır ve raporlanır. `m_diverged` olan pencere geçerli bir ATE taşıyorsa içeride kalır;
  büyük hata gerçek bozulmadır.

**Payların hesabı:**
- Artış: Δ_c = ATE_c − ATE_yok, c ∈ {tilt, blur, ikisi}.
- Birincil analiz D/r₀ = 16'da, D/r₀ = 8 ikincildir.
- pay_blur = medyan_u(Δ_blur) / medyan_u(Δ_ikisi); pay_tilt = medyan_u(Δ_tilt) / medyan_u(Δ_ikisi).
  - Medyanların oranı kullanılır, çünkü birim başına oran Δ_ikisi ≈ 0 iken kararsız olur.
  - %95 GA: birimler yeniden örneklenir (10 000 bootstrap, tohum 0).
- Paylar 1'e eklenmek zorunda değildir: etkileşim (alt ya da üst toplamsallık) ayrıca raporlanır,
  pay_blur + pay_tilt − 1.

**Karar** (ON_KAYIT §3):
- pay_blur ≥ 0.75 ve pay_tilt ≤ 0.25 → öngörü tuttu.
- pay_tilt > 0.50 → çürüdü.
- Öteki her durum → kısmen / kararsız; olduğu gibi raporlanır.

**İkincil analiz:**
- Birim başına oranlar ATE_c / ATE_yok: medyan, Wilcoxon tilt ile blur (eşli, iki yönlü).
- Sekans bazında paylar.
- PSNR ve kare başına eşleşme sayısı (betimsel).

**Uygulama ve zamanlama:**
- Betik: `pyeco_turb/turb/h2_analiz.py`.
- `ikisi` ve `yok` ana ızgaradan gelir. Bu yüzden karar, ana ızgarada D/r₀ ∈ {0, 8, 16} A0 hücreleri tamamlanınca verilir.
  O zamana kadar betik yalnız eksik sayısını basar, sayı basmaz.
- Betik sha256: d0b1e386…e353.
  - Öz-test geçti: sentetik payların geri kazanımı, çürütme dalı, kayıt eşlemesi (A0 dışı kayıtlar yok sayılıyor).
  - `--durum` çıktısı: FAKT'ın `tilt`/`blur` hücrelerinde eksik yok. Ana ızgarada D/r₀ 16 için 45 birimden 5'i tam;
    `yok` ve `ikisi` kayıtları bekleniyor.
- Ana ızgara bitince `python -m turb.h2_analiz` koşulacak (koordinatör).
  - Bekçi betiği koşarken düzenlenmedi: bash, çalışan betiği okumaya devam eder.

### 8.1 H6 parça 1 spekülatif başlatıldı (2026-11-10 12:45; koordinatör)

**Ne oldu:**
- H6 parça 0, 2026-11-10 04:52'de bitti: 210/210 hücre, 630 kayıt, hepsi `ok`, 3090, 08e28b83407c.
- Bekçide "parça 0, EuRoC'tan önce biter" olayı yoktu. Bu yüzden makine yaklaşık 8 saat %29 boşta kaldı.

**Kırılma durumu:**
- Ana ızgarada D/r₀ = 0 ve 4, 40/45 sekans-tohum çiftinde tam. Medyan ATE oranı 2.66, kırılma şimdilik 4.
- Bu oran 2'den kesin olarak büyük değil. En kötü durumda, eksik 15 pencere taban için +∞ ve D/r₀ = 4 için 0 alınırsa, oran 1.80'e iner.
- Bu nedenle başlatma, §8'deki gibi yalnız sıralama amaçlı spekülatif bir başlatmadır. Kesin K, ızgara tamamlanınca
  `kirilma_hesap` ile belirlenir.

**Başlatma:**
- Komut: `faz2_kos.sh baslat 8 --deney H6 --kirilma 4 --parca 1/2 --gpu_jobs 2`, nice 5, PID 4118627. Kuru denetim geçti:
  195 hücre, G1–G4.
- `--gpu_jobs 2`: H6 tek kollu (A0) ve ön-uçla sınırlı. R2/§7'ye göre GPU paylaşımı bit-aynılığı bozmuyor.

**Bekçi düzeltmesi:**
- Bekçi durduruldu, yeniden başlatıldı (PID 4119124). Bayrak dosyaları korundu. Yedek: `sefer_bekci.sh.yedek_20261110`.
- Olay 3b eklendi: parça 0 bitmiş ve EuRoC sürüyorsa parça 1 K=4 ile başlatılır.
- K=4 dalı parça 1'i ikinci kez başlatmaz. K≠4 dalı parça 1'i de durdurur, iki parçayı doğru K ile yeniden başlatır.
- EuRoC bitince `turb.h2_analiz` kendiliğinden koşar.

### 8.2 Duraklatma (2026-11-10 13:14; kullanıcının ≈ 2.5 sa hız testi)

- EuRoC ana (3862389) ve H6 parça 1 (4118627) süreç ağaçları SIGSTOP ile duraklatıldı: 29 süreç, PID listesi
  `evidence/bekci/duraklatilan_pidler.txt`.
  - Hiçbir süreç öldürülmedi. Yarım hücreler korunur; runner'da duvar saati zaman aşımı yok.
- Bekçi duraklatma süresince kapalı. Açık kalsaydı EuRoC'u bitmiş sayıp yeniden başlatırdı.
- **Etki:** duraklatma anında süren hücrelerin `t_backend` / `t_frontend` süreleri bu aralığı içerir.
  - Süre analizlerinde bu kayıtlar duraklatma aralığına göre ayıklanacak.
  - Metrikler (ATE, NEES) etkilenmez.
- Devam: aynı PID'lere SIGCONT, ardından bekçi yeniden başlatılır.
- **Devam (2026-11-10 14:56:49):** 29 sürecin 29'u canlıydı; SIGCONT sonrası hepsi R/S durumunda. Bekçi yeniden
  başlatıldı (PID 74314, 14:57:09); bayraklar korunduğu için olaylar kaldığı yerden işlenir.
  - Duraklatma aralığı: 13:14:15–14:56:49 (1 sa 42 dk). Bu aralıkla kesişen kayıtların süre alanları ayıklanacak.
  - Devamda EuRoC 285/315 hücre; H6 parça 1 henüz kayıt yazmamıştı.

### 8.3 EuRoC ana ızgara bitti (2026-11-10 18:32); kırılma, H2 ve H6 parça 1'in yeniden başlatılması

- **Izgara:** 315/315 hücre, 5880 kayıt, hepsi `ok`, hepsi RTX 3090, code_fp `08e28b83407c`; D/r0=0'da A1==A0 135/135
  bit-aynı.
- **Denetimin iki uyarısı incelendi; ikisi de beklenen davranış, yeniden koşum yok:**
  - *Temiz ATE > 0.15 m:* tek pencere, V2_02_medium / tohum 1 / pencere 100, 0.170 m (aynı pencere tohum 2 ve 4'te
    0.112 / 0.109 m). Bu pencere kalibrasyonda da en zor pencereydi (§4, 9.8 cm). Türbülans yok (D/r0=0), A1==A0
    bit-aynı; kestiricinin hızlı V2_02 hareketindeki tohum saçılımı.
  - *Iraksama (max_err > 5 m) 383 kayıt:* D/r0 0/4/8'de bütün kollarda 0; yalnız D/r0 ≥ 12'de var. A0 sayıları
    kirilma_hesap'in ıraksama oranlarıyla birebir uyuşuyor (12: 3/135, 16: 5/135, 20: 19/135, 24: 26/135).
- **Kırılma (ON_KAYIT §2, `turb.kirilma_hesap`):** K = 4 (ATE>2×: medyan 0.0599 / 0.0266 m = 2.25). Spekülatif
  K=4 başlatması doğru çıktı; H6 parça 0 kayıtları geçerli, parça 1 dokunulmadan sürdü.
- **H2 (`turb.h2_analiz`, §9 tanımı):** D/r0=16'da pay_tilt 0.836 [0.64, 1.02] > 0.50 → **ÇÜRÜDÜ**. Ayrıntı
  `evidence/FAZ2_FAKT/h2_analiz.json`; makaleye giden sayılar bu dosyadan alınacak.
- **H6 parça 1 yeniden başlatıldı (18:35):** EuRoC bitince yalnız parça 1'in 8 işçisi kaldı, işlemci ~%40 boştaydı.
  - Bekçi durduruldu; süreç TERM ile, yanıt vermeyen 8 arka-uç işçisi KILL ile kapatıldı; o anda süren hücreler
    (ikinci yığın, kayıt yazılmamış) kaybedildi. 24 tamam hücre (72 kayıt) korunur, devam kaydından atlanır.
  - Yeni ayar: `--jobs 20 --gpu_jobs 3 --chunk 8` (önceki 8 / 2 / 24). Yığın boyu yalnız sırayı değiştirir;
    hücre sonucu hücre içi tohumla belirlenir; GPU paylaşımının bit-aynılığı R2'de gösterildi.
  - Parça 1: 195 hücre (parça 0: 210; toplam 405 = 9 × 3 × 3 × 5), kalan 171. Bekçi yeniden başlatıldı.
- **Özetler (sha256):** `FAZ2/records.jsonl` 511b9a6f…b0fa; `FAZ2/kirilma.json` ec1bc637…c820;
  `FAZ2_FAKT/h2_analiz.json` 9c921449…73af.
  - **Düzeltme (18:39):** `--chunk 8` arka-uç eşzamanlılığını 8'e sınırlıyordu (runner yığın içinde `imap` yapar;
    eşzamanlı hücre ≤ yığın boyu). Yalnız ilk yığının ön-uç işi kaybedilerek `--jobs 20 --gpu_jobs 3 --chunk 40` ile
    yeniden başlatıldı.
- **H6 parça 1 ikiye bölündü (21:07):** `run_experiment_pipelined` her yığının bitmesini bekler. 40'lık ilk yığının
  son 3 hücresi uzarken 20 işçinin 17'si boştu (işlemci %85 boş); swap da dolmuştu.
  - Süreç durduruldu (o an süren 3 hücre ve ikinci yığının ön-uç işi kaybedildi). 61 tamam hücre (183 kayıt) korundu.
  - Kod değiştirilmeden, var olan `--parca` bölmesiyle ikiye ayrıldı: sha1 % 4 ∈ {1, 3} ⇔ sha1 % 2 = 1, yani 1/4 ∪ 3/4 =
    parça 1/2, kesişim boş.
    - `records_p1.jsonl` bölündü: 1/4 hücreleri yerinde kaldı (30 hücre, 90 kayıt), 3/4 hücreleri
      `records_p3.jsonl`'a taşındı (31 hücre, 93 kayıt). Bölme öncesi yedek: `records_p1.jsonl.yedek_bolme_20261110`
      (sha256 77e1f566…).
    - Kuru denetim: 1/4 89 hücre (kalan 59), 3/4 106 hücre (kalan 75); 89 + 106 = 195.
  - Ayar: 1/4 `--jobs 8 --gpu_jobs 2 --chunk 16`; 3/4 `--jobs 10 --gpu_jobs 2 --chunk 20`. İki sürecin yığın sonları
    birbirine denk gelmediği için biri beklerken öteki işlemciyi doldurur.
  - **H6 kayıtları artık üç dosyada:** `records_p0.jsonl` (0/2), `records_p1.jsonl` (1/4), `records_p3.jsonl` (3/4).
    Analiz üçünü birleştirir; beklenen 405 hücre × 3 pencere = 1215 tekil kayıt.
  - Bekçi kapalı bırakıldı: H6 desenleri `--parca 1/2`'ye bağlı, açık kalsaydı H6'yı bitmiş sayardı. Kalan tek olay
    olan H6 bitişini ayrı bir bekleyici izliyor (iki süreç bitince üç parçanın denetimi ve birleşik sayım).

### 8.4 H6 bitti (2026-11-11 01:23) — bütün Faz 2 koşumları tamam

- İki çeyrek süreç 01:00 (1/4) ve 01:20 (3/4) bitti. Birleşik: 405/405 hücre, 1215 kayıt, 1215 tekil
  (key, pencere, kol); hepsi `ok`, hepsi RTX 3090, code_fp `08e28b83407c`, ıraksama 0.
- Denetim uyarısı "temiz ATE > 0.15 m" (3/4 parçası): yine V2_02_medium / pencere 100 (tohum 0, r0.2_kp512: 0.173 m;
  tohum 2, r0.4_kp256: 0.111 m). Ana ızgaradaki §8.3 bulgusuyla aynı pencere; kestiricinin bu penceredeki saçılımı,
  türbülans yok (H6'da D/r0 = K − 4 = 0). Yeniden koşum yok.
- Özetler (sha256): `records_p0.jsonl` dd3c5958…c793; `records_p1.jsonl` 65811a0e…3a97; `records_p3.jsonl` 9b8a2cd8…0a14.
- Bekçi kapalı; kuyrukta koşum kalmadı.

## 10. Ana ızgara + H6 analizinin işlemsel tanımı (2026-11-11, analiz betiğinden önce yazıldı)

**Görülenler beyanı:** koşum sırasında ara tablolar görüldü: kırılma tablosu (§8.3), kısmi kol ANEES oranları
(A0 temiz ≈ 150, A3/AG 0.70–0.95), H6 ham kırılma tablosu (§8.4 sonrası tek bakış). Bu yüzden aşağıdaki tanımlar
ON_KAYIT §2–3'ün metnine harfiyen bağlıdır; serbest seçim yalnız metnin açık bırakmadığı yerlerde ve gerekçesiyle.

- **Veri:** `evidence/FAZ2/records.jsonl` (5880), H6 için `evidence/FAZ2_H6/records_p{0,1,3}.jsonl` (1215) + ana
  ızgaranın A0 hücreleri (varsayılan yapılandırma r0.2_kp256).
- **Metrikler:** ATE = `m_ate_se3` (gerçek zamanlı yörünge, §2 birincil); ANEES_poz = `n_nees_pos_mean`, ANEES_dön =
  `n_nees_rot_mean` (yalnız `n_valid`); bant χ²(3n)/n %95, n = pencere NEES örnek sayısı `n_n`.
- **Birim:** (sekans, tohum), n = 45; birim değeri = pencere değerlerinin ortalaması (ANEES için geometrik ortalama;
  ANEES oranları ölçek-tipi olduğu için). Bir birimin bir kol/düzeyde geçerli penceresi yoksa o birim o
  karşılaştırmadan düşer (sayısı raporlanır).
- **Karşılaştırma:** eşli log-oran log(kol / A0) birim düzeyinde; iki yönlü Wilcoxon işaretli sıra; etki = medyan oran
  (exp medyan log-oran), %95 GA birim bootstrap 10 000 (tohum 0). Holm ailesi = (kollar {A1, A3, AG, A0C, A0C+A3} ×
  D/r₀ {0, 4, …, 24}) aynı metrik için (ATE ailesi, ANEES_poz ailesi). A5 yalnız MH02+MH04 (n = 10), betimsel.
- **H1:** (i) D/r₀ = 8'de A0 ATE oranı 8/0, (sekans, tohum, pencere) eşli (§2 v0.3), medyan + birim bootstrap GA.
  Öngörü ≤ 1.2; çürütme > 1.5. (ii) Kırılma K = `kirilma.json` (4); öngörü 12–20.
- **H3 (EuRoC tarafı, betimsel):** A0 temiz ANEES_poz medyanı ve bant içi/üstü pencere oranı. Karar AirSim'de verildi
  (SONUC_AIRSIM_KOLLAR.md).
- **H4:** (i) D/r₀ ≥ 8'de A3, A1, AG için ANEES_poz log(kol/A0) < 0 ve Holm-anlamlı mı. (ii) A3–AG: D/r₀ ∈ {8, …, 24}
  her düzeyde eşli log(A3/AG), iki yönlü Wilcoxon, ayrı Holm ailesi (5 test); "A3 AG'yi geçer" = medyan oran < 1 ve
  Holm-anlamlı. Bütün kollarda ANEES bandın üstündeyse düşük ANEES = daha tutarlı. Bandın altına inen birim varsa
  ikincil ölçüt |log(ANEES/3)| ile aynı test raporlanır. (iii) "Fark %25'ten küçük" = medyan A3/AG oranı ∈ [0.75, 1.25].
- **H5:** ATE ailesi içinde D/r₀ ≤ 8'de herhangi bir kol için medyan oran < 0.90 ve Holm-anlamlı → çürütme.
- **H6:** her (r, kp) yapılandırması için kırılma = kendi D/r₀ = 0 medyanına göre §2 tanımı (medyan ATE > 2× ya da
  ıraksama > %20), düzeyler {0, 4, 8}; varsayılan yapılandırma ana ızgaranın A0'ı. Kayma = K_yapı − 4. Birincil
  yapılandırma (0.4, 512). 8'de de kırılma yoksa kayma "> 4". Öngörü ≥ 4; çürütme < 2. Ek: D/r₀ = 4'te oran
  (yapı / kendi temizi) birim bootstrap GA'sıyla, betimsel.
- **Çıktı:** `evidence/FAZ2/analiz.json`; betik `turb/faz2_analiz.py` (`--oz_test` sentetik denetim).
