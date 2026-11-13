# AirSim ana ızgarası — doğrulayıcı sonuçlar (ON_KAYIT_FAZ2 v0.3 §4)

- **Tarih:** 2026-11-08. Izgara 120/120 hücre, hepsi `tamam`, hepsi RTX 3090 (GPU 1), ATE eşleşme dışlaması 0.
- **Analiz:** `airsim_sim/analiz_v03.py` (sha256 45224044…ed98a), girdi `cikti_izgara/ana_v03/w*.jsonl`.
- **Çıktı:** `cikti_izgara/analiz_v03_ana_v03.json` (sha256 0822d780…3137). Sayılar bu dosyadan alındı.
- **Ön kayıt:** `ON_KAYIT_FAZ2.md` v0.3; sapma yok. Sapma 1 yalnız AirSim kollarını ilgilendiriyor.
- **Etki ölçüsü:** eşli oran RPE1s(C0) / RPE1s(temiz), aynı (sekans, pencere, tohum).

## HA3 — F5 (5°): öngörü TUTTU

| C0 | n | medyan oran | boot %95 | GA90 (TOST) | beklenen | p_Holm | >1 |
|---|---|---|---|---|---|---|---|
| 1e-13 | 10 | 1.048 | 1.01–1.08 | 1.005–1.072 (eşdeğer ±%15) | eşdeğer | 0.042 | 8/10 |
| 1e-12 | 10 | 1.129 | 1.10–1.18 | 1.108–1.173 | 1.05–1.30 | 0.0039 | 10/10 |
| 3e-12 | 10 | 1.397 | 1.34–1.44 | 1.332–1.436 | 1.20–1.60 | 0.0039 | 10/10 |
| 1e-11 | 10 | 2.056 | 1.90–2.28 | 1.912–2.239 | 1.60–2.60 | 0.0039 | 10/10 |

- Spearman ρ = 0.957, tek yönlü p < 1e-4 (n = 40).
- Üç bozulma düzeyi Holm sonrası anlamlı ve medyanları öngörülen aralıklarda. 1e-13 ±%15 içinde eşdeğer.
- Not: 1e-13'te küçük (%5) ama istatistiksel olarak saptanabilir bir artış da var (p_Holm 0.042). Eşdeğerlik
  sınırının içinde olduğu için öngörüyle çelişmiyor; makalede ayrıca belirtilecek.

## HA2 — F15 (15°): 3/4 düzey eşdeğer, 1e-11 KARARSIZ (çürütülmedi)

| C0 | n | medyan oran | boot %95 | GA90 (TOST) | eşdeğer | p_Holm | >1 |
|---|---|---|---|---|---|---|---|
| 1e-13 | 14 | 1.010 | 0.93–1.02 | 0.894–1.011 | evet | 0.90 | 8/14 |
| 1e-12 | 14 | 1.010 | 0.88–1.08 | 0.914–1.120 | evet | 0.90 | 7/14 |
| 3e-12 | 14 | 1.025 | 0.96–1.10 | 0.975–1.067 | evet | 0.59 | 9/14 |
| 1e-11 | 14 | 1.098 | 0.97–1.24 | 1.000–1.438 | hayır | 0.18 | 10/14 |

- **Karar:** öngörü 1e-13, 1e-12 ve 3e-12'de tuttu. 1e-11'de GA90'ın üst ucu 1.15'i aşıyor, ama Holm sonrası
  anlamlı bozulma yok (p_Holm 0.18).
  - Çürütme ölçütü (anlamlı bozulma VE medyan > 1.15) sağlanmadı.
  - §4 kuralına göre bu düzey "kararsız" olarak raporlanır. HA2 tam desteklenmedi.
- İkincil bilgi: F15'te Spearman ρ = 0.304, tek yönlü p = 0.011 (n = 56). Bu zayıf bir doz-yanıt eğilimi. 1e-11'de
  medyan ATE oranı 1.56.

## HA5 — Metrik ölçek: TUTTU

- 120 hücrede en büyük |Sim(3) ölçek − 1| = 0.0091 (< 0.02).

## Keşif (hipotez değil)

Benzer σ_tilt'te etki benzer:
- F15 @ 1e-11 (σ_tilt ≈ 0.205 px): medyan oran 1.098.
- F5 @ 1e-12 (σ_tilt ≈ 0.220 px): medyan oran 1.129.

Pilotta görülen fark (6/6 ile 3/6) büyük örneklemde küçüldü. Piksel cinsinden σ_tilt'in etkiyi büyük ölçüde açıkladığı
izlenimi veriyor. Bu keşifsel bir gözlem; σ_tilt eksenindeki iki eğri makalede keşif olarak çizilecek.
