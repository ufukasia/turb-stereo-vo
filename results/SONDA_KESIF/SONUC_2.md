# Keşif sondası 2 — sonuç (KEŞİF; ön kayıt dışı, makalede doğrulayıcı kanıt DEĞİL)

Yazım: 2026-11-12 11:00. Kaynak: `ANALIZ_2_CIKTI.txt` (`kapi_sonda_analiz.py`, tam veriyle tek kez).
Girdi: `kapi2.jsonl` 87/87 kayıt, hepsi ok, ıraksayan 0; code_fp 08e28b83407c; özet `analiz2_girdi_sha256.txt`.
Koşul: R0, D/r₀ 8, tohum 0, 9 sekans × 3 pencere. Etki = sekanslarda eşli medyan oran [en küçük, en büyük].

| Ön beyan | Büyüklük | Beklenen | Ölçülen | Sonuç |
|---|---|---|---|---|
| B1 | K ≡ ilk sondanın A0g kayıtları | 6/6 bit-aynı | 6/6 | tuttu |
| B3 (Q1) | Gd / A0 ATE | 0.50–0.80, ≥ 7/9 | 0.492 [0.36, 0.89], 9/9 | kural tuttu (beklenen aralığın hemen altında) |
| B3 | Gd / A0g ATE | 1.0–1.5 | 0.988 [0.85, 1.13] | kayıp yok |
| B4 (Q2) | Gn / A0 ATE | 0.55–0.90 | 0.483 [0.37, 1.00], 8/9 | kazancı TAŞIYOR (beklenenden güçlü) |
| B4 | Gp / A0 ATE | 0.85–1.15 | 1.051 [0.77, 1.32], 3/9 | taşımıyor |
| B5 | Gn / A0 konum ANEES'i | < 0.3 | 0.049 | tuttu |
| B5 | Gd / A0 konum ANEES'i | < 0.3 | 0.067 | tuttu |
| B5 | Gp / A0 konum ANEES'i | 0.5–1.5 | 0.886 | tuttu |
| B6 | ıraksayan pencere | ≤ 1/27 | 0/27 (üç varyant) | tuttu |

Okuma:
- Q1: kapı ölçeği oracle yerine kestirim düzeyine (DIMM σ̂) indirildiğinde ATE kazancı korunuyor
  (Gd/A0 0.492; Gd/A0g 0.988). Tutarlılıkta bedel var: Gd/A0g konum ANEES'i 1.869 (0/9 < 1).
- Q2: kazancı gürültü grubu (px_noise + stereo_sigma; χ² kapıları dahil) taşıyor; PnP eşiği tek başına etkisiz
  (Gp/A0 1.051; Gp/A0g ATE 2.252).
- Düzey ortancaları (ATE m / konum ANEES'i): A0 0.1339 / 909.8; A0g 0.1052 / 39.2; Gd 0.0860 / 64.0;
  Gn 0.1116 / 41.7; Gp 0.1483 / 809.5; A1 0.1974 / 90.3.
- Sınır: tek tohum, 9 birim, test ve güven aralığı yok; Gd'de σ̂ nedensel değil (önceki koşumun medyanı);
  0.9–1.2 arasındaki oranlar sayısal gürültüden ayırt edilemez. Doğrulayıcı kanıt yalnız Sapma 6 ızgarasından gelir.
