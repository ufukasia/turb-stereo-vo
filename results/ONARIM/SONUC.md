# Onarma kolu (Sapma 3) — sonuç özeti

Yazım: 2026-11-12 10:39. Kaynak: `ANALIZ_CIKTI.txt` (`onarim_analiz.py`, tam veriyle TEK kez çalıştırıldı).
Girdi: `records.jsonl` 405 + `ikiz.jsonl` 108 = 513 kayıt, hepsi ok; code_fp 08e28b83407c;
girdi özetleri `analiz_girdi_sha256.txt`. n = 9 sekans; ortanca eşli oran [%95 GA].

| Hipotez | Büyüklük | Ön kayıtlı beklenti | Ölçülen | Sonuç |
|---|---|---|---|---|
| HD1 | R0 ATE O/N (A0) | [0.8, 1.25] | 1.001 [0.708, 1.376] | tuttu; sayısal gürültüden ayırt edilemiyor |
| HD2a | R5 ATE O/N (A0) | ≥ 0.7 | 0.406 [0.218, 0.833], 9/9 < 1 | TUTMADI (onarma öngörülenden çok daha yararlı) |
| HD2b | R5 ATE O/temiz (A0) | ≥ 2 | 4.704 [2.724, 11.949] | tuttu |
| HD3 | ANEES A1/A0 | ≤ 0.33 | R0-O 0.080, R5-O 0.102, R5-C 0.077 | tuttu |
| HD4 | R5 ATE C/N (A0) | ≥ 0.8 | 0.755 [0.482, 1.034] | tutmadı (sınırda; GA 0.8'i içeriyor) |

Sayısal ikiz (N′ = fx, fy × (1+1e-12)): R0 ATE A0 N′/N 0.970 [0.702, 1.028], tipik çarpan ×/÷ 1.264,
birim aralığı [0.525, 1.748]; R5 ATE A0 1.074, ×/÷ 1.189, birim aralığı [0.841, 1.227].
HD2a ve HD4 ortancaları R5 ikiz aralığının dışında; HD1 içinde.

Betimleyici: R5 düzey ortancaları A0 ATE N 0.2543 m → O 0.0814 m; R0 N 0.2260 → O 0.1894 m.
R5 ATE O/N: A1 0.533, A3 0.461. R0 ATE A1/A0 (N) 1.004, A3/A0 0.971. N/temiz ≈ 13.
Dosya yolu ile kanca yolu arasındaki ATE farkı ortanca %15–24 (sayısal duyarlılık düzeyinde).

Makale için okuma:
- Hızlı rejimde (R5) derin onarma ATE'yi ≈ 0.41 katına indiriyor, ama temiz düzeyin hâlâ ≈ 4.7 katı:
  onarma tek başına yetmiyor. Tutarlılık (ANEES) kazancı onarmadan bağımsız sürüyor (HD3).
- Yavaş rejimde (R0) onarmanın ATE etkisi yok (sayısal gürültü içinde).
- HD2a ve HD4 öngörüleri tutmadı; makalede ön kayıtlı beklentiyle birlikte açıkça raporlanacak.
