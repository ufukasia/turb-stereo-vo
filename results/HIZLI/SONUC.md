# Hızlı rejim doğrulayıcı ızgarası (Sapma 5) — sonuç özeti

Yazım: 2026-11-12 19:03. Kaynak: `ANALIZ_CIKTI.txt`, `analiz.json` (`hizli_analiz.py`, tam veriyle TEK kez).
Girdi: 1296/1296 kayıt (72 hücre), hepsi ok, boş alan 0; code_fp 08e28b83407c; özetler `analiz_girdi_sha256.txt`.
Birim (sekans, tohum), n = 36 (tohum 1–4); medyan eşli oran [bootstrap %95 GA]; Wilcoxon; Holm aile içinde.
R5 = hızlı rejim (τ₀/Δt = 0.0314), R0 = yarı-statik ana ızgara.

| Hipotez | Öngörü | Ölçülen | Karar |
|---|---|---|---|
| HZ1 | A0 ATE R5@8 / temiz > 1.5 (medyan 8–25); kırılma K = 4 | 14.54 [9.59, 21.09]; K = 4; R5@4 2.96 [2.11, 3.36] | TUTTU |
| HZ2 | A0 ATE R5/R0: D/r₀ 8'de > 1, anlamlı (1.3–2.2); 4'te anlamsız (0.85–1.20) | 8: 1.353 [1.186, 1.532], p_H 1.2e-4; 4: 1.037 [0.863, 1.255], p_H 1 | TUTTU |
| HZ3 (i) | A1, A3, AG konum ANEES'i / A0 anlamlı düşer (iki düzey) | A1 0.165 / 0.078; A3 0.337 / 0.227; AG 0.493 / 0.340 (hepsi p_H ≤ 1.6e-9) | TUTTU |
| HZ3 (ii) | A3/AG ANEES ≥ 1 (A3 üstünlüğünü yitirir) | 4: 0.738 [0.637, 0.824]; 8: 0.671 [0.576, 0.770]; p_H 2.3e-6 | TUTMADI (A3 hızlı rejimde de AG'yi geçiyor) |
| HZ4 | hiçbir kovaryans kolu ATE'yi anlamlı ve > %10 düşürmez | en iyi: A3@8 0.913 [0.872, 0.981], p_H 0.030 (< %10); A1@4 0.898, p_H 0.69 | TUTTU |
| HZ5 | konum ANEES'i R5/R0: 8'de > 1 (A1 1.3–2.3, A0 1.6–3.0), 4'te < 1 | 8: A0 1.665, A1 1.630 (anlamlı), A3 1.159 (anlamsız), AG 1.409; 4: A0 0.810, A1 0.753, A3 0.750, AG 0.738 (anlamlı) | A0/A1 için TUTTU; A3@8 anlamsız |

Okuma:
- Türbülans bozulması gerçekçi enine hızda da sürüyor ve D/r₀ 8'de yarı-statik rejime göre %35 daha kötü
  (ana ızgara D/r₀ 8'de iyimser).
- Kovaryans kolları tek başına doğruluğu iyileştirmiyor (≤ %10); tutarlılığı güçlü biçimde iyileştiriyor.
  Doğruluk kazancı kapıdan geliyor (bkz. pyeco_turb_kapi/evidence/KAPI/SONUC.md).
- Öngörünün tersine fizik-tabanlı kestirim (A3) artık-tabanlı AG'yi hızlı rejimde de geçiyor (0.74 / 0.67).
