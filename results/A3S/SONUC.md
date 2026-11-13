# A3S ayırma kolu (Sapma 4) — sonuç özeti

Yazım: 2026-11-12 15:58. Kaynak: `ANALIZ_CIKTI.txt`, `analiz.json` (`a3s_analiz.py`, tam veriyle TEK kez).
Girdi: `records.jsonl` 702 kayıt, 225 hücre, hepsi geçerli, hücre hatası 0; A0 bit-aynı 27/27; code_fp 08e28b83407c.
Birim (sekans, tohum), n = 45; medyan eşli oran [bootstrap %95 GA]; Holm ailesi beş seviye; eşdeğerlik t-TOST ±%15.

| Hipotez | Büyüklük | D/r₀ 8 | 12 | 16 | 20 | 24 | Karar |
|---|---|---|---|---|---|---|---|
| HS1 | ANEES A3S/A3 (eşdeğerlik 8, 12, 16) | 1.032 [0.952, 1.304] | 1.085 [1.000, 1.192] | 1.119 [1.029, 1.181] | 1.035 (eşd.) | 1.212 [1.053, 1.233] | TUTMADI (eşdeğerlik yalnız 20'de) |
| HS2 | ATE A3S/A3 (eşdeğerlik, beş seviye) | 1.018 | 1.021 | 0.976 | 1.034 | 1.012 | TUTTU (5/5 eşdeğer) |
| HS3 | ANEES A3S/AG (< 1, Holm p < 0.05; 8, 12, 16) | 0.895, p_H 0.045 | 0.830, p_H 0.008 | 0.925, p_H 0.355 | 0.883 | 0.918 | TUTMADI (16'da anlamlı değil) |

Yorum kuralına göre sonuç: **kararsız** — seviye seviye GA'larıyla raporlanır.
- A3S/A3 ANEES medyanı hiçbir seviyede > 1.15 eşiğini (8, 12, 16) aşmıyor (1.03 / 1.09 / 1.12; Holm p 0.229 / 0.052 /
  0.039): yama korelasyonunun tutarlılığa katkısı küçük (%3–12), ama eşdeğerlik de gösterilemedi.
- Doğruluk aynı (HS2): yama modeli ATE'yi değiştirmiyor.
- A3S, AG'ye karşı 8 ve 12'de anlamlı biçimde iyi; 16'da yön aynı ama anlamlı değil.
- Betimsel ANEES (kol/A0): A3S 0.351 / 0.301 / 0.268 / 0.247 / 0.177; A3 0.337 / 0.271 / 0.236 / 0.225 / 0.175.
- σ̂ oranı (A3S/A3) medyan 0.998 / 1.001 / 0.999 / 0.993 / 0.946; ±%5 içinde 42 / 37 / 27 / 14 / 4 (45'te).
- Arka-uç süresi bu koşumda makine yüküyle karışık (833–3556 s); süre karşılaştırması makalede KULLANILMAZ
  (eşit yük altında ayrı ölçüm gerekir).
