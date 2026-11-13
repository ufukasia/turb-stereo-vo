# Uyarlanır kapı kolu — doğrulayıcı ızgara sonucu (Sapma 6)

Yazım: 2026-11-12 19:02. Kaynak: `ANALIZ_CIKTI.txt`, `dogrula_analiz.json` (`kapi_dogrula_analiz.py` sha256 c7bd583a55d720a2…,
tam veriyle TEK kez). Girdi özetleri `analiz_girdi_sha256.txt`. Kayıt: 1674 yeni + 1188 başvuru = 2862/2862, hepsi ok;
yeni kod parmak izi bedca2849bb2. Doğrulamalar: V1 A0y ≡ A0 378/378 bit-aynı; V2 216; V3 972; V4 702; V5 1674; tamam.
Birim (sekans, tohum), n = 36; medyan eşli oran [bootstrap %95 GA]; iki yönlü Wilcoxon; Holm (m = 3).
KAPI = A0C+A3, cal 1.0, dimm ölçeği 1.0, kapı 3 (çevrimiçi σ̂ ile); A1g = oracle σ ile kapılı A1.

## Birincil hipotezler (ATE, KAPI / A0)
| Hipotez | Koşul | Beklenti | Ölçülen | azalan | p_Holm | Karar |
|---|---|---|---|---|---|---|
| HK1 | R0, D/r₀ 8 | < 1, anlamlı; ≤ 0.85 | 0.549 [0.415, 0.638] | 35/36 | 1.2e-10 | TUTTU |
| HK2 | R5, D/r₀ 8 | < 1, anlamlı; ≤ 0.85 | 0.362 [0.347, 0.529] | 36/36 | 8.7e-11 | TUTTU |
| HK3 | R5, D/r₀ 4 | < 1, anlamlı; ≤ 0.85 | 0.613 [0.564, 0.712] | 33/36 | 4.0e-09 | TUTTU |
| HK4 | temiz (D/r₀ 0) | GA üst sınırı < 1.10 | 1.020 [0.976, 1.047] | 17/36 | — | GEÇTİ (temizde bedel yok) |

## İkincil
- R0, D/r₀ 16: KAPI/A0 ATE 0.794 [0.752, 0.849], 34/36.
- Kapının yalın etkisi (ATE): KAPI/A3S 0.450 / 0.368 / 0.681 / 0.813 (R0-8 / R5-8 / R5-4 / R0-16);
  A1g/A1 0.406 / 0.299 / 0.671 / 0.554. Kapısız kollar A0'dan ayırt edilemiyor (A1/A0 0.90–1.04; A3S/A0 0.94–1.03).
- Çevrimiçi / oracle (KAPI/A1g ATE): 1.244 / 1.236 / 1.046 / 1.423 → oracle σ ile %5–42 daha iyi; çevrimiçi σ̂ oracle'ın
  altında kalıyor (kapı sigması medyanı 2.01 / 2.23 / 1.70 / 2.31 px; oracle 4.52 / 4.52 / 2.67 / 7.93 px).
- Tutarlılık (konum ANEES'i KAPI/A0): 0.156 / 0.075 / 0.401 / 0.239; KAPI/A3S 0.455 / 0.281 / 0.999 / 0.964.
  Mutlak ANEES hâlâ 1'den çok büyük (KAPI medyan 226 / 182 / 66 / 1866): tutarlılık iyileşiyor ama sağlanmıyor.
- Temiz: KAPI/A0 konum ANEES'i 0.683; ATE değişmiyor.
- Birim medyan ATE [m]: R0-8 A0 0.2647 → KAPI 0.1801 (A1g 0.1489); R5-8 0.4673 → 0.2481 (0.1746);
  R5-4 0.0781 → 0.0478 (0.0430); R0-16 0.8181 → 0.6157 (0.4874); temiz 0.0262 → 0.0263.

## Yürütme sapmaları (makalede raporlanacak)
- II. parçada D/r₀ 16 hücreleri D/r₀ 8'den önce koştu (yürütme sırası; içerik ve örneklem aynı).
- I. parça bellek/yük nedeniyle durdurulup aynı komutla sürdürüldü (son yeniden başlatma 2026-11-12 10:28, 8 işçi);
  koşucu devam ettirilebilir, tamamlanan kayıtlar yeniden koşulmadı.
- Ara bakışlar yalnız bütünlük içindi; ATE/ANEES oranlarına bu çözümlemeden önce bakılmadı.
