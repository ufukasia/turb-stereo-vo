# Sapma 2 rakip koşumları — sonuç (HR1, HI1)

Üretim: `rakip_ek_analiz.py` → `evidence/FAZ2_rakip_ek/rakip_ek_analiz.json` (sha256 75b244352b9743da…), 2026-11-11 17:03. Sayılar bu JSON'dan otomatik yazıldı.
Kayıt kimliği: {'ana': ['f578efb6ff85'], 'ek': ['f578efb6ff85'], 'si': ['2e0701051eb2']}. Altyapı hatası 0. Betimsel (Holm ailesi yok; ON_KAYIT §5 Sapma 2).

## HR1 — erken kırılma PyEco'ya özgü mü (D/r₀ 4 / 0, nedensel ATE, ikisinde de başarılı pencereler)

| Rakip | pencere çifti medyanı [%95 GA] (n) | birim, log-oran ort. [GA] | birim, ATE ort. (PyEco H1 ile aynı) [GA] | Wilcoxon p (iki yönlü, birim ATE ort.) | başarısız |
|---|---|---|---|---|---|
| ORB-SLAM3 stereo | 2.23 [1.95, 2.48] (78) | 2.38 [1.87, 3.03] | 2.65 [1.91, 3.15] | 1.5e-08 | 3/81 |
| AirSLAM | 1.66 [1.52, 1.84] (78) | 1.7 [1.46, 2.08] | 1.7 [1.43, 2.12] | 2.8e-07 | 0/81 |
| DROID-SLAM | 1.89 [1.72, 2.40] (81) | 2.12 [1.97, 2.46] | 2.02 [1.81, 2.42] | 1.5e-08 | 0/81 |
| PyEco A0 (ana ızgara, 45 birim) | — | — | 2.95 [2.23, 3.28] | — | — |

- **Ön kayıtlı ölçüt (en az iki rakipte oran ≥ 2):** pencere çifti düzeyinde 1/3 → çürüdü; birim (sekans, tohum) düzeyinde 2/3 → tuttu. Ön kayıt metni eşlemeyi pencere düzeyinde, istatistik birimini (§1) (sekans, tohum) olarak tanımlıyor; hangi toplamanın eşiğe uygulanacağı yazılmamıştı. İki okuma da raporlanır; sonuç görüldükten sonra seçim yapılmaz.
- Üç rakipte de bozulma var ve anlamlı; ORB-SLAM3 açıkça ≥ 2, AirSLAM açıkça < 2, DROID ≈ 2 (GA 2'yi kapsıyor). PyEco'nun oranı dört sistemin en yükseği.
- AirSLAM V2_01 / pencere 100 duyarlılığı etkisiz: o pencere temizde üç tohumda da başarısız, çiftlere girmiyor.

## HI1 — ataletsel destek (ORB-SLAM3 stereo-ataletsel, son yörünge; Sapma 2 eki)

| D/r₀ | SI başarısız | stereo başarısız | McNemar p* | SI oran [GA] (n) | stereo oran [GA] (n) | SI ATE medyan [m] | stereo ATE medyan [m] |
|---|---|---|---|---|---|---|---|
| 0 | 9/81 | 0/81 | 0.0039 | 1 [1.00, 1.00] (72) | 1 [1.00, 1.00] (81) | 0.0155 | 0.0295 |
| 4 | 9/81 | 3/81 | 0.15 | 1.46 [1.33, 1.76] (72) | 2.14 [1.74, 2.63] (78) | 0.0247 | 0.0569 |
| 8 | 12/81 | 5/81 | 0.12 | 2.24 [1.98, 2.73] (69) | 4.93 [3.92, 6.29] (76) | 0.0388 | 0.124 |
| 16 | 31/81 | 36/81 | 0.53 | 6.19 [4.57, 7.56] (49) | 19.4 [12.92, 25.34] (45) | 0.107 | 0.406 |
| 24 | 58/81 | 64/81 | 0.29 | 10.4 [8.38, 12.39] (23) | 40.8 [23.52, 51.94] (17) | 0.253 | 0.38 |

\* Kesin McNemar (uyumsuz pencereler üzerinde binom); sonradan eklenen betimsel sınama, ön kayıtta yok.

- **Karar (ön kayıtlı ölçüt): TUTTU.** D/r₀ 16'da SI başarısızlık oranı 38.3%, stereo 44.4% (SI düşük); SI eşli medyan ATE oranı 16/0 = 6.19 [4.57, 7.56] ≥ 2.
- **Çekince:** başarısızlık farkı küçük ve eşli sınamada anlamlı değil (yalnız SI 18, yalnız stereo 23, ikisi 13 pencere; p = 0.53). Makalede "ataletsel destek ıraksamayı azaltır" denmez; "başarısızlık oranı benzer, doğruluk kaybı belirgin biçimde küçük ama sürüyor" denir.
- **SI'nin temizdeki yapısal başarısızlıkları:** MH03/w162, MH04/w132, MH05/w102 bütün tohum ve tekrarlarda kapsam ölçütünden (IMU ilklemesi bitmeden yörünge başlamıyor) başarısız (9/81). Bu pencereler hariç D/r₀ 16: SI 23/72, stereo 28/72 (p = 0.52).
- Nedensel (gerçek zamanlı) SI metrikleri ikincil olarak JSON'da (`si_nedensel (ikincil)`).
