# Hızlı rejim doğrulayıcı ızgarası — ön test (ON_KAYIT_FAZ2.md §5 Sapma 5)

Yazım: 2026-11-11 15:48, **ön test koşulmadan, keşif sondası tamamlanmadan ve hızlı rejimde (R5) hiçbir sonuç (ATE, ANEES)
görülmeden önce.** Kural: uzun koşum ön-test kuralı (yedi madde). Bu belgenin kapıları geçmeden 72 hücrelik
ızgara başlatılmaz. Ön testin tohum 1 hücreleri doğrulayıcı kümenin parçasıdır; bu yüzden ön test, Sapma 5 metni
(öngörüleriyle) ON_KAYIT §5'e yazılıp özeti alındıktan **sonra** koşulur.

## 1. Kimlik (madde 1)

| Öğe | Değer |
|---|---|
| Kod parmak izi | `08e28b83407c` (ana ızgarayla aynı; parmak izli hiçbir dosya değişmedi) |
| Ön kayıt | `ON_KAYIT_FAZ2.md` — Sapma 5 eklendikten sonraki sha256 `kimlik_sha256.txt` dosyasına ön testten önce yazılır |
| Koşucu / analiz / denetim | `hizli_kos.py` `b27d9664…`, `hizli_analiz.py` `18dfa6a4…`, `hizli_ontest_denetim.py` (tam özetler `kimlik_sha256.txt`) |
| Yapılandırma | `config.yaml` `92e1dc7f…` (ana ızgaradaki dosya) |
| Pencereler | `evidence/pencereler_faz2_kapsam1m.json` `cd788381…` |
| Eşleme kaynağı | `evidence/FAZ2/records.jsonl` `511b9a6f…` (5880 kayıt), `evidence/FAZ2/analiz.json` `7bd332cf…` |
| Gerileme başvurusu | `evidence/SONDA_KESIF/records.jsonl` (V1_01_easy, tohum 0, R5, D/r₀ 8, fe0; kollar A0, A1) |
| Veri | `SLAM_DATA_ROOT=<data-disk>/sdata` (EuRoC; türbülans kanca yoluyla ön-uçta uygulanır, disk verisi üretilmez) |
| Ortam | `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=4` (RTX 3090; koşucu başka kartta durur) |
| Hücre tanımı | `turb.experiments.cells_FAZ2` ile aynı alanlar; fark: `tau0_over_dt` (R5: v = 5 m/s) ve `tag` (`hizli_R5`; ön testte `hizli_ontest_R5`) |

## 2. Ön test hücreleri (madde 3)

V1_01_easy, üç pencere (100, 500, 900); üç hücre, 14 kol koşumu, 42 kayıt:

| Hücre | Tohum | D/r₀ | τ₀/Δt | Kollar | Kayıt | Amaç |
|---|---|---|---|---|---|---|
| gerileme | 0 | 8 | 0.0157 | A0, A1 | 6 | sonda kaydıyla bit-aynılık (madde 2) |
| duman 8 | 1 | 8 | 0.0157 | A0, A0C, A1, A3, AG, A0C+A3 | 18 | bütün kollar uçtan uca (madde 3) |
| duman 4 | 1 | 4 | 0.0314 | A0, A0C, A1, A3, AG, A0C+A3 | 18 | ana etkenin ikinci düzeyi |

"Temiz" (D/r₀ = 0) düzeyi bu ızgarada yoktur (rejim temiz görüntüde tanımsız); temiz başvuru ana ızgaranın kaydıdır.

## 3. Beklenen değerler (madde 2, 4, 5) — sonuçtan önce yazıldı

Ana ızgara başvuruları (V1_01_easy, tohum 1, yarı-statik τ₀/Δt = 1; bilinen veri). Konum ANEES'i üç pencerenin
geometrik ortalaması / ATE üç pencerenin aritmetik ortalaması [m]:

- D/r₀ 0: A0 2139.4 / 0.0246; A0C 20.4 / 0.0240.
- D/r₀ 4: A0 1474.5 / 0.0408; A1 220.3 / 0.0379; A3 620.7 / 0.0363; AG 693.2 / 0.0341; A0C 15.3 / 0.0450; A0C+A3 12.8 / 0.0454.
- D/r₀ 8: A0 847.3 / 0.0823; A1 48.4 / 0.1236; A3 242.7 / 0.1037; AG 384.0 / 0.1193; A0C 10.8 / 0.1204; A0C+A3 10.7 / 0.0900.
- Ön-uç (D/r₀ 4 / 8): zamansal eşleşme 160.53 / 134.59; stereo 4.35 / 11.22; anahtar nokta 253.9 / 235.3.
- Kol istatistikleri (D/r₀ 4 / 8): A1 σ 2.476291 / 4.412249 px; A3 σ̂ medyanı 1.0266 / 1.3143 px; AG artık RMS medyanı
  0.8999 / 1.1319 px; karşılaştırılan poz sayısı 400 / pencere.
- Süre (ana ızgara, yük altında): ön-uç 126–136 s; arka-uç A0, A1, AG, A0C, A0C+A3 ≈ 270–290 s, A3 ≈ 1350–1630 s.

| Kapı | Ölçü | Beklenen | DUR ölçütü |
|---|---|---|---|
| G1 kimlik | kayıtlardaki `code_fp` | hepsi `08e28b83407c` | başka değer |
| G2 gerileme (bit-aynılık) | gerileme hücresi: A0 ve A1, üç pencere; `m_ate_se3`, `n_nees_pos_mean`, `n_nees_rot_mean`, `m_rpe_trans` sonda kaydıyla `==`; ön-uç istatistikleri (süre dışında) aynı | 6 / 6 bit-aynı; ön-uç aynı | bir kayıt bile farklı |
| G3 duman | iki duman hücresi | 36 kayıt; hücre hatası 0; bozuk satır 0; "ok" kayıtlarda dört metrik alanı dolu; "ok" ≥ 30 / 36 | herhangi biri tutmazsa |
| G4a birim | `tau0_over_dt` | 0.0157 (D/r₀ 8), 0.0314 (D/r₀ 4) = 0.314·(D/(D/r₀))/(v·Δt), D = 0.10 m, v = 5 m/s, Δt = 0.05 s | başka değer |
| G4b etki uygulandı (rejim) | ön-uç istatistikleri ana ızgaranın aynı (sekans, tohum, D/r₀) R0 hücresinden farklı; ortalama anahtar nokta sayısı R5 / R0 | farklı; oran 0.90–1.10 (bulanıklık rejimden bağımsız, yalnız tilt'in zamansal korelasyonu değişir) | aynıysa (rejim uygulanmamış) ya da oran bandın dışındaysa |
| G4c kol ayarları | A1'in kullandığı σ; A0C `cal_factor`; A3 yama faktörü sayısı | 4.412249 / 2.476291 px (ana ızgarayla `==`: oracle, rejimden bağımsız); 12.0; > 0 | farklıysa |
| G4d DIMM kestirimi | A3 σ̂ medyanı R5 / R0 (aynı sekans, tohum, D/r₀) | 0.75–1.33 (anlık sol-sağ farkından kestirilir; rejime duyarsız olmalı) | bandın dışında |
| G4e artık kestirimi | AG artık RMS medyanı R5 / R0 | > 1 (öngörü: hızlı tilt artıklarda görünür). Kapı değil; yönü raporlanır | — |
| G4f etki uygulandı (kollar) | konum ANEES'i (pencere geometrik ort.) A1 / A0 ve A0C / A0 | ikisi de < 1, iki düzeyde de | biri ≥ 1 |
| G4f′ | A3 / A0 ve AG / A0 | < 1 beklenir (HZ3 (i)); kapı değil | — |
| G4g ölçek | A0 ATE [m], pencere başına | sonlu ve 0.005–3 m (temiz 0.025 m; yarı-statik D/r₀ 8'de 0.08 m) | aralığın dışında |
| G4h pencere tanımı | karşılaştırılan poz sayısı, R5 / R0 aynı pencere | ±%5 içinde (400) | dışında |
| G5 çıktı ve toplama | `hizli_analiz.py --hizli evidence/HIZLI/on_test.jsonl --tohumlar 1 --sekanslar V1_01_easy --etiket hizli_ontest_R5 --butunluk` | "36/36 kayıt", tam = True, eksik 0; `--butunluk` olmadan "DUR: 1 birim" (test tablosu basılmaz) | tutmazsa |
| G6 süre ve disk | `t_frontend`, kol başına `t_backend`; /tmp boş alanı | bilgi: toplam = 72 × (ön-uç + altı kol); /tmp ≥ 20 GB boş (şimdi 77 GB; hücre dökümü ≈ 0.3–0.7 GB, hücre bitince silinir) | /tmp < 20 GB |

Yön öngörüsü (kapı değil, raporlanır): zamansal eşleşme sayısı R5'te R0'dakinden büyük olmaz; A1'in konum ANEES'i
R5'te R0 değerinden (48.4 / 220.3) küçüktür (HZ5 ile aynı yön).

Tablo `hizli_ontest_denetim.py` ile basılır (betik ön testten önce yazıldı; boş veride bütün kapılar "TUTMADI" veriyor).
Kapılardan biri tutmazsa ızgara başlatılmaz, neden araştırılır. Izgarada V1_01_easy tohum 1 hücreleri yeniden koşulur ve
ön test kayıtlarıyla bit-aynı olmalıdır (`hizli_kos.py` denetimi; tutmazsa koşucu DUR ile çıkar): bu, koşum içi erken
denetimdir.

## 4. Sonuç (ön testten sonra doldurulur)

Ön test 2026-11-11 18:44–19:45 arasında koşuldu (`evidence/HIZLI/on_test_1844.log`; 3 hücre, 14 kol koşumu, 42 kayıt, sorunlu
hücre 0, çıkış kodu 0). Bu bölüm 2026-11-11 19:48'de yazıldı. Kayıt dosyası `evidence/HIZLI/on_test.jsonl` sha256
`fac80184f225362135102930f18512b2eb63eb19494ee9c3d668f05dd41933ee`.

`hizli_ontest_denetim.py` çıktısı (olduğu gibi):

| no | denetim | ölçülen | beklenen | sonuç |
|---|---|---|---|---|
| G2 | tohum 0, R5, D/r0 8, A0 + A1 ↔ sonda | 6/6 kayıt bit-aynı; ön-uç istatistikleri aynı: True | 6/6 bit-aynı; ön-uç aynı | geçti |
| G3 | tohum 1, D/r0 8 ve 4, altı kol | 36 kayıt, durum {'ok': 36}, hücre hatası 0, bozuk satır 0, 'ok' kayıtlarda alanlar dolu: True | 36 kayıt; hücre hatası 0; bozuk 0; alanlar dolu; ok ≥ 30 | geçti |
| G1 | kod parmak izi (kayıtlarda) | {'08e28b83407c'} | {'08e28b83407c'} | geçti |
| G4a | D/r0 8: tau0/dt | {0.0157} | {0.0157} (v = 5 m/s) | geçti |
| G4b | D/r0 8: rejim uygulandı (ön-uç istatistikleri R0'dan farklı; zamansal eşleşme R5 / R0) | aynı: False; avg_temporal_matches 122.0 / 134.6; kp1 239 / 235 | farklı; anahtar nokta sayısı ±%10 içinde | geçti |
| G4c | D/r0 8: A1 oracle σ (kullanılan) | 4.412249159365505 | 4.412249159365505 (ana ızgarayla aynı) | geçti |
| G4c | D/r0 8: A0C cal_factor; A3 yama faktörü | 12.0; 9706 | 12.0; > 0 | geçti |
| G4d | D/r0 8: A3 σ̂ medyanı R5 / R0 (aynı hücre) | 1.5190310153796762 / 1.3142868954143399 = 1.16 | 0.75–1.33 | geçti |
| G4e | D/r0 8: AG artık RMS medyanı R5 / R0 (öngörü: büyür) | 1.3196715717860106 / 1.131922118774412 = 1.17 | > 1 (öngörü; kapı değil) | geçti |
| G4f | D/r0 8: ANEES_poz (pencere geometrik ort.) A1, A0C < A0 | A1/A0 0.086; A0C/A0 0.026 | ikisi de < 1 (etki uygulandı) | geçti |
| G4f | D/r0 8: ANEES_poz A3, AG < A0 | A3/A0 0.450; AG/A0 0.298 | ikisi de < 1 (beklenen; kapı değil) | geçti |
| G4g | D/r0 8: A0 ATE [m] (pencereler); temiz ana ızgara medyanı | [0.2626, 0.1853, 0.1635]; 0.0271 | sonlu, 0.005–3 m aralığında | geçti |
| G4h | D/r0 8: pencere tanımı ana ızgarayla aynı (karşılaştırılan poz sayısı) | [(400, 400), (400, 400), (400, 400)] | çiftler ±%5 içinde | geçti |
| G6 | D/r0 8: süre [s] ön-uç; arka-uç kol başına | 151; A0 328, A0C 308, A1 324, A3 1841, AG 380, A0C+A3 328 | bilgi (toplam süre kestirimi için) | — |
| G4a | D/r0 4: tau0/dt | {0.0314} | {0.0314} (v = 5 m/s) | geçti |
| G4b | D/r0 4: rejim uygulandı (ön-uç istatistikleri R0'dan farklı; zamansal eşleşme R5 / R0) | aynı: False; avg_temporal_matches 159.6 / 160.5; kp1 254 / 254 | farklı; anahtar nokta sayısı ±%10 içinde | geçti |
| G4c | D/r0 4: A1 oracle σ (kullanılan) | 2.476291114536357 | 2.476291114536357 (ana ızgarayla aynı) | geçti |
| G4c | D/r0 4: A0C cal_factor; A3 yama faktörü | 12.0; 7751 | 12.0; > 0 | geçti |
| G4d | D/r0 4: A3 σ̂ medyanı R5 / R0 (aynı hücre) | 1.075477652246771 / 1.0266457035657626 = 1.05 | 0.75–1.33 | geçti |
| G4e | D/r0 4: AG artık RMS medyanı R5 / R0 (öngörü: büyür) | 0.9423908254624146 / 0.8999371808661043 = 1.05 | > 1 (öngörü; kapı değil) | geçti |
| G4f | D/r0 4: ANEES_poz (pencere geometrik ort.) A1, A0C < A0 | A1/A0 0.166; A0C/A0 0.014 | ikisi de < 1 (etki uygulandı) | geçti |
| G4f | D/r0 4: ANEES_poz A3, AG < A0 | A3/A0 0.377; AG/A0 0.409 | ikisi de < 1 (beklenen; kapı değil) | geçti |
| G4g | D/r0 4: A0 ATE [m] (pencereler); temiz ana ızgara medyanı | [0.0449, 0.0514, 0.0507]; 0.0271 | sonlu, 0.005–3 m aralığında | geçti |
| G4h | D/r0 4: pencere tanımı ana ızgarayla aynı (karşılaştırılan poz sayısı) | [(400, 400), (400, 400), (400, 400)] | çiftler ±%5 içinde | geçti |
| G6 | D/r0 4: süre [s] ön-uç; arka-uç kol başına | 115; A0 206, A0C 184, A1 220, A3 1137, AG 192, A0C+A3 177 | bilgi (toplam süre kestirimi için) | — |

ÖN TEST: GEÇTİ

G5 (çıktı ve toplama): `hizli_analiz.py … --butunluk` → "bütünlük: 36/36 kayıt; eksik 0; durum {'ok': 36}; boş alanlı 'ok' 0;
code_fp {'08e28b83407c': 36}; tam = True"; `--butunluk` olmadan → "DUR: 1 birim; test tablosu en az 6 birimle basılır" (test
tablosu basılmadı). **Geçti.**

G6 (süre ve disk): ön-uç 151 s (D/r₀ 8) / 115 s (D/r₀ 4); hücre başına altı kolun toplamı 3509 s (D/r₀ 8) / 2116 s (D/r₀ 4)
(makine yükü ≈ 60 iken ölçüldü; A3 tek başına 1841 / 1137 s). Toplam ≈ 36 × (151 + 3509) + 36 × (115 + 2116) ≈ 212 000 işçi-saniyesi
≈ 59 işçi-saati; 6 arka-uç işçisiyle ≈ 9.5–10 saat (yük sürerse daha uzun). /tmp boş alanı 80 GB (≥ 20 GB). **Geçti.**

Yön öngörüleri (kapı değil; §3'ün son paragrafı):

- Zamansal eşleşme sayısı R5'te R0'dakinden büyük olmaz: 122.0 ≤ 134.6 (D/r₀ 8), 159.6 ≤ 160.5 (D/r₀ 4). **Tuttu.**
- A1'in konum ANEES'i R5'te R0 değerinden küçüktür: D/r₀ 4'te 153.4 < 220.3 (**tuttu**); D/r₀ 8'de 90.5 > 48.4 (**tutmadı**).
  Bu not sondadan önce yazılmıştı; keşif sondası (tohum 0, dokuz sekans) da D/r₀ 8'de aynı yönde çelişmişti. Tek birimlik
  gözlemdir, test değildir; Sapma 5'in HZ hipotezleri ön kayıttaki hâliyle kalır ve yalnız tam ızgarayla değerlendirilir.

Karar: G1–G6 kapılarının hepsi geçti, hiçbir DUR ölçütü tetiklenmedi. 72 hücrelik ızgara `scratchpad/hizli_baslat.sh` ile
başlatılır (1 ön-uç + 6 arka-uç işçisi; bellek payı öbür koşumlarla paylaşıldığı için ön-uç işçisi 1). Koşum içi erken denetim:
V1_01_easy tohum 1 hücreleri ön test kayıtlarıyla bit-aynı olmalı (tutmazsa koşucu DUR ile çıkar).
