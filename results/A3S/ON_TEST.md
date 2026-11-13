# A3S ayırma kolu — ön test (ON_KAYIT_FAZ2.md §5 Sapma 4)

Yazım: 2026-11-11 13:32 (dosya 13:32:42, koşum günlüğü 13:32:45), **ön test koşulmadan ve hiçbir A3S (k = 1) sonucu görülmeden önce.**
Kural: uzun koşum ön-test kuralı (yedi madde). Bu belge geçmeden 225 hücrelik ızgara başlatılmaz.

## 1. Kimlik (madde 1)

| Öğe | Değer |
|---|---|
| Kod parmak izi | `08e28b83407c` (ana ızgarayla aynı; parmak izli hiçbir dosya değişmedi) |
| Ön kayıt | `ON_KAYIT_FAZ2.md` sha256 `140472d6…bf83` (Sapma 4 dahil) |
| Koşucu / analiz | `a3s_kos.py` `dd9af6e1…`, `a3s_analiz.py` `91b7f120…` (tam özetler `kimlik_sha256.txt`) |
| Yapılandırma | `config.yaml` `92e1dc7f…` (ana ızgaradaki dosya) |
| Pencereler | `evidence/pencereler_faz2_kapsam1m.json` `cd788381…` |
| Eşleme kaynağı | `evidence/FAZ2/records.jsonl` `511b9a6f…` (5880 kayıt) |
| Veri | `SLAM_DATA_ROOT=<data-disk>/sdata` (EuRoC) |
| Ortam | `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=4` (RTX 3090; koşucu başka kartta durur) |
| Hücre tanımı | `turb.experiments.cells_FAZ2` ile aynı alanlar; fark: `arms`, `backend_argv=("--turb_cal_factor","1.0")`, `tag="a3s"` |

## 2. Ön test hücreleri (madde 3)

V1_01_easy, tohum 0, üç pencere (100, 500, 900), D/r₀ ∈ {8, 0}; dört hücre, 21 kayıt:

| Hücre | Etiket | k | Kollar | Kayıt |
|---|---|---|---|---|
| denetim, D/r₀ 8 | `a3s_ontest_k12` | 12 | A0, A0C+A3 | 6 |
| A3S, D/r₀ 8 | `a3s_ontest` | 1 | A0, A0C+A3 (= A3S) | 6 |
| denetim, D/r₀ 0 | `a3s_ontest_k12` | 12 | A0, A0C+A3 | 6 |
| A3S, D/r₀ 0 | `a3s_ontest` | 1 | A0C+A3 (= A3S) | 3 |

Aynı (sekans, tohum, D/r₀) için ön-uç iki kez koşulur; bu, ön-ucun tekrar üretilebilirliğini de sınar.

## 3. Beklenen değerler (madde 2, 4, 5) — sonuçtan önce yazıldı

Ana ızgara referansları (V1_01_easy, tohum 0; `evidence/FAZ2/records.jsonl`):

- D/r₀ 8: sekans düzeyi konum ANEES'i (üç pencerenin geometrik ortalaması) A0 413.7, A1 73.8, A3 386.0, AG 493.9,
  A0C+A3 (k = 12) 16.8; A3'ün ortalama ATE'si 0.12725 m; A3 `sigma_turb_used_mean` 1.2889 px.
- D/r₀ 0: ANEES A0 2024.0, A3 1421.1, AG 1587.2, A0C+A3 (k = 12) 21.2; A3 ortalama ATE 0.02353 m; A3 σ̂ 0.5785 px.
- Ön-uç (D/r₀ 8 / 0): stereo eşleşme 11.339 / 4.245, üçgenlenen 8.702 / 4.217, zamansal 136.755 / 162.445,
  etkin derinlik izi 75.699 / 95.518.

| Kapı | Ölçü | Beklenen | DUR ölçütü |
|---|---|---|---|
| G1 bit-aynılık (özellik kapalı = k 12) | denetim hücrelerinde A0 ve A0C+A3, ayrıca A3S hücresindeki A0 (k = 1 argümanıyla): `m_ate_se3`, `n_nees_pos_mean`, `n_nees_rot_mean`, `m_rpe_trans` ana ızgara kaydıyla `==` | 15 / 15 bit-aynı (ör. A0 D/r₀ 8 pencere 100 ATE 0.0695479719992662) | bir kayıt bile farklıysa |
| G2 ön-uç | `frontend_stats` (süre alanı hariç) ana ızgarayla aynı; `fe_gpu` RTX 3090; `code_fp` 08e28b83407c | 4 / 4 hücrede aynı | biri farklıysa |
| G3 kol kimliği | A3S kayıtlarında `cal_factor` = 1.0, `n_patch_factors` = 0, `dimm_sigma_scale` = 1.0; denetim kayıtlarında `cal_factor` = 12.0 | hepsi | biri tutmazsa |
| G4 σ̂ | A3S `sigma_turb_used_mean` / ana ızgara A3 | D/r₀ 8: 1.29–1.31 px (oran 1.00–1.02); D/r₀ 0: 0.578–0.583 px | oran [0.95, 1.05] dışında (1.2245–1.3533; 0.5496–0.6074) |
| G5 etki uygulandı, makul | A3S sekans düzeyi konum ANEES'i | D/r₀ 8: A3'ün [0.6, 1.6] katı = [231.6, 617.6], A1'den (73.8) büyük, k = 12'ninkinin (16.8) en az 5 katı; D/r₀ 0: [852.6, 2273.7], A0'dan (2024.0) küçük | A3'ün [0.33, 3] katı dışında; ya da A3S kaydı A0 ya da k = 12 kaydıyla bit-aynı (etki uygulanmamış) |
| G6 doğruluk | A3S ortalama ATE / A3 ortalama ATE | [0.6, 1.6]: D/r₀ 8 [0.0763, 0.2036] m; D/r₀ 0 [0.0141, 0.0376] m | bir pencerede ATE, ana ızgara kollarının (A0, A1, A3, AG, A0C+A3) o penceredeki [0.5 × en küçük, 2 × en büyük] aralığının dışında; ıraksama |
| G7 çıktı | `on_test.jsonl` 21 kayıt; `status` ok, `m_valid`, `n_valid`; ATE, RPE, ANEES, `n_n`, süre, ön-uç, kart, parmak izi alanları dolu; hücre hatası 0 | hepsi | biri eksikse |
| G8 analiz | `a3s_analiz.py oz_test`: ana ızgaranın A3 kayıtları "A3S" diye verilince A3S / A3 tam 1 ve eşdeğer, A3S / AG tablosu `analiz.json` H4 ii ile sayı sayı aynı; sentetik TOST doğru; `a3s_analiz.py ontest` dosyayı okur | geçer | geçmezse |
| G9 süre | A3S kol koşumu `t_backend` (yük altında) ve ön-uç süresi | arka-uç 150–800 s (yüksüz ana ızgara 170–312 s; şimdiki yük ≈ 1.85×); ön-uç 50–450 s | arka-uç > 1600 s (bütçe iki katını aşar: plan yeniden yapılır) |

## 4. İşlemsel tanımlar (sonuçtan önce; Sapma 4'ün uygulanışı)

- **"`frontend_stats` aynı"** koşulu `total_execution_time_seconds` dışındaki bütün alanlar için `==` ile uygulanır
  (o alan duvar saati ölçümüdür; keşif sondasının ana ızgarayla bit-aynı çıkan hücrelerinde de yalnız o alan farklıydı).
- **Bit-aynılık alanları:** `m_ate_se3`, `n_nees_pos_mean`, `n_nees_rot_mean`, `m_rpe_trans`.
- **Hücre geçerliliği (eşlemeye girme):** ön-uç aynı, kart RTX 3090, parmak izi aynı, üç pencerenin A3S kaydı var,
  `cal_factor` = 1, yama faktörü 0, `dimm_sigma_scale` = 1; A0 koşulan dokuz hücrede A0 bit-aynı. Koşucu her biten
  hücrede ön-uç ve A0 denetimini yapar; tutmazsa yeni hücre göndermez (çıkış kodu 2).
- **Yorum kuralının üçüncü maddesi** ("A3S / A3 medyanı < 0.87 ve anlamlı") ikinci maddeyle aynı ölçütle uygulanır:
  D/r₀ 8, 12, 16'nın en az ikisinde, Holm sonrası.
- **HS3 kararı:** medyan < 1 ve Holm p < 0.05 (D/r₀ 8, 12, 16'nın üçünde de). "Beklenen medyan 0.78–0.92" ayrıca
  içinde / dışında diye raporlanır; karara girmez.
- **σ̂ oranı:** hücre başına `sigma_turb_used_mean` (A3S) / (ana ızgara A3).
- **Karar yalnız tam veriyle verilir** (beş seviye, 225 hücre). Ara bakışlar yalnız denetim içindir (kayıt sayısı, hata,
  bit-aynılık, süre); ara bakışta test tablosu üretilmez, örneklem değiştirilmez.
- **Holm ailesi:** karşılaştırma ve metrik başına beş seviye (HS1: A3S/A3 ANEES; HS2: A3S/A3 ATE; HS3: A3S/AG ANEES).

## 5. Sonuçlar

Ön test: 2026-11-11 13:32:45 → 14:02:11 (29.4 dk), çıkış kodu 0, "4 hücre, sorunlu 0"; `on_test.jsonl` 21 kayıt,
günlük `on_test_1332.log`. Bu bölüm 14:14'te, `a3s_analiz.py ontest` çıktısı ve kayıtlar üzerinden yazıldı.

| Kapı | Ölçülen | Beklenen | Sonuç |
|---|---|---|---|
| G1 bit-aynılık | 15 / 15 kayıt ana ızgarayla bit-aynı (A0 D/r₀ 8 pencere 100 ATE 0.0695479719992662; k = 12 A0C+A3 kayıtları da aynı) | 15 / 15 | geçti |
| G2 ön-uç | 4 / 4 hücrede `frontend_stats` aynı (stereo 11.339 / 4.245, üçgenlenen 8.702 / 4.217, zamansal 136.755 / 162.445, etkin iz 75.699 / 95.518); RTX 3090; `08e28b83407c` | 4 / 4 | geçti |
| G3 kol kimliği | A3S: `cal_factor` 1.0, `n_patch_factors` 0, `dimm_sigma_scale` 1.0; denetim: `cal_factor` 12.0 | hepsi | geçti |
| G4 σ̂ | D/r₀ 8: 1.2983 px (σ̂ 1.3082), oran 1.007; D/r₀ 0: 0.5763 px (σ̂ 0.5783), oran 0.996 | 1.29–1.31 px (1.00–1.02); 0.578–0.583 px | D/r₀ 8 beklenen aralıkta. **D/r₀ 0 beklenen aralığın 0.0017 px altında** (beklenti ıskalandı, aralık değiştirilmedi); DUR sınırının [0.5496, 0.6074] içinde → geçti |
| G5 etki | sekans düzeyi ANEES D/r₀ 8: 472.2 (A3'ün 1.223 katı; A1 73.8'den büyük; 16.8'in 28 katı); D/r₀ 0: 1484.7 (A3'ün 1.045 katı; A0 2024.0'dan küçük); A0 ve k = 12 kayıtlarıyla bit-aynı değil | [231.6, 617.6]; [852.6, 2273.7] | geçti |
| G6 doğruluk | ortalama ATE A3S / A3: D/r₀ 8 0.1001 / 0.1272 m = 0.786; D/r₀ 0 0.0242 / 0.0235 m = 1.030; altı pencere ATE'si de ana ızgara kollarının [0.5 × en küçük, 2 × en büyük] aralığında; ıraksama yok | [0.6, 1.6] | geçti |
| G7 çıktı | 21 kayıt; hepsi `status` ok, `m_valid`, `n_valid`; boş alan 0; hücre hatası 0 | hepsi | geçti |
| G8 analiz | `oz_test` geçti (ön testten önce); `ontest` dosyayı okudu ve tabloyu üretti | geçer | geçti |
| G9 süre | A3S `t_backend` 493 s (D/r₀ 8), 546 s (D/r₀ 0); ön-uç 104–321 s (yük ortalaması 60–90) | 150–800 s; 50–450 s | geçti |

A3S pencere kayıtları (ATE m / konum ANEES / dönme ANEES / n):

| D/r₀ | pencere 100 | pencere 500 | pencere 900 |
|---|---|---|---|
| 8 | 0.06244 / 607.87 / 245.74 / 29 | 0.10835 / 366.19 / 597.54 / 58 | 0.12943 / 473.13 / 283.89 / 57 |
| 0 | 0.02891 / 3184.25 / 648.01 / 29 | 0.02661 / 1190.17 / 915.44 / 52 | 0.01722 / 863.55 / 444.67 / 46 |

Notlar:

- Ön test hücresi (V1_01_easy, tohum 0, D/r₀ 8) ızgaranın 45 biriminden biridir; buradaki sayılar karar üretmez, karar
  kuralları (Sapma 4, §4) bu sayılardan önce yazılmıştır. Izgaranın `a3s` etiketli aynı hücresi bu değerleri bit-aynı
  üretmelidir (ek tutarlılık denetimi; ilk biten hücrelerde bakılır).
- D/r₀ 0 hücresi ızgarada yoktur (yalnız süre kalırsa, betimsel).
- **Kaynak ve süre (madde 6):** 225 ön-uç koşumu × ≈ 250 s (yük altında) ve 234 kol koşumu × ≈ 500 s ≈ 32 işlemci-saat.
  Ekran kartı belleği kısıtlıdır (14:14'te 17.1 / 24.6 GB dolu: rakip DROID ≈ 8.7 GB, onarım hazırlığı 2 × 5 GB,
  keşif sondası 3 × 1 GB; tepe toplamı ≈ 22.7 GB) → ızgara **tek ön-uç işçisiyle** (`--gpu_jobs 1 --jobs 4`) başlar;
  rakip kuyruğu (≈ 15:35), sonda ön-ucu ve onarım hazırlığı (≈ 17:00) bitince koşucu durdurulup üç işçiyle yeniden
  başlatılır (koşucu kaldığı yerden sürer; yalnız o an süren hücreler yeniden koşulur). Beklenen bitiş ≈ 22:00–23:00.
  Disk: `/` 83 GB boş; spool geri basınçlı (`/tmp/turb_spool_a3s`).
- **Karar: bütün kapılar geçti → 225 hücrelik ızgara ön kayıttaki sırayla (8, 16, 24, 12, 20) başlatılır.**
