# Derin onarma kolu — ön test (ON_KAYIT_FAZ2.md Sapma 3)

Tarih: 2026-11-11. Beklenen değerler sonuçlardan ÖNCE yazıldı (aşağıda "önceden görüldü" denenler dışında).

## 1. Beklenenler

| No | Denetim | Beklenen (önceden yazıldı) | Dur kuralı |
|---|---|---|---|
| K1 | Kimlik | code_fp 08e28b83407c; config.yaml sha 92e1dc7f…; MambaTM commit 79aa41f45d2f, ağırlık sha d207a6cd…; GPU RTX 3090; `onar_sekans.py`, `onarim_kos.py`, `uret_sekans.py` sha256 kaydedilir | fp farklı → dur |
| G1 | `--ad` ve `--tau0_over_dt 1.0` açık verilince üretim | otomatik adlı üretimle (aynı parametreler) bütün PNG'ler bayt-aynı | fark → dur |
| G2 | Onarım tekrarı (60 kare, O) | iki koşum bayt-aynı. **Önceden görüldü:** aynı (0f4d173d…) | — |
| G3 | Dosya yolu = kanca yolu (V1_01_easy, R0, D/r₀ 8, N; A0, A1, A3; 3 pencere) | ana ızgara tohum 0 kaydıyla pencere ATE farkı medyanı < %1 (kalibrasyon yaml'ı 12 haneye yuvarlı → bit-aynılık beklenmez); ANEES farkı < %10; ön-uç stereo eşleşme ortalaması farkı < %1 | ATE farkı medyanı > %10 → dur |
| G4 | Temiz dosya yolu (V1_01_easy, D/r₀ 0, σ = 1 DN, N; A0) | ana ızgara D/r₀ 0 tohum 0 A0 ile pencere ATE farkı medyanı < %1 | > %10 → dur |
| G5 | Temiz görüntünün onarımı VO'yu bozmaz (V1_01_easy, D/r₀ 0, O; A0) | sekans ATE oranı (O / N) [0.67, 1.5]; stereo eşleşme ortalaması ≥ 0.8 × N; ıraksama yok | oran > 3 ya da ıraksama → dur, incele |
| G6 | Duman: her kol × her yol | G3: A0, A1, A3 (N yolu); G5: A0 (O yolu); C kipi 60 karede çalışır, çıktısı O'dan farklı; bütün kayıtlarda status ok, m_/ms_/n_ alanları dolu | hata → dur |
| G7 | Etki gerçekten uygulandı | O çıktısı girdiden farklı (ortalama mutlak fark > 1 DN); A1'in `sigma_turb_used` değeri N ve O hücrelerinde ana ızgaradakiyle aynı (4.412 px) | — |
| G8 | Süre | onarım O ≈ 0.06–0.13 s/kare/kamera (GPU paylaşımlı), C ≈ 5 × O; ön-uç 180–400 s; kol başına arka-uç 6–25 dk | toplam > 10 saat → ızgarayı küçült |
| G9 | Çıktı | `on_test.jsonl` satır sayısı = hücre × pencere × kol (3×3 + 3×1 + 3×1 = 15); `onarim` alanı dolu | eksik → dur |

İlk ızgara hücrelerinde erken denetim (V1_01_easy R5 N/O/C, R0 O): görüntü düzeyi ölçü (`onarim_kos.py olc`), beklenen
temize göre Farneback RMS oranı R0'da O/N [0.9, 1.1], R5'te O/N [0.6, 0.95], C/N [0.75, 1.0]; PSNR değişimi [−0.5, +1.5] dB.
R5 N hücresi sondanın R5/fe0 kaydıyla (aynı sekans, tohum 0) G3'teki ölçüyle uyumlu olmalı.

### 1b. Veri aşaması (üretim + onarım) için görüntü düzeyi beklentiler — 13:10'da, `olc --ontest` koşulmadan ÖNCE yazıldı

Ölçü: cam0, her 20 karede bir, başvuru = aynı sekansın gürültüsüz temiz rektifiye karesi (D/r₀ 0, σ = 0).
Hazırlık (üretim + onarım) ızgarası bu aşama geçince, PyEco koşumları G3–G9 geçince başlatılır.

| No | Denetim | Beklenen | Dur kuralı |
|---|---|---|---|
| G10a | Temiz + gürültü (T, N) | PSNR 47.8 ± 0.5 dB (σ = 1 DN + nicemleme); RMS akış < 0.1 px; \|dOrt\| < 0.1 DN | PSNR < 45 dB → dur |
| G10b | Temiz + gürültü, onarılmış (T, O) | dOrt +1…+3 DN (ajanın ölçtüğü +2 DN kayma); PSNR 30–42 dB; SSIM ≥ 0.90; RMS akış < 0.3 px (ağ temiz görüntüye geometrik bozulma eklememeli) | RMS akış > 1 px → dur |
| G10c | R0, D/r₀ 8, N | RMS akış 1.5–4.5 px (oracle σ_tilt 2.19 px/eksen → 2B RMS ≈ 3.1 px; Farneback düzgünleştirir); PSNR 18–26 dB | RMS akış < 0.5 px → bozulma uygulanmamış, dur |
| G11 | Tam sekans onarımı (1300 kare × 2 kamera, orkestra yoluyla) | hatasız biter; `onarim_meta.json` iki kamerayı da içerir; kare sayısı girdiyle aynı; s/kare 0.06–0.20 (yük altında) | hata → dur |

## 2. Sonuçlar

### 2a. Veri aşaması (13:12'de tamamlandı; hepsi geçti)

| No | Beklenen | Ölçülen | Sonuç |
|---|---|---|---|
| K1 (kısmi) | code_fp 08e28b83407c; ağırlık d207a6cd…; commit 79aa41f45d2f; RTX 3090 | `kos` çıktısı: code_fp 08e28b83407c; `onarim_meta.json`: ağırlık d207a6cdea95…, commit 79aa41f45d2f…, GPU RTX 3090 | GEÇTİ |
| G1 | 2600 PNG bayt-aynı | `V1_01_easy__R0_dr8_s0` ile otomatik adlı üretim: 2600/2600 PNG bayt-aynı; data.csv, GT, turb_meta.json aynı; sensor.yaml yalnız adı taşıyan yorum satırında farklı | GEÇTİ |
| G2 | bayt-aynı | bayt-aynı (0f4d173d…; önceden görüldü) | GEÇTİ |
| G10a | PSNR 47.8 ± 0.5 dB; RMS akış < 0.1 px; \|dOrt\| < 0.1 DN | PSNR 47.99 dB; SSIM 0.990; RMS akış 0.031 px; dOrt −0.03 DN | GEÇTİ |
| G10b | dOrt +1…+3 DN; PSNR 30–42 dB; SSIM ≥ 0.90; RMS akış < 0.3 px | dOrt +2.14 DN; PSNR 36.04 dB; SSIM 0.978; RMS akış 0.294 px (sınıra yakın) | GEÇTİ |
| G10c | RMS akış 1.5–4.5 px; PSNR 18–26 dB | RMS akış 3.51 px; PSNR 22.92 dB; SSIM 0.737 | GEÇTİ |
| G11 | hatasız; iki kamera; kare sayısı aynı; 0.06–0.20 s/kare | 1300 × 2 kare; 0.078 / 0.100 s/kare; toplam 455.8 s (ileri 231 s, oku 71 s, yaz 64 s); tepe 3652 MiB; ort. \|fark\| 2.88 DN | GEÇTİ |

Kaynak: `goruntu_olcu_on_test.jsonl` (n_kare 65), `hazirlik.jsonl`, `on_test_1257.log`.

Not (13:14): `onarim_kos.py` dosyasına yalnız `hazirla` için `--kipler` süzgeci eklendi (önce N ve O, sonra C hazırlansın
diye). Üç yerde değişiklik: argparse seçeneği, `hazirla` çağrısındaki süzgeç, kullanım satırı. `plan` çıktısı değişiklikten
önce ve sonra bayt-aynı (45 hücre); `hucre()`, `kos()`, `uret()`, `onar()` değişmedi. Eski sha fbaa431f… (ön-test bu
sürümle başladı), yeni sha 99270271….

Karar (13:16): veri aşaması geçti → ızgara verisinin hazırlığı (üretim + onarım) başlatıldı. PyEco koşumları G3–G9
sonuçlanana kadar başlatılmaz.

### 2c. Erken denetim, ızgara verisi (V1_01_easy; 13:43'te ölçüldü, beklenti §1'de önceden yazılıydı)

| Ölçü | Beklenen | Ölçülen | Sonuç |
|---|---|---|---|
| R0, RMS akış O / N | [0.9, 1.1] | 3.148 / 3.507 = 0.898 | **sınırın hemen altında** (0.002 dışarıda) |
| R5, RMS akış O / N | [0.6, 0.95] | 2.465 / 3.283 = 0.751 | beklenen içinde |
| R0, PSNR değişimi O − N | [−0.5, +1.5] dB | 22.69 − 22.92 = −0.22 dB | beklenen içinde |
| R5, PSNR değişimi O − N | [−0.5, +1.5] dB | 24.78 − 23.44 = +1.35 dB | beklenen içinde |
| R5, C / N | [0.75, 1.0] | C henüz üretilmedi (ikinci geçiş) | bekliyor |
| R5 N ile sondanın R5/fe0 kaydı | uyumlu | sondada V1_01_easy R5/fe0 kaydı henüz yok | bekliyor |

Not: R0'da onarma, yarı-statik tilti beklediğimden biraz fazla (%10) azalttı; aralık bir beklentiydi, DUR ölçütü
değildi; sonuç olduğu gibi kaydedildi, aralık değiştirilmedi. Ortalama parlaklık kayması O'da +2.0…+2.3 DN (bilinen).
Kaynak: `goruntu_olcu_erken.jsonl` (n_kare 65; tam ızgara ölçüsü ayrıca `goruntu_olcu.jsonl`'e yazılacak).

### 2b. PyEco ön testi (G3–G9) — sonuçlar (14:20; `onarim_analiz.py ontest`, `on_test.jsonl` 15 kayıt, 56.2 dk)

| No | Ölçülen | Beklenen | Sonuç |
|---|---|---|---|
| G3 | dosya yolu (R0, D/r₀ 8, N) / ana ızgara tohum 0: pencere ATE farkı medyanı **%42.2** (en çok %261), ANEES farkı medyanı **%42.6** (en çok %474; A0'da üç pencerede de dosya yolu yüksek: 1628 / 1370 / 1560, ana ızgara 284 / 267 / 933); ön-uç istatistikleri (süre dışında bütün alanlar) **bit-aynı** (stereo eşleşme farkı %0) | ATE farkı medyanı < %1; ANEES < %10; stereo eşleşme < %1 | **GEÇMEDİ — dur kuralı (> %10) tetiklendi.** PyEco ızgarası başlatılmadı; inceleme aşağıda |
| G4 | temiz dosya yolu / ana ızgara D/r₀ 0 A0: ATE farkı medyanı %0.04 (en çok %0.68), ANEES farkı en çok %0.05; ön-uç istatistikleri bit-aynı | < %1 | geçti |
| G5 | temiz + O / temiz N: sekans ATE oranı 1.069; stereo eşleşme 4.35 / 4.24 (1.02 ×); ıraksama yok | [0.67, 1.5]; ≥ 0.8 × | geçti |
| G6 | A0, A1, A3 (N yolu) ve A0 (O yolu) koştu; 15 / 15 kayıt `status` ok, alanlar dolu; C kipi 60 karede veri aşamasında çalıştı (`_test_c60`) | hata yok | geçti |
| G7 | O çıktısı girdiden farklı: ortalama mutlak fark 3.2 / 6.3 / 3.2 DN (R0, kare 100 / 600 / 1100); A1 `sigma_turb_used` 4.412249 px (N hücresi; ön testte türbülanslı O hücresi yok) | > 1 DN; 4.412 px | geçti |
| G8 | onarım O 0.078–0.130 s/kare; ön-uç 153–204 s; arka-uç A0 527 s, A1 445 s, **A3 2167 s (36 dk)** (yük ortalaması 60–90) | 0.06–0.13 s/kare; 180–400 s; 6–25 dk | A3 süresi beklenenin üstünde (dur kuralı: toplam > 10 saat; 45 hücre × ≈ 52 dk ≈ 39 işlemci-saat, 6–8 işle ≈ 5–6.5 saat → sınırın altında) |
| G9 | 15 satır, `onarim` alanı dolu, eksik alan 0 | 15 | geçti |

#### G3 incelemesi (dur kuralı tetiklendikten sonra; PyEco ızgarası beklemede)

Eldeki olgular (yeni koşum olmadan, 14:20):

1. Ön-uç istatistikleri (anahtar nokta, stereo / zamansal eşleşme, üçgenlenen, etkin iz, `avg_s_d1`; ortalama ve standart
   sapmalar) iki yolda son haneye kadar aynı → ön-uç iki yolda aynı görüntüleri görüyor.
2. Arka-uca geçen kalibrasyon farklı (işlemcide, `StereoMatcher.load_calibration` iki dizinde çağrılarak ölçüldü):
   `rect_K` kanca yolu (436.244295647102, 364.44123458862305, 256.95167541503906), dosya yolu (436.244308600007,
   364.44122314453125, 256.9516830444336) → göreli fark 3.0e-8 (fx'te 1.3e-5 px, cx'te 1.1e-5 px, cy'de 7.6e-6 px);
   taban 2.2e-12; T_B_C 5e-13. Neden: ön-uç zaten rektifiye çifte yeniden `stereoRectify` (alpha = 0) uyguluyor, odak
   1e-8 düzeyinde yeniden ölçekleniyor. Dosya yolundaki rektifikasyon haritası özdeşlikten en çok 2.3e-5 px sapıyor
   (`cv2.remap` sabit noktalı 1/32 px → görüntü değişmez; olgu 1 ile uyumlu).
3. Dosya yolu sonuçları aynı sekansın ana ızgaradaki tohumlar arası aralığının içinde ya da hemen üstünde (A0 ANEES
   pencere 100: tohum 0–4 284 / 1262 / 461 / 800 / 1814, dosya 1628; pencere 500: 267–1250, dosya 1370; pencere 900:
   711–1463, dosya 1560; A1 pencere 900 ATE: 0.057–0.146 m, dosya 0.418 m).

**Hipotez:** arka-uç, 1e-8 düzeyindeki kalibrasyon farkını ayrık kararlar (kapılama, landmark ayıklama, RANSAC)
yoluyla kaotik olarak büyütüyor; temiz görüntüde zayıf (G4: %0.04), D/r₀ 8'de güçlü. G3'ün "fark < %1" beklentisi
arka-ucun küçük pertürbasyona düzgün yanıt verdiği varsayımına dayanıyordu.

**Sınama (yalnız arka-uç, ekran kartı yok; `g3_duyarlilik.py` → `evidence/ONARIM/g3_duyarlilik.jsonl`). Beklenenler
sonuçtan ÖNCE (14:20) yazıldı:**

| Sınama | Girdi | Beklenen | Yorum kuralı |
|---|---|---|---|
| T-K (denetim) | MH01, D/r₀ 8, tohum 0 kanca dökümü (A3S ızgarasından kopya), değiştirilmeden, A0 | ana ızgara A0 kaydıyla bit-aynı (üç pencere; ATE, RPE, iki ANEES) | değilse sınama düzeneği yanlış → dur |
| T-D1 (yeniden üretim) | V1_01_easy, D/r₀ 8, tohum 0 kanca dökümü + dosya yolu kalibrasyonu (`rect_K`, taban, `T_BS`, `R1`), hücre = ön testin R0 N hücresi, A0 ve A1 | ön testin R0 N kayıtlarıyla **bit-aynı** | bit-aynıysa G3 farkının tamamı kalibrasyon pertürbasyonundan gelir. Değilse gözlemler de farklıdır → dosya yolu ön-ucu yeniden koşulup dökümler karşılaştırılır; ızgara başlamaz |
| T-D2 | MH01 kanca dökümü + MH01 dosya yolu kalibrasyonu, A0 | ana ızgarayla bit-aynı değil; en az bir pencerede ATE farkı ≥ %5 | üç pencerede de < %1 ise duyarlılık sekansa özgüdür (hipotez zayıflar) |
| T-E | MH01 kanca dökümü, yalnız fx ve fy × (1 + 1e-12), A0 | bit-aynı değil; fark T-D2 ile aynı mertebede (en az bir pencerede ≥ %5) | bit-aynı ya da < %0.1 ise duyarlılık eşiği 1e-12 ile 3e-8 arasındadır |

**Karar kuralı (önceden):** T-K geçer ve T-D1 bit-aynı çıkarsa G3 "veri / boru hattı hatası" değil, "kapının varsayımı
yanlış" diye sınıflanır; dosya yolu N hücresi kanca yolunun sayısal olarak pertürbe edilmiş ikizi sayılır; ön kayda
"Sapma 3 eki" yazılır (onarım ızgarası sonucu görülmeden) ve ızgara başlatılır. Onarım karşılaştırmaları zaten yalnız
dosya yolu içindedir (N / O / C); ana ızgara hücreleriyle sayı sayı karşılaştırma yapılmaz. Duyarlılık makalede sınır
olarak sayılarıyla yazılır. T-D1 bit-aynı çıkmazsa ızgara başlamaz, kullanıcıya sorulur.

**Düzeltme (14:26; hiçbir sınama koşulmadan önce).** `stereo_front_end.landmarks_to_obs` okununca görüldü: dökümdeki
gözlem satırı [id, u, v, has_depth, fresh_depth, X, Y, Z] biçimindedir ve X, Y, Z ön-uçta ön-ucun kendi
kalibrasyonuyla üçgenlenir. Dosya yolunda bu üç sütun da 1e-8 düzeyinde farklı olmalıdır; yalnız kalibrasyon alanlarını
değiştirmek ön test kayıtlarını bit-aynı üretemez. Yukarıdaki T-D1 beklentisi ("bit-aynı") bu yüzden **geri alındı**;
yerine aşağıdaki üç sınama geçer (T-K, T-D2, T-E değişmedi):

| Sınama | Girdi | Beklenen | Yorum kuralı |
|---|---|---|---|
| T-F (döküm karşılaştırması; belirleyici) | V1_01_easy, D/r₀ 8, tohum 0: kanca dökümü (A3S ızgarasından kopya) ve dosya yolu dökümü (`V1_01_easy__R0_dr8_s0` üzerinde ön-uç yeniden koşulur, döküm saklanır) | kare sayısı, zaman damgaları, kare başına gözlem sayısı aynı; id, u, v, has_depth, fresh_depth sütunları bütün karelerde bit-aynı; NaN deseni aynı; X, Y, Z farkı (‖Δp‖ / ‖p‖) en çok < 1e-6 (beklenen ≈ 1e-8–1e-7); kalibrasyon farkı olgu 2'deki gibi; yeniden koşulan dosya yolu ön-ucunun istatistikleri ön testle bit-aynı | ayrık fark (kare, gözlem sayısı, id, u, v, bayrak) varsa iki yol aynı gözlemi üretmiyor → boru hattı farkı; ızgara başlamaz |
| T-D1 (yeniden tanımlandı: tekrarlanabilirlik) | yeniden koşulan dosya yolu dökümü, hücre = ön testin R0 N hücresi, A0 | ön testin R0 N A0 kayıtlarıyla bit-aynı (üç pencere) | değilse dosya yolu belirlenimci değil ya da düzenek yanlış → dur |
| T-D1k (ek duyarlılık örneği) | V1_01_easy kanca dökümü, yalnız `rect_K`, taban, `T_BS`, `R1` dosya yolundakiyle değiştirilir (gözlemler aynı), ana ızgara hücresi, A0 | ana ızgarayla bit-aynı değil; en az bir pencerede ATE farkı ≥ %5 | — |

**Karar kuralı (güncel, sonuçtan önce):** (i) T-K bit-aynı, (ii) T-F'de yalnız sürekli büyüklüklerde ≤ 1e-6 göreli fark,
(iii) T-D1 bit-aynı ve (iv) T-D2, T-E, T-D1k'den en az biri en az bir pencerede ≥ %5 ATE farkı gösterirse G3 "kapının
varsayımı yanlış" diye sınıflanır (arka-uç belirlenimci + girdiler yalnız 1e-8 düzeyinde farklı ⇒ çıktı farkı o farkın
büyümesidir) ve yukarıdaki yol izlenir (ön kayda Sapma 3 eki, sonra ızgara). (i)–(iii)'ten biri tutmazsa ızgara başlamaz.
(iv) tutmaz ama (i)–(iii) tutarsa sonuç kullanıcıya sorulur.

#### G3 incelemesi — sonuçlar (15:04; `evidence/ONARIM/g3_duyarlilik.jsonl`)

| Sınama | Ölçülen | Beklenen | Sonuç |
|---|---|---|---|
| T-K | MH01, D/r₀ 8, tohum 0 kanca dökümü değişmeden, A0: üç pencerede ATE, RPE ve iki ANEES ana ızgara kaydıyla bit-aynı (1184 s) | bit-aynı | geçti |
| T-F | dosya yolu ön-ucu yeniden koşuldu (190 s; ön-uç istatistikleri ön testle aynı). İki döküm: 1300 kare, zaman damgaları, kare başına gözlem sayısı, id / u / v / has_depth / fresh_depth sütunları, NaN deseni, anahtar noktalar ve tanımlayıcılar bit-aynı (193 466 gözlem). Fark yalnız X, Y, Z'de: 98 409 satırın 4 493'ü bit-aynı, göreli fark medyan 4.0e-8, %99 1.6e-7, en büyük 2.4e-7. Kalibrasyon: `rect_K` göreli 3.0e-8, taban 2.2e-12 | ayrık içerik bit-aynı; sürekli büyüklüklerde ≤ 1e-6 (beklenen düzey ≈ 1e-8–1e-7) | geçti (en büyük fark beklenen düzeyin biraz üstünde, 1e-6 eşiğinin altında) |
| T-D1 | yeniden koşulan dosya yolu dökümü, hücre = ön testin R0 N hücresi, A0 (435 s): üç pencerede ATE, RPE ve iki ANEES ön test kaydıyla bit-aynı (ATE 0.058958 / 0.176120 / 0.214404 m) | ön testin R0 N A0 kayıtlarıyla bit-aynı | geçti |
| T-D2 | MH01 kanca dökümü + dosya yolu kalibrasyonu, A0: bit-aynı değil; pencere ATE farkı +%79.5, +%10.6, −%18.4; ANEES −%42.9, +%45.2, +%20.1 | bit-aynı değil; en az bir pencerede ≥ %5 | beklendiği gibi |
| T-E | MH01 kanca dökümü, yalnız fx ve fy × (1 + 1e-12) (436.244295647102 → 436.24429564753825), A0: bit-aynı değil; ATE +%5.1, −%17.3, +%22.5; ANEES −%40.4, −%20.3, +%4.1 | bit-aynı değil; T-D2 ile aynı mertebe | beklendiği gibi |
| T-D1k | V1_01_easy kanca dökümü + dosya yolu kalibrasyonu (gözlemler aynı), A0: bit-aynı değil; ATE −%26.3, +%471.3, +%128.9; ANEES +%322, +%1820, +%144 | bit-aynı değil; en az bir pencerede ≥ %5 | beklendiği gibi |

Notlar:

- 1e-12 göreli pertürbasyonun etkisi (T-E: en çok %22) 3e-8'lik pertürbasyonunkiyle (T-D2: en çok %79) aynı
  mertebededir; yanıt pertürbasyonun büyüklüğüyle ölçeklenmiyor. Bu, ayrık kararlar üzerinden kaotik büyüme hipoteziyle
  uyumlu, düzgün (türevlenebilir) bir duyarlılıkla uyumsuzdur.
- V1_01_easy'de iki pertürbe koşum da (ön testin dosya yolu N'si ve T-D1k) pencere 500 ve 900'de ana ızgara tohum 0
  kaydından kötü çıktı. İki örnekle yön hakkında sonuç çıkarılmaz; sistematik yön olup olmadığını sayısal ikiz denetimi
  (aşağıda) dokuz sekansta ölçer.
- Arka-uç belirlenimcidir: aynı girdi aynı çıktıyı verir (T-K, T-D1). Tekrarlanabilirlik sorunu yoktur; sorun, girdideki
  son hane farkının pencere düzeyi sonucu değiştirmesidir.

**Karar (önceden yazılan kurala göre):** (i) T-K bit-aynı: tuttu. (ii) T-F'de yalnız sürekli
büyüklüklerde ≤ 1e-6 göreli fark: tuttu (en büyük 2.4e-7). (iii) T-D1 bit-aynı: tuttu. (iv) T-D2, T-E ve T-D1k'nin üçü
de en az bir pencerede ≥ %5 ATE farkı gösterdi: tuttu. G3 "kapının varsayımı yanlış" diye sınıflanır; veri ya da boru
hattı hatası bulunmadı. Ön kayıttaki dur kuralı tetiklendiği için ızgaranın başlatılması bir sapmadır: ön kayda
"Sapma 3 eki" yazılır (onarım ızgarasından hiçbir sonuç görülmeden; sayısal ikiz denetimi N′ ile birlikte), sonra
ızgara başlar.

#### Koşum planı (sonuçtan önce)

- Izgara: `onarim_kos.py kos --bekle --gpu_jobs 1 --jobs 4`, nice 5 (bellek sınırlı: boş bellek ≈ 11 GB, ön-uç işçisi
  ≈ 2.1 GB, arka-uç işçisi ≈ 1.3 GB; ön-uç tek işçiyle hücre başına ≈ 5 dk, dört arka-uç işini besler). 45 hücre ×
  ≈ 52 dk ≈ 39 işçi-saat → ≈ 10 saat (dur sınırı: toplam > 10 saat; sınırda. Sonda bitince (≈ 20:30) işlemci boşalır).
- N dökümleri: `ikiz_kos.py sakla` (spool dizininden sabit bağlantı; ızgara koşucusu değişmez).
- Sayısal ikiz: önce denetim (`--eps 0 --yalniz V1_01_easy --jobs 2`; 12 kayıt ızgaranın N kayıtlarıyla bit-aynı
  olmalı), geçerse `--eps 1e-12 --jobs 2` (18 hücre × A0, A1).
- Erken denetim (ilk hücreler): ızgaranın V1_01_easy R0 N kayıtları (A0, A1, A3; üç pencere) bu belgedeki ön test
  kayıtlarıyla bit-aynı; ön-uç istatistikleri ön testle aynı; hata kaydı yok. Tutmazsa ızgara durdurulur.
- R5 C verisi hazır olunca görüntü düzeyi denetim (`onarim_kos.py olc`): C'nin N'ye göre konumu Sapma 3'teki beklentiyle.
- Karar tabloları (`onarim_analiz.py`) yalnız 45 hücre ve ikiz tamamlanınca üretilir; ara bakış yalnız bütünlük içindir.

**Koşum planı güncellemesi (15:07, başlatmadan önce; yalnız yürütme düzeni, hücreler ve kod aynı).** Takas alanı
dolu (8.15 / 8.19 GB), boş bellek ≈ 11.8 GB; sonda ≈ 20:30'a kadar ≈ 18 GB tutuyor. Izgara bu yüzden iki parçada koşulur:
parça A şimdi (`--yalniz V1_01_easy,MH01,MH02,V1_02_medium,MH03`, 25 hücre, `--gpu_jobs 1 --jobs 4`; ≈ 5.5 saat),
parça B bellek boşalınca (`--yalniz V2_01_easy,MH04,V2_02_medium,MH05`, 20 hücre, daha çok işle). İki parça aynı
`records.jsonl` dosyasına ekler; bitişte `onarim_kos.py plan` ile 45 hücrede kayıt = 9/9 ve bozuk satır = 0 denetlenir.
Sayısal ikiz de parça parça koşulur (denetim yalnız A'da; `ikiz.jsonl`).

### 2d. Görüntü düzeyi denetim, ızgara verisi — altı sekans (17:30; beklenti §1'de önceden yazılıydı)

Kaynak: `goruntu_olcu.jsonl` (30 satır: V1_01_easy, MH02, MH03, MH04, MH01, V1_02_medium × R5 N/O/C, R0 N/O; her 20
karede bir, cam0, temiz gürültüsüz kareye göre). Kalan üç sekans (V2_01_easy, V2_02_medium, MH05) C verisi bitince ölçülür.

| Sekans | R0 RMS akış O/N [0.9, 1.1] | R5 O/N [0.6, 0.95] | R5 C/N [0.75, 1.0] | PSNR R0 O−N | R5 O−N | R5 C−N [−0.5, +1.5] dB |
|---|---|---|---|---|---|---|
| V1_01_easy | 0.898 | 0.751 | 0.852 | −0.22 | +1.35 | +0.38 |
| MH02 | 0.831 | 0.721 | 0.837 | +0.23 | +1.23 | +0.23 |
| MH03 | 0.864 | 0.731 | 0.808 | +0.13 | +0.91 | +0.36 |
| MH04 | 0.864 | 0.655 | 0.772 | +0.10 | +1.15 | +0.33 |
| MH01 | 0.853 | 0.707 | 0.816 | +0.25 | +1.45 | +0.44 |
| V1_02_medium | 0.847 | 0.788 | 0.850 | −0.26 | +0.28 | −0.04 |

- R5 O/N (6/6), R5 C/N (6/6) ve PSNR değişimleri (18/18) beklenen aralıkta. Nedensel kip (C) çevrimdışı kipin (O) kazancının
  bir bölümünü veriyor (C/N 0.77–0.85, O/N 0.66–0.79): beklenen sıra O < C < N altı sekansta da tutuyor.
- **R0 O/N altı sekansın altısında beklenen aralığın altında (0.83–0.90; alt sınır 0.9).** Onarma, yarı-statik rejimde de
  görüntü kaymasını %10–17 azaltıyor; beklentim "≈ 1" idi, yanlış çıktı. Aralık bir beklentiydi, DUR ölçütü değildi; sonuç
  olduğu gibi kaydedildi, aralık değiştirilmedi. Bu, HD kararlarında R0'ın "onarmanın etkisiz olduğu denetim rejimi" diye
  yorumlanamayacağı anlamına gelir; karar tablosunda R0 O/N sonucu bu ölçüyle birlikte okunur.
- Bütünlük: 30/30 hücrede dOrt N'de |·| < 0.4 DN, O ve C'de +1.75…+2.39 DN (bilinen parlaklık kayması); eksik hücre yok.

**2d eki (17:55): kalan üç sekans ölçüldü; `goruntu_olcu.jsonl` 45 satır, 45 tekil hücre (ızgaranın tamamı).**

| Sekans | R0 O/N [0.9, 1.1] | R5 O/N [0.6, 0.95] | R5 C/N [0.75, 1.0] | PSNR R0 O−N | R5 O−N | R5 C−N [−0.5, +1.5] dB |
|---|---|---|---|---|---|---|
| V2_01_easy | 0.866 | 0.765 | 0.860 | −0.24 | +1.33 | +0.36 |
| V2_02_medium | 0.780 | 0.698 | 0.750 | −0.32 | +0.52 | +0.02 |
| MH05 | 0.838 | 0.661 | 0.765 | +0.06 | +1.30 | +0.50 |

Dokuz sekansın özeti: R5 O/N 0.66–0.79 (9/9 aralıkta); R5 C/N 0.750–0.860 (8/9 aralıkta; V2_02_medium aşağıdaki sayıyla
sınırda); PSNR değişimleri 27/27 aralıkta; **R0 O/N 0.78–0.90, 9/9 beklenen aralığın altında** (yukarıdaki not geçerli).
V2_02_medium R5 C/N = 2.8016 / 3.7358 = 0.74994: alt sınırın 0.00006 altında (aralık dışı sayılır, DUR ölçütü değil; değer
olduğu gibi kaydedildi). Yön beklendiği gibi (C, N'den iyi; O, C'den iyi) dokuz sekansın dokuzunda tutuyor.
