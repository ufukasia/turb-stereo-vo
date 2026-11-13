# Faz 2 rakip yöntemler (ORB-SLAM3 stereo, AirSLAM, DROID-SLAM) — ön-testler; beklenen değerler ÖNCE yazıldı

Tarih: 2026-11-08 (başlangıç 12:12, yük ortalaması 1 dk 72.9 / 5 dk 154 / 15 dk 130; 24 çekirdek).
Ön kayıt: `ON_KAYIT_FAZ2.md` v0.3 (sha256 `e119c0f7…be27`, `evidence/on_kayit_sha256.txt` ile aynı; DEĞİŞTİRİLMEDİ).
İlgili madde §1: "Rakipler: ORB-SLAM3 stereo, AirSLAM, DROID-SLAM. Aynı bozulmuş sekanslar üzerinde çalışır
(`turb/uret_sekans.py`, deterministik, sha256 kayıtlı). D/r₀ ∈ {0, 8, 16, 24}; tohum 0–2." v0.3: her hücrede σ = 1 DN,
gürültü tohumu = hücre tohumu; pencere listesi `evidence/pencereler_faz2_kapsam1m.json` (sha256 `cd788381…9e95`);
tek kart RTX 3090; eşleme (sekans, tohum, pencere).

Bu dosyanın §0–§3 bölümleri **hiçbir rakip ön-test sonucu görülmeden** yazıldı. Sonuçlar §4'e eklenir; §1–§3'te sonradan
değişiklik yapılmaz (gerekirse §4'te "sapma" olarak yazılır).

## 0. Ortam, mevcut altyapı ve sabitler

- Python (üretici, değerlendirme): `turbulans vio/kod/turbvio/bin/python`. Veri kökü: `<data-disk>/sdata`
  (`turb/experiments.py::EUROC_FAZ2`, `faz2_kos.sh`).
- Rakip koşucuları önceki makaleden (pyeco_v2 rakip kıyası) **değiştirilmeden** yeniden kullanılır:
  `<home>/Desktop/makaleler/slam/slam_hakem_kosum/rakip/{orbslam3/run_orbslam3_asl.py, airslam/run_airslam_asl.py,
  droid/run_droid_asl.py}` (sözleşme `rakip/SOZLESME.md`; çıktı gövde çerçevesinde ASL CSV + `.json` yan dosya).
  Rakip kaynak kodlarına ve bu koşuculara dokunulmaz; yeni orkestrasyon `pyeco_turb/turb/rakip_*.py` dosyalarındadır.
- Seçilen yapılandırmalar (önceki makalenin "yazar ayarları" varyantları):
  - **ORB-SLAM3** = `orbslam3_ned` varyantı: `rakipler/ORB_SLAM3` (UZ-SLAMLab v1.0 + yalnız başsız örnek
    `Examples/Stereo/stereo_asl_headless`), System::STEREO, EuRoC.yaml ORB ayarları (nFeatures 1200), döngü kapama
    kapalı (PyEco gibi), kamera hızında besleme (stereo_euroc gibi), `--timeout 30000`.
  - **AirSLAM** = yazarların kodu (sair-lab d1265771) + `belirlenimci=2` (yalnız PLNet çizgi ağındaki GPU yazma
    yarışının düzeltmesi ve kapanış senkronu; yazarların iki iş parçacıklı hattı; öteki bütün `ek` bayrakları kapalı),
    kanonik derleme `rakipler/airslam_ws_ek11c` (kaynak `AirSLAM_ek11`), docker `xukuanhit/air_slam:v4`,
    `vo_euroc.yaml` (400 nokta), stereo VO (use_imu 0), `--no-scale-thresholds`. Gerekçe: yazarların yapısı aynı
    girdide koşumdan koşuma farklı yörünge veriyor (`rakip/airslam/AIRSLAM_KUSURLAR.md` #2); eşli tasarım için
    belirlenimli sürüm gerekli. `belirlenimci=1` (önceki `sa_taban`) aynı sonucu ~4× yavaş veriyor.
  - **DROID-SLAM** = `droid_cevrimici` varyantı: `rakipler/DROID-SLAM` (2dfd39f), `droid.pth`, stereo, test_euroc.py
    varsayılanları, `--causal-out` (global BA öncesi çevrimiçi anahtar kare yörüngesi).
- **Birincil yörünge = gerçek zamanlı (nedensel)**, PyEco'nun `m_*` ölçütüyle aynı tanım: ORB-SLAM3 `<out>_nedensel.csv`,
  AirSLAM `<out>_airslam/frames_raw_causal.csv`, DROID `--causal-out`. İkincil = son yörünge (ORB-SLAM3
  SaveTrajectoryEuRoC, AirSLAM anahtar kareye göre birleştirilmiş, DROID global BA + doldurma), PyEco `ms_*` karşılığı.
- **Pencere metriği:** PyEco ile birebir aynı fonksiyon `turb/metrics.py::trajectory_metrics(est, gt, t_min, t_max)`
  (pencere içindeki kestirim damgalarında GT enterpolasyonu, SE(3) Umeyama ATE, RPE 1 s, Sim(3) ölçeği, max_err,
  ıraksama > 5 m). Pencere sınırları PyEco ile aynı: `t_min = ts[w]`, `t_max = ts[w+399]`, `ts` = kaynak sekansın cam0
  damgalarının sayısal sıralı listesi (`runner._run_arms`).
- **Pencere kapsamı** `kapsam_p` = pencerenin 400 cam0 damgasından, kestirimin kapsadığı (tam damga ya da aralarındaki
  boşluk ≤ `max_gap` olan iki kestirim arasında kalan) damgaların oranı. `max_gap` = 0.5 s (ORB-SLAM3, AirSLAM: kare
  başına poz), 2.0 s (DROID çevrimiçi: yalnız anahtar kareler; önceki makaledeki `droid_cevrimici` kuralı).
- **Başarısızlık kodu (koşumdan önce sabitlendi):** pencere "başarısız" (`durum_p = basarisiz`) ⇔ pencerede < 5 kestirim
  (yörünge yok) **veya** `kapsam_p < 0.90` (**X = %90**) **veya** yayılım (kestirim/GT ağırlık merkezli RMS yarıçap,
  SE(3) eşleşmiş örneklerde) [0.5, 2.0] dışı **veya** `max_err > 5 m` (PyEco `m_diverged`). Neden alanı ayrıca yazılır.
- Bozulmuş sekanslar: `<data-disk>/sdata_turb_rakip/` (2 TB diskte ~800 GB boş; `/` 136 GB boş, kullanılmaz).
- Kaynak paylaşımı: aynı anda en çok 1 rakip süreci; GPU yalnız RTX 3090 (DROID `CUDA_DEVICE_ORDER=PCI_BUS_ID
  CUDA_VISIBLE_DEVICES=1`; AirSLAM docker `--gpus device=<3090 UUID>`); 4060 Ti'ye dokunulmaz.

## 1. Denetimler ve beklenen değerler

| # | denetim | beklenen (koşumdan önce) | geçme ölçütü |
|---|---|---|---|
| G0a | `uret_sekans.py`'ye sensör gürültüsü eklenmesi — gürültü kapalı bit-aynılık | Değişiklikten ÖNCE ve SONRA üretilen MH02 kare 600:620 (D/r₀ 0 ve 8, tohum 0): bütün PNG'ler, `data.csv`, `sensor.yaml`, GT/IMU kopyaları, `zamanlar*.txt` bayt bayt aynı. `turb_meta.json` yalnız `kod_parmak_izi.dosyalar["uret_sekans.py"]` alanında farklı. Çıktı klasör adı aynı. | fark 0 (meta: yalnız o alan) |
| G0b | gürültü açık (σ=1) — referans uygulamayla bayt-aynılık | Her PNG == `hooks.sensor_noise(σ=0 çıktısı, 1.0, tohum, kare, kamera)` bayt-aynı; ayrıca == `TurbulenceInjector(sensor_noise_dn=1).apply_stereo(rektifiye temiz cam0, cam1, kare)` (gerçek kanca yolu) bayt-aynı, D/r₀ 0 ve 8. Fark istatistiği: std ≈ 1.04 DN (önceki N2: 1.040), ortalama ≈ 0, cam0–cam1 korelasyonu < 0.02. | bayt farkı 0; std 0.95–1.10 |
| G0c | Türbülans parametreleri PyEco Faz 2 hücresiyle aynı | `TurbParams` (üretici: `kod/config_euroc.yaml` + CLI) == PyEco hücresinin etkin `turbulence` bölümünden kurulan `TurbParams` (`runner._build_cfg` + `hooks.params_from_cfg`); D/r₀ ∈ {0, 8, 16, 24}, tohum 0–2 → `p.hash()` aynı. | 12/12 aynı |
| G0d | Üretim belirlenimi | Aynı komut iki kez (ve `--isci 4` ile `--isci 1`) → görüntü dosyalarının sha256'ları aynı (σ=1, D/r₀=8). | fark 0 |
| G0e | Tüketici rektifikasyonu kimlik | Üretilen rektifiye `sensor.yaml` ile üç rakibin kendi rektifikasyonu (ORB-SLAM3: float32, alpha=-1; DROID: alpha=0; AirSLAM: alpha=0, yaml %.12g) birim dönüşüme çok yakın: harita sapması < 0.01 px; remap edilen görüntü bayt-aynı (ya da < %0.1 piksel 1 DN farklı). | sapma < 0.01 px |
| R1 | Kimlik | Her rakip için commit/parmak izi, yapılandırma dosyası sha256, komut satırı kaydedilir (§4.2). | eksiksiz |
| R2 | Pencere mi bütün sekans mı (temiz pencereler, σ=1, D/r₀=0, tohum 0, V1_01 ve MH02) | **(b) önek koşumu** (kare 0 → son pencere sonu, PyEco'nun aynı protokolü: "arka-uç her zaman kare 0'dan başlar", pencere kesilir): üç rakip de kare 0'dan sonraki ilk 60 karede ilklenir (DROID MH02'nin başındaki durağan bölümde gecikebilir; o zaman ilk pencere başlangıcından önce ilklenmesi yeter); 6 pencerenin hepsi `kapsam_p ≥ 0.90`. **(a) pencere-başına kısa sekans** (yalnız o 400 kare): ORB-SLAM3 ve AirSLAM ilk ≤ 5 karede ilklenir (stereo ilklenme tek karede); DROID ≤ 40 karede (15 anahtar karelik ısınma) → DROID'de `kapsam_p` 0.90–1.00. (a)/(b) pencere ATE oranı 0.5–2 (sistemler kaotik, tek koşum). Maliyet: (b) kare sayısı (a)'nın 1.08–1.5 katı, ama süreç sayısı 1/3. | karar kuralı aşağıda |
| R3a | Temiz regresyon (ham V1_01_easy, bütün sekans, özgün görüntü + özgün kalibrasyon = önceki protokol) | Önceki pyeco_v2 sonuçları (bütün sekans, SE(3) ATE): ORB-SLAM3 `orbslam3_ned` son 0.0508 / nedensel 0.0520 m; AirSLAM yazar yapısı (3090) 0.0831, `sa_taban` (4060 Ti) 0.0882 m; DROID `droid_cevrimici` 0.0677 m. Yayımlanmış: ORB-SLAM3 stereo V101 0.035 m (döngü kapama AÇIK; burada kapalı, büyük olması beklenir). Beklenen bantlar: ORB-SLAM3 [0.030, 0.075]; AirSLAM [0.060, 0.115]; DROID çevrimiçi [0.045, 0.095]. AirSLAM için bayt aynılığı beklenmez (önceki koşum farklı kart/derleme); ORB/DROID belirlenimsiz. | üçü de bandın içinde, kırık yok |
| R3b | Rektifiye girdi yolu | Aynı V1_01'in üretici çıktısı (D/r₀=0, σ=0, kare 0:1300) üzerinde nedensel pencere ATE'leri, R3a koşumunun aynı karelerden (0:1300) kesilen nedensel pencere ATE'leriyle aynı mertebede: oran 0.5–2 (yalnız rektifikasyon ve belirlenimsizlik farkı). Poz dönüşümü (rektifiye T_BS) hatası olsaydı Sim(3) ölçeği ≠ 1 ya da sistematik büyük ATE görülürdü: ölçek 0.97–1.03. | oran 0.5–2, ölçek 0.97–1.03 |
| R4 | Duman (her rakip × D/r₀ {0, 8, 24} × tohum 0 × σ=1 × {MH02, V1_01}, önek koşumu, 18 koşum, 54 pencere kaydı) | Çökme yok: 18/18 koşumda `.csv` + `.json` yazılır, `durum` ∈ {ok, ilklenemedi} (D/r₀=24'te ilklenemedi/kayıp olabilir; bu bir sonuçtur, altyapı hatası değil). 54/54 pencere kaydında metrik alanları dolu ya da `basarisiz` + neden. D/r₀=0: pencere ATE ≤ 0.10 m (ORB, DROID), ≤ 0.15 m (AirSLAM); başarısız pencere yok. D/r₀=8: medyan oran 8/0 > 1 (PyEco 2–6×; rakipler için 1.2–10× beklenir; DROID'in yoğun akışı daha dayanıklı olabilir: 1–5×). D/r₀=24: ORB-SLAM3 ve AirSLAM'de medyan oran > 3 ya da ≥ 1/3 pencere başarısız; DROID'de medyan oran > 2. Etki uygulandı: PSNR(üretilen, rektifiye temiz): D/r₀=0+σ=1 ≈ 47.8 dB (20·log10(255/1.04)), D/r₀=8 < 25 dB, D/r₀=24 < D/r₀=8. | çökme 0, boş alan 0, yön denetimleri |
| R5 | Belirlenim (V1_01, D/r₀=8, σ=1, tohum 0; her rakip toplam 3 koşum) | AirSLAM (`belirlenimci=2`): 3 koşumda nedensel yörünge dosyası bayt-aynı (sha256). DROID: CUDA atomik işlemleri nedeniyle bayt-aynı DEĞİL; pencere ATE tekrar CV'si < %10. ORB-SLAM3: çok iş parçacıklı (LocalMapping zamanlaması) → bayt-aynı değil; pencere ATE tekrar CV'si %5–30. | ölçülür; öneri kuralı aşağıda |
| R6 | Süre ve disk | Üretim (D/r₀ > 0, `--isci 4`, yüklü makinede): 0.03–0.1 s/kare; D/r₀=0+σ=1: < 0.02 s/kare. Kare başı disk ≈ 0.30 MB (iki kamera PNG). Rakip süresi (önek ~1300–1800 kare): ORB-SLAM3 ≈ gerçek zaman (n/20 s, +%10–50 yüklü makinede), AirSLAM 15–30 fps + 20–60 s docker/TensorRT açılışı, DROID 8–17 fps + 15–30 s model yükleme. Tam rakip ızgarası (9 sekans × 4 D/r₀ × 3 tohum × 3 rakip, önek koşumu): toplam 108 üretim + 324 rakip koşumu ≈ 10–16 saat (tek rakip süreci; ORB-SLAM3 için tekrar gerekirse +2–5 saat). Disk: hepsi tutulursa ≈ 46 GB; kuyruk her hücreden sonra sildiği için anlık < 2 GB. | ölçülür, §4'te tahmin |
| R7 | Analiz betiği | `turb/rakip_analiz.py` duman kayıtlarından eşli (sekans, tohum, pencere) D/r₀ / temiz oranlarını ve başarısızlık oranlarını JSON olarak yazar; elle hesaplanan bir oranla aynı. | elle denetim eşit |

**R2 karar kuralı (koşumdan önce):**
1. Tercih (b) önek koşumudur, çünkü PyEco'nun ızgarasıyla birebir aynı protokoldür (aynı kareler 0…son pencere sonu,
   ortak ısınma geçmişi, aynı pencere kesimi) ve ön kayıttaki "aynı bozulmuş sekanslar" ifadesini harfiyen karşılar.
2. (b) seçilir ⇔ temiz testte üç rakip × 6 pencerenin hepsinde `durum_p = ok` (kapsam ≥ 0.90, yayılım bandında) VE
   (b)'nin hücre başına rakip süresi (a)'nın üç penceresinin toplamının ≤ 1.5 katı.
3. (b) bu ölçütü geçemez, (a) geçerse (a) seçilir ve sapma olarak yazılır. İkisi de geçemezse başlatılmaz, rapor edilir.
4. (a) yalnız ek bilgi olarak raporlanır (pencere başı soğuk başlangıcın etkisi).

**R5 öneri kuralı (koşumdan önce):** tekrar CV'si c_r (pencere başına 3 koşumun CV'si, pencereler üzerinden medyan).
PyEco'nun temiz tohumlar-arası CV'si ≈ %9 (ON_TEST_FAZ2_EUROC §5.3) vekil c_s olarak alınır.
- c_r ≤ %5 (≤ c_s/2; tekrar gürültüsü varyansı ≤ %25 şişirir) → tohum başına 1 koşum.
- c_r > %5 → tohum başına 3 koşum, pencere ATE'si üç koşumun medyanı (başarısızlık: 3 koşumdan ≥ 2'si başarısızsa).
- Bayt-aynı ise (AirSLAM beklenen) → 1 koşum.

## 2. Koşum düzeni

1. G0a öncesi: değişiklikten ÖNCE `uret_sekans.py` ile referans çıktılar (MH02 600:620, D/r₀ 0 ve 8, tohum 0) →
   scratch. Sonra gürültü eklenir (`--gurultu_dn`, `--gurultu_tohum`; varsayılan 0 = kapalı), G0a–G0e.
2. R3a: ham V1_01 bütün sekans, üç rakip (sırayla; ORB → AirSLAM → DROID).
3. Duman sekansları üretilir (σ=1, tohum 0): V1_01 0:1300 ve MH02 0:1800, D/r₀ {0, 8, 24}; ayrıca V1_01 0:1300 D/r₀=0 σ=0 (R3b).
4. R3b (rektifiye temiz, σ=0) ve R2 (σ=1 D/r₀=0 önek + pencere görünümleri (a)).
5. R4 duman ızgarası (geri kalan hücreler), R5 tekrarları.
6. R6 tahmin, R7 analiz betiği, kuyruk betiği (`turb/rakip_kuyruk.py`; çalıştırılmaz).

## 3. Engelleyici durumlar (tam rakip ızgarası başlatılmaz)

- G0a–G0d'den herhangi birinin bozulması (gürültü uygulaması PyEco ile aynı değil ya da eski çıktıyı değiştiriyor).
- G0e: tüketici rektifikasyonu görüntüyü yeniden örnekliyor (> 0.01 px).
- R3a: bir rakip bandın dışında ya da kırık (yapılandırma/derleme bozulmuş olabilir).
- R4: altyapı kaynaklı çökme (GPU belleği, TensorRT, docker), boş metrik alanı, D/r₀=0'da başarısız pencere.
- R5: AirSLAM belirlenimli değilse (o zaman yapılandırma yanlış).

---

## 4. Sonuçlar (koşumdan SONRA eklendi; §0–§3 değiştirilmedi)

### 4.1 G0a–G0e — sensör gürültüsü ve üretici kapıları (`pyeco_turb/evidence/FAZ2_rakip/on_test/g0.json`; betik `turb/rakip_g0.py`)

`turb/uret_sekans.py`'ye `--gurultu_dn` / `--gurultu_tohum` eklendi (yedek: `uret_sekans.py.yedek_20261108_gurultu_oncesi`,
sha256 `985b9484…c011`). Gürültü, PyEco kancasının kullandığı fonksiyonun kendisiyle (`turb/hooks.py::sensor_noise`)
türbülanstan sonra uygulanıyor; tohum varsayılanı türbülans tohumu (kanca: `noise_seed None → params.seed`); kapalıyken
işçi durumu, klasör adı ve meta alanları eskisiyle aynı (açıkken ada `_n1_ns<tohum>`, metaya `sensor_gurultusu` ve
`hooks.py` parmak izi eklenir).

| # | ölçülen | sonuç |
|---|---|---|
| G0a | MH02 600:620, D/r₀ 0 ve 8: 54'er dosyanın hepsi bayt-aynı, klasör adı aynı; `turb_meta.json`'da yalnız `kod_parmak_izi.dosyalar.uret_sekans.py` farklı | **GEÇTİ** |
| G0b | D/r₀ 0 ve 8, 20 kare × 2 kamera: `hooks.sensor_noise(σ=0 çıktısı)` ile 0 piksel farkı; gerçek kanca yolu `TurbulenceInjector(sensor_noise_dn=1).apply_stereo` ile 0 piksel farkı; fark std 1.0405 / 1.0404 DN, ortalama 0.0002 / 0.0003, cam0–cam1 korelasyonu −0.0005 | **GEÇTİ** |
| G0c | D/r₀ {0, 8, 16, 24} × tohum {0, 1, 2}: üretici `TurbParams.hash()` = PyEco Faz 2 hücresi (`runner._build_cfg` + `hooks.params_from_cfg`, `pyeco_turb/config.yaml`) 12/12; PyEco hücresinde `sensor_noise_dn=1.0`, `noise_seed=None` | **GEÇTİ** |
| G0d | σ=1, D/r₀=8: iki ayrı üretim ve `--isci 4` ↔ `--isci 1`, 40 PNG'nin hepsi aynı sha256 | **GEÇTİ** |
| G0e | Rektifiye sensor.yaml üç rakibin kendi kuralıyla yeniden rektifiye edilince harita sapması ORB-SLAM3 7.6e-6 px, DROID 9.1e-6 px, AirSLAM 9.4e-6 px; remap edilen görüntüde 0 farklı piksel (iki kamera). Yani rakipler bozulmuş görüntüyü yeniden örneklemiyor (DROID'in kendi 320×512'ye küçültmesi hariç, tasarım gereği) | **GEÇTİ** |

### 4.2 R1 — Kimlik (`evidence/FAZ2_rakip/on_test/kimlik.json`; `python -m turb.rakip_kos kimlik`)

| rakip | kod / sürüm | yapılandırma (sha256) | komut |
|---|---|---|---|
| ORB-SLAM3 | `rakipler/ORB_SLAM3` (UZ-SLAMLab v1.0 kopyası, git yok): `libORB_SLAM3.so` `f66b1bc4…f1a2`, `stereo_asl_headless` `4720d029…c5f0`, `Tracking.cc` `3a3654bc…3264` (yazarların 500 ilklenme eşiği), ORBvoc `f8dd027f…3570` | yazar `Examples/Stereo/EuRoC.yaml` `72160296…b442`; koşucu `settings.yaml`'ı sensor.yaml'dan üretir (EuRoC.yaml ORB değerleri, nFeatures 1200, ThDepth 60, `loopClosing: 0`); her koşumun yan dizininde | `turbvio/bin/python rakip/orbslam3/run_orbslam3_asl.py --timeout 3600 --seq S --calib S --out O.csv` (koşucu `1caa3992…b317`) |
| AirSLAM | yazar commit `d1265771`; kaynak `AirSLAM_ek11` git `c633b516`; derleme `airslam_ws_ek11c`, `visual_odometry` `c79c5199…f807`; docker `xukuanhit/air_slam:v4` (`sha256:9df819fd…8e7`); 3090 TensorRT motorları (plnet_s0 `86dd79ad…`, plnet_s1 `f26c8e2b…`, SP+LG `958889f0…`, SP `2af2710b…`) | `configs/visual_odometry/vo_euroc.yaml` `cd3903a4…db80`; yalnız `ek: {belirlenimci: 2}`; kamera yaml'ı koşucu tarafından sensor.yaml'dan (EuRoC profili: derinlik 0.1–10 m, max_y_diff 1) | `turbvio/bin/python rakip/airslam/run_airslam_asl.py --ws airslam_ws_ek11c --no-scale-thresholds --ek belirlenimci=2 --gpu <3090 UUID> --timeout 3600 --seq S --calib S --out O.csv` (koşucu `709a6845…5087`) |
| DROID-SLAM | `rakipler/DROID-SLAM` git `2dfd39f0`, `droid.pth` `46476ef6…4526` | test_euroc.py varsayılanları (koşucu argümanları; stereo, taban ölçeği B/0.1) | `CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 envs/droid/bin/python rakip/droid/run_droid_asl.py --no-viz --timeout 3600 --seq S --calib S --out O_global.csv --causal-out O.csv` (koşucu `1a75fb92…e644`) |

Değerlendirme: `turb/metrics.py` `923c06a6…5405` (PyEco ile aynı dosya). Hat kimlik parmak izi (`kimlik_fp`) her kayıtta.

### 4.3 R3a — Temiz regresyon (ham V1_01_easy, bütün sekans, 2912 kare; özgün görüntü + özgün kalibrasyon)

Bütün sekans ATE önceki makalenin değerlendiricisiyle (`rakip/degerlendir.py`, GT damgalarına enterpolasyon, SE(3)):

| rakip | beklenen bant | önceki (pyeco_v2) | ölçülen son / nedensel [m] | Sim(3) ölçek | kapsam | süre (yük 1 dk) | sonuç |
|---|---|---|---|---|---|---|---|
| ORB-SLAM3 | 0.030–0.075 | 0.0508 / 0.0520 | **0.0367 / 0.0392** | 1.008 | 1.000, tek harita, kayıp 0 | 616 s = 4.7 fps (50→158) | **GEÇTİ** |
| AirSLAM (belirlenimci=2) | 0.060–0.115 | 0.0831 (yazar yapısı, 3090) | **0.0815 / 0.0821** | 1.014 | 1.000 | 246 s; döngü 187 s = 15.6 fps (158) | **GEÇTİ** |
| DROID-SLAM | 0.045–0.095 (çevrimiçi) | 0.0677 (çevrimiçi) | global **0.0404** / çevrimiçi **0.0677** | 1.009 / 1.018 | 1.000 (global), 586 anahtar kare pozu | 481 s = 6.3 fps, tepe GPU 11.5 GB (191) | **GEÇTİ** |

Aynı koşumdan PyEco protokolüyle kesilen nedensel pencere ATE'leri (pencere 100 / 500 / 900): ORB-SLAM3 0.0375 / 0.0256 /
0.0154; AirSLAM 0.0262 / 0.0182 / 0.0176; DROID 0.0272 / 0.0206 / 0.0189 m (DROID w100 kapsam 0.980). PyEco A0 aynı
pencerelerde (gürültüsüz, 4060 Ti) 0.0288 / 0.0271 / 0.0185 m — üç rakip PyEco ile aynı mertebede.

Not (yük): makinede AirSim ana ızgarası ve EuRoC ön-testleri koşarken yük ortalaması 50–200 arasında dalgalandı;
ORB-SLAM3 bu yükte gerçek zamanın ~¼'ünde (4.7 fps) çalıştı (kamera hızında besleme yalnız beklemeyi ekler, kare düşürmez).

### 4.4 SAPMA 1 (engelleyici bulgu, düzeltildi): üreticinin stereo dış parametresi PyEco'nunkiyle aynı değildi

**Belirti (ilk duman koşumları, 12:49–13:00):** üretilmiş V1_01 (σ=1, D/r₀=0, kare 0:1300) üzerinde ORB-SLAM3 temiz
pencere ATE'si 0.45 / 0.50 / 0.31 m (ham girdide aynı pencereler 0.038 / 0.026 / 0.015 m), Sim(3) ölçeği 1.12–1.69,
yayılım 0.50–0.87, 1300 karede 509 anahtar kare (ham: 2912 karede 88). AirSLAM aynı sekansta **ilklenemedi**.

**Kök neden (ölçüldü):** `turb/uret_sekans.py::Rektifikator`, cam1/sensor.yaml'da `T_cn_cnm1` yoksa eski `kod/` ön-ucunun
**gömülü** cam0→cam1 matrisini (`EUROC_T_C1_C0`) kullanıyordu. Bu matris EuRoC'un T_BS'lerinden türetilen dış parametreden
x ekseni etrafında ~0.8° (0.0141 rad) farklı. Rektifiye görüntülerde SIFT eşleşmelerinin dikey kayması: gömülü matrisle
medyan **−6.7 px** (V1_01) / **−7.3 px** (MH02), |dy| ≤ 1 px oranı %0.1 / %0; T_BS'den türetilen matrisle medyan 0.05 / −0.05 px,
|dy| ≤ 1 px oranı %97 / %96. AirSLAM stereo kabulü |dy| ≤ 1 px (max_y_diff 1) → hiç stereo nokta yok → ilklenemedi; ORB-SLAM3
stereo araması ±2 px'lik satır bandında → stereo derinlik çok az ve yanlı → ölçek kayması.
**PyEco Faz 2 ızgarası bu matrisi kullanmıyor:** `pyeco_turb/config.yaml` (ve `config.yaml.pyeco_orijinal`) `frontend.
stereo_extrinsics: sensor_yaml` → `stereo_front_end.py` T_c1_c0 = inv(T_B_c1)·T_B_c0. Yani üreticinin "kanca eşdeğerliği"
eski `kod/` ön-ucuna (gömülü matris) göreydi; Faz 2 PyEco'sunun gördüğü görüntüyle aynı değildi. Ön kayıttaki "aynı bozulmuş
sekanslar" koşulu bu hâliyle sağlanmıyordu.

**Düzeltme:** `uret_sekans.py`'ye `--stereo_dis {gomulu, sensor_yaml}` eklendi (varsayılan `gomulu` = eski davranış, G0a
bit-aynılığı korunur; `sensor_yaml` seçilince klasör adına `_sy` eklenir). Rakip kuyruğu ve ön-testler `sensor_yaml` kullanır.
Yedek: `uret_sekans.py.yedek_20261108_stereo_dis_oncesi`. Gömülü matrisle üretilen ön-test sekansları silindi, kayıtları
`on_test_gomulu_GECERSIZ.jsonl` / `uret_duman_gomulu_GECERSIZ.log` / `on_test_zincir_gomulu_GECERSIZ.log` adıyla geçersiz
olarak saklandı (analize girmez).

**Yeni kapı G0f (beklenen: PyEco ön-ucuyla bayt-aynılık; koşumdan önce yazıldı, sonuç `g0.json`):**

| denetim | ölçülen | sonuç |
|---|---|---|
| σ=0, D/r₀=0 çıktısı == PyEco `StereoMatcher.load_calibration` (`stereo_extrinsics: sensor_yaml`) + `rectify_images` (MH02 600:620, 2 kamera) | 0 piksel farkı | **GEÇTİ** |
| σ=1, D/r₀ 0 ve 8 çıktısı == `TurbulenceInjector(sensor_noise_dn=1).apply_stereo(PyEco rektifiye kare)` | 0 / 0 piksel farkı | **GEÇTİ** |
| dikey kayma (SIFT, 3 kare) | sensor_yaml: medyan |dy| 0.13 px (%96 ≤ 1 px); gömülü: 7.34 px (%0) | **GEÇTİ** |
| G0a–G0e yeniden (yeni kodla): gömülü varsayılan eski çıktıyla bayt-aynı; G0e sensor_yaml çıktısında da harita sapması ≤ 2.3e-5 px, 0 farklı piksel | hepsi | **GEÇTİ** |

Not (koordinatöre): `turb/test_uret_sekans.py`'deki "kanca eşdeğerliği" testi `kod/config_euroc.yaml` (gömülü matris) ile
yazılmış; Faz 2 PyEco'su için geçerli eşdeğerlik G0f'tir. EuRoC için üreticiyi kullanan başka bir iş varsa
`--stereo_dis sensor_yaml` ile yeniden üretilmeli (AirSim etkilenmez: `T_cn_cnm1` var, kip `kimlik`).

### 4.5 SAPMA 2 (türbülanslı hiçbir sonuç görülmeden, 13:25): DROID çevrimiçi yörüngesinde kapsam tanımı

**Belirti:** temiz MH02 (σ=1, D/r₀=0) önek koşumunda DROID'in pencere 600'ü §0'daki kurala göre "başarısız" çıktı:
`kapsam_p` = 0.560 (< 0.90), ATE 0.0077 m (iyi). Pencerenin ilk 8.8 s'sinde hiç anahtar kare yok (zaman kapsamı 0.551);
öteki boşluklar ≤ 1.25 s. MH02/600 penceresinin başında araç neredeyse durağan (pencere GT yolu 5.03 m, eşiğin hemen üstü).

**Neden kural yanlıştı:** DROID'in çevrimiçi çıktısı yalnız anahtar karelerdir ve DROID, son anahtar kareye ortalama optik akış
eşiğin altındayken (kamera durağan) yeni kare eklemez (`MotionFilter` 2.4 px, `keyframe_thresh` 3.0). Bu bir izleme kaybı değil,
çıktı biçimidir. §0'da "önceki makalenin `droid_cevrimici` kuralı" diye 2.0 s boşluk sınırı yazıldı; önceki makalede bu sınır yalnız
ATE enterpolasyonu içindi, kırık ölçütü zaman kapsamıydı (aktarım hatası).

**Yeni kural (yalnız DROID; ORB-SLAM3 ve AirSLAM değişmedi):** `kapsam_p` sıfırıncı derece tutmayla: pencere karesi, kendinden
önce bir anahtar kare varsa kapsanır (koşum yarıda kaldıysa son anahtar kareden sonraki kareler kapsanmaz). Iraksama
yayılım [0.5, 2] ve max_err > 5 m ile, çökme / ilklenememe yörünge yokluğuyla yakalanmaya devam eder. Eski tanım
`m_kapsam_p_bosluk2s` alanında bilgi için tutulur. ATE yine yalnız anahtar karelerde (PyEco fonksiyonu, değişmedi).
Uygulama: `turb/rakip_kos.py::kapsam_tutma`; bütün ön-test kayıtları `python -m turb.rakip_kos yeniden` ile yeni kurala göre
yeniden değerlendirildi (eski dosyalar `.eski_<zaman>`).

### 4.6 R2 — Pencere protokolü: (a) pencere-başına kısa sekans ↔ (b) önek koşumu (σ=1, D/r₀=0, tohum 0; `sensor_yaml` sekansları)

Nedensel pencere ATE [m] (son kurallarla; `evidence/FAZ2_rakip/on_test/on_test_ozet.json`). (a) = yalnız o pencerenin 400 karesi
(aynı baytlar, görünüm dizini); (b) = kare 0 → son pencere sonu.

| sekans/pencere | ORB-SLAM3 (a) / (b) | AirSLAM (a) / (b) | DROID (a) / (b) |
|---|---|---|---|
| V1_01/100 | 0.0326 / 0.0330 | 0.0261 / 0.0260 | 0.0277 / 0.0277 |
| V1_01/500 | 0.0219 / 0.0256 | 0.0195 / 0.0194 | 0.0204 / 0.0206 |
| V1_01/900 | 0.0183 / 0.0197 | 0.0189 / 0.0186 | 0.0209 / 0.0199 |
| MH02/600 | 0.0076 / 0.0088 | 0.0124 / 0.0160 | 0.0075 / 0.0077 |
| MH02/1000 | 0.0096 / 0.0171 | 0.0116 / 0.0119 | 0.0098 / 0.0124 |
| MH02/1400 | 0.0148 / 0.0165 | 0.0524 / 0.0442 | 0.0119 / 0.0148 |

- İlklenme: (a)'da üç rakip de pencerenin **ilk karesinde** poz veriyor (ORB-SLAM3 ve AirSLAM stereo tek-kare ilklenme;
  DROID ısınma anahtar karelerini ilklenme anındaki pozla yazıyor); (b)'de üçü de kare 0'dan itibaren poz veriyor.
  Kısa sekansta ilklenme sorunu **görülmedi** (beklenen: ORB/AirSLAM ≤ 5 kare, DROID ≤ 40 kare → geçti).
- Kapsam: (b)'de 18/18 pencere `ok` (DROID MH02/600 yalnız §0'daki ilk kapsam tanımıyla "başarısız" idi → SAPMA 2).
  (a)'da 18/18 `ok`.
- (a)/(b) oranı: 0.56–1.18, medyan 0.975 (beklenen bant 0.5–2 → geçti). (a) biraz daha küçük; (b)'de pencere öncesi biriken
  harita/sürüklenme geçmişi var (PyEco'da da aynı).
- Maliyet (rakip süresi toplamı, hücre başına): V1_01 (b) 494 s ↔ (a) 594 s (oran 0.83); MH02 (b) 542 s ↔ (a) 258 s (oran **2.10**);
  iki hücre toplamı 1036 ↔ 852 s (oran 1.22). MH02'de oranı DROID büyütüyor: 1800 karelik önekte izleme 198 s + global BA/doldurma
  118 s (koşucu her koşumda yapıyor; birincil metrik bunu kullanmıyor), pencerede 22 + 14 s; ayrıca koşumlar farklı yük/GPU
  paylaşımı anlarına denk geldi (yük 25–150). Izgara ölçeğinde kare oranı (b)/(a) = 12 917 / 10 800 = 1.20.
- **Karar kuralının uygulanması:** (1) koşulu (18/18 `ok`) sağlandı. (2) maliyet koşulu toplamda (1.22 ≤ 1.5) sağlandı, MH02 hücresinde
  tek başına (2.10) sağlanmadı; kural "hücre başına" ifadesinin ortalama mı tek tek mi olduğunu belirtmiyordu. Izgara toplamında
  beklenen fark ≈ %20–30 (≈ 2–4 saat). **Öneri: (b) önek koşumu** — PyEco ile birebir aynı protokol (aynı kareler, ortak ısınma
  geçmişi, aynı pencere kesimi, ön kayıttaki "aynı bozulmuş sekanslar"), temiz ATE (a) ile aynı mertebede. Koordinatör (a)'yı
  seçerse kuyruk hazır: `python -m turb.rakip_kuyruk --protokol pencere` (yalnız pencere kareleri üretilir, her pencere ayrı koşum).

### 4.7 R3b — Rektifiye girdi yolu (V1_01, D/r₀=0, σ=0, kare 0:1300, `sensor_yaml`) ↔ R3a'nın aynı kareleri

| rakip | pencere 100 / 500 / 900: rektifiye ↔ ham [m] | oran | Sim(3) ölçek (rektifiye) | sonuç |
|---|---|---|---|---|
| ORB-SLAM3 | 0.0245↔0.0375 / 0.0316↔0.0256 / 0.0170↔0.0154 | 0.65 / 1.24 / 1.10 | 1.002 / 1.020 / 1.004 | **GEÇTİ** |
| AirSLAM | 0.0284↔0.0262 / 0.0168↔0.0182 / 0.0148↔0.0176 | 1.08 / 0.92 / 0.84 | 1.010 / 1.010 / 1.005 | **GEÇTİ** |
| DROID | 0.0272↔0.0272 / 0.0206↔0.0206 / 0.0189↔0.0189 | 1.00 / 1.00 / 1.00 | 1.006 / 1.010 / 1.008 | **GEÇTİ** |

DROID'in iki yoldan aynı sonucu vermesi beklenen: koşucusu ham görüntüyü T_BS'den aynı rektifikasyonla düzeltiyor, girdi pratikte aynı.

### 4.8 R4 — Duman ızgarası (3 rakip × D/r₀ {0, 8, 24} × tohum 0 × σ=1 × {V1_01, MH02}, önek; 18 koşum, 54 pencere kaydı)

PSNR (üretilen cam0 ↔ rektifiye temiz, her 50. kare, `etki_psnr.json`): D/r₀=0+σ=1 47.9–48.0 dB (beklenen 47.8); D/r₀=8 V1_01 22.4,
MH02 20.9 dB (< 25; PyEco'nun aynı D/r₀'daki kanca PSNR'ı 22.9 / 20.8 dB); D/r₀=24 19.8 / 16.6 dB (< D/r₀=8) → **etki uygulandı**.

Nedensel pencere ATE [m] (`ok` / **B** = başarısız, neden) ve eşli oran:

| rakip | sekans/pencere | D/r₀=0 | D/r₀=8 | D/r₀=24 | oran 8/0 | oran 24/0 |
|---|---|---|---|---|---|---|
| ORB-SLAM3 | V1_01/100 | 0.0330 | 0.0669 | B (kapsam 0.69) | 2.03 | — |
| | V1_01/500 | 0.0256 | 0.1273 | B (kapsam 0.62) | 4.98 | — |
| | V1_01/900 | 0.0197 | 0.1490 | B (yörünge yok) | 7.55 | — |
| | MH02/600 | 0.0088 | 0.0204 | 0.3802 | 2.32 | 43.4 |
| | MH02/1000 | 0.0171 | 0.0481 | 0.5839 | 2.81 | 34.2 |
| | MH02/1400 | 0.0165 | 0.1002 | 0.8950 | 6.07 | 54.3 |
| AirSLAM | V1_01/100 | 0.0260 | 0.0268 | B (ilklenemedi) | 1.03 | — |
| | V1_01/500 | 0.0194 | 0.0421 | B (ilklenemedi) | 2.17 | — |
| | V1_01/900 | 0.0186 | 0.0459 | B (ilklenemedi) | 2.47 | — |
| | MH02/600 | 0.0160 | 0.0421 | B (ilklenemedi) | 2.64 | — |
| | MH02/1000 | 0.0119 | 0.0752 | B (ilklenemedi) | 6.34 | — |
| | MH02/1400 | 0.0442 | 0.0542 | B (ilklenemedi) | 1.22 | — |
| DROID | V1_01/100 | 0.0277 | 0.0511 | 0.0808 | 1.84 | 2.91 |
| | V1_01/500 | 0.0206 | 0.0231 | 0.1427 | 1.12 | 6.92 |
| | V1_01/900 | 0.0199 | 0.0304 | 0.1463 | 1.53 | 7.34 |
| | MH02/600 | 0.0077 | 0.0370 | 0.1052 | 4.79 | 13.6 |
| | MH02/1000 | 0.0124 | 0.1367 | 0.1149 | 11.1 | 9.30 |
| | MH02/1400 | 0.0148 | 0.1375 | 0.2863 | 9.28 | 19.3 |

- **Çökme yok:** 18/18 koşumda `.csv` + `.json` yazıldı; `durum` = ok (16) / ilklenemedi (2: AirSLAM D/r₀=24, iki sekans).
  AirSLAM D/r₀=24 bir sonuçtur, altyapı hatası değil: 1300/1800 karenin hepsi işlendi, TensorRT hata satırı 0, ilklenme yok
  (yazarların 90 stereo eşleşme eşiği + 1 px dikey tolerans; iki kameranın türbülans eğimi farklı). ORB-SLAM3 V1_01 D/r₀=24:
  ilk OK kare 224, 2 harita, ilklenme sonrası 264 kayıp kare; en büyük harita dışındaki kareler çıktıda yok (SaveTrajectoryEuRoC).
  ORB-SLAM3 MH02 D/r₀=24: tek harita, kapsam 1.0, ama Sim(3) ölçeği 1.35–1.42, yayılım 0.70–0.73, max_err 0.9–1.5 m (ıraksama
  eşiği 5 m'nin altında → "ok" ama çok bozuk).
- **Alanlar:** 54/54 pencere kaydında metrikler dolu ya da `basarisiz` + neden → **GEÇTİ**.
- **D/r₀=0:** 18/18 pencere `ok`, ATE 0.008–0.044 m (≤ 0.10 / 0.15) → **GEÇTİ**.
- **D/r₀=8 medyan oran** (`rakip_analiz.py`, tohum 0): ORB-SLAM3 3.88, AirSLAM 2.32, DROID 3.32 (beklenen 1.2–10; DROID 1–5) → **GEÇTİ**.
  Hiçbir pencere başarısız değil. (PyEco A0 aynı pencerelerde 2–6×.)
- **D/r₀=24:** ORB-SLAM3 3/6 başarısız + kalan 3'te oran 34–54; AirSLAM 6/6 başarısız; DROID 0/6 başarısız, medyan oran 8.3 (> 2) → **GEÇTİ**.
  Bilgi: "başarısızlık dahil" medyan oran ORB-SLAM3 ve AirSLAM için +∞.

### 4.9 R5 — Belirlenim (V1_01, D/r₀=8, σ=1, tohum 0; her rakip 3 koşum)

| rakip | pencere 100 / 500 / 900 ATE, koşum 0 / 1 / 2 [m] | nedensel yörünge sha256 | tekrar CV (pencere) | medyan CV |
|---|---|---|---|---|
| AirSLAM (belirlenimci=2) | 0.0268 ×3 / 0.0421 ×3 / 0.0459 ×3 | 3/3 aynı (`70e397f8…`) | 0 / 0 / 0 | **0 (bayt-aynı)** |
| DROID | 0.0511, 0.0510, 0.0511 / 0.0231 ×3 / 0.0304 ×3 | 3 farklı | 0.12 % / 0.01 % / 0.02 % | **0.02 %** |
| ORB-SLAM3 | 0.0669, 0.0829, 0.0752 / 0.1273, 0.1177, 0.1263 / 0.1490, 0.1029, 0.1197 | 3 farklı | 10.7 % / 4.3 % / 18.9 % | **10.7 %** |

- Beklentiler: AirSLAM bayt-aynı → geçti; DROID bayt-aynı değil, CV < %10 → geçti (CV pratikte 0); ORB-SLAM3 CV %5–30 → geçti.
- **Öneri (§1 kuralıyla):** AirSLAM ve DROID tohum başına **1** koşum; ORB-SLAM3 c_r = %10.7 > %5 → tohum başına **3** koşum, pencere
  değeri üç koşumun medyanı (başarısızlık: 3'ten ≥ 2'si başarısızsa). Kuyrukta `TEKRAR = {orbslam3: 3, airslam: 1, droid: 1}`.
  ORB-SLAM3'ün tekrar değişkenliği D/r₀=8'de PyEco'nun temiz tohumlar-arası CV'sine (%9) eşit; tek koşumla eşli oranın gürültüsü
  ~iki katına çıkardı.

### 4.10 R6 — Süre ve disk (ölçülen; yük ortalaması 1 dk her koşumun başında kayıtlı)

Kare başına rakip süresi [s/kare] (önek koşumları ~1300–1900 kare; docker/model açılışı dahil):

| | hafif yük (yük 13–60) | ağır yük (yük 100–200; AirSim ızgarası + EuRoC ön-testleri) |
|---|---|---|
| ORB-SLAM3 | 0.066–0.081 (≈ gerçek zaman + açılış) | 0.12–0.21 |
| AirSLAM | 0.040–0.061 (D/r₀ ≤ 16); 0.10–0.12 (ilklenemeyen D/r₀=24: her kare stereo) | 0.061–0.095 |
| DROID | 0.079–0.18 (GPU'yu AirSim ızgarası ile paylaşırken üst uç; global BA/doldurma süresin ~%30'u) | 0.16–0.29 |
| Üretim (`--isci 4`) | D/r₀=0: 0.019 ; D/r₀ 8–24: 0.055–0.08 | D/r₀=0: 0.05 ; D/r₀ 8: 0.08–0.12 (isci 3'te 0.30'a kadar) |

Izgara: 9 sekans × 4 D/r₀ × 3 tohum = 108 hücre; önek kare toplamı (D/r₀, tohum) başına 12 917 → rakip tekrarı başına 155 004 kare.
Tekrar: ORB-SLAM3 3, AirSLAM 1, DROID 1 → 540 rakip koşumu.

| senaryo | ORB-SLAM3 (×3) | AirSLAM | DROID | toplam (üretim önden, örtüşür) |
|---|---|---|---|---|
| hafif yük | 9.3 sa | 2.8 sa | 6.0 sa | **≈ 18–19 sa** |
| ağır yük | 22 sa | 3.9 sa | 9.5 sa | **≈ 35 sa** |
| (ORB-SLAM3 ×1 seçilirse) | 3.1 / 7.3 sa | | | ≈ 12 / 21 sa |

Üretim toplamı 3.4–8.6 sa, ama kuyruk bir sonraki hücrenin sekansını rakipler koşarken üretir (hücre başı rakip süresi ≈ 430–1000 s
> üretim ≈ 30–430 s) → toplam süreye neredeyse eklenmez. Uçtan uca kuyruk testinde (§4.11) bir hücre (V1_02, D/r₀=16, tohum 1;
5 koşum) 71 s üretim + 82+81+84 (ORB×3) + 92 (AirSLAM) + 255 (DROID) = 665 s sürdü (yük 13–24).

Disk: üretilen sekans kare başına 0.31–0.43 MB görüntü (σ=1 gürültü PNG sıkıştırmasını büyütür; D/r₀=0 en büyük) + ~10 MB GT/IMU
kopyası → hücre başı 0.4–0.85 GB. Kuyruk hücre bitince siler, bir hücre önden üretir → 2 TB diskte anlık ≤ ~1.7 GB. Hepsi
tutulsaydı ≈ 54 GB. `/` üzerinde yalnız rakip çıktıları (`pyeco_turb/evidence/FAZ2_rakip/kosum`): ORB-SLAM3 ~1.4 MB, AirSLAM ~1.1 MB,
DROID ~0.25 MB / koşum → ≈ 0.6 GB; JSONL ≈ 10 MB. GPU belleği: DROID 5.6 GB (≤ 1500 kare) / 10.2–11.5 GB (MH01, MH02 önekleri,
tampon 1024), AirSLAM 2.5–3.1 GB; kuyruk DROID'den önce 3090'da ≥ 12 GB, AirSLAM'den önce ≥ 3.5 GB boş bellek bekler.

### 4.11 R7 — Analiz betiği ve uçtan uca kuyruk testi

- `turb/rakip_analiz.py` duman kayıtlarında (`on_test_analiz.json`): ORB-SLAM3 V1_01/100 D/r₀=8 eşli oranı 2.2785 = medyan(0.0669,
  0.0829, 0.0752)/0.0330 (elle 2.277; fark yalnız yuvarlamadan) → **GEÇTİ**. Tekrar medyanı, başarısızlık çoğunluk kuralı, pencere ve
  (sekans, tohum) birimi oranları, bootstrap GA, tek yönlü Wilcoxon, sekans başına döküm ve "başarısızlık dahil" medyan JSON'da.
- Uçtan uca kuyruk testi (gerçek uzun-koşum yolu, ayrı kayıt dizini): `python -m turb.rakip_kuyruk --ev evidence/FAZ2_rakip/on_test/kuyruk_testi
  --sekanslar V1_02_medium --dr 16,0 --tohumlar 1 --en_cok_hucre 1 --sil` → manifest yazıldı (kimlik_fp `f578efb6ff85`), sekans üretildi
  (sha256 kaydı `uretim.jsonl`), 5 koşum, 15 pencere kaydı (`uretim_sha256` alanı dahil), ara denetim "TAMAM", önden üretilen ikinci
  hücre (D/r₀=0) bitince silindi, üretim dizini boş kaldı → **GEÇTİ**. (Bu hücre sonucu: V1_02 D/r₀=16'da ORB-SLAM3 3 pencerenin 3'ü
  çoğunlukla başarısız, AirSLAM 3/3 başarısız (geç ilklendi), DROID 3/3 ok; ATE 0.12–0.51 m. Izgara kaydına girmez.)

### 4.12 Ön-test tablosu (beklenen ↔ ölçülen)

| # | beklenen | ölçülen | sonuç |
|---|---|---|---|
| G0a | gürültü kapalıyken eski çıktıyla bayt-aynı | 54/54 dosya aynı, meta yalnız parmak izinde farklı | GEÇTİ |
| G0b | `hooks.sensor_noise` ve gerçek kanca yoluyla bayt-aynı; std ≈ 1.04 DN | 0 piksel farkı; std 1.0404 DN | GEÇTİ |
| G0c | TurbParams PyEco hücresiyle aynı (12 hücre) | 12/12 hash aynı | GEÇTİ |
| G0d | iki üretim ve isci 4 ↔ 1 aynı | 40/40 PNG aynı | GEÇTİ |
| G0e | tüketici rektifikasyonu kimlik (< 0.01 px) | ≤ 2.3e-5 px, 0 piksel farkı | GEÇTİ |
| G0f (sapma 1) | PyEco ön-ucunun (`sensor_yaml`) gördüğü görüntüyle bayt-aynı | 0 piksel farkı (σ 0/1, D/r₀ 0/8); dikey kayma 0.13 px (gömülü: 7.3 px) | GEÇTİ (düzeltmeden sonra) |
| R1 | kimlikler kayıtlı | `kimlik.json`, §4.2 | GEÇTİ |
| R2 | (b) 18/18 ok; (a) ilklenme ≤ 5/40 kare; a/b 0.5–2 | (b) 18/18 ok (sapma 2 ile), (a) 18/18 ok, ilk karede poz; a/b 0.56–1.18 (medyan 0.975); maliyet toplam 1.22×, MH02 hücresi 2.10× | (b) önerildi |
| R3a | ORB 0.030–0.075, AirSLAM 0.060–0.115, DROID 0.045–0.095 | 0.0367, 0.0815, 0.0677 (çevrimiçi) | GEÇTİ |
| R3b | rektifiye/ham 0.5–2; ölçek 0.97–1.03 | 0.65–1.24; ölçek 1.002–1.020 | GEÇTİ |
| R4 | çökme 0; D/r₀=0 ok, ≤ 0.10/0.15 m; 8/0 > 1; 24'te bozulma; PSNR yönleri | 18/18 çıktı; D/r₀=0 18/18 ok (≤ 0.044 m); 8/0 medyan 3.9 / 2.3 / 3.3; 24'te ORB 3/6 + AirSLAM 6/6 başarısız, DROID 8.3×; PSNR 48 → 21–22 → 17–20 dB | GEÇTİ |
| R5 | AirSLAM bayt-aynı; DROID CV < %10; ORB CV %5–30 | bayt-aynı; %0.02; %10.7 | GEÇTİ → ORB×3 |
| R6 | 10–16 sa (ORB tekrarsız varsayımı); disk anlık < 2 GB | 18–19 sa (hafif yük) / ≈ 35 sa (ağır yük) ORB×3 ile; anlık ≤ 1.7 GB | ölçüldü (ORB×3 yüzünden beklenenden uzun) |
| R7 | analiz elle hesapla aynı; uçtan uca kuyruk | aynı; kuyruk testi tamam | GEÇTİ |

### 4.13 Engelleyici durum değerlendirmesi ve başlatma planı

- §3'teki engelleyicilerden yalnız biri görüldü (G0'ın PyEco ile aynı görüntüyü üretmemesi = sapma 1) ve düzeltildi (G0f geçti).
  Sapma 2 (DROID kapsam tanımı) türbülanslı sonuç görülmeden yapıldı. Başka engelleyici yok.
- **Koordinatörün karar vermesi gerekenler:** (i) protokol (b) önek (öneri) / (a) pencere; (ii) ORB-SLAM3 tohum başına 3 tekrar
  (öneri; ~+6–15 sa) ya da 1; (iii) EuRoC tam ızgarasıyla eşzamanlı başlatılırsa ORB-SLAM3 ve DROID süreleri ağır yük satırına yaklaşır
  (sonuçlar etkilenmez: ORB kare düşürmez, AirSLAM bayt-aynı, DROID CV %0.02; ama ORB tekrar değişkenliği yükle değişebilir).
- **Başlatma komutu (çalıştırılmadı):**
  `cd "<project>/pyeco_turb" && ./rakip_kuyruk.sh baslat`
  (önce `./rakip_kuyruk.sh kuru`; koşum sırasında `./rakip_kuyruk.sh denetle`; bitince `./rakip_kuyruk.sh analiz`).
  (a) seçilirse: `CUDA_DEVICE_ORDER=PCI_BUS_ID nohup "$PY" -m turb.rakip_kuyruk --sil --protokol pencere > evidence/FAZ2_rakip/kuyruk_pencere.log 2>&1 &`.
- **Erken denetim (koşum sırasında):** kuyruk 1., 4., 12. hücreden ve her 36 hücrede bir `denetle` basar (tek kimlik_fp, GPU = 3090,
  D/r₀=0'da başarısız pencere yok, temiz ATE ≤ 0.15 m). İlk dört hücre MH01 tohum 0'ın dört D/r₀'ı: D/r₀=0'da üç rakip `ok` ve
  ATE 0.005–0.05 m beklenir; değilse kuyruk durdurulup incelenmeli.

## 5. Koordinatör kararları (2026-11-08 14:39, tam ızgaranın türbülanslı sonuçları görülmeden)
- Protokol: (b) önek koşumu (PyEco protokolüyle aynı). Maliyet kuralı toplam üzerinden yorumlandı (1.22× ≤ 1.5×); MH02 hücresi tek başına 2.10× — bu yorum burada kayıtlıdır.
- ORB-SLAM3: önceden yazılmış kurala uygun olarak tohum başına 3 koşum, medyan.
- Eşzamanlılık: AirSim ana ızgarası (41/120) ile aynı anda başlatıldı; sonuçlar yükten etkilenmez (R5 belirlenim), yalnız süre uzar.
- Başlatma: PID 2908763, manifest evidence/FAZ2_rakip/manifest_20261108_143920.json, kimlik_fp=f578efb6ff85.

## 6. Tam ızgara bitti (2026-11-09 14:51; koordinatör incelemesi)

- 540/540 koşum; 1620 pencere kaydı.
  - kimlik_fp f578efb6ff85 (tek); hepsi RTX 3090; altyapı hatası 0.
- Kayıtlar: `evidence/FAZ2_rakip/rakip_kayitlari.jsonl` (sha256 0393121a…31a3).
- Özet: `rakip_kayitlari_ozet.json` (sha256 62514d16…4ee).
- Denetim "SORUN VAR" dedi. Tek tetikleyici: AirSLAM D/r₀=0'da 3 başarısız pencere (R4 kuralı). İnceleme sonucu aşağıda.

**AirSLAM, V2_01_easy, pencere 100, tohum 0/1/2: yöntem davranışı, altyapı hatası değil.**
- Üç pencerede de `kapsam_p` 0.36–0.37, ATE iyi (0.013–0.022 m). Aynı koşumun 500 ve 900 pencereleri `ok`.
- Neden (`airslam.log`): AirSLAM ilklenmek için ≥ 90 iyi stereo nokta ister (`min_init_stereo_feature: 90`,
  yazarların `vo_euroc.yaml` değeri).
  - V2_01'in ilk 354 karesinde en çok 84 nokta bulunuyor; 354 kez "Not enough stereo points to initialize!" yazıyor.
  - İlk poz kare ≈ 354'te (17.7 s), bu yüzden pencere 100'ün yalnız son ≈ %36'sı kapsanıyor.
- ORB-SLAM3 ve DROID aynı pencerede `ok` (kapsam 1.0).
- R4 dumanı yalnız V1_01 ve MH02'deydi, bu yüzden V2_01'in başlangıcı ön testte görülmedi.
- Belirlenimci yapı (`belirlenimci=2`) yeniden koşumda aynı sonucu verir; yeniden koşulmadı.

**Yüksek türbülansta AirSLAM başarısızlığı da aynı kapıdan geliyor.**
- D/r₀=24'te 81/81 pencere "< 5 kestirim": AirSLAM hiç ilklenmiyor (MH01: 1921 karede en çok 59 nokta; V1_01: en çok 54).
- D/r₀=16'da 42/81 başarısız.
- Makalede "izleme kaybı" değil, "ilklenemedi (stereo nokta kapısı, varsayılan eşik 90)" diye yazılacak.

**Analiz kararı:**
- Ön-kayıtlı başarısızlık kodu (§0, koşumdan önce sabitlendi) aynen uygulanır. V2_01/100, AirSLAM için bütün D/r₀'larda başarısız
  sayılır. Eşli oranlar temiz referansı `ok` olan pencerelerde hesaplanır; temiz medyan 78/81 pencereden.
- Duyarlılık analizi (sonradan eklendi, açıkça böyle etiketlenecek): V2_01/100'ü üç rakipten ve PyEco'dan birlikte çıkaran
  ortak pencere kümesi. Ana sonuç değişmiyorsa bu bir dipnotla raporlanır.

## 7. Sapma 2 (a): rakipler D/r₀ = 4'te — ön test (2026-11-11, sonuçlardan önce yazıldı)

- **Kimlik:** aynı koşucu ve rakip ikilileri (`kimlik_fp` §6'daki f578efb6ff85 olmalı); ön kayıt sha256 2d7daaa8…1504
  (Sapma 2 sonrası, `evidence/on_kayit_sha256.txt`'ye eklendi).
- **Plan:** `python -m turb.rakip_kuyruk --sil --dr 4 --ev evidence/FAZ2_rakip_ek`; kuru denetim: 27 hücre,
  135 rakip koşumu, kapılar temiz (RTX 3090, çıktı diski 881 GB boş).
- **Duman testi:** ilk hücre (`--en_cok_hucre 1`: MH01, D/r₀ 4, tohum 0; üç rakip, ORB-SLAM3 3 tekrar).
- **Beklenen (sonuçtan önce):** üç rakip de `ok`; her pencerenin ATE'si aynı rakibin aynı (sekans, tohum, pencere)
  D/r₀ 0 ve 8 değerleri arasında, saçılım payıyla [0.5 × ATE₀, 1.5 × ATE₈]:
  - ORB-SLAM3 (w1/w2/w3): 0.003–0.036 / 0.020–0.229 / 0.006–0.178 m.
  - DROID: 0.002–0.069 / 0.013–0.136 / 0.008–0.143 m.
  - AirSLAM: 0.008–0.048 / 0.014–0.225 / 0.010–0.094 m.
  - `kimlik_fp` §6 ile aynı; `uretim_sha256` kayıtlı; üretilen sekans silinmiş.

## 8. Sapma 2 (b): ORB-SLAM3 stereo-ataletsel (SI) — ön test (2026-11-11)

- **İkili:** `ORB_SLAM3/Examples/Stereo-Inertial/stereo_inertial_asl_headless` (kaynak `ORB_SLAM3_ek2`'den aynen
  kopyalandı; CMakeLists'e hedef eklendi, yedek `CMakeLists.txt.yedek_20261111_si`), sha256 09dbd264…; yalnız bu
  hedef derlendi. `lib/libORB_SLAM3.so` (f66b1bc4…) ve `stereo_asl_headless` (4720d029…) derlemeden önce ve sonra
  aynı, yani stereo rakibin kimliği değişmedi.
- **Koşucu:** `run_orbslam3_asl.py --imu` (var olan kip; IMU bloğu Stereo-Inertial/EuRoC.yaml değerleri, gürültü
  imu0/sensor.yaml'dan; T_b_c1 = inv(T_B_I) T_B_C0).
- **Temiz özgün EuRoC (beklenen: tam ATE 2–6 cm, ilkleme başarılı, kare ≥ %95 OK):**

| sekans | tam ATE se3 / sim3 (son) | ölçek | OK kare | pencereler (son) | pencereler (nedensel) |
|---|---|---|---|---|---|
| MH01 | 0.0417 / 0.0219 m | 1.008 | 3676/3682 | 0.010 / 0.029 / 0.009 | **0.302** / 0.039 / 0.016 |
| V2_01 | 0.0344 / 0.0238 m | 1.011 | — | 0.010 / 0.012 / 0.012 | **0.184** / 0.017 / 0.011 |

  → **GEÇTİ** (son yörünge). Nedensel yörünge IMU ilklemesini (VIBA) içeren ilk pencerede bozuk; ilkleme haritayı
  döndürüp ölçekler. Karar (türbülanslı SI sonucu görülmeden): HI1 son yörüngeyle (ms_) değerlendirilir, stereo
  ORB-SLAM3 karşılaştırması da ms_ ile; ON_KAYIT Sapma 2 eki (sha256 `evidence/on_kayit_sha256.txt` son satır).

### 7.1 Duman testi sonucu (MH01, D/r₀ 4, tohum 0) → GEÇTİ; tam koşum başladı

| rakip | w721 | w1121 | w1521 | beklenen aralık |
|---|---|---|---|---|
| ORB-SLAM3 r0/r1/r2 | 0.016 / 0.015 / 0.011 | 0.078 / 0.075 / 0.084 | 0.090 / 0.102 / 0.091 | 0.003–0.036 / 0.020–0.229 / 0.006–0.178 |
| DROID | 0.013 | 0.036 | 0.030 | 0.002–0.069 / 0.013–0.136 / 0.008–0.143 |
| AirSLAM | 0.027 | 0.101 | 0.060 | 0.008–0.048 / 0.014–0.225 / 0.010–0.094 |

- 15/15 pencere `ok`, kapsam 1.000; `kimlik_fp` f578efb6ff85 (§6 ile aynı); üretim sha256 2418dc49…; denetim TAMAM.
- Tam koşum: aynı komut `--en_cok_hucre` olmadan (devam kaydı ilk hücreyi atlar), 26 hücre kaldı.
- **Koşucu entegrasyonu:** `rakip_kos.py` / `rakip_kuyruk.py`'ye `orbslam3_si` eklendi (yedekler `*.yedek_20261111_si`);
  ön kayıtlı üç rakibin komutu ve listesi değişmedi (yalnız `--rakipler orbslam3_si` ile seçilir). Modül sha256'sı
  kimliğe girdiği için SI kayıtlarının `kimlik_fp`'si yeni değer alır; D/r₀ = 4 koşumu başlangıçta sabitlenen
  f578efb6ff85 ile sürer.
- **Izgara ve paralellik:** D/r₀ {0, 4, 8, 16, 24} × 9 sekans × tohum 0–2 × 3 tekrar = 405 koşum; her tohum ayrı
  kuyruk örneği ve kayıt dizini (`evidence/FAZ2_rakip_si_t{0,1,2}/`; kilit dizin başına).
- **Kuyruk duman testi (MH01, D/r₀ 0, tohum 0) — beklenen (sonuçtan önce):** 3 tekrarın hepsi `ok`; son yörünge (ms_)
  pencere ATE'leri özgün temiz SI ön testiyle aynı düzeyde, [0.5×, 3×] payla: w721 0.005–0.030, w1121 0.015–0.090,
  w1521 0.005–0.027 m; nedensel w721 bozuk olabilir (VIBA, Sapma 2 eki).

### 8.1 SI kuyruk duman testi (MH01, D/r₀ 0, tohum 0) → GEÇTİ; ızgara başladı

| tekrar | son yörünge w721 / w1121 / w1521 (m) | nedensel w721 / w1121 / w1521 (m) |
|---|---|---|
| r0 | 0.0101 / 0.0327 / 0.0089 | 0.0191 / 0.0417 / 0.0157 |
| r1 | 0.0125 / 0.0274 / 0.0090 | 0.0405 / 0.0362 / 0.0165 |
| r2 | 0.0114 / 0.0288 / 0.0096 | 0.0170 / 0.0370 / 0.0169 |

- Beklenen (son): 0.005–0.030 / 0.015–0.090 / 0.005–0.027 → hepsi içinde; 9/9 `ok`, kimlik_fp 2e0701051eb2; denetim
  TAMAM. Bu üretilmiş sekansta nedensel w721 de bozuk değil (önek koşumunda IMU ilklemesi pencereden önce bitiyor).
- Tam ızgara: üç kuyruk örneği (tohum 0/1/2), her biri 45 hücre × 3 tekrar.

### 8.2 Makine yükü denetimi (2026-11-11 11:20; koşum sırasında, sonuç analizinden önce)

- Ek koşumlar eski rakip ızgarasından farklı yükte koşuyor: eski ızgara yük medyanı 47 (p90 68; EuRoC/H6 ile eş
  zamanlı), yeni D/r₀ = 4 ve SI koşumları ≈ 10–13 (ilk hücreler; dört kuyruk eş zamanlı olduğunda ~35'e çıkabilir).
- ORB-SLAM3 kamera hızında beslendiği için duyarlılık eski veride ölçüldü: aynı (sekans, D/r₀, tohum, pencere)'nin
  3 tekrarı farklı yüklerde koştu; grup içi yük – log ATE ilişkisi Spearman ρ = 0.034 (p = 0.38, n = 642); eğim
  yük +10 → ATE ×0.997 (%95 GA 0.989–1.005); D/r₀ başına 0.991 / 1.024 / 0.999 / 1.011 (hepsi p ≥ 0.10).
- Sonuç: yük farkının ATE etkisi ≲ %1–3, türbülans etkisinin (×2–10) çok altında. Koşumlara müdahale edilmedi;
  analizde yük ortak değişken olarak raporlanır.
- İlk türbülanslı SI hücreleri (MH01, D/r₀ 4, tohum 1–2): son yörünge 0.013–0.067 m, hepsi `ok`; bir tekrarda nedensel
  yörünge 9.01 m (VIBA artefaktı; Sapma 2 ekinin öngördüğü durum).
