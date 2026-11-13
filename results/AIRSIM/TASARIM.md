# AirSim (Colosseum) ile Türbülanslı Ortamda Drone Stereo VO — Tasarım Belgesi

**Sürüm:** 1.0 · 2026-11-07
**Kararlar (kullanıcı):** dil Türkçe; hesap yerel makine (RTX 3090 + RTX 4060 Ti, 24 çekirdek, tam performans); TartanAir çıkarıldı; sistem **stereo VO** (IMU yok); kamera **aşağı bakar**; mesafeler (metrik ölçek) doğru olmalı.

---

## 1. Amaç ve senaryo (v1.1 — yangın senaryosu çıkarıldı, 2026-11-07)

**Karar:** Yangın senaryosu çıkarıldı. Duman Cn² değerleri için güvenilir ölçüm bulunamadı ve senaryo niş kalıyor. Senaryo literatüre ve fiziğe göre seçildi.

### 1.1 Fiziksel mantık: "shower-curtain" etkisi
Gündüz yüzey tabakasında Cn² ∝ h^(-4/3) (Wyngaard vd. JOSA 1971). Sıcak hava yere yakın birikir. Küresel dalga ağırlığı (z/L)^(5/3) ile, nesneye (zemine) yakın türbülans görüntüyü en az bozar.

Bu yüzden **aşağı bakan drone kamerası türbülansı neredeyse görmez:**
- Tilt ≤ 0.003 px.
- Drone→yer r₀ değeri, yer→drone r₀ değerinin 3.6–5.4 katı.

Piksel ölçeğinde etki için üç koşul birlikte gerekir:
1. Uzun görüş hattı (100 m – km),
2. Sıcak zemine **yakından sıyıran** bakış,
3. Cn²(2 m) ≥ 1e-13…1e-12 (çöl/asfalt gündüz; Walters & Kunkel JOSA 1981: gündüz r₀ ≈ 34 mm).

### 1.2 Senaryolar (rejim haritası)

| Kod | Platform / bakış | Rol | Beklenen tilt |
|---|---|---|---|
| **N** | 15–30 m irtifa, **nadir** stereo | Negatif kontrol + temiz SLAM doğrulaması | ≤ 0.003 px (etkisiz) |
| **O** | 3–10 m irtifa, **eğik −20…−35°** stereo | **Ana senaryo**: alt yarı yakın zemin (stereo derinlik iyi), üst yarı uzak (türbülans etkili) | 0.02–0.5 px, derinliğe bağlı |
| **F** | 2–5 m irtifa, ileri (−5°), dar FOV | Üst sınır / yer aracı benzeri | 0.1–2 px |

Senaryo O, kullanıcının "uzak noktalar sorunu" kaygısıyla fiziği uzlaştırır. VO'yu taşıyan yakın zemin noktaları temiz kalır; türbülans uzak noktaları ve stereo derinliğini bozar. Hardie tipi yola bağlı modelin gerekliliği tam olarak budur.

### 1.3 Ölçülecekler
1. Türbülansın stereo VO'yu **hangi rejimde** bozduğu (N/O/F × Cn² × irtifa) — rejim haritası.
2. Bozulmanın **hangi bileşenden** ve **hangi derinlik aralığından** geldiği (tilt/blur; yakın/uzak noktalar).
3. Derinliğe bağlı (öznitelik-başına) ağırlıklandırmanın skaler ve jenerik uyarlamalı kestiricilere katkısı. Türbülans burada derinlikle heterojen olduğu için C3'ün yeniden sınanması anlamlı.

**Literatür boşluğu:** Türbülans altında VO/SLAM/stereo için doğrudan çalışma bulunamadı. En yakın çalışmalar: Oreifej vd. TPAMI 2013 (türbülanslı videoda stabilizasyon) ve BRIAR (UAV/uzun menzil türbülans kovaryatı, 1 km).

**Not:** Aşağıdaki §5.4'teki yangın Cn² tablosu artık kullanılmıyor. Senaryo O/F için hesaplar `scratchpad/tilt.py` ve bu belgenin bir sonraki sürümünde verilecek.

---

## 2. Simülatör seçimi — Colosseum (AirSim halefi) + UE 5.6.1

| Seçenek | Durum | Karar |
|---|---|---|
| Microsoft AirSim | Arşivlendi, yalnız UE 4.27 | ✗ |
| Cosys-AirSim 3.3 | UE 5.5 ister (makinede 5.6.1 var) | Yedek |
| **Colosseum** | Makinede UE 5.6.1 ile derlenmiş paket var (`BlocksV2_Linux`). 2026-09-05 logu ve bugünkü test: Ubuntu 24.04 + Vulkan + RTX 3090 üzerinde çalışıyor | **✓ Ana** |

**Duman testi (2026-11-07) — geçti:**
- Ortam `~/sim_envs/BlocksV2_Linux` altına kopyalandı. NTFS diskte çalıştırma izni olmadığı için yerel diske alındı.
- Başlatıcı: `airsim_sim/sim_baslat.sh <ortam.sh> <gpu>`. Görüntüsüz (`-RenderOffscreen`), `-graphicsadapter=1` (3090).
- Python istemcisi: Colosseum'un kendi istemcisi (`airsim_sim/lib/airsim`), conda ortamı `colosseum` (py3.10).
  - Microsoft istemcisi veri tipi uyumsuzluğu verdi; kullanılmamalı.
- Stereo (cam0/cam1) + DepthPlanar + kamera pozu alındı.
- **Pozlama sorunu:** UE5'in fiziksel aydınlatmasında görüntüler bembeyaz çıkıyor. `simRunConsoleCommand("r.DefaultFeature.AutoExposure.Bias -10")` ile düzeldi. Kalıcı çözüm: kendi seviyemizde manuel pozlamalı Post Process Volume.
- **Kamera yönelimi doğrulandı:** settings'te `Pitch: -90` → API kamera kuaterniyonu (w=0.707, y=−0.707), yani nadir (tam aşağı).

---

## 3. Ortam (sahne)

| | E-A: Cesium + Google Photorealistic 3D Tiles (mevcut paket, Manhattan) | **E-B: Kendi UE 5.6 seviyemiz** |
|---|---|---|
| Doku gerçekçiliği | Gerçek fotogrametri dokusu (VO için iyi) | Megascans/Fab varlıkları (foto-taramalı) |
| Kontrol | Yok; internetten akış, ışık dokuya gömülü | Tam: yangın yerleri, rüzgâr, güneş, manuel pozlama |
| Bilinen sorun | Aşağı bakan yakalama kameralarında zemin parçaları yüklenmiyor (Cesium yalnız ana kamerayı izliyor; derinlik %95 = ∞). Drone gökdelen tepesinde doğuyor. | Colosseum eklentisi UE 5.6 için yeniden derlenmeli (kaynak `~/Desktop/Colosseum` silinmiş) |
| Lisans | Google Map Tiles koşulları türetilmiş veri setinin dağıtımını kısıtlar → **açık benchmark olarak yayımlanamaz** | Fab/Megascans: UE ile üretilen render'ların yayımlanması serbest |
| Kullanım | Yalnız iç ön-doğrulama | **Ana deney ortamı** |

**E-B sahne tasarımı:**
- 600 × 600 m kırsal/sanayi alanı: tarla, orman kenarı, toprak yol, depo çatıları, asfalt saha.
- Sahne zengin ama tekrarsız dokulu olmalı (TartanAir dersi: sentetik dokuda iz ömrü ≈ 2 çıktı).
- 2–4 yangın odağı: Niagara alev + duman. Duman görsel opaklığı ayrı bir eksen olarak başlangıçta **kapalı** tutulur, böylece türbülans etkisi izole edilir.
- Sabit güneş, manuel pozlama.

**Ortam kapısı** (koşumlara geçmeden önce, temiz veride):
- 9-kare iz hayatta kalma > %35
- Temiz ATE/yol < %1
- Sim(3) ölçek hatası < %2

Kapıyı geçmeyen sahne kullanılmaz.

---

## 4. Kamera ve uçuş geometrisi

### 4.1 Neden aşağı bakan (nadir)

Kullanıcının sezgisi doğru: ufka bakan kamerada noktaların çoğu 100 m+ uzakta kalır, stereo derinlik çöker. Nadir kamerada derinlik ≈ irtifa olur ve görüntü boyunca homojendir.

İkincil ablasyon: **oblik −45°** (öne-aşağı). Önü görerek yol bulma senaryosu için.

### 4.2 Stereo geometri — mesafe doğruluğu için sayısal tasarım

Derinlik hatası δZ/Z = Z·δd / (f·B). Alt-piksel eşleşmede δd ≈ 0.25 px alındı.
- Çözünürlük 1024×768, FOV 70° → f = 731 px.
- Gerçek sensör: piksel 3.45 µm → f = 2.52 mm. f/2 diyafram → **D = 1.26 mm** (türbülans hesabında kullanılan gerçek açıklık).

| İrtifa Z | B = 0.30 m | B = **0.40 m** | Z/B (0.40) |
|---|---|---|---|
| 10 m | %1.14 | %0.85 | 25 |
| 15 m | %1.71 | %1.28 | 38 |
| **20 m** | %2.28 | **%1.71** | 50 |
| 30 m | %3.42 | %2.56 | 75 |

**Seçim:**
- B = 0.40 m (ağır sınıf multirotor kol açıklığıyla uyumlu).
- İrtifa 10–20 m ana aralık; 30 m zorlama koşulu.
- Z/B ≤ 50: ORB-SLAM2/pyslam'in "yakın nokta" ölçütü 35–40·B ile aynı mertebe.

### 4.3 Hat uyarlamaları (EuRoC ayarlarının KITTI'de yaptığı hatayı tekrarlamamak için)
- `max_init_depth`: sabit 19 m kullanılırsa 20 m irtifada **bütün** noktalar kesilir. Taban çizgisi katına çevrilmeli (≈ 75·B = 30 m) ya da ≥ 1.5 × irtifa (`KITTI_KOK_NEDEN.md` dersi).
- `jump_max_trans`: hız/kare hızına göre ölçeklenmeli (5 m/s, 20 Hz → 0.25 m/kare; kapı ≥ 0.6 m).
- Kare hızı 20 Hz; hız 2–5 m/s.
- Bu değişiklikler EuRoC yapılandırmasını bozmaz: ayrı bir `config_airsim.yaml` kullanılır.

### 4.4 Çerçeveler ve EuRoC ASL çıktısı
- AirSim dünya çerçevesi NED; gövde FRD.
- Kamera dış parametresi T_BS:
  - nadir: optik eksen = gövde +z (aşağı), görüntü-x = gövde +y, görüntü-y = gövde −x
  - cam1, cam0'a göre +y'de B kadar
- İç parametreler: fx = fy = 731, cx = 512, cy = 384, bozulmasız; rektifiye çift üretilir.
- Çıktı EuRoC ASL biçiminde: `mav0/cam0`, `cam1`, `state_groundtruth_estimate0`, `sensor.yaml`. Hattımız değişmeden okur.

**Zorunlu doğrulamalar** (TartanAir/KITTI'de yaşanan kusurlara karşı, her kayıtta otomatik):
1. Epipolar test: GT pozlarıyla dönme hatası < 0.3°, öteleme kosinüsü > 0.999.
2. Taban çizgisi: rektifiye P2'den −P2[0,3]/fx = B (±1e-4).
3. Derinlik–disparite: render derinliğinden beklenen disparite ile eşleşme disparitesi arasındaki medyan hata < 0.3 px.
4. Zaman damgaları sabit 19 haneli ns; sayısal sıralama = sözlüksel sıralama.

### 4.5 Uçuş senaryoları
Her senaryo 90–120 s, 6 tohum (rüzgâr ve duman gerçeklemesi). Tohum değişse de temiz render aynıdır.

| Kod | Senaryo | Amaç |
|---|---|---|
| S1 | Biçerdöver taraması (lawnmower), 15 m | Yangın alanı haritalama — ana senaryo |
| S2 | Yangın çevresinde orbit, 20 m | Duman kenarında heterojen türbülans |
| S3 | Duman sütununun üstünden geçiş, 10→20 m | Türbülans kameraya yakın — en zor |
| S4 | İrtifa merdiveni 10/15/20/30 m | Derinlik ve yol boyu ölçeklemesi |
| S5 | Oblik −45° ile ileri uçuş | Kamera yönelimi ablasyonu |

### 4.6 Senkron kayıt protokolü
Deterministik ve senkron kayıt için gerçek zamanlı kayıt yerine şu döngü kullanılır:

`simPause(True)` → her adımda `simContinueForTime(0.05)` → tek `simGetImages` çağrısında cam0 + cam1 + cam0 derinlik → aynı anda `simGetVehiclePose`

Hareket API'si: `moveOnPathAsync` (yol noktaları) ya da yörüngeyi önceden hesaplayıp `simSetVehiclePose` ile sürmek. İkincisi tam tekrarlanabilir; GT tanım gereği tam olur. Ana yol bu.

---

## 5. Türbülans fiziği — derinliğe bağlı, çok katmanlı model

### 5.1 İlke
Unreal'ın "heat haze" efekti Kolmogorov istatistiği taşımaz; kullanılmaz. AirSim **temiz** görüntü + **piksel başına derinlik** verir. Türbülansı `turb/` modülünün yeni sürümü sonradan uygular.

Bir temiz render'dan çok sayıda türbülans koşulu üretilir. Render bir kez yapılır; türbülans koşulları CPU/GPU'da paralel uygulanır. Bu, iki GPU'nun verimli kullanımı demektir.

### 5.2 Model (Hardie vd. 2017 OE 56(7):071502; Chimitt & Chan 2020 OE 59(8):083101)
- **Faz ekranları yatay katmanlardır** (nadir bakış için doğal): yükseklikler h₁…h_M (M ≈ 8–12; dumanın yakınında sık).
- **Cn²(x, y, h)** iki bileşenden oluşur:
  1. Arka plan sınır tabakası: yere yakın konvektif profil Cn² ∝ h^(−4/3), sıcak yüzey haritasıyla ölçeklenir.
  2. Yangın dumanı: rüzgârla eğilen, yükseldikçe genişleyen 3B Gauss sütunları. Merkezde Cn²_duman, yarıçap ve yükseklikle azalır.
- **Piksel başına ışın:** pikselin derinliğinden yer noktası bulunur; kamera–yer noktası doğrusu her ekranı keser. Kesişim noktasında ekranın açıklık-ortalamalı gradyan tilt'i alınır. Ayak izi küresel dalga için D·(z/L), z kaynaktan (yer noktası) ölçülür. Tilt'ler küresel-dalga ağırlığıyla toplanır.
- Bu yapı **anizoplanatizmi** (komşu piksellerin ışınları ekranları yakın noktalarda keser → korelasyon) ve **derinlik bağımlılığını** (çatıdaki noktanın yolu kısa) fiziksel olarak üretir. EuRoC modelindeki tek katman ve derinlikten bağımsızlık sınırlılığı ortadan kalkar.
- **Kısa-poz bulanıklık:** yerel yol-integralli r₀'dan tilt-çıkarılmış kısa-poz PSF. **Doğru piksel ölçeği**: PSF örneği λ/(4D)·f px olarak yeniden örneklenir (EuRoC'taki B1 kusurunun düzeltilmiş hâli). Karo bazlı uzaysal değişken PSF uygulanır.
- **Zaman:** katman başına donmuş akış (rüzgâr) + dumanda yukarı akış ve AR(1) "kaynama".
- **Sintilasyon:** opsiyonel, başlangıçta kapalı.

### 5.3 Fizik doğrulama kuralları (koşumlardan önce kapı)

| # | Denetim | Beklenen |
|---|---|---|
| S1 | Tek katman, homojen Cn²: tilt varyansı | 0.170 (D/r₀)^(5/3) (λ/D)² ± %10 |
| S2 | Ölçekleme üssü | 5/3 ± 0.2 |
| S3 | Diferansiyel tilt ↔ DIMM (Sarazin–Roddier / Tokovinin 2002) | < %15 sapma |
| S4 | Yol boyu: homojen Cn²'de σ² ∝ L | eğim 1 ± 0.1 |
| S5 | Konum ağırlığı: aynı Cn² kameraya yakın vs zemine yakın | (z/L)^(5/3) ağırlığıyla uyum |
| S6 | PSF: kırınım FWHM ≈ λ/D·f | 0.32 px ± %15 |
| S7 | PSF: kısa-poz genişliği vs λ/r₀ (Fried kısa-poz bağıntısı) | ± %15 |
| S8 | Anizoplanatik korelasyon: piksel ayrımına göre tilt korelasyonu | Kapalı-form tensörle (Basu 2015 / Hardie 2021) |

### 5.4 Beklenen piksel etkileri (gerçek kamera: D = 1.26 mm, f = 731 px)

| Yol-integralli r₀ | D/r₀ | Tilt σ | Bulanıklık ölçeği λ/r₀ |
|---|---|---|---|
| 10 mm | 0.13 | 0.02 px | 0.04 px |
| 3 mm | 0.42 | 0.06 px | 0.13 px |
| 1 mm | 1.26 | 0.16 px | 0.40 px |
| 0.3 mm | 4.2 | 0.44 px | 1.3 px |
| 0.1 mm | 12.6 | 1.09 px | 4.0 px |

Homojen Cn² ve küresel dalga ağırlığıyla hesaplanan r₀:

| Cn² (m^-2/3) | L = 5 m | L = 20 m |
|---|---|---|
| 1e-11 | 15.5 mm | 6.8 mm |
| 1e-10 | 3.9 mm | 1.7 mm |
| 1e-9 | 0.98 mm | 0.43 mm |
| 1e-8 | 0.25 mm | 0.11 mm |
| 1e-7 | 0.06 mm | 0.03 mm |

**Yorum (önsel beklenti):**
- Normal atmosferik koşulda (Cn² ≤ 1e-12) aşağı bakan geniş açılı drone kamerası türbülansı **görmez**. Bu bir bulgudur; makalede açıkça yazılmalı.
- Piksel ölçeğinde etki ancak **dumanın içinde ya da hemen üstünde** (Cn² ~ 1e-9 – 1e-7, yol 5–20 m) oluşur.
- Bulanıklık tilt'ten ≈ 3–4 kat büyüktür. EuRoC'taki "blur baskın" gözlemi fiziksel ölçekte de beklenir, ama artık doğru büyüklükte.

**Açık iş — literatür değeri:** Yangın dumanı Cn² ölçümleri için *Atmospheric and Oceanic Optics* 15(4) 2002 (yanan orman malzemesi üzerinde konvektif duman, optik ve sıcaklık-gradyanı ölçümleri) ve 37(7) 2024 (model yangın üzerinde 10 m'de lidar Cn²) tam metinlerinden sayısal aralık alınacak. Prof. Ata'nın erişimi/uzmanlığı burada önemli.

Alternatif türetme: Cn² = (79·10⁻⁶ P/T²)² C_T², C_T² duman sıcaklık dalgalanması ölçümlerinden. Süpürme aralığı bu değerler gelene kadar 1e-11 – 1e-7.

---

## 6. Deney tasarımı

### 6.1 Eksenler
- Cn²_duman ∈ {0, 1e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7}
- Dumanın kameraya göre konumu: zemin yakını / orta / kamera yakını
- İrtifa (S4)
- Bileşen: {yalnız tilt, yalnız blur, ikisi} — B1 bulgusunun nedensel testi
- Senaryo S1–S5

### 6.2 Kollar (EuRoC denetimi §4 Faz 1 ile aynı aile)
- A0
- A0C (örneklem dışı kalibre)
- **A1-harita:** oracle, **piksel başına** σ haritası. Lokalize dumanda artık anlamlı.
- **A2-harita:** öznitelik-başına oracle σ (heteroskedastik)
- A3 (çevrimiçi r₀)
- A0C+A3
- AG (jenerik artık-tabanlı uyarlama: VB / Sage–Husa ya da artık RMS)

### 6.3 Önsel beklentiler (koşumdan önce `ON_KAYIT_AIRSIM.md`'ye yazılır)
- **H1:** Cn² ≤ 1e-10'da temiz veriye göre ATE farkı istatistiksel olarak ayırt edilemez.
- **H2:** Diz noktası blur'lu koşulda tilt-only'den daha düşük Cn²'de çıkar.
- **H3 (C3'ün ikinci şansı):** Lokalize dumanda öznitelik-başına σ (A2-harita), skaler şişirmeyi (A1-skaler) tutarlılıkta en az 1.2× geçer. EuRoC'ta tilt alanı homojendi ve yapı işe yaramadı. Burada türbülans uzayda **heterojen**; kuram farkın burada çıkmasını öngörür.
- **H4:** A3 ≈ AG (± %20). Tutarsa C4 bir "ölçüm" katkısı olarak kurulur.
- **H5:** Kamera dumana yakınken bozulma, aynı Cn²'nin zemine yakın olduğu duruma göre belirgin biçimde büyüktür (S5 ağırlığının navigasyon karşılığı).

### 6.4 Ölçütler ve istatistik
- ATE (SE3), RPE, **Sim(3) ölçek hatası** (mesafe doğruluğu), ıraksama oranı, ANEES (N bağımsız birim bandı), kapsama.
- Bağımsız birim: (senaryo, tohum).
- D/r₀ ya da Cn² başına ayrı rapor, Holm düzeltmesi, %95 bootstrap GA.
- Bütün sayılar tek betikten (`make_tables.py`) gelir.

---

## 7. Hesap düzeni (yerel, tam performans)

| İş | Kaynak |
|---|---|
| UE render (sahne başına bir kez) | RTX 3090 (`-graphicsadapter=1`) |
| Türbülans uygulaması (ışın izleme + PSF) | CuPy/PyTorch, RTX 4060 Ti (CUDA 0) + CPU havuzu |
| ALIKED + LightGlue ön-uç | İki GPU, iş kuyruğuyla paylaşılır |
| GTSAM arka-uç | 20 CPU işçi (`OMP_NUM_THREADS=1`) |
| Ön-uçuş kapısı | Her kod değişikliğinden sonra (mevcut altyapı) |

---

## 8. İş planı ve kapılar

| Adım | İş | Süre | Kapı |
|---|---|---|---|
| 1 | Colosseum kaynağını UE 5.6 için yeniden kur, eklentiyi derle; boş seviyede API testi | 1–2 gün | API + stereo + derinlik |
| 2 | E-B sahnesini kur (Fab/Megascans, yangın odakları, manuel pozlama); seviye Python ile yerleştirilebilir | 3–5 gün | — |
| 3 | `kayit.py`: senkron kayıt + EuRoC ASL + 4 otomatik doğrulama | 1–2 gün | Epipolar / taban / disparite / sıralama |
| 4 | Ortam kapısı: temiz S1–S5'te iz ömrü, ATE/yol, ölçek | 1 gün | §3 kapısı |
| 5 | `turb/` derinliğe bağlı çok katmanlı modül + S1–S8 | 4–6 gün | S1–S8 yeşil |
| 6 | Önsel beklentileri yaz, küçük ızgarada ön-uçuş | 1 gün | P/E maddeleri |
| 7 | Tam koşum | 2–4 gün makine | — |
| 8 | Makaleye "AirSim yangın senaryosu" bölümü | 3–4 gün | Orkestra denetimi |

Paralelde: EuRoC tarafında PSF düzeltmesi ve kolların yeniden koşumu (denetim raporu Faz 1–2).

## 9. Riskler
1. Sentetik dokuda iz ömrü kısa kalabilir (TartanAir dersi). Önlem: Megascans foto-tarama dokular, §3 kapısı, gerekirse kare hızını 30 Hz'e çıkarmak.
2. Colosseum'un UE 5.6'da yeniden derlenmesi sorun çıkarabilir. Daha önce başarıldığı log ile sabit; olmazsa Cosys-AirSim + UE 5.5.
3. Yangın Cn² değerleri literatürde dağınık. Önlem: süpürme + C_T² üzerinden türetme + Prof. Ata.
4. Hâlâ sentetik: gerçek veri itirazı kalır (denetim raporu Faz 3).

---

## 10. Makale üretim orkestrasıyla uyum (`~/Desktop/skillerim/yzmakaleuretimskill`, `/makale-uret`)

Deneyler bittiğinde makale bu orkestrayla üretilecek. Deney çıktıları şimdiden onun sözleşmelerine göre üretilir:

| Orkestra kuralı | Bizim deney tarafında karşılığı |
|---|---|
| Kanıt-bağlılık: her sayı `evidence/evidence.json`'da id ile (`evidence-contract.md`) | `make_tables.py` tabloların yanında **doğrudan `evidence.json` kayıtları** üretir (id, metrik, yöntem, veri_seti, değer, kaynak_dosya, kaynak_konum, durum) |
| Terk edilen kanıt silinmez, `durum: terk` + paradigma no | PSF düzeltmesi öncesi EuRoC kanıtları (E1a/E1b/E1c/E3/E4/E8/E9, blur'lu) → **P-001 "PSF ölçek kusuru"**; TartanAir kayıtları → **P-002 "sıralama kusuru"**. `notes/paradigma-kaydi.md`'ye yazılır. |
| Figürler kod-tabanlı (TikZ/pgfplots), veriden okunur; raster yasak | Her figür için `figures/data/*.dat` (pgfplots okuyacak düz tablo) üretilir. matplotlib PNG ana figür olmaz. |
| **Nitel gerçek-çıktı paneli zorunlu** (7b, §A2); kare seçimi ölçütle | AirSim eğik senaryo: temiz / türbülanslı kare + öznitelik izleri + derinliğe göre bozulma haritası. Kare seçim ölçütü önceden: ör. "Cn² orta düzeyde, uzak bölgedeki iz kaybı medyana en yakın kare". |
| Sonuç tahlili S/D/I kuralları, kaynaklı (`sonuc-tahlili.md`) | `turb/sanity.py`'ye S6–S8 (PSF, yol boyu, anizoplanatik korelasyon) eklenir; her biri künyeli. |
| Ön-kayıt / ≥3 tohum / beklenti dosyası (`deney-sozlesmesi.md`) | `ON_KAYIT_FAZ2.md` + `ON_KAYIT_AIRSIM.md` koşumdan önce; tohum 6. |
| Dil Türkçe, şablon varsayılan Springer (IEEE istenirse açıkça söylenir) | Q1 hedef dergi seçildiğinde şablon kararı verilecek (IEEE Sensors J. / T-IM ise IEEEtran). |
