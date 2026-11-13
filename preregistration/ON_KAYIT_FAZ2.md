# Faz 2 ön-kayıt (koşumlardan ÖNCE yazıldı) — v0.3, 2026-11-08 (v0.1: 2026-11-07; v0.2 taslağı: 2026-11-08 gece)

Platform: `pyeco_turb` (PyEco 52fe017 + türbülans modülü; E0 bit-özdeşlik kapısı geçti).
Bu dosya dondurulduktan sonra hipotez, öngörü veya çürütme ölçütü **değiştirilmez**; değişiklik
gerekirse altta tarihli "sapma" kaydı açılır ve makalede bildirilir.

---

## 0. Faz 2 öncesi tespit edilen olgu (pilot, hipotez değil)

Temiz MH01'de (k=16 A0C dökümünden A0 kovaryansına geri ölçeklenerek, n=125 anahtar kare örneği):

| | medyan | p90 |
|---|---|---|
| göreli poz hatası, konum | 11.2 mm | 20.5 mm |
| A0'ın beyan ettiği σ, konum | 1.21 mm | 1.69 mm |
| göreli poz hatası, dönme | 0.30° | 0.44° |
| A0'ın beyan ettiği σ, dönme | 0.062° | 0.077° |

Pencere ≈ 8 s, göreli öteleme ≈ 0.27 m. A0 hatayı ~9× (konum) / ~5× (dönme) küçük beyan ediyor
→ NEES ≈ 367 / 150. Buna rağmen gerçek artık RMS (0.76 px) beyan edilen σ₀'dan (1 px) küçük.
İki açıklama birbirinden ayrılamıyor:

- **(i) Kestirici-içi aşırı güven:** gözlemler bağımsız kabul ediliyor. Aynı noktanın kareler
  arası hata korelasyonu ve ön-uç seçim yanlılığı bu varsayımı bozuyor. BA literatüründe bilinen
  bir durum.
- **(ii) Değerlendirme tabanı:** EuRoC MH'de GT, Leica yalnız-konum ile füzyonlanmış yönelimden
  oluşuyor. Buna kamera–gövde dış kalibrasyonu ve zaman senkronu hataları ekleniyor. Bunlar
  mm/0.1° düzeyinde ve A0'ın σ'sından büyük.

**Pilot olgular (v0.3'e eklendi; hipotez değil, H1–H6 metni değiştirilmedi):**

- **ATE kaotik, tek temiz koşum geçerli referans değil.** PyEco aynı GPU'da bit-belirlenimci. Ama RTX 4060 Ti ↔ RTX 3090
  arasında temiz AirSim ATE'si %47 oynadı; yalnız 1 DN sensör gürültüsü tohumu değişince bütün-sekans ATE'si 2.6 kat oynadı.
  EuRoC'ta 1 DN gürültüyle 5 tohumun temiz pencere ATE CV'si %3–24 (medyan %9). Kaynak: `ON_TEST_FAZ2_EUROC.md` §4.11, §5;
  `airsim_sim/africa/ON_TEST_GURULTU_TABANI.md`.
- **EuRoC D/r₀=8 ön-bakışı:** MH02 + V1_01, 3'er pencere, tohum 0–4, eşli ATE oranı 8/0 30/30 çiftte > 1.5 (medyan 5.7).
  H1'in ilk yarısı (D/r₀ ≤ 8'de ≤ 1.2×) büyük olasılıkla çürüyecek. H1 bilerek değiştirilmeden donduruldu.
- **AirSim pilotu (tohum 0–2):** F5 eşli RPE1s oranı 1e-12 → 1.15, 3e-12 → 1.38, 1e-11 → 2.07 (tekdüze); F15'te 1e-12 → 0.93,
  1e-11 → 1.02. §4 v0.3 öngörüleri bu pilottan sonra yazıldı; bu yüzden doğrulayıcı AirSim ızgarası pilotta kullanılmamış
  YENİ tohumlarla koşulur (§4).

**Sonuç:** EuRoC NEES'i mutlak tutarlılık ölçüsü olarak kullanılamaz, yalnız **kollar arası
göreli** karşılaştırma için kullanılır. Mutlak tutarlılık iddiası AirSim'de (kusursuz GT,
kalibrasyon, senkron) ve Vicon odalarında (V1/V2, 6-SD GT) test edilir. Bkz. H3.

---

## 1. Tasarım (sabit)

- **Sekanslar:**
  - EuRoC MH01–MH05: ana analiz.
  - V1_01, V1_02, V2_01, V2_02: Vicon 6-SD GT, NEES dönme analizi ve A0C kalibrasyonu.
- **Pencereler:** 400 kare. GT yol uzunluğu ≥ 5 m olan, örtüşmeyen, sekans başına en fazla
  3 pencere, hareket ölçütüyle seçilir.
  - Eski sabit başlangıçlar (500/900/1300/1600) terk edilir; MH01/500 penceresi 0.075 m ile
    neredeyse hareketsiz (denetim B6).
  - Seçim yalnız GT hareketine bakar, hiçbir koşum sonucuna bakmaz; seçim betiği ve listesi
    koşumdan önce kaydedilir.
- **Pencere kapsam ek ölçütü (v0.3):** GT konum kapsamı ≥ 1 m (yerinde sallama pencereleri dışlanır; MH01/100 ve
  MH02/100). Liste: `evidence/pencereler_faz2_kapsam1m.json` (yalnız GT hareketinden seçildi).
- **Tohumlar:** 0–4 (türbülans alanı). İstatistik birimi (sekans, tohum); pencereler birim
  içinde ortalanır.
- **Sensör gürültüsü (v0.3):** D/r₀=0 dahil bütün hücrelerde σ = 1 DN Gauss sensör gürültüsü (türbülanstan sonra,
  `rint` + kırpma); gürültü tohumu = hücre tohumu. D/r₀=0'da da tohum 0–4 koşulur; temiz referans tek koşum değil,
  aynı tohumlu temiz koşumdur (eşli) ve "temiz medyan" 5 tohumun medyanıdır.
- **Kovaryans örnekleme (v0.3):** `cov_every = 2` (NEES örnek sayısı pencere başına ≈ 40–80).
- **Donanım (v0.3):** bütün EuRoC ve AirSim hücreleri tek kartta, RTX 3090 (`CUDA_DEVICE_ORDER=PCI_BUS_ID`,
  `CUDA_VISIBLE_DEVICES=1`). Her kayıt GPU adını taşır; analiz betiği aile içinde tek GPU olduğunu denetler.
  A0C k kalibrasyonu da aynı kartta, `cov_every=2` ve σ=1 ile yinelenir (ölçüt aynı).
- **D/r₀ ekseni:** {0, 4, 8, 12, 16, 20, 24}. PSF fiziksel ölçekte (`psf_legacy_scale: false`).
- **Kollar:**
  - Hepsinde: A0, A1 (oracle), A3 (DIMM), AG (genel artık-uyarlamalı), A0C+A3.
  - A0C: k, aşağıdaki kalibrasyondan sonra.
  - A5 (oracle yama, yavaş): yalnız MH02 + MH04 alt kümesi.
  - A2 ve A4: kapsam dışı (ek malzeme).
- **Robust rejim:** bütün kollarda PyEco'nun kendi Huber'i (mono 2.4477, stereo 2.7955).
- **A0C k kalibrasyonu (örneklem dışı):**
  - k, yalnız V-serisi temiz (D/r₀=0) koşumlarında seçilir: medyan konum NEES'i χ²₃ medyanına
    (2.37) en yakın k ∈ {1, 2, 4, 8, 12, 16, 24}.
  - MH sekansları kalibrasyonda kullanılmaz.
- **Faktöriyel (tilt × blur):** D/r₀ ∈ {8, 16}; {yok, yalnız tilt, yalnız blur, ikisi}; yalnız A0.
- **Nedensel müdahale:** D/r₀ ∈ {kırılma noktası −4, kırılma, +4}.
  `replenish_ratio_kp` ∈ {0.2 (varsayılan), 0.4}, `max_keypoints` ∈ {256, 512}; yalnız A0.
- **Rakipler:** ORB-SLAM3 stereo, AirSLAM, DROID-SLAM. Aynı bozulmuş sekanslar üzerinde
  çalışır (`turb/uret_sekans.py`, deterministik, sha256 kayıtlı). D/r₀ ∈ {0, 8, 16, 24};
  tohum 0–2.

## 2. Metrikler ve istatistik (sabit)

- **Birincil:** ATE (SE(3) hizalı, gerçek zamanlı yörünge) ve pencere-göreli NEES (konum ve
  dönme ayrı; ANEES).
  - ANEES bandı örneklem sayısına göre kurulur: χ²(3n)/n, %95.
- **İkincil:**
  - RPE.
  - İzleme kaybı / ıraksama oranı.
  - Kare başına eşleşme sayısı.
  - Stereo eşleşme oranı.
- **Karşılaştırma:** her D/r₀ seviyesinde kol − A0, eşli (aynı sekans, tohum, pencere).
  - D/r₀ etkisi (H1, kırılma): türbülanslı / temiz oranı aynı (sekans, tohum, pencere) üzerinden eşli (v0.3).
  - Wilcoxon işaretli sıra testi.
  - Aile içi Holm düzeltmesi: aile = (kollar × D/r₀) aynı metrik için.
  - Etki büyüklüğü: medyan oran, 10 000 örnekli bootstrap %95 GA.
- **"Kırılma noktası":** A0'ın medyan ATE'sinin, D/r₀=0 medyanının 2 katını ilk aştığı D/r₀,
  ya da ıraksama oranının %20'yi ilk aştığı D/r₀; hangisi önce gelirse.

## 3. Hipotezler, öngörüler, çürütme ölçütleri

**H1 — Dayanıklılık ve kırılma.**
- Öngörü: D/r₀ ≤ 8'de A0 ATE'si temizin en fazla 1.2 katı. Kırılma noktası 12–20 aralığında.
- Çürütme:
  - D/r₀ = 8'de oran > 1.5 ise H1'in ilk yarısı çürür.
  - 24'e kadar kırılma yoksa makale "test edilen aralıkta eşik yok" diye raporlar.
  - Eski "D/r₀≈5 eşiği" P-001 ile terk edilmiştir; hangi durumda olursa olsun geri gelmez.

**H2 — Mekanizma blur, tilt değil.**
- Öngörü: D/r₀ = 16'da ATE artışının ≥ %75'i yalnız-blur koşulunda da görülür. Yalnız-tilt,
  artışın ≤ %25'ini açıklar.
- Çürütme: yalnız-tilt payı > %50.

**H3 — Aşırı güvenin kaynağı.**
- Öngörü: AirSim'de temiz veride A0 ANEES'i yine bandın çok üstünde (> 10). Yani aşırı güven
  kestiriciye içkin, GT'den gelmiyor.
- Çürütme: AirSim temiz A0 ANEES'i bant içinde. Bu durumda EuRoC'taki aşırı güven
  GT/kalibrasyon kaynaklıdır ve makalede EuRoC NEES'i yalnız göreli kullanılır.

**H4 — Fizik-tabanlı kollar, genel uyarlamalı kola karşı.** Makalenin ana iddiası bu karşılaştırmaya bağlı.
- Öngörü (dürüst beklenti):
  - D/r₀ ≥ 8'de A3 ve A1, A0'a göre ANEES'i anlamlı düşürür.
  - AG de düşürür; A3 ile AG arasındaki fark %25'ten küçük kalır.
- Çürütme / yeniden çerçeveleme:
  - A3, AG'yi Holm sonrası anlamlı biçimde geçemezse "türbülans fiziği gerekli" iddiası
    düşürülür. Katkı "türbülans altında tutarlılık bozulmasının nicelenmesi ve basit
    uyarlamalı düzeltme" olarak yeniden yazılır.
  - A3 AG'yi geçerse fizik-tabanlı katkı korunur.

**H5 — Doğruluk.**
- Öngörü: hiçbir kol D/r₀ ≤ 8'de ATE'yi A0'a göre > %10 iyileştirmez. Kollar tutarlılığı
  değiştirir, doğruluğu değiştirmez.
- Çürütme: Holm sonrası anlamlı > %10 iyileşme. Bu da raporlanır.

**H6 — Kırılma özellik açlığından.**
- Öngörü: `max_keypoints` 512 ve `replenish_ratio_kp` 0.4 ile kırılma noktası ≥ 4 birim yukarı
  kayar.
- Çürütme: kayma < 2 birim. Bu durumda kırılma, eşleştirme kalitesi kaynaklıdır, sayı
  kaynaklı değildir.

---

## 4. AirSim kolu — v0.3 (2026-11-08; kullanıcı onayı 2026-11-08)

v0.1'deki "eğik geniş açılı kamera ana senaryo" öngörüsü, koşum yapılmadan fizik hesabıyla
çürütüldü: 70° kamerada tilt 0.01–0.07 px. Bkz. `airsim_sim/rejim_hesap.py`. v0.1 metni
`ON_KAYIT_FAZ2.md.yedek_v01`, v0.2 taslağı `ON_KAYIT_FAZ2.md.yedek_v02` dosyasında duruyor.
v0.2 → v0.3 değişikliklerinin gerekçesi pilot ön testlerdir (`airsim_sim/africa/ON_TEST_AIRSIM_PILOT.md`,
`ON_TEST_GURULTU_TABANI.md`).

**Ortam:** Africa_Savannah.
- Su, hayvanlar ve insan figürleri her kayıtta silinir (`kayit.py --temizle`).
- Sahne dondurulur (`--dondur 1`).
- Alan derinliği ve hareket bulanıklığı kapatılır (`--render_temiz 1`).
- Işık gün batımına sabit. Sahne geometri ve doku içindir; türbülans gücü Cn² parametresidir.

**Model:** `pyeco_turb/turb/derinlik_model.py`, doğrulama `turb/test_derinlik_model.py`; 12/12 test geçti.
- 12 kameraya dik faz ekranı; piksel başına hat integrali, DTM'den yükseklik.
- Cn²(h) = C0 (h/2 m)^(-4/3); küresel dalga ağırlığı.
- 4 px üstü tilt yapısı çarpıtma; 4 px altı tilt ile kısa-poz fazlası bulanıklık.
- Donmuş akış, rüzgâr 2 m/s.
- Sensör gürültüsü σ = 1 DN, türbülanstan sonra (`uret_derinlik.py --gurultu 1`), gürültü tohumu = türbülans tohumu.

**Kameralar:** 1024×768, B = 0.40 m, pitch −3°, 5 m AGL arazi takibi, 3 m/s, 20 Hz. Rota 150 m gidiş-dönüş
(`africa/rota_F1.json`). Kayıtlar D3/D4/D5 doğrulamasından geçti.

| Ad | FOV | fx | açıklık D | rol |
|---|---|---|---|---|
| F15 | 15° | 3889 px | 4.8 mm | tele; ana |
| F5 | 5° | 11700 px | 14.4 mm | zoom; üst sınır |

- **N15 (nadir negatif kontrol) düşürüldü (v0.3).** Pilot: veri doğrulaması geçti, ama PyEco temiz veride Sim(3) ölçeğini
  %3.2–3.3 kaçırıyor. Neden: nadir bakışta ≈ 49 px/kare yanal akışta fırsatçı stereo tetikleyicisi stereo eşleşmeyi atlıyor
  (zorlanınca %1.0). Bu türbülanstan bağımsız bir etken. "Türbülans ≈ 0" rolünü F15/F5'in C0 = 1e-13 düzeyi üstlenir.
  HA1 bu nedenle test edilmez.

**Pencereler (v0.3):** her sekansta iki düz şerit, W1 = kare 20:857, W2 = kare 1157:1994 (her biri 41.8 s, ≈ 125 m).
Gerekçe: F5 bütün sekansta U dönüşünde kopuyor (temiz ATE 42 m). Pencereler yalnız rotanın geometrisinden (düz bacaklar)
seçildi, türbülanslı sonuçlara bakılmadı. Her pencere ayrı bir PyEco koşumudur.

**Cn² seviyeleri (C0 = Cn² @ 2 m):** 0, 1e-13, 1e-12, 3e-12, 1e-11 m^(-2/3).
- 1e-14 düşürüldü (v0.3): σ_tilt < 0.01 px, 1e-13 ile aynı rol.
- 1e-11 sıcak çöl öğleni için üst uçtur; literatür aralığı cite1 ile doğrulanacak.

**Tohumlar (v0.3):** doğrulayıcı ızgara pilotta kullanılmamış tohumlarla koşulur: F5 tohum 10–14 (5), F15 tohum 10–16 (7).
Pilot tohumları (0–2) doğrulayıcı analize girmez; pilot olarak ayrıca raporlanır. F15'te daha çok tohum, eşdeğerlik
testinin gücü içindir (pilot log-oran SD ≈ 0.17 → ±%15 sınırında güç 0.8 için ≈ 14 çift).

**Izgara:** F15 5 C0 × 7 tohum × 2 pencere = 70 hücre; F5 5 × 5 × 2 = 50 hücre; toplam 120 hücre, hepsi RTX 3090.

**Modelin önceden hesapladığı büyüklükler** (σ_tilt toplam, baskın derinlik dilimi; SLAM koşumundan ÖNCE):

| C0 | F15 | F5 |
|---|---|---|
| 1e-13 | 0.018 px | ≈ 0.06 px (ölçekten) |
| 1e-12 | 0.065 px | 0.220 px |
| 3e-12 | 0.112 px | 0.381 px |
| 1e-11 | 0.205 px | 0.695 px (D/r0 ≈ 4.9) |

**Metrikler (v0.3):**
- **Birincil:** pencere RPE1s. Yörünge GT'ye SE(3) hizalanır, Δ = 20 kare (1 s) öteleme farkı hatasının RMS'i
  (`airsim_sim/rpe_ozet.py`). Etki ölçüsü: eşli log-oran ln(RPE1s_C0 / RPE1s_temiz), aynı (sekans, pencere, tohum).
- **İkincil:** pencere ATE (SE(3)), RPE10s, Sim(3) ölçeği.
- Gerekçe: pilotta temiz RPE1s CV %4–19, ATE CV %10–61. Eşli RPE1s etkisi GPU değişince korundu (1.87 ↔ 1.81), ATE korunmadı.

**İstatistik:** birim (pencere, tohum) çifti; F15 n = 14, F5 n = 10 her C0 düzeyinde.
- Eşdeğerlik: TOST, sınır ±ln(1.15), α = 0.05 (log-oran ortalamasının %90 GA'sı sınır içinde). t ve bootstrap birlikte
  raporlanır; karar t-TOST ile.
- Bozulma: tek yönlü Wilcoxon işaretli sıra (log-oran > 0), Holm düzeltmesi aile içinde (aile = sekans × C0 > 0).
- Etki büyüklüğü: medyan oran, 10 000 örnekli bootstrap %95 GA. Pencere başına ayrıca raporlanır.
- Doz-yanıt: log-oranın log10(C0) ile Spearman korelasyonu (C0 > 0 düzeyleri, bütün çiftler), tek yönlü.

**HA2 — F15, gerçekçi tele.**
- Öngörü: C0 ∈ {1e-13, 1e-12, 3e-12, 1e-11} düzeylerinin hepsinde RPE1s temizle eşdeğer (TOST ±%15).
  Medyan oran 0.9–1.1.
- Çürütme: herhangi bir düzeyde Holm sonrası anlamlı bozulma ve medyan oran > 1.15. GA sınırı aşıp anlamlı olmayan
  durum "kararsız" diye raporlanır.

**HA3 — F5, zoom.**
- Öngörü: C0 ∈ {1e-12, 3e-12, 1e-11} düzeylerinin üçünde de RPE1s anlamlı artar (Holm). Beklenen medyan oranlar
  1e-12: 1.05–1.30, 3e-12: 1.20–1.60, 1e-11: 1.60–2.60. 1e-13 temizle eşdeğerdir (TOST ±%15).
  Doz-yanıt tekdüze (Spearman ρ > 0, p < 0.05).
- Çürütme: 1e-11'de anlamlı artış yok ya da ρ ≤ 0.

**HA5 — Metrik ölçek.**
- Öngörü: bütün hücrelerde Sim(3) ölçek hatası < %2 (pilot en çok %0.3).

**Keşif (hipotez değil):** aynı σ_tilt (~0.2 px) pilotta F5 1e-12'de 6/6, F15 1e-11'de 3/6 çiftte oran > 1 verdi.
Piksel cinsinden σ_tilt'in tek başına açıklayıcı olup olmadığı (D/r0, bulanıklık oranı, açısal çözünürlük) iki
kameranın doz-yanıt eğrileri σ_tilt ekseninde karşılaştırılarak keşifsel olarak incelenir.

**Bu dondurmanın kapsamı dışında (ayrı ön-test sonrası, tarihli sapma kaydıyla eklenecek):**
- **HA4 (derinlik bağımlılığı)** ve **H3'ün AirSim testi (kollar A0, AG, A3, A0C+A3).** İkisi de PyEco kol koşucusunun
  (`turb/runner.py`) AirSim verisine bağlanmasını ve nokta başına artık dökümünü gerektiriyor. Bu hat henüz ön-testten
  geçmedi. Öngörüleri v0.2'deki gibidir (HA4: F5 @ 1e-11'de ek yeniden-izdüşüm artığı derinlik dilimiyle artar,
  Spearman ρ > 0.8; H3: temiz A0 ANEES > 10).

**Analiz kodu:** `airsim_sim/analiz_cift.py` (eşli oranlar; yörünge–kayıt ATE eşleşme denetimi) ve dondurmadan sonra
yazılacak `airsim_sim/analiz_v03.py` (yukarıdaki testlerin birebir uygulaması; önce pilot verisinde sınanır).

## 5. Dondurma

- [x] Kullanıcı onayı: 2026-11-08 ("önerilerine katılıyorum, v0.3'e işle ve başlat").
- Dosya sha256'sı `evidence/on_kayit_sha256.txt`'ye yazılır, sonra koşum başlar.
- EuRoC A0C k değeri §1'deki yordamla 3090'da yeniden seçilir; seçilen değer sapma değildir, `ON_TEST_FAZ2_EUROC.md`'ye
  yazılır.
- Sapma kayıtları: aşağıda (Sapma 1, …).

### Sapma 1 — 2026-11-08: AirSim kolları (H3'ün AirSim testi + HA4) eklendi
- §4 "kapsam dışı" maddesinin uygulanışı. Öngörüler v0.2'den aynen; işlemsel tanımlar, ızgara ve kararlar
  `airsim_sim/africa/SAPMA_TASLAK_KOLLAR.md` (§1–4, §7) dosyasında; ön testler `airsim_sim/africa/ON_TEST_KOLLAR.md`.
- Hepsi türbülanslı kol ve artık sonuçlarından önce yazıldı. Pilot görüntüleme beyanı: temiz A0 ANEES_poz 1.09 (F15) /
  1.81 (F5), bandın altında; tek pilot çiftinde HA4 ρ = 0.976 (HA4 tanımı bundan önce yazılmıştı).
- H3: medyan ANEES_poz bant altıysa "öngörü tutmadı" (çürütmeyle aynı sonuç). A0C+A3 AirSim'de kapsam dışı
  (A0 aşırı güvenli olmadığı için k ≥ 1 kalibrasyonunun hedefi yok).
- Izgara: F5 × C0 {0, 1e-12, 1e-11} × tohum 10–14 × {W1, W2} = 30 hücre × {A0, AG, A3} + 20 artık dökümü.

### Sapma 2 — 2026-11-11: hedef dergi (IEEE TAES) için ek rakip koşumları
- **Zaman ve bilinenler beyanı:** bütün ön kayıtlı koşumlar bittikten sonra yazıldı. Bilinen sonuçlar: ana ızgara
  (K = 4), H2 (çürüdü), H6 ham tablosu, rakiplerin D/r₀ {0, 8, 16, 24} sonuçları, AirSim kolları. Aşağıdaki ek
  hücrelerin hiçbir sonucu görülmedi.
- **Gerekçe:** kırılma D/r₀ = 4'te çıktı, ama rakipler 4'te koşulmadı. Hedef dergi okurunun (havacılık, seyrüsefer)
  "ataletsel destek bunu çözmez mi" sorusuna veri yok.
- **(a) Rakipler D/r₀ = 4'te:** ORB-SLAM3 stereo, AirSLAM, DROID-SLAM; 9 sekans × tohum 0–2; §1 rakip protokolü aynen
  (önek koşumu, σ = 1 DN, ORB-SLAM3 tohum başına 3 tekrar). Kayıtlar ayrı dizinde (`evidence/FAZ2_rakip_ek/`).
  - **HR1 — Erken kırılma PyEco'ya özgü değil.** Öngörü: üç rakipten en az ikisinde D/r₀ = 4'te eşli medyan ATE oranı
    (4 / 0; aynı sekans, tohum, pencere; ikisinde de başarılı pencereler) ≥ 2. Çürütme: en fazla bir rakipte ≥ 2 →
    makale erken kırılmayı PyEco'ya özgü diye raporlar.
- **(b) ORB-SLAM3 stereo-ataletsel (SI) rakibi:** aynı bozulmuş görüntüler + EuRoC'un özgün `imu0` verisi
  (değiştirilmeden; türbülans IMU'yu etkilemez). D/r₀ ∈ {0, 4, 8, 16, 24}; 9 sekans × tohum 0–2 × 3 tekrar. Ayar:
  yazarın EuRoC stereo-ataletsel ayarı (ORB değerleri ORB-SLAM3 stereo koşumuyla aynı, nFeatures 1200; IMU gürültü
  değerleri yazar dosyasından), döngü kapama kapalı. Koşucu: `stereo_asl_headless`'in IMU beslemeli eşleniği (yeni
  ikili; kimliği ve sha256'sı kaydedilir). Ön test: temiz EuRoC'ta SI ATE'si yayımlanmış düzeyde (MH için birkaç cm).
  - **HI1 — Ataletsel destek ıraksamayı azaltır, türbülansın doğruluk etkisini gidermez.** Öngörü: D/r₀ = 16'da SI'nin
    başarısız pencere oranı stereo ORB-SLAM3'ünkinden düşük; ama SI'nin başarılı pencerelerde eşli medyan ATE oranı
    (16 / 0) ≥ 2. Çürütme: SI oranı < 2 → makale "ataletsel destek türbülansın doğruluk etkisini büyük ölçüde
    gideriyor" diye raporlar. İkinci yarı (başarısızlık) tutmazsa o da raporlanır.
- **(c) Derin onarma kolu:** model kurulumu sürüyor. Tanım, ızgara ve öngörü model seçildikten sonra, koşumdan önce
  Sapma 3 olarak eklenecek.
- (a) ve (b) betimseldir (Holm ailesi yok). Metrikler rakip analizindekiyle aynı: başarısız pencere oranı, ATE
  medyanı, eşli oranlar (bootstrap %95 GA, 10 000).
- **Sapma 2 eki (2026-11-11, SI ön testinden sonra, türbülanslı SI sonucu görülmeden):** temiz ön testte SI'nin gerçek
  zamanlı (nedensel) yörüngesi IMU ilklemesi (VIBA) süren pencerede bozuk çıktı (MH01 w721: nedensel 0.302 m, son
  yörünge 0.010 m; V2_01 w100: 0.184 / 0.010 m). Nedeni: ilklemede harita yerçekimine döndürülür ve ölçeklenir,
  nedensel pozlar pencere ortasında sıçrar. Bu yüzden HI1 **son yörüngeyle** (ms_ metrikleri, başarısızlık ölçütü
  ms_durum_p) değerlendirilir, karşılaştırılan stereo ORB-SLAM3 için de aynı tür (ms_) kullanılır. Nedensel metrikler
  ikincil olarak raporlanır.

### Sapma 3 — 2026-11-11: derin onarma kolu (öğrenilmiş türbülans giderme ağı ön-işlemci olarak)
- **Zaman ve bilinenler beyanı:** bütün ön kayıtlı koşumlar bittikten sonra, Sapma 2 koşumları sürerken yazıldı.
  Bilinenler: ana ızgara, H1–H6, rakip sonuçlarının bir bölümü, keşif sondasının (`evidence/SONDA_KESIF/`) tek
  pencerelik duman testi (V1_01, D/r₀ 8: hızlı rejimde A0 ATE'si 2.2 kat kötü). Sondanın kendisi sürüyor, kaydı
  görülmedi. Onarım ağlarının yalnız görüntü düzeyi ölçümleri biliniyor: bağımsız sentetik bozulmada (eğim 1.5 px,
  bulanıklık, gürültü) MambaTM çevrimdışı PSNR 23.20 → 23.88 dB, temize göre Farneback RMS kayması 1.93 → 1.66 px;
  nedensel L = 5: 23.43 dB, 1.78 px. Benzetici görüntüsünde yalnız 60 karelik süre denemesi yapıldı (ortalama mutlak
  değişim 2.9 DN, parlaklık +2 DN; temize göre ölçüm yapılmadı). **Onarılmış görüntüde hiçbir VO sonucu görülmedi.**
- **Gerekçe:** hedef dergi okurunun "görüntüyü önce bir türbülans giderme ağıyla onarmak yetmez mi" sorusuna veri yok.
- **Model:** MambaTM (xg416/MambaTM, commit 79aa41f45d2f; yazarların yayımladığı `MambaTM_dynamic.pth`, sha256
  d207a6cd…1a45c; yeniden eğitim ya da ince ayar yok; fp32; 752×480 kare iki 480×480 yamaya bölünür, örtüşme
  ortalanır). Seçim gerekçesi: sentetik testte üç adayın (DATUM, MambaTM, TMT) en iyisi ve en hızlısı. Sıra kancadaki
  onarım aşamasıyla aynı: türbülans → sensör gürültüsü (σ = 1 DN) → onarım → VO ön-ucu. İki kamera bağımsız onarılır.
  - **O (çevrimdışı, nedensel değil):** 20 karelik parçalar, parça kenarında 2 kare atılır; ağ gelecek kareleri de
    görür. Çevrimiçi bir sistem için üst sınırdır.
  - **C (nedensel):** her kare için yalnız [t−4, t] penceresi verilir, son çıktı alınır.
  - **N (onarım yok):** aynı dosya yolundan geçen denetim hücresi.
- **Hat (parmak izli kodda değişiklik yok; code_fp 08e28b83407c):** `turb/uret_sekans.py` (rektifiye, `sensor_yaml`,
  kare 0 → son pencere sonu; görüntüler kancanınkiyle bayt-aynı, rakip ön testi G0b/G0f) →
  `kod/derin_onarim/betikler/onar_sekans.py` → PyEco, kanca kapalı (`turbulence.enabled = False`, sensör gürültüsü 0;
  ikisi de görüntüde zaten var). Arka-uç kolları gerçek D/r₀'ı alır (A1 kâhin değeri ham türbülansınkidir, onarımdan
  sonra değiştirilmez). Koşucu `onarim_kos.py`, kayıtlar `evidence/ONARIM/`.
- **Izgara:** 9 sekans × tohum 0 × D/r₀ = 8 × {R0: τ₀/Δt = 1 (ana ızgara rejimi) → N, O; R5: enine hız 5 m/s
  (τ₀/Δt = 0.0157) → N, O, C} × kollar {A0, A1, A3}; pencereler ana ızgarayla aynı. 45 hücre. C yalnız R5'te koşulur
  (R0'da eğim kareden kareye 0.09 px değişir, zamansal ağın ortalayacağı bir şey yok; GPU süresi 5 kat).
- **Birim ve istatistik:** birim = sekans (n = 9); ATE pencerelerde aritmetik, ANEES geometrik ortalama; eşli medyan
  oran, bootstrap %95 GA (10 000). Betimseldir (Holm ailesi yok, tek tohum). Temiz başvuru: ana ızgara D/r₀ = 0,
  tohum 0, A0.
- **Öngörüler (sonuçtan önce):**
  - **HD1 — Donmuş eğimde onarım VO'yu değiştirmez.** R0: eşli medyan ATE oranı (O / N, A0) [0.8, 1.25] içinde.
    Dışındaysa yönüyle birlikte raporlanır.
  - **HD2 — Hızlı rejimde çevrimdışı onarım doğruluk kaybının küçük bölümünü geri alır.** R5: (a) eşli medyan ATE
    oranı (O / N, A0) ≥ 0.7; (b) onarılmış görüntüde A0, temizden hâlâ en az 2 kat kötü: eşli medyan ATE oranı
    (O, D/r₀ 8 / temiz) ≥ 2. Çürütme: (a) < 0.7 → makale "onarım hızlı rejimde doğruluğun önemli bölümünü geri alıyor"
    diye raporlar; (b) < 2 → "onarım temiz doğruluğa yaklaştırıyor".
  - **HD3 — Onarımdan sonra da kovaryans kolu gerekir.** İki rejimde de onarılmış görüntüde eşli medyan ANEES oranı
    (A1 / A0) ≤ 0.33. Çürütme: > 0.33 → "onarım tutarlılık açığını büyük ölçüde kapatıyor".
  - **HD4 — Nedensel kullanım çevrimdışından iyi değildir.** R5: eşli medyan ATE oranı (C / N, A0) ≥ 0.8.
- **İkincil (betimsel):** ön-uç verimi ve stereo eşleşme sayısı, A3'ün çevrimiçi σ̂ değeri (N'ye göre), görüntü düzeyi
  ölçüler (temize göre PSNR, SSIM, Farneback RMS), kare başına süre, A3 / AG karşılaştırması yok.
- **Ön test (koşumdan önce, `evidence/ONARIM/ON_TEST.md`):** dosya yolu N hücresi kanca yolundaki ana ızgara
  kaydıyla uyumlu (kalibrasyon yaml'ı 12 haneye yuvarlandığı için bit-aynılık beklenmez; pencere ATE farkı medyanı
  < %1; > %10 ise dur); temiz görüntünün (D/r₀ = 0) onarımı VO'yu bozmaz (ATE oranı [0.67, 1.5]); onarım tekrarı
  bayt-aynı; süre ölçümü.
- **Sınır (makalede yazılacak):** ağ bu benzeticinin dağılımında eğitilmedi (yayımlanmış ağırlıklar); sonuç "hazır
  ağ ön-işlemci olarak" durumunu ölçer, eğitilebilecek en iyi ağı değil.

### Sapma 4 — 2026-11-11: A3S ayırma kolu (A3'ün marjinal varyansı, yama korelasyonu olmadan)
- **Zaman ve bilinenler beyanı:** 13:22'de, Sapma 2 ve Sapma 3 koşumları ile keşif sondası sürerken yazıldı. Bilinenler:
  ana ızgaranın bütün sonuçları (A0C+A3, k = 12 dahil), H4 kararı (A3 / AG ANEES medyan oranı 0.83 / 0.79 / 0.85 / 0.87 /
  0.83, D/r₀ 8–24) ve ön kayıt sonrası bulgu F2: D/r₀ = 8'de bu oran yalnız skaler büyüklük farkından öngörülen değere
  eşit, (1 + σ_AG²) / (1 + σ̂_A3²) = 0.83 (12 ve 16'da öngörü 0.82 ve 0.90, ölçülen 0.79 ve 0.85). Yani A3'ün AG'ye
  üstünlüğünün yama korelasyon yapısından mı, yalnız daha büyük σ̂'dan mı geldiği eldeki veriyle ayrılamıyor.
  **k = 1 ile koşulmuş hiçbir A0C+A3 (A3S) sonucu görülmedi.**
- **Gerekçe:** H4'ün "fizik-tabanlı katkı korunur" kararının mekanizması (korelasyon yapısı mı, DIMM büyüklüğü mü)
  doğrudan sınanmadı; hedef dergi hakemi bu ayırmayı ister.
- **Kol (parmak izli kodda değişiklik yok; code_fp 08e28b83407c):** A3S = mevcut "A0C+A3" kolu, `--turb_cal_factor 1.0`
  ile. Ölçüm σ'sı σ² = σ₀² + ŝ², ŝ = A3'ün kullandığı çevrimiçi DIMM kestirimi (`dimm_sigma_scale` = 1.0, aynı
  kestirici); her gözlem ayrı faktördür, yama faktörü kurulmaz. A3'te tek başına kalan gözlemin marjinal σ'sı da budur
  (`_add_single`). A3 ile A3S arasındaki fark yalnız ortak yama faktörüdür (korelasyonlu Σ ve onun dayanıklı çekirdeği).
- **Izgara:** ana ızgaranın aynısı (9 sekans × tohum 0–4, τ₀/Δt = 1, σ = 1 DN, aynı pencereler, tek kart RTX 3090),
  D/r₀ ∈ {8, 12, 16, 20, 24} → 225 hücre, tek kol. Eşleme ana ızgaranın A3, AG, A0, A1 kayıtlarıyla (sekans, tohum,
  pencere). Ön-uç yeniden koşulur (döküm saklanmadı); ön-uç aynı kartta bayt-belirlenimlidir. Denetim: (a) her hücrenin
  `frontend_stats` alanı ana ızgara kaydıyla aynı olmalı; (b) dokuz sekansın (D/r₀ 8, tohum 0) hücresinde A0 da koşulur,
  ana ızgara A0 kaydıyla bit-aynı olmalı. Biri tutmazsa o hücre eşlemeye girmez ve koşum durdurulup incelenir.
  Koşum sırası 8, 16, 24, 12, 20. D/r₀ 0 ve 4 (90 hücre) yalnız betimseldir ve yalnız süre kalırsa koşulur; bu karar
  sonuçlara bakılmadan, süreye göre verilir.
- **Birim ve istatistik:** ana analizle aynı (`turb/faz2_analiz.py`): birim (sekans, tohum), n = 45; ATE pencerelerde
  aritmetik, ANEES geometrik ortalama; medyan oran, bootstrap %95 GA (10 000, tohum 0), iki yönlü Wilcoxon, Holm
  ailesi = beş D/r₀ seviyesi (karşılaştırma ve metrik başına). Eşdeğerlik: §4'teki kural (TOST, sınır ±ln 1.15,
  α = 0.05; birim log-oranlarının ortalamasının %90 t-GA'sı sınır içinde; bootstrap da raporlanır, karar t-TOST ile).
  Örneklem sabit (ana ızgaranın 45 birimi); tohum eklenmez.
- **Öngörüler (sonuçtan önce):**
  - **HS1 — Yama korelasyonu tutarlılığı ölçülebilir biçimde değiştirmez.** Konum ANEES'i oranı (A3S / A3) D/r₀ 8, 12
    ve 16'da eşdeğerdir (TOST ±%15). 20 ve 24'te de eşdeğerlik beklenir ama güven düşüktür (kırılma bölgesi; ana
    ızgarada A3 ile A0C+A3'ün σ̂ medyanları 2.21'e 2.66 ve 3.80'e 4.09 diye ayrışıyor).
  - **HS2 — Doğruluk değişmez.** ATE oranı (A3S / A3) beş seviyede de eşdeğerdir (TOST ±%15).
  - **HS3 — A3S de AG'yi geçer.** Konum ANEES'i oranı (A3S / AG) D/r₀ 8, 12, 16'da < 1 ve Holm sonrası anlamlıdır;
    beklenen medyan 0.78–0.92. 20 ve 24 için yön öngörüsü yok (AG'nin artık RMS'i bu seviyelerde bazı hücrelerde
    patlıyor: medyan 3.1 ve 18.7 px).
- **Yorum kuralı (sonuçtan önce):**
  - HS1 üç seviyede de tutarsa: H4'ün mekanizması "A3'ün AG'ye üstünlüğü yama korelasyon modelinden değil, DIMM'in
    verdiği büyüklükten gelir" diye yazılır; yarı-statik rejim için önerilen kol A3S olur (arka-uç süresi A3'ün
    yaklaşık beşte biri) ve yama modeli "bu rejimde ölçülebilir katkı vermedi" diye raporlanır.
  - A3S / A3 medyanı > 1.15 ve Holm sonrası anlamlı (8, 12, 16'nın en az ikisinde): korelasyon yapısı tutarlılığa
    katkı veriyor; A3 önerilen kol olarak kalır, A3S ucuz seçenek olarak raporlanır.
  - A3S / A3 medyanı < 0.87 ve anlamlı: skaler kol daha tutarlı; öyle raporlanır.
  - Diğer durumlar seviye seviye, GA'larıyla "kararsız" diye raporlanır.
- **İkincil (betimsel):** dönme ANEES'i; σ̂(A3S) / σ̂(A3) (beklenen: D/r₀ ≤ 16'da ±%5 içinde); A3S / A0 ANEES oranı
  (A3, AG, A1 ile aynı tabloda); arka-uç süresi oranı (aynı koşumdan ölçülen A0C+A3 / A3 oranı ana ızgarada 0.19; yük
  farklı olduğu için yeni ölçüm yalnız betimsel).
- **Ön test (koşumdan önce, `evidence/A3S/ON_TEST.md`):** V1_01_easy, tohum 0, D/r₀ ∈ {0, 8}: denetim hücresi
  (A0 ve k = 12 ile A0C+A3) ana ızgarayla bit-aynı; A3S hücresinde `cal_factor` = 1, yama faktörü sayısı 0, σ̂ ana
  ızgaranın A3 değerinin ±%5'i içinde; bütün alanlar dolu; süre ölçümü. Koşucu `a3s_kos.py`, analiz `a3s_analiz.py`
  (ikisi de parmak izi kümesinin dışında), kayıtlar `evidence/A3S/`.

### Sapma 3 eki — 2026-11-11 15:05: onarım kolu ön testinde tutmayan kapı (G3) ve sayısal ikiz denetimi
- **Zaman ve bilinenler beyanı:** PyEco ön testinden sonra, onarım ızgarasından hiçbir sonuç görülmeden yazıldı.
  Bilinenler: ön testin 15 kaydı (V1_01_easy; `evidence/ONARIM/on_test.jsonl`) ve G3 incelemesinin altı sınaması.
- **Olgu:** ön testteki "dosya yolu N ↔ kanca yolu ana ızgara: pencere ATE farkı medyanı < %1; > %10 ise dur"
  kapısı tutmadı: V1_01_easy, D/r₀ 8, tohum 0'da fark medyanı %42.2 (ANEES %42.6); temiz görüntüde %0.04. Dur kuralı
  uygulandı, ızgara başlatılmadı ve neden araştırıldı (`evidence/ONARIM/ON_TEST.md` §2b; beklenenler ve karar kuralı
  sınamalardan önce yazıldı; sonuçlar `evidence/ONARIM/g3_duyarlilik.jsonl`):
  1. İki yol ön-uca aynı görüntüyü veriyor: ön-uç istatistikleri bit-aynı; dökümlerde kare sayısı, zaman damgaları,
     gözlem kimlikleri, piksel konumları, bayraklar, anahtar noktalar ve tanımlayıcılar bit-aynı (T-F; 1300 kare,
     193 466 gözlem).
  2. Fark yalnız ön-ucun arka-uca verdiği kalibrasyonda (`rect_K` göreli 3.0e-8; ön-uç zaten rektifiye çifte yeniden
     `stereoRectify` uyguluyor) ve o kalibrasyonla üçgenlenen X, Y, Z sütunlarında (göreli medyan 4.0e-8, en çok 2.4e-7).
  3. Arka-uç belirlenimci (T-K: değişmemiş dökümle ana ızgara kaydı bit-aynı; T-D1: dosya yolu ön-ucu ve arka-ucu
     yeniden koşulunca ön test kaydı bit-aynı), ama küçük girdi farkını büyütüyor. MH01, D/r₀ 8, A0: fx ve fy'yi 1e-12
     göreli değiştirmek üç pencerenin ATE'sini +%5, −%17, +%22 değiştiriyor (T-E); dosya yolu kalibrasyonu +%79, +%11,
     −%18 (T-D2). V1_01_easy'de dosya yolu kalibrasyonu (gözlemler aynı) −%26, +%471, +%129 (T-D1k). Yanıt pertürbasyonun
     büyüklüğüyle ölçeklenmiyor (1e-12 ile 3e-8 aynı mertebede etki veriyor).
  Sonuç: kapının "< %1" beklentisi arka-ucun küçük pertürbasyona düzgün yanıt verdiği varsayımına dayanıyordu; varsayım
  yanlış. Veri ya da boru hattı hatası bulunmadı. **Bu bir sapmadır:** ön kayıttaki dur kuralı tetiklendiği halde ızgara
  aşağıdaki eklerle başlatılır.
  - (a) Onarım karşılaştırmaları yalnız dosya yolu içinde kalır (N, O, C aynı kalibrasyonu kullanır); ana ızgara
    hücreleriyle sayı sayı karşılaştırma yapılmaz. HD2 (b)'deki temiz başvuru (ana ızgara D/r₀ = 0, tohum 0, A0) kalır:
    temiz görüntüde iki yolun farkı %0.04 (G4).
  - (b) **Sayısal ikiz denetimi (yeni, sonuçtan önce tanımlandı; `ikiz_kos.py`, parmak izi kümesinin dışında):** her N
    hücresinin (18 hücre: 9 sekans × R5, R0) ön-uç dökümü saklanır (spool dizininden sabit bağlantı; ızgara koşucusu
    değişmez); aynı dökümde yalnız `rect_K`'nın fx ve fy değerleri (1 + 1e-12) ile çarpılarak A0 ve A1 yeniden koşulur
    (N′; çıktı `evidence/ONARIM/ikiz.jsonl`). Birim başına eşli oran N′ / N (ATE, ANEES) O / N ile aynı istatistikle
    (eşli medyan, bootstrap %95 GA) raporlanır. HD1–HD4 kararları ön kayıttaki eşiklerle verilir, eşikler değişmez; her
    kararın yanına N′ / N medyanı ve birim oranlarının [en küçük, en büyük] aralığı yazılır. Karşılaştırılan medyan oran
    (O / N ya da C / N) N′ / N birim oranlarının aralığının içindeyse karar "sayısal gürültüden ayırt edilemiyor" notuyla
    raporlanır. A3 ikizde koşulmaz (süre: kol başına ≈ 36 dk); HD3 için A1 / A0 ANEES oranı N′ üzerinde de verilir.
  - (c) Makalede sınır olarak yazılır: türbülans altında tek bir koşumun pencere düzeyi sonucu 1e-12 düzeyindeki girdi
    farkına duyarlıdır; sonuçlar yalnız eşli ve çok birimli istatistikle yorumlanır.
  - (d) **Beklenenler (sonuçtan önce):** (1) ikiz denetimi (eps = 0; V1_01_easy'nin iki N hücresi, A0 ve A1, 12 kayıt)
    ızgaranın N kayıtlarıyla bit-aynı; değilse ikiz koşulmaz. (2) Izgaranın V1_01_easy R0 N kayıtları (A0, A1, A3; üç
    pencere) ön testin aynı hücredeki kayıtlarıyla bit-aynı; değilse ızgara durdurulur. (3) N′ / N'de sistematik yön
    beklenmez: A0 ATE eşli medyan oranı [0.8, 1.25] içinde, birim oranları 1'in iki yanında. Medyan bu aralığın dışına
    çıkarsa HD1 kararı "tolerans bandı sayısal gürültü tabanından dar" notuyla raporlanır. (4) A1 / A0 ANEES oranının
    HD3 eşiğine (≤ 0.33) göre konumu N ve N′ üzerinde aynı; farklıysa HD3 "sayısal gürültüye duyarlı" notuyla raporlanır.
  - (e) Ek maliyet: 18 hücre × (A0 + A1) ≈ 5 işlemci-saat; ekran kartı kullanmaz.

### Sapma 5 — 2026-11-11: hızlı rejim doğrulayıcı ızgarası (sabit enine hız, τ₀ ≪ Δt)
- **Zaman ve bilinenler beyanı:** 18:43'te, keşif sondası (`evidence/SONDA_KESIF/`, ön kayıt dışı) tamamlandıktan
  (18:33) sonra, Sapma 3–4 koşumları sürerken yazıldı. Metnin taslağı 15:42'de, sonda sürerken ve R5 rejiminde hiçbir
  sonuç görülmeden yazılmıştı; sonda bittikten sonra değişen yerler aşağıda "sonda sonrası" diye işaretlidir.
  Bilinenler: ana ızgaranın bütün sonuçları (τ₀/Δt = 1); Sapma 2 rakip koşumlarının sonuçları (yarı-statik rejim, başka
  sistemler); sondanın sonuçları (yalnız tohum 0, 9 sekans, birim = sekans, eşli medyan oran; kollar A0, A1 ve oracle
  ölçekli geniş kapılı A0 = A0g; R5'te D/r₀ 4 / 8; `evidence/SONDA_KESIF/SONUC.md`): A0 ATE R5 / temiz 2.81 / 14.5;
  A0 ATE R5 / R0 1.00 / 1.66; konum ANEES'i R5 / R0: A0 0.63 / 2.24, A1 0.85 / 1.72; R5'te A1 / A0: ANEES 0.19 / 0.07,
  ATE 0.90 / 1.02; ıraksayan pencere (A0, fe0) 0 / 27 ve 0 / 27; A0g / A0 ATE 0.65 / 0.31 (9 / 9 sekans). Sondanın ön
  beyanındaki S1 (kapı) ve S2 (R5 / R0 ATE > 1.3) tuttu, S3 (A1 ANEES'i R5'te R0'ın yarısından küçük) tutmadı (1.72).
  Onarım ızgarasının (Sapma 3) R5 hücreleri (tohum 0, D/r₀ 8; kollar A0, A1, A3) koşuyor; kayıtlarına yalnız bütünlük
  için bakıldı (satır sayısı, durum alanı), metrik değerleri okunmadı, test tablosu üretilmedi. **Görülmeyenler:** R5
  rejiminde tohum 1–4'ün hiçbir sonucu; R5 rejiminde A3, AG, A0C ve A0C+A3 kollarının hiçbir sonucu (hiçbir tohumda).
  A0g bu ızgarada yoktur: uyarlamalı kapı kolu ayrı bir sapma ve kullanıcı onayı gerektirir (sonda ön beyanı S1).
- **Gerekçe:** ön kayıt sonrası bulgu F6: Faz 2'nin bütün hücreleri yarı-statik rejimde koşuldu (τ₀ = Δt; D/r₀ = 8'de
  0.079 m/s enine akışa denk; tilt kareden kareye 0.09 px değişiyor). Gerçek enine hız (rüzgâr + platform) m/s
  düzeyindedir; o rejimde tilt kareden kareye bağımsıza yakındır. H1, H4 ve H5 kararlarının bu rejimde geçerli olup
  olmadığı sınanmadı; hedef dergi hakemi bunu ister. Sonda bunu yalnız tohum 0'da ve yalnız A0, A1 için gösterdi.
- **Izgara (parmak izli kodda değişiklik yok; code_fp 08e28b83407c):** ana ızgaranın hücre tanımı (9 sekans, σ = 1 DN,
  aynı pencereler, cov_every = 2, tek kart RTX 3090), yalnız `tau0_over_dt` farklı: rejim **R5**, sabit enine hız
  v = 5 m/s, τ₀/Δt = 0.314·(D/(D/r₀))/(v·Δt), D = 0.10 m, Δt = 0.05 s → D/r₀ = 8'de 0.0157, 4'te 0.0314. Düzeyler
  D/r₀ ∈ {8, 4} (ana ızgarada karşılığı olan ve sondada koşulan iki düzey); tohumlar 1–4 (tohum 0 sondada görüldüğü
  için doğrulayıcı kümeye girmez); kollar A0, A0C (k = 12), A1, A3, AG, A0C+A3 → 72 hücre, 432 kol koşumu, 1296
  kayıt. Koşum sırası: düzey 8 → 4; düzey içinde tohum → sekans. Örneklem sabit: 36 birim / düzey; tohum eklenmez.
- **Eşleme:** (a) rejim içi: kollar aynı ön-uç dökümünü paylaşır (kol / A0, A3 / AG). (b) rejimler arası: R5 hücresi,
  ana ızgaranın aynı (sekans, tohum, D/r₀, pencere) R0 kaydıyla; ön-uç koşumları farklıdır. (c) bozulma: A0, R5@düzey
  ile ana ızgaranın aynı (sekans, tohum, pencere) temiz (D/r₀ = 0) kaydı.
- **Birim ve istatistik:** ana analizle aynı (`turb/faz2_analiz.py` işlevleri; betik `hizli_analiz.py`, ızgaradan
  önce yazıldı, öz-testi ana ızgaranın aile değerlerini sayı sayı yeniden üretiyor): birim (sekans, tohum), n = 36;
  ATE pencerelerde aritmetik, ANEES geometrik ortalama; eşli medyan oran, bootstrap %95 GA (10 000, tohum 0), iki
  yönlü Wilcoxon işaretli sıra; Holm aileleri metrik başına: **F1** rejim içi kol / A0 = 5 kol × 2 düzey; **F2**
  A3 / AG konum ANEES'i = 2 düzey; **F3** rejim etkisi R5 / R0 = 6 kol × 2 düzey. Geçersiz pencere eşlemeden düşer;
  düşen birim sayısı ve başarısız pencere oranı (rejim, kol, düzey başına) ayrıca raporlanır.
- **Tasarım duyarlılığı (sonuçtan önce; ana ızgaranın birim log-oran saçılımından, hızlı rejimde saçılımın benzer
  olduğu varsayımıyla; normal yaklaşım, Wilcoxon ARE 0.955, Holm'un en sıkı basamağı α / m, güç %80):** n = 36 ile
  saptanabilir en küçük eşli medyan oran: kol / A0 konum ANEES'i 0.85 (D/r₀ 4), 0.76–0.80 (D/r₀ 8); kol / A0 ATE
  0.86–0.88 (A0C: 0.80 / 0.90); A3 / AG ANEES 0.925 (D/r₀ 4), 0.83 (D/r₀ 8); bozulma (A0 / temiz) 1.21 / 1.40.
  Sonuç: (1) HZ4'ün %10 eşiği saptama sınırındadır; "anlamlı iyileşme yok" kararı %10–14 arası gerçek bir
  iyileşmeyi dışlamaz, güven aralığı kararla birlikte verilir. (2) A3 / AG farkı hızlı rejimde yarı-statikteki kadar
  kalırsa (0.68 / 0.83) D/r₀ 4'te güç ≈ 1, D/r₀ 8'de ≈ 0.80'dir; D/r₀ 8'de anlamsız sonuç "fark yok" diye
  yorumlanmaz, güven aralığıyla raporlanır. Örneklem bu hesaba göre büyütülmez (sabit: 36).
- **Öngörüler (sonuçtan önce):**
  - **HZ1 — Bozulma hızlı rejimde de sürer.** A0'ın ATE oranı (R5, D/r₀ = 8) / (temiz) > 1.5 ve kırılma ölçütü
    (medyan ATE > 2 × temiz ya da ıraksama > %20; ana ızgaradaki tanım) D/r₀ ≤ 8'de sağlanır. *Sonda sonrası sayısal
    beklenti:* tohum 0'da oran 14.5 (D/r₀ 8; sekans aralığı 6.8–44) ve 2.81 (D/r₀ 4; 1.7–5.2); tohum 1–4'te medyan
    oran D/r₀ 8'de 8–25, D/r₀ 4'te 2.0–4.0; kırılma düzeyi K = 4 (ana ızgarayla aynı; D/r₀ 4'te eşiğe yakın).
  - **HZ2 — Rejim etkisi (doğruluk).** A0'ın ATE oranı R5 / R0 (eşleme b). Taslakta yön yazılmamıştı (iki karşıt ön
    bilgi: eski hat τ₀/Δt = 0.1'de 1.6 kat kötü; "beyaz tilt düzleştiricide ortalanır" → daha iyi). *Sonda sonrası:*
    tohum 0'da D/r₀ 8'de 1.66 (9 sekansın 9'unda > 1), D/r₀ 4'te 1.00 (4 / 9 sekansta < 1). Öngörü: D/r₀ 8'de
    oran > 1, Holm sonrası (F3) anlamlı, medyan 1.3–2.2; D/r₀ 4'te yön öngörülmez, medyan 0.85–1.20, anlamsız.
  - **HZ3 — Tutarlılık kolları.** (i) A1, A3 ve AG'nin konum ANEES'i iki düzeyde de A0'a göre Holm sonrası anlamlı
    düşer. (ii) **A3 / AG:** yarı-statik rejimde 0.68 (D/r₀ = 4) ve 0.83 (D/r₀ = 8) idi. Beklenti: hızlı rejimde A3 üstünlüğünü yitirir,
    oran ≥ 1 olur. Gerekçe: AG artık RMS'ten kestirir; yarı-statik rejimde yavaş tilt duruma soğrulduğu için artıklar
    küçük kalır (ana ızgarada σ_AG medyanı D/r₀ 4 / 8'de 1.01 / 1.32 px; oracle tilt 1.23 / 2.19 px), hızlı rejimde tilt
    artıklarda tam görünür ve σ_AG büyür; A3'ün DIMM kestirimi (1.16 / 1.51 px) anlık sol-sağ farkından gelir, rejimden
    bağımsız kalması beklenir. *Sonda sonrası:* (i) için A1 / A0 konum ANEES'i tohum 0'da 0.19 (D/r₀ 4) ve 0.07
    (D/r₀ 8); beklenen aralık 0.12–0.28 ve 0.04–0.12. A3 ve AG sondada koşulmadı; (ii)'deki beklenti sondadan
    bağımsızdır ve değiştirilmedi (öznel olasılık, A3 / AG ≥ 1: D/r₀ 8'de ≈ %65, D/r₀ 4'te ≈ %55).
  - **HZ4 — Kovaryans kolları doğruluğu iyileştirmez (H5'in karşılığı).** F1'de hiçbir kol ATE'yi A0'a göre Holm
    sonrası anlamlı ve %10'dan fazla düşürmez (oran < 0.90). *Sonda sonrası:* A1 / A0 ATE tohum 0'da 0.90 (D/r₀ 4;
    7 / 9 sekansta < 1) ve 1.02 (D/r₀ 8; 4 / 9). Beklenen: bütün kol / A0 ATE medyanları 0.88–1.12 içinde; eşiğe en
    yakın aday A1, D/r₀ 4 (saptama sınırı 0.86–0.88; gerçek etki 0.90 ise büyük olasılıkla anlamsız çıkar).
  - **HZ5 — Rejim etkisi (tutarlılık).** Taslaktaki (15:42, sonda görülmeden) öngörü: "A1'in konum ANEES'i R5'te
    R0'dakinden düşüktür (gürültü beyaza yaklaşınca skaler şişirme daha tutarlı): oran R5 / R0 < 1, Holm sonrası
    anlamlı; A3 ve AG için yön < 1." **Sonda tohum 0 bunu D/r₀ 8'de çürüttü:** A1'de oran 1.72 (9 sekansın
    8'inde > 1), A0'da 2.24 (9 / 9); D/r₀ 4'te yön taslaktaki gibi: A1 0.85 (7 / 9 sekansta < 1), A0 0.63 (7 / 9). Öngörü bu
    yüzden *sonda yönünün bağımsız tohumlarda yinelenmesi* olarak, düzeye göre yeniden yazıldı: D/r₀ 8'de A1
    oranı > 1 (medyan 1.3–2.3) ve A0 oranı > 1 (1.6–3.0), ikisi de Holm sonrası (F3) anlamlı; D/r₀ 4'te A1 oranı < 1
    (0.65–1.05) ve A0 oranı < 1 (0.5–0.9), anlamlılık öngörülmez. A3 ve AG için aynı oran betimsel raporlanır
    (yön öngörüsü: D/r₀ 8'de A3 > 1, çünkü σ̂ rejimden bağımsız kalırken hata büyür; AG için öngörü yok).
- **Yorum kuralı (sonuçtan önce):**
  - HZ3 (ii)'de A3 / AG < 1 ve Holm sonrası anlamlı (en az bir düzeyde): "fizik-tabanlı kestirimin artık-tabanlı
    kestirime üstünlüğü iki rejimde de geçerli" diye yazılır.
  - A3 / AG ≥ 1 (anlamlı olsun olmasın, iki düzeyde de): ana ızgaranın H4 kararı "yarı-statik rejime özgü" diye
    sınırlandırılır; makalenin katkı cümlesi rejime bağlı öneriye çevrilir (yavaş tilt: DIMM tabanlı kol; hızlı
    tilt: artık tabanlı kol yeterli) ve bu sınır özet ile sonuç bölümünde açıkça yazılır.
  - HZ1 tutmazsa (oran ≤ 1.5): kırılma eşiğinin rejime bağlı olduğu yazılır; H1 kararı "yarı-statik rejimde" diye
    sınırlandırılır.
  - HZ2 ve HZ5 D/r₀ 8'de tutarsa: makalede "yarı-statik ana ızgara D/r₀ 8'de iyimserdir; gerçekçi enine hızda hata ve
    tutarsızlık daha büyüktür" diye yazılır ve ana ızgaranın mutlak düzeyleri "yarı-statik rejimde" diye
    sınırlandırılır. D/r₀ 8'de oran ≤ 1.15 ya da Holm sonrası anlamsız çıkarsa "sondadaki rejim etkisi bağımsız
    tohumlarda yinelenmedi" diye yazılır; sonda bulgusu keşif notu olarak kalır.
  - HZ4 çürürse (bir kol ATE'yi anlamlı ve > %10 iyileştirirse): H5'in "kovaryans modeli doğruluğu değiştirmez"
    kararı rejime bağlıdır; o kol ve düzey açıkça raporlanır.
  - Her karar ana ızgaradaki karşılığıyla aynı tabloda, iki rejim yan yana verilir. Sapma 3 ekinin sınırı geçerlidir:
    tek hücre düzeyinde sayısal duyarlılık büyük olduğu için yalnız eşli, çok birimli istatistik yorumlanır.
- **İkincil (betimsel):** dönme ANEES'i; σ̂ (A3) ve σ_AG medyanlarının R5 / R0 oranı (beklenen: σ̂ ±%20 içinde,
  σ_AG > 1); ANEES'in %95 ki-kare bandına göre konumu (kol, düzey, rejim); ön-uç verimi (zamansal eşleşme sayısı);
  arka-uç süresi.
- **Ön test (koşumdan önce, `evidence/HIZLI/ON_TEST.md`; beklenenler sonuçtan önce yazılır):** V1_01_easy: (1) tohum
  0, R5, D/r₀ 8, kollar A0 ve A1 → sonda kaydıyla bit-aynı ve ön-uç istatistikleri aynı (gerileme kapısı); (2) tohum
  1, R5, D/r₀ 8 ve 4, altı kol (duman, alan doluluğu, süre). Izgarada aynı hücreler yeniden koşulur ve ön test
  kaydıyla bit-aynı olmalıdır. Koşucu `hizli_kos.py`, analiz `hizli_analiz.py` (ikisi de parmak izi kümesinin
  dışında), kayıtlar `evidence/HIZLI/`.
- **Maliyet:** hücre başına ön-uç ≈ 4–10 dk (yüke göre), arka-uç ≈ 75 dk (altı kol; A3 ≈ 36 dk) → ≈ 100 işlemci-saat;
  A3S ve onarım ızgarasıyla birlikte gece boyunca koşar.
- **Döküm saklama (kararlara etkisi yok):** hücrelerin ön-uç dökümleri spool dizininden sabit bağlantıyla saklanır
  (koşucu değişmez). Amaç: kullanıcı uyarlamalı kapı kolunu onaylarsa o kol, ayrı bir sapma olarak, bu ızgaranın A0
  kayıtlarıyla aynı ön-uç çıktısında koşabilsin.

### Sapma 6 — 2026-11-11: uyarlanır izleme kapısı kolu (çevrimiçi yöntem kolu) ve doğrulayıcı ızgarası
- **Zaman ve bilinenler beyanı:** saat 21:20; kapı kolunun geliştirme koşumu (tohum 0) bitip (20:42) yapılandırma
  dondurulduktan (20:46) ve doğrulayıcı koşucunun ön testi geçtikten (21:14) sonra, Sapma 3–5 koşumları sürerken yazıldı. Tasarım bölümlerinin (ızgara, eşleme, istatistik,
  hipotez aileleri) taslağı 20:03'te, geliştirme koşumu sürerken ve oran tablosu görülmeden yazıldı; geliştirme sonucu
  görüldükten sonra eklenen yerler "geliştirme sonrası" diye işaretlidir. **Bilinenler:** ana ızgaranın bütün
  sonuçları; Sapma 2 rakip koşumları; ilk keşif sondası (`evidence/SONDA_KESIF/SONUC.md`; tohum 0, birim = sekans;
  oracle ölçekli geniş kapılı A0 = A0g): A0g / A0 ATE R0-8'de 0.517 (9 / 9 sekans), R5-4'te 0.65, R5-8'de 0.31 (9 / 9);
  A0g / A0 konum ANEES'i R0-8'de 0.044; keşif sondası 2 (kapı ayırma; tohum 0, R0-8): sürüyor, 29 koşumun 8'i bitti (K kontrolünün bit-aynılık kapısı 6 / 6 geçti); oran tablosu üretilmedi, görülmedi (betik veri tam değilken tablo basmaz); kapı kolunun ön
  testi (`pyeco_turb_kapi/evidence/KAPI/on_test_1930.log`; tohum 0, V1_01_easy, R0-8; ATE oranı kol / A0): D1 0.896,
  D2 0.749, A1g 0.798, yalnız gürültü grubu (D1n) 0.951, yalnız PnP (D1p) 1.094; geliştirme koşumu (tohum 0, R0-8,
  9 sekans): D1 / A0 ATE 0.512 [0.42, 0.90] (9 / 9 sekans < 1), D2 / A0 0.522 [0.35, 1.19] (8 / 9), D1 / A3S 0.437, D1 / A0g 1.082, D1 / A1 0.519; konum ANEES'i D1 / A0 0.174, D2 / A0 0.062; başarısız pencere 0 / 27; D1 kapı sigması medyanı 1.98 px (`SONUC_GELISTIRME.md`; başvuru olarak ana ızgaranın ve A3S ızgarasının **tohum 0**, D/r₀ 8 kayıtları okundu: A1 / A0 1.214, A3S / A0 1.214); doğrulayıcı koşucunun ön testi (tohum 0, V1_01_easy; R0-8, R0-0, R5-8, R0-16):
  13 / 13 koşum, sert kapılar geçti; ATE oranı kol / A0: R0-0 A3S 0.997, KAPI 1.012; R5-8 KAPI 0.437, A3S 0.779, A1g 0.542; R0-16 KAPI 0.660, A1g 0.530 (R0-8 ilk ön testle bit-aynı: 0.896 / 0.798); hızlı rejim ön testi (Sapma 5; **tohum 1**, V1_01_easy, R5, D/r₀ 8 ve 4; eski kodun altı kolu): bu
  hücrelerin A0 ve A1 kayıtları HK2 / HK3'ün başvurusudur ve 36 birimden birinin başvuru değeri görülmüştür (kapı
  kolu sonucu görülmedi). Hızlı rejim (Sapma 5), A3S (Sapma 4) ve onarım (Sapma 3) ızgaraları koşuyor; kayıtlarına
  yalnız bütünlük için bakıldı (satır sayısı, durum alanı, bit-aynılık kapıları), metrik değerleri okunmadı, test
  tablosu üretilmedi (istisnalar yukarıda sayıldı: A3S ızgarasının tohum 0, D/r₀ 8 kayıtları ve hızlı rejim ön testi). **Görülmeyenler:** tohum 1–4'te kapılı hiçbir kolun (KAPI, A1g, A0g ve türevleri) hiçbir
  sonucu (hiçbir rejimde); tohum 1–4'te A3S'in R5 ve temiz sonuçları; KAPI'nin R5-4 sonucu (hiçbir tohumda).
- **Gerekçe:** ön kayıtlı H5 kararı: kovaryans kolları (oracle A1 dahil) tutarlılığı düzeltiyor ama doğruluğu
  iyileştirmiyor. Keşif sondaları nedeni gösterdi: izleme cephesinin kapıları (motion-only BA sigması, χ² aykırı
  eleme, graf kapısı, yerel harita süzgeci, PnP-RANSAC eşiği) σ₀ = 1 px'e göre sabittir; türbülans tilt'i geçerli
  eşleşmeleri bu kapıların dışına iter, kovaryans şişirmesi graf faktörlerine ulaşmadan önce veri atılmış olur.
  Kapılar oracle sigmayla genişletilince ATE yarıya iniyor (tohum 0). Kullanıcı 2026-11-11'de çevrimiçi kapı kolunu
  ve doğrulayıcı ızgarasını onayladı. Makalede yöntem iddiası olarak yazılabilmesi için (i) oracle bilgi kullanmayan
  bir kol, (ii) geliştirmede görülmemiş tohumlarda, (iii) sonuçtan önce yazılmış hipotezlerle sınanmalıdır.
- **Kod değişikliği (ilk kez parmak izli kod değişiyor):** yeni bayrak `--turb_gate` (varsayılan 0 = kapalı; bit 1
  gürültü grubu, bit 2 PnP-RANSAC eşiği). Kapılar, kolun **son graf karesinde** kullandığı bireysel-faktör
  sigmasıyla ölçeklenir (nedensel; kestirici ek çağrı almaz; sabit piksel pencereleri ölçeklenmez). Değişen dosyalar
  `config.py`, `gtsam_backend.py`, `turb/backend.py`; geliştirme kopyası `pyeco_turb_kapi`, kod parmak izi
  `bedca2849bb2`. Koşan ağaç `pyeco_turb` (`08e28b83407c`) değişmedi; Sapma 3–5 ızgaraları eski kodla sürüyor. Yeni
  kod, bütün koşumlar bittikten sonra ana ağaca alınır. Kanıtlar (ON_BEYAN_3.md; hepsi geçti): **T1** kapı kapalıyken
  A0, A0C, A1, A3S, A3 eski kodun kayıtlarıyla bit-aynı (15 / 15); **T2** A0C(k = s) + kapı ≡ A0 + üç bayrak ve
  A0C(k = 1) + kapı ≡ A0C(k = 1) (metrik ve iki yörünge dosyasının sha256'sı aynı); **T3** duman (kapı sigması
  alanları, etki uygulandı).
- **Dondurulmuş yapılandırma (geliştirme sonrası):** KAPI = `A0C+A3 --turb_cal_factor 1.0 --turb_dimm_sigma_scale
  1.0 --turb_gate 3` (aday D1): σ² = σ₀² + (c·σ̂_DIMM)², skaler; aynı sigma kapıları ölçekler. Seçim, sonuçtan
  önce yazılmış kuralla (ON_BEYAN_3.md §3; sha256 `598ebfac…`) geliştirme koşumundan yapıldı: D1 / A0 ATE 0.512 (9 / 9 sekans < 1; başarısız pencere 0 / 27), D2 / A0 0.522 (8 / 9; 0 / 27): ikisi de "geçer"; iki medyan oranın oranı 1.019 ∈ [0.90, 1.10] → kural 2 gereği D1 (c = 1, ek sabit yok). Beklenti dökümü: G1, G4–G7 tuttu; G2 tutmadı (D1, beklenen 0.60–0.95 bandının altında: beklenenden iyi); G3'te yön tuttu, D2 / A3S 0.430 beklenen 0.5–0.9 bandının altında. Tutarlılıkta D2 daha iyiydi (ANEES D2 / D1 0.333); tutarlılık ön beyan gereği seçim ölçütü değildir.
  `evidence/KAPI/SECIM.json` sha256 `cd4f4520183e0ee4b4eb72458c116181dd0fefe2da21f21e7440d25796bcae2d`. Bu yapılandırma değiştirilmez; tohum 1–4'te başka kapı varyantı
  koşulmaz (oracle üst sınır A1g dışında). Hiçbir eşik, ölçek ya da kural doğrulayıcı veriye bakılarak ayarlanmaz.
- **Izgara:** 9 sekans, tohum 1–4 (tohum 0 keşifte görüldü, girmez), ana ızgaranın hücre tanımı (σ = 1 DN, aynı
  pencereler, cov_every = 2, RTX 3090). Örneklem sabit: koşul başına 36 birim; tohum eklenmez. Kollar: **KAPI**;
  **A1g** = A1 + `--turb_gate 3` (oracle kovaryans + oracle ölçekli kapı; üst sınır); **A3S** = `A0C+A3
  --turb_cal_factor 1.0` (aynı çevrimiçi kovaryans, kapı kapalı); **A0y** = yeni kodla A0 (yalnız doğrulama).
  - **Parça I — hızlı rejim R5 (v = 5 m/s), D/r₀ ∈ {8, 4}:** Sapma 5 ızgarasının sakladığı ön-uç dökümleri üstünde
    (aynı ön-uç çıktısı): KAPI, A3S, A1g; A0y yalnız tohum 1'de (döküm eşlemesinin denetimi). 234 arka-uç koşumu.
  - **Parça II — yarı-statik R0, D/r₀ ∈ {8, 16}:** A0y, KAPI, A1g. Döküm: A3S ızgarasının sakladığı döküm varsa o
    (D/r₀ 16), yoksa ön-uç aynı ayarla yeniden üretilir (D/r₀ 8). 216 koşum.
  - **Parça III — temiz, D/r₀ = 0:** ön-uç yeniden üretilir; A0y, A3S, KAPI. 108 koşum.
  Toplam 558 arka-uç koşumu, 1674 kayıt. Koşucu `kapi_dogrula.py`, kayıtlar `pyeco_turb_kapi/evidence/KAPI/
  dogrula_{I,II,III}.jsonl`.
- **Başvurular ve eşleme:** birim (sekans, tohum); bütün karşılaştırmalar aynı (sekans, tohum, D/r₀, rejim, pencere)
  içinde ve **aynı ön-uç çıktısı** üstünde eşlidir. A0 ve A1: eski kodun kayıtları (R0: ana ızgara; R5: Sapma 5
  ızgarası). A3S: R0-8 ve R0-16'da Sapma 4 ızgarası (eski kod), R5 ve temizde bu ızgara (yeni kod, kapı kapalı; T1
  gereği eski kodla bit-aynı). Aynı ön-uç çıktısının kanıtı: her A0y kaydı başvuru A0 kaydıyla dört alanda bit-aynı
  olmalıdır (parça II ve III: 108 hücrenin hepsi; parça I: tohum 1'in 18 hücresi — kalan hücrelerde döküm, başvuru
  A0'ın koştuğu dosyanın sabit bağlantısıdır); yeniden üretilen ön-ucun istatistikleri başvuru kaydınınkiyle aynı
  olmalıdır. Biri tutmazsa koşucu durur (çıkış 2) ve o parça yorumlanmaz.
- **İstatistik:** ana analizle aynı (`turb/faz2_analiz.py`; betik `kapi_dogrula_analiz.py`, ızgaradan önce yazıldı,
  öz-testi ana ızgaranın A1 / A0 aile değerlerini sayı sayı yeniden üretiyor): ATE pencerelerde aritmetik, ANEES
  geometrik ortalama; eşli medyan oran, birim bootstrap %95 GA (10 000, tohum 0); iki yönlü Wilcoxon işaretli sıra
  (log oran); Holm aile içinde; α = 0.05. Geçersiz pencere eşlemeden düşer; düşen birim ve başarısız pencere
  sayıları raporlanır. Karar tabloları yalnız veri tam ve doğrulamalar (V1–V5) tutuyorken üretilir.
- **Birincil hipotezler (aile B, m = 3, Holm):** çevrimiçi kapı kolu türbülans altında doğruluğu iyileştirir —
  KAPI / A0 ATE eşli medyan oranı < 1 ve Holm sonrası anlamlı:
  **HK1** R0, D/r₀ 8; **HK2** R5, D/r₀ 8; **HK3** R5, D/r₀ 4.
  Karar: p_Holm < 0.05 ve oran < 1 → "tuttu"; p_Holm < 0.05 ve oran > 1 → "ters yön"; aksi → "tutmadı".
- **HK4 (temiz veride bedel yok; eşdeğer-altı olmama):** D/r₀ = 0'da KAPI / A0 ATE eşli medyan oranının %95 GA üst
  sınırı < 1.10. (Temizde σ̂_DIMM sıfır değildir — ana ızgarada medyan 0.74 px, 0.57–1.16 — kapı ≈ %15–55 genişler;
  bu yüzden gerçek bir sınamadır.)
- **İkincil aileler (Holm aile içinde; kararları birincil iddiayı değiştirmez):**
  **I1** KAPI / A0 ATE, R0-16. **I2 kapının katkısı:** KAPI / A3S (çevrimiçi) ve A1g / A1 (oracle), ATE, dört
  türbülans koşulu (m = 8). **I3 çevrimiçi ↔ oracle:** KAPI / A1g, ATE (m = 4). **I4** A1g / A0, A1 / A0, A3S / A0, ATE
  (m = 12). **I5 tutarlılık:** konum ANEES'i (ve ayrı aile olarak dönme ANEES'i): KAPI / A0, A1g / A0, A3S / A0,
  A1 / A0, KAPI / A3S (m = 20). **T2 temiz, betimsel:** A3S / A0 ve KAPI / A3S (ATE); KAPI / A0 ve A3S / A0 (ANEES).
  Betimsel: kol × koşul birim medyanları (2 × 2 tablo: {oracle, çevrimiçi} × {yalnız kovaryans, kovaryans + kapı} =
  A1, A1g, A3S, KAPI), başarısız pencere sayıları, kapı sigması ve σ̂ dağılımı, ANEES'in χ² bandına göre konumu.
- **Sayısal beklentiler (sonuçtan önce; geliştirme sonrası):** (20:49'da, doğrulayıcı koşucunun ön testi sürerken, hiçbir çıktısı okunmadan yazıldı.) Birim = (sekans, tohum), koşul başına 36; ölçü = eşli medyan ATE oranı.
  - KAPI / A0: **HK1** (R0-8) 0.45–0.75 [geliştirme 0.512; seçim aynı veride yapıldı ve tohum 0'da A1 / A0 1.214 iken tohum 1–4'te 1.027 → iyimserlik payı bırakıldı]; **HK2** (R5-8) 0.35–0.70 [ilk sonda A0g / A0 0.31; σ̂ R0-8 ile aynı düzeyde]; **HK3** (R5-4) 0.60–0.92 [ilk sonda A0g / A0 0.65; kapı ≈ 1.5 px, oracle 2.67 px — en zayıf koşul, "tutmadı" olasılığı gerçek]; **HK4** (temiz) medyan 0.95–1.05, üst sınır < 1.10 [riskli: aşağıdaki duyarlılık]; **I1** (R0-16) 0.50–0.95.
  - **I2:** KAPI / A3S dört koşulda < 1 (R0-8'de 0.40–0.75); A1g / A1 < 1 (0.40–0.80). **I3:** KAPI / A1g R0-8'de 0.95–1.25; R5-8 ve R0-16'da 1.0–2.0 (kestirim oracle'ın altında). **I4:** A1g / A0 R0-8 0.45–0.75, R5-8 0.30–0.55, R5-4 0.55–0.85, R0-16 0.40–0.75; A1 / A0 ve A3S / A0 0.90–1.25 (H5 ile tutarlı; A1 / A0 R0'da ana ızgaradan biliniyor: 1.027 ve 1.039).
  - **I5:** konum ANEES'i KAPI / A0 R0-8'de 0.10–0.40, A1g / A0 0.02–0.10; KAPI'nin ANEES'i χ² beklenen değerinin (3) çok üstünde kalır (geliştirmede medyan 160): tutarlılık iyileşir, sağlanmaz.
  - Kapı sigması medyanı (KAPI): R0-8 1.8–2.2, R5-8 1.6–2.4, R5-4 1.3–1.9, R0-16 2.0–2.8, temiz 1.1–1.6 px. Başarısız pencere: KAPI'de koşul başına ≤ 3 / 108.
- **Tasarım duyarlılığı (sonuçtan önce):** Sapma 5'teki formül: saptanabilir en küçük etki (log oran) = (z_{1−α/(2m)} + z_{0.80})·SS / √(0.955·n); m = 3, n = 36 → 0.552·SS. SS için iki dayanak: geliştirme koşumunda D1 / A0 sekans log oranlarının SS'i 0.320 (tohum 0, n = 9; sekanslar arası etki farkı) ve ana ızgarada kovaryans kollarının A0'a göre birim log oranı SS'i 0.15–0.20 (tohum 1–4, n = 36; sayısal duyarlılık tabanı). SS 0.32 → oran ≤ 0.84; SS 0.45 → ≤ 0.78; SS 0.60 → ≤ 0.72 (%80 güç). Geliştirmedeki etki (0.512; log −0.67) bu sınırların çok ötesinde: etki görülmemiş tohumlarda 0.80'e kadar küçülse de tasarım yakalar; 0.85–0.95 arası gerçek etkilerde güç düşüktür ve "tutmadı" sonucu "etki yok" diye okunmaz (güven aralığı verilir). **HK4:** medyanın %95 GA yarı genişliği ≈ 0.41·SS (log); temizde ana ızgara SS'i 0.15–0.20 → üst sınır ≈ medyan × 1.06–1.08; gerçek oran 1.00 ise geçme olasılığı kabaca %75–80, 1.03 ise ≈ %50. Bu yüzden HK4'ün kalması tek başına "bedel var" kanıtı sayılmaz; medyan ve GA ile raporlanır. Marj 1.10 sonuçtan önce (20:00) analiz betiğinde sabitlendi; değiştirilmez.
- **Yorum kuralları (sonuçtan önce):**
  - HK1–HK3'ün üçü de tutarsa: "çevrimiçi uyarlanır kapı, iki zamansal rejimde ve iki türbülans düzeyinde ATE'yi
    düşürür" diye yazılır; etki büyüklüğü olarak yalnız bu ızgaranın eşli medyan oranları ve güven aralıkları
    verilir (sonda ve geliştirme sayıları değil).
  - Yalnız bir kısmı tutarsa: iddia tutan koşullarla sınırlanır; tutmayan koşul aynı tabloda, güven aralığıyla
    yazılır. Hiçbiri tutmazsa: kol olumsuz sonuç olarak raporlanır; sonda bulgusu "tek tohumda görülen, bağımsız
    tohumlarda yinelenmeyen keşif" diye kalır; özet ve katkı listesinde yöntem iddiası yer almaz.
  - HK4 kalırsa (üst sınır ≥ 1.10): "kol temiz veride doğruluk bedeli getirir" diye sınır olarak yazılır; çözüm
    önerisi (ör. σ̂ eşiğiyle devreye alma) yalnız gelecek iş olarak anılır, bu veriyle yeni kol denenmez.
  - I2: KAPI / A3S < 1 ve anlamlı → "kazancı kovaryans değil kapı taşıyor" (H5 ile tutarlı); A1g / A1 için aynı
    okuma oracle tarafında. I3: KAPI / A1g > 1 ve anlamlı → "çevrimiçi kestirim oracle'ın gerisinde; kestirimin
    iyileştirilmesi gelecek iş"; anlamsız → güven aralığıyla "ayırt edilemedi" (eşdeğerlik iddiası kurulmaz).
  - I5: KAPI'nin ANEES'i A0'a göre düşer ama χ² bandının dışında kalırsa "tutarlılık iyileşir, tam sağlanmaz" diye
    yazılır; KAPI / A3S ANEES oranı > 1 çıkarsa (geniş kapı daha çok aykırı gözlem alır) açıkça raporlanır.
  - R0-16 (I1) birincil değildir: σ̂ orada oracle toplam sigmanın ≈ dörtte biridir (1.86 / 7.86 px); sonuç
    "kestirimin doyduğu düzey" başlığı altında, yönü ne olursa olsun raporlanır.
  - Sapma 3 ekinin sınırı geçerlidir: tek hücre düzeyinde sayısal duyarlılık büyüktür; yalnız eşli, çok birimli
    istatistik yorumlanır. A1g oracle σ_tilt kullanır; gerçek sistemde uygulanamaz, yalnız üst sınırdır.
- **Ara bakış kuralı:** koşum sürerken yalnız `kapi_dogrula_analiz.py --butunluk` (sayım, durum, parmak izi, V1–V5)
  ve koşucu günlüğü (durum, A0y bit-aynılığı, kapı sigması istatistiği; ATE / ANEES yazmaz) okunur. ATE / ANEES
  oranları yalnız 1674 kayıt + başvurular tam olduğunda, bir kez üretilir. Erken durdurma ya da uzatma yoktur.
- **Ön test (`pyeco_turb_kapi/evidence/KAPI/ON_TEST_DOGRULA.md`; beklenenler 20:00'de, sonuçtan önce yazıldı):**
  tohum 0, V1_01_easy, dört hücre (R0-8, R0-0, R5-8 ön-uç yeniden üretimiyle; R0-16 saklı dökümle), 13 arka-uç
  koşumu: K1 A0y ≡ eski A0 (12 / 12), K2 ön-uç istatistikleri, K3 iki döküm yolunun eşdeğerliği (geliştirme
  kayıtlarıyla bit-aynı + yörünge sha256), K4 alanlar ve kapı istatistiği, K4b analiz betiğinin okuması, yumuşak
  beklentiler ve süre. **Sonuç (21:16; §4):** K1 12 / 12, K2 3 / 3, K3 KAPI ve A1g 3 / 3 + yörünge sha256, K4 39 / 39, K4b tam = True (V1 12, V2 9, V3 21, V4 18, V5 39): **hepsi geçti.** Yumuşak beklentilerden ikisi tutmadı: R5-8'de A3S / A0 0.779 (beklenen 0.80–1.20) ve süre (R0-16 koşumları 566–919 s, beklenen 250–700 s; temiz ön-uç 64 s, beklenen 150–400 s).
- **Maliyet:** ön testte ölçülen arka-uç süreleri (V1_01_easy; makine yükü ≈ 65–72, Sapma 3–5 ızgaralarıyla birlikte): R0-8 331–445 s, temiz 325–334 s, R5-8 440–676 s, R0-16 566–919 s; dokuz sekansın ortalama uzunluğu × 1.10 → parça I ≈ 40, parça II ≈ 38, parça III ≈ 11 işçi-saat (toplam ≈ 89). Ön-uç yeniden üretimi 72 hücre (≈ 170 s ve ≈ 70 s; yalnız RTX 3090). Parça I, Sapma 5 ızgarasının döküm üretme hızına bağlıdır (o ızgara bu yükte ≈ 19 saat). Üç parça ayrı süreçlerle, bellek kabul denetimiyle ve kaldığı yerden sürdürülebilir biçimde koşar; parça II önce D/r₀ 8'i, sonra 16'yı koşar. Disk: saklı dökümler ≈ 20 GB (boş alan 78 GB).
- **Makalede kullanım:** doğrulayıcı kanıt yalnız bu ızgaradan gelir. Sonda 1–2, kapı kolunun ön testi ve
  geliştirme koşumu "keşif (ön kayıt dışı), tohum 0" etiketiyle, yöntemin nasıl bulunduğunu anlatan bölümde anılır.

### Sapma 7 — 2026-11-12: kapılı genel uyarlamalı taban çizgisi (AGg)
- **Zaman ve bilinenler beyanı:** saat 22:42; Sapma 6 ızgarası bitip (18:00) çözümlendikten (19:02) ve TAES taslağı
  yazıldıktan sonra, AGg kolunun hiçbir koşumu (ön test dahil) yapılmadan yazıldı. **Bilinenler:** önceki bütün
  ızgaraların sonuçları; özellikle Sapma 6: KAPI / A0 ATE 0.549 (R0-8), 0.362 (R5-8), 0.613 (R5-4), 0.794 (R0-16),
  temiz 1.020 [0.976, 1.047]; A1g / A1 0.406 / 0.299 / 0.671 / 0.554; KAPI / A1g 1.244 / 1.236 / 1.046 / 1.423; KAPI
  kapı sigması medyanı 2.01 / 2.23 / 1.70 / 2.31 px (oracle 4.52 / 4.52 / 2.67 / 7.93 px); kapısız AG kolunun (ana ızgara ve
  hızlı ızgara) ATE'si A0'dan ayırt edilemiyor (≈ 0.96–1.06). Bu sapma için bugün okunan tek yeni nicelik: kapısız
  AG kayıtlarındaki `ag_rms` alanının medyanı (tohum 1–4, ilk pencere): temiz 0.639, R0-8 1.318, R0-16 1.598,
  R5-8 1.517, R5-4 1.096 px. **Görülmeyenler:** AGg'nin hiçbir sonucu (hiçbir tohum, hiçbir koşul).
- **Gerekçe:** taslağın sınırlılık (vii) maddesi: uyarlanır kapı, kapılı bir genel (türbülans fiziği kullanmayan)
  uyarlamalı gürültü yöntemiyle karşılaştırılmadı. Q1 hakemi "kazancı fiziksel kestirim mi, yoksa herhangi bir
  uyarlamalı ölçekle kapıyı genişletmek mi taşıyor" diye sorar. Kullanıcı 2026-11-12'de gerekli ek koşumu onayladı.
- **Kol:** **AGg** = `AG --turb_cal_factor 12.0 --turb_gate 3`: AG'nin ana ızgaradaki yapılandırması (σ² = σ₀² + r̂²;
  r̂ = arka-uç yeniden izdüşüm artıklarının kayan RMS'i, pencere 60, 10 karede bir güncelleme, en az 50 bileşen;
  ilk kestirime kadar r̂ = 0) + KAPI ile aynı kapı kuralı (son graf karesinin bireysel-faktör sigması motion-only
  BA sigmasını, üç χ² sınamasını ve PnP eşiğini ölçekler). **Kod değişikliği yok:** `--turb_gate` A2 dışındaki her
  satır-içi kolda çalışır; parmak izi `bedca2849bb2` aynı kalır. AG'nin hiçbir parametresi ayarlanmaz; başka kapılı
  genel varyant (ör. uyarlamalı çekirdek) koşulmaz.
- **Izgara:** Sapma 6 ile aynı: 9 sekans, tohum 1–4, R0 D/r₀ {8, 16}, R5 D/r₀ {8, 4}, temiz; koşul başına 36 birim,
  örneklem sabit. Parça I (R5): Sapma 5'in saklı dökümleri, yalnız AGg (72 koşum). Parça II (R0): D/r₀ 16 saklı
  döküm, D/r₀ 8 ön-uç yeniden üretimi; A0 + AGg (144). Parça III (temiz): ön-uç yeniden üretimi; A0 + AGg (72).
  Toplam 288 arka-uç koşumu. Koşucu `pyeco_turb_kapi/sapma7_kos.py` (`kapi_dogrula.py` kopyası; yalnız kol kümesi,
  çıktı dizini, AGg denetimi değişti), kayıtlar `evidence/SAPMA7/s7_{I,II,III}.jsonl`.
- **Başvurular ve eşleme:** birim (sekans, tohum); aynı (sekans, tohum, D/r₀, rejim, pencere) ve aynı ön-uç
  çıktısı. A0: eski kayıtlar (R0 ana ızgara, R5 hızlı ızgara). KAPI ve A1g: Sapma 6 kayıtları. Kapısız AG: ana ızgara
  (R0, temiz) ve hızlı ızgara (R5). Aynı ön-uç kanıtı: parça II ve III'te her hücrenin A0 koşumu başvuru A0 ile
  bit-aynı olmalı (108 hücre); parça I'de döküm Sapma 6'da doğrulanan dosyanın kendisidir. Tutmazsa koşucu durur.
- **İstatistik:** Sapma 6 ile aynı (eşli medyan oran, birim bootstrap %95 GA 10 000, iki yönlü Wilcoxon, Holm aile
  içinde, TOST ±%15, α = 0.05). Betik `sapma7_analiz.py`, ızgara sonuçlarından önce yazılır; tam veriyle bir kez.
- **Hipotezler:**
  **HG1 (m = 4, Holm):** AGg / A0 ATE < 1 ve anlamlı: R0-8, R5-8, R5-4, R0-16. (Kapı genişletmenin kendisi, ölçeğin
  kaynağından bağımsız olarak doğruluğu iyileştirir.)
  **HG2 (m = 4, Holm; iki yönlü, yön öngörüsü YOK):** KAPI / AGg ATE, aynı dört koşul; ayrıca TOST ±%15.
  **HG3 (temiz):** AGg / A0 ATE %95 GA üst sınırı < 1.10.
  **İkincil:** AGg / AG (kapısız) ATE (m = 4); konum ANEES'i AGg / A0, KAPI / AGg, AGg / AG (m = 12); AGg / A1g ATE
  (m = 4); temizde KAPI / AGg ATE ve AGg / A0 ANEES (betimsel); kapı sigması dağılımı; başarısız pencere sayısı.
- **Sayısal beklentiler (sonuçtan önce):**
  - AGg kapı sigması medyanı: kapı geri beslemesi yoksa √(1 + r̂²) = 1.65 (R0-8), 1.82 (R5-8), 1.48 (R5-4),
    1.89 (R0-16), 1.19 (temiz) px. Geniş kapı daha büyük artıkları içeri aldığı için r̂'nin büyümesi beklenir
    (pozitif geri besleme): öngörü R0-8 1.7–2.8, R5-8 1.9–3.2, R5-4 1.5–2.2, R0-16 1.9–3.6, temiz 1.15–1.35 px.
    En küçük kapı sigması her koşumda tam 1.0 (ilk kestirimden önce).
  - **HG1:** AGg / A0 R0-8 0.45–0.85; R5-8 0.33–0.75; R5-4 0.58–0.95; R0-16 0.55–0.95. Dördünün de < 1 çıkması
    beklenir (A1g ve KAPI ikisi de kazanıyor; ölçek kaynağı ikincil görünüyor). R5-4 ve R0-16 en zayıf.
  - **HG2:** KAPI / AGg 0.80–1.25, yön belirsiz. İki karşıt mekanizma: (a) AGg'nin artık RMS'i kapıyla büyüyüp
    oracle'a yaklaşabilir (o zaman AGg KAPI'den iyi olur; A1g KAPI'den %5–42 iyi); (b) artık RMS'i poz hatasını ve
    aykırı eşleşmeleri de içerir, kapı kontrolsüz genişleyip yanlış eşleşmeleri alabilir (o zaman KAPI iyi olur).
    Benim nokta tahminim: 0.95–1.10 (ayırt edilemez ya da AGg hafif iyi).
  - **HG3:** medyan 0.97–1.05; üst sınır < 1.10 geçme olasılığı ≈ %75 (Sapma 6'daki duyarlılık).
  - ANEES: AGg / A0 R0-8'de 0.05–0.40; KAPI / AGg 0.6–2.0 (yön belirsiz). Başarısız pencere ≤ 3 / 108 koşul başına.
- **Yorum kuralları (sonuçtan önce):**
  - HK1–HK4 kararları (Sapma 6) bu sapmayla değişmez; AGg ek bir karşılaştırma koludur.
  - HG1 tutar ve HG2'de fark anlamsızsa (ya da TOST eşdeğer): makalede "doğruluk kazancını kapının ölçeklenmesi
    taşır; ölçeğin yapı-fonksiyonu kestiriminden ya da genel artık RMS'inden gelmesi doğrulukta ayırt edilemedi"
    yazılır. Fiziksel kestirimin katkısı doğruluk üstünlüğü olarak SUNULMAZ; katkı "r₀/σ_t'yi fiziksel birimle
    veren, gerçek sintilometreyle doğrulanan kestirim + kapı ölçekleme mekanizmasının gösterilmesi" olarak kalır.
    Özet ve katkı listesi buna göre düzeltilir.
  - KAPI / AGg < 1 ve anlamlıysa: "fizik tabanlı ölçek genel ölçekten daha doğru" etki büyüklüğüyle yazılır.
  - KAPI / AGg > 1 ve anlamlıysa: açıkça raporlanır ("genel artık ölçeğiyle kapı daha doğru"); yöntem bölümünde
    kapı kuralı ölçek kaynağından bağımsız sunulur, fiziksel kestirimin doğruluk üstünlüğü iddiası çıkarılır ve bu
    sonuç özet düzeyinde belirtilir.
  - HG1 tutmazsa (AGg kazanmaz): "kapının kazancı ölçeğin kalitesine bağlıdır" diye, GA ile yazılır.
  - HG3 kalırsa: AGg'nin temizde bedeli olduğu yazılır.
  - Sonuç ne olursa olsun AGg tabloya girer; sınırlılık (vii) "genel kapılı kol koşuldu; uyarlamalı çekirdek ve
    dağılım-tabanlı eleme koşulmadı" diye daraltılır.
- **Ara bakış kuralı:** koşum sürerken yalnız koşucu günlüğü (durum, A0 bit-aynılığı, kapı sigması istatistiği) ve
  kayıt sayısı; ATE / ANEES oranları yalnız 288 kayıt-kümesi + başvurular tamken, bir kez. Erken durdurma yok.
- **Ön test (tohum 0, V1_01_easy; dört hücre R0-8, temiz, R5-8 ön-uç yeniden üretimiyle; R0-16 saklı dökümle; 9
  arka-uç koşumu):** K1 A0 ≡ başvuru A0 (4 hücre × 3 pencere); K2 ön-uç istatistikleri başvuruyla aynı; K3 R0-8'de
  KAPI yeniden koşumu Sapma 6 geliştirme kaydıyla bit-aynı + yörünge sha256 (koşucu kopyasının ve kodun değişmediği);
  K4 AGg: `gate == 3`, kapı sigması sayısı = kare − 1, en küçük = 1.0, en büyük > 1.0 (etki uygulandı), bütün metrik
  alanları dolu; K5 analiz betiği ön test kayıtlarını okuyor. Yumuşak beklentiler: AGg kapı sigması medyanı R0-8
  1.5–2.8, temiz 1.1–1.4 px; AGg arka-uç süresi 300–900 s. (Tohum 0 ön test oranları tek hücredir, yorumlanmaz.)
- **Maliyet:** Sapma 6'da ölçülen medyan süreler (A0 566–746 s, kapılı kol 686–1218 s) → ≈ 65 işçi-saat; 8–10 işçiyle
  ≈ 7–9 saat. Ön-uç yeniden üretimi 72 hücre ≈ 2.4 GPU-saat (yalnız RTX 3090). Disk: yeni döküm ≈ 6 GB (boş 67 GB).
- **Makalede kullanım:** AGg sonuçları Sapma 7 etiketiyle, doğrulayıcı ızgaranın eki olarak raporlanır.
