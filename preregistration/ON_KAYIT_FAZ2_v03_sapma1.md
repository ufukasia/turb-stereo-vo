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
