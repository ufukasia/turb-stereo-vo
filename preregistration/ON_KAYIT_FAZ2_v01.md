# Faz 2 ön-kayıt (koşumlardan ÖNCE yazıldı) — TASLAK v0.1, 2026-11-07

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
- **Tohumlar:** 0–4 (türbülans alanı). İstatistik birimi (sekans, tohum); pencereler birim
  içinde ortalanır.
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

## 4. AirSim kolu (önkoşul: derinliğe bağlı türbülans modu)

- **Model:**
  - Her piksel, render derinliğinden elde edilen eğik görüş hattı boyunca hesaplanır.
  - Cn²(h) ∝ h^(-4/3) profili, küresel dalga r₀ (z nesneden ölçülür).
  - σ_tilt piksel başına değişir; açısal izoplanatik bölge derinlik ve yükseklikle değişir.
  - Stereo taban çizgisi 0.40 m. İki kameranın tilt alanları, örtüşen hacim oranına göre
    kısmi korelasyonlu.
- **Rejimler (TASARIM.md §1):**
  - N: nadir, h = 5–10 m. Negatif kontrol.
  - O: h = 3–10 m, eğik −20…−35°. Ana senaryo.
  - F: ileri bakış, dar FOV. Üst sınır.
- **Cn² seviyeleri:** 10⁻¹⁴, 10⁻¹³, 10⁻¹² m^(-2/3) (zayıf, orta, güçlü yer-yakını gündüz).
  Literatür aralığı kaynakla birlikte verilir.

**HA1 — Negatif kontrol.**
- Öngörü: N rejiminde en güçlü Cn²'de bile medyan tilt < 0.1 px ve ATE farkı anlamsız
  (eşdeğerlik testi, TOST ±%10).
- Çürütme: N rejiminde anlamlı bozulma. Bu durumda model veya hat hatası aranır.

**HA2 — Eğik geometri.**
- Öngörü: O rejiminde etki derinliğe göre monoton artar. Uzak bölgedeki (Z > 20 m) eşleşmeler
  bozulmanın çoğunu taşır. A3/AG'nin ANEES kazancı N'de yoktur, O'da vardır.
- Çürütme: O ile N arasında fark yok.

**HA3 — Metrik ölçek.**
- Öngörü: temiz hatta ölçek hatası < %0.5 (pilot: Sim(3) ölçeği ≈ 1). Türbülans ölçeği değil,
  gürültüyü etkiler; en güçlü Cn²'de ölçek hatası < %2.
- Çürütme: ölçek hatası > %2.

---

## 5. Dondurma

- [ ] Kullanıcı onayı → dosya sha256'sı `evidence/on_kayit_sha256.txt`'ye yazılır, sonra koşum başlar.
- Sapma kayıtları: (boş)
