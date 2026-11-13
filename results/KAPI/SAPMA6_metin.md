
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
