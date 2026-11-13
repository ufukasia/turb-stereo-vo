
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
