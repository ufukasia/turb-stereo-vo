# Ön-uçuş kapısı

Üretim: 2026-11-08 01:15 · sonuç: **KALDI**

**P** maddeleri yazılımın MEKANİĞİNİ denetler (her biri yaşanmış bir kusura karşılık gelir). **E** maddeleri ÖN-KAYITLI BİLİMSEL BEKLENTİLERdir: koşum mekanik olarak kusursuz olup tezi çürütebilir ve bu, günler süren kuyruk başlamadan ÖNCE bilinmelidir. Bir E maddesi düşerse sonuç 'hata' değil **DUR-VE-İNCELE**dir — beklenti de yanlış olabilir, yöntem de.

| # | durum | denetim | ölçülen | beklenen | hangi olaydan doğdu |
|---|---|---|---|---|---|
| P1 | ✅ | çekirdek kolların hiçbiri çökmüyor (A4 keşifsel — raporlanır) | A0:%0; A0C:%0; A0C+A3:%0; A1:%0; A3:%0; A5:%0; AG:%0  [A4 keşifsel: %0] | çekirdek kollar < %10 | A4 %47-79 ıraksıyordu |
| P2 | 🛑 | kovaryans verimi (pencere başına örnek) | medyan 11 / beklenen 40 (126 pencere) | >= 30 | %89 sessiz başarısızlık (Cholesky) |
| P3 | ✅ | kovaryans PD kapısı arka-uçta etkin | gtsam_backend._pose_marginal_schur PD kontrolü | PD olmayan matris döndürülmez | 6979'un 138'i bozuktu (A3) |
| P4 | ✅ | NEES sonlu ve pozitif | negatif 0, sonsuz/NaN 0 / 126 | 0 / 0 | tabloya -5977 düşmüştü |
| P5 | ✅ | her kol GERÇEKTEN farklı bir yöntem mi | 10 ikili, hepsi farklı | her ikilide medyan bağıl fark > 1e-3 | A0C, cal_factor=1.0 yüzünden A0 ile aynıydı (plasebo) |
| P6 | ✅ | DIMM theta_c kestirimi fiziksel aralıkta | 113-113 px (görüntü 752 px) | [8, 1504] px | 58-1460 px arası salınıyordu |
| P7 | ✅ | her kayıt güncel kodla damgalı | {'65dbe4761b2b': 126} (güncel 65dbe4761b2b) | hepsi 65dbe4761b2b | code_fp hatlı kod yoluna eklenmemişti |
| P8 | ✅ | hiçbir kol saçma ATE üretmiyor | en büyük ATE 0.63 m | < 10 m | A3, D/r0=8'de 17.9 m veriyordu |
| P9 | ✅ | NEES hatası kovaryansla AYNI teğet uzayında | retract == compose(Expmap) | retract == compose(Expmap) | hata ham ötelenme alınıyordu; |w|=1 rad'da %40 sapma |
| P10 | ✅ | A1 ile A5 AYNI köşegeni taşıyor (fark yalnız yapıda) | köşegen bağıl fark 6.2e-11; köşegen-dışı max |rho| 0.362 | köşegen farkı < 1e-6 VE köşegen-dışı yapı var (|rho| > 0.01) | nugget iki kez sayılınca A5'in köşegeni A1'den ayrışmıştı |
| P11 | 🛑 | BÜTÜN kollar aynı aykırı-değer çekirdeği rejiminde | hiçbir kol robust sayacı bildirmedi (sayaç eklenmemiş?) | her kolda robust faktör oranı AYNI | yama ve tekil yollar Huber'ı atlıyordu; A1'de vardı, A5'te yoktu |
| P12 | 🛑 | Schur tümleyeni BAĞLI landmark'larda da doğru | sınanamadı: TypeError: 'gtsam.gtsam.Ordering' object is not iterable | yoğun tersle bağıl fark < 1e-9 | Hll blok-köşegen VARSAYILIYORDU; yama faktörü onu bağlıyor, yapılı kolların kovaryansı fazla küçük çıkıyordu |
| P13 | ✅ | D/r0=0'da A1 ile A5 BİREBİR aynı (yapı terimi yok) | 6 temiz hücrede A1 vs A5 ATE, max bağıl fark 8.19e-14 | bağıl fark < 1e-3 | Huber granülerliği A5'i temiz veride 2.26 kat kötüleştiriyordu |
| P14 | ✅ | çevrimiçi SIGMA bir ÖLÇÜM mü (sabit değil, oracle mertebesinde) | 12 kayıt · sigma_hat/oracle = 0.53x · bağıl saçılım 0.128 | saçılım > 0.02 VE sigma_hat/oracle 0.5-2.0 arasında | theta_c sınıra çakılıp sabit değer raporluyordu; artık sabit nominal kullanılıyor ve yalnız sigma kestiriliyor |
| E1 | ✅ | Türbülans-farkında kovaryans, temel çizgiden TUTARLI olmalı (NEES bandın merkezine YAKLAŞMALI) | A0/A3=5.28x  A0/A5=24.20x  (A0=1130 A3=204 A5=37)  [4 eşleşmiş hücre] | eşleştirilmiş oranlar > 1 | ön-kayıtlı beklenti (2026-08-19; eşleştirme 2026-08-21) |
| E2 | ✅ | YAPI, skaler şişirmeden İYİ olmalı — C3'ün ANA İDDİASI (aynı köşegen, fark YALNIZ köşegen-dışı yapıdan) | A1/A5 = 1.30x  (A1=46  A5=37)  [4 eşleşmiş hücre] | eşleştirilmiş A1/A5 >= 0.85 (teorik beklenti 1.2-1.5) | ön-kayıtlı beklenti (2026-08-19; eşleştirme 2026-08-21) |
| E3 | ✅ | Kalibre temel çizgi (A0C) tutarlılığı satın almalı ama DOĞRULUKLA ödemeli — yapının gerekçesi bu ödünleşimdir | A0/A0C = 57.36x, ATE 0.173 -> 0.210  [4 eşleşmiş hücre] | eşleştirilmiş A0/A0C > 1 | ön-kayıtlı beklenti (2026-08-19; eşleştirme 2026-08-21) |
| E4 | ✅ | Türbülans-farkında kollar DOĞRULUĞU bozmamalı (tutarlılık doğruluk pahasına satın alınmamalı) | A5/A0 ATE = 1.07x  (A0=0.173 m  A5=0.190 m)  [4 eşleşmiş hücre] | eşleştirilmiş A5/A0 ATE <= 1.5 | ön-kayıtlı beklenti (2026-08-19; eşleştirme 2026-08-21) |

## Kol özeti (bu küçük ızgarada)

| kol | n | ATE ort. | NEES medyan | kapsama | başarısız |
|---|---|---|---|---|---|
| A0 | 18 | 0.1565 | 734.9 | 0.002 | 0 |
| A0C | 18 | 0.1329 | 13.5 | 0.457 | 0 |
| A1 | 18 | 0.1436 | 78.5 | 0.134 | 0 |
| A2 | 0 | nan | nan | nan | 0 |
| A3 | 18 | 0.1430 | 145.8 | 0.066 | 0 |
| A4 | 0 | nan | nan | nan | 0 |
| A5 | 18 | 0.1396 | 65.4 | 0.167 | 0 |
