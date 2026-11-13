# Sapma 7 — kapılı genel uyarlamalı kol (AGg): sonuç

- **Izgara:** parça I 13:24:38, parça II 13:29:45, parça III 13:08:23 (2026-11-13) bitti; 288/288 arka-uç koşumu,
  864 pencere kaydı, sorun 0, ön-uç hatası 0. Kod parmak izi `bedca2849bb2`.
- **Analiz:** `sapma7_analiz.py`, tam veriyle bir kez, 2026-11-13 13:29 (`s7_analiz.json`, günlük `s7_analiz.log`).
  Ara bakışlar yalnız bütünlük (durum, A0 bit-aynılığı, kapı sigması) ile sınırlı kaldı; ATE/ANEES oranına
  analizden önce bakılmadı.
- **Bütünlük:** 2916/2916 kayıt (eksik 0). V1 A0 yeniden koşumu ≡ başvuru A0: 324/324 pencere bit-aynı. V2 ön-uç
  istatistikleri: 216/216 aynı. V3 AGg kapı alanları (gate = 3, en küçük sigma 1.0, en büyük > 1.0): 540/540.
  V5 bayraklar 864/864.
- **Yürütmedeki sapma:** ön testin ilk iki denemesi K1/K3'te (yalnız NEES ortalama/medyan, göreli ~1e-13) kaldı;
  neden OMP_NUM_THREADS'in başvuru koşumlarındaki 4 değerine ayarlanmamasıydı. Üçüncü deneme ayarla, kapı
  gevşetilmeden geçti (`ON_TEST.md`). Başarısız denemeler saklı.

## Kararlar (eşli medyan oran, birim bootstrap %95 GA, Wilcoxon, Holm aile içinde; n = 36 birim)

| Hipotez | Koşul | Oran [%95 GA] | p_Holm | TOST ±%15 | Karar |
|---|---|---|---|---|---|
| HG1 AGg / A0 | R0-8 | 0.542 [0.485, 0.700] | 2.62e-10 | – | TUTTU |
| HG1 AGg / A0 | R5-8 | 0.440 [0.414, 0.542] | 1.16e-10 | – | TUTTU |
| HG1 AGg / A0 | R5-4 | 0.672 [0.611, 0.709] | 7.36e-09 | – | TUTTU |
| HG1 AGg / A0 | R0-16 | 0.849 [0.809, 0.891] | 8.15e-10 | – | TUTTU |
| HG2 KAPI / AGg | R0-8 | 0.890 [0.842, 1.001] | 0.0314 | ort 0.909 GA90 [0.858, 0.963], p 0.101, eşdeğerlik yok | KAPI daha doğru |
| HG2 KAPI / AGg | R5-8 | 0.895 [0.824, 0.951] | 0.0314 | ort 0.917 GA90 [0.871, 0.966], p 0.047, eşdeğer | KAPI daha doğru |
| HG2 KAPI / AGg | R5-4 | 0.952 [0.923, 1.023] | 0.108 | ort 0.955 GA90 [0.914, 0.997], p 4.9e-4, eşdeğer | ayırt edilemedi; TOST eşdeğer |
| HG2 KAPI / AGg | R0-16 | 0.953 [0.903, 0.988] | 0.0111 | ort 0.953 GA90 [0.928, 0.978], p 4.6e-7, eşdeğer | KAPI daha doğru |
| HG3 AGg / A0 | temiz | 0.989 [0.956, 1.023] | 0.822 | – | GEÇTİ (üst sınır 1.023 < 1.10) |

Not: R0-8'de HG2 Wilcoxon'la anlamlı, ancak medyan oranın bootstrap GA'sı 1.001'e uzanıyor; bu ayrım metinde
belirtilir. R5-8 ve R0-16'da fark hem anlamlı hem ±%15 içinde eşdeğer: küçük ama sistematik.

## İkincil

- AGg / AG (kapısız) ATE: 0.540 / 0.449 / 0.636 / 0.820 (R0-8 / R5-8 / R5-4 / R0-16), hepsi p_Holm < 1e-9.
- AGg / A1g ATE: 1.317 / 1.340 / 1.050 / 1.447 (oracle kapı hâlâ daha iyi).
- Konum ANEES: AGg / A0 0.209 / 0.109 / 0.472 / 0.286; KAPI / AGg 0.819 / 0.743 / 0.843 / 0.794 (hepsi anlamlı);
  AGg / AG 0.540 / 0.332 / 1.037 / 1.150 (son ikisi anlamsız).
- Temiz: KAPI / AGg ATE 1.019 [0.999, 1.037] (p_Holm 0.304; TOST eşdeğerlik gösterilemedi, GA90 [0.970, 1.179]);
  AGg / AG ATE 1.006 (eşdeğer); ANEES AGg / A0 0.746, KAPI / AGg 0.915.
- Kapı sigması medyanı (AGg): 1.792 / 1.974 / 1.549 / 2.023 / temiz 1.195 px; KAPI 2.011 / 2.232 / 1.702 / 2.311 /
  1.253 px. Başarısız pencere AGg: R0-16'da 5/108 (A0 5/108, KAPI 5/108, AG 11/108), diğerlerinde 0.

## Sonuçtan önce yazılan beklentilerle karşılaştırma

| Beklenti | Değer | Sonuç |
|---|---|---|
| HG1 bantları (0.45–0.85 / 0.33–0.75 / 0.58–0.95 / 0.55–0.95) | 0.542 / 0.440 / 0.672 / 0.849 | dördü de bant içinde |
| HG2 bandı 0.80–1.25 | 0.890 / 0.895 / 0.952 / 0.953 | bant içinde |
| HG2 nokta tahmini 0.95–1.10 (ayırt edilemez ya da AGg hafif iyi) | R0-8 0.890, R5-8 0.895 | **tutmadı** (KAPI daha iyi); R5-4, R0-16 bant içinde |
| HG3 medyan 0.97–1.05, üst sınır < 1.10 | 0.989, 1.023 | tuttu |
| Kapı sigması (1.7–2.8 / 1.9–3.2 / 1.5–2.2 / 1.9–3.6 / 1.15–1.35) | 1.792 / 1.974 / 1.549 / 2.023 / 1.195 | hepsi bant içinde (alt uca yakın) |
| ANEES AGg / A0 R0-8 0.05–0.40 | 0.209 | tuttu |
| ANEES KAPI / AGg 0.6–2.0 | 0.743–0.843 | tuttu |
| Başarısız pencere ≤ 3/108 koşul başına | R0-16: 5/108 | **tutmadı** (A0 ve KAPI de 5/108) |

## Ön-yazılı yorum kuralına göre okuma

HG1 dört koşulda tuttu: kapıyı genel bir artık ölçeğiyle büyütmek de doğruluğu iyileştiriyor; kazancın büyük
kısmını kapının ölçeklenmesi taşıyor. HG2'de KAPI / AGg üç koşulda < 1 ve anlamlı: kural gereği "fizik tabanlı
ölçek genel ölçekten daha doğru" etki büyüklüğüyle yazılır: medyan %5–11 (R0-8 0.890, R5-8 0.895, R0-16 0.953),
R5-4'te ayırt edilemedi (eşdeğer). Fark küçük; R5-8 ve R0-16'da ±%15 eşdeğerlik sınırı içinde. KAPI ayrıca
tutarlılıkta açık biçimde daha iyi (konum ANEES %16–26 düşük). HK1–HK4 kararları değişmez. Sınırlılık (vii):
"genel kapılı kol koşuldu; uyarlamalı çekirdek ve dağılım tabanlı eleme koşulmadı" diye daraltılır.
