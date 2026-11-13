# AirSim kolları: H3 (AirSim testi) ve HA4 sonuçları (ON_KAYIT_FAZ2 §5 Sapma 1)

- **Tarih:** 2026-11-09 06:16. Izgara 110/110 kayıt, hepsi `ok`: 90 kol koşumu ve 20 artık dökümü.
  - Hepsi RTX 3090'da koştu; code_fp `08e28b83407c`.
  - A0, ana ızgarayla bayt bayt aynı (30/30).
  - Koşum başına NEES örnek sayısı en az 378, medyan 401.
- **Analiz:** `airsim_sim/analiz_kollar.py hepsi --etiket izgara` (sha256 c6d471c2…a261).
- **Çıktı:** `cikti_kollar/izgara/analiz.json` (sha256 28b6896b…89ce). Analiz yeniden koşturuldu, çıktı bayt bayt aynı.
  - Aşağıdaki sayılar bu dosyadan alındı. İstisna: §3'teki "A0 mutlak ANEES" sütunu `kayit.jsonl`'dan doğrudan
    hesaplandı ve betimseldir.
- **Tanımlar:** `africa/SAPMA_TASLAK_KOLLAR.md` §2. Koordinatör kararları aynı dosyanın §7'sinde; sonuçlar görülmeden
  önce verildi.

## 1. H3, AirSim testi: öngörü TUTMADI (bandın altında, eksik güven)

| birim | n | medyan ANEES_poz | boot %95 | bant (n medyan 404) | bant üstü | bant içi |
|---|---|---|---|---|---|---|
| A0, C0 = 0, F5, tohum 10–14 × W1/W2 | 10 | 1.490 | 1.12–1.87 | 2.77–3.24 | 0/10 | 0/10 |

- **Sonuç:** 10 birimin hepsi bandın altında; medyan ANEES_dön 1.195.
  - W1 (20:857): 1.76–2.07.
  - W2 (1157:1994): 1.04–1.22.
- **Karar (§7, sonuç öncesi):** "öngörü tutmadı". Aşırı güven AirSim'de yok, kestiriciye içkin değil.
  - Bu, çürütmeyle aynı sonuca bağlanır: EuRoC'taki aşırı güven veriye/GT'ye özgüdür.
- Pilot (1.09 / 1.81) ile tutarlı.

## 2. HA4, derinlik bağımlılığı: öngörü TUTTU

- **Birincil ölçüt:** ek artık Δσ² (F5, 1e-11 ile temiz karşılaştırması), GT pozlarıyla üçgenleme, eşit sayılı 8 dilim.
  - Bütün dilimler 10/10 çiftte geçerli.

| dilim (m) | 25–65 | 65–80 | 80–90 | 90–105 | 105–115 | 115–130 | 130–145 | 145–400 |
|---|---|---|---|---|---|---|---|---|
| medyan ek artık (px) | 0.365 | 0.420 | 0.430 | 0.455 | 0.505 | 0.545 | 0.595 | 0.656 |
| σ_rob temiz (px) | 0.231 | 0.231 | 0.238 | 0.258 | 0.270 | 0.293 | 0.322 | 0.302 |
| σ_rob türbülanslı (px) | 0.427 | 0.479 | 0.490 | 0.525 | 0.569 | 0.623 | 0.680 | 0.721 |

- **Sonuç:** ρ = 1.000 (> 0.8, öngörü tuttu). Çift başına ρ'nun medyanı 1.000; 10 çiftin 10'unda ρ > 0.8.
- **İkincil sonuçlar:**
  - Model dilimleri (25/50/100/200 m): dolu olan 3 dilimde ρ = 1.000; çift başına ρ > 0.8 olan çift 7/10.
  - A0 kestirim pozlarıyla artık: ρ = 0.976; çift başına ρ'nun medyanı 0.869 (6/10).
  - Model σ_tot: 0.40 / 0.50 / 0.59 / 0.91 px (25–50 / 50–100 / 100–200 / 200–400 m), bilgi amaçlı.
    Ölçülen ek artık aynı büyüklükte ve aynı yönde.

## 3. Kollar (betimsel; v0.2'de yönlü öngörü yok)

| C0 | A0 mutlak ANEES* | AG ANEES | AG/A0 oranı | A3 ANEES | A3/A0 oranı | A3/AG oranı | RPE1s oranı (kol/A0) |
|---|---|---|---|---|---|---|---|
| 0 | 1.490 | 1.193 | 0.800 | 0.943 | 0.653 | 0.807 | 1.000 / 1.007 |
| 1e-12 | 1.745 | 1.375 | 0.779 | 1.084 | 0.596 | 0.759 | 1.014 / 1.013 |
| 1e-11 | 3.958 | 2.794 | 0.694 | 1.994 | 0.497 | 0.720 | 0.998 / 1.000 |

\* `kayit.jsonl`'dan doğrudan hesaplandı (betimsel); öteki sütunlar `analiz.json`'dan.

- **Testler:** her kol × C0 için iki yönlü Wilcoxon; Holm sonrası p = 0.0117. Bu, n = 10'da elde edilebilecek en küçük
  değer: 10/10 birim aynı yönde.
- **Gözlemler** (keşifsel, hipotez değil):
  - Temiz ve 1e-12'de A0 zaten eksik güvenli. AG ve A3 kovaryansı daha da büyütüyor, ANEES'i banttan uzaklaştırıyor.
  - 1e-11'de A0 aşırı güvenli (3.96, bandın üstünde). AG ANEES'i 2.79'a, yani bandın alt kenarına çekiyor.
    A3 ise fazla düzeltiyor (1.99).
  - Kolların doğruluğa etkisi yok: RPE1s oranı 0.998–1.014.
  - Yorum: kollar türbülansın eklediği güven açığını kapatıyor. Temiz veride bulunan eksik güvenli taban hâli ise
    düzeltmenin üstüne biniyor. EuRoC'ta tabanın aşırı güvenli olduğu düşünülürse bu fark makalede tartışılacak.
