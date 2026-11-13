# Keşif sondası — ön-test kaydı (2026-11-11, koşum öncesi)

Ön beyan: `ON_BEYAN.md` (sha256 ffad5722…, `ON_BEYAN.sha256`). Bu dosya zorunlu ön-test denetimlerinin sonucudur
(kural: uzun koşum ön-test kuralı). Sonda sonuçlarına bakılmadan yazıldı.

| # | denetim | beklenen | ölçülen | sonuç |
|---|---|---|---|---|
| 1 | kod parmak izi | 08e28b83407c (ana ızgara) | 08e28b83407c | GEÇTİ |
| 1 | config.yaml sha256 | 92e1dc7f…948aa1 (ana ızgara) | 92e1dc7f…948aa1 | GEÇTİ |
| 1 | GPU | RTX 3090 (`fe_gpu`) | RTX 3090 | GEÇTİ |
| 2 | bit-aynılık kapısı (E1): V1_01_easy, tohum 0, D/r0 = 8, R0/fe0, pencereler 100/500/900, A0 + A1 | ana ızgarayla bit-aynı | 6/6 kayıt, 35/35 metrik alanı bit-aynı (`e1_test.jsonl`, etiket `e1tam`) | GEÇTİ |
| 3 | duman: 1 sekans × 1 pencere × {R0, R5} × {fe0, fe1} × {A0, A1, A0g} | çökme yok, bütün alanlar dolu (E5) | 12/12 kayıt `ok`, NEES örnek 28–59 | GEÇTİ |
| 4 | E2: fe1 üçgenleme verimi (D/r0 = 8, R0) | ≥ 0.88 | 0.839 → 0.999 | GEÇTİ |
| 4 | E3a: R5 enine hız | 5.000 m/s | 5.000 m/s (kuru koşum) | GEÇTİ |
| 4 | E3b: R5/fe0 zamansal eşleşme (D/r0 = 8) | < 100 (9 sekans ortalaması; R0 112.6) | tek sekans: 147.3 → 130.0 (oran 0.88) | sondada karara bağlanacak |
| 4 | E4: ANEES(A0g)/ANEES(A0), R0, D/r0 = 8 | < 0.3 ve A1 oranının 3 katı içinde | tek pencere: 0.73; A0g/A1 = 3.1 | TUTMADI (tek pencere); sondada karara bağlanacak |
| 5 | çıktı / analiz betiği | `sonda_analiz.py` kayıtları okur, tabloları basar | 12 kayıt okundu, bütün tablolar dolu | GEÇTİ |
| 6 | süre | ON_BEYAN: 2–2.5 sa | yeniden tahmin: ön-uç 103 336 kare × ≈ 0.15–0.2 s / 3 eşzamanlı ≈ 1.5–2 sa; arka-uç ≈ 38 CPU-sa / 12 işçi ≈ 3–3.5 sa | toplam ≈ 3.5 sa |
| 6 | disk / bellek | yeterli | /tmp 93 GB boş (en kötü durum ≈ 36 GB önbellek + spool); RAM 36 GB kullanılabilir | GEÇTİ |

## E1 farkının açıklaması (duman testinde 4. anlamlı basamak farkı)

Duman testi yalnız ilk pencereyle (500 kare) koştu ve `m_ate_se3` ana ızgaradan farklı çıktı
(A0: 0.06956956 ↔ 0.06954797). Ayrıştırma testi (`e1_test.jsonl`):

- **Tam pencere kümesi** (1300 kare, GPU'da tek başına): ana ızgarayla 6/6 kayıt bit-aynı.
- **Kırpılmış** (500 kare, tek başına): duman testiyle (3 eşzamanlı ön-uç) bit-aynı → GPU paylaşımı bit-aynılığı bozmuyor.
- Ön-uç dökümü önek-aynı: tam dökümün ilk 500 karesi kırpılmış dökümle bayt-aynı.
- Neden: arka-uç koşumun **son karesini** her zaman eşzamanlı optimize eder (`f == n_frames - 1`; `be_aralik = 3`,
  `be_paralel = 2`). Kırpılmış koşumda kare 499 son karedir; tam koşumda ara karedir. Fark yalnız koşumun son
  karesinin nedensel pozundadır; nedensellik ihlali değildir.

Sonuç: E1 yalnız tam pencere kümesiyle sınanır. Sonda tam pencere kümesiyle koşar; ilk biten R0/fe0 hücrelerinde
E1 yeniden denetlenecek, fark çıkarsa sonda durdurulacak.

## E4 notu

A0g ile A1 aynı faktör kovaryansını kullanır (σ = 4.524 px, D/r0 = 8); fark yalnız kapılardadır. Geniş kapı daha çok
ölçüm kabul eder → sonsal kovaryans küçülür → ANEES A1'den yüksek kalır. Beklenti yanlış kurulmuş olabilir; 9 sekanslık
sondada değerlendirilecek. Bu bir durdurma kuralı değildir (yalnız E1 durdurma kuralıdır).

## Koşum komutu

`nice -n 5 ../kod/turbvio/bin/python sonda_kesif.py --gpu_jobs 3 --jobs 12` (ortam `faz2_kos.sh` ile aynı;
`TURB_DUMP_CACHE=/tmp/turb_cache_sonda`). Çıktı `records.jsonl`, günlük `run_*.log`.
