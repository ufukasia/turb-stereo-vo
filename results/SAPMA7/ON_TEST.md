# Sapma 7 ön testi (AGg kolu) — 2026-11-13 07:24:00

Kimlik: code_fp bedca2849bb2 (ana ağaç 08e28b83407c); koşucu sapma7_kos.py sha256 f5fe8f902b338b01;
SECIM.json cd4f4520183e0ee4; ortam CUDA_DEVICE_ORDER=PCI_BUS_ID, CUDA_VISIBLE_DEVICES=1, OMP_NUM_THREADS=4.
Girdi: V1_01_easy, tohum 0 (keşif; kanıt değildir), dört koşul, 9 arka-uç koşumu, 27 kayıt.

| Kapı | Beklenen (koşumdan önce yazıldı) | Ölçülen (deneme 3) | Sonuç |
|---|---|---|---|
| K1 A0 başvuruyla bit-aynı (ATE, RPE, NEES poz/rot) | 12/12 pencere | 12/12 | geçti |
| K2 ön-uç istatistikleri başvuruyla aynı | 3/3 yeniden üretilen hücre | 3/3 | geçti |
| K3 KAPI geliştirme kaydıyla bit-aynı + yörünge sha256 | 3/3 | 3/3, sha256 ve gerçek-zaman aynı | geçti |
| K4 AGg alanlar dolu, gate=3, kapı sigması tabanı 1.0, büyüyor | 12/12 ok | 12/12 ok; medyan 1.677 / 1.952 / 1.755 / 1.120 (R0-8 / R0-16 / R5-8 / temiz) | geçti |
| K5 (yumuşak, karar kuralı değil) AGg/A0 ATE | temizde ~1, türbülansta < 1 | 1.012 / 0.946 / 0.741 / 0.597 | bilgi |

Süre: arka-uç 266–558 s (medyan 333 s), ön test toplam 14.8 dk. Izgara: 288 arka-uç koşumu, 72 hücrede ön-uç
yeniden üretimi (kuru çıktı: II 36/72 saklı, III 0/36 saklı, I 72/72 saklı); Sapma 6 aynı düzenle ~13 saat sürdü.
Disk / 66 GB boş; bellek 41 GB kullanılabilir.

## Başarısız denemeler (saklandı)
Deneme 1 (jobs 9) ve deneme 2 (jobs 5) OMP_NUM_THREADS ayarsız başlatıldı (başlatma hatası; Sapma 6 koşumlarının tümü
OMP_NUM_THREADS=4 idi). Sonuç: K1/K3 yalnız NEES ortalama/ortanca alanlarında bağıl ~1e-13 fark (ATE, RPE, yörünge
sha256 aynı); iki deneme birbirine bit-aynı. Arka-uç 7175–19137 s. Kayıtlar s7_ontest_deneme{1_jobs9,2_jobs5}.{jsonl,log}.
Kapı gevşetilmedi; ortam düzeltilip yeniden koşuldu.
