# Keşif sondası — ön beyan (2026-11-11)

**Statü:** KEŞİF. Ön kayıtlı (ON_KAYIT_FAZ2.md) doğrulayıcı analizlerin parçası DEĞİL; makalede kullanılırsa
"keşif (ön kayıt dışı)" diye etiketlenir. Bu belge sondanın hiçbir sonucu görülmeden yazıldı; özeti
`ON_BEYAN.sha256` dosyasında. Bilinen sonuçlar: Faz 2'nin bütün ön kayıtlı sonuçları (ana ızgara, H2, H6,
rakipler D/r₀ {0, 8, 16, 24}, AirSim kolları).

## 1. Neden (mevcut kanıttan çıkan iki soru)
1. **Sabit kapılar.** Ön-uç stereo satır kapısı 2 px sabit; sol-sağ diferansiyel dikey tilt D/r₀ = 8'de 1.4 px,
   16'da 2.3–2.6 px. Ana ızgarada üçgenleme verimi 0.965 (temiz) → 0.73 (8) → 0.47 (16) → 0.36 (24) düşüyor.
   Arka-uç kapıları (PnP 2.45 px, χ² kapıları σ₀ = 1 px) da sabit. Kollar yalnız düzleştirici kovaryansını
   değiştiriyor. Soru: kapılar türbülans ölçeğiyle (oracle) genişletilirse doğruluk ve sağlamlık değişir mi?
2. **Zamansal rejim.** Bütün Faz 2 hücreleri τ₀ = Δt ile koşuldu; bu, D/r₀ = 8'de 0.079 m/s, 16'da 0.039 m/s enine
   akış demek (tilt kareden kareye 0.09 px değişiyor, 1 s gecikmede korelasyon 0.65–0.90). Gerçekçi enine hız
   (rüzgâr + platform) m/s düzeyinde. Soru: hızlı rejimde kırılma ve kolların tutarlılığı değişir mi?

## 2. Tasarım (kod değişikliği YOK; yalnız `cfg_overrides`, `backend_argv`, `tau0_over_dt`)
- 9 sekans × tohum 0 × kapsam1m pencereleri; σ = 1 DN; NEES açık, cov_every = 2 (ana ızgarayla aynı).
- Rejim: **R0** τ₀/Δt = 1 (ana ızgara), D/r₀ ∈ {8, 16}; **R5** sabit enine hız v = 5 m/s
  (τ₀/Δt = 0.314·(D/(D/r₀))/(v·Δt) = 0.0314 / 0.0157 → D/r₀ = 4 / 8), D/r₀ ∈ {4, 8}.
- Ön-uç: **fe0** varsayılan; **fe1** `y_axis_threshold = 2·sqrt(1+σ²)` ve (yalnız R5)
  `temporal_ransac_threshold = 1.96·sqrt(1+σ²)`; σ = oracle σ_tilt(D/r₀) [px]. R0'da kareden kareye tilt değişimi
  0.09 px olduğu için zamansal eşik değiştirilmez.
- Arka-uç: **A0**, **A1** (oracle skaler kovaryans, kapılar sabit), **A0g** = A0 + `--px_noise s --stereo_sigma s
  --pnp_reproj 2.4477·s`, s = sqrt(1 + s_tot²), s_tot = hypot(1, 1.75)·σ (A1'in kovaryansı + aynı ölçekte genişlemiş
  izleme kapıları).
- Hücre başına: fe0 {A0, A1, A0g}, fe1 {A0, A1, A0g}. 72 ön-uç koşumu, 216 kol koşumu.
- Kayıt: `evidence/SONDA_KESIF/records.jsonl`; betik `sonda_kesif.py` (parmak izi kümesinin dışında).

## 3. Beklenen değerler (sonuçtan ÖNCE yazıldı)
| no | denetim | beklenen | tutmazsa |
|---|---|---|---|
| E1 | R0/fe0 taban hücrelerinde A0 ve A1 `m_ate_se3`, ana ızgaranın aynı (sekans, tohum 0, D/r₀, pencere) kaydıyla | bit-aynı (aynı kart, code_fp 08e28b83407c, aynı config) | DUR, incele |
| E2 | fe1 üçgenleme verimi (triangulated / stereo_matches), R0 | D/r₀ = 8: 0.73 → ≥ 0.88; 16: 0.47 → ≥ 0.85 | kapı uygulanmamış: DUR |
| E3 | R5 hücrelerinde enine hız | 5.000 m/s (kuru denetimde); R5/fe0 `avg_temporal_matches` D/r₀ = 8'de < 100 (R0: 112.6) | rejim uygulanmamış: DUR |
| E4 | ANEES(A0g) / ANEES(A0), R0 D/r₀ = 8 | < 0.3 ve ANEES(A1)/ANEES(A0)'ın 3 katı içinde | bayraklar uygulanmamış: DUR |
| E5 | bütün kayıtlar | çökme yok, metrik alanları dolu; ıraksama yalnız türbülanslı hücrelerde | incele |

## 4. Sondanın soruları ve karar kuralı (önceden)
- **S1 (kapı):** A0g+fe1'in A0+fe0'a göre eşli medyan ATE oranı (birim: sekans; pencereler sekans içinde log-ortalama)
  en az bir (rejim, D/r₀) hücresinde ≤ 0.85 ve 9 sekansın ≥ 7'sinde < 1 — ya da başarısız pencere sayısı yarıya
  iner. Tutarsa: çevrimiçi kestirimli (stereo diferansiyel dikey kayma = DIMM) uyarlamalı kapı kolu yazılır;
  doğrulayıcı ızgara görülmemiş tohumlarda, ON_KAYIT §5'e sapma olarak ve kullanıcı onayıyla. Tutmazsa: yöntem
  yükseltmesi bırakılır, sonuç keşif notu olarak raporlanır. Öznel ön olasılık: R0 ~%40, R5 ~%60.
- **S2 (rejim):** ATE(A0, R5) / ATE(A0, R0), D/r₀ = 8'de > 1.3 (eski hat, D/r₀ = 5: τ₀/Δt 0.1'e karşı 1 → 1.6 kat).
- **S3 (tutarlılık):** ANEES(A1, R5) < 0.5 × ANEES(A1, R0), D/r₀ = 8 (gürültü beyaza yaklaşınca skaler şişirme
  daha tutarlı).
- S2/S3 ne çıkarsa çıksın hızlı rejim ızgarası ayrı bir sapma olarak önerilecek; sonda yalnız düzey seçimine ve
  öngörülerin yazımına yardım eder.
- fe0 ile fe1 farklı ön-uç koşumlarıdır (zorunlu); A0 / A1 / A0g aynı dökümü paylaşır.

## 5. Kaynak tahmini
Ön-uç 72 × ~240 s (3 eşzamanlı) ≈ 1.6 sa; arka-uç 216 × ~450 s (12 işçi) ≈ 2.3 sa; boru hatlı toplam ≈ 2–2.5 sa.
Disk: döküm önbelleği en çok ~20 GB (hücre bitince silinir); boş 94 GB.
