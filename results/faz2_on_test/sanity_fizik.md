# Tahlil raporu — fizik
Üretim: 2026-11-07 23:34 · genel durum: **SARI**

## S — Simülatör fiziği (yazılım doğruluğu)

### ⚠️ `S7b-kisa-poz-fwhm` — Aynı karşılaştırma FWHM ile (bilgi amaçlı; çekirdek+hale yapısında kararsız)
- **Ölçülen:** D/r0: FWHM simülasyon/teori = 1:0.97; 3:0.94; 5:0.83; 8:0.58; 12:0.58
- **Beklenen:** 1.00 ± 0.20
- **Kaynak:** Fried 1966 JOSA 56:1372
- **Yapılacak:** FWHM oranı [8.0, 12.0] düzeylerinde bandın dışında: simülasyonun ortalama kısa-poz PSF çekirdeği Fried'in kısa-poz yaklaşımından dar (EE50 uyumlu).

### ✅ `S1-tilt-varyansi` — Hattın ürettiği açıklık-ortalamalı G-tilt varyansı kapalı-form değere uyuyor mu?
- **Ölçülen:** sigma_hat/sigma_teori = 0.950  (D/r0=5.0, N=1024, ppd=8, kalibrasyon=1.085, n=4)
- **Beklenen:** 1.00 ± 0.10;  sigma_a^2 = 0.17 (D/r0)^(5/3) (lam/D)^2
- **Kaynak:** Tyler 1994 JOSA A 11:358; Lofdahl 2010 (K_G=0.340/2)

### ✅ `S1b-kalibrasyon-mesruiyeti` — Eğim kalibrasyon sabiti, ölçülen ayrıklaştırma açığıyla açıklanıyor mu?
- **Ölçülen:** ham oran = 0.891;  1/ham = 1.122;  kullanılan kalibrasyon = 1.085;  çarpım = 0.967
- **Beklenen:** kalibrasyon in [0.80, 1.25] ve ham x kalibrasyon = 1.00 ± 0.12
- **Kaynak:** Schmidt 2010 'Numerical Simulation of Optical Wave Propagation' §9 (ayrık faz ekranlarında düşük/yüksek frekans açığı ve düzeltmesi)

### ✅ `S2-olcekleme-usteli` — Eğim varyansı (D/r0)^{5/3} ile ölçekleniyor mu?
- **Ölçülen:** ölçülen üstel = 1.736  (D/r0 = 2, 4, 8)
- **Beklenen:** 5/3 = 1.667 ± 0.20
- **Kaynak:** Kolmogorov spektrumu, f^{-11/3}

### ✅ `S3-yapi-fonksiyonu` — Faz yapı fonksiyonu D_phi(r) = 6.88 (r/r0)^{5/3} — ekran güç açığı
- **Ölçülen:** ince ızgara (D/32): 0.957 · kaba ızgara (D/8): 0.860  (0.15 r0 < r < 1.5 r0, ekran genişliği 64 D)
- **Beklenen:** 1.00 ± 0.10 (ayrık ekranlarda %10'a kadar açık literatürde beklenir; > %25 ise normalizasyon hatası)
- **Kaynak:** Fried 1965; Roddier 1981 Prog. Opt. 19; Schmidt 2010 §9 (ayrık ekran güç açığı)

### ✅ `S4-korelasyon-tensoru` — Eğim korelasyon tensörü BAĞIMSIZ DIMM kapalı-formuyla uyuşuyor mu?
- **Ölçülen:** s/D=1:(0.589|0.814), 2:(0.427|0.653), 4:(0.321|0.509), 6:(0.272|0.437), 8:(0.241|0.392)  — en büyük bağıl sapma 11.0%; rho(0)=1 True; tekdüze True
- **Beklenen:** rho_boy = 0.541 (D/s)^{1/3}, rho_en = 0.810 (D/s)^{1/3}, sapma < %15
- **Kaynak:** Sarazin & Roddier 1990 A&A 227:294 (DIMM); Fried 1975 JOSA 65:1427

### ✅ `S5-anizotropi-orani` — Eğim korelasyonunun anizotropisi (boyuna/enine) doğru mu? Bu oran, köşegen-olmayan Sigma'nın VAROLUŞ GEREKÇESİDİR
- **Ölçülen:** rho_boyuna/rho_enine medyanı = 0.631
- **Beklenen:** 2/3 = 0.667 ± 0.10 (Kolmogorov, s >> D)
- **Kaynak:** Sarazin & Roddier 1990: 0.0968/0.1450 = 2/3

### ✅ `S6-kirinim-fwhm` — Kırınım sınırlı çekirdeğin FWHM'i görüntü pikselinde lam/D*f/p mi?
- **Ölçülen:** FWHM = 0.949 px; lam/D*f/p = 0.939 px; oran = 1.011 (Airy 1.029); çekirdek toplamı 1 = True; psf_legacy_scale=False; eski kip oranı = 4.78
- **Beklenen:** oran 1.00 ± 0.15 ve çekirdek toplamı 1
- **Kaynak:** Born & Wolf, Principles of Optics §8.5.2 (Airy FWHM = 1.029 lam/D)

### ✅ `S7-kisa-poz-olcekleme` — Tilt çıkarılmış kısa-poz PSF genişliği (EE50 çapı) D/r0 ile Fried kısa-poz bağıntısına göre ölçekleniyor mu?
- **Ölçülen:** D/r0: EE50 simülasyon/teori [px] = 1:1.12/1.14=0.98; 3:2.37/2.36=1.01; 5:3.91/3.88=1.01; 8:6.33/6.37=0.99; 12:9.45/9.86=0.96
- **Beklenen:** her D/r0'da oran 1.00 ± 0.20 (alpha = 1, yakın alan)
- **Kaynak:** Fried 1966 JOSA 56:1372 (kısa-poz OTF); Goodman Statistical Optics §8.7
