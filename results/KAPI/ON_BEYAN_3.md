# Uyarlanır izleme kapısı kolu — ön testler ve geliştirme koşumu: ön beyan (KEŞİF; ön kayıt dışı)

Yazım: 2026-11-11 19:30 (hiçbir T1–T3 / geliştirme sonucu görülmeden). Ağaç: `pyeco_turb_kapi` (geliştirme kopyası;
kod parmak izi `bedca2849bb2`; koşan ana ağaç `pyeco_turb` `08e28b83407c` değişmedi). Koşucu: `kapi_kol.py`.
Veri: saklı ön-uç dökümleri, **tohum 0**, R0 (yarı-statik), D/r₀ 8, 9 sekans (`scratchpad/a3s_dokum_s0/`).

**Statü.** Tohum 0 hücreleri ilk sondada görüldü; buradaki hiçbir sayı doğrulayıcı kanıt değildir, makalede yalnız
"keşif / geliştirme" etiketiyle anılabilir. Bu koşumun tek amacı (i) yeni kodun doğruluğunu kanıtlamak, (ii) doğrulayıcı
ızgarada (ON_KAYIT §5 Sapma 6; görülmemiş tohumlar 1–4) dondurulacak tek yapılandırmayı seçmektir.

## 1. Değişiklik
`--turb_gate` (varsayılan 0 = kapalı): izleme cephesinin sabit ölçekli kapıları, kolun **son graf karesinde**
kullandığı bireysel-faktör sigmasıyla ölçeklenir (nedensel; kestirici ek çağrı almaz).
bit 1: motion-only BA sigması + χ² aykırı eleme + graf kapısı + yerel harita χ² süzgeci; bit 2: PnP-RANSAC eşiği
(`pnp_reproj · σ_mono / px_noise`). Sabit piksel pencereleri (`moba_pencere_px`, harita arama yarıçapı) ölçeklenmez.
Çevrimiçi kol = `A0C+A3`, `--turb_cal_factor 1` (A3S: σ² = σ₀² + (c·σ̂_DIMM)², skaler) + `--turb_gate 3`.
Adaylar: **D1** c = 1 (σ̂ toplam eksen sigması sayılır; ek sabit yok), **D2** c = hypot(1, 1.75) = 2.0156 (A1'in çarpanı).
Oracle üst sınır: **A1g** = A1 + `--turb_gate 3` (ilk sondanın A0g'siyle aynı ölçek, s = 4.5241510).

## 2. Ön testler (V1_01_easy; kapı tutmazsa geliştirme koşumu başlamaz)
**T1 gerileme (kapı kapalı).** Yeni kodda A0, A0C (k = 12), A1, A3S (A0C+A3, k = 1), A3: `m_ate_se3`,
`n_nees_pos_mean`, `n_nees_rot_mean`, `m_rpe_trans` × 3 pencere, eski kodun kayıtlarıyla (ana ızgara / A3S ızgarası,
tohum 0) **bit-aynı**. Beklenen: 15 / 15. (A3 yavaş; ayrı aşamada koşar, doğrulayıcı ızgaradan önce tamamlanır.)

**T2 eşdeğerlik.** (a) A0C(k = 4.524151) + gate 3 ≡ A0 + `--px_noise 4.524151 --stereo_sigma 4.524151 --pnp_reproj
11.0737644027`; (b) A0C(k = 1) + gate 3 ≡ A0C(k = 1). Beklenen: metrik 3 / 3 bit-aynı **ve** iki yörünge dosyasının
(pürüzsüz + gerçek zamanlı) sha256'sı aynı. Bilgi (kapı değil): A0C(k = 1) ile A0; yeni A0g ile ilk sondanın A0g kaydı
(PnP eşiği 1e-6 yuvarlama farkı var; bit-aynı olması beklenmez, G3 duyarlılığı).

**T3 duman** (D1, D2, D1n = gate 1, D1p = gate 2, A1g). Sert:
- H1: 3 / 3 pencere `status` ok, `m_valid`, `n_valid`, `m_n_compare` 400 ± 20, dört metrik sonlu;
- H2: `arm_stats.gate` beyan edilen değer; `gate_sigma_n` = 1299 (= 1300 kare − 1; kare 0'da izleme yok);
- H3: çevrimiçi kollarda `gate_sigma_min` = √(1 + (c·0.05)²) (ilk kestirimden önce taban: 1.0012492 / 1.0050653) ve
  `gate_sigma_max` > 1.3; A1g'de min = maks = 4.5241510 (± 1e-9);
- H4: etki uygulandı — üç pencerede ATE kapısız karşılığından (A3S / A1 kaydı) farklı.

Yumuşak (yalnız raporlanır; tek sekans, G3 duyarlılığı):
- S1: kapı sigması medyanı D1 ∈ [1.6, 4.0] px; D2 ∈ [2.8, 9.0] px. Gerekçe: kapısız A3S'te σ̂ medyanı 1.28–1.62 px
  (s ≥ 1.6 / 2.8); kapı genişleyince kesilme yanlılığı azalır, σ̂ yükselir; oracle toplam 4.41 px.
- S2: D1n, D1p ∈ [1.5, 4.0] px.
- S3: V1_01_easy'de ATE (3 pencere ortalaması) / A0: D1, D2, D1n ∈ [0.3, 1.0]; A1g ∈ [0.3, 0.8]; D1p ∈ [0.7, 1.3].

## 3. Geliştirme koşumu (9 sekans × {D1, D2}; birim = sekans; test yok)
Ölçü: ATE pencerelerde aritmetik, konum ANEES'i geometrik ortalama; etki = sekanslarda eşli medyan oran
(ilk sondayla aynı kural). Karşılaştırma kayıtları eski koddan (T1 bit-aynılığı gerekçesiyle): A0, A1 (ana ızgara),
A3S (A3S ızgarası), A0g (ilk sonda).

Beklenen değerler (sonuç görülmeden):
| no | ölçü | beklenen |
|---|---|---|
| G1 | D2 / A0 ATE | 0.50–0.85, ≥ 7 / 9 sekansta < 1 |
| G2 | D1 / A0 ATE | 0.60–0.95 (kapı daha dar; kazanç daha küçük) |
| G3 | D / A3S ATE (aynı kovaryans modeli, kapı farkı) | < 1 (D2: 0.5–0.9) |
| G4 | D / A0g ATE (oracle'a göre bedel) | 1.0–1.6 |
| G5 | konum ANEES'i D / A0 | < 0.5; D2 < D1 |
| G6 | ıraksayan / başarısız pencere | her aday ≤ 1 / 27 |
| G7 | kapı sigması medyanı (sekans medyanı) | D1 1.6–4.0 px; D2 2.8–9.0 px; D2'de 9 px üstü = kaçak geri besleme işareti |

**Seçim kuralı (dondurulacak yapılandırma).**
"Geçer" = aday / A0 eşli medyan ATE oranı ≤ 0.85 **ve** 9 sekansın ≥ 7'sinde < 1 **ve** başarısız pencere ≤ 1 / 27.
1. Yalnız biri geçerse o seçilir.
2. İkisi de geçerse medyan oranı küçük olan seçilir; iki medyan oranın oranı [0.90, 1.10] içindeyse **D1** seçilir
   (ek sabit yok, öz-tutarlı).
3. Hiçbiri geçmezse çevrimiçi kol bu hâliyle doğrulayıcı ızgaraya götürülmez; sonuç olduğu gibi raporlanır. Yeni aday
   (ör. başka c, kesilme düzeltmeli σ̂) yalnız tohum 0'da ve bu dosyaya koşumdan **önce** eklenerek denenebilir.
Tutarlılık (ANEES) seçim ölçütü değildir; ikincil sonuç olarak raporlanır.
Doğrulayıcı ızgara tek, dondurulmuş yapılandırmayla, bir kez, görülmemiş tohumlarda koşulur (Sapma 6'da ön kayıt).

## 4. Sınırlar
- Tek tohum, 9 birim; güven aralığı ve test yok. 0.9–1.1 arası oranlar sayısal duyarlılıktan (G3) ayırt edilemez.
- A0g / A1g oracle σ_tilt kullanır (gerçek sistemde bilinmez); çevrimiçi kollar yalnız kendi artıklarını kullanır.
- Keşif sondası 2 (kapı ayırma; aynı dökümler) bu yazım anında sürüyor; oran tabloları görülmedi.
