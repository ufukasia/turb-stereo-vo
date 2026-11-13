# Keşif sondası 2 — kapı ayırma ve kestirimli kapı: ön beyan (2026-11-11 18:53:00)

**Statü:** KEŞİF. Ön kayıtlı (ON_KAYIT_FAZ2.md) doğrulayıcı analizlerin parçası DEĞİL. Yalnız tohum 0; bu hücrelerin
taban ve geniş kapı sonuçları zaten görüldü (`SONUC.md`). Burada yazılan hiçbir sayı makalede doğrulayıcı kanıt
olamaz; amaç, kullanıcının "çevrimiçi uyarlanır kapı kolu yazılsın mı" kararına iki bilgi eklemek. Doğrulayıcı kapı
ızgarası bu sondanın parçası değildir ve kullanıcı onayı olmadan koşulmaz. Bu belge sonda 2'nin hiçbir sonucu
görülmeden yazıldı; özeti `ON_BEYAN_2.sha256` dosyasında.

**Bilinenler:** Faz 2'nin bütün ön kayıtlı sonuçları; keşif sondasının bütün sonuçları (`SONUC.md`: R0, D/r₀ 8, fe0'da
A0g / A0 ATE oranı 0.517, 9 / 9; A1 / A0 1.214; A0g / A0 konum ANEES'i 0.044); ana ızgaranın tohum 0 A3 kayıtlarındaki
DIMM σ̂ medyanları.

## 1. Sorular
- **Q1 (kestirim):** A0g'nin kazancı oracle σ_tilt ile ölçüldü (gerçek sistemde bilinmez). Kapı ölçeği, sistemin
  kendi çevrimiçi kestirimi düzeyine (DIMM σ̂; D/r₀ 8'de oracle değerin 0.59–0.74 katı) indirilirse kazanç korunur mu?
- **Q2 (ayırma):** A0g üç bayrağı birlikte değiştiriyordu. Kazancı hangisi taşıyor: gürültü grubu (`--px_noise`,
  `--stereo_sigma`: düzleştirici gürültü modeli + bütün χ² kapıları, çünkü χ² eşikleri σ₀² ile ölçeklenir) mi,
  PnP RANSAC eşiği (`--pnp_reproj`) mi?

## 2. Tasarım (parmak izli kodda değişiklik YOK; `kapi_sonda.py` parmak izi kümesinin dışında)
- Koşul: R0 (τ₀/Δt = 1), D/r₀ = 8, tohum 0, 9 sekans, kapsam1m pencereleri (3 pencere), fe0. Kol: A0;
  `--turb_cal_factor 12.0` (sondayla aynı; A0'da etkisiz).
- Ön-uç çıktısı: A3S ızgarasının sakladığı tohum 0 dökümleri (`scratchpad/a3s_dokum_s0/{SEKANS}__dr8__s0.pkl`; aynı
  ön-uç ayarı; sondanın ön-uç çıktısıyla aynı olması beklenir, B1 bunu sınar). Ekran kartı kullanılmaz.
- s(σ) = sqrt(1 + (hypot(1, 1.75)·σ)²) (ilk sondayla aynı formül). Oracle: σ = 2.1891 px → s = 4.524151.
- Varyantlar:

| ad | bayraklar | sekans |
|---|---|---|
| K (kontrol) | `--px_noise 4.524151 --stereo_sigma 4.524151 --pnp_reproj 11.073765` (ilk sondanın A0g hücresiyle aynı) | V1_01_easy, V2_01_easy |
| Gd (kestirimli) | üç bayrak, s = s(σ̂); σ̂ = ana ızgara, tohum 0, D/r₀ 8, aynı sekansın A3 kaydındaki `sigma_hat_median` (1.284–1.623 px → s = 2.775–3.420) | 9 |
| Gn (yalnız gürültü grubu) | `--px_noise 4.524151 --stereo_sigma 4.524151`; PnP eşiği varsayılan (2.4477) | 9 |
| Gp (yalnız PnP) | `--pnp_reproj 11.073765`; gürültü varsayılan (1.0) | 9 |

- Başvurular (yeniden koşulmaz): ilk sondanın R0 / D/r₀ 8 / fe0 kayıtları (A0, A1, A0g).
- Birim = sekans; ATE pencerelerde aritmetik, konum ANEES'i geometrik ortalama; etki = sekanslarda eşli medyan oran
  ve oranı < 1 olan sekans sayısı (`sonda_analiz.py` ile aynı kural). Test ve güven aralığı yok.
- Koşum sırası: önce ön-test kümesi (K × 2 sekans ve V1_01_easy'de Gd, Gn, Gp); kapılar tutarsa kalanlar.

## 3. Beklenen değerler (sonuçlardan ÖNCE yazıldı)
- **B1 (regresyon kapısı):** K kayıtları ilk sondanın A0g kayıtlarıyla bit-aynı (6 / 6 kayıt; `m_ate_se3`,
  `n_nees_pos_mean`, `n_nees_rot_mean`, `m_rpe_trans`). Tutmazsa DUR (döküm ya da zincir farklı); kalanlar koşulmaz.
- **B2 (duman):** ön-test kümesindeki bütün kayıtlar `status = ok`, `m_valid` ve `n_valid` dolu, karşılaştırılan poz
  sayısı 400 ± %5; Gd, Gn, Gp kayıtlarının `m_ate_se3` değeri aynı penceredeki A0 kaydından farklı (bayrak uygulandı).
  Tutmazsa DUR.
- **B3 (Q1):** Gd / A0 ATE medyan oranı **0.50–0.80**, 9 sekansın ≥ 7'sinde < 1. Karar kuralı (ilk sondanın S1
  kuralıyla aynı eşik): oran ≤ 0.85 ve ≥ 7 / 9 sekansta < 1 ise "kestirim düzeyindeki ölçekle kazanç korunuyor";
  değilse "kazanç oracle ölçeğe bağlı". Gd / A0g: 1.0–1.5 (kestirim oracle'dan dar; biraz kayıp beklenir).
- **B4 (Q2):** Gn / A0 ATE **0.55–0.90**; Gp / A0 ATE **0.85–1.15**. Bir grup "kazancı taşıyor" sayılır: oran ≤ 0.85
  ve ≥ 7 / 9 sekansta < 1. Öznel olasılıklar: yalnız gürültü grubu taşır ≈ %55; yalnız PnP ≈ %10; ikisi ayrı ayrı
  taşır ≈ %10; hiçbiri tek başına taşımaz (birlikte gerekir) ≈ %25.
- **B5 (tutarlılık, etkinin uygulandığı):** konum ANEES'i Gn / A0 < 0.3 ve Gd / A0 < 0.3 (gürültü modeli şişti);
  Gp / A0 0.5–1.5 (gürültü modeli değişmedi).
- **B6:** çökme yok; ıraksayan pencere varyant başına ≤ 1 / 27 (ilk sondada R0 D/r₀ 8'de hiç yoktu).

## 4. Sonuçların kullanımı ve sınırlar
- Q1 tutarsa ve kullanıcı onaylarsa: çevrimiçi kol tasarımında kapı ölçeği DIMM σ̂'dan alınır; Q2, hangi eşiklerin
  uyarlanacağını belirler. İkisi de tutmazsa kullanıcıya "kazanç oracle bilgiye bağlı" diye raporlanır.
- Gd'de σ̂ aynı sekansın önceki bir koşumunun medyanıdır (nedensel değil, sabit); bu, kestirimin **düzeyini** sınar,
  çevrimiçi bir uygulamayı değil. Çevrimiçi kol parmak izli kodda değişiklik ister; koşumlar bittikten sonra,
  regresyon kapısıyla ve yalnız kullanıcı onayıyla yazılır.
- Tek tohum, 9 birim; arka-uç sayısal olarak duyarlı (Sapma 3 eki): 0.9–1.2 arasındaki oranlar gürültüden ayırt
  edilemez.
