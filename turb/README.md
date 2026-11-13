# `turb/` — Türbülans-Farkında Görsel-Ataletsel Odometri

Bu paket, `turbulans_farkinda_vio_proje_raporu.md` yürütme planının kod karşılığıdır.
Mevcut GTSAM/iSAM2 tabanlı stereo VIO hattına **minimal kanca** ile bağlanır:
türbülans kapalıyken hattın davranışı bit düzeyinde değişmez (rapor E0 kabul ölçütü).

## Kurulum

```bash
./turbvio/bin/python -c "import turb, gtsam, torch; print('ok')"
```

`turbvio/` bu proje için kurulmuş sanal ortamdır (py3.12, torch 2.5.1+cu124,
gtsam 4.2.1, opencv 4.9, kornia 0.7.3). Eski `vio/` ortamı bozuktur (silinebilir).

## Mimari

```
config.yaml (turbulence: bölümü)
   │
   ├── turb.hooks.injector_from_cfg ──► stereo_front_end.StereoMatcher.turb
   │        (rektifikasyon sonrası: tilt warp + kısa-poz PSF + sintilasyon
   │         + isteğe bağlı restorasyon kolu)
   │
   └── config.backend_args --turb_arm ──► gtsam_backend.turb_builder
            (A0..A5 ölçüm modeli kolları; yama-gruplu tam kovaryans)
```

### Fizik çekirdeği

| Modül | Ne yapar | Doğrulaması |
|---|---|---|
| `optics.py` | kapalı-form Kolmogorov bağıntıları | bağımsız referans |
| `screen.py` | faz ekranı (FFT + subharmonik) | `D_φ(r) = 6.88 (r/r₀)^{5/3}` ile %1 |
| `fields.py` | açıklık-ortalamalı G-tilt, kısa-poz PSF, sintilasyon | σ_α oranı 0.92 → `tilt_calibration = 1.085` |
| `correlation.py` | tilt korelasyon tensörü (boyuna/enine) + rapor §4.2 modeli | analitik katsayı **0.1706** ≈ literatür 0.170 |
| `simulator.py` | kare-tutarlı bozulma (donmuş akış) | `validate_simulator.py` |

Kritik sayısal ayrıntılar (hepsi ölçümle gerekçelendirildi, kodda yorumlu):

* **Subharmonik şart:** düz FFT ekranı r = 1 m'de yapı fonksiyonunu %58 eksik üretir.
* **Tilt spektral süzgeçle alınır:** sonlu fark + ayrık disk, iç ölçek olmadan ıraksayan
  gradyan varyansı yüzünden ızgaraya bağımlıdır. Doğru yol tam analitik transfer
  fonksiyonudur: `H_x(f) = 2πi f_x · A(f)/k`, `A(f) = 2J₁(πDf)/(πDf)`.
* **PSF birinci momenti sıfırlanır:** görüntü kayması zaten geometrik warp'tadır;
  aksi halde PSF'in artık ağırlık merkezi tüm alana ORTAK bir kayma ekler ve
  ölçülen korelasyonu bozar (ampirik: A > 1, üs ≈ 0.2).
* **Yer değiştirme alanı kaba ızgarada:** korelasyon ölçeği ≥ 37 px olduğundan
  8 px adımlı ızgara + bilinear yükseltme hatasızdır ve ~30× hızlıdır.

### Beyan edilen optik geometri

EuRoC/TUM VI iç mekân, geniş FOV kayıtlarıdır; türbülans o geometride ihmal
edilebilirdir (rapor R3). Görüntüler **beyan edilen** uzun menzilli tele-objektif
geometrisinden görülüyormuş gibi ele alınır:

```
D = 0.10 m,  λ = 550 nm,  f/p = 170648 px/rad  (IFOV 5.86 µrad/px)
FOV(752 px) = 4.41 mrad = 0.25°
σ_px = √0.170 (D/r₀)^{5/6} (λ/D)(f/p):  D/r₀=1 → 0.39 px,  5 → 1.48 px,  12 → 3.07 px
```

İddialar makalede **"denetimli bozulma benchmark'ı"** (EuRoC/TUM VI) ve
**"fiziksel olarak makul rejim"** (TartanAir dış mekân) olarak ayrılarak yazılır.

## Deneyler

```bash
./turbvio/bin/python -m turb.validate_simulator            # simülatör doğrulaması
./turbvio/bin/python -m turb.experiments --exp E1a --jobs 8 --gpu_jobs 2
./turbvio/bin/python -m turb.exp_e2 --blur none            # model doğrulaması (H-2)
./turbvio/bin/python -m turb.analyze --exp E1a             # tablo + şekil + testler
bash run_all_experiments.sh                                # tüm kuyruk
```

| Deney | Kayıt | Katkı |
|---|---|---|
| E0 | `evidence/E0/` | temel çizgi doğrulaması |
| E1a/b/c | `evidence/E1a,E1b,E1c/` | bozulma haritası (C1) |
| E2 | `evidence/E2/` | ρ̂(Δθ) model doğrulaması (C3 temeli) |
| E3 | `evidence/E3/` | restorasyon boru hattı (C2) |
| E4/E5 | `evidence/E4/` | ablasyon + tutarlılık (C3) |
| E6 | E4 kayıtları içinde (`arm_stats`) | çevrimiçi r₀ (C4) |
| TARTAN | `evidence/TARTAN/` | dış geçerlilik |

### Protokol (rapor §6.0'dan devralındı)

* 400 karelik değerlendirme penceresi, başlangıçlar 500 / 900 / 1300 / 1600
* Kabul ölçütü üçlüsü: bootstrap %95 GA sıfırı kesmeyecek **ve** Wilcoxon p < 0.05
  **ve** işaret testi aynı yön
* Hücre başına Monte-Carlo tohumları; hat belirlenimci olduğundan **belirsizliğin
  tek kaynağı türbülans gerçeklemesidir**

> **Uygulama notu (ölçümle bulundu):** adaptif ön-uçta taze stereo kareleri seyrektir
> (MH01'de kare 400–900 arasında hiç yok). Bu yüzden arka-uç pencere ortasından
> başlatılamaz — hiç landmark init edilemez ve graf boş kalır. Protokol bu nedenle
> "sekansın başından tek koşum + pencere bazlı değerlendirme" olarak uygulanır;
> bu, raporun "ortak ısınma geçmişi" şartının birebir karşılığıdır.

## Ölçüm modeli kolları (E4)

| Kol | Model | Amaç |
|---|---|---|
| A0 | izotropik köşegen R | temel çizgi |
| A1 | skaler şişirme, **oracle σ** | naif alternatif (en iyi hâliyle) |
| A2 | öznitelik-başına köşegen σ_j (çevrimiçi artık varyansı) | ara model |
| A3 | blok-korelasyonlu Σ_yama + **DIMM ile kestirilen r₀** | önerilen |
| A4 | A3 + izoplanatik-yama farkındalıklı seçim | tam sistem |
| A5 | A3 + oracle r₀/θ_c | üst sınır |

`--turb_corr_model physical|iso53` ile korelasyon modeli değiştirilebilir
(makalede ayrı bir ablasyon).

## Bulgular (kod tarafı)

1. **Rapor Ek B/5 kapandı:** 0.170 katsayısı **G-tilt**'e aittir (Z-tilt 0.182);
   sayısal integral 0.1706 verir.

2. **Rapor §4.2 modeli yetersiz** (E2, 144 hücre, `evidence/E2/e2_report.md`).
   Ampirik ρ̂ ağır kuyrukludur: 660 px ayrımda hâlâ ρ ≈ 0.25 iken iso53 modeli
   250 px'te sıfıra iner. Serbest üs uydurmasının ortalaması **0.98**, rapordaki
   5/3 = 1.67 değil.

   | D/r₀ | σ_ölçülen/σ_teorik | R² fizik (oracle θc) | R² iso53 (oracle θc) | R² fizik (θc uydurulmuş) | R² iso53 (θc uydurulmuş) |
   |---|---|---|---|---|---|
   | 3 | 1.08 | 0.51 | −4.6 | **0.82** | 0.72 |
   | 5 | 1.18 | −0.01 | −8.4 | **0.77** | 0.73 |
   | 8 | 1.20 | −1.76 | −12.9 | **0.58** | 0.50 |
   | 12 | 1.18 | −1.78 | −13.1 | **0.67** | 0.58 |

   H-2'nin ön-kayıtlı eşiği (R² ≥ 0.85) hücrelerin %35'inde fizik tensörüyle,
   %8'inde iso53 ile sağlanır. **Eşik genel olarak SAĞLANMAMIŞTIR** ve rapor §6.0
   gereği gevşetilmez; makalede bu haliyle raporlanır.

3. **Ölçülen tilt genliği teoriyle uyumlu:** σ̂/σ_teorik = 1.08–1.20 (dedektör
   gürültüsü "nugget" düzeltmesiyle ayrıştırıldı, korelasyonlu pay A ≈ 0.92–1.10).
   Ön-uç GT eğim alanını 0.36–0.44 px hatayla ölçüyor; tekrarlanabilirlik
   D/r₀ = 3 → 12 arasında 0.94 → 0.92.
