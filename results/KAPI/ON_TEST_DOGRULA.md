# Sapma 6 doğrulayıcı ızgarası — koşucu ön testi (tohum 0, V1_01_easy): ön beyan ve sonuç

**§1–§3 yazım: 2026-11-11 20:00** — geliştirme koşumu sürerken (16 koşumun 1'i bitmişti; oran tablosu görülmedi),
doğrulayıcı koşucunun hiçbir çıktısı üretilmeden. §4 koşumdan sonra eklenir (yalnız ekleme).

**Statü.** Tohum 0, tek sekans: doğrulayıcı kanıt değildir (tohum 0 hücreleri sondalarda ve geliştirme koşumunda
görüldü). Amaç yalnız koşucu + analiz zincirinin doğruluğunu ve süresini uzun koşumdan önce kanıtlamaktır.

## 1. Kimlik
- Ağaç `pyeco_turb_kapi`, kod parmak izi `bedca2849bb2` (ana ağaç `08e28b83407c`; koşan ağaçta değişiklik yok).
- Koşucu `kapi_dogrula.py` sha256 `55a91eb69839ea54…`; analiz `kapi_dogrula_analiz.py` sha256 `c7bd583a55d720a2…`
  (öz-testi geçti: sentetik veri + ana ızgaranın A1 / A0 aile değerleri 4 / 4 sayı sayı aynı).
- Yapılandırma: `evidence/KAPI/SECIM.json` (geliştirme koşumu bitince ON_BEYAN_3.md kuralıyla dondurulur; ön test
  ondan sonra başlar). Aşağıda aday D1 ve D2 için beklenenler ayrı ayrı yazıldı.
- Kuru koşum (19:4x, `--kuru --aday D1`): 184 hücre, 571 arka-uç koşumu (1713 kayıt; ön test dahil); hücre
  tanımları ana ızgara (parça II, III) ve hızlı rejim ızgarasıyla (parça I) kol / bayrak / etiket dışında aynı.

## 2. Plan (4 hücre, 13 arka-uç koşumu, 39 kayıt; 3 ön-uç yeniden üretimi)
| parça | koşul | döküm | kollar |
|---|---|---|---|
| II | R0, D/r₀ 8 | ön-uç yeniden üretilir | A0, KAPI, A1g |
| III | R0, D/r₀ 0 (temiz) | ön-uç yeniden üretilir | A0, A3S, KAPI |
| I | R5, D/r₀ 8 | ön-uç yeniden üretilir | A0, KAPI, A3S, A1g |
| II | R0, D/r₀ 16 | saklı döküm (`a3s_dokum_s0`) | A0, KAPI, A1g |

## 3. Beklenen değerler (sonuçtan önce)
**Sert kapılar (biri tutmazsa ızgara başlamaz):**
- **K1 gerileme / bit-aynılık:** yeni kodun A0 kayıtları eski kodun kayıtlarıyla dört alanda (`m_ate_se3`,
  `n_nees_pos_mean`, `n_nees_rot_mean`, `m_rpe_trans`) bit-aynı: R0-8, R0-16, R0-0 → ana ızgara tohum 0 (3 / 3'er);
  R5-8 → ilk sonda `sonda_R5_fe0` (3 / 3). Toplam 12 / 12.
- **K2 ön-uç:** yeniden üretilen üç ön-ucun istatistikleri (süre alanı dışında) başvuru A0 kaydınınkiyle aynı (3 / 3).
- **K3 iki döküm yolunun eşdeğerliği:** R0-8'de KAPI ve A1g (yeniden üretilen ön-uç) geliştirme / ilk ön test
  kayıtlarıyla (saklı döküm) bit-aynı: metrik 3 / 3 + pürüzsüz ve gerçek-zamanlı yörünge dosyalarının sha256'sı aynı.
- **K4 alanlar ve kapı istatistiği:** 39 / 39 kayıt `status` ok, `m_valid`, `n_valid`, `m_n_compare` 400 ± 20, dört
  metrik sonlu; kapılı kollarda `gate` = 3, `gate_sigma_n` = 1299; A1g kapı sigması min = maks = 4.5241510 (D/r₀ 8),
  7.9250780 (D/r₀ 16) (± 1e-9); KAPI en küçük kapı sigması = taban (D1 1.0012492, D2 1.0050653).
- **K4b analiz okuması:** `kapi_dogrula_analiz.py --ontest --butunluk` → tam = True; yeni kayıt 39 / 39; başvuru
  kayıtları R0|A0 9, R0|A1 6, R0|A3S 6, R5|A0 3, R5|A1 3; V1 (A0y ≡ A0) 12 / 12 aynı; V2 9 / 9 aynı (3 hücre × 3
  pencere); V3 (kapı istatistiği) 21 / 21; V4 (kapısız kollarda kapı kapalı) 18 / 18; V5 (bayraklar, aday) 39 / 39.

**Yumuşak (yalnız raporlanır; tek sekans, sayısal duyarlılık ±%10–20):** ATE = 3 pencere ortalaması, oran = kol / A0.
| koşul | ölçü | beklenen | dayanak |
|---|---|---|---|
| R0-8 | KAPI / A0 | D1 0.896, D2 0.749 (K3 gereği tam bu) | ilk ön test |
| R0-8 | A1g / A0 | 0.798 (K3 gereği) | ilk ön test |
| R0-0 | A3S / A0 | 0.90–1.10 | H5: kovaryans doğruluğu değiştirmiyor |
| R0-0 | KAPI / A0 | 0.88–1.12 | temizde σ̂ medyanı bu sekansta 0.57 px → kapı ≈ 1.15 (D1), 1.5 (D2) |
| R0-0 | KAPI kapı sigması medyanı | D1 1.10–1.45; D2 1.4–2.6 px | aynı |
| R5-8 | A1g / A0 | 0.30–0.60 | ilk sonda: A0g / A0 bu sekansta 0.415 |
| R5-8 | KAPI / A0 | 0.35–0.90 | kestirim oracle'ın altında (σ̂ ≈ 1.5 px; R5 ön testi) |
| R5-8 | A3S / A0 | 0.80–1.20 | kapısız; ilk sonda A1 / A0 0.84 |
| R5-8 | KAPI kapı sigması medyanı | D1 1.6–4.0; D2 2.8–9.0 px | R0-8 ile aynı σ̂ düzeyi |
| R0-16 | A1g / A0 | 0.35–0.80 | ilk sonda: A0g / A0 bu sekansta 0.53 |
| R0-16 | KAPI / A0 | 0.40–1.00 | σ̂ ≈ 1.74 px (kapısız), oracle toplam 7.86 px: kapı oracle'ın çok altında |
| R0-16 | KAPI kapı sigması medyanı | D1 1.9–5.0; D2 3.5–12 px | |
| hepsi | etki uygulandı | KAPI ATE'si aynı hücredeki A3S / A0 kaydından üç pencerede farklı | |
| süre | arka-uç koşumu | 250–700 s (yük altında); ön-uç 150–400 s | ilk ön test, geliştirme koşumu |

## 4. Sonuç
**Yazım: 2026-11-11 21:16.** Koşum 20:47:25 – 21:14:29 (27.0 dk; 5 arka-uç + 1 ön-uç işçisi; makine yükü ≈ 65–72), çıkış kodu 0.
Yapılandırma: `SECIM.json` sha256 `cd4f4520183e0ee4…` (aday **D1**, c = 1.0, kapı 3; 20:46'da donduruldu). Günlük
`dogrula_ontest_2047.log`; kayıtlar `dogrula_ontest.jsonl` (39 satır; sha256 `26b585db15616917…`). Koşucu ve analiz
betiği §1'deki özetlerle aynı (değişmedi).

**Sert kapılar: hepsi GEÇTİ → "[kapı] Sapma 6 ön testi (K1-K4): GEÇTİ".**
| kapı | beklenen | ölçülen | sonuç |
|---|---|---|---|
| K1 A0y ≡ eski A0 (4 alan) | 12 / 12 | 12 / 12 (R0-8, R0-16, R0-0 ana ızgara; R5-8 ilk sonda) | geçti |
| K2 ön-uç istatistikleri | 3 / 3 | 3 / 3 ("istatistikler başvuruyla aynı") | geçti |
| K3 iki döküm yolu (R0-8 KAPI, A1g) | metrik 3 / 3 + yörünge sha256 | KAPI 3 / 3, A1g 3 / 3; pürüzsüz ve gerçek-zaman yörüngeleri aynı | geçti |
| K4 alanlar | 39 / 39 ok; `m_n_compare` 400 ± 20 | 39 / 39; hepsi 400 | geçti |
| K4 kapı istatistiği | `gate` 3; n 1299; A1g 4.5241510 / 7.9250780; KAPI min 1.0012492 | 3; 1299; 4.5241510413 / 7.9250782499; 1.0012492197 | geçti |
| K4b analiz okuması (`--ontest --butunluk`) | tam = True; yeni 39 / 39; başvuru 9 / 6 / 6 / 3 / 3; V1 12, V2 9, V3 21, V4 18, V5 39 | tam = True; 39 / 39; R0\|A0 9, R0\|A1 6, R0\|A3S 6, R5\|A0 3, R5\|A1 3; V1 12, V2 9, V3 21, V4 18, V5 39 | geçti |

**Yumuşak beklentiler (tek sekans, tohum 0; yalnız rapor):**
| koşul | ölçü | beklenen | ölçülen | sonuç |
|---|---|---|---|---|
| R0-8 | KAPI / A0 | 0.896 (D1) | 0.896 | tuttu (K3) |
| R0-8 | A1g / A0 | 0.798 | 0.798 | tuttu (K3) |
| R0-0 | A3S / A0 | 0.90–1.10 | 0.997 | tuttu |
| R0-0 | KAPI / A0 | 0.88–1.12 | 1.012 | tuttu |
| R0-0 | kapı sigması medyanı | 1.10–1.45 px | 1.154 (maks 1.317; σ̂ medyanı 0.577) | tuttu |
| R5-8 | A1g / A0 | 0.30–0.60 | 0.542 | tuttu |
| R5-8 | KAPI / A0 | 0.35–0.90 | 0.437 | tuttu |
| R5-8 | A3S / A0 | 0.80–1.20 | 0.779 | **tutmadı** (bandın 0.02 altında) |
| R5-8 | kapı sigması medyanı | 1.6–4.0 px | 1.980 (maks 2.252; σ̂ 1.849) | tuttu |
| R0-16 | A1g / A0 | 0.35–0.80 | 0.530 | tuttu |
| R0-16 | KAPI / A0 | 0.40–1.00 | 0.660 | tuttu |
| R0-16 | kapı sigması medyanı | 1.9–5.0 px | 2.235 (maks 2.418; σ̂ 1.980) | tuttu |
| hepsi | etki uygulandı | KAPI üç pencerede A0'dan / A3S'ten farklı | 3 / 3 (dört hücre); A3S'ten 3 / 3 (iki hücre) | tuttu |
| süre | arka-uç 250–700 s; ön-uç 150–400 s | | arka-uç 325–919 s (R0-16: 566 / 762 / 919 s bandın üstünde); ön-uç 156 / 171 s, temiz 64 s (bandın altında: türbülans benzetimi yok) | **kısmen tutmadı** |

Notlar (keşif; tohum 0, tek sekans):
- R5-8'de KAPI (0.437) oracle kapılı A1g'den (0.542) daha düşük ATE verdi; tek hücre, sayısal duyarlılık içinde
  yorumlanmaz. R0-8'de bu sekansta KAPI'nin konum ANEES'i A0'dan yüksek (559 / 414): geliştirme koşumunda görülen tek
  ters sekans (oran 1.35) budur.
- Ölçülen süreler (V1_01_easy, 1300 kare): R0-8 331 / 415 / 445 s (A0 / KAPI / A1g); temiz 325 / 334 / 333 s (A0 / A3S /
  KAPI); R5-8 465 / 586 / 440 / 676 s (A0 / KAPI / A3S / A1g); R0-16 566 / 762 / 919 s. Dokuz sekansın ortalama kare
  sayısı 1434 (× 1.10). Izgara tahmini: parça I ≈ 40, parça II ≈ 38, parça III ≈ 11 işçi-saat (toplam ≈ 89), bu yük
  düzeyinde; ön-uç yeniden üretimi 36 × ≈ 170 s (R0-8) + 36 × ≈ 70 s (temiz) GPU.

**Karar:** sert kapıların hepsi geçti; doğrulayıcı ızgara, Sapma 6 ON_KAYIT §5'e eklendikten sonra başlatılır.
