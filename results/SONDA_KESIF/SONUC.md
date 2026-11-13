# Keşif sondası — sonuç (KEŞİF; ön kayıt dışı)

Yazım: 2026-11-11 18:42. Koşum 12:38–18:33 (354.4 dk), 144 / 144 hücre, 648 kayıt.
Kaynak: `records.jsonl` (sha256 `9e8b920a…`), `python sonda_analiz.py records.jsonl`; ön beyan `ON_BEYAN.md`
(sha256 `ffad5722…`, 12:24:24, koşumdan önce). Kod parmak izi bütün kayıtlarda `08e28b83407c`.

**Statü:** yalnız tohum 0, 9 sekans, birim = sekans. Hiçbir hipotez testi yapılmadı; sayılar eşli medyan oran ve
"oranı < 1 olan sekans sayısı"dır. Geniş kapı kolu (A0g) ve fe1 **oracle** σ_tilt kullanır (gerçek sistemde bilinmez).
Bu sonuçlar makalede doğrulayıcı kanıt olarak kullanılamaz; yalnız sonraki ızgaraların tasarımına girer.

## 1. Ön beyandaki beklenen değerler (E1–E5)

| no | beklenen | ölçülen | sonuç |
|---|---|---|---|
| E1 | R0/fe0 taban hücreleri (A0, A1) ana ızgara tohum 0 kaydıyla bit-aynı | 108 bit-aynı, 0 farklı | tuttu |
| E2 | fe1 üçgenleme verimi R0: D/r₀ 8 ≥ 0.88; 16 ≥ 0.85 | 0.732 → 0.973; 0.473 → 0.947 | tuttu |
| E3 | R5/fe0 zamansal eşleşme, D/r₀ 8: < 100 (R0: 112.6) | 101.9 (R0: 113.2; −%10) | **sayısal eşik tutmadı** (yön beklenen) |
| E4 | ANEES(A0g) / ANEES(A0), R0 D/r₀ 8: < 0.3 ve A1 / A0'ın 3 katı içinde | 0.044 (A1 / A0: 0.101) | tuttu |
| E5 | çökme yok, metrik alanları dolu; ıraksama yalnız türbülanslı hücrelerde | 648 / 648 `status` ok, `m_valid` ve `n_valid` 648 / 648; ıraksayan 8 pencere (hepsi türbülanslı) | tuttu |

E3 notu: eşik 1.9 eşleşme farkla tutmadı. Kuralın amacı rejimin uygulandığını doğrulamaktı; bu ayrıca kayıtlardaki
`tau0_over_dt` alanıyla (R5: 0.0157 / 0.0314; R0: 1.0) ve R5 ile R0 sonuçlarının 9 / 9 sekansta farklı olmasıyla
doğrulanıyor. Eşik yanlış kalibre edilmişti; değiştirilmedi, tutmadığı raporlanıyor.

Iraksayan pencereler (8): R0 D/r₀ 16: A1/fe0 MH03-562, MH04-532; A0g/fe0 MH03-562, MH04-532; A1/fe1 MH04-532;
A0g/fe1 MH03-562. R5: A1/fe1 D/r₀ 4 V1_02_medium-100; A0/fe1 D/r₀ 8 V1_02_medium-100. A0/fe0 hiçbir hücrede ıraksamadı.

## 2. Sonda soruları (S1–S3)

### S1 (kapı) — TUTTU
Kural: (A0g, fe1) / (A0, fe0) eşli medyan ATE oranı en az bir (rejim, D/r₀) hücresinde ≤ 0.85 ve 9 sekansın ≥ 7'sinde < 1.

| rejim, D/r₀ | (A0g,fe1) / (A0,fe0) | A0g / A0 [fe0] | A0g / A1 [fe0] | A1 / A0 [fe0] | fe1 / fe0 [A0] | fe1 / fe0 [A0g] |
|---|---|---|---|---|---|---|
| R0, 8 | 0.462 (8/9) | 0.517 (9/9) | 0.500 (9/9) | 1.214 (3/9) | 1.241 (3/9) | 0.971 (5/9) |
| R0, 16 | 0.529 (9/9) | 0.590 (9/9) | 0.628 (9/9) | 0.971 (5/9) | 0.988 (6/9) | 0.921 (7/9) |
| R5, 4 | 0.620 (8/9) | 0.647 (9/9) | 0.710 (9/9) | 0.902 (7/9) | 0.984 (5/9) | 1.021 (4/9) |
| R5, 8 | 0.247 (8/9) | 0.312 (9/9) | 0.341 (9/9) | 1.020 (4/9) | 0.880 (7/9) | 0.678 (8/9) |

(ATE; parantez: oranı < 1 olan sekans / 9.) Kural dört hücrenin dördünde sağlandı.

Okuma:
- Kazancın kaynağı **arka-uç ayarı**: aynı ön-uçta (fe0) A0g / A0 = 0.31–0.65, 9 / 9 sekans. Ön-uç kapısı (fe1) tek
  başına tutarlı değil (R0 D/r₀ 8'de A0 ile 1.24: kötüleştiriyor; yalnız R5 D/r₀ 8'de A0g ile birlikte 0.68).
- **Kovaryansı şişirmek tek başına doğruluğu değiştirmiyor** (A1 / A0 = 0.90–1.21; ana ızgaranın H5 kararıyla aynı), aynı
  ölçekle kapıları da genişletmek hatayı 0.34–0.71 katına indiriyor (A0g / A1, 9 / 9). Yani doğruluk kaybının baskın
  nedeni, türbülansla kaymış ama geçerli gözlemlerin sabit kapılarca atılması.
- A0g üç bayrağı birlikte değiştirir (`--px_noise`, `--stereo_sigma`, `--pnp_reproj`); hangisinin taşıdığı bu sondayla
  ayrılamaz (ayırma koşumu gerekir).
- Tutarlılık: A0g'nin konum ANEES'i A0'ın 0.02–0.21 katı; medyan düzey 27–112 (hâlâ iyimser). A1'e göre
  0.27–0.56 (R5 D/r₀ 4 hariç: 1.05).
- Başarısız pencere: A0/fe0 hiçbir hücrede ıraksamadı; R0 D/r₀ 16'da A0g/fe0 ve A1/fe0 2 / 27 ıraksadı (aynı iki
  pencere). Şişirilmiş gürültü modeli D/r₀ 16'da sağlamlığı artırmıyor (iki pencerede bozuyor), kalan pencerelerde
  doğruluğu artırıyor. S1'in "başarısız pencere yarıya iner" kolu sağlanmadı; S1 ATE koluyla tuttu.

Düzeyler (medyan ATE [m] / medyan konum ANEES'i / ıraksayan pencere, n = 27):

| rejim, D/r₀, ön-uç | A0 | A0g | A1 |
|---|---|---|---|
| R0, 8, fe0 | 0.1339 / 909.8 / 0 | 0.1052 / 39.2 / 0 | 0.1974 / 90.3 / 0 |
| R0, 8, fe1 | 0.2031 / 941.3 / 0 | 0.1070 / 36.4 / 0 | 0.1410 / 63.6 / 0 |
| R0, 16, fe0 | 0.5767 / 7593.0 / 0 | 0.2711 / 106.8 / 2 | 0.4700 / 203.5 / 2 |
| R0, 16, fe1 | 0.5227 / 6458.6 / 0 | 0.2427 / 111.7 / 1 | 0.5277 / 234.0 / 1 |
| R5, 4, fe0 | 0.0550 / 176.6 / 0 | 0.0355 / 33.3 / 0 | 0.0455 / 35.0 / 0 |
| R5, 4, fe1 | 0.0491 / 151.2 / 0 | 0.0347 / 41.7 / 0 | 0.0498 / 31.6 / 1 |
| R5, 8, fe0 | 0.2609 / 1572.4 / 0 | 0.1219 / 41.0 / 0 | 0.2403 / 132.5 / 0 |
| R5, 8, fe1 | 0.2087 / 1053.3 / 1 | 0.0761 / 27.1 / 0 | 0.2204 / 83.8 / 0 |

(Temiz başvuru, ana ızgara tohum 0, A0: medyan ATE 0.0257 m.)

S1'in sonucu (ön beyan): "çevrimiçi kestirimli (stereo diferansiyel dikey kayma = DIMM) uyarlamalı kapı kolu yazılır;
doğrulayıcı ızgara görülmemiş tohumlarda, ON_KAYIT §5'e sapma olarak **ve kullanıcı onayıyla**." Karar kullanıcıda.

### S2 (rejim, doğruluk) — TUTTU
Kural: ATE(A0, R5) / ATE(A0, R0), D/r₀ 8'de > 1.3. Ölçülen: **1.658** (9 sekansın 9'unda > 1; aralık 1.34–2.83).
A1 için 1.439 (9 / 9 > 1).

### S3 (rejim, tutarlılık) — TUTMADI (ters yön)
Kural: ANEES(A1, R5) < 0.5 × ANEES(A1, R0), D/r₀ 8. Ölçülen: **1.717** (9 sekansın 8'inde > 1). A0 için 2.238 (9 / 9),
A0g için 1.013 (4 / 9 < 1). "Gürültü beyaza yaklaşınca skaler şişirme daha tutarlı olur" beklentisi D/r₀ 8'de yanlış:
hızlı rejimde hata, kovaryanstan daha çok büyüyor.

### Ek (ön beyanda yoktu; Sapma 5 öngörülerini yazmak için; ana ızgara tohum 0 ile eşli)
`scratchpad/sonda_ek.py` (medyan [en küçük, en büyük sekans oranı], oranı < 1 olan sekans):

| ölçü | D/r₀ 4 | D/r₀ 8 |
|---|---|---|
| A0 ATE, R5 / temiz | 2.81 [1.66, 5.16] (0/9) | 14.46 [6.75, 44.28] (0/9) |
| A0 ATE, R0 / temiz (ana ızgara) | 3.20 [1.60, 4.94] | 10.52 [3.27, 15.89] |
| A0 ATE, R5 / R0 | 1.00 [0.58, 1.61] (4/9) | 1.66 [1.34, 2.83] (0/9) |
| A0 konum ANEES'i, R5 / R0 | 0.63 [0.49, 1.82] (7/9) | 2.24 [1.62, 4.50] (0/9) |
| A1 ATE, R5 / R0 | 0.80 [0.57, 1.42] (6/9) | 1.44 [1.14, 2.36] (0/9) |
| A1 konum ANEES'i, R5 / R0 | 0.85 [0.60, 1.73] (7/9) | 1.72 [0.65, 2.78] (1/9) |
| A1 / A0 konum ANEES'i, R5 | 0.19 [0.12, 0.27] (9/9) | 0.07 [0.04, 0.09] (9/9) |
| A1 / A0 ATE, R5 | 0.90 [0.63, 1.44] (7/9) | 1.02 [0.78, 1.31] (4/9) |

Rejim etkisi düzeye bağlı: D/r₀ 4'te doğruluk değişmiyor, tutarlılık biraz iyileşiyor; D/r₀ 8'de ikisi de belirgin
kötüleşiyor. Ana ızgaranın yarı-statik rejimi D/r₀ 8'de iyimser görünüyor (tohum 0; doğrulanması Sapma 5'in işi).

## 3. Sınırlar
- Tek tohum, 9 birim; güven aralığı ve test yok. Sapma 3 ekindeki sayısal duyarlılık (tek koşumun pencere sonucu
  1e-12 girdi farkına duyarlı) burada da geçerli; 9 / 9 sekansta aynı yön ve 2–4 katlık etki bunun üstünde, ama
  0.9–1.2 arasındaki oranlar (fe1 / fe0, A1 / A0) gürültüden ayırt edilemez.
- A0g ve fe1 oracle σ kullanır; çevrimiçi kestirimle (DIMM σ̂ ana ızgarada oracle değerin 0.69–0.94 katı) kazancın
  ne kadarının korunacağı bilinmiyor.
- Bu sonuçlar görüldü: aynı hücrelerde (tohum 0) yapılacak hiçbir test doğrulayıcı sayılamaz.
