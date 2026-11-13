# Uyarlanır izleme kapısı kolu — geliştirme koşumu sonucu ve seçim (KEŞİF; ön kayıt dışı)

Yazım: 2026-11-11 20:46. Ön beyan: `ON_BEYAN_3.md` (19:30, sonuçtan önce; sha256 `598ebfac7aba46d8…`). Analiz betiği
`kapi_kol_analiz.py`; ham çıktısı `gelistirme_analiz.txt` (sha256 `bf34585782bd61e5…`).
Veri: tohum 0, R0 (yarı-statik), D/r₀ 8, 9 sekans × 3 pencere; kod parmak izi `bedca2849bb2`.
Kayıtlar: `gelistirme.jsonl` (16 koşum; sha256 `31aef8be4847d2f1…`) + `on_test.jsonl` içindeki V1_01_easy D1 / D2
koşumları (sha256 `a3f5dc4cb64dab71…`). Koşum 19:48:22 – 20:42:28, çıkış kodu 0.

**Statü.** Tohum 0, birim = sekans (9), test ve güven aralığı yok. Buradaki hiçbir sayı doğrulayıcı kanıt değildir;
makalede yalnız "keşif / geliştirme (ön kayıt dışı)" etiketiyle anılabilir. Tek işlevi doğrulayıcı ızgarada
(ON_KAYIT §5 Sapma 6; tohum 1–4) dondurulacak yapılandırmayı seçmektir.

## 1. Bütünlük
54 / 54 kayıt (D1 27, D2 27); hepsi `status` ok, `m_valid`, `n_valid`; ıraksayan 0; eksik koşum 0; tek parmak izi.
Kapı sayacı `gate_sigma_n` = (son pencere başı + 400 − 1) dokuz sekansın hepsinde tuttu (MH01 1920, MH02 1799,
MH03 1361, MH04 1331, MH05 1301, V1 / V2 1299).

## 2. Sonuç (eşli medyan oran; köşeli parantez: en küçük ve en büyük sekans; parantez: oranı < 1 olan sekans)
| oran | ATE | konum ANEES'i |
|---|---|---|
| D1 / A0 | 0.512 [0.42, 0.90] (9 / 9) | 0.174 [0.09, 1.35] (8 / 9) |
| D2 / A0 | 0.522 [0.35, 1.19] (8 / 9) | 0.062 [0.02, 0.51] (9 / 9) |
| D1 / A3S (aynı kovaryans, kapı farkı) | 0.437 [0.28, 0.92] (9 / 9) | 0.514 [0.20, 1.18] (8 / 9) |
| D2 / A3S | 0.430 [0.26, 1.28] (8 / 9) | 0.160 [0.07, 0.44] (9 / 9) |
| D1 / A1 (oracle kovaryans, kapısız) | 0.519 [0.29, 1.17] (8 / 9) | 1.909 [0.86, 7.57] (2 / 9) |
| D2 / A1 | 0.490 [0.26, 1.12] (8 / 9) | 0.648 [0.20, 2.83] (6 / 9) |
| D1 / A0g (oracle ölçekli kapı) | 1.082 [0.91, 1.24] (2 / 9) | 3.993 [3.41, 5.32] (0 / 9) |
| D2 / A0g | 1.009 [0.56, 1.65] (4 / 9) | 1.322 [0.80, 2.62] (1 / 9) |
| D2 / D1 | 0.913 [0.45, 1.60] (6 / 9) | 0.333 [0.18, 0.66] (9 / 9) |
| başvurular: A0g / A0 | 0.517 (9 / 9) | 0.044 |
| A1 / A0 | 1.214 (3 / 9) | 0.101 |
| A3S / A0 | 1.214 (2 / 9) | 0.351 |

Sekans bazında ATE oranı D1 / A0: MH01 0.498, MH02 0.512, MH03 0.447, MH04 0.454, MH05 0.415, V1_01_easy 0.896,
V1_02_medium 0.893, V2_01_easy 0.858, V2_02_medium 0.664 (log oran ortalaması −0.514, SS 0.320).
D2 / A0: MH01 0.795, MH02 0.522, MH03 0.408, MH04 0.352, MH05 0.372, V1_01_easy 0.749, V1_02_medium 0.406,
V2_01_easy 1.190, V2_02_medium 0.660.

Düzeyler (27 pencere; medyan ATE [m] / medyan konum ANEES'i / başarısız pencere): A0 0.1339 / 909.8 / 0;
A1 0.1974 / 90.3 / 0; A3S 0.1729 / 366.2 / 0; A0g 0.1052 / 39.2 / 0; D1 0.1132 / 160.2 / 0; D2 0.1338 / 50.8 / 0.

Kapı sigması medyanı (sekanslar): D1 1.90–2.04 px (sekans medyanı 1.983; en büyük maks 2.282); D2 3.72–4.09 px
(3.922; en büyük maks 4.859). σ̂_DIMM medyanı: D1 1.63–1.79, D2 1.79–1.99 px (kapısız A3S 1.30–1.64; oracle σ_tilt
2.189, oracle toplam 4.412 px). Kapı genişleyince σ̂ yükseliyor (kesilme yanlılığı azalıyor) ama kaçak geri besleme yok.

## 3. Ön beyan edilen beklentilerin sonucu (dürüst döküm)
| no | ölçü | beklenen | ölçülen | sonuç |
|---|---|---|---|---|
| G1 | D2 / A0 ATE | 0.50–0.85, ≥ 7 / 9 | 0.522, 8 / 9 | tuttu |
| G2 | D1 / A0 ATE | 0.60–0.95 | 0.512, 9 / 9 | **tutmadı** — bandın altında (beklenenden iyi); "dar kapı daha küçük kazanç verir" öngörüsü yanlış çıktı |
| G3 | D / A3S ATE | < 1 (D2: 0.5–0.9) | D1 0.437; D2 0.430 | yön tuttu; D2 **bandın altında** (beklenenden büyük etki) |
| G4 | D / A0g ATE | 1.0–1.6 | D1 1.082; D2 1.009 | tuttu |
| G5 | konum ANEES'i D / A0 | < 0.5; D2 < D1 | D1 0.174; D2 0.062 | tuttu |
| G6 | başarısız pencere | ≤ 1 / 27 | 0 / 27; 0 / 27 | tuttu |
| G7 | kapı sigması medyanı | D1 1.6–4.0; D2 2.8–9.0 px; 9 px üstü yok | 1.983; 3.922; maks 4.859 | tuttu |

Okuma (keşif): (i) kazancı kapı taşıyor — aynı kovaryansla kapısız A3S, A0'dan kötü (1.214) iken kapı eklenince
0.512; (ii) çevrimiçi kestirim oracle ölçekli kapıya yakın (D1 / A0g 1.082; 0.9–1.1 bandı sayısal duyarlılıktan
ayırt edilemez); (iii) D1'in kapısı oracle toplam sigmanın yarısından dar (1.98 / 4.41 px) olduğu hâlde ATE kazancı
aynı düzeyde, yani kazanç için kapının "tam" oracle genişliğine ulaşması gerekmiyor; (iv) tutarlılıkta D2 belirgin
iyi (ANEES D2 / D1 0.333, 9 / 9): D1'in kovaryansı hâlâ iyimser (medyan ANEES 160; D2 51; A0 910). Tutarlılık seçim
ölçütü değildir (ön beyan); doğrulayıcı ızgarada ikincil aile I5 olarak raporlanır ve makalede sınır olarak yazılır.

## 4. Seçim (ON_BEYAN_3.md §3 kuralı, harfiyen)
"Geçer": D1 evet (0.512 ≤ 0.85; 9 / 9 ≥ 7; 0 / 27 ≤ 1), D2 evet (0.522; 8 / 9; 0 / 27). İkisi de geçti → kural 2: iki
medyan oranın oranı D2 / D1 = 0.522 / 0.512 = 1.019 ∈ [0.90, 1.10] → **D1** (c = 1; ek sabit yok).
Dondurulan yapılandırma: `A0C+A3 --turb_cal_factor 1.0 --turb_dimm_sigma_scale 1.0 --turb_gate 3`
(`evidence/KAPI/SECIM.json`). D2 doğrulayıcı ızgarada koşulmaz; D2'nin tutarlılık üstünlüğü yalnız bu keşif
dosyasında kalır.

## 5. Sınırlar
- Tek tohum; 9 birim; test yok. MH sekanslarında oran 0.42–0.51, V sekanslarında 0.66–0.90: etki sekansa bağlı.
- Tohum 0 ilk sondada da görülmüştü; seçim aynı veride yapıldı → iyimserlik yanlılığı olası; doğrulayıcı ızgaranın
  (tohum 1–4) etki büyüklüğü daha küçük çıkabilir.
- Karşılaştırma kayıtları (A0, A1, A3S, A0g) eski koddan; T1 bit-aynılığı (15 / 15) gerekçesiyle.
