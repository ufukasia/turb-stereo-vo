# Gerçek türbülans verisiyle kestirici denetimi — sonuç

Yazım: 2026-11-12 00:54. Ön bildirim: ON_BILDIRIM.md (c6ce9c58…), ek: DEGISIKLIK_1.md (9faaf6f3…).
Sayıların kaynağı: saha_ozet.json (saha_temmuz.jsonl + saha_ekim.jsonl). Kestirici kodu
değişmedi (code_fp 08e28b83407c, koşum öncesi ve sonrası aynı).
Statü: dış geçerlilik denetimi; VO doğruluk kanıtı değil, ön kayıtlı aileye girmez.

## Birincil: Temmuz 2021, L = 724 m, ortalama-görüntü başvurusu (119/139 dakika)

| No | Ölçüt | Beklenen | Ölçülen [%95 GA] | Sonuç |
|----|-------|----------|------------------|-------|
| G1 | Spearman ρ(σ̂², Cn²_sint) | 0.6–0.9 | 0.863 [0.585, 0.903] | GEÇTİ |
| G2 | log–log eğim (kuram 1.0) | 0.6–1.0 | 0.833 [0.665, 0.883] | GEÇTİ |
| G3 | ortanca Cn²_kest/Cn²_sint (küresel dalga kuramı 0.375) | 0.15–1.0 | 0.399 [0.366, 0.413] | GEÇTİ |
| G4 | uydurma R² ortancası; yapışma; θ_c | ≥ 0.9; ≤ %30; 5–60 px | 0.615; %2.5; 232 px | KALDI |

σ̂ / σ_kuram(küresel dalga), Cn² üçte birliklerine göre: 1.08 / 1.00 / 0.98
(σ̂ ortanca 0.57 / 0.81 / 1.06 px).

## G4 neden kaldı

Gerçek yapı fonksiyonu 180 px'e kadar doygunlaşmadan büyüyor (ortanca D(sep)/D(son):
0.17, 0.35, 0.48, 0.62, 0.72, 0.78, 0.83, 0.87, 0.90, 0.93, 0.97, 1.00). Benzetimde
ayarlanan çekirdek (A = 0.75, p = 0.40) bu biçimi iyi tutmuyor ve korelasyon ölçeği
beklenenden (5–60 px) çok uzun (232 px). Yani kestiricinin GENLİĞİ gerçek veride
doğrulanıyor, uzamsal korelasyon BİÇİMİ doğrulanmıyor. Olası nedenler (sınanmadı):
yatay yol boyunca dağılmış türbülans (benzetim tek katman), yazarların özniteliğe
dayalı ikinci sabitleme adımı, dış ölçek. Makalede sınır olarak yazılacak.

## Duyarlılık (betimsel)

- D = 0.03 m: G3 oranı 0.316. LK penceresi 21: ρ 0.869, eğim 0.823, oran 0.350, R² 0.568.
- Sınıra yapışan dakikalar hariç (116): ρ 0.863, eğim 0.837, oran 0.397.
- Kare-başvurulu kol, Temmuz (127 dakika): ρ 0.861, eğim 0.835, oran 0.443, R² 0.637.
- Doğrudan RMS artık ile Cn²: ρ 0.580 (fark ölçümü daha iyi izliyor).

## İkincil: Ekim 2021, L = 1450 m

- Ön bildirimdeki yöntem: 5/166 dakika → kullanılamadı.
- Kare-başvurulu kol (Değişiklik 1): 51/166 dakika; ρ 0.665 [0.071, 0.907], eğim 0.750
  [0.066, 0.964], oran 0.438 [0.341, 0.512], R² 0.663 (G4 yine kalıyor). Dakikaların
  %69'u dışlandığı için zayıf, destekleyici kanıt.

## Sınırlar

- Temmuz'da dışlanan 20 dakika rastgele değil: hepsi 11:11–11:48 arasında, Cn² ortancası
  3.3e-14 (tutulanlar 1.1e-14). Yani aralığın en güçlü ucu (Cn² > ~3e-14) sınanmadı;
  son dakikalarda panonun önünde insanlar var.
- Tek kameralı sabit telefoto; D/r₀ ≈ 0.9–4; stereo ve yörünge yok.
- Açıklık (0.03/0.06 m) ve kırpım ölçeği belirsiz; G3 oranı buna %25 düzeyinde duyarlı.
- Dakika başına 1–2 s'lik patlama; sintilometre aygıt belirsizliği.
- OTIS çözümlemesi (betimsel) yapılmadı; Cn² içermediği için ek kanıt değeri düşük.
