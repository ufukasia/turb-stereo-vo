# Değişiklik 1 — kare-başvurulu izleme kolu (ön bildirime ek)

Yazım zamanı: 2026-11-12 00:02:39. ON_BILDIRIM.md (sha256 c6ce9c58…) değiştirilmedi.

## Bu ana kadar görülenler (dürüst kayıt)

Duman testi, ay başına eşit aralıklı 3 dakika (toplam 6):

- Temmuz 09:30 (Cn² 2.9e-15): σ̂ 0.485 px, θ_c 259 px (sınıra yapışık), R² 0.57.
- Temmuz 10:39 (Cn² 6.9e-15): σ̂ 0.629 px, θ_c 250 px, R² 0.65.
- Temmuz 11:48: 0 öznitelik → dışlandı. Neden: karelerde panonun önünde yürüyen
  insanlar var (deney sonu); dışlama yerinde.
- Ekim'in 3 dakikası: 0 / 1 / 19 öznitelik → dışlandı. Neden: dakika içindeki kareler
  birbirine hizalı değil (pano kırpım içinde yer değiştiriyor); ortalama görüntü
  bulanık çıkıyor ve ortalamadan kareye izleme başarısız (ileri–geri hata 1.4–7 px).

Yani ön bildirimdeki yöntem Ekim kümesinde veri yapısı yüzünden çalışmıyor. Bu bir
sonuç değil, bir yöntem uygunsuzluğu; dışlama Cn² ile değil hizasızlıkla ilgili.

## Değişiklik (dondurulmuş; tam koşumdan önce)

- **Birincil çözümleme değişmez:** Temmuz, ortalama-görüntü başvurulu yöntem,
  ölçütler G1–G4 aynen.
- **Ek kol "kare":** başvuru = dakikanın sıralı karelerinin ortadaki karesi; köşeler
  bu karede bulunur (aynı parametreler); LK başvurudan her kareye, piramit 3 düzey,
  pencere 15; ileri–geri < 0.5 px ve ≥ %90 kuralı aynı; artık = konum − özniteliğin
  dakika ortalaması. Ortak kayma çıkarılmaz (kestirici fark ölçer, ortak kaymaya duyarsız).
- Ek kol iki aya da uygulanır: Ekim için ikincil çözümlemenin yöntemi olur; Temmuz için
  duyarlılık (betimsel) olarak raporlanır. Ölçütler ve bantlar aynı.
- Beklenti (ek kol, Temmuz): birincil ile aynı işaret ve benzer büyüklük (G1 farkı < 0.15).
- Görülen iki Temmuz dakikası G4'ün (R² ≥ 0.8) kalabileceğine işaret ediyor; ölçüt
  gevşetilmez, sonuç olduğu gibi raporlanır.
