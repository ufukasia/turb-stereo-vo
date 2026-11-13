# Gerçek türbülans verisiyle doğrulama — ön bildirim

Yazım zamanı: 2026-11-11 23:56:36. Bu belge, hiçbir görüntü istatistiği (eğim, yapı fonksiyonu, kestirim)
hesaplanmadan ÖNCE yazıldı. O ana kadar bakılanlar: dosya yapısı, ReadMe, Cn² tablosu
(bağımsız değişken), iki dakikadan dörder örnek kare (içerik türü için) ve altı kare
arasında faz korelasyonu (karelerin aynı sahne olup olmadığı için).

Statü: benzetici ve çevrimiçi kestirici için **dış geçerlilik denetimi**. VO doğruluk
kanıtı DEĞİLDİR; ön kayıtlı hipotez ailesine (ON_KAYIT_FAZ2.md) girmez, p-değeri
ailesini değiştirmez. Makalede "gerçek veriyle kestirici doğrulaması" olarak ve
sınırlarıyla birlikte raporlanır; sonuç ne çıkarsa çıksın raporlanır.

## Veri

- Saha ve ark., "Turbulence strength Cn² estimation from video using physics-based deep
  learning", Optics Express 30, 40854 (2022); arXiv 2408.16623. Lisans belirtilmemiş:
  yalnız çözümlenir ve atıf yapılır, depoya konmaz.
- Temmuz 2021 (birincil): L = 724 m, tek hedef panosu, 139 dakika, Cn² 2.4e-15 – 6.2e-14.
- Ekim 2021 (ikincil): L = 1450 m, 166 dakika, Cn² 2.1e-15 – 7.3e-14. Örnek karede
  doku az; yeterli öznitelik çıkmazsa dakika dışlanır ve sayı raporlanır.
- Birim: bir dakika = bir sintilometre kaydı (Kipp & Zonen LAS MkII) + 120–240 kare
  (120 fps patlama), 256×256 kırpım, yazarlarca iki aşamalı sabitleme uygulanmış.
- Optik (makaleden): piksel görüş açısı 5.95 µrad, D = 0.06 m. ReadMe f/18 ve 539 mm
  diyor (D ≈ 0.03 m); çelişki var. Birincil D = 0.06 m, duyarlılık D = 0.03 m.
  λ = 0.55 µm. Kırpımın yeniden ölçeklenmediği VARSAYILIR (doğrulanamıyor).

## Yöntem (dondurulmuş)

1. Dakikanın kareleri DSCN numarasına göre sıralanır; ortalama görüntü alınır.
2. Ortalama görüntüde Shi–Tomasi köşeleri (en çok 300, kalite 0.02, en az aralık 5 px).
3. Ortalama görüntüden her kareye piramitli LK (pencere 15, 2 düzey); ileri–geri hata
   < 0.5 px; karelerin ≥ %90'ında izlenen öznitelikler tutulur.
4. Artık = konum − özniteliğin dakika ortalaması. Kare başına `add_frame(uv, artık)`.
5. Kestirici: makaledeki `turb.dimm.DifferentialMotionMonitor`, kod DEĞİŞMEDEN;
   yalnız kurucu parametreleri: px_per_rad = 1/5.95e-6, D, λ, n_bins = 12,
   max_sep_px = 180, window_frames = dakikadaki kare sayısı.
6. Cn²_kest = `turb.optics.cn2_from_r0(r0_kest, L, λ)` (kodun düzlem dalga kuralı).
7. Dışlama: < 30 öznitelik, < 60 kare ya da kestirim NaN. Sayılar raporlanır.
8. Güven aralığı: dakikalar üzerinde 10 dakikalık bloklarla özyükleme (2000 tekrar).

## Beklenen değerler ve geçme ölçütleri (Temmuz birincil)

Kuram (sintilometre Cn²'sinden, D = 0.06 m): tek eksen eğim sigması düzlem dalga
0.60 / 1.37 / 3.03 px (en düşük / ortanca / en yüksek Cn²), küresel dalga (hedef
noktası için doğrusu; varyans 3/8) 0.37 / 0.84 / 1.86 px. D/r₀ ≈ 0.9 – 6.5.

| No | Ölçüt | Beklenen | Geçer |
|----|-------|----------|-------|
| G1 | σ̂² ile sintilometre Cn² arasında Spearman ρ | 0.6 – 0.9 | ρ ≥ 0.5 ve %95 GA alt ucu > 0 |
| G2 | log σ̂² – log Cn² eğimi (kuram 1.0; gürültü tabanı düşürür) | 0.6 – 1.0 | nokta kestirim 0.5 – 1.3 |
| G3 | ortanca Cn²_kest / Cn²_sint (kuram, küresel dalga: 0.375) | 0.15 – 1.0 | 0.12 – 1.2 |
| G4 | yapı fonksiyonu doygunlaşan biçimde: uydurma R² ortancası; sınıra yapışma oranı; θ_c | R² ≥ 0.9; ≤ %30; 5 – 60 px | R² ≥ 0.8 ve yapışma ≤ %50 |

Ekim: aynı ölçütler, ikincil; ayrı raporlanır (kuramsal ortanca sigma küresel ≈ 1.6 px).
OTIS: Cn² yok; yalnız betimsel (yapı fonksiyonu biçimi), geçer/kalır yok.

Duyarlılık (betimsel): D = 0.03 m; LK penceresi 21; sınıra yapışan dakikalar hariç;
doğrudan RMS artık (sabitleme tek noktadaki ortak hareketi siler, bu yüzden yalnız betimsel).

## Bilinen sınırlar (sonuçtan bağımsız yazılır)

- Tek kameralı sabit telefoto; stereo yok, yörünge yok: VO doğruluğu sınanamaz.
- İzleme gürültüsü ve bulanıklık Cn² ile birlikte artar; σ̂'yı aynı yönde şişirebilir.
- 1–2 saniyelik patlama: dakika başına bağımsız eğim örneği az (rüzgâr 2 m/s için ~30).
- Açıklık ve ölçek belirsiz (D 0.03/0.06; kırpım ölçeği varsayım); G3 bu yüzden geniş bantlı.
- Sintilometrenin kendisi de %35'e varan aygıt-arası sapma taşır (makalenin belirttiği).
- Çekirdek sabitleri (A = 0.75, p = 0.40) benzetimde ayarlandı; gerçek yatay yol için
  ayrıca ayarlanmadı (ayarlanmayacak; bu denetimin amacı zaten bu).
