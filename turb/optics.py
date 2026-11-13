"""
Atmosferik türbülansın optik ilişkileri (Kolmogorov).

Bu modül SADECE kapalı-form fiziksel bağıntıları içerir; rastgele alan üretimi
`screen.py` ve `tiltfield.py`'dedir. Buradaki bağıntılar, simülatörün ürettiği
istatistiğin DOĞRULANMASI için bağımsız bir referans oluşturur (bkz.
`validate_simulator.py`) — yani simülatör bu formüllerle üretilmez, bu
formüllerle SINANIR. E2'deki model doğrulamasının dairesel olmamasının nedeni budur.

Semboller (proje raporu §2.1):
    r0    Fried parametresi [m]
    D     açıklık çapı [m]
    theta0 izoplanatik açı [rad]
    tau0  koherens süresi [s]
    lam   dalga boyu [m]
"""
from __future__ import annotations

import numpy as np

# --- Tilt varyansı katsayıları (tek eksen) -----------------------------------
# sigma_alpha^2 = C * (D/r0)^(5/3) * (lambda/D)^2      [rad^2]
#
# Z-TILT (Zernike/tam açıklık en-küçük-kareler eğimi): Noll'un artık varyanslarından
#   Delta_1 - Delta_3 = 1.0299 - 0.134 = 0.8959 (D/r0)^{5/3}, iki eksene eşit bölünür
#   => sigma_{a2}^2 = 0.448 (D/r0)^{5/3} [rad^2 faz]
#   a2 -> açı dönüşümü: W = a2*Z2/k, Z2 = 2 rho cos(phi) => alpha = 4 a2 / (k D)
#   => sigma_alpha^2 = 16*0.448/(4 pi^2) (D/r0)^{5/3} (lam/D)^2 = 0.1816 (...)
# G-TILT (açıklık üzerinde ORTALAMA dalga cephesi eğimi, "gradient tilt"):
#   Sasiela/Tyler: 0.170 (D/r0)^{5/3} (lam/D)^2
#
# Proje raporu Ek B/5 bu ayrımın netleştirilmesini istiyordu: simülatörümüz
# açıklık üzerinde gradyan ORTALAMASI aldığı için ürettiği büyüklük G-tilt'tir ve
# doğrulama referansı C_GTILT'tir. Rapordaki 0.170 katsayısı da G-tilt'e aittir.
C_GTILT = 0.170
C_ZTILT = 0.1816


def sigma_tilt_rad(D_over_r0: float, lam: float, D: float, kind: str = "gtilt") -> float:
    """Tek eksen görüntü-hareketi (tilt) standart sapması [rad].

    sigma = sqrt(C) * (D/r0)^(5/6) * (lam/D)
    """
    C = C_GTILT if kind == "gtilt" else C_ZTILT
    return float(np.sqrt(C) * D_over_r0 ** (5.0 / 6.0) * (lam / D))


def sigma_tilt_px(D_over_r0: float, lam: float, D: float, px_per_rad: float,
                  kind: str = "gtilt") -> float:
    """Tek eksen tilt jitteri [piksel] — px_per_rad = f/p (odak uzaklığı / piksel adımı)."""
    return sigma_tilt_rad(D_over_r0, lam, D, kind) * px_per_rad


def r0_from_D_over_r0(D_over_r0: float, D: float) -> float:
    """r0 [m]."""
    if D_over_r0 <= 0:
        return np.inf
    return D / D_over_r0


def r0_from_cn2(cn2: float, L: float, lam: float) -> float:
    """Düzgün (homojen) yatay yol için Fried parametresi.

    r0 = (0.423 k^2 Cn^2 L)^(-3/5),  k = 2 pi / lam
    """
    k = 2.0 * np.pi / lam
    return float((0.423 * k * k * cn2 * L) ** (-3.0 / 5.0))


def cn2_from_r0(r0: float, L: float, lam: float) -> float:
    """r0 -> Cn^2 (düzgün yatay yol)."""
    k = 2.0 * np.pi / lam
    return float(r0 ** (-5.0 / 3.0) / (0.423 * k * k * L))


def isoplanatic_angle(r0: float, z_layer: float) -> float:
    """Tek eşdeğer katman için izoplanatik açı [rad]: theta0 = 2.914^(-3/5) r0 / z ≈ 0.314 r0/z.

    (Homojen yol için theta0 = 2.914 k^2 Cn^2 L^{8/3})^{-3/5}; tek ince katman
    yaklaşımında bu theta0 = 0.314 r0 / z_layer'a indirgenir.)
    """
    if z_layer <= 0:
        return np.inf
    return float(0.314 * r0 / z_layer)


def tilt_correlation_angle(D: float, z_layer: float) -> float:
    """TILT için açısal korelasyon ölçeği [rad].

    Yüksek-mertebe faz (theta0) ile TILT'in korelasyon açısı FARKLIDIR: tilt,
    açıklık üzerinde ortalandığı için korelasyonunu yitirmesi ışın ayrımı
    z*dtheta'nın AÇIKLIK çapı D ile kıyaslanabilir olmasına bağlıdır:
        theta_c ≈ D / z_layer
    Anizoplanatizmin tilt-özel ölçeği budur ve E1-c'de taranan eksendir.
    """
    if z_layer <= 0:
        return np.inf
    return float(D / z_layer)


def layer_distance_for_theta_c(D: float, theta_c: float) -> float:
    """Verilen tilt korelasyon açısını üreten eşdeğer katman mesafesi [m]."""
    if theta_c <= 0:
        return np.inf
    return float(D / theta_c)


def tau0(r0: float, v_wind: float) -> float:
    """Koherens süresi [s]: tau0 = 0.314 r0 / v."""
    if v_wind <= 0:
        return np.inf
    return float(0.314 * r0 / v_wind)


def wind_for_tau0(r0: float, tau0_target: float) -> float:
    """Hedef tau0'ı veren enine rüzgâr hızı [m/s]."""
    if tau0_target <= 0:
        return np.inf
    return float(0.314 * r0 / tau0_target)


def rho_theoretical(dtheta: np.ndarray, theta_c: float, exponent: float = 5.0 / 3.0) -> np.ndarray:
    """Mühendislik korelasyon modeli (proje raporu §4.2):

        rho(dtheta) = exp[ -(dtheta/theta_c)^(5/3) ]

    E2'de bu model, fizik-temelli faz ekranından ölçülen ampirik rho_hat'e
    UYDURULUR; katsayı uyumu (R^2) katkının doğrulamasıdır.
    """
    dtheta = np.asarray(dtheta, dtype=np.float64)
    with np.errstate(over="ignore"):
        return np.exp(-np.power(np.abs(dtheta) / max(theta_c, 1e-30), exponent))


def short_exposure_blur_sigma_rad(r0: float, lam: float, D: float) -> float:
    """Kısa-poz (tilt çıkarılmış) bulanıklık için eşdeğer Gauss sigma [rad].

    Uzun-poz seeing FWHM ≈ 0.98 lam/r0. Tilt çıkarıldığında kalan kısa-poz
    genişliği yaklaşık olarak açıklık kırınımı ile türbülans arasındaki
    kareler-toplamıdır; D >> r0 rejiminde
        FWHM_SE ≈ 0.98 (lam/r0) * [1 - 0.37 (r0/D)^{1/3}]^{1/2}
    (Fried 1966). sigma = FWHM / 2.3548.
    """
    if not np.isfinite(r0) or r0 <= 0:
        return 0.0
    fwhm_le = 0.98 * lam / r0
    corr = max(1.0 - 0.37 * (r0 / D) ** (1.0 / 3.0), 0.0)
    fwhm_se = fwhm_le * np.sqrt(corr)
    return float(fwhm_se / 2.3548)


def rytov_variance(cn2: float, L: float, lam: float) -> float:
    """Düzlem dalga Rytov varyansı: sigma_R^2 = 1.23 Cn^2 k^{7/6} L^{11/6}."""
    k = 2.0 * np.pi / lam
    return float(1.23 * cn2 * k ** (7.0 / 6.0) * L ** (11.0 / 6.0))


def scintillation_index_weak(sigma_R2: float) -> float:
    """Zayıf türbülansta sintilasyon indeksi ≈ Rytov varyansı (doygunluk yok)."""
    return float(sigma_R2)


def fresnel_scale(lam: float, L: float) -> float:
    """Fresnel ölçeği sqrt(lam L) [m] — sintilasyon benek boyutu."""
    return float(np.sqrt(max(lam * L, 0.0)))
