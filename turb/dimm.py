"""
VIO izlerinden çevrimiçi türbülans kestirimi — DIMM ilkesi (C4, rapor §3/E6).

Astronomideki Differential Image Motion Monitor (DIMM): r0, açısal olarak ayrık iki
noktanın DİFERANSİYEL görüntü hareketinin varyansından ölçülür. Diferansiyel almak,
teleskop titremesi / izleme hatası gibi ORTAK bileşeni yok eder.

VIO karşılığı: bir karedeki gözlemlerin yeniden-projeksiyon ARTIKLARI r_j alınır.
Bu artıklar (i) türbülans eğimi d(theta_j), (ii) dedektör gürültüsü n_j ve
(iii) kestiricinin kendi poz hatası e_poz'u içerir. Poz hatası TÜM özniteliklerde
ortaktır (birinci mertebeden düzgün bir alandır); bu yüzden ÇİFT FARKLARI alınca
düşer:

    D(dtheta) = < |r_a - r_b|^2 >  =  2 sigma^2 [1 - rho(dtheta)]  +  2 sigma_n^2

Bu yapı fonksiyonuna sigma (tilt jitteri) ve theta_c (korelasyon açısı) uydurulur;
oradan r0 kapalı formdan çözülür:

    sigma = sqrt(0.170) (D/r0)^{5/6} (lambda/D) (f/p)   ->   r0 = D (sigma_ref/sigma)^{6/5}

KAPANAN DÖNGÜ: VIO izleri -> r0 kestirimi -> Sigma_yama parametresi -> daha iyi VIO.
"""
from __future__ import annotations

import numpy as np

from . import correlation as _corr
from .covariance import EMP_A, EMP_P
from .optics import C_GTILT


class DifferentialMotionMonitor:
    """Kayan pencereli DIMM kestiricisi.

    Kullanım (kare başına):
        dimm.add_frame(theta_xy, residuals)      # (N,2) rad, (N,2) px
        s, tc, r0 = dimm.estimate()
    """

    def __init__(self, px_per_rad: float, aperture_D: float, wavelength: float,
                 z_layer: float | None = None, n_bins: int = 12,
                 max_sep_px: float = 400.0, window_frames: int = 60,
                 min_pairs: int = 400, max_pairs_per_frame: int = 4000,
                 base_sigma_px: float = 0.5, seed: int = 0):
        self.px_per_rad = float(px_per_rad)
        self.D = float(aperture_D)
        self.lam = float(wavelength)
        self.z_layer = z_layer
        self.base_sigma_px = float(base_sigma_px)
        self.window = int(window_frames)
        self.min_pairs = int(min_pairs)
        self.max_pairs = int(max_pairs_per_frame)
        self.rng = np.random.default_rng(seed)
        # ayrım kutuları [piksel cinsinden açısal ayrım]
        self.edges = np.linspace(0.0, max_sep_px, n_bins + 1)
        self.centers = 0.5 * (self.edges[1:] + self.edges[:-1])
        self._hist = []          # [(sum_D, count)] kare başına
        self.n_frames = 0

    # ------------------------------------------------------------------ ekleme
    def add_frame(self, uv: np.ndarray, residuals: np.ndarray):
        """Bir karenin artıklarını ekle. uv: (N,2) piksel, residuals: (N,2) piksel."""
        uv = np.asarray(uv, dtype=np.float64)
        r = np.asarray(residuals, dtype=np.float64)
        n = len(uv)
        if n < 4:
            return
        # çift örnekleme (N^2 yerine rastgele alt küme)
        npair = min(self.max_pairs, n * (n - 1) // 2)
        ia = self.rng.integers(0, n, size=npair)
        ib = self.rng.integers(0, n, size=npair)
        m = ia != ib
        ia, ib = ia[m], ib[m]
        sep = np.hypot(uv[ia, 0] - uv[ib, 0], uv[ia, 1] - uv[ib, 1])
        dd = np.sum((r[ia] - r[ib]) ** 2, axis=1)          # |r_a - r_b|^2 (2 eksen)
        idx = np.digitize(sep, self.edges) - 1
        ok = (idx >= 0) & (idx < len(self.centers))
        s = np.zeros(len(self.centers))
        c = np.zeros(len(self.centers))
        np.add.at(s, idx[ok], dd[ok])
        np.add.at(c, idx[ok], 1.0)
        self._hist.append((s, c))
        if len(self._hist) > self.window:
            self._hist.pop(0)
        self.n_frames += 1

    # ---------------------------------------------------------------- kestirim
    def structure_function(self):
        if not self._hist:
            return self.centers, np.full(len(self.centers), np.nan), np.zeros(len(self.centers))
        S = np.sum([h[0] for h in self._hist], axis=0)
        C = np.sum([h[1] for h in self._hist], axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            Dv = S / C
        return self.centers, Dv, C

    # v2 (2026-08-25): ampirik cekirdek + bilesen sayimi duzeltmesi.
    # Gerekce ve Monte Carlo dogrulamasi: evidence/DIMM_KUSURU.md
    last_pinned = False
    last_bounds = (float("nan"), float("nan"))

    def estimate(self, theta_c_px: float | None = None):
        """(sigma_px, theta_c_px, r0_m). Yetersiz veri -> (nan, nan, nan)."""
        sep, Dv, C = self.structure_function()
        ok = np.isfinite(Dv) & (C > 10)
        if ok.sum() < 3 or C.sum() < self.min_pairs:
            return np.nan, np.nan, np.nan

        # KOVARYANS YOLUYLA AYNI ÇEKİRDEK (turb/covariance.py EMP_A, EMP_P).
        # rho(0) = A < 1: nugget (dedektör gürültüsü + korelasyonsuz kalıntı)
        # genliğin İÇİNDEDİR, bu yüzden AYRI bir gürültü terimi EKLENMEZ —
        # eklemek nugget'i iki kez sayardı (2026-08-20'de kovaryans yolunda
        # bulunan ve düzeltilen hatanın aynısı).
        # BILESEN SAYIMI (2026-08-23'te bulunan IKINCI kusur, ilkinden bagimsiz).
        # add_frame `dd = sum((r_a - r_b)**2, axis=1)` ile IKI EKSENI toplar:
        #     D(sep) = E|dr|^2 = 2 * [2 sigma^2 (1-rho)] = 4 sigma^2 (1-rho)
        # Eski model katsayisi 2.0 idi -> s2 = 2 sigma^2 uyduruyor ->
        # sigma_hat = sqrt(2) * sigma_gercek, yani SISTEMATIK %41 fazla kestirim.
        # E4'te olculen sigma_hat/sigma_oracle = 1.25-2.01 araliginin bir kismi
        # tam olarak budur (sqrt(2) = 1.414); kalani cekirdek uyusmazligindandi.
        # Bu, r0 kestirimini de (sigma^-6/5 ile) yanli kilar -> E6 eksenini etkiler.
        def model(sep_px, s2, tc_px):
            x = np.maximum(sep_px, 0.0) / max(tc_px, 1e-6)
            return 4.0 * s2 * (1.0 - EMP_A * np.exp(-(x ** EMP_P)))

        from scipy.optimize import curve_fit
        x, y, w = sep[ok], Dv[ok], np.sqrt(C[ok])
        try:
            if theta_c_px is None:
                # ÖLÇÜLEBİLİRLİK SINIRI korunur ama gevşetilir: alt sınır artık
                # 2x değil 1x ilk ayrım (0.5x güvenlik payıyla), üst sınır
                # 4x en büyük ayrım. Gerçek theta_c (~113 px) bu aralığın
                # ORTASINDA kalmalı; sınıra dayanırsa P14 kırmızı verir.
                # OLCULEBILIRLIK: en uzun ayrimdan UZUN bir korelasyon boyu
                # olculemez -- yapi fonksiyonu o araligin otesinde duz olur ve
                # uydurma gurultuye oturur (olculdu: tc_hi = 4*max_sep iken
                # gercek 200 px icin 1533 px donuyordu, %667 hata).
                # Ust sinir 1.5 * en uzun ayrim.
                tc_lo = max(0.5 * float(np.min(x)), 2.0)
                tc_hi = max(1.5 * float(np.max(x)), tc_lo * 4.0)
                # COK BASLANGICLI UYDURMA. Tek p0 ile yuksek sigma'da yerel
                # minimuma dusuyordu (olculdu: theta_c=112.8 gercekte iken 63.9,
                # theta_c=200 iken 282 veriyordu). Birkac baslangic denenip
                # agirlikli kalinti karesi EN KUCUK olan secilir.
                best = None
                _s0 = max(np.nanmax(y) / 4.0, 0.1)
                for _f in (0.25, 0.5, 1.0, 2.0, 4.0):
                    _tc0 = float(np.clip(np.nanmedian(x) * _f, tc_lo, tc_hi))
                    try:
                        _po, _ = curve_fit(model, x, y, p0=[_s0, _tc0],
                                           sigma=1.0 / w,
                                           bounds=([1e-6, tc_lo], [1e4, tc_hi]),
                                           maxfev=20000)
                    except Exception:
                        continue
                    _r = float(np.sum(w * (model(x, *_po) - y) ** 2))
                    if best is None or _r < best[0]:
                        best = (_r, _po)
                if best is None:
                    return np.nan, np.nan, np.nan
                popt = best[1]
                s2, tc = float(popt[0]), float(popt[1])
                self.last_bounds = (tc_lo, tc_hi)
                # SINIRA YAPISMA TESPITI — v1'in ASIL KUSURU BUYDU.
                # v1 her hucrede alt siniri (33.3 px) donduruyordu ve bunu
                # OLCUM gibi raporluyordu; hicbir sey uyarmiyordu. Sinira
                # yapisan bir kestirim bir OLCUM DEGIL, "olculemedi"dir.
                self.last_pinned = bool(tc <= tc_lo * 1.01 or tc >= tc_hi * 0.99)
            else:
                f = lambda xx, s2_: model(xx, s2_, theta_c_px)     # noqa: E731
                popt, _ = curve_fit(f, x, y, p0=[max(np.nanmax(y) / 2.0, 0.1)],
                                    sigma=1.0 / w, bounds=([1e-6], [1e4]),
                                    maxfev=20000)
                s2, tc = float(popt[0]), float(theta_c_px)
                self.last_bounds = (np.nan, np.nan)
                self.last_pinned = False
        except Exception:
            return np.nan, np.nan, np.nan

        sigma_px = float(np.sqrt(max(s2, 0.0)))
        r0 = r0_from_sigma_px(sigma_px, self.D, self.lam, self.px_per_rad)
        return sigma_px, tc, r0


def r0_from_sigma_px(sigma_px: float, D: float, lam: float, px_per_rad: float) -> float:
    """sigma_px -> r0 [m] (G-tilt kapalı formunun tersi)."""
    if not np.isfinite(sigma_px) or sigma_px <= 0:
        return np.nan
    sigma_rad = sigma_px / px_per_rad
    # sigma = sqrt(C) (D/r0)^{5/6} (lam/D)  ->  (D/r0) = [sigma D/(sqrt(C) lam)]^{6/5}
    dr = (sigma_rad * D / (np.sqrt(C_GTILT) * lam)) ** (6.0 / 5.0)
    return float(D / dr) if dr > 0 else np.nan


def sigma_px_from_r0(r0: float, D: float, lam: float, px_per_rad: float) -> float:
    if not np.isfinite(r0) or r0 <= 0:
        return np.nan
    return float(np.sqrt(C_GTILT) * (D / r0) ** (5.0 / 6.0) * (lam / D) * px_per_rad)
