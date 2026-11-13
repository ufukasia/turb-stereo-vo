"""
Anizoplanatik türbülans alanları: tilt (geometrik bozulma), yüksek-mertebe PSF
(bulanıklık) ve sintilasyon.

TASARIM İLKESİ — dairesellikten kaçınma
---------------------------------------
Bozulma alanı, rapor §4.2'deki rho(dtheta) = exp[-(dtheta/theta_c)^{5/3}] MODELİNDEN
üretilmez. Üretim tamamen fizikseldir:

    1. Kolmogorov faz ekranı (screen.py, D_phi = 6.88 (r/r0)^{5/3} ile doğrulanmış)
    2. Görüş hattı theta yönü için ekran üzerinde AÇIKLIK AYAK İZİ merkezi c = z*theta
    3. O ayak izi üzerinde dalga cephesi GRADYAN ORTALAMASI  ->  G-tilt(theta)
    4. Görüntü kayması [px] = G-tilt [rad] * (f/p)

Anizoplanatizm buradan KENDİLİĞİNDEN doğar: iki yönün ayak izleri z*dtheta kadar
ayrıldığında ortak açıklık alanı azalır, korelasyon düşer. rho(dtheta) bu sürecin
bir SONUCUDUR; E2 bunu ölçüp mühendislik modeline uydurur (R^2 testi).

Katman geometrisi: tek eşdeğer ince katman, z = D / theta_c. Bu, tilt korelasyon
açısını doğrudan denetlemeyi sağlar (bkz. optics.tilt_correlation_angle).

Zamansal evrim: Taylor "donmuş akış" — ekran v = 0.314 r0 / tau0 hızıyla ötelenir.
Ekran periyodik olduğundan öteleme sarmalıdır; tau0 << pencere olduğu hızlı rejimde
bu tekrar istatistiksel olarak etkisizdir (ardışık kareler zaten korelasyonsuzdur).
"""
from __future__ import annotations

import numpy as np

from . import screen as _screen


# =============================================================================
#  Tilt alanı — açıklık-ortalamalı dalga cephesi gradyanı (G-tilt)
# =============================================================================
class TiltField:
    """Ekran üzerinde önceden hesaplanmış (alpha_x, alpha_y) açı alanları [rad].

    Bir kez kurulur (bir tohum = bir atmosfer gerçeklemesi), sonra her kare için
    yalnızca ÖRNEKLEME yapılır — donmuş akış ötelemesi örnekleme koordinatına
    eklenir. Bu, kare başına maliyeti bir bilinear interpolasyona indirir.
    """

    def __init__(self, N: int, delta: float, r0: float, D: float, lam: float,
                 rng: np.random.Generator, n_subharmonics: int = 8,
                 calibration: float = 1.0, L0: float | None = None,
                 keep_phase: bool = False):
        self.N, self.delta, self.r0, self.D, self.lam = N, delta, r0, D, lam
        real = _screen.KolmogorovRealization(N, delta, r0, rng,
                                             n_subharmonics=n_subharmonics, L0=L0)
        self.alpha_x, self.alpha_y = self._aperture_averaged_tilt(real, D, lam)
        self.phz = real.phase() if keep_phase else None
        if calibration != 1.0:
            self.alpha_x *= calibration
            self.alpha_y *= calibration

    @staticmethod
    def _aperture_averaged_tilt(real: "_screen.KolmogorovRealization", D: float, lam: float):
        """G-tilt alanı: alpha_x(c) = (1/(k A)) * integral_{|x-c|<D/2} d(phi)/dx dA.

        Frekans uzayında TAM analitik transfer fonksiyonu ile:
            H_x(f) = 2 pi i f_x * A(f) / k ,  A(f) = 2 J1(pi D f) / (pi D f)
        A(f) normalize dairesel açıklığın Fourier dönüşümüdür (açıklık ortalaması),
        2 pi i f_x ise türevdir. Bu çarpımın integrali yakınsaktır — sonlu fark +
        ayrık disk yaklaşımının ızgaraya bağımlılığı ve Nyquist şişmesi ortadan kalkar.
        Subharmonik (düzlem dalga) modlara aynı H analitik olarak uygulanır.
        """
        k = 2.0 * np.pi / lam
        FX, FY = real.FX, real.FY
        F = np.sqrt(FX ** 2 + FY ** 2)
        A = _airy_aperture(F, D)
        Hx = (2j * np.pi * FX) * A / k
        Hy = (2j * np.pi * FY) * A / k

        def hx(fx, fy):
            return (2j * np.pi * fx) * _airy_aperture(np.hypot(fx, fy), D) / k

        def hy(fx, fy):
            return (2j * np.pi * fy) * _airy_aperture(np.hypot(fx, fy), D) / k

        ax = real.filtered(Hx, hx)
        ay = real.filtered(Hy, hy)
        return ax, ay

    def sample(self, X: np.ndarray, Y: np.ndarray):
        """(X, Y) [m] konumlarında bilinear örnekleme (periyodik sarmalı)."""
        gx = X / self.delta + self.N / 2.0
        gy = Y / self.delta + self.N / 2.0
        return (_bilinear_periodic(self.alpha_x, gx, gy),
                _bilinear_periodic(self.alpha_y, gx, gy))

    def tilt_std(self) -> float:
        """Alanın tek eksen tilt standart sapması [rad] (doğrulama için)."""
        return float(0.5 * (self.alpha_x.std() + self.alpha_y.std()))


def _airy_aperture(f: np.ndarray, D: float) -> np.ndarray:
    """Normalize dairesel açıklığın Fourier dönüşümü: A(f) = 2 J1(pi D f)/(pi D f)."""
    from scipy.special import j1
    x = np.pi * D * np.asarray(f, dtype=np.float64)
    out = np.ones_like(x)
    nz = np.abs(x) > 1e-12
    out[nz] = 2.0 * j1(x[nz]) / x[nz]
    return out


def _bilinear_periodic(F: np.ndarray, gx: np.ndarray, gy: np.ndarray) -> np.ndarray:
    """Periyodik sınırlı bilinear örnekleme (grid koordinatlarında)."""
    N0, N1 = F.shape
    x0 = np.floor(gx).astype(np.int64)
    y0 = np.floor(gy).astype(np.int64)
    fx = gx - x0
    fy = gy - y0
    x0 %= N1; y0 %= N0
    x1 = (x0 + 1) % N1; y1 = (y0 + 1) % N0
    return (F[y0, x0] * (1 - fx) * (1 - fy) + F[y0, x1] * fx * (1 - fy)
            + F[y1, x0] * (1 - fx) * fy + F[y1, x1] * fx * fy)


# =============================================================================
#  Yüksek-mertebe (bulanıklık) alanı — kısa-poz PSF
# =============================================================================
class SpeckleBlurField:
    """Tilt'i çıkarılmış ince pupil ekranından kısa-poz PSF üretimi.

    Neden ikinci (ince) ekran: tilt alanının ekranı geniş kapsama için kaba
    örneklenir (delta = D/8); PSF hesabı ise pupili en az ~32 örnekle çözmek ister.
    İki ölçek AYRIŞTIRILIR: kaba ekran tilt'i, ince ekran yalnız YÜKSEK MERTEBE
    aberasyonu taşır (her pupil yamasından düzlem eğim en-küçük-kareler ile
    çıkarılır). Böylece çift-sayım olmaz ve her iki bileşen de kendi ölçeğinde
    doğru istatistiğe sahiptir.
    """

    #: Sifir dolgu katsayisi: odak duzlemi ornek boyu = lam / (PAD_FACTOR * pupil_px * delta)
    #: = lam / (4 D) rad (pupil_px * delta = D).
    PAD_FACTOR = 4

    def __init__(self, N: int, delta: float, r0: float, D: float, lam: float,
                 rng: np.random.Generator, pupil_px: int = 32, psf_px: int = 21,
                 n_subharmonics: int = 4, px_per_rad: float | None = None,
                 max_px: int | None = None):
        """px_per_rad None ise ESKI davranis: PSF ornekleri dogrudan goruntu pikseli
        sayilir (bir ornek = lam/(4D) rad, yani px_per_rad=170648'de 0.235 px; cekirdek
        ~4.3 kat genis uygulanir). px_per_rad verilirse PSF, ornek boyu
        lam/(4D) * px_per_rad piksel olacak sekilde goruntu izgarasina enerji korunumlu
        (kutu-ortusme) yeniden orneklenir; max_px cekirdek kenarina ust sinirdir."""
        self.N, self.delta, self.D, self.lam = N, delta, D, lam
        self.pupil_px = int(pupil_px)
        self.psf_px = int(psf_px) | 1               # tek sayı
        self.px_per_rad = None if px_per_rad is None else float(px_per_rad)
        self.max_px = None if max_px is None else int(max_px)
        pad = self.PAD_FACTOR * self.pupil_px
        #: bir PSF ornegi [rad] ve [goruntu pikseli]
        self.sample_rad = lam / (pad * delta)
        self.sample_px = (None if self.px_per_rad is None
                          else self.sample_rad * self.px_per_rad)
        self.phz = _screen.phase_screen(N, delta, r0, rng, n_subharmonics=n_subharmonics)
        # Pupil maskesi ve tilt-çıkarma tabanı
        p = self.pupil_px
        yy, xx = np.mgrid[0:p, 0:p] - (p - 1) / 2.0
        self.mask = (xx ** 2 + yy ** 2) <= (p / 2.0) ** 2
        A = np.stack([np.ones(self.mask.sum()), xx[self.mask], yy[self.mask]], axis=1)
        self.pinv = np.linalg.pinv(A)
        self.basis = (np.ones_like(xx), xx, yy)

    def _focal_psf(self, X: float, Y: float):
        """Tilt'i cikarilmis pupil yamasinin sifir-dolgulu FFT yeginligi (pad x pad)."""
        p = self.pupil_px
        i0 = int(round(Y / self.delta + self.N / 2.0))
        j0 = int(round(X / self.delta + self.N / 2.0))
        rows = (np.arange(p) - p // 2 + i0) % self.N
        cols = (np.arange(p) - p // 2 + j0) % self.N
        patch = self.phz[np.ix_(rows, cols)]
        # düzlem (piston + tilt) çıkar — tilt kaba ekrandan geliyor
        coef = self.pinv @ patch[self.mask]
        patch = patch - (coef[0] * self.basis[0] + coef[1] * self.basis[1]
                         + coef[2] * self.basis[2])
        pupil = self.mask * np.exp(1j * patch)
        # sıfır-dolgulu FFT -> odak düzlemi yeğinliği
        pad = self.PAD_FACTOR * p
        big = np.zeros((pad, pad), dtype=np.complex128)
        big[:p, :p] = pupil
        psf = np.abs(np.fft.fftshift(np.fft.fft2(big))) ** 2
        return psf, pad

    def psf_at(self, X: float, Y: float, oversample: int = 1) -> np.ndarray:
        """(X, Y) [m] merkezli açıklık için kısa-poz PSF (toplamı 1).

        Eski kip (px_per_rad None): psf_px x psf_px, bir eleman = bir PSF ornegi.
        Yeni kip: goruntu izgarasinda, eleman boyu 1/oversample piksel (oversample>1
        yalniz dogrulama icindir; filtrelemede 1 kullanilir)."""
        psf, pad = self._focal_psf(X, Y)
        if self.sample_px is not None:
            return self._resampled(psf, pad, max(1, int(oversample)))
        h = self.psf_px // 2
        cpix = pad // 2
        # Ağırlık merkezi TAM PSF üzerinde (kırpmadan ÖNCE) hesaplanır: kırpma,
        # geniş kanatlı PSF'lerde (büyük D/r0) merkezi yanlı kestirir.
        s_full = psf.sum()
        if s_full > 0:
            yyf, xxf = np.mgrid[0:pad, 0:pad]
            cy_f = float((psf * yyf).sum()) / s_full - cpix
            cx_f = float((psf * xxf).sum()) / s_full - cpix
            sy, sx = int(round(cy_f)), int(round(cx_f))
            cpy, cpx = cpix + sy, cpix + sx
        else:
            cpy = cpx = cpix
        out = psf[max(cpy - h, 0):cpy + h + 1, max(cpx - h, 0):cpx + h + 1]
        if out.shape != (self.psf_px, self.psf_px):
            out = psf[cpix - h:cpix + h + 1, cpix - h:cpix + h + 1]
        s = out.sum()
        if s <= 0:
            k = np.zeros((self.psf_px, self.psf_px)); k[h, h] = 1.0
            return k
        out = out / s
        # BİRİNCİ MOMENT SIFIRLAMA (kritik): görüntü kayması = dalga cephesi EĞİMİdir ve
        # geometrik warp tarafından ZATEN uygulanır. Pupil yamasından düzlem eğim
        # çıkarılsa bile ayrık PSF'in kalan ağırlık merkezi sıfırdan farklı olur; bu
        # kalıntı, alan boyunca AYNI olduğu için tüm özniteliklere ortak (mükemmel
        # korelasyonlu) bir kayma ekler ve E2'de ölçülen korelasyonu bozar
        # (ampirik: A > 1, üs ~ 0.2). Bu yüzden PSF ağırlık merkezine göre yeniden
        # merkezlenir -> bulanıklık kanalı SAF radyometriktir, geometrisi yoktur.
        import cv2 as _cv2
        yy, xx = np.mgrid[0:self.psf_px, 0:self.psf_px]
        for _ in range(3):                      # kırpma sonrası kalıntıyı yinele
            cy = float((out * yy).sum()) - h
            cx = float((out * xx).sum()) - h
            if abs(cx) < 1e-3 and abs(cy) < 1e-3:
                break
            M = np.array([[1.0, 0.0, -cx], [0.0, 1.0, -cy]], dtype=np.float32)
            out = _cv2.warpAffine(out.astype(np.float32), M,
                                  (self.psf_px, self.psf_px),
                                  flags=_cv2.INTER_LINEAR,
                                  borderMode=_cv2.BORDER_CONSTANT, borderValue=0.0)
            s2 = out.sum()
            if s2 > 0:
                out = out / s2
        return out.astype(np.float64)

    def _resampled(self, psf: np.ndarray, pad: int, oversample: int) -> np.ndarray:
        """Odak duzlemi PSF'ini goruntu izgarasina enerji korunumlu yeniden ornekle.

        Her PSF ornegi, genisligi sample_px olan duzgun bir kutu kabul edilir ve
        kutlesi cikis piksellerine kutu-ortusme uzunlugu oraninda dagitilir
        (ayrilabilir: K = R_y P R_x^T). Bir ornegin agirliklari toplami 1'dir; yani
        cikis penceresine dusen enerji birebir korunur, kirpma kaybi sonda
        normalizasyonla giderilir. Izgara, PSF agirlik merkezi cekirdek merkezine
        dusecek bicimde kurulur (alt-ornek hassasiyetinde)."""
        s_full = float(psf.sum())
        os_ = int(oversample)
        w = 1.0 / os_
        n = _odd(int(np.ceil(pad * self.sample_px * os_)))
        if os_ == 1 and self.max_px is not None:
            n = min(n, _odd(self.max_px))
        h = n // 2
        if s_full <= 0:
            k = np.zeros((n, n)); k[h, h] = 1.0
            return k
        P = psf / s_full
        idx = np.arange(pad, dtype=np.float64)
        cy = float((P.sum(axis=1) * idx).sum())
        cx = float((P.sum(axis=0) * idx).sum())
        Ry = _box_rebin_matrix(pad, cy, self.sample_px, n, w)
        Rx = _box_rebin_matrix(pad, cx, self.sample_px, n, w)
        out = Ry @ P @ Rx.T
        s = out.sum()
        if s <= 0:
            k = np.zeros((n, n)); k[h, h] = 1.0
            return k
        out = out / s
        import cv2 as _cv2
        yy, xx = np.mgrid[0:n, 0:n]
        for _ in range(3):
            cyo = float((out * yy).sum()) - h
            cxo = float((out * xx).sum()) - h
            if abs(cxo) < 1e-3 and abs(cyo) < 1e-3:
                break
            M = np.array([[1.0, 0.0, -cxo], [0.0, 1.0, -cyo]], dtype=np.float32)
            out = _cv2.warpAffine(out.astype(np.float32), M, (n, n),
                                  flags=_cv2.INTER_LINEAR,
                                  borderMode=_cv2.BORDER_CONSTANT, borderValue=0.0)
            s2 = out.sum()
            if s2 > 0:
                out = out / s2
        return out.astype(np.float64)


def _box_rebin_matrix(n_in: int, c_in: float, s_in: float, n_out: int,
                      w_out: float) -> np.ndarray:
    """(n_out x n_in) kutu-ortusme yeniden ornekleme matrisi.

    Giris ornegi i, [u_i - s_in/2, u_i + s_in/2] araligini kaplar, u_i = (i - c_in) s_in;
    cikis pikseli j, merkezi (j - n_out//2) w_out ve genisligi w_out olan araliktir.
    Agirlik = ortusme uzunlugu / s_in; tamamen pencere icindeki bir ornegin
    agirliklari toplami 1'dir (kutle korunumu)."""
    u = (np.arange(n_in, dtype=np.float64) - c_in) * s_in
    lo, hi = u - 0.5 * s_in, u + 0.5 * s_in
    c = (np.arange(n_out, dtype=np.float64) - n_out // 2) * w_out
    L, R = c - 0.5 * w_out, c + 0.5 * w_out
    ov = np.minimum(hi[None, :], R[:, None]) - np.maximum(lo[None, :], L[:, None])
    return np.clip(ov, 0.0, None) / s_in


def _odd(n: int) -> int:
    n = max(int(n), 1)
    return n if n % 2 == 1 else n + 1


# =============================================================================
#  Sintilasyon alanı — log-normal yeğinlik dalgalanması
# =============================================================================
class ScintillationField:
    """Fresnel ölçeğinde korelasyonlu log-normal yeğinlik modülasyonu.

    I' = I * exp(2*chi - 2*var(chi)) ,  var(2 chi) = ln(1 + sigma_I^2)
    (ortalama yeğinliği koruyacak biçimde normalize edilir).
    """

    def __init__(self, N: int, delta: float, corr_len: float, sigma_I2: float,
                 rng: np.random.Generator):
        self.N, self.delta = N, delta
        if sigma_I2 <= 0:
            self.chi = None
            return
        f = np.fft.fftfreq(N, d=delta)
        FX, FY = np.meshgrid(f, f, indexing="xy")
        F2 = FX ** 2 + FY ** 2
        # Gauss korelasyonlu alan: PSD ~ exp(-(pi corr_len f)^2)
        psd = np.exp(-((np.pi * corr_len) ** 2) * F2)
        w = rng.standard_normal((N, N))
        g = np.fft.ifft2(np.fft.fft2(w) * np.sqrt(psd)).real
        g /= g.std()
        sigma_ln = np.sqrt(np.log1p(sigma_I2))
        self.chi = (0.5 * sigma_ln) * g          # 2*chi ~ N(0, sigma_ln^2)

    def sample(self, X: np.ndarray, Y: np.ndarray) -> np.ndarray:
        """Çarpımsal yeğinlik kazancı (ortalaması 1)."""
        if self.chi is None:
            return np.ones_like(X)
        gx = X / self.delta + self.N / 2.0
        gy = Y / self.delta + self.N / 2.0
        c = _bilinear_periodic(self.chi, gx, gy)
        v = float(np.var(self.chi)) * 4.0
        return np.exp(2.0 * c - 0.5 * v)
