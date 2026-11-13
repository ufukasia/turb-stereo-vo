"""
Derinliğe bağlı, yol-integralli, çok-ekranlı türbülans modeli (AirSim sekansları için).

Fizik (TASARIM.md §5.2; Hardie vd. 2017 OE 56(7):071502, Chimitt & Chan 2020 OE 59(8):083101):

  * Piksel başına görüş hattı: kamera merkezinden, pikselin render derinliğinden bulunan
    nesne noktasına. Hattın zemine göre yüksekliği h(s) arazi modelinden (DTM) alınır.
  * Cn²(h) = C0 · (max(h, h_min) / 2 m)^(-4/3)   (Wyngaard 1971, gündüz yüzey tabakası).
  * Küresel dalga, z NESNEDEN ölçülür:
        r0^(-5/3) = 0.423 k² ∫_0^L Cn²(s) (z/L)^(5/3) ds          (z = L − s, s kameradan)
    Tek eksen G-tilt varyansı (toplam):  σ²_tot = 0.170 (D/r0)^(5/3) (λ/D)²
  * Uzaysal yapı (anizoplanatizm): kameraya dik, uzaklıkla geometrik aralıklı M faz ekranı.
    Ekran m'nin dilimindeki türbülans, pikselin kendi hattı üzerindeki ∫_dilim Cn² ds ile
    ölçeklenen ORTAK bir Kolmogorov/von Kármán gerçeklemesiyle temsil edilir. Ekranın
    açıklık-ortalamalı gradyanı açıyı (z/L) ağırlığıyla saptırır (perde etkisi).
    Ekran alanı açısal çözünürlükte (varsayılan 4 px) bant-sınırlıdır: bu "çarpıtma" (warp)
    bileşenidir. Bant dışında kalan ince ölçekli tilt varyansı (σ²_tot − σ²_düşük) bir
    karenin içinde ayrışamaz; piksel ölçeğinde BULANIKLIK olarak eklenir.
  * Kısa-poz (tilt çıkarılmış) türbülans bulanıklığı: Fried (1966) yakın-alan kısa-poz MTF'i
    × kırınım MTF'inden sayısal PSF; ek Gauss varyansı FWHM eşlemesiyle
    σ²_ek = (FWHM_top² − FWHM_kırınım²) / (8 ln 2). Kırınımın kendisi EKLENMEZ: render, kameranın
    türbülanssız görüntüsü sayılır, yalnız türbülansın getirdiği fazlalık eklenir
    (C0 = 0 ⇒ çıktı = girdi, bit düzeyinde).
  * Zaman: ekranlar havaya bağlıdır (donmuş akış). Kamera dönmesi/yanal ötelemesi ve rüzgâr
    ekran koordinatını kaydırır; eksenel ilerleme dilim kalınlığı boyunca iki bağımsız
    gerçekleme arasında enerji-koruyan karışımla içeriği yeniler.
  * Stereo: cam1 hattı ekranları B(1 − d/L) kadar kaymış konumda keser; cam0–cam1 tilt
    korelasyonu bu geometriden kendiliğinden doğar.

Sınırlılıklar (makalede yazılacak): sintilasyon yok; dilim içinde Cn² yatay olarak yalnız
piksel hattı üzerinden değişir (ortak desen genliği modüle edilir); 4 px altı tilt yapısı
Gauss bulanıklığına indirgenir; DTM kapsamı dışında arazi kenar değeriyle uzatılır.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict

import numpy as np

from . import screen as _screen

C_GTILT = 0.170


# =============================================================================
#  Parametreler
# =============================================================================
@dataclass
class DerinlikParam:
    C0: float = 1e-13            # Cn²(2 m) [m^-2/3]
    h_min: float = 0.5           # [m] bunun altında Cn² sabit
    D: float = 4.8e-3            # [m] açıklık (tele: f=13.4 mm, F/2.8)
    lam: float = 5.5e-7          # [m]
    L0: float = 5.0              # [m] dış ölçek (yer yakını, ~yükseklik)
    M: int = 12                  # ekran sayısı
    d_min: float = 0.5           # [m] ilk dilim alt sınırı
    d_max: float = 3000.0        # [m] son dilim üst sınırı
    alt_ornek: int = 4           # dilim başına hat örneği
    warp_px: float = 4.0         # çarpıtma bileşeninin açısal çözünürlüğü [px]
    N_ekran: int = 1024
    alt_harmonik: int = 8
    ruzgar: float = 2.0          # [m/s] yatay rüzgâr
    ruzgar_yon_deg: float = 30.0 # NED x'ten (kuzey) saat yönünde
    L_gok: float = 2000.0        # [m] isabetsiz (gök) hatların uzunluğu
    tohum: int = 0

    def anahtar(self) -> str:
        return f"C0{self.C0:.0e}_s{self.tohum}".replace("+", "")


# =============================================================================
#  Arazi
# =============================================================================
class Arazi:
    """Zemin yüksekliği H(x, y) [m, yukarı +] — DTM ızgarası ya da düz düzlem."""

    def __init__(self, dtm_npz: str | None = None, z_zemin_ned: float = 0.0):
        if dtm_npz:
            d = np.load(dtm_npz)
            self.H = (-d["dtm"]).astype(np.float32)          # NED z -> yükseklik
            self.x0, self.y0, self.res = float(d["x0"]), float(d["y0"]), float(d["res"])
        else:
            self.H = None
            self.h0 = -float(z_zemin_ned)

    def yukseklik(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        if self.H is None:
            return np.full(np.shape(x), self.h0, np.float32)
        H = self.H
        fi = np.clip((x - self.x0) / self.res, 0, H.shape[0] - 1.001)
        fj = np.clip((y - self.y0) / self.res, 0, H.shape[1] - 1.001)
        i0 = fi.astype(np.int32); j0 = fj.astype(np.int32)
        di = fi - i0; dj = fj - j0
        return (H[i0, j0] * (1 - di) * (1 - dj) + H[i0 + 1, j0] * di * (1 - dj)
                + H[i0, j0 + 1] * (1 - di) * dj + H[i0 + 1, j0 + 1] * di * dj)


# =============================================================================
#  Kısa-poz türbülans bulanıklığı (Gauss eşdeğeri, kırınım fazlası)
# =============================================================================
def _fwhm(mtf, nu, r):
    """Dairesel simetrik MTF(ν) -> PSF(r) (Hankel) -> FWHM; birimler λ/D."""
    from scipy.special import j0
    psf = np.trapz(mtf[None, :] * j0(2 * np.pi * nu[None, :] * r[:, None]) * nu[None, :], nu, axis=1)
    psf = psf / psf[0]
    i = int(np.argmax(psf < 0.5))
    r1, r2, p1, p2 = r[i - 1], r[i], psf[i - 1], psf[i]
    return 2 * (r1 + (0.5 - p1) * (r2 - r1) / (p2 - p1))


def _se_tablo():
    """Kısa-poz türbülans bulanıklığının, kırınıma EK Gauss varyansı [(λ/D)²]:
    FWHM eşlemesi:  σ²_ek = (FWHM_top² − FWHM_kırınım²) / (8 ln 2).
    Toplam MTF = kırınım MTF × Fried yakın-alan kısa-poz MTF
                 exp{−3.44 (D/r0)^{5/3} ν^{5/3} [1 − ν^{1/3}]},  ν = λf/D."""
    x = np.logspace(-3, math.log10(30.0), 120)
    nu = np.linspace(0.0, 1.0, 4001)
    mtf_d = (2 / np.pi) * (np.arccos(nu) - nu * np.sqrt(np.clip(1 - nu * nu, 0, 1)))
    r = np.linspace(0, 40.0, 16001)
    f_d = _fwhm(mtf_d, nu, r[:2000])
    s2 = []
    for xi in x:
        mt = mtf_d * np.exp(-3.44 * xi ** (5 / 3) * nu ** (5 / 3) * (1 - nu ** (1 / 3)))
        rr = r if xi > 3 else r[:2000]
        f = _fwhm(mt, nu, rr)
        s2.append(max(f * f - f_d * f_d, 0.0) / (8 * math.log(2)))
    return x, np.array(s2), f_d


_SE_X, _SE_S2, _FWHM_KIRINIM = _se_tablo()


def se_bulaniklik_var_rad2(D_over_r0: np.ndarray, lam: float, D: float) -> np.ndarray:
    """Kısa-poz türbülans MTF'inin Gauss eşdeğer varyansı [rad²] (kırınıma EK)."""
    q = np.asarray(D_over_r0, np.float64)
    s2 = np.interp(np.log(np.maximum(q, 1e-6)), np.log(_SE_X), _SE_S2, left=0.0)
    s2 = np.where(q > _SE_X[-1], _SE_S2[-1] * (q / _SE_X[-1]) ** (5 / 3), s2)
    return s2 * (lam / D) ** 2


# =============================================================================
#  Model
# =============================================================================
def _airy(f, D):
    from scipy.special import j1
    x = np.pi * D * np.asarray(f, np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        A = np.where(x > 1e-12, 2 * j1(x) / np.where(x > 1e-12, x, 1), 1.0)
    return A


class DerinlikModeli:
    def __init__(self, p: DerinlikParam, fx: float, arazi: Arazi):
        self.p, self.fx, self.arazi = p, float(fx), arazi
        self.k = 2 * np.pi / p.lam
        b = np.geomspace(p.d_min, p.d_max, p.M + 1)
        b[0] = 0.0
        self.sinir = b                                        # dilim sınırları [m]
        self.d_ekran = np.sqrt(np.maximum(b[:-1], p.d_min * 0.5) * b[1:])
        self.dlp = p.warp_px * self.d_ekran / self.fx         # warp bant sınırı (açıklık eşdeğeri) [m]
        self.delta = self.dlp / 2.0                           # ekran ızgara adımı [m] (bant sınırının yarısı)
        self._ekranlar = None
        self._var_dusuk = None
        # Testler için: Cn²/C0 profilini (h, eksenel uzaklık a, Z) -> göreli Cn² ile değiştirir.
        self.cn_fn = None

    # ------------------------------------------------------------------ ekranlar
    def _ekran_kur(self):
        """Her dilim için 2 bağımsız gerçekleme (A, B) × (αx, αy), birim güç (r0 = 1 m)."""
        p = self.p
        ss = np.random.SeedSequence([p.tohum, 0xD3E7])
        ekr, var = [], []
        for m in range(p.M):
            dm, Dlp = self.delta[m], self.dlp[m]
            alanlar = []
            for g in range(2):
                rng = np.random.default_rng(ss.spawn(1)[0])
                real = _screen.KolmogorovRealization(p.N_ekran, dm, 1.0, rng,
                                                     n_subharmonics=p.alt_harmonik, L0=p.L0)
                F = np.hypot(real.FX, real.FY)
                A = _airy(F, p.D) * _airy(F, Dlp)
                Hx = (2j * np.pi * real.FX) * A / self.k
                Hy = (2j * np.pi * real.FY) * A / self.k
                hx = lambda a, b_, D_=p.D, L_=Dlp: (2j * np.pi * a) * _airy(np.hypot(a, b_), D_) * _airy(np.hypot(a, b_), L_) / self.k
                hy = lambda a, b_, D_=p.D, L_=Dlp: (2j * np.pi * b_) * _airy(np.hypot(a, b_), D_) * _airy(np.hypot(a, b_), L_) / self.k
                alanlar.append((real.filtered(Hx, hx).astype(np.float32),
                                real.filtered(Hy, hy).astype(np.float32)))
            ekr.append(alanlar)
            var.append(self._var_dusuk_izgara(real, Dlp))
        self._ekranlar, self._var_dusuk = ekr, np.array(var)

    def _var_dusuk_izgara(self, real, Dlp: float) -> float:
        """Gerçeklenen sürecin BEKLENEN tek eksen varyansı [rad²], birim güç (r0 = 1 m):
        FFT ızgarası modları + alt-harmonik modlar üzerinden Σ Φ |H|² Δf² (ekranın temsil
        ettiği bant tam olarak; Nyquist üstü kısım σ²_tot − σ²_düşük ile bulanıklığa gider)."""
        p = self.p
        N, dl = p.N_ekran, real.delta
        df = 1.0 / (N * dl)
        F = np.hypot(real.FX, real.FY)
        psd = _screen.kolmogorov_psd(F, 1.0, p.L0); psd[N // 2, N // 2] = 0.0
        A2 = (_airy(F, p.D) * _airy(F, Dlp)) ** 2
        g = psd * A2 * (2 * np.pi / self.k) ** 2 * 0.5 * (real.FX ** 2 + real.FY ** 2)
        v = float(g.sum() * df * df)
        Ds = N * dl
        for q in range(1, p.alt_harmonik + 1):
            dfs = 1.0 / (3.0 ** q * Ds)
            fx_ = np.array([-1.0, 0.0, 1.0]) * dfs
            SX, SY = np.meshgrid(fx_, fx_); SF = np.hypot(SX, SY)
            sp = _screen.kolmogorov_psd(SF, 1.0, p.L0); sp[1, 1] = 0.0
            v += float(np.sum(sp * (_airy(SF, p.D) * _airy(SF, Dlp)) ** 2
                              * (2 * np.pi / self.k) ** 2 * 0.5 * (SX ** 2 + SY ** 2)) * dfs * dfs)
        return v

    def _var_dusuk_analitik(self, Dlp: float) -> float:
        """Birim güç (r0 = 1 m) bant-sınırlı ekranın tek eksen tilt varyansı [rad²]:
        ∫ Φ(f) (2π f_x / k)² A(fD)² A(f Dlp)² d²f  (von Kármán Φ, açısal ortalama f_x² = f²/2)."""
        p = self.p
        f = np.geomspace(1e-4, 50.0 / min(p.D, Dlp), 20000)
        psd = _screen.kolmogorov_psd(f, 1.0, p.L0)
        g = psd * (2 * np.pi / self.k) ** 2 * (f ** 2 / 2) * (_airy(f, p.D) * _airy(f, Dlp)) ** 2
        return float(np.trapz(g * 2 * np.pi * f, f))

    # ------------------------------------------------------------------ hat integrali
    def hat_integralleri(self, Z: np.ndarray, egim_x: np.ndarray, egim_y: np.ndarray,
                         R_wc: np.ndarray, t_wc: np.ndarray, ofset_x: float = 0.0):
        """Piksel hatları boyunca Cn² integralleri.

        Z: (H, W) düzlemsel derinlik [m] (inf = gök); egim_x/y: (H, W) (u−cx)/fx, (v−cy)/fy;
        R_wc, t_wc: kamera (görüntü çerçevesi) -> dünya (NED); ofset_x: kamera merkezinin
        görüntü-x yönündeki kayması (cam1 için +B).
        Dönüş: W_m (M, H, W) = 0.423 k² ∫_dilim Cn² ds  [m^-5/3],
               zL_m (M, H, W) dilimin etkin ekranında z/L,  d_eff (M, H, W),
               r0_kure^-5/3 (H, W).
        """
        p = self.p
        Zt = np.where(np.isfinite(Z), Z, p.L_gok).astype(np.float64)
        sec = np.sqrt(1 + egim_x ** 2 + egim_y ** 2)          # ds / d(eksenel)
        H_, W_ = Z.shape
        M, K = p.M, p.alt_ornek
        Wm = np.zeros((M, H_, W_)); zL = np.zeros((M, H_, W_)); deff = np.zeros((M, H_, W_))
        r0s = np.zeros((H_, W_))
        # dünya: P(a) = t + R [ofset + a*egim_x, a*egim_y, a]   (a = eksenel uzaklık)
        R = np.asarray(R_wc, np.float64); t = np.asarray(t_wc, np.float64)
        dirw = (R[:, 0:1, None] * egim_x[None] + R[:, 1:2, None] * egim_y[None] + R[:, 2:3, None])  # (3,H,W)
        orig = t[:, None, None] + R[:, 0:1, None] * ofset_x
        c0 = 0.423 * self.k ** 2 * p.C0
        for m in range(M):
            lo = self.sinir[m]; hi = np.minimum(self.sinir[m + 1], Zt)
            gecerli = hi > lo
            if not gecerli.any():
                continue
            hi_ = np.where(gecerli, hi, lo)
            # dilim içinde K orta-nokta örneği
            fr = (np.arange(K) + 0.5) / K
            acc = np.zeros((H_, W_)); acc_s = np.zeros((H_, W_))
            for f_ in fr:
                a = lo + f_ * (hi_ - lo)
                Pw = orig + dirw * a[None]
                h = -Pw[2] - self.arazi.yukseklik(Pw[0], Pw[1])
                cn = (self.cn_fn(h, a, Zt) if self.cn_fn is not None
                      else (np.maximum(h, p.h_min) / 2.0) ** (-4 / 3))
                ds = (hi_ - lo) / K * sec
                acc += cn * ds
                acc_s += cn * ds * np.clip(1 - a / Zt, 0, 1) ** (5 / 3)
            Wm[m] = np.where(gecerli, c0 * acc, 0.0)
            r0s += np.where(gecerli, c0 * acc_s, 0.0)
            dm = 0.5 * (lo + hi_)
            deff[m] = dm
            zL[m] = np.where(gecerli, np.clip(1 - dm / Zt, 0, 1), 0.0)
        return Wm, zL, deff, r0s

    # ------------------------------------------------------------------ kare
    def kare(self, Z, egim_x, egim_y, R_wc, t_wc, ofset_m, karisim_m, ofset_x=0.0):
        """Tek kamera, tek kare: (disp_x, disp_y) [px], blur_sigma [px], tanı sözlüğü.
        ofset_m: (M, 2) ekran kayması [m]; karisim_m: (M,) A/B karışım açısı [rad]."""
        p = self.p
        if self._ekranlar is None:
            self._ekran_kur()
        Wm, zL, deff, r0s = self.hat_integralleri(Z, egim_x, egim_y, R_wc, t_wc, ofset_x)
        ax = np.zeros(Z.shape); ay = np.zeros(Z.shape); var_lo = np.zeros(Z.shape)
        for m in range(p.M):
            if not Wm[m].any():
                continue
            X = ofset_x + deff[m] * egim_x + ofset_m[m, 0]
            Y = deff[m] * egim_y + ofset_m[m, 1]
            gx = X / self.delta[m] + p.N_ekran / 2.0
            gy = Y / self.delta[m] + p.N_ekran / 2.0
            c, s = math.cos(karisim_m[m]), math.sin(karisim_m[m])
            (Ax, Ay), (Bx, By) = self._ekranlar[m]
            sx = c * _bilin(Ax, gx, gy) + s * _bilin(Bx, gx, gy)
            sy = c * _bilin(Ay, gx, gy) + s * _bilin(By, gx, gy)
            amp = zL[m] * np.sqrt(Wm[m])
            ax += amp * sx; ay += amp * sy
            var_lo += (amp ** 2) * self._var_dusuk[m]
        var_tot = C_GTILT * p.lam ** 2 * p.D ** (-1 / 3) * r0s       # r0s = r0_kure^(-5/3)
        D_r0 = p.D * r0s ** (3 / 5)
        var_hi = np.maximum(var_tot - var_lo, 0.0)
        var_se = se_bulaniklik_var_rad2(D_r0, p.lam, p.D)
        blur = np.sqrt(var_hi + var_se) * self.fx
        tani = {"sigma_tot_px": np.sqrt(var_tot) * self.fx, "sigma_lo_px": np.sqrt(var_lo) * self.fx,
                "D_r0": D_r0}
        return ax * self.fx, ay * self.fx, blur, tani

    def info(self) -> dict:
        return {"param": asdict(self.p), "fx": self.fx, "dilim_sinir_m": self.sinir.tolist(),
                "ekran_uzaklik_m": self.d_ekran.tolist(), "ekran_adim_m": self.delta.tolist()}


def _bilin(F, gx, gy):
    """Periyodik bilinear örnekleme; satır = y, sütun = x (screen.py ile aynı)."""
    N = F.shape[0]
    x0 = np.floor(gx).astype(np.int64); y0 = np.floor(gy).astype(np.int64)
    fx = gx - x0; fy = gy - y0
    x0 %= N; y0 %= N; x1 = (x0 + 1) % N; y1 = (y0 + 1) % N
    return (F[y0, x0] * (1 - fx) * (1 - fy) + F[y0, x1] * fx * (1 - fy)
            + F[y1, x0] * (1 - fx) * fy + F[y1, x1] * fx * fy)


# =============================================================================
#  Görüntüye uygulama
# =============================================================================
def uygula(img: np.ndarray, disp_x: np.ndarray, disp_y: np.ndarray, blur: np.ndarray,
           n_seviye: int = 10) -> np.ndarray:
    """Düşük çözünürlüklü alanları görüntü boyuna büyütüp çarpıtma + uzaysal değişken Gauss
    bulanıklığı uygular. img uint8 (H, W)."""
    import cv2
    H, W = img.shape[:2]
    dx = cv2.resize(disp_x.astype(np.float32), (W, H), interpolation=cv2.INTER_CUBIC)
    dy = cv2.resize(disp_y.astype(np.float32), (W, H), interpolation=cv2.INTER_CUBIC)
    sb = cv2.resize(blur.astype(np.float32), (W, H), interpolation=cv2.INTER_LINEAR)
    u, v = np.meshgrid(np.arange(W, dtype=np.float32), np.arange(H, dtype=np.float32))
    w = cv2.remap(img.astype(np.float32), u - dx, v - dy, cv2.INTER_CUBIC, borderMode=cv2.BORDER_REFLECT)
    smax = float(sb.max())
    if smax < 0.05:
        out = w
    else:
        lv = np.linspace(0.0, smax, n_seviye)
        stk = [w] + [cv2.GaussianBlur(w, (0, 0), s) for s in lv[1:]]
        pos = np.clip(sb / (lv[1] - lv[0]), 0, n_seviye - 1 - 1e-6)
        i0 = pos.astype(np.int32); fr = pos - i0
        out = np.zeros_like(w)
        for k in range(n_seviye - 1):
            m = i0 == k
            if m.any():
                out[m] = stk[k][m] * (1 - fr[m]) + stk[k + 1][m] * fr[m]
    return np.clip(np.rint(out), 0, 255).astype(np.uint8)
