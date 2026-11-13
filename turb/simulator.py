"""
Türbülans simülatörü — kare-tutarlı, tohumlu, tekrarlanabilir.

Bir `TurbulenceSimulator` örneği TEK bir atmosfer gerçeklemesini temsil eder
(tohum + D/r0 + theta_c + tau0). Kurulumda ekranlar bir kez üretilir; her kare
yalnızca örnekleme + warp maliyetindedir.

BEYAN EDİLEN OPTİK GEOMETRİ (makalede açıkça yazılacak)
------------------------------------------------------
EuRoC/TUM VI görüntüleri geniş görüş açılı, kısa menzilli iç mekân kayıtlarıdır;
türbülans bu geometride fiziksel olarak ihmal edilebilirdir (rapor R3). Bu yüzden
görüntüler, BEYAN EDİLEN bir uzun menzilli tele-objektif geometrisinden görülüyormuş
gibi ele alınır:

    D        = 0.10 m      açıklık
    lambda   = 550 nm
    f/p      = 170648 px/rad   (f = 1.0 m, p = 5.86 um)  -> IFOV = 5.86 urad/px
    FOV(752 px) = 4.41 mrad = 0.25 deg

Bu geometride D/r0, tilt jitterini sigma_px = sqrt(0.170) (D/r0)^{5/6} (lambda/D) (f/p)
olarak belirler: D/r0=1 -> 0.39 px, D/r0=5 -> 1.48 px, D/r0=12 -> 3.05 px.
Sahne içeriği EuRoC'tan, bozulma fiziği bu geometriden gelir. İddialar makalede
"denetimli bozulma benchmark'ı" (EuRoC/TUM VI) ve "fiziksel olarak makul rejim"
(TartanAir dış mekân) olarak AYRILARAK yazılır.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field, asdict

import cv2
import numpy as np

from . import optics
from .fields import TiltField, SpeckleBlurField, ScintillationField


# =============================================================================
@dataclass
class TurbParams:
    """Türbülans denetim parametreleri (config.yaml `turbulence:` bölümü)."""
    enabled: bool = False
    seed: int = 0

    # --- ana denetim eksenleri (rapor §5.3) ---
    D_over_r0: float = 0.0             # bozulma şiddeti
    theta_c_over_fov: float = 0.15     # uzaysal korelasyon ölçeği (E1-c)
    tau0_over_dt: float = 1.0          # zamansal rejim (E1-b)

    # --- beyan edilen optik geometri ---
    aperture_D: float = 0.10           # [m]
    wavelength: float = 550e-9         # [m]
    px_per_rad: float = 170648.0       # f/p
    frame_dt: float = 0.05             # [s] (EuRoC 20 Hz)
    stereo_baseline: float = 0.11      # [m] — cam1 ayak izi ötelemesi

    # --- bozulma kanalları ---
    tilt: bool = True
    blur: str = "speckle_global"       # none | gaussian | speckle_global | speckle_grid
    blur_grid: tuple = (3, 4)          # speckle_grid için (satır, sütun)
    scint_sigma_I2: float = 0.0        # sintilasyon indeksi (0 = kapalı)

    # --- sayısal ayarlar ---
    warp_grid_step: int = 8            # yer değiştirme alanı kaba ızgarada hesaplanır
    psf_max_px: int = 41               # PSF çekirdek üst sınırı [px]
    # PSF ölçeği. False (varsayılan): kısa-poz PSF, örnek boyu lam/(4D)*px_per_rad piksel
    # olacak şekilde görüntü ızgarasına enerji korunumlu yeniden örneklenir (kırınım
    # FWHM ≈ 1.03 lam/D*px_per_rad ≈ 0.97 px). True: 2026-11-07 öncesi davranış — PSF
    # örnekleri doğrudan piksel sayılır, bulanıklık ~4.3 kat geniş (E1a-E9 kanıtları
    # bu kiple üretildi; yeniden üretim için korunur).
    psf_legacy_scale: bool = False
    screen_N: int = 2048
    screen_px_per_D: int = 8           # kaba ekran çözünürlüğü (delta = D/8)
    subharmonics: int = 8
    fine_N: int = 1024
    fine_px_per_D: int = 32
    tilt_calibration: float = 1.085      # ölçülen/teorik sapmayı kapatan katsayı
    wind_angle_deg: float = 30.0       # rüzgâr yönü (ekran düzleminde)

    def key(self) -> str:
        """Bu parametre kümesinin kısa özeti (dosya adı / önbellek anahtarı)."""
        s = (f"dr{self.D_over_r0:g}_tc{self.theta_c_over_fov:g}"
             f"_t0{self.tau0_over_dt:g}_s{self.seed}_{self.blur}"
             f"_sc{self.scint_sigma_I2:g}")
        return s

    def hash(self) -> str:
        return hashlib.sha1(repr(sorted(asdict(self).items())).encode()).hexdigest()[:10]


# =============================================================================
class TurbulenceSimulator:
    """Tek atmosfer gerçeklemesi; kare kare bozulma uygular ve GT alanını verir."""

    def __init__(self, p: TurbParams, image_size=(752, 480)):
        self.p = p
        self.W, self.H = int(image_size[0]), int(image_size[1])
        self.active = bool(p.enabled) and float(p.D_over_r0) > 0.0
        if not self.active:
            return

        D, lam = p.aperture_D, p.wavelength
        self.r0 = optics.r0_from_D_over_r0(p.D_over_r0, D)

        # --- katman geometrisi: theta_c hedefinden z ---
        self.fov_x = self.W / p.px_per_rad                       # [rad]
        self.theta_c = p.theta_c_over_fov * self.fov_x
        self.z_layer = optics.layer_distance_for_theta_c(D, self.theta_c)
        self.theta0 = optics.isoplanatic_angle(self.r0, self.z_layer)

        # --- donmuş akış hızı: tau0 = 0.314 r0 / v ---
        self.tau0 = p.tau0_over_dt * p.frame_dt
        self.v_wind = optics.wind_for_tau0(self.r0, self.tau0)
        a = np.deg2rad(p.wind_angle_deg)
        self.wind_vec = np.array([np.cos(a), np.sin(a)]) * self.v_wind

        # --- ekranlar ---
        delta = D / p.screen_px_per_D
        rng = np.random.default_rng(_mix_seed(p.seed, 0xA1))
        self.tiltfield = TiltField(p.screen_N, delta, self.r0, D, lam, rng,
                                   n_subharmonics=p.subharmonics,
                                   calibration=p.tilt_calibration)
        self.screen_extent = p.screen_N * delta

        self.blurfield = None
        if p.blur in ("speckle_global", "speckle_grid"):
            fdelta = D / p.fine_px_per_D
            rng_f = np.random.default_rng(_mix_seed(p.seed, 0xB2))
            psf_px = min(int(p.psf_max_px),
                         _odd(int(6 * max(1.0, p.D_over_r0)) + 5))
            if p.psf_legacy_scale:
                self.blurfield = SpeckleBlurField(p.fine_N, fdelta, self.r0, D, lam, rng_f,
                                                  pupil_px=p.fine_px_per_D,
                                                  psf_px=_odd(psf_px))
            else:
                self.blurfield = SpeckleBlurField(p.fine_N, fdelta, self.r0, D, lam, rng_f,
                                                  pupil_px=p.fine_px_per_D,
                                                  psf_px=_odd(psf_px),
                                                  px_per_rad=p.px_per_rad,
                                                  max_px=p.psf_max_px)
            self.fine_extent = p.fine_N * fdelta
        self.gauss_sigma_px = (optics.short_exposure_blur_sigma_rad(self.r0, lam, D)
                               * p.px_per_rad)

        self.scint = None
        if p.scint_sigma_I2 > 0:
            lF = optics.fresnel_scale(lam, self.z_layer)
            rng_s = np.random.default_rng(_mix_seed(p.seed, 0xC3))
            self.scint = ScintillationField(1024, D / 8.0, lF, p.scint_sigma_I2, rng_s)

        # --- piksel -> ayak izi taban ızgarası ---
        # Tilt alanının uzaysal korelasyon ölçeği theta_c * (f/p) pikseldir (en dar
        # ayarda ~37 px). Bu yüzden yer değiştirme alanı `warp_grid_step` adımlı KABA
        # ızgarada hesaplanıp bilinear yükseltilir: alan zaten bu ölçekte düzgün olduğu
        # için hata ihmal edilebilir, maliyet ~step^2 kat düşer. Seyrek GT sorguları
        # (`displacement_at`) her zaman TAM çözünürlükte ve interpolasyonsuzdur.
        step = max(1, int(p.warp_grid_step))
        self.step = step
        # Kaba ızgara konumları, cv2.resize'ın bilinear örnekleme kuralıyla TAM uyumlu
        # seçilir: dst(i) <- src((i+0.5)*W_c/W - 0.5). Böylece yükseltme, alanın doğru
        # piksel konumlarındaki bilinear interpolasyonudur (yarım piksel kayması yok).
        Wc = max(2, int(np.ceil(self.W / step)))
        Hc = max(2, int(np.ceil(self.H / step)))
        self._coarse_wh = (Wc, Hc)
        gu = (np.arange(Wc) + 0.5) * (self.W / Wc) - 0.5
        gv = (np.arange(Hc) + 0.5) * (self.H / Hc) - 0.5
        GU, GV = np.meshgrid(gu, gv)
        self._grid_shape = GU.shape
        self.X0 = self.z_layer * (GU - (self.W - 1) / 2.0) / p.px_per_rad
        self.Y0 = self.z_layer * (GV - (self.H - 1) / 2.0) / p.px_per_rad
        # tam çözünürlüklü piksel ızgaraları (remap ve harmanlama için)
        self.U, self.V = np.meshgrid(np.arange(self.W, dtype=np.float32),
                                     np.arange(self.H, dtype=np.float32))

        self.sigma_px_theory = optics.sigma_tilt_px(p.D_over_r0, lam, D, p.px_per_rad)

    # ------------------------------------------------------------------ bilgi
    def info(self) -> dict:
        if not self.active:
            return {"active": False}
        return {
            "active": True,
            "D_over_r0": self.p.D_over_r0,
            "r0_m": self.r0,
            "z_layer_m": self.z_layer,
            "theta_c_rad": self.theta_c,
            "theta_c_px": self.theta_c * self.p.px_per_rad,
            "theta0_rad": self.theta0,
            "fov_x_rad": self.fov_x,
            "tau0_s": self.tau0,
            "v_wind_ms": self.v_wind,
            "sigma_px_theory": self.sigma_px_theory,
            "sigma_px_field": self.tiltfield.tilt_std() * self.p.px_per_rad,
            "blur_sigma_px_gauss": self.gauss_sigma_px,
            "screen_extent_m": self.screen_extent,
        }

    # ------------------------------------------------------- yer değiştirme
    def _footprint(self, frame_idx: int, cam: int):
        t = frame_idx * self.p.frame_dt
        dx = self.wind_vec[0] * t + (self.p.stereo_baseline if cam == 1 else 0.0)
        dy = self.wind_vec[1] * t
        return self.X0 + dx, self.Y0 + dy

    def displacement(self, frame_idx: int, cam: int = 0, full: bool = True):
        """Kare için (du, dv) yer değiştirme alanı [piksel].

        Tanım: temiz görüntüde (u,v) konumundaki sahne noktası, bozuk görüntüde
        (u + du, v + dv) konumunda görünür (ileri yönlü tanım).
        full=False ise kaba ızgara çözünürlüğünde döner.
        """
        if not self.active or not self.p.tilt:
            z = np.zeros((self.H, self.W), dtype=np.float32)
            return z, z.copy()
        X, Y = self._footprint(frame_idx, cam)
        ax, ay = self.tiltfield.sample(X, Y)
        s = self.p.px_per_rad
        du = (ax * s).astype(np.float32)
        dv = (ay * s).astype(np.float32)
        if not full or self.step == 1:
            return du, dv
        du = cv2.resize(du, (self.W, self.H), interpolation=cv2.INTER_LINEAR)
        dv = cv2.resize(dv, (self.W, self.H), interpolation=cv2.INTER_LINEAR)
        return du, dv

    def displacement_at(self, frame_idx: int, cam: int, u: np.ndarray, v: np.ndarray):
        """Seyrek (u, v) TEMİZ konumlarında ileri yer değiştirme [piksel] (E2 GT'si)."""
        if not self.active or not self.p.tilt:
            z = np.zeros_like(np.asarray(u, dtype=np.float64))
            return z, z.copy()
        p = self.p
        t = frame_idx * p.frame_dt
        dx = self.wind_vec[0] * t + (p.stereo_baseline if cam == 1 else 0.0)
        dy = self.wind_vec[1] * t
        th_x = (np.asarray(u, dtype=np.float64) - (self.W - 1) / 2.0) / p.px_per_rad
        th_y = (np.asarray(v, dtype=np.float64) - (self.H - 1) / 2.0) / p.px_per_rad
        X = self.z_layer * th_x + dx
        Y = self.z_layer * th_y + dy
        ax, ay = self.tiltfield.sample(X, Y)
        return ax * p.px_per_rad, ay * p.px_per_rad

    # -------------------------------------------------------------- uygulama
    def apply(self, img: np.ndarray, frame_idx: int, cam: int = 0) -> np.ndarray:
        """Tek kareye bozulmayı uygula (gri seviye, uint8)."""
        if not self.active:
            return img
        out = img
        if self.p.tilt:
            du, dv = self.displacement(frame_idx, cam)
            # Geri-yönlü harita: out(u) = in(u - d)  ->  içerik +d kadar kayar.
            map_x = (self.U - du).astype(np.float32)
            map_y = (self.V - dv).astype(np.float32)
            out = cv2.remap(out, map_x, map_y, cv2.INTER_LINEAR,
                            borderMode=cv2.BORDER_REFLECT_101)
        out = self._blur(out, frame_idx, cam)
        if self.scint is not None:
            X, Y = self._footprint(frame_idx, cam)
            g = self.scint.sample(X, Y)
            out = np.clip(out.astype(np.float32) * g, 0, 255).astype(np.uint8)
        return out

    def _blur(self, img: np.ndarray, frame_idx: int, cam: int) -> np.ndarray:
        mode = self.p.blur
        if mode == "none" or not self.active:
            return img
        if mode == "gaussian":
            s = self.gauss_sigma_px
            if s < 0.3:
                return img
            k = _odd(int(6 * s) + 1)
            return cv2.GaussianBlur(img, (k, k), s)
        if self.blurfield is None:
            return img
        X, Y = self._footprint(frame_idx, cam)          # kaba ızgara boyutunda
        Hc, Wc = X.shape
        if mode == "speckle_global":
            psf = self.blurfield.psf_at(float(X[Hc // 2, Wc // 2]),
                                        float(Y[Hc // 2, Wc // 2]))
            return cv2.filter2D(img, -1, psf, borderType=cv2.BORDER_REFLECT_101)
        # speckle_grid: karo bazlı uzaysal değişken PSF + bilinear harmanlama
        gr, gc = self.p.blur_grid
        acc = np.zeros((self.H, self.W), dtype=np.float32)
        wsum = np.zeros((self.H, self.W), dtype=np.float32)
        ys = np.linspace(0, self.H - 1, gr + 2)[1:-1]
        xs = np.linspace(0, self.W - 1, gc + 2)[1:-1]
        hy = self.H / gr
        hx = self.W / gc
        for yi in ys:
            for xi in xs:
                ci = min(Hc - 1, max(0, int(round((yi + 0.5) * Hc / self.H - 0.5))))
                cj = min(Wc - 1, max(0, int(round((xi + 0.5) * Wc / self.W - 0.5))))
                psf = self.blurfield.psf_at(float(X[ci, cj]), float(Y[ci, cj]))
                f = cv2.filter2D(img, cv2.CV_32F, psf, borderType=cv2.BORDER_REFLECT_101)
                w = (np.clip(1 - np.abs(self.U - xi) / hx, 0, None)
                     * np.clip(1 - np.abs(self.V - yi) / hy, 0, None)).astype(np.float32)
                acc += f * w
                wsum += w
        wsum = np.maximum(wsum, 1e-6)
        return np.clip(acc / wsum, 0, 255).astype(np.uint8)


# =============================================================================
def _mix_seed(seed: int, salt: int) -> int:
    return int((np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) ^ np.uint64(salt))
               % np.uint64(2 ** 63 - 1))


def _odd(n: int) -> int:
    n = max(int(n), 1)
    return n if n % 2 == 1 else n + 1
