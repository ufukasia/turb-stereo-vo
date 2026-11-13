"""
Mevcut ön-uca minimal türbülans enjeksiyonu.

TASARIM KISITI (rapor E0 kabul ölçütü): türbülans KAPALIYKEN hattın davranışı
bit düzeyinde değişmemelidir. Bu yüzden enjeksiyon tek bir `if self.turb is None`
kapısının arkasındadır ve kapalıyken hiçbir dizi kopyalanmaz, hiçbir hesap yapılmaz.

Enjeksiyon noktası: REKTİFİKASYONDAN SONRA, ön-işlemeden ÖNCE.
Gerekçesi (makalede beyan edilir): bozulma, rektifiye edilmiş görüntü düzleminde
tanımlı bir açısal tilt alanıdır; rektifikasyon sabit bir yeniden örneklemedir ve
türbülansı ondan önce uygulamak, ölçülen yer değiştirmeyi rektifikasyon jakobiyeni
ile karıştırıp E2'nin ground-truth eğim alanını bulanıklaştırırdı. Rektifiye
düzlemde uygulamak, GT yer değiştirmesini ANALİTİK olarak bilinir kılar.
"""
from __future__ import annotations

import numpy as np

from .simulator import TurbParams, TurbulenceSimulator


def sensor_noise(img: np.ndarray, sigma: float, seed: int, frame_idx: int, cam: int) -> np.ndarray:
    """Sensör gürültüsü (Gauss, sigma DN), türbülanstan SONRA (Faz 2, 2026-11-08).

    `turb/uret_derinlik.py::gurultu_ekle` ile BİREBİR aynı matematik (AirSim yolu): üreteç
    default_rng([gürültü tohumu, kare, kamera]) -> aynı tohum aynı baytlar; kamera ve kare
    bağımsız; float32 + N(0, sigma), rint, [0, 255], uint8. sigma <= 0 -> girdi aynen döner.
    """
    if sigma <= 0:
        return img
    rng = np.random.default_rng([int(seed), int(frame_idx), int(cam)])
    return np.clip(np.rint(img.astype(np.float32) + rng.normal(0.0, sigma, img.shape)),
                   0, 255).astype(np.uint8)


class TurbulenceInjector:
    """cam0/cam1 için ortak atmosfer gerçeklemesi; stereo ayak izi ötelemesiyle."""

    def __init__(self, params: TurbParams, restoration: str = "none",
                 restoration_window: int = 7, restoration_sharpen: float = 0.6,
                 score_quality: bool = True, sensor_noise_dn: float = 0.0,
                 noise_seed: int | None = None):
        self.p = params
        # Sensör gürültüsü (Faz 2): tohumlu temiz dağılım + eşli tasarım için. 0 = kapalı (eski yol).
        self.noise_dn = float(sensor_noise_dn or 0.0)
        self.noise_seed = int(params.seed if noise_seed is None else noise_seed)
        self._sim = None
        self._shape = None
        self.frame_idx = 0
        self.n_applied = 0
        # [E3] restorasyon kolu — türbülanstan SONRA, ön-işlemeden ÖNCE uygulanır.
        self.restorer = None
        self.score_quality = bool(score_quality)
        # restoration="none" + score_quality=True -> ozdes donusum ama PSNR/SSIM olculur.
        # E3'un referans kolu ("bozuk -> VIO") icin gereklidir: algisal kalite ekseni
        # ile ATE ekseni AYNI kosumda kaydedilir.
        if restoration and restoration != "none":
            from .restoration import RestorationHook
            self.restorer = RestorationHook(restoration, restoration_window,
                                            restoration_sharpen)
        elif self.score_quality and params.enabled:
            from .restoration import RestorationHook
            self.restorer = RestorationHook("none", restoration_window,
                                            restoration_sharpen)

    # -------------------------------------------------------------- kurulum
    def _ensure(self, img: np.ndarray):
        shp = img.shape[:2]
        if self._sim is None or self._shape != shp:
            self._shape = shp
            self._sim = TurbulenceSimulator(self.p, image_size=(shp[1], shp[0]))
        return self._sim

    @property
    def sim(self) -> TurbulenceSimulator:
        if self._sim is None:
            raise RuntimeError("TurbulenceInjector henüz bir kare görmedi")
        return self._sim

    def info(self) -> dict:
        return self._sim.info() if self._sim is not None else {"active": False}

    # ------------------------------------------------------------ uygulama
    def apply_stereo(self, rect1: np.ndarray, rect2, frame_idx: int):
        """(cam0, cam1) rektifiye görüntülerine bozulmayı uygula."""
        sim = self._ensure(rect1)
        noisy = self.noise_dn > 0
        if not sim.active and self.restorer is None and not noisy:
            return rect1, rect2
        self.frame_idx = int(frame_idx)
        self.n_applied += 1
        out1 = sim.apply(rect1, self.frame_idx, cam=0) if sim.active else rect1
        out2 = (sim.apply(rect2, self.frame_idx, cam=1)
                if (sim.active and rect2 is not None) else rect2)
        if noisy:                       # sensör gürültüsü: optik bozulmadan SONRA, restorasyondan ÖNCE
            out1 = sensor_noise(out1, self.noise_dn, self.noise_seed, self.frame_idx, 0)
            if out2 is not None:
                out2 = sensor_noise(out2, self.noise_dn, self.noise_seed, self.frame_idx, 1)
        if self.restorer is not None:
            res1 = self.restorer.apply(out1, 0)
            if self.score_quality:
                # PSNR/SSIM referansı TEMİZ rektifiye karedir (C2'nin "algısal" ekseni)
                self.restorer.score(rect1, res1)
            out1 = res1
            if out2 is not None:
                out2 = self.restorer.apply(out2, 1)
        return out1, out2

    def quality_summary(self):
        return self.restorer.summary() if self.restorer is not None else None

    def gt_displacement(self, frame_idx: int, cam: int, u, v):
        """Temiz (u,v) konumlarındaki GT yer değiştirme [px] (E2)."""
        return self._ensure_shape_free().displacement_at(frame_idx, cam, u, v)

    def _ensure_shape_free(self):
        if self._sim is None:
            raise RuntimeError("Simülatör henüz kurulmadı (bir kare işleyin)")
        return self._sim


# =============================================================================
def params_from_cfg(cfg: dict) -> TurbParams | None:
    """config.yaml'ın `turbulence:` bölümünden TurbParams üret (yoksa None)."""
    t = (cfg or {}).get("turbulence") or {}
    if (not t.get("enabled", False) and str(t.get("restoration", "none")) == "none"
            and not float(t.get("sensor_noise_dn", 0.0) or 0.0) > 0):
        return None
    known = set(TurbParams.__dataclass_fields__.keys())
    kw = {k: v for k, v in t.items() if k in known}
    if "blur_grid" in kw and kw["blur_grid"] is not None:
        kw["blur_grid"] = tuple(kw["blur_grid"])
    return TurbParams(**kw)


def injector_from_cfg(cfg: dict) -> TurbulenceInjector | None:
    p = params_from_cfg(cfg)
    if p is None:
        return None
    t = (cfg or {}).get("turbulence") or {}
    return TurbulenceInjector(
        p,
        restoration=str(t.get("restoration", "none")),
        restoration_window=int(t.get("restoration_window", 7)),
        restoration_sharpen=float(t.get("restoration_sharpen", 0.6)),
        score_quality=bool(t.get("score_quality", True)),
        sensor_noise_dn=float(t.get("sensor_noise_dn", 0.0) or 0.0),
        noise_seed=t.get("noise_seed"))
