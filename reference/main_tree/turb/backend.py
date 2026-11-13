"""
Arka-uç ölçüm modeli kolları — PyEco-SLAM uyarlaması (pyeco_turb, 2026-11-07).

    A0      değişiklik yok: kol nesnesi None, PyEco kod yolu birebir                 (temel çizgi)
    A0C     kalibre taban: sigma = k * sigma_0  (k = cal_factor)                     (kontrol)
    A1      oracle skaler şişirme: sigma^2 = sigma_0^2 + s_tot^2,
            s_tot = hypot(1, ind_over_tilt) * sigma_tilt(D/r0)                        (naif alternatif)
    A2      öznitelik-başına (iz artığından) skaler şişirme                           (ara model)
    A3      yama-korelasyonlu Sigma_yama + çevrimiçi DIMM sigma (turb/dimm_v2)        (önerilen)
    A4      A3 + izoplanatik-yama farkındalıklı seçim
    A5      yama-korelasyonlu Sigma_yama + oracle sigma/theta_c                       (üst sınır)
    AG      JENERİK artık-tabanlı skaler ölçek: sigma^2 = sigma_0^2 + r_hat^2,
            r_hat = son `ag_window` graf karesindeki yeniden-izdüşüm artıklarının eksen
            başına RMS'i. Türbülans fiziği YOK (yapı fonksiyonu/DIMM uydurması yok).
            Sage-Husa / VB uyarlamalı R kestiriminin en yalın hâli (denetim B3).
    A0C+A3  kalibre taban x çevrimiçi kestirim: sigma^2 = k^2 sigma_0^2 + s_hat^2,
            s_hat = dimm_sigma_scale * sigma_hat_DIMM (A3 ile aynı kestirici), skaler.

sigma_0 = PyEco'nun px_noise (mono) / stereo_sigma (stereo). Stereo modelde aynı ek
varyans kullanılır (oran korunmaz, eklenir) — eski hatla aynı.

ROBUST REJİM (denetim B4'ün düzeltmesi). Eski hatta A0 mono faktörlerde Huber
kullanıyor, kollar `robust: 0` ile hiç kullanmıyordu; yani kollar yalnız kovaryansta
farklı değildi. Burada PyEco'nun varsayılan rejimi BÜTÜN kollara aynen uygulanır:
    mono   : Robust(Huber(huber_mono_k = 2.4477)),  stereo: Robust(Huber(huber_stereo_k = 2.7955))
eşikler beyazlatılmış birimdedir, yani sigma ne olursa olsun aynı kural. Kol başına
robust bayrağı YOKTUR; rejim backend'in kendi argümanlarından okunur. Yama faktöründe
(A3/A4/A5) Huber GÖZLEM BAŞINA, aynı eşiklerle uygulanır (turb/factors.py belgesi).
Korelasyonsuz kalan yamalar (tek gözlem ya da köşegen Sigma, ör. A5'te D/r0=0)
bireysel GTSAM faktörleri olarak eklenir.

KAPSAM: kollar yalnız PENCERE (smoother) yeniden-izdüşüm faktörlerinin kovaryansını
değiştirir: ana gözlem döngüsü, yerel harita eşleşmeleri, search-by-projection ve
çok-görüşlü ilklemenin geçmiş-poz faktörleri. İzleme (PnP-RANSAC, motion-only BA,
chi2 kapıları) bütün kollarda aynıdır. Kestiriciler (DIMM, AG, A2) yalnız graf
karelerinin (anahtar kare kipinde anahtar karelerin) gözlemlerinden beslenir ve kare f'in
gürültüsü, f'ten ÖNCEKİ graf karelerinin artıklarıyla belirlenir (nedensel).

Bireysel-faktör kolları (A0C/A1/A2/AG/A0C+A3) faktörü PyEco'nun kendi ekleme
noktasında, yalnız gürültü modelini değiştirerek ekler; böylece faktör sırası A0 ile
aynı kalır ve A1, D/r0=0'da A0 ile BİT DÜZEYİNDE aynıdır (P13). Yama kolları
gözlemleri kare sonunda gruplayıp `flush` ile ekler.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import gtsam

from . import covariance as _cov
from . import optics as _optics
from .factors import make_patch_factor


ARMS = ("A0", "A0C", "A1", "A2", "A3", "A4", "A5", "AG", "A0C+A3")
INLINE_ARMS = ("A0C", "A1", "A2", "AG", "A0C+A3")
PATCH_ARMS = ("A3", "A4", "A5")
DIMM_ARMS = ("A3", "A4", "A0C+A3")


@dataclass
class ArmConfig:
    arm: str = "A0"
    corr_model: str = "empirical"
    D_over_r0: float = 0.0
    theta_c_over_fov: float = 0.15
    theta_c_nominal_over_fov: float = 0.15
    aperture_D: float = 0.10
    wavelength: float = 550e-9
    px_per_rad: float = 170648.0
    image_width: int = 752
    max_per_patch: int = 5
    a4_keep_per_patch: int = 5
    cal_factor: float = 16.0
    shrinkage: float = 0.05
    ind_over_tilt: float = 1.75
    dimm_version: int = 2
    dimm_window: int = 60
    dimm_update_every: int = 10
    dimm_sigma_scale: float = 1.0
    ag_window: int = 60
    sigma_floor_px: float = 0.05
    huber_mono_k: float = 2.4477
    huber_stereo_k: float = 2.7955

    @property
    def theta_c(self) -> float:
        return self.theta_c_over_fov * (self.image_width / self.px_per_rad)

    @property
    def theta_c_px(self) -> float:
        return self.theta_c * self.px_per_rad

    @property
    def theta_c_nominal_px(self) -> float:
        return self.theta_c_nominal_over_fov * self.image_width

    @property
    def z_layer(self) -> float:
        return _optics.layer_distance_for_theta_c(self.aperture_D, self.theta_c)

    @property
    def sigma_px_oracle(self) -> float:
        if self.D_over_r0 <= 0:
            return 0.0
        return _optics.sigma_tilt_px(self.D_over_r0, self.wavelength,
                                     self.aperture_D, self.px_per_rad)


def _L(i):
    return gtsam.symbol("l", int(i))


class TurbMeasurementBuilder:
    """Kare bazında kolun gürültü modellerini üretir; yama kollarında faktörleri kurar."""

    def __init__(self, cfg: ArmConfig, fx, fy, cx, cy, baseline, px_noise, stereo_sigma,
                 K_gtsam, K_stereo):
        self.c = cfg
        self.fx, self.fy, self.cx, self.cy = fx, fy, cx, cy
        self.baseline = float(baseline)
        self.px_noise = float(px_noise)
        self.stereo_sigma = float(stereo_sigma)
        self.K_gtsam, self.K_stereo = K_gtsam, K_stereo
        self.inline = cfg.arm in INLINE_ARMS
        self.dimm = None
        if cfg.arm in DIMM_ARMS:
            if int(cfg.dimm_version) == 2:
                from .dimm_v2 import DifferentialMotionMonitorV2 as _DM
            else:
                from .dimm import DifferentialMotionMonitor as _DM
            self.dimm = _DM(cfg.px_per_rad, cfg.aperture_D, cfg.wavelength,
                            window_frames=cfg.dimm_window, base_sigma_px=px_noise)
        self._needs_residuals = cfg.arm in DIMM_ARMS or cfg.arm in ("AG", "A2")
        self._track_res = {}          # A2: lid -> (n, toplam yarı-kare artık)
        self._ag_hist = []            # AG: [(toplam kare artık, bileşen sayısı)] graf karesi başına
        self._frame_uv, self._frame_res = [], []
        self._obs = []                # yama kolları: (lid, u, v, Z|None)
        self._must_keep = set()
        self._frame_key = None
        self._n_frames = 0
        self._sigma_hat = np.nan      # DIMM kestirimi (eksen başına, px)
        self._theta_c_hat_px = np.nan
        self._ag_rms = np.nan
        self._models = {}             # (mono_sigma, stereo_sigma) -> (mono_model, stereo_model)
        self.cur_sigma_turb = 0.0
        self.stats = {"arm": cfg.arm, "n_patch_factors": 0, "n_patch_single": 0,
                      "n_patch_uncorrelated": 0, "n_obs_patch": 0, "n_fac_inline": 0,
                      "sigma_hat": [], "theta_c_hat": [], "r0_hat": [], "dimm_pinned": 0,
                      "ag_rms": [], "sigma_turb_used": [],
                      "robust": {"huber_mono_k": cfg.huber_mono_k,
                                 "huber_stereo_k": cfg.huber_stereo_k,
                                 "granularity": "per-observation (all arms)"}}

    # ----------------------------------------------------------------- modeller
    def _model_pair(self, s_mono: float, s_stereo: float):
        """PyEco'nun robust rejimiyle sarılmış (mono, stereo) gürültü modeli."""
        key = (float(s_mono), float(s_stereo))
        if key not in self._models:
            mono = gtsam.noiseModel.Robust.Create(
                gtsam.noiseModel.mEstimator.Huber.Create(float(self.c.huber_mono_k)),
                gtsam.noiseModel.Isotropic.Sigma(2, key[0]))
            st = gtsam.noiseModel.Isotropic.Sigma(3, key[1])
            if float(self.c.huber_stereo_k or 0.0) > 0:
                st = gtsam.noiseModel.Robust.Create(
                    gtsam.noiseModel.mEstimator.Huber.Create(float(self.c.huber_stereo_k)), st)
            self._models[key] = (mono, st)
        return self._models[key]

    def _dimm_sigma(self):
        """A3/A4/A0C+A3: DIMM sigma kestirimi (theta_c sabit nominal; bkz. eski belge)."""
        c = self.c
        if (self._n_frames % max(1, c.dimm_update_every)) != 0 and np.isfinite(self._sigma_hat):
            return float(self._sigma_hat), float(self._theta_c_hat_px)
        tc_nom = c.theta_c_nominal_px
        s, tc, r0 = self.dimm.estimate(theta_c_px=tc_nom)
        if not np.isfinite(s):
            s = self._sigma_hat if np.isfinite(self._sigma_hat) else c.sigma_floor_px
            tc = tc_nom
        else:
            self._sigma_hat, self._theta_c_hat_px = s, tc
            self.stats["sigma_hat"].append(float(s))
            self.stats["theta_c_hat"].append(float(tc))
            self.stats["r0_hat"].append(float(r0))
            if bool(getattr(self.dimm, "last_pinned", False)):
                self.stats["dimm_pinned"] += 1
        return float(s), float(tc)

    def _ag_sigma(self):
        if (self._n_frames % max(1, self.c.dimm_update_every)) != 0 and np.isfinite(self._ag_rms):
            return float(self._ag_rms)
        S = sum(h[0] for h in self._ag_hist)
        n = sum(h[1] for h in self._ag_hist)
        if n >= 50:
            self._ag_rms = float(np.sqrt(S / n))
            self.stats["ag_rms"].append(self._ag_rms)
        return float(self._ag_rms) if np.isfinite(self._ag_rms) else 0.0

    def current_sigma(self):
        """(sigma_turb [px] eklenen, theta_c [px]) — kolun kullandığı ek ölçek."""
        c = self.c
        if c.arm in ("A1", "A5"):
            s = c.sigma_px_oracle
            return float(np.hypot(s, c.ind_over_tilt * s)), c.theta_c_px
        if c.arm in DIMM_ARMS:
            s, tc = self._dimm_sigma()
            return float(c.dimm_sigma_scale) * s, tc
        if c.arm == "AG":
            return self._ag_sigma(), np.nan
        return 0.0, np.nan

    # ------------------------------------------------------------------- kare
    def begin_frame(self, frame_key):
        """Graf karesinin başı. Dönüş: (mono_model, stereo_model) — bireysel faktör modeli."""
        self._frame_key = frame_key
        self._obs = []
        self._must_keep = set()
        self._frame_uv, self._frame_res = [], []
        c = self.c
        if c.arm == "A0C":
            s_t, self._theta_cur = 0.0, np.nan
            sm, ss = c.cal_factor * self.px_noise, c.cal_factor * self.stereo_sigma
        else:
            s_t, self._theta_cur = self.current_sigma()
            k = c.cal_factor if c.arm == "A0C+A3" else 1.0
            sm = float(np.sqrt((k * self.px_noise) ** 2 + s_t ** 2))
            ss = float(np.sqrt((k * self.stereo_sigma) ** 2 + s_t ** 2))
        self.cur_sigma_turb = float(s_t)
        self.stats["sigma_turb_used"].append(float(s_t))
        self._cur = self._model_pair(sm, ss)
        return self._cur

    def model_for(self, lid, stereo: bool):
        """A2: iz artığından öznitelik başına model; diğer kollarda karenin modeli."""
        if self.c.arm != "A2":
            return self._cur[1] if stereo else self._cur[0]
        n, sq = self._track_res.get(int(lid), (0, 0.0))
        extra = np.sqrt(max(sq / n - self.px_noise ** 2, 0.0)) if n >= 3 else 0.0
        m = self._model_pair(_cov.scalar_inflation(self.px_noise, extra),
                             _cov.scalar_inflation(self.stereo_sigma, extra))
        return m[1] if stereo else m[0]

    def observe(self, lid, u, v, pred):
        """Grafa giren bir gözlemin artığını kestiricilere kaydet (pred: (u,v) ya da None)."""
        if not self._needs_residuals or pred is None:
            return
        r = (float(u) - float(pred[0]), float(v) - float(pred[1]))
        self._frame_uv.append((float(u), float(v)))
        self._frame_res.append(r)
        if self.c.arm == "A2":
            n, s = self._track_res.get(int(lid), (0, 0.0))
            self._track_res[int(lid)] = (n + 1, s + 0.5 * (r[0] ** 2 + r[1] ** 2))

    def add(self, lid, u, v, Z=None, pred=None, must_keep=False):
        """Yama kolları: gözlemi kare sonuna ertele. must_keep: bu karede grafa giren yeni nokta
        (A4 seçimi atamaz; atılırsa değişken faktörsüz kalır ve sistem tekilleşir)."""
        self._obs.append((int(lid), float(u), float(v), None if Z is None else float(Z)))
        if must_keep:
            self._must_keep.add(int(lid))
        self.observe(lid, u, v, pred)

    def flush(self, graph):
        """Graf karesinin sonu: yama faktörlerini ekle, artıkları kestiricilere it."""
        if self.c.arm in PATCH_ARMS and self._obs:
            self._flush_patch(graph)
        self._obs = []
        if self._frame_res:
            res = np.asarray(self._frame_res, dtype=np.float64)
            if self.dimm is not None and len(res) >= 6:
                self.dimm.add_frame(np.asarray(self._frame_uv, dtype=np.float64), res)
            if self.c.arm == "AG":
                self._ag_hist.append((float(np.sum(res ** 2)), int(res.size)))
                if len(self._ag_hist) > int(self.c.ag_window):
                    self._ag_hist.pop(0)
        self._frame_uv, self._frame_res = [], []
        self._n_frames += 1

    # ------------------------------------------------------------ A3/A4/A5
    def _add_single(self, graph, lid, u, v, Z, s_tot):
        m, st = self._model_pair(_cov.scalar_inflation(self.px_noise, s_tot),
                                 _cov.scalar_inflation(self.stereo_sigma, s_tot))
        if Z is None:
            graph.add(gtsam.GenericProjectionFactorCal3_S2(
                np.array([u, v]), m, self._frame_key, _L(lid), self.K_gtsam, False, False))
        else:
            uR = u - (self.fx * self.baseline / Z)
            graph.add(gtsam.GenericStereoFactor3D(
                gtsam.StereoPoint2(u, uR, v), st, self._frame_key, _L(lid), self.K_stereo))

    def _flush_patch(self, graph):
        c = self.c
        s_tot, theta_c_px = self.cur_sigma_turb, self._theta_cur
        obs = self._obs
        uv = np.array([[o[1], o[2]] for o in obs])
        if c.arm == "A4":
            pri = np.array([1.0 if o[3] is not None else 0.0 for o in obs])
            keep = _cov.select_patch_aware(uv, pri, theta_c_px,
                                           max_per_patch=c.a4_keep_per_patch, sigma_px=s_tot)
            keep |= np.array([o[0] in self._must_keep for o in obs], dtype=bool)
            obs = [o for o, k in zip(obs, keep) if k]
            uv = uv[keep]
            if not obs:
                return
        if c.arm == "A5" or not np.isfinite(theta_c_px):
            theta_c_use, z_layer = c.theta_c, c.z_layer
            theta_c_px = c.theta_c_px
        else:
            theta_c_use = float(theta_c_px) / c.px_per_rad
            z_layer = _optics.layer_distance_for_theta_c(c.aperture_D, theta_c_use)
        patches = _cov.group_isoplanatic(uv, theta_c_px, max_per_patch=c.max_per_patch)
        # s_tot: korelasyonlu + bağımsız TOPLAM ölçek (A5: hypot(1, ind_over_tilt)*sigma_tilt;
        # A3: dimm_sigma_scale*sigma_hat). empirical çekirdekte nugget genliğin içindedir.
        sigma_tilt = (s_tot / float(np.hypot(1.0, c.ind_over_tilt))) if c.arm == "A5" else s_tot
        for idx in patches:
            grp = [obs[i] for i in idx]
            if len(grp) == 1:
                lid, u, v, Z = grp[0]
                self._add_single(graph, lid, u, v, Z, s_tot)
                self.stats["n_patch_single"] += 1
                continue
            is_st = np.array([g[3] is not None for g in grp])
            theta = np.stack([[(g[1] - self.cx) / c.px_per_rad,
                               (g[2] - self.cy) / c.px_per_rad] for g in grp])
            meas = []
            for lid, u, v, Z in grp:
                meas += [u, v] if Z is None else [u, u - (self.fx * self.baseline / Z), v]
            if c.corr_model == "empirical":
                sig_kernel, sig_ind = s_tot, 0.0
            else:
                sig_kernel, sig_ind = sigma_tilt, float(np.sqrt(max(s_tot ** 2 - sigma_tilt ** 2, 0.0)))
            Sigma = _cov.patch_covariance(
                theta, is_st, sig_kernel, self.px_noise, self.stereo_sigma,
                z_layer, c.aperture_D, self.baseline,
                model=c.corr_model, theta_c=theta_c_use,
                shrinkage=c.shrinkage, sigma_ind_px=sig_ind)
            off = Sigma - np.diag(np.diag(Sigma))
            if not np.any(off):
                # korelasyonsuz yama (ör. sigma=0): bireysel faktörler, aynı robust rejim
                for lid, u, v, Z in grp:
                    self._add_single(graph, lid, u, v, Z, s_tot)
                self.stats["n_patch_uncorrelated"] += 1
                continue
            f, _ = make_patch_factor(self._frame_key, [_L(g[0]) for g in grp],
                                     np.array(meas), is_st, Sigma,
                                     self.fx, self.fy, self.cx, self.cy, self.baseline,
                                     robust=False, huber_k_mono=c.huber_mono_k,
                                     huber_k_stereo=(c.huber_stereo_k or None))
            graph.add(f)
            self.stats["n_patch_factors"] += 1
        self.stats["n_obs_patch"] += len(obs)

    # ------------------------------------------------------------- raporlama
    def summary(self) -> dict:
        s = dict(self.stats)
        for k in ("sigma_hat", "theta_c_hat", "r0_hat", "ag_rms", "sigma_turb_used"):
            v = np.array(s.pop(k), dtype=float)
            v = v[np.isfinite(v)]
            s[k + "_mean"] = float(v.mean()) if v.size else float("nan")
            s[k + "_median"] = float(np.median(v)) if v.size else float("nan")
            s[k + "_n"] = int(v.size)
        s["sigma_px_oracle"] = self.c.sigma_px_oracle
        s["theta_c_px_oracle"] = self.c.theta_c_px
        s["cal_factor"] = self.c.cal_factor
        s["dimm_version"] = int(self.c.dimm_version)
        s["dimm_sigma_scale"] = float(self.c.dimm_sigma_scale)
        return s


def build_from_args(args, fx, fy, cx, cy, baseline, K_gtsam, K_stereo):
    """gtsam_backend argümanlarından kol nesnesi (A0 -> None: PyEco kod yolu birebir)."""
    arm = str(getattr(args, "turb_arm", "A0") or "A0").upper()
    if arm == "A0":
        return None
    if arm not in ARMS:
        raise ValueError(f"bilinmeyen kol: {arm} (geçerli: {', '.join(ARMS)})")
    if int(getattr(args, "karanlik_sigma", 0) or 0):
        raise ValueError("karanlik_sigma ile türbülans kolları birlikte kullanılamaz")
    cfg = ArmConfig(
        arm=arm,
        corr_model=str(getattr(args, "turb_corr_model", "empirical")),
        D_over_r0=float(getattr(args, "turb_D_over_r0", 0.0)),
        theta_c_over_fov=float(getattr(args, "turb_theta_c_over_fov", 0.15)),
        theta_c_nominal_over_fov=float(getattr(args, "turb_theta_c_nominal_over_fov", 0.15)),
        aperture_D=float(getattr(args, "turb_aperture_D", 0.10)),
        wavelength=float(getattr(args, "turb_wavelength", 550e-9)),
        px_per_rad=float(getattr(args, "turb_px_per_rad", 170648.0)),
        image_width=int(getattr(args, "turb_image_width", 752)),
        max_per_patch=int(getattr(args, "turb_max_per_patch", 5)),
        a4_keep_per_patch=int(getattr(args, "turb_a4_keep", 5)),
        shrinkage=float(getattr(args, "turb_shrinkage", 0.05)),
        cal_factor=float(getattr(args, "turb_cal_factor", 16.0)),
        ind_over_tilt=float(getattr(args, "turb_ind_over_tilt", 1.75)),
        dimm_version=int(getattr(args, "turb_dimm_version", 2)),
        dimm_window=int(getattr(args, "turb_dimm_window", 60)),
        dimm_update_every=int(getattr(args, "turb_dimm_every", 10)),
        dimm_sigma_scale=float(getattr(args, "turb_dimm_sigma_scale", 1.0)),
        ag_window=int(getattr(args, "turb_ag_window", 60)),
        huber_mono_k=float(getattr(args, "huber_mono_k", 1.345)),
        huber_stereo_k=float(getattr(args, "huber_stereo_k", 0.0) or 0.0),
    )
    return TurbMeasurementBuilder(cfg, fx, fy, cx, cy, baseline,
                                  float(args.px_noise), float(args.stereo_sigma),
                                  K_gtsam, K_stereo)
