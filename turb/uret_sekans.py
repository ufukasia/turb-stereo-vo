"""
Önceden bozulmuş (türbülanslı) EuRoC-ASL sekans üreticisi.

Türbülansı görüntülere DİSKTE uygular; çıktı aynı ASL düzeninde yeni bir sekanstır
(`mav0/cam0`, `mav0/cam1`, sensor.yaml'lar, GT/IMU kopyası). Böylece aynı bozulmuş
veri PyEco dışında DEĞİŞTİRİLMEMİŞ sistemlere (ORB-SLAM3, AirSLAM, DROID-SLAM ...)
verilebilir.

REKTİFİKASYON KARARI
--------------------
Varsayılan (a) `rektifiye`: ham görüntü -> stereo rektifikasyon (eski/PyEco kancasının
BİREBİR aynısı: cv2.stereoRectify CALIB_ZERO_DISPARITY, alpha=0, aynı taban matrisi,
INTER_LINEAR remap) -> türbülans -> rektifiye pinhole sekans. sensor.yaml'lar rektifiye
K, sıfır bozulma, döndürülmüş T_BS ve saf x-ötelemeli `T_cn_cnm1` ile yazılır.
Gerekçe: (1) kancanın ön-uca verdiği görüntüyle bayt bayt aynı görüntü üretir;
(2) simülatörün yer değiştirme alanı rektifiye düzlemde tanımlıdır, oracle/GT alanı
analitik kalır (rektifikasyon jakobiyeni karışmaz); (3) tüketici sistem kendi
rektifikasyonunu yaptığında kimlik dönüşüm olur, bozulma YENİDEN örneklenmez.
(b) `ham`: türbülans ham (bozulmalı) görüntüye uygulanır, kalibrasyon kopyalanır.
Fiziksel sıra (atmosfer optikten önce) budur ama ham görüntü düzlemi lens
bozulmasıyla çarpık olduğundan "atmosfer -> lens" sırasını da TAM temsil etmez
(alan bozulmuş piksel ızgarasında üretilir). Yalnız ablasyon için.
`kimlik`: girdi zaten rektifiye pinhole (AirSim); görüntü dokunulmadan geçer.
`otomatik`: bozulma sıfır + K'lar eşit + T_cn_cnm1 saf x-ötelemesiyse kimlik,
aksi halde rektifiye.

KARE İNDEKSİ: türbülans zaman ekseni, cam0 zaman damgalarının SAYISAL sıralı
listesindeki MUTLAK indekstir (kanca: `_seq_frame_idx = i`). `--kareler 400:900`
yalnız o kareleri yazar ama türbülans o karelerin mutlak indeksleriyle üretilir;
yani pencere çıktısı tam sekans koşumundaki kancayla aynıdır.

KULLANIM
    python -m turb.uret_sekans --girdi "/data/sdata/MH01" \
        --cikti_kok "/data/sdata_turb" --D_over_r0 5 --tohum 0 \
        --kareler 400:900
    # çok pencere + ısınma:  --kareler 400:900,800:1300,...
    # AirSim:  --girdi ~/sim_envs/sdata_airsim/NH_egik_h5_v1 --kareler 0:300
    # Faz 2 sensör gürültüsü (ON_KAYIT v0.3: sigma = 1 DN, gürültü tohumu = hücre tohumu):
    #          --gurultu_dn 1        (varsayılan 0 = kapalı: çıktı baytları eskisiyle aynı)

STEREO DIŞ PARAMETRESİ (2026-11-08): `--stereo_dis sensor_yaml` rektifikasyonda cam0->cam1 dönüşümünü
iki kameranın T_BS'sinden alır (inv(T_B_c1) T_B_c0) — pyeco_turb/config.yaml `frontend.stereo_extrinsics:
sensor_yaml` ile PyEco'nun Faz 2'de gördüğü rektifiye görüntünün AYNISI. Varsayılan `gomulu` eski kanca sabit
matrisidir (kod/ ön-ucu); EuRoC'ta ~0.8 derece farklı ve rektifiye görüntüde ~6.7 px dikey kayma bırakır
(ON_TEST_RAKIPLER.md §4). Faz 2 rakip sekansları `sensor_yaml` ile üretilir (ada `_sy` eklenir).

SENSÖR GÜRÜLTÜSÜ (2026-11-08): `--gurultu_dn > 0` ise türbülanstan SONRA her kareye
`turb/hooks.py::sensor_noise` (PyEco kancasının kullandığı fonksiyonun KENDİSİ) uygulanır:
Gauss sigma DN, rint + [0, 255] kırpma + uint8, üreteç default_rng([gürültü tohumu, kare, kamera]);
kare = yukarıdaki MUTLAK indeks. Gürültü tohumu varsayılan olarak türbülans tohumudur (kanca:
noise_seed None -> params.seed). D/r0 = 0'da da uygulanır. Kapalıyken kod yolu ve çıktı aynıdır.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from dataclasses import asdict
from pathlib import Path

import cv2
import numpy as np
import yaml

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    __package__ = "turb"

from turb import optics  # noqa: E402
from turb.simulator import TurbParams, TurbulenceSimulator  # noqa: E402
from turb.hooks import sensor_noise  # noqa: E402  (PyEco kancasıyla aynı fonksiyon)

MODUL_DIZINI = Path(__file__).resolve().parent
META_SURUM = "turb_meta/1"
#: Kancanın (kod/stereo_front_end.py::load_calibration) varsayılan cam0->cam1 matrisi.
EUROC_T_C1_C0 = np.array([
    [0.99999609, 0.0023145, -0.00136141, -0.110073],
    [-0.00231508, 0.99999693, -0.00042232, -0.000399],
    [0.00136043, 0.00042547, 0.99999868, -0.000284],
    [0., 0., 0., 1.]])
#: A1 kolunun bağımsız/tilt oranı (backend.ArmConfig.ind_over_tilt).
A1_IND_OVER_TILT = 1.75
VARSAYILAN_YAPILANDIRMA = (MODUL_DIZINI.parent.parent / "kod" / "config_euroc.yaml")
PARMAK_IZI_DOSYALARI = ("simulator.py", "fields.py", "optics.py", "screen.py",
                        "uret_sekans.py")


# =============================================================================
#  Yardımcılar
# =============================================================================
def sha256_dosya(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def kareleri_ayristir(spec: str | None, n: int) -> tuple[list[int], list[tuple[int, int]]]:
    """'400:900,1200:1500' -> (sıralı tekil indeksler, aralıklar). Yarı açık [a, b).

    Boş/None/'tum' -> tüm kareler. 'a:' sona kadar, ':b' baştan. Tek sayı = tek kare.
    """
    if spec is None or str(spec).strip() in ("", "tum", "all", ":"):
        return list(range(n)), [(0, n)]
    araliklar = []
    for parca in str(spec).split(","):
        parca = parca.strip()
        if not parca:
            continue
        if ":" in parca:
            a, b = parca.split(":", 1)
            a = int(a) if a.strip() else 0
            b = int(b) if b.strip() else n
        else:
            a = int(parca)
            b = a + 1
        a, b = max(0, a), min(n, b)
        if b <= a:
            raise ValueError(f"boş/sekans dışı kare aralığı: {parca!r} (n={n})")
        araliklar.append((a, b))
    idx = sorted({i for a, b in araliklar for i in range(a, b)})
    return idx, araliklar


def _csv_oku(p: Path) -> dict[int, str]:
    out = {}
    with open(p) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            ts, fn = [s.strip() for s in ln.split(",")[:2]]
            out[int(ts)] = fn
    return out


def asl_oku(seq: Path) -> dict:
    """ASL sekansının kamera kalibrasyonu ve kare listesi."""
    m = seq / "mav0"
    with open(m / "cam0/sensor.yaml") as f:
        c0 = yaml.safe_load(f)
    with open(m / "cam1/sensor.yaml") as f:
        c1 = yaml.safe_load(f)
    f0 = _csv_oku(m / "cam0/data.csv")
    f1 = _csv_oku(m / "cam1/data.csv")
    # Kanca ile aynı: cam0 damgalarının SAYISAL sıralı tekil listesi.
    ts = sorted(f0.keys())
    return {"cam0": c0, "cam1": c1, "f0": f0, "f1": f1, "ts": ts}


def _K(c) -> np.ndarray:
    i = c["intrinsics"]
    return np.array([[i[0], 0, i[2]], [0, i[1], i[3]], [0, 0, 1]], dtype=np.float64)


def _T_BS(c) -> np.ndarray:
    return np.asarray(c["T_BS"]["data"], dtype=np.float64).reshape(4, 4)


# =============================================================================
#  Rektifikasyon (kancanın birebir kopyası)
# =============================================================================
class Rektifikator:
    def __init__(self, cal: dict, kip: str = "otomatik", stereo_dis: str = "gomulu"):
        c0, c1 = cal["cam0"], cal["cam1"]
        self.K1, self.K2 = _K(c0), _K(c1)
        self.D1 = np.array(c0["distortion_coefficients"], dtype=np.float64)
        self.D2 = np.array(c1["distortion_coefficients"], dtype=np.float64)
        res = c0.get("resolution")
        self.size = (int(res[0]), int(res[1])) if res else (752, 480)
        if "T_cn_cnm1" in c1:
            self.T_c1_c0 = np.array(c1["T_cn_cnm1"], dtype=np.float64).reshape(4, 4)
            self.taban_kaynagi = "cam1/sensor.yaml:T_cn_cnm1"
        elif stereo_dis == "sensor_yaml":
            # PyEco frontend.stereo_extrinsics: sensor_yaml (pyeco_turb/config.yaml; Faz 2 ızgarası):
            # T_c1_c0 = inv(T_B_c1) @ T_B_c0  (stereo_front_end.py::load_calibration ile aynı ifade)
            self.T_c1_c0 = np.linalg.inv(_T_BS(c1)) @ _T_BS(c0)
            self.taban_kaynagi = "sensor.yaml T_BS: inv(T_B_c1) T_B_c0 (PyEco stereo_extrinsics: sensor_yaml)"
        elif stereo_dis == "gomulu":
            self.T_c1_c0 = EUROC_T_C1_C0.copy()
            self.taban_kaynagi = "kanca varsayilani (EuRoC sabit matrisi)"
        else:
            raise ValueError(f"bilinmeyen stereo_dis: {stereo_dis}")
        self.stereo_dis = stereo_dis
        self.c0, self.c1 = c0, c1
        if kip == "otomatik":
            kip = "kimlik" if self._zaten_rektifiye() else "rektifiye"
        if kip not in ("rektifiye", "kimlik", "ham"):
            raise ValueError(f"bilinmeyen rektifikasyon kipi: {kip}")
        self.kip = kip
        self.maps = None
        R, T = self.T_c1_c0[:3, :3], self.T_c1_c0[:3, 3]
        # Rektifiye geometri her kipte hesaplanır (taban çizgisi bilgisi için).
        (self.R1, self.R2, self.P1, self.P2, self.Q, _, _) = cv2.stereoRectify(
            self.K1, self.D1, self.K2, self.D2, self.size, R, T,
            flags=cv2.CALIB_ZERO_DISPARITY, alpha=0)
        self.baseline = float(-self.P2[0, 3] / self.P2[0, 0])
        if kip == "rektifiye":
            m1 = cv2.initUndistortRectifyMap(self.K1, self.D1, self.R1, self.P1,
                                             self.size, cv2.CV_32FC1)
            m2 = cv2.initUndistortRectifyMap(self.K2, self.D2, self.R2, self.P2,
                                             self.size, cv2.CV_32FC1)
            self.maps = (m1, m2)
        else:
            self.baseline = float(abs(self.T_c1_c0[0, 3])) if kip == "kimlik" \
                else float(np.linalg.norm(T))

    def _zaten_rektifiye(self) -> bool:
        T = self.T_c1_c0
        return (np.allclose(self.D1, 0) and np.allclose(self.D2, 0)
                and np.allclose(self.K1, self.K2)
                and np.allclose(T[:3, :3], np.eye(3), atol=1e-9)
                and np.allclose(T[1:3, 3], 0, atol=1e-9))

    def uygula(self, i0, i1):
        if self.kip != "rektifiye":
            return i0, i1
        (m1x, m1y), (m2x, m2y) = self.maps
        return (cv2.remap(i0, m1x, m1y, cv2.INTER_LINEAR),
                cv2.remap(i1, m2x, m2y, cv2.INTER_LINEAR))

    def bilgi(self) -> dict:
        d = {"kip": self.kip, "taban_matrisi_kaynagi": self.taban_kaynagi,
             "T_c1_c0_kullanilan": self.T_c1_c0.tolist(), "baseline_m": self.baseline,
             "goruntu_boyutu": list(self.size)}
        if self.kip == "rektifiye":
            d.update(R1=self.R1.tolist(), R2=self.R2.tolist(), P1=self.P1.tolist(),
                     P2=self.P2.tolist(),
                     K_rekt=[float(self.P1[0, 0]), float(self.P1[1, 1]),
                             float(self.P1[0, 2]), float(self.P1[1, 2])])
            # Bilgi: kanca sabit matrisi ile T_BS'den türetilen cam0->cam1 farkı
            try:
                T_bs = np.linalg.inv(_T_BS(self.c1)) @ _T_BS(self.c0)
                d["T_BS_turetilen_fark_maks"] = float(np.abs(T_bs - self.T_c1_c0).max())
            except Exception:  # noqa: BLE001
                pass
        return d

    def sensor_yaml_yaz(self, hedef_mav0: Path, ad: str):
        """Rektifiye kipte yeni pinhole sensor.yaml'lar; diğer kiplerde kaynak kopya."""
        if self.kip != "rektifiye":
            return False
        fx, fy, cx, cy = self.P1[0, 0], self.P1[1, 1], self.P1[0, 2], self.P1[1, 2]
        B = self.baseline
        R1 = np.eye(4); R1[:3, :3] = self.R1
        T_B_r0 = _T_BS(self.c0) @ R1.T                 # T_B_cam0 * T_cam0_rect0
        T_r1_r0 = np.eye(4); T_r1_r0[0, 3] = -B        # rect1 = rect0 + B x
        T_B_r1 = T_B_r0 @ np.linalg.inv(T_r1_r0)
        rate = self.c0.get("rate_hz", 20)
        w, h = self.size
        for cam, T_BS, T_cn in ((0, T_B_r0, None), (1, T_B_r1, T_r1_r0)):
            src = self.c0 if cam == 0 else self.c1
            satir = [
                "# Turb-VIO uret_sekans: REKTIFIYE pinhole (bozulma 0).",
                "# Kaynak: " + str(src.get("comment", "")).replace("\n", " "),
                "sensor_type: camera",
                f"comment: rektifiye cam{cam} ({ad})",
                "T_BS:", "  cols: 4", "  rows: 4",
                "  data: [" + ", ".join(f"{v:.12f}" for v in T_BS.ravel()) + "]",
                f"rate_hz: {rate}",
                f"resolution: [{w}, {h}]",
                "camera_model: pinhole",
                f"intrinsics: [{fx:.9f}, {fy:.9f}, {cx:.9f}, {cy:.9f}]",
                "distortion_model: radial-tangential",
                "distortion_coefficients: [0.0, 0.0, 0.0, 0.0]"]
            if T_cn is not None:
                satir.append("T_cn_cnm1: [" + ", ".join(f"{v:.12f}" for v in T_cn.ravel())
                             + "]")
            (hedef_mav0 / f"cam{cam}/sensor.yaml").write_text("\n".join(satir) + "\n")
        return True


# =============================================================================
#  Derinliğe bağlı kip — ARAYÜZ (uygulama bir sonraki adım)
# =============================================================================
class DerinlikOkuyucu:
    """AirSim `mav0/depth0/data/<ts>.npy` okuyucu (float16 düzlemsel derinlik [m]).

    65500 (float16 tavanı) gibi değerler "boş/sonsuz" kabul edilir ve +inf döner.
    """
    GECERSIZ_ESIK = 6.0e4

    def __init__(self, seq: Path):
        self.dizin = Path(seq) / "mav0/depth0/data"

    def var_mi(self) -> bool:
        return self.dizin.is_dir()

    def yol(self, ts: int) -> Path:
        return self.dizin / f"{int(ts)}.npy"

    def oku(self, ts: int) -> np.ndarray:
        d = np.load(self.yol(ts)).astype(np.float32)
        d[~np.isfinite(d) | (d >= self.GECERSIZ_ESIK) | (d <= 0)] = np.inf
        return d


def derinlik_buyut(d: np.ndarray, hedef_wh=(1024, 768), kenar_orani: float = 1.15
                   ) -> np.ndarray:
    """Düzlemsel derinliği (ör. 256x192) görüntü çözünürlüğüne (1024x768) büyüt.

    Yöntem: TERS derinlikte bilinear (yarım-piksel merkez kuralı, cv2.resize ile
    aynı). Düzlemsel derinlikte bir düzlemin 1/Z'si piksel koordinatlarında afindir;
    yani düzlemler (zemin, duvar) TAM korunur ve +inf (gökyüzü) 1/Z = 0 olarak
    doğal biçimde interpole edilir. Derinlik süreksizliklerinde (2x2 destek içinde
    maks/min derinlik > kenar_orani ya da sonlu/sonsuz karışımı) en-yakın komşuya
    geçilir: ön/arka plan arasında "uçan" ara derinlik üretilmez.
    Dönüş: float32, (H, W); sonsuz pikseller +inf.
    """
    d = np.asarray(d, dtype=np.float32)
    h, w = d.shape
    W, H = int(hedef_wh[0]), int(hedef_wh[1])
    with np.errstate(divide="ignore"):
        inv = np.where(np.isfinite(d) & (d > 0), 1.0 / d, 0.0).astype(np.float32)
    bil = cv2.resize(inv, (W, H), interpolation=cv2.INTER_LINEAR)
    # en yakın komşu: kaynak piksel = floor((x+0.5)*w/W)
    xs = np.minimum(((np.arange(W) + 0.5) * w / W).astype(np.int64), w - 1)
    ys = np.minimum(((np.arange(H) + 0.5) * h / H).astype(np.int64), h - 1)
    yak = inv[np.ix_(ys, xs)]
    # bilinear destek hücresi (cv2 kuralı: xs = (x+0.5)*w/W - 0.5, sınırda kenetli)
    fx = (np.arange(W) + 0.5) * w / W - 0.5
    fy = (np.arange(H) + 0.5) * h / H - 0.5
    x0 = np.clip(np.floor(fx).astype(np.int64), 0, w - 2)
    y0 = np.clip(np.floor(fy).astype(np.int64), 0, h - 2)
    a, b = inv[:-1, :-1], inv[:-1, 1:]
    c, e = inv[1:, :-1], inv[1:, 1:]
    mx = np.maximum(np.maximum(a, b), np.maximum(c, e))
    mn = np.minimum(np.minimum(a, b), np.minimum(c, e))
    sonsuz_karisik = (mn == 0) & (mx > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        kenar = sonsuz_karisik | ((mn > 0) & (mx / np.where(mn > 0, mn, 1) > kenar_orani))
    k = kenar[np.ix_(y0, x0)]
    out_inv = np.where(k, yak, bil)
    with np.errstate(divide="ignore"):
        out = np.where(out_inv > 0, 1.0 / out_inv, np.inf).astype(np.float32)
    return out


def derinlik_ic_parametre(K_goruntu, goruntu_wh, derinlik_wh):
    """Derinlik haritasının iç parametreleri (aynı FOV, yarım-piksel kuralı)."""
    sx = derinlik_wh[0] / goruntu_wh[0]
    sy = derinlik_wh[1] / goruntu_wh[1]
    fx, fy, cx, cy = K_goruntu
    return [fx * sx, fy * sy, (cx + 0.5) * sx - 0.5, (cy + 0.5) * sy - 0.5]


DERINLIK_KIPLERI = ("piksel", "derinlik")


def derinlik_kipi_denetle(kip: str, seq: Path):
    if kip == "piksel":
        return None
    if kip == "derinlik":
        okuyucu = DerinlikOkuyucu(seq)
        if not okuyucu.var_mi():
            raise SystemExit(f"--derinlik_kipi derinlik: {okuyucu.dizin} yok")
        raise NotImplementedError(
            "Derinliğe bağlı (yol-integralli, çok katmanlı) kip henüz uygulanmadı "
            "(airsim_sim/TASARIM.md §5.2). Arayüz: DerinlikOkuyucu + derinlik_buyut.")
    raise ValueError(kip)


# =============================================================================
#  Parametreler
# =============================================================================
def parametreleri_kur(a, B_seq: float) -> tuple[TurbParams, dict]:
    """Yapılandırmanın `turbulence:` bölümü + CLI geçersiz kılmaları."""
    kw = {}
    yp = a.yapilandirma
    if yp and Path(yp).is_file():
        with open(yp) as f:
            t = (yaml.safe_load(f) or {}).get("turbulence") or {}
        known = set(TurbParams.__dataclass_fields__)
        kw = {k: v for k, v in t.items() if k in known}
        if kw.get("blur_grid") is not None:
            kw["blur_grid"] = tuple(kw["blur_grid"])
    kw.update(enabled=bool(a.D_over_r0 > 0), seed=int(a.tohum),
              D_over_r0=float(a.D_over_r0))
    for cli, alan in (("theta_c_over_fov", "theta_c_over_fov"),
                      ("tau0_over_dt", "tau0_over_dt"), ("blur", "blur"),
                      ("scint", "scint_sigma_I2"), ("frame_dt", "frame_dt"),
                      ("px_per_rad", "px_per_rad")):
        v = getattr(a, cli)
        if v is not None:
            kw[alan] = v
    if a.tilt is not None:
        kw["tilt"] = a.tilt == "acik"
    if a.psf_eski:
        kw["psf_legacy_scale"] = True
    # Stereo ayak izi ötelemesi
    cfg_B = float(kw.get("stereo_baseline", TurbParams.stereo_baseline))
    if a.stereo_taban in (None, "otomatik"):
        if abs(abs(B_seq) - cfg_B) <= 0.01 * cfg_B:
            B, neden = cfg_B, (f"yapilandirma degeri ({cfg_B}); sekans tabani "
                               f"{abs(B_seq):.6f} m ile %1 icinde (kanca esdegerligi)")
        else:
            B, neden = round(abs(B_seq), 4), (f"sekans tabani ({abs(B_seq):.6f} m, 4 haneye "
                                             f"yuvarlandi); yapilandirma {cfg_B} m'den farkli")
    elif a.stereo_taban == "sekans":
        B, neden = round(abs(B_seq), 4), "sekans tabani (CLI)"
    else:
        B, neden = float(a.stereo_taban), "CLI"
    kw["stereo_baseline"] = B
    return TurbParams(**kw), {"stereo_baseline_m": B, "neden": neden,
                              "yapilandirma": str(yp) if yp else None}


def cikti_adi(kaynak: str, p: TurbParams, araliklar, gurultu_dn: float = 0.0,
              gurultu_tohum: int | None = None, stereo_dis: str = "gomulu") -> str:
    ad = f"{kaynak}__{p.key()}"
    if stereo_dis == "sensor_yaml":    # yalnız açıkken: gömülü matriste ad eskisiyle aynı
        ad += "_sy"
    if p.psf_legacy_scale:
        ad += "_psfeski"
    if not p.tilt:
        ad += "_tiltsiz"
    if gurultu_dn > 0:                 # yalnız açıkken: kapalıyken ad eskisiyle aynı
        ad += f"_n{gurultu_dn:g}_ns{p.seed if gurultu_tohum is None else int(gurultu_tohum)}"
    if araliklar:
        ad += "__k" + "_".join(f"{a}-{b}" for a, b in araliklar)
    return ad


# =============================================================================
#  Kare işleme (işçi süreçlerinde de çalışır)
# =============================================================================
_IS = {}


def _isci_kur(durum: dict):
    cv2.setNumThreads(int(durum.get("cv_threads", 0)))
    cal = asl_oku(Path(durum["girdi"]))
    _IS.clear()
    _IS.update(durum)
    _IS["cal"] = cal
    _IS["rekt"] = Rektifikator(cal, durum["rekt_kip"], durum.get("stereo_dis", "gomulu"))
    p = TurbParams(**durum["params"])
    _IS["sim"] = TurbulenceSimulator(p, image_size=_IS["rekt"].size)


def _psf_ozet(psf: np.ndarray) -> tuple[float, float]:
    n = psf.shape[0]
    yy, xx = np.mgrid[0:n, 0:n] - n // 2
    s = float(psf.sum())
    cy = float((psf * yy).sum()) / s
    cx = float((psf * xx).sum()) / s
    var = float((psf * ((xx - cx) ** 2 + (yy - cy) ** 2)).sum()) / s / 2.0
    return float(np.sqrt(var)), float(psf.max())


def _kare_ozet(sim: TurbulenceSimulator, i: int, cam: int) -> list:
    """[du_ort, dv_ort, d_rms_eksen, d_maks, psf_sigma_px, psf_tepe] (kaba ızgara)."""
    if not sim.active:
        return [0.0] * 6
    du, dv = sim.displacement(i, cam, full=False)
    du = du.astype(np.float64); dv = dv.astype(np.float64)
    rms = float(np.sqrt(0.5 * np.mean(du ** 2 + dv ** 2)))
    dmax = float(np.sqrt(du ** 2 + dv ** 2).max())
    ps, pk = 0.0, 0.0
    if sim.p.blur == "speckle_global" and sim.blurfield is not None:
        X, Y = sim._footprint(i, cam)
        Hc, Wc = X.shape
        ps, pk = _psf_ozet(sim.blurfield.psf_at(float(X[Hc // 2, Wc // 2]),
                                                float(Y[Hc // 2, Wc // 2])))
    elif sim.p.blur == "gaussian":
        ps = float(sim.gauss_sigma_px)
    return [float(du.mean()), float(dv.mean()), rms, dmax, ps, pk]


def _kareleri_isle(gorev: list[int]) -> list[dict]:
    cal, rekt, sim = _IS["cal"], _IS["rekt"], _IS["sim"]
    src = Path(_IS["girdi"]) / "mav0"
    dst = Path(_IS["hedef"]) / "mav0"
    png = [cv2.IMWRITE_PNG_COMPRESSION, int(_IS["png_sikistirma"])]
    alan_her = int(_IS.get("alan_her", 0))
    out = []
    for i in gorev:
        ts = cal["ts"][i]
        n0, n1 = cal["f0"].get(ts), cal["f1"].get(ts)
        if n0 is None or n1 is None:
            out.append({"i": i, "ts": ts, "atlandi": "cam1 damgasi yok"})
            continue
        i0 = cv2.imread(str(src / "cam0/data" / n0), cv2.IMREAD_GRAYSCALE)
        i1 = cv2.imread(str(src / "cam1/data" / n1), cv2.IMREAD_GRAYSCALE)
        if i0 is None or i1 is None:
            out.append({"i": i, "ts": ts, "atlandi": "okunamadi"})
            continue
        r0, r1 = rekt.uygula(i0, i1)
        d0 = sim.apply(r0, i, cam=0)
        d1 = sim.apply(r1, i, cam=1)
        g_dn = float(_IS.get("gur_dn", 0.0) or 0.0)
        if g_dn > 0:                   # sensör gürültüsü: türbülanstan SONRA (kanca ile aynı sıra)
            d0 = sensor_noise(d0, g_dn, _IS["gur_tohum"], i, 0)
            d1 = sensor_noise(d1, g_dn, _IS["gur_tohum"], i, 1)
        ad = f"{ts}.png"
        cv2.imwrite(str(dst / "cam0/data" / ad), d0, png)
        cv2.imwrite(str(dst / "cam1/data" / ad), d1, png)
        rec = {"i": i, "ts": ts, "ad": ad, "c0": _kare_ozet(sim, i, 0),
               "c1": _kare_ozet(sim, i, 1)}
        if alan_her and sim.active and (i % alan_her == 0):
            a0 = np.stack(sim.displacement(i, 0, full=False)).astype(np.float32)
            a1 = np.stack(sim.displacement(i, 1, full=False)).astype(np.float32)
            np.save(Path(_IS["hedef"]) / "turb_alanlar" / f"{i:06d}.npy",
                    np.stack([a0, a1]))
        out.append(rec)
    return out


# =============================================================================
#  Ana üretim
# =============================================================================
def _kopyala_yardimci(src_seq: Path, dst_seq: Path, ts_secili: list[int],
                      derinlik_kopyala: bool):
    """GT/IMU vb. kamera dışı her şeyi kopyala; depth0'dan yalnız seçili kareler."""
    for f in sorted(src_seq.iterdir()):
        if f.is_file():
            shutil.copyfile(f, dst_seq / f.name)
    sm, dm = src_seq / "mav0", dst_seq / "mav0"
    for e in sorted(sm.iterdir()):
        if e.name in ("cam0", "cam1", "depth0"):
            continue
        if e.is_dir():
            shutil.copytree(e, dm / e.name, copy_function=shutil.copyfile)
        else:
            shutil.copyfile(e, dm / e.name)
    dep = sm / "depth0"
    if derinlik_kopyala and dep.is_dir():
        (dm / "depth0/data").mkdir(parents=True, exist_ok=True)
        for f in sorted(dep.iterdir()):
            if f.is_file():
                shutil.copyfile(f, dm / "depth0" / f.name)
        for ts in ts_secili:
            s = dep / "data" / f"{ts}.npy"
            if s.exists():
                shutil.copyfile(s, dm / "depth0/data" / s.name)


def uret(a) -> Path:
    t_bas = time.time()
    girdi = Path(a.girdi).expanduser().resolve()
    cal = asl_oku(girdi)
    n = len(cal["ts"])
    idx, araliklar = kareleri_ayristir(a.kareler, n)
    tum = (araliklar == [(0, n)])
    derinlik_kipi_denetle(a.derinlik_kipi, girdi)
    rekt = Rektifikator(cal, a.rektifikasyon, a.stereo_dis)
    p, taban_karari = parametreleri_kur(a, rekt.baseline)
    gur_dn = float(a.gurultu_dn or 0.0)
    gur_tohum = int(p.seed if a.gurultu_tohum is None else a.gurultu_tohum)
    ad = a.ad or cikti_adi(girdi.name, p, None if tum else araliklar, gur_dn, a.gurultu_tohum,
                           a.stereo_dis if "T_cn_cnm1" not in cal["cam1"] else "gomulu")
    hedef = Path(a.cikti_kok).expanduser() / ad
    if hedef.exists():
        if not a.uzerine_yaz:
            raise SystemExit(f"çıktı var: {hedef} (--uzerine_yaz)")
        shutil.rmtree(hedef)
    gecici = hedef.with_name(hedef.name + ".yaziliyor")
    if gecici.exists():
        shutil.rmtree(gecici)
    for s in ("mav0/cam0/data", "mav0/cam1/data"):
        (gecici / s).mkdir(parents=True)
    if a.alan_her:
        (gecici / "turb_alanlar").mkdir()

    durum = {"girdi": str(girdi), "hedef": str(gecici), "rekt_kip": rekt.kip,
             "params": asdict(p), "png_sikistirma": a.png_sikistirma,
             "alan_her": a.alan_her, "cv_threads": 0 if a.isci <= 1 else 1}
    if gur_dn > 0:                     # yalnız açıkken: kapalıyken işçi durumu eskisiyle aynı
        durum.update(gur_dn=gur_dn, gur_tohum=gur_tohum)
    if a.stereo_dis != "gomulu":
        durum.update(stereo_dis=a.stereo_dis)
    t_kur = time.time()
    kayitlar = []
    if a.isci <= 1:
        _isci_kur(durum)
        t_kur = time.time() - t_kur
        kayitlar = _kareleri_isle(idx)
        sim = _IS["sim"]
    else:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        parca = max(1, int(np.ceil(len(idx) / (a.isci * 4))))
        gorevler = [idx[k:k + parca] for k in range(0, len(idx), parca)]
        with ProcessPoolExecutor(a.isci, mp_context=mp.get_context("spawn"),
                                 initializer=_isci_kur, initargs=(durum,)) as ex:
            for r in ex.map(_kareleri_isle, gorevler):
                kayitlar.extend(r)
        t_kur = float("nan")
        sim = TurbulenceSimulator(p, image_size=rekt.size)   # meta bilgisi için
    kayitlar.sort(key=lambda r: r["i"])
    yazilan = [r for r in kayitlar if "ad" in r]
    atlanan = [{"i": r["i"], "ts": r["ts"], "neden": r["atlandi"]}
               for r in kayitlar if "ad" not in r]

    # --- data.csv, zaman dosyaları ---
    for cam in (0, 1):
        satir = ["#timestamp [ns],filename"] + [f"{r['ts']},{r['ad']}" for r in yazilan]
        (gecici / f"mav0/cam{cam}/data.csv").write_text("\n".join(satir) + "\n")
    (gecici / "zamanlar.txt").write_text("".join(f"{r['ts']}\n" for r in yazilan))
    if not tum:
        for (x, y) in araliklar:
            (gecici / f"zamanlar_k{x:05d}-{y:05d}.txt").write_text(
                "".join(f"{r['ts']}\n" for r in yazilan if x <= r["i"] < y))

    # --- kalibrasyon ---
    if not rekt.sensor_yaml_yaz(gecici / "mav0", ad):
        for cam in (0, 1):
            shutil.copyfile(girdi / f"mav0/cam{cam}/sensor.yaml",
                            gecici / f"mav0/cam{cam}/sensor.yaml")
    _kopyala_yardimci(girdi, gecici, [r["ts"] for r in yazilan], a.derinlik_kopyala)

    # --- meta ---
    s_or = optics.sigma_tilt_px(p.D_over_r0, p.wavelength, p.aperture_D, p.px_per_rad) \
        if p.D_over_r0 > 0 else 0.0
    sutun = ["du_ort_px", "dv_ort_px", "d_rms_eksen_px", "d_maks_px", "psf_sigma_px",
             "psf_tepe"]
    meta = {
        "bicim": META_SURUM,
        "kaynak": {"yol": str(girdi), "ad": girdi.name, "n_kare": n,
                   "cam0_sensor_sha256": sha256_dosya(girdi / "mav0/cam0/sensor.yaml"),
                   "cam1_sensor_sha256": sha256_dosya(girdi / "mav0/cam1/sensor.yaml")},
        "parametreler": asdict(p),
        "parametre_anahtari": p.key(),
        "parametre_hash": p.hash(),
        "stereo_taban_karari": taban_karari,
        "derinlik_kipi": a.derinlik_kipi,
        "rektifikasyon": rekt.bilgi(),
        "kareler": {"belirtim": a.kareler, "araliklar": araliklar,
                    "indeks_tanimi": "cam0 damgalarinin sayisal sirali listesinde mutlak "
                                     "indeks (kanca _seq_frame_idx); yari acik [a,b)",
                    "n_yazilan": len(yazilan), "atlanan": atlanan},
        "oracle": {
            "sigma_px": s_or,
            "sigma_px_tanimi": "optics.sigma_tilt_px(D/r0, lam, D, px_per_rad) = "
                               "backend.ArmConfig.sigma_px_oracle (A1/A5)",
            "A1_ind_over_tilt": A1_IND_OVER_TILT,
            "A1_sigma_toplam_px": float(np.hypot(s_or, A1_IND_OVER_TILT * s_or)),
            "theta_c_px": float(p.theta_c_over_fov * rekt.size[0]),
            "simulator_info": sim.info(),
        },
        "kare_bilgisi": {
            "aciklama": "cam0/cam1 icin kaba izgara (warp_grid_step) yer degistirme "
                        "alani ozeti [px] ve merkez PSF (speckle_global) ozeti",
            "sutunlar": sutun,
            "i": [r["i"] for r in yazilan], "ts": [r["ts"] for r in yazilan],
            "cam0": [r["c0"] for r in yazilan], "cam1": [r["c1"] for r in yazilan]},
        "alanlar": ({"dizin": "turb_alanlar", "her": a.alan_her,
                     "bicim": "[cam, (du,dv), Hc, Wc] float32 kaba izgara; tam cozunurluk "
                              "= cv2.resize(INTER_LINEAR) (simulator.displacement)"}
                    if a.alan_her else None),
        "cikti": {"png_sikistirma": a.png_sikistirma, "goruntu_boyutu": list(rekt.size),
                  "dtype": "uint8 gri (kanca da uint8 calisir: nicemleme farki 0)"},
        "kod_parmak_izi": {
            "dosyalar": {f: sha256_dosya(MODUL_DIZINI / f) for f in
                         (PARMAK_IZI_DOSYALARI + (("hooks.py",) if gur_dn > 0 else ()))},
            "opencv": cv2.__version__, "numpy": np.__version__,
            "scipy": __import__("scipy").__version__},
    }
    if gur_dn > 0:                     # yalnız açıkken: kapalıyken meta eskisiyle aynı alanları taşır
        meta["sensor_gurultusu"] = {
            "sigma_dn": gur_dn, "tohum": gur_tohum,
            "tohum_kaynagi": "CLI --gurultu_tohum" if a.gurultu_tohum is not None
                             else "turbulans tohumu (kanca: noise_seed None -> params.seed)",
            "fonksiyon": "turb/hooks.py::sensor_noise (PyEco kancasiyla ayni)",
            "tanim": "turbulanstan SONRA; float32 + N(0, sigma), rint, [0,255], uint8; "
                     "default_rng([tohum, mutlak kare indeksi, kamera])"}
    (gecici / "turb_meta.json").write_text(json.dumps(meta, indent=1, ensure_ascii=False))
    gecici.rename(hedef)
    sure = time.time() - t_bas
    boyut = sum(f.stat().st_size for f in hedef.rglob("*") if f.is_file())
    log = {"uretim_zamani_unix": t_bas, "sure_s": sure, "simulator_kurulum_s": t_kur,
           "isci": a.isci, "boyut_bayt": boyut, "n_kare": len(yazilan),
           "komut": sys.argv}
    (hedef / "turb_uretim.json").write_text(json.dumps(log, indent=1))
    print(f"[uret_sekans] {hedef}  kare={len(yazilan)}  atlanan={len(atlanan)}  "
          f"sure={sure:.1f}s  boyut={boyut / 1e6:.1f} MB  rekt={rekt.kip}  "
          f"sigma_or={s_or:.3f}px  B_turb={p.stereo_baseline}")
    return hedef


def arg_ayristirici() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Türbülanslı EuRoC-ASL sekans üreticisi")
    ap.add_argument("--girdi", required=True, help="kaynak ASL sekans kökü (mav0'ın üstü)")
    ap.add_argument("--cikti_kok", default="/data/sdata_turb")
    ap.add_argument("--ad", default=None, help="çıktı klasör adı (vars.: otomatik)")
    ap.add_argument("--D_over_r0", type=float, required=True)
    ap.add_argument("--tohum", type=int, default=0)
    ap.add_argument("--theta_c_over_fov", type=float, default=None)
    ap.add_argument("--tau0_over_dt", type=float, default=None)
    ap.add_argument("--blur", default=None,
                    choices=["none", "gaussian", "speckle_global", "speckle_grid"])
    ap.add_argument("--tilt", default=None, choices=["acik", "kapali"])
    ap.add_argument("--scint", type=float, default=None, help="sintilasyon sigma_I^2")
    ap.add_argument("--psf_eski", action="store_true",
                    help="psf_legacy_scale=True (E1a-E9 kanıtlarının eski PSF ölçeği)")
    ap.add_argument("--frame_dt", type=float, default=None)
    ap.add_argument("--px_per_rad", type=float, default=None)
    ap.add_argument("--stereo_taban", default="otomatik",
                    help="otomatik | sekans | <metre> (cam1 ayak izi ötelemesi)")
    ap.add_argument("--kareler", default=None,
                    help="yarı açık aralıklar, ör. 400:900 veya 400:900,800:1300")
    ap.add_argument("--rektifikasyon", default="otomatik",
                    choices=["otomatik", "rektifiye", "kimlik", "ham"])
    ap.add_argument("--derinlik_kipi", default="piksel", choices=DERINLIK_KIPLERI,
                    help="piksel: mevcut piksel-eşdeğer kip; derinlik: (henüz yok)")
    ap.add_argument("--yapilandirma", default=str(VARSAYILAN_YAPILANDIRMA),
                    help="turbulence: bölümü taban alınacak yaml (boş = TurbParams vars.)")
    ap.add_argument("--isci", type=int, default=1, help="paralel süreç sayısı")
    ap.add_argument("--png_sikistirma", type=int, default=3)
    ap.add_argument("--alan_her", type=int, default=0,
                    help=">0: her K. karede kaba yer değiştirme alanını kaydet")
    ap.add_argument("--derinlik_kopyala", type=int, default=1)
    ap.add_argument("--uzerine_yaz", action="store_true")
    ap.add_argument("--stereo_dis", default="gomulu", choices=["gomulu", "sensor_yaml"],
                    help="cam0->cam1 dış parametresi (sensor.yaml'da T_cn_cnm1 yoksa): gomulu = eski kanca "
                         "sabit matrisi (kod/stereo_front_end.py varsayılanı); sensor_yaml = inv(T_B_c1) T_B_c0, "
                         "pyeco_turb/config.yaml frontend.stereo_extrinsics ile AYNI (Faz 2 için zorunlu)")
    ap.add_argument("--gurultu_dn", type=float, default=0.0,
                    help="sensör gürültüsü sigma [DN], türbülanstan sonra (0 = kapalı)")
    ap.add_argument("--gurultu_tohum", type=int, default=None,
                    help="gürültü tohumu (vars.: türbülans tohumu = kanca davranışı)")
    return ap


def main(argv=None):
    return uret(arg_ayristirici().parse_args(argv))


if __name__ == "__main__":
    main()
