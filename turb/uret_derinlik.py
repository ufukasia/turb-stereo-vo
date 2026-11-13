"""
Derinliğe bağlı türbülanslı sekans üreticisi (AirSim ASL kayıtları için).

    python -m turb.uret_derinlik --girdi /data/sdata_airsim/AF_F1_h5 \
        --cikti_kok /data/sdata_turb_airsim --C0 1e-12 --tohum 0 \
        --dtm "../airsim_sim/africa/harita_africa.npz" --D 4.8e-3 --isci 12

Çıktı: <cikti_kok>/<seq>__C0<..>_s<tohum>/mav0/... (cam0/cam1 türbülanslı; GT, sensor.yaml,
data.csv kopya; depth0 sembolik bağ) + turb_meta.json (parametreler, kod özeti, kare başına
derinlik dilimli istatistik) + ornek_alanlar.npz (her 100. karenin düşük çözünürlüklü alanları).
C0 = 0: görüntüler bayt bayt kopyalanır (kontrol kolu).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import sys
import time
from dataclasses import asdict
from pathlib import Path

import cv2
import numpy as np
import yaml

from .derinlik_model import Arazi, DerinlikModeli, DerinlikParam, uygula

DERIN_DILIM = [0, 25, 50, 100, 200, 400, 1e9]
_G = {}


def _R(q):
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def _rotvec(R):
    c = np.clip((np.trace(R) - 1) / 2, -1, 1); a = math.acos(c)
    if a < 1e-9:
        return np.zeros(3)
    v = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    return v * a / (2 * math.sin(a))


def sekans_oku(seq: Path):
    m = seq / "mav0"
    c0 = yaml.safe_load(open(m / "cam0/sensor.yaml")); c1 = yaml.safe_load(open(m / "cam1/sensor.yaml"))
    fx, fy, cx, cy = c0["intrinsics"]; W, H = c0["resolution"]
    B = -float(np.array(c1["T_cn_cnm1"]).reshape(4, 4)[0, 3])
    ts = [l.split(",")[0].strip() for l in open(m / "cam0/data.csv") if not l.startswith("#")]
    gt = np.loadtxt(m / "state_groundtruth_estimate0/data.csv", delimiter=",", comments="#")
    gts = [l.split(",")[0].strip() for l in open(m / "state_groundtruth_estimate0/data.csv") if not l.startswith("#")]
    assert len(gts) == len(ts) and all(int(t) == int(g) for t, g in zip(ts, gts)), "GT/kare damgası uyuşmuyor"
    return dict(m=m, fx=fx, fy=fy, cx=cx, cy=cy, W=int(W), H=int(H), B=B, ts=ts, gt=gt)


def ofsetler(model: DerinlikModeli, gt: np.ndarray, dt: float):
    """Donmuş akış: kare başına ekran kayması (M, 2) ve A/B karışım açısı (M,)."""
    p = model.p
    n, M = len(gt), p.M
    of = np.zeros((n, M, 2)); km = np.zeros((n, M))
    kal = np.diff(model.sinir)
    yon = math.radians(p.ruzgar_yon_deg)
    vw = p.ruzgar * np.array([math.cos(yon), math.sin(yon), 0.0])
    for i in range(1, n):
        R0, R1 = _R(gt[i - 1, 4:8]), _R(gt[i, 4:8])
        w = _rotvec(R0.T @ R1)
        dtc = R0.T @ (gt[i, 1:4] - gt[i - 1, 1:4])
        vc = R0.T @ vw
        d = model.d_ekran
        of[i, :, 0] = of[i - 1, :, 0] + d * w[1] + dtc[0] - vc[0] * dt
        of[i, :, 1] = of[i - 1, :, 1] - d * w[0] + dtc[1] - vc[1] * dt
        km[i] = km[i - 1] + (abs(dtc[2]) + abs(vc[2]) * dt) / kal * (math.pi / 2)
    return of, km


def _egim(fx, fy, cx, cy, w, h, W, H):
    """Derinlik ızgarası (w×h) piksel merkezlerinin eğimleri (görüntü iç parametreleriyle)."""
    sx, sy = w / W, h / H
    fxl, fyl = fx * sx, fy * sy
    cxl, cyl = (cx + 0.5) * sx - 0.5, (cy + 0.5) * sy - 0.5
    u, v = np.meshgrid(np.arange(w, dtype=np.float64), np.arange(h, dtype=np.float64))
    return (u - cxl) / fxl, (v - cyl) / fyl, fxl


def _derinlik(path):
    d = np.load(path).astype(np.float32)
    d[~np.isfinite(d) | (d >= 6.0e4) | (d <= 0)] = np.inf
    return d


def _cam1_derinlik(Z0, fxl, B):
    """cam1 pikselinin nesne derinliği: Z1(u) = Z0(u + fx·B/Z) (2 sabit-nokta yinelemesi)."""
    h, w = Z0.shape
    u = np.arange(w)[None, :].repeat(h, 0); v = np.arange(h)[:, None].repeat(w, 1)
    Z = Z0.copy()
    for _ in range(2):
        disp = np.where(np.isfinite(Z), fxl * B / Z, 0.0)
        uu = np.clip(np.rint(u + disp).astype(int), 0, w - 1)
        Z = Z0[v, uu]
    return Z


def _istat(Z, tani, blur, dx, dy):
    out = {}
    ok = np.isfinite(Z)
    for a, b in zip(DERIN_DILIM[:-1], DERIN_DILIM[1:]):
        m = ok & (Z >= a) & (Z < b)
        if m.sum() < 20:
            continue
        out[f"{a:g}-{b:g}"] = {"n": int(m.sum()),
                               "sigma_tot_px": float(np.median(tani["sigma_tot_px"][m])),
                               "sigma_lo_px": float(np.median(tani["sigma_lo_px"][m])),
                               "disp_rms_px": float(np.sqrt(np.mean(dx[m] ** 2 + dy[m] ** 2) / 2)),
                               "blur_px": float(np.median(blur[m])),
                               "D_r0": float(np.median(tani["D_r0"][m]))}
    return out


def gurultu_ekle(img, sigma, tohum, i, cam):
    """Sensör gürültüsü (Gauss, sigma DN). Tohum (gürültü tohumu, kare, kamera) ile belirlenimli."""
    if sigma <= 0:
        return img
    rng = np.random.default_rng([int(tohum), int(i), int(cam)])
    return np.clip(np.rint(img.astype(np.float32) + rng.normal(0.0, sigma, img.shape)), 0, 255).astype(np.uint8)


def _kare_isle(i):
    g = _G
    S, model, of, km = g["S"], g["model"], g["of"], g["km"]
    ts = S["ts"][i]
    Z0 = _derinlik(S["m"] / "depth0/data" / f"{ts}.npy")
    h, w = Z0.shape
    ex, ey, fxl = _egim(S["fx"], S["fy"], S["cx"], S["cy"], w, h, S["W"], S["H"])
    R = _R(S["gt"][i, 4:8]); t = S["gt"][i, 1:4]
    sonuc = {"i": i}
    for cam in (0, 1):
        Z = Z0 if cam == 0 else _cam1_derinlik(Z0, fxl, S["B"])
        dx, dy, blur, tani = model.kare(Z, ex, ey, R, t, of[i], km[i], ofset_x=(0.0 if cam == 0 else S["B"]))
        img = cv2.imread(str(S["m"] / f"cam{cam}/data/{ts}.png"), cv2.IMREAD_UNCHANGED)
        cv2.imwrite(str(g["out"] / f"mav0/cam{cam}/data/{ts}.png"),
                    gurultu_ekle(uygula(img, dx, dy, blur), g["gur"], g["gur_t"], i, cam))
        if cam == 0:
            sonuc["istat"] = _istat(Z, tani, blur, dx, dy)
            if i % 100 == 0:
                sonuc["alan"] = {"dx": dx.astype(np.float32), "dy": dy.astype(np.float32),
                                 "blur": blur.astype(np.float32),
                                 "sigma_tot": tani["sigma_tot_px"].astype(np.float32), "Z": Z0}
        else:
            sonuc["dx1"] = float(np.sqrt(np.mean(dx[np.isfinite(Z)] ** 2)))
            if i % 100 == 0:
                sonuc["alan"]["dx1"] = dx.astype(np.float32)
    return sonuc


def kod_ozeti() -> str:
    h = hashlib.sha1()
    for f in ("derinlik_model.py", "uret_derinlik.py", "screen.py"):
        h.update(open(Path(__file__).parent / f, "rb").read())
    return h.hexdigest()[:12]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--girdi", required=True)
    ap.add_argument("--cikti_kok", required=True)
    ap.add_argument("--C0", type=float, required=True)
    ap.add_argument("--tohum", type=int, default=0)
    ap.add_argument("--dtm", default=None, help="harita.py npz; yoksa düz zemin (--z_zemin)")
    ap.add_argument("--z_zemin", type=float, default=0.0)
    ap.add_argument("--D", type=float, default=4.8e-3)
    ap.add_argument("--L0", type=float, default=5.0)
    ap.add_argument("--ruzgar", type=float, default=2.0)
    ap.add_argument("--warp_px", type=float, default=4.0)
    ap.add_argument("--kareler", default=None, help="a:b (yalnız bu aralık; GT/data.csv de kırpılır)")
    ap.add_argument("--isci", type=int, default=12)
    ap.add_argument("--ad", default=None)
    ap.add_argument("--gurultu", type=float, default=0.0, help="sensör gürültüsü sigma (DN); 0 = yok")
    ap.add_argument("--sabit_kayma", type=float, default=0.0,
                    help="KONTROL: türbülans yok; iki kameraya aynı sabit alt-piksel kayma (px, x ve y) uygulanır. "
                         "Yalnız ara-değerleme yumuşatmasının etkisini ölçer. C0=0 ile kullanılır.")
    ap.add_argument("--gurultu_tohum", type=int, default=None, help="varsayılan: --tohum")
    a = ap.parse_args(argv)

    seq = Path(a.girdi); S = sekans_oku(seq)
    n = len(S["ts"])
    i0, i1 = (0, n) if not a.kareler else tuple(int(x) for x in a.kareler.split(":"))
    p = DerinlikParam(C0=a.C0, tohum=a.tohum, D=a.D, L0=a.L0, ruzgar=a.ruzgar, warp_px=a.warp_px)
    gur_t = a.tohum if a.gurultu_tohum is None else a.gurultu_tohum
    ad = a.ad or (f"{seq.name}__{p.anahtar()}" + (f"_g{a.gurultu:g}t{gur_t}" if a.gurultu > 0 else "")
                  + (f"_ks{a.sabit_kayma:g}" if a.sabit_kayma > 0 else "")
                  + (f"_k{i0}-{i1}" if a.kareler else ""))
    out = Path(a.cikti_kok) / ad
    if out.exists():
        shutil.rmtree(out)
    for c in ("cam0/data", "cam1/data", "state_groundtruth_estimate0"):
        (out / "mav0" / c).mkdir(parents=True, exist_ok=True)
    sec = list(range(i0, i1))
    hdr = "#timestamp [ns],filename\n"
    for c in ("cam0", "cam1"):
        shutil.copy(S["m"] / c / "sensor.yaml", out / "mav0" / c / "sensor.yaml")
        open(out / "mav0" / c / "data.csv", "w").write(hdr + "".join(f"{S['ts'][i]},{S['ts'][i]}.png\n" for i in sec))
    gl = [l for l in open(S["m"] / "state_groundtruth_estimate0/data.csv")]
    gh = [l for l in gl if l.startswith("#")]; gd = [l for l in gl if not l.startswith("#")]
    open(out / "mav0/state_groundtruth_estimate0/data.csv", "w").write("".join(gh + [gd[i] for i in sec]))
    if (S["m"] / "depth0").is_dir():
        os.symlink(S["m"] / "depth0", out / "mav0/depth0")
    t0 = time.time()
    model = DerinlikModeli(p, S["fx"], Arazi(a.dtm, a.z_zemin))
    istat, ornek, dx1 = {}, {}, {}
    if a.sabit_kayma > 0:
        assert a.C0 == 0, "--sabit_kayma yalnız C0=0 ile"
        F = np.full((4, 4), a.sabit_kayma, np.float32); Z0 = np.zeros((4, 4), np.float32)
        for i in sec:
            for k, c in enumerate(("cam0", "cam1")):
                img = cv2.imread(str(S["m"] / f"{c}/data/{S['ts'][i]}.png"), cv2.IMREAD_UNCHANGED)
                cv2.imwrite(str(out / f"mav0/{c}/data/{S['ts'][i]}.png"),
                            gurultu_ekle(uygula(img, F, F, Z0), a.gurultu, gur_t, i, k))
    elif a.C0 == 0 and a.gurultu > 0:
        for i in sec:
            for k, c in enumerate(("cam0", "cam1")):
                img = cv2.imread(str(S["m"] / f"{c}/data/{S['ts'][i]}.png"), cv2.IMREAD_UNCHANGED)
                cv2.imwrite(str(out / f"mav0/{c}/data/{S['ts'][i]}.png"), gurultu_ekle(img, a.gurultu, gur_t, i, k))
    elif a.C0 == 0:
        for i in sec:
            for c in ("cam0", "cam1"):
                shutil.copyfile(S["m"] / f"{c}/data/{S['ts'][i]}.png", out / f"mav0/{c}/data/{S['ts'][i]}.png")
    else:
        model._ekran_kur()
        dt = (int(S["ts"][1]) - int(S["ts"][0])) * 1e-9
        of, km = ofsetler(model, S["gt"], dt)
        _G.update(S=S, model=model, of=of, km=km, out=out, gur=a.gurultu, gur_t=gur_t)
        import multiprocessing as mp
        with mp.get_context("fork").Pool(a.isci) as pool:
            for r in pool.imap_unordered(_kare_isle, sec, chunksize=4):
                istat[r["i"]] = r["istat"]; dx1[r["i"]] = r["dx1"]
                if "alan" in r:
                    ornek[r["i"]] = r["alan"]
        if ornek:
            np.savez_compressed(out / "ornek_alanlar.npz",
                                **{f"{k}_{i}": v for i, al in ornek.items() for k, v in al.items()})
    meta = {"kaynak": str(seq), "kareler": [i0, i1], "model": model.info(), "kod_ozeti": kod_ozeti(),
            "dtm": a.dtm, "gurultu_dn": a.gurultu, "sabit_kayma_px": a.sabit_kayma, "gurultu_tohum": gur_t, "sure_s": time.time() - t0,
            "kare_istat": {str(i): istat[i] for i in sorted(istat)},
            "cam1_disp_rms_px": {str(i): dx1[i] for i in sorted(dx1)}}
    json.dump(meta, open(out / "turb_meta.json", "w"), indent=1)
    print(json.dumps({"cikti": str(out), "sure_s": round(meta["sure_s"], 1), "kod": meta["kod_ozeti"]}))
    return out


if __name__ == "__main__":
    main()
