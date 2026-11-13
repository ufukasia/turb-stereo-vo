"""
Rakip hattı ön-test kapıları G0a–G0e (ON_TEST_RAKIPLER.md §1). Yalnız denetim; hiçbir koşum sonucunu değiştirmez.

    python -m turb.rakip_g0 --once <değişiklik öncesi çıktı kökü> --is <scratch> --out <json>

G0a  gürültü kapalı: değişiklik sonrası çıktı == değişiklik öncesi (meta yalnız uret_sekans.py parmak izinde farklı)
G0b  gürültü açık: PNG == hooks.sensor_noise(σ=0 çıktısı) == TurbulenceInjector.apply_stereo(rektifiye temiz)
G0c  TurbParams(üretici) == PyEco Faz 2 hücresinin etkin turbulence bölümü (D/r0 {0,8,16,24} x tohum 0-2)
G0d  belirlenim: aynı komut iki kez ve --isci 4 == --isci 1 (görüntü sha256)
G0e  tüketici rektifikasyonu (ORB-SLAM3 / DROID / AirSLAM kuralları) üretilen rektifiye sensor.yaml'da kimlik
G0f  (2026-11-08, sapma) --stereo_dis sensor_yaml: σ=0, D/r0=0 çıktısı == PyEco ön-ucunun (pyeco_turb
     stereo_front_end.py, stereo_extrinsics: sensor_yaml) rectify_images çıktısı bayt-aynı; σ=1, D/r0 0 ve 8 ==
     TurbulenceInjector.apply_stereo(PyEco rektifiye kare); dikey kayma (SIFT eşleşmeleri) gömülü ↔ sensor_yaml
G0b/G0e "sy" sürümleri: aynı denetimler sensor_yaml çıktısında
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import cv2
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
PY = sys.executable
SEQ = Path("/data/sdata/MH02")
KARE = "600:620"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def ozet(d: Path, haric=("turb_uretim.json", "turb_meta.json")) -> dict:
    return {str(f.relative_to(d)): sha(f) for f in sorted(d.rglob("*"))
            if f.is_file() and f.name not in haric}


def uret(kok: Path, dr: float, gur: float = 0.0, isci: int = 1, ek=(), sy: bool = False) -> Path:
    cmd = [PY, "-m", "turb.uret_sekans", "--girdi", str(SEQ), "--cikti_kok", str(kok),
           "--D_over_r0", str(dr), "--tohum", "0", "--kareler", KARE, "--isci", str(isci), *ek]
    if gur > 0:
        cmd += ["--gurultu_dn", str(gur)]
    if sy:
        cmd += ["--stereo_dis", "sensor_yaml"]
    out = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                         env={**os.environ, "OMP_NUM_THREADS": "4"})
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-2000:])
    line = [l for l in out.stdout.splitlines() if l.startswith("[uret_sekans]")][-1]
    return Path(line.split()[1])


def meta_fark(a: dict, b: dict, yol="") -> list:
    out = []
    for k in sorted(set(a) | set(b)):
        y = f"{yol}.{k}" if yol else k
        if k not in a or k not in b:
            out.append(y)
        elif isinstance(a[k], dict) and isinstance(b[k], dict):
            out += meta_fark(a[k], b[k], y)
        elif a[k] != b[k]:
            out.append(y)
    return out


def g0a(once: Path, kok: Path) -> dict:
    res = {}
    for dr in (0, 8):
        eski = next(once.glob(f"MH02__dr{dr}_*__k600-620"))
        yeni = uret(kok / "g0a", dr)
        oe, oy = ozet(eski), ozet(yeni)
        farkli = [k for k in sorted(set(oe) | set(oy)) if oe.get(k) != oy.get(k)]
        me = json.loads((eski / "turb_meta.json").read_text())
        my = json.loads((yeni / "turb_meta.json").read_text())
        mf = meta_fark(me, my)
        res[f"dr{dr}"] = {"ad_ayni": eski.name == yeni.name, "n_dosya": len(oe), "farkli_dosya": farkli,
                          "meta_farkli_alanlar": mf,
                          "gecti": eski.name == yeni.name and not farkli
                          and mf == ["kod_parmak_izi.dosyalar.uret_sekans.py"]}
    res["gecti"] = all(v["gecti"] for v in res.values())
    return res


def _cal_rect(seq: Path):
    from turb.uret_sekans import asl_oku, Rektifikator
    cal = asl_oku(seq)
    return cal, Rektifikator(cal, "rektifiye")


def g0b(kok: Path) -> dict:
    from turb.hooks import sensor_noise, TurbulenceInjector
    from turb.uret_sekans import arg_ayristirici, parametreleri_kur
    cal, rekt = _cal_rect(SEQ)
    res = {}
    for dr in (0, 8):
        temiz = next((kok / "g0a").glob(f"MH02__dr{dr}_*__k600-620"))     # σ=0 çıktısı (değişiklik sonrası)
        gur = uret(kok / "g0b", dr, 1.0)
        a = arg_ayristirici().parse_args(["--girdi", str(SEQ), "--D_over_r0", str(dr), "--tohum", "0"])
        p, _ = parametreleri_kur(a, rekt.baseline)
        inj = TurbulenceInjector(p, sensor_noise_dn=1.0, noise_seed=None)
        n_kare, n_fark_ref, n_fark_kanca, farklar = 0, 0, 0, []
        c01 = []
        for i in range(600, 620):
            ts = cal["ts"][i]
            g = [cv2.imread(str(gur / f"mav0/cam{c}/data/{ts}.png"), cv2.IMREAD_UNCHANGED) for c in (0, 1)]
            t = [cv2.imread(str(temiz / f"mav0/cam{c}/data/{ts}.png"), cv2.IMREAD_UNCHANGED) for c in (0, 1)]
            ref = [sensor_noise(t[c], 1.0, 0, i, c) for c in (0, 1)]
            i0 = cv2.imread(str(SEQ / "mav0/cam0/data" / cal["f0"][ts]), cv2.IMREAD_GRAYSCALE)
            i1 = cv2.imread(str(SEQ / "mav0/cam1/data" / cal["f1"][ts]), cv2.IMREAD_GRAYSCALE)
            r0, r1 = rekt.uygula(i0, i1)
            k0, k1 = inj.apply_stereo(r0, r1, i)
            for c, kk in ((0, k0), (1, k1)):
                n_fark_ref += int(np.count_nonzero(g[c] != ref[c]))
                n_fark_kanca += int(np.count_nonzero(g[c] != kk))
            d = [g[c].astype(np.float64) - t[c].astype(np.float64) for c in (0, 1)]
            m = [(t[c] > 3) & (t[c] < 252) for c in (0, 1)]
            farklar.append(d[0][m[0]])
            mm = m[0] & m[1]
            c01.append(np.corrcoef(d[0][mm], d[1][mm])[0, 1])
            n_kare += 1
        f = np.concatenate(farklar)
        res[f"dr{dr}"] = {"ad": gur.name, "n_kare": n_kare, "piksel_farki_sensor_noise": n_fark_ref,
                          "piksel_farki_kanca": n_fark_kanca, "std_dn": float(f.std()),
                          "ort_dn": float(f.mean()), "cam01_korelasyon": float(np.mean(c01)),
                          "gecti": n_fark_ref == 0 and n_fark_kanca == 0 and 0.95 <= f.std() <= 1.10}
    res["gecti"] = all(v["gecti"] for v in res.values())
    return res


def g0c() -> dict:
    """Üretici TurbParams == PyEco Faz 2 hücresinin TurbParams'ı (aynı hash)."""
    os.environ.setdefault("SLAM_DATA_ROOT", "/data/sdata")
    from turb.experiments import Cell, FAZ2_NOISE_DN, EUROC_FAZ2
    from turb.runner import _build_cfg
    from turb.hooks import params_from_cfg
    from turb.uret_sekans import arg_ayristirici, parametreleri_kur
    import config as _cfg
    cal, rekt = _cal_rect(SEQ)
    out, ok = [], True
    for dr in (0, 8, 16, 24):
        for sd in (0, 1, 2):
            cell = Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence="MH02", seed=sd, D_over_r0=float(dr),
                        windows=(600, 1000, 1400), window_len=400, psf_legacy_scale=False,
                        sensor_noise_dn=FAZ2_NOISE_DN, noise_seed=None)
            cfg = _build_cfg(cell, "/tmp/_g0c_yok")
            pp = params_from_cfg(cfg)
            a = arg_ayristirici().parse_args(["--girdi", str(SEQ), "--D_over_r0", str(dr), "--tohum", str(sd),
                                              "--gurultu_dn", "1"])
            pu, _ = parametreleri_kur(a, rekt.baseline)
            from dataclasses import asdict
            fark = {k: (v, asdict(pu)[k]) for k, v in asdict(pp).items() if asdict(pu)[k] != v}
            t = cfg["turbulence"]
            gur_ok = (float(t.get("sensor_noise_dn", 0)) == 1.0 and t.get("noise_seed") is None)
            e = pp.hash() == pu.hash() and gur_ok
            ok &= e
            out.append({"D_over_r0": dr, "tohum": sd, "hash_pyeco": pp.hash(), "hash_uretici": pu.hash(),
                        "fark": fark, "pyeco_gurultu": [t.get("sensor_noise_dn"), t.get("noise_seed")], "ayni": e})
    return {"hucreler": out, "config": str(_cfg.DEFAULT_PATH), "gecti": ok}


def g0d(kok: Path) -> dict:
    a = uret(kok / "g0d1", 8, 1.0)
    b = uret(kok / "g0d2", 8, 1.0)
    c = uret(kok / "g0d4", 8, 1.0, isci=4)
    oa, ob, oc = (ozet(x) for x in (a, b, c))
    img = [k for k in oa if k.endswith(".png")]
    r = {"n_png": len(img), "iki_kez_farkli": [k for k in img if oa[k] != ob.get(k)],
         "isci4_farkli": [k for k in img if oa[k] != oc.get(k)]}
    r["gecti"] = not r["iki_kez_farkli"] and not r["isci4_farkli"] and len(img) == 40
    return r


def g0e(kok: Path, desen="MH02__dr8_*n1_ns0__k600-620", alt="g0b") -> dict:
    """Üretilen rektifiye sensor.yaml'ı her rakibin kendi kuralıyla yeniden rektifiye et; sapma ve görüntü farkı."""
    d = next((kok / alt).glob(desen))
    y0 = yaml.safe_load(open(d / "mav0/cam0/sensor.yaml"))
    y1 = yaml.safe_load(open(d / "mav0/cam1/sensor.yaml"))
    W, H = y0["resolution"]
    K0 = np.array([[y0["intrinsics"][0], 0, y0["intrinsics"][2]], [0, y0["intrinsics"][1], y0["intrinsics"][3]], [0, 0, 1.]])
    K1 = np.array([[y1["intrinsics"][0], 0, y1["intrinsics"][2]], [0, y1["intrinsics"][1], y1["intrinsics"][3]], [0, 0, 1.]])
    D = np.zeros(5)
    Tb0 = np.array(y0["T_BS"]["data"]).reshape(4, 4)
    Tb1 = np.array(y1["T_BS"]["data"]).reshape(4, 4)
    Tcn = np.array(y1["T_cn_cnm1"]).reshape(4, 4)
    ts = sorted(int(l.split(",")[0]) for l in open(d / "mav0/cam0/data.csv") if l[0].isdigit())[0]
    im = [cv2.imread(str(d / f"mav0/cam{c}/data/{ts}.png"), cv2.IMREAD_GRAYSCALE) for c in (0, 1)]
    U, V = np.meshgrid(np.arange(W, dtype=np.float64), np.arange(H, dtype=np.float64))
    out = {}

    def olc(ad, KA, KB, DA, DB, R, T, alpha, f32=False, yeni_boyut=True):
        args = dict(flags=cv2.CALIB_ZERO_DISPARITY, alpha=alpha)
        if yeni_boyut:
            args["newImageSize"] = (W, H)
        R1, R2, P1, P2, *_ = cv2.stereoRectify(KA, DA, KB, DB, (W, H), R, T, **args)
        r = {}
        for c, (K, Dd, Rr, P) in enumerate(((KA, DA, R1, P1), (KB, DB, R2, P2))):
            mx, my = cv2.initUndistortRectifyMap(K, Dd, Rr, P[:3, :3], (W, H), cv2.CV_32FC1)
            sap = float(np.hypot(mx - U, my - V).max())
            rm = cv2.remap(im[c], mx, my, cv2.INTER_LINEAR)
            r[f"cam{c}"] = {"harita_sapma_maks_px": sap, "farkli_piksel": int(np.count_nonzero(rm != im[c])),
                            "farkli_oran": float(np.count_nonzero(rm != im[c]) / rm.size)}
        r["gecti"] = all(v["harita_sapma_maks_px"] < 0.01 for v in r.values() if isinstance(v, dict))
        out[ad] = r

    # ORB-SLAM3 (Settings::precomputeRectificationMaps): T_c1_c2 = inv(T_cn_cnm1); float32; alpha -1
    T01 = np.linalg.inv(Tcn)
    T21 = np.linalg.inv(T01.astype(np.float32).astype(float))
    f = lambda K: K.astype(np.float32).astype(float)
    olc("orbslam3", f(K0), f(K1), D[:4].reshape(1, -1), D[:4].reshape(1, -1), T21[:3, :3].copy(), T21[:3, 3].reshape(3, 1).copy(), -1)
    # DROID koşucusu: T_C0_C1 = inv(T_cn); R,T = inv(T_C0_C1) = T_cn; alpha 0
    olc("droid", K0, K1, D, D, Tcn[:3, :3].copy(), Tcn[:3, 3].copy(), 0, yeni_boyut=False)
    # AirSLAM (camera.cc): koşucu T_B_C1 = T_B_C0 inv(T_cn) yazar (%.12g), R10/t10 = inv(T_B_C1) T_B_C0; alpha 0
    r12 = lambda M: np.array([[float(f"{v:.12g}") for v in row] for row in M])
    Tb1a = r12(Tb0 @ np.linalg.inv(Tcn))
    T10 = np.linalg.inv(Tb1a) @ r12(Tb0)
    Ki = lambda y: np.array([[float(f"{y['intrinsics'][0]:.10g}"), 0, float(f"{y['intrinsics'][2]:.10g}")],
                             [0, float(f"{y['intrinsics'][1]:.10g}"), float(f"{y['intrinsics'][3]:.10g}")], [0, 0, 1.]])
    olc("airslam", Ki(y0), Ki(y1), D, D, T10[:3, :3].copy(), T10[:3, 3].copy(), 0)
    out["T_BS_fark_cam1_maks"] = float(np.abs(Tb1a - Tb1).max())
    out["gecti"] = all(v["gecti"] for k, v in out.items() if isinstance(v, dict))
    return out


def _pyeco_on_uc():
    """pyeco_turb ön-ucunun kalibrasyon/rektifikasyon kodu (modeller yüklenmeden; yalnız load_calibration)."""
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
    import config as _cfg
    from stereo_front_end import StereoMatcher as StereoFrontEnd
    cfg = _cfg.load()
    fe = object.__new__(StereoFrontEnd)
    fe.base_path = SEQ.parent
    fe.sequence = SEQ.name
    fe.STEREO_EXTRINSICS = str(cfg["frontend"].get("stereo_extrinsics", "gomulu"))
    fe.load_calibration()
    return fe, cfg["frontend"].get("stereo_extrinsics")


def _dy(r0s, r1s):
    sift = cv2.SIFT_create(2000)
    bf = cv2.BFMatcher()
    out = []
    for r0, r1 in zip(r0s, r1s):
        k0, d0 = sift.detectAndCompute(r0, None)
        k1, d1 = sift.detectAndCompute(r1, None)
        m = [a for a, b in bf.knnMatch(d0, d1, k=2) if a.distance < 0.7 * b.distance]
        dy = np.array([k0[a.queryIdx].pt[1] - k1[a.trainIdx].pt[1] for a in m])
        dx = np.array([k0[a.queryIdx].pt[0] - k1[a.trainIdx].pt[0] for a in m])
        out.append(dy[(dx > 0) & (dx < 120)])
    d = np.concatenate(out)
    return {"medyan_dy_px": float(np.median(d)), "medyan_abs_dy_px": float(np.median(np.abs(d))),
            "abs_dy_le_1px_orani": float(np.mean(np.abs(d) <= 1)), "n_eslesme": int(len(d))}


def g0f(kok: Path) -> dict:
    from turb.hooks import TurbulenceInjector
    from turb.uret_sekans import arg_ayristirici, parametreleri_kur, asl_oku
    fe, sx = _pyeco_on_uc()
    cal = asl_oku(SEQ)
    res = {"pyeco_stereo_extrinsics": sx, "pyeco_baseline_matrisi": fe.R.tolist() + [fe.T.tolist()]}
    temiz = uret(kok / "g0f", 0, 0.0, sy=True)
    res["ad_temiz"] = temiz.name
    n_fark, r0s, r1s, g0s, g1s = 0, [], [], [], []
    for i in range(600, 620):
        ts = cal["ts"][i]
        i0 = cv2.imread(str(SEQ / "mav0/cam0/data" / cal["f0"][ts]), cv2.IMREAD_GRAYSCALE)
        i1 = cv2.imread(str(SEQ / "mav0/cam1/data" / cal["f1"][ts]), cv2.IMREAD_GRAYSCALE)
        r0, r1 = fe.rectify_images(i0, i1)
        g = [cv2.imread(str(temiz / f"mav0/cam{c}/data/{ts}.png"), cv2.IMREAD_UNCHANGED) for c in (0, 1)]
        n_fark += int(np.count_nonzero(g[0] != r0)) + int(np.count_nonzero(g[1] != r1))
        r0s.append(r0); r1s.append(r1)
    res["sigma0_dr0_piksel_farki_pyeco_rectify"] = n_fark
    for dr in (0, 8):
        gur = uret(kok / "g0f", dr, 1.0, sy=True)
        a = arg_ayristirici().parse_args(["--girdi", str(SEQ), "--D_over_r0", str(dr), "--tohum", "0"])
        from turb.uret_sekans import Rektifikator
        p, _ = parametreleri_kur(a, Rektifikator(cal, "rektifiye", "sensor_yaml").baseline)
        inj = TurbulenceInjector(p, sensor_noise_dn=1.0, noise_seed=None)
        nf = 0
        for k, i in enumerate(range(600, 620)):
            ts = cal["ts"][i]
            k0, k1 = inj.apply_stereo(r0s[k], r1s[k], i)
            g = [cv2.imread(str(gur / f"mav0/cam{c}/data/{ts}.png"), cv2.IMREAD_UNCHANGED) for c in (0, 1)]
            nf += int(np.count_nonzero(g[0] != k0)) + int(np.count_nonzero(g[1] != k1))
        res[f"sigma1_dr{dr}_piksel_farki_pyeco_kanca"] = nf
        res[f"ad_dr{dr}"] = gur.name
    # dikey kayma: sensor_yaml (PyEco) ↔ gömülü (eski üretici varsayılanı), aynı kareler
    eski = next((kok / "g0a").glob("MH02__dr0_*__k600-620"))
    e0 = [cv2.imread(str(eski / f"mav0/cam0/data/{cal['ts'][i]}.png"), 0) for i in (600, 610, 619)]
    e1 = [cv2.imread(str(eski / f"mav0/cam1/data/{cal['ts'][i]}.png"), 0) for i in (600, 610, 619)]
    res["dikey_kayma_sensor_yaml"] = _dy([r0s[0], r0s[10], r0s[19]], [r1s[0], r1s[10], r1s[19]])
    res["dikey_kayma_gomulu"] = _dy(e0, e1)
    res["gecti"] = (sx == "sensor_yaml" and n_fark == 0 and res["sigma1_dr0_piksel_farki_pyeco_kanca"] == 0
                    and res["sigma1_dr8_piksel_farki_pyeco_kanca"] == 0
                    and res["dikey_kayma_sensor_yaml"]["medyan_abs_dy_px"] < 0.5)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--once", required=True)
    ap.add_argument("--is", dest="isdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--yalniz", default="a,b,c,d,e")
    a = ap.parse_args()
    kok = Path(a.isdir)
    kok.mkdir(parents=True, exist_ok=True)
    R = json.loads(Path(a.out).read_text()) if Path(a.out).exists() else {}
    for k in a.yalniz.split(","):
        f = {"a": lambda: g0a(Path(a.once), kok), "b": lambda: g0b(kok), "c": g0c,
             "d": lambda: g0d(kok), "e": lambda: g0e(kok), "f": lambda: g0f(kok),
             "e_sy": lambda: g0e(kok, "MH02__dr8_*_sy_n1_ns0__k600-620", "g0f")}[k]
        R[f"G0{k}"] = f()
        print(f"G0{k}:", "GEÇTİ" if R[f"G0{k}"]["gecti"] else "KALDI", flush=True)
        Path(a.out).write_text(json.dumps(R, indent=1, ensure_ascii=False, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


if __name__ == "__main__":
    main()
