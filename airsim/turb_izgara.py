"""
AirSim türbülans ızgarası: (sekans, C0, tohum) -> türbülanslı sekans üret -> PyEco (A0) koş -> metrik.

    python turb_izgara.py --sekanslar AF_F15_h5 AF_F5_h5 --C0 0 1e-12 1e-11 --tohumlar 0 \
        --cikti cikti_izgara/pilot.jsonl

Sonuç satırı (jsonl): sekans, C0, tohum, ATE_SE3, Sim3 ölçeği, Sim3 ATE, yol, son kare hatası,
türbülans istatistikleri (derinlik dilimli medyanlar, sekans ortalaması), süreler, kod özetleri.
Kaldığı yerden devam eder (aynı anahtar varsa atlar). Türbülanslı sekans, PyEco bitince silinir
(--sakla ile tutulur), disk yükü sınırlı kalır.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PE = HERE.parent / "pyeco_turb"
PY = str(HERE.parent / "kod/turbvio/bin/python")
VERI = Path("/data/sdata_airsim")
TURB_KOK = Path("/data/sdata_turb_airsim")
D_FOV = {15: 4.8e-3, 5: 14.4e-3}


def _umeyama(src, dst, with_scale):
    mu_s, mu_d = src.mean(0), dst.mean(0)
    xs, xd = src - mu_s, dst - mu_d
    C = xd.T @ xs / len(src)
    U, D, Vt = np.linalg.svd(C)
    S = np.eye(3)
    if np.linalg.det(U) * np.linalg.det(Vt) < 0:
        S[2, 2] = -1
    R = U @ S @ Vt
    s = np.trace(np.diag(D) @ S) / xs.var(0).sum() if with_scale else 1.0
    return s, R, mu_d - s * R @ mu_s


def metrik(traj_csv, gt_csv):
    gt = np.loadtxt(gt_csv, delimiter=",", comments="#")
    e = np.loadtxt(traj_csv, delimiter=",", comments="#")
    ex = np.stack([np.interp(gt[:, 0], e[:, 0], e[:, k]) for k in (1, 2, 3)], 1)
    out = {}
    for ad, ws in (("se3", False), ("sim3", True)):
        s, R, t = _umeyama(ex, gt[:, 1:4], ws)
        al = (s * (R @ ex.T)).T + t
        err = np.linalg.norm(al - gt[:, 1:4], axis=1)
        out[f"ate_{ad}"] = float(np.sqrt(np.mean(err ** 2)))
        if ws:
            out["sim3_olcek"] = float(s)
        else:
            out["son_hata"] = float(err[-1])
    out["yol_m"] = float(np.linalg.norm(np.diff(gt[:, 1:4], axis=0), axis=1).sum())
    out["n_poz"] = int(len(e))
    return out


def turb_ozet(meta_json):
    if not Path(meta_json).exists():
        return {}
    m = json.load(open(meta_json))
    ki = m.get("kare_istat", {})
    agg = {}
    for st in ki.values():
        for dil, v in st.items():
            a = agg.setdefault(dil, {k: [] for k in v})
            for k, x in v.items():
                a[k].append(x)
    return {"kod": m.get("kod_ozeti"), "uretim_s": m.get("sure_s"),
            "dilim": {d: {k: float(np.median(x)) for k, x in a.items()} for d, a in agg.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sekanslar", nargs="+", required=True)
    ap.add_argument("--C0", nargs="+", type=float, required=True)
    ap.add_argument("--tohumlar", nargs="+", type=int, default=[0])
    ap.add_argument("--dtm", default=str(HERE / "africa/harita_africa.npz"))
    ap.add_argument("--cikti", required=True)
    ap.add_argument("--isci", type=int, default=16)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--sakla", action="store_true")
    ap.add_argument("--pencereler", nargs="+", default=[None], help="a:b kare pencereleri (yoksa bütün sekans)")
    ap.add_argument("--gurultu", type=float, default=0.0, help="sensör gürültüsü sigma DN (tüm kollara)")
    ap.add_argument("--sabit_kayma", type=float, default=0.0, help="KONTROL kolu (C0=0): sabit alt-piksel kayma px")
    a = ap.parse_args()
    out = Path(a.cikti); out.parent.mkdir(parents=True, exist_ok=True)
    yapildi = set()
    if out.exists():
        for l in open(out):
            r = json.loads(l); yapildi.add((r["sekans"], r["C0"], r["tohum"], r.get("pencere"), r.get("gurultu", 0.0), r.get("sabit_kayma", 0.0)))
    # GPU başına ayrı üretim ve çıktı dizini: aynı turb_sekans adının başka GPU koşusunca ezilmesini önler
    turb_kok = TURB_KOK / f"gpu{a.gpu}"
    turb_kok.mkdir(parents=True, exist_ok=True)
    cikti_kok = HERE / "cikti_izgara" / f"gpu{a.gpu}"
    for seq in a.sekanslar:
        meta = json.load(open(VERI / seq / "kayit_meta.json"))
        fov = int(round(meta["arguman"]["fov"]))
        D = D_FOV.get(fov, 4.8e-3)
        for pen in a.pencereler:
          for C0 in a.C0:
            for th in a.tohumlar:
                if C0 == 0 and th != a.tohumlar[0] and a.gurultu == 0:
                    continue                                   # gürültüsüz temiz kol tohumdan bağımsız
                if (seq, C0, th, pen, a.gurultu, a.sabit_kayma) in yapildi:
                    continue
                t0 = time.time()
                if C0 == 0 and pen is None and a.gurultu == 0 and a.sabit_kayma == 0:
                    tseq_root, tseq_name = VERI, seq
                else:
                    ek = (["--kareler", pen] if pen else []) + (["--gurultu", str(a.gurultu)] if a.gurultu else []) \
                        + (["--sabit_kayma", str(a.sabit_kayma)] if a.sabit_kayma else [])
                    r = subprocess.run([PY, "-m", "turb.uret_derinlik", "--girdi", str(VERI / seq),
                                        "--cikti_kok", str(turb_kok), "--C0", str(C0), "--tohum", str(th),
                                        "--dtm", a.dtm, "--D", str(D), "--isci", str(a.isci)] + ek,
                                       cwd=PE, capture_output=True, text=True)
                    if r.returncode != 0:
                        print(r.stderr[-2000:]); sys.exit(f"üretim hatası {seq} {C0} {th}")
                    tseq_root = turb_kok; tseq_name = json.loads(r.stdout.strip().splitlines()[-1])["cikti"].split("/")[-1]
                t1 = time.time()
                od = cikti_kok / f"{tseq_name}"
                env = dict(os.environ, GPU=a.gpu)
                r = subprocess.run([str(HERE / "pyeco_kos.sh"), str(tseq_root), tseq_name, str(od)],
                                   env=env, capture_output=True, text=True)
                t2 = time.time()
                traj = od / f"gtsam_traj_{tseq_name}.csv"
                rec = {"sekans": seq, "C0": C0, "tohum": th, "pencere": pen, "gurultu": a.gurultu, "sabit_kayma": a.sabit_kayma, "fov": fov, "D": D, "turb_sekans": tseq_name, "cikti_dizin": str(od.relative_to(HERE)), "gpu": a.gpu,
                       "uretim_s": t1 - t0, "pyeco_s": t2 - t1}
                if traj.exists():
                    rec.update(metrik(traj, tseq_root / tseq_name / "mav0/state_groundtruth_estimate0/data.csv"))
                    rec["durum"] = "tamam"
                else:
                    rec["durum"] = "pyeco_hata"; rec["log"] = (r.stdout + r.stderr)[-500:]
                rec["turb"] = turb_ozet(tseq_root / tseq_name / "turb_meta.json")
                open(out, "a").write(json.dumps(rec) + "\n")
                print(json.dumps({k: rec.get(k) for k in ("sekans", "C0", "tohum", "durum", "ate_se3", "sim3_olcek",
                                                          "uretim_s", "pyeco_s")}), flush=True)
                if tseq_root == turb_kok and not a.sakla:
                    shutil.rmtree(tseq_root / tseq_name, ignore_errors=True)


if __name__ == "__main__":
    main()
