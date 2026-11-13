"""
AirSim'de arka-uç kolları + mutlak NEES (GT kusursuz): (sekans, C0, tohum) × kollar.

    python kol_izgara.py --sekanslar AF_F5_h5 --C0 0 1e-11 --kollar A0 AG A3 --cikti cikti_izgara/kollar.jsonl

- Türbülanslı sekans bir kez üretilir (saklanır, sonda silinir).
- A0 ön-ucu koşar; diğer kollar aynı ön-uç dökümünü kullanır (DUMP_KAYNAK). Kollar için
  turbulence.enabled=true, D_over_r0=0 (ön-uç kancası kimlik), arm=<kol>.
- NEES: pencere-göreli kovaryans dökümü (COV_OUT), turb.metrics.nees_from_dump; bütün sekans.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from turb_izgara import PE, PY, VERI, TURB_KOK, D_FOV, metrik  # noqa: E402


def nees(cov, gt):
    code = ("import sys,json; sys.path.insert(0,'.'); from turb.metrics import nees_from_dump; "
            f"print(json.dumps(nees_from_dump({str(cov)!r}, {str(gt)!r})))")
    r = subprocess.run([PY, "-c", code], cwd=PE, capture_output=True, text=True)
    try:
        return json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return {"hata": r.stderr[-400:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sekanslar", nargs="+", required=True)
    ap.add_argument("--C0", nargs="+", type=float, required=True)
    ap.add_argument("--tohumlar", nargs="+", type=int, default=[0])
    ap.add_argument("--kollar", nargs="+", default=["A0", "AG", "A3"])
    ap.add_argument("--cov_every", type=int, default=3)
    ap.add_argument("--dtm", default=str(HERE / "africa/harita_africa.npz"))
    ap.add_argument("--cikti", required=True)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--isci", type=int, default=14)
    a = ap.parse_args()
    out = Path(a.cikti); out.parent.mkdir(parents=True, exist_ok=True)
    yapildi = set()
    if out.exists():
        for l in open(out):
            r = json.loads(l); yapildi.add((r["sekans"], r["C0"], r["tohum"], r["kol"]))
    for seq in a.sekanslar:
        fov = int(round(json.load(open(VERI / seq / "kayit_meta.json"))["arguman"]["fov"]))
        D = D_FOV.get(fov, 4.8e-3)
        for C0 in a.C0:
            for th in a.tohumlar:
                if C0 == 0 and th != a.tohumlar[0]:
                    continue
                if all((seq, C0, th, k) in yapildi for k in a.kollar):
                    continue
                if C0 == 0:
                    root, name = VERI, seq
                else:
                    r = subprocess.run([PY, "-m", "turb.uret_derinlik", "--girdi", str(VERI / seq), "--cikti_kok",
                                        str(TURB_KOK), "--C0", str(C0), "--tohum", str(th), "--dtm", a.dtm,
                                        "--D", str(D), "--isci", str(a.isci)], cwd=PE, capture_output=True, text=True)
                    if r.returncode != 0:
                        sys.exit(r.stderr[-1500:])
                    root, name = TURB_KOK, json.loads(r.stdout.strip().splitlines()[-1])["cikti"].split("/")[-1]
                gt = root / name / "mav0/state_groundtruth_estimate0/data.csv"
                dump = HERE / "cikti_kollar" / f"{name}__A0" / f"landmarks_{name}.pkl"
                for kol in a.kollar:
                    if (seq, C0, th, kol) in yapildi:
                        continue
                    od = HERE / "cikti_kollar" / f"{name}__{kol}"
                    od.mkdir(parents=True, exist_ok=True)
                    cov = od / "cov.npz"
                    env = dict(os.environ, GPU=a.gpu, COV_OUT=str(cov), COV_EVERY=str(a.cov_every))
                    ek = []
                    if kol != "A0":
                        ek = ["turbulence.enabled=true", "turbulence.D_over_r0=0.0", f"turbulence.arm={kol}"]
                        env["DUMP_KAYNAK"] = str(dump)
                    t0 = time.time()
                    subprocess.run([str(HERE / "pyeco_kos.sh"), str(root), name, str(od)] + ek, env=env,
                                   capture_output=True, text=True)
                    rec = {"sekans": seq, "C0": C0, "tohum": th, "kol": kol, "fov": fov, "sure_s": time.time() - t0}
                    traj = od / f"gtsam_traj_{name}.csv"
                    if traj.exists():
                        rec.update(metrik(traj, gt)); rec["durum"] = "tamam"
                        rec["nees"] = nees(cov, gt) if cov.exists() else {"hata": "kovaryans yok"}
                    else:
                        rec["durum"] = "hata"
                    open(out, "a").write(json.dumps(rec) + "\n")
                    print(json.dumps({k: rec.get(k) for k in ("sekans", "C0", "kol", "durum", "ate_se3", "sure_s")}
                                     | {"nees_poz_med": rec.get("nees", {}).get("nees_pos_median"),
                                        "nees_don_med": rec.get("nees", {}).get("nees_rot_median")}), flush=True)
                if C0 != 0:
                    shutil.rmtree(root / name, ignore_errors=True)


if __name__ == "__main__":
    main()
