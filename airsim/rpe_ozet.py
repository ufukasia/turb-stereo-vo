"""Mevcut tüm AirSim yörüngeleri için ATE + RPE (Δ=1 s ve 10 s; öteleme RMSE, m) özeti.
GT kaynak sekanstan alınır (türbülanslı sekanslar silinmiş olabilir; GT aynıdır, pencere kırpılır)."""
import glob, json, re, sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from turb_izgara import _umeyama, VERI

def rpe(est, gt, d):
    # est, gt: (N,3) aynı zaman damgalarında; SE3 hizalı. Öteleme farkının farkı (yön-bağımsız kaba RPE)
    e = (est[d:] - est[:-d]) - (gt[d:] - gt[:-d])
    return float(np.sqrt(np.mean(np.sum(e ** 2, 1))))

def ozet(traj, seq, pen):
    g = np.loadtxt(VERI / seq / "mav0/state_groundtruth_estimate0/data.csv", delimiter=",", comments="#")
    if pen:
        a, b = map(int, pen.split("-")); g = g[a:b]
    e = np.loadtxt(traj, delimiter=",", comments="#")
    ex = np.stack([np.interp(g[:, 0], e[:, 0], e[:, k]) for k in (1, 2, 3)], 1)
    s, R, t = _umeyama(ex, g[:, 1:4], False); al = (R @ ex.T).T + t
    return {"ate": float(np.sqrt(np.mean(np.sum((al - g[:, 1:4]) ** 2, 1)))),
            "rpe1s": rpe(al, g[:, 1:4], 20), "rpe10s": rpe(al, g[:, 1:4], 200)}

def main():
    rows = []
    for traj in sorted(glob.glob("cikti_izgara/*/gtsam_traj_*.csv") + glob.glob("cikti_izgara/gpu*/*/gtsam_traj_*.csv") + glob.glob("cikti_kontrol/*/gtsam_traj_*.csv")):
        d = Path(traj).parent.name
        name = Path(traj).stem.replace("gtsam_traj_", ""); seq = name.split("__")[0]
        m = re.search(r"_k(\d+-\d+)", d); pen = m.group(1) if m else None
        if not (VERI / seq).exists():
            continue
        o = ozet(traj, seq, pen); o["kosum"] = Path(traj).parent.parent.name + "/" + d
        cfg = Path(traj).parent / "fe.log"
        rows.append(o)
    for o in rows:
        print(f"{o['kosum']:<55} ATE {o['ate']:.3f}  RPE1s {o['rpe1s']:.4f}  RPE10s {o['rpe10s']:.4f}")
    json.dump(rows, open("cikti_izgara/rpe_ozet.json", "w"), indent=1)


if __name__ == "__main__":
    main()
