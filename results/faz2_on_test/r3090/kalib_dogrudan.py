"""V-serisi kol dizinlerinden (sürücü kayıtları sıra beklerken) doğrudan kalibrasyon: runner ile aynı metrik fonksiyonları."""
import sys, os, json, glob, numpy as np
sys.path.insert(0, "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb")
from turb.metrics import nees_from_dump, trajectory_metrics
from turb.pencere_sec import kare_zamanlari
F = "/tmp/faz2_on_test/r3090"
BASE = "/media/masa24/2tb ana/sdata"
W = json.load(open("/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb/evidence/pencereler_faz2_kapsam1m.json"))["sequences"]
rows = {}
for s in ["V1_01_easy", "V1_02_medium", "V2_01_easy", "V2_02_medium"]:
    ts = kare_zamanlari(BASE, s)
    gt = f"{BASE}/{s}/mav0/state_groundtruth_estimate0/data.csv"
    for d in sorted(glob.glob(f"{F}/fe/{s}_dr0_s0_n1/kol_*")):
        tag = os.path.basename(d)[4:]
        cov = glob.glob(f"{d}/cov_*.npz"); rt = glob.glob(f"{d}/traj_*_realtime.csv")
        if not cov or not rt:
            continue
        for w in W[s]["starts"]:
            t0, t1 = float(ts[w]), float(ts[min(w + 400, len(ts)) - 1])
            n = nees_from_dump(cov[0], gt, t_min=t0, t_max=t1, rel_key="rel")
            m = trajectory_metrics(rt[0], gt, t_min=t0, t_max=t1)
            rows[(s, w, tag)] = (n["nees_pos_median"], n["n"], m["ate_se3"], n["nees_rot_median"])
tags = sorted({k[2] for k in rows})
print("etiketler:", tags)
a0 = {k[:2]: v for k, v in rows.items() if k[2] == "A0"}
print("A0 n/pencere:", {f"{k[0][:5]}/{k[1]}": v[1] for k, v in sorted(a0.items())})
print("A0 NEES_pos:", {f"{k[0][:5]}/{k[1]}": round(v[0], 1) for k, v in sorted(a0.items())})
med_a0 = np.median([v[0] for v in a0.values()])
print(f"A0 pencere-medyanlarının medyanı {med_a0:.2f} (n pencere {len(a0)})")
print("| k | pencere | gerçek A0C medyan NEES_pos | A0/k² | NEES_rot | ATE A0C/A0 (medyan) |")
best = None
for k in (1, 2, 4, 8, 12, 16, 24):
    t = {kk[:2]: v for kk, v in rows.items() if kk[2] == f"A0C_k{k}"}
    if not t:
        continue
    med = np.median([v[0] for v in t.values()]); rot = np.median([v[3] for v in t.values()])
    ate = np.median([t[c][2] / a0[c][2] for c in t if c in a0])
    print(f"| {k} | {len(t)} | {med:.2f} | {med_a0/k**2:.2f} | {rot:.2f} | {ate:.3f} |")
    if len(t) == 12 and (best is None or abs(med - 2.37) < abs(best[1] - 2.37)):
        best = (k, med)
print("SEÇİLEN k:", best)
c10 = {kk[:2]: v for kk, v in rows.items() if kk[2] == "A0_cov10"}
if c10:
    print("cov10 n/pencere:", {f"{k[0][:5]}/{k[1]}": v[1] for k, v in sorted(c10.items())},
          "NEES_pos cov10/cov2 medyan oranı:", round(float(np.median([c10[c][0] / a0[c][0] for c in c10 if c in a0])), 3),
          "ATE aynı:", all(abs(c10[c][2] - a0[c][2]) < 1e-12 for c in c10 if c in a0))
