"""Çiftli analiz: aynı GPU ailesi içinde (sekans, pencere, tohum) başına C0>0 / C0=0 oranları (ATE, RPE1s, RPE10s).
Kaynak jsonl'ler GPU'ya göre verilir:  python analiz_cift.py gpu1:cikti_izgara/f5_cift.jsonl,cikti_izgara/gece2_gpu1.jsonl gpu0:cikti_izgara/gece2_gpu0.jsonl"""
import json, sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from rpe_ozet import ozet

HERE = Path(__file__).parent
satir = []
for arg in sys.argv[1:]:
    gpu, files = arg.split(":")
    for f in files.split(","):
        for l in open(HERE / f):
            r = json.loads(l)
            if r.get("durum") != "tamam" or not r.get("gurultu"):
                continue
            od = HERE / r["cikti_dizin"] if r.get("cikti_dizin") else HERE / "cikti_izgara" / r["turb_sekans"]
            traj = od / f"gtsam_traj_{r['turb_sekans']}.csv"
            pen = r["pencere"].replace(":", "-") if r.get("pencere") else None
            o = ozet(traj, r["sekans"], pen)
            # yörünge dosyası bu kayda mı ait? (aynı turb_sekans adı başka GPU koşusunca ezilmiş olabilir)
            if abs(o["ate"] - r["ate_se3"]) > 1e-6 * max(1.0, r["ate_se3"]):
                print(f"UYARI ezilmiş yörünge: {gpu} {r['turb_sekans']} jsonl ATE {r['ate_se3']:.4f} dosya {o['ate']:.4f} -> RPE yok")
                o = dict(o, ate=r["ate_se3"], rpe1s=float("nan"), rpe10s=float("nan"))
            satir.append(dict(gpu=gpu, sekans=r["sekans"], pencere=r.get("pencere"), C0=r["C0"], tohum=r["tohum"], **o))
fam = {}
for s in satir:
    fam.setdefault((s["gpu"], s["sekans"], s["pencere"]), {}).setdefault(s["C0"], {})[s["tohum"]] = s
out = []
for (g, sq, pen), d in sorted(fam.items(), key=lambda x: str(x[0])):
    if 0.0 not in d:
        continue
    for C0 in sorted(c for c in d if c > 0):
        for th, t in sorted(d[C0].items()):
            c = d[0.0].get(th)
            if c:
                out.append(dict(gpu=g, sekans=sq, pencere=pen, C0=C0, tohum=th,
                                ate_c=c["ate"], ate_t=t["ate"], r_ate=t["ate"] / c["ate"],
                                rpe_c=c["rpe1s"], rpe_t=t["rpe1s"], r_rpe=t["rpe1s"] / c["rpe1s"]))
print(f"{'gpu':5}{'sekans':11}{'pencere':11}{'C0':>7}{'th':>3}  ATE c/t     oran   RPE1s c/t     oran")
for o in out:
    print(f"{o['gpu']:5}{o['sekans']:11}{str(o['pencere']):11}{o['C0']:7.0e}{o['tohum']:3d}  {o['ate_c']:.3f}/{o['ate_t']:.3f} {o['r_ate']:5.2f}   {o['rpe_c']:.4f}/{o['rpe_t']:.4f} {o['r_rpe']:5.2f}")
print("\nözet (sekans, C0): medyan oran ATE / RPE1s, n, temiz RPE1s CV")
grp = {}
for o in out:
    grp.setdefault((o["sekans"], o["C0"]), []).append(o)
for k, v in sorted(grp.items()):
    rr = [x["r_rpe"] for x in v if np.isfinite(x["r_rpe"])]
    print(k, "ATE", round(float(np.median([x["r_ate"] for x in v])), 2), "n", len(v),
          "| RPE1s", round(float(np.median(rr)), 2) if rr else None, "n", len(rr), "RPE1s>1.3:", sum(x > 1.3 for x in rr))
for (g, sq, pen), d in sorted(fam.items(), key=lambda x: str(x[0])):
    if 0.0 in d and len(d[0.0]) > 1:
        r = np.array([x["rpe1s"] for x in d[0.0].values()])
        if not np.all(np.isfinite(r)): r = np.array([np.nan])
        a = np.array([x["ate"] for x in d[0.0].values()])
        print("temiz", g, sq, pen, "RPE1s CV %.0f%%" % (100 * r.std(ddof=1) / r.mean()), "ATE CV %.0f%%" % (100 * a.std(ddof=1) / a.mean()))
json.dump(out, open(HERE / "cikti_izgara/analiz_cift.json", "w"), indent=1)
