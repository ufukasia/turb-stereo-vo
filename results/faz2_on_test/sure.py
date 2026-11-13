"""Tam ızgara süre/disk tahmini: ölçülen hücre sürelerinden."""
import json, sys, numpy as np
from collections import defaultdict
F = "/tmp/faz2_on_test"
R = [json.loads(l) for l in open(f"{F}/duman.jsonl")] + [json.loads(l) for l in open(f"{F}/kalib.jsonl")]
# kol-çalışması başına bir değer (pencere kayıtları aynı t_backend'i taşır)
be = defaultdict(list)       # (kol, dr) -> [s / kare]
seen = set()
for r in R:
    if r.get("status") != "ok":
        continue
    tag = r["kol_etiket"]
    if "_cov" in tag:
        continue
    c = r["cell"]
    k = (c["sequence"], c["D_over_r0"], tag)
    if k in seen:
        continue
    seen.add(k)
    arm = r["arm"]
    nfr = r["n_dump_frames"]
    be[(arm, c["D_over_r0"])].append(r["t_backend"] / nfr)
fe = {}
for r in R:
    c = r["cell"]
    fe[(c["sequence"], c["D_over_r0"])] = r["t_frontend"] / r["n_dump_frames"]
print("arka-uç s/kare (6 paralel işçi, OMP=4, makinede başka iş de var):")
for k in sorted(be):
    print(f"  {k}: medyan {np.median(be[k]):.3f}  maks {max(be[k]):.3f}  n={len(be[k])}")
print("ön-uç s/kare:", {k: round(v, 3) for k, v in sorted(fe.items())})

W = json.load(open("/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb/evidence/pencereler_faz2.json"))["sequences"]
DR = (0, 4, 8, 12, 16, 20, 24)
ARMS = ("A0", "A0C", "A1", "A3", "AG", "A0C+A3")
A5S = ("MH02", "MH04")


def fe_rate(dr):
    # ölçülen: 0, 8, 16 -> doğrusal ara/dış değerleme
    xs = [0, 8, 16]
    ys = [np.median([v for (s, d), v in fe.items() if d == x]) for x in xs]
    return float(np.interp(dr, xs, ys)) if dr <= 16 else ys[2] + (ys[2] - ys[1]) / 8 * (dr - 16)


def be_rate(arm, dr):
    xs = sorted({d for (a, d) in be if a == arm})
    ys = [np.median(be[(arm, x)]) for x in xs]
    if dr <= xs[-1]:
        return float(np.interp(dr, xs, ys))
    sl = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]) if len(xs) > 1 else 0
    return ys[-1] + max(sl, 0) * (dr - xs[-1])


seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 5
fe_tot = be_tot = 0.0
ncell = nrun = 0
for s, v in W.items():
    nf = v["windows"][-1]["end_excl"]
    for dr in DR:
        for sd in (range(seeds) if dr > 0 else [0]):
            ncell += 1
            fe_tot += fe_rate(dr) * nf
            for a in ARMS + (("A5",) if s in A5S else ()):
                nrun += 1
                be_tot += be_rate(a, dr) * nf
print(f"hücre {ncell}, kol koşumu {nrun}")
print(f"ön-uç toplam (tek akış, ölçülen yük altında): {fe_tot/3600:.1f} sa")
print(f"arka-uç toplam (tek işçi-süresi): {be_tot/3600:.1f} sa  -> 6 işçiyle {be_tot/6/3600:.1f} sa, "
      f"10 işçiyle (ölçeklenirse) {be_tot/10/3600:.1f} sa")
for g, nfe in (("GPU1 (4060 Ti) tek", 1), ("4060 Ti, 2 eşzamanlı ön-uç", 2), ("iki GPU, 2+2 ön-uç", 4)):
    for w in (6, 10):
        t = max(fe_tot / nfe, be_tot / w)
        print(f"  {g}, {w} arka-uç işçisi: ~{t/3600:.1f} sa (darboğaz {'ön-uç' if fe_tot/nfe > be_tot/w else 'arka-uç'})")
