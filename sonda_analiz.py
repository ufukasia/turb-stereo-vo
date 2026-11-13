#!/usr/bin/env python
"""Keşif sondası analizi (ON_BEYAN.md §3-§4). Birim: sekans; ATE pencerelerde aritmetik, ANEES geometrik ortalama
(turb/faz2_analiz.py ile aynı kural); etki = birimlerde medyan oran.

    python sonda_analiz.py [duman.jsonl | records.jsonl]
"""
import collections
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(ROOT, "evidence", "SONDA_KESIF")
dosya = os.path.join(EV, sys.argv[1] if len(sys.argv) > 1 else "records.jsonl")


def oku(p):
    out = {}
    with open(p) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "arm" not in r:
                print("  !! hücre hatası:", r.get("key"), r.get("error"))
                continue
            c = r["cell"]
            rejim = "R0" if abs(c["tau0_over_dt"] - 1.0) < 1e-9 else "R5"
            tag = c.get("tag", "")
            fe = "fe1" if "_fe1" in tag else "fe0"
            kol = "A0g" if tag.endswith("_bg") else r["arm"]
            out[(rejim, float(c["D_over_r0"]), fe, kol, c["sequence"], int(r["window"]))] = r      # son kayıt geçerli
    return out


S = oku(dosya)
print(f"[analiz] {dosya}: {len(S)} kayıt")

# --- ana ızgara (E1) ---
ANA = {}
with open(os.path.join(ROOT, "evidence", "FAZ2", "records.jsonl")) as f:
    for line in f:
        r = json.loads(line)
        c = r.get("cell") or {}
        if c.get("seed") == 0 and r.get("arm") in ("A0", "A1"):
            ANA[(float(c["D_over_r0"]), r["arm"], c["sequence"], int(r["window"]))] = r


def basarili(r):
    return r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged")


def ate(r):
    v = r.get("m_ate_se3")
    return float(v) if (r.get("status") == "ok" and r.get("m_valid") and v is not None and math.isfinite(v) and v > 0) else None


def anees(r):
    v = r.get("n_nees_pos_mean")
    return float(v) if (r.get("status") == "ok" and r.get("n_valid") and v is not None and math.isfinite(v) and v > 0) else None


print("\n== E1: R0/fe0 taban hücreleri ana ızgarayla (m_ate_se3) ==")
n_ayni = n_fark = 0
farklar = []
for (rejim, dr, fe, kol, s, w), r in sorted(S.items()):
    if rejim == "R0" and fe == "fe0" and kol in ("A0", "A1"):
        m = ANA.get((dr, kol, s, w))
        if m is None:
            continue
        a, b = r.get("m_ate_se3"), m.get("m_ate_se3")
        if a == b:
            n_ayni += 1
        else:
            n_fark += 1
            farklar.append((s, dr, w, kol, a, b))
print(f"  bit-aynı {n_ayni}, farklı {n_fark}")
for x in farklar[:12]:
    print("   fark:", x)

print("\n== ön-uç istatistikleri (taban hücreleri, sekans ortalaması) ==")
print("  rejim D/r0 fe   n  stereo_es  ucgen  verim  zamansal  aktif_derinlik  stereo_oran")
FE = collections.defaultdict(dict)
for (rejim, dr, fe, kol, s, w), r in S.items():
    if kol == "A0" and r.get("frontend_stats"):
        FE[(rejim, dr, fe)][s] = r["frontend_stats"]
for k in sorted(FE):
    v = list(FE[k].values())
    g = lambda n: float(np.mean([x.get(n, np.nan) for x in v]))      # noqa: E731
    print(f"  {k[0]}  {k[1]:4.0f}  {k[2]}  {len(v):2d}  {g('avg_stereo_matches'):8.1f}  "
          f"{g('avg_current_stereo_triangulated_points'):6.1f}  "
          f"{g('avg_current_stereo_triangulated_points') / g('avg_stereo_matches'):5.3f}  "
          f"{g('avg_temporal_matches'):8.1f}  {g('avg_active_depth_tracks'):10.1f}  {g('stereo_run_ratio'):8.2f}")


def birim(rejim, dr, fe, kol, f, tur):
    """sekans -> pencere ortalaması (yalnız geçerli pencereler)."""
    g = collections.defaultdict(list)
    for (rj, d, fe_, k, s, w), r in S.items():
        if (rj, d, fe_, k) == (rejim, dr, fe, kol):
            v = f(r)
            if v is not None:
                g[s].append((w, v))
    return g


def oran(pay, payda, f, tur):
    """Eşli (sekans, pencere) kesişiminde birim ortalaması; birimlerde medyan oran."""
    A, B = birim(*pay, f, tur), birim(*payda, f, tur)
    lr = []
    for s in sorted(set(A) & set(B)):
        a, b = dict(A[s]), dict(B[s])
        ws = sorted(set(a) & set(b))
        if not ws:
            continue
        xa, xb = np.array([a[w] for w in ws]), np.array([b[w] for w in ws])
        if tur == "geo":
            lr.append(float(np.mean(np.log(xa)) - np.mean(np.log(xb))))
        else:
            lr.append(float(np.log(np.mean(xa) / np.mean(xb))))
    if not lr:
        return "      -      "
    lr = np.array(lr)
    return f"{math.exp(np.median(lr)):6.3f} ({int(np.sum(lr < 0))}/{len(lr)})"


def tablo(baslik, f, tur):
    print(f"\n== {baslik}: medyan eşli oran (oran < 1 olan sekans / sekans) ==")
    print("  rejim D/r0 |  A1/A0 [fe0]    A0g/A0 [fe0]   A0g/A1 [fe0]   fe1/fe0 [A0]   fe1/fe0 [A0g]  (A0g,fe1)/(A0,fe0)  (A0g,fe1)/(A1,fe0)")
    for rejim, dr in sorted(set((k[0], k[1]) for k in S)):
        o = lambda p, q: oran((rejim, dr) + p, (rejim, dr) + q, f, tur)      # noqa: E731
        print(f"  {rejim}  {dr:4.0f}  | {o(('fe0', 'A1'), ('fe0', 'A0'))}  {o(('fe0', 'A0g'), ('fe0', 'A0'))}  "
              f"{o(('fe0', 'A0g'), ('fe0', 'A1'))}  {o(('fe1', 'A0'), ('fe0', 'A0'))}  {o(('fe1', 'A0g'), ('fe0', 'A0g'))}  "
              f"   {o(('fe1', 'A0g'), ('fe0', 'A0'))}       {o(('fe1', 'A0g'), ('fe0', 'A1'))}")


tablo("ATE (m_ate_se3)", ate, "ari")
tablo("ANEES (n_nees_pos_mean)", anees, "geo")

print("\n== düzeyler: medyan ATE [m] / medyan ANEES / başarısız pencere ==")
print("  rejim D/r0 fe  kol    n   ATE_med   ANEES_med  başarısız")
G = collections.defaultdict(list)
for (rejim, dr, fe, kol, s, w), r in S.items():
    G[(rejim, dr, fe, kol)].append(r)
for k in sorted(G):
    v = G[k]
    a = [ate(r) for r in v if ate(r) is not None]
    n = [anees(r) for r in v if anees(r) is not None]
    print(f"  {k[0]}  {k[1]:4.0f}  {k[2]} {k[3]:4s} {len(v):3d}  {np.median(a) if a else float('nan'):8.4f}  "
          f"{np.median(n) if n else float('nan'):9.1f}  {sum(not basarili(r) for r in v):3d}")

print("\n== S2/S3: rejim etkisi (R5 / R0, D/r0 = 8, fe0) ==")
for kol, f, tur, ad in (("A0", ate, "ari", "ATE A0"), ("A1", ate, "ari", "ATE A1"), ("A0", anees, "geo", "ANEES A0"),
                        ("A1", anees, "geo", "ANEES A1"), ("A0g", anees, "geo", "ANEES A0g")):
    print(f"  {ad:10s}: {oran(('R5', 8.0, 'fe0', kol), ('R0', 8.0, 'fe0', kol), f, tur)}")
