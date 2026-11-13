#!/usr/bin/env python
"""Derin onarma kolu analizi (ON_KAYIT_FAZ2.md Sapma 3). Birim: sekans; ATE pencerelerde aritmetik, ANEES geometrik
ortalama (turb/faz2_analiz.py ile aynı kural); etki = birimlerde medyan eşli oran, bootstrap %95 GA (10 000).

    python onarim_analiz.py ontest      # evidence/ONARIM/on_test.jsonl  (G3-G9)
    python onarim_analiz.py             # evidence/ONARIM/records.jsonl  (HD1-HD4 + betimsel tablolar)

Sayısal ikiz denetimi (Sapma 3 eki, madde b): evidence/ONARIM/ikiz.jsonl varsa N′ kayıtları ("Ni" kipi) okunur;
N′ / N oranları aynı istatistikle raporlanır ve O / N, C / N kararlarının yanına N′ / N birim aralığı yazılır.
"""
import collections
import json
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(ROOT, "evidence", "ONARIM")
ONTEST = len(sys.argv) > 1 and sys.argv[1] == "ontest"
dosya = os.path.join(EV, "on_test.jsonl" if ONTEST else "records.jsonl")


def oku(p):
    out = {}
    with open(p) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            o = r.get("onarim") or {}
            if "arm" not in r:
                print("  !! hücre hatası:", r.get("key"), r.get("error"))
                continue
            out[(o["rejim"], o["mod"], r["arm"], o["sekans"], int(r["window"]))] = r        # son kayıt geçerli
    return out


def ana_izgara():
    A = {}
    with open(os.path.join(ROOT, "evidence", "FAZ2", "records.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            c = r.get("cell") or {}
            if c.get("seed") == 0 and r.get("arm") in ("A0", "A1", "A3") and float(c.get("D_over_r0", -1)) in (0.0, 8.0):
                A[(float(c["D_over_r0"]), r["arm"], c["sequence"], int(r["window"]))] = r
    return A


def sonda():
    S = {}
    p = os.path.join(ROOT, "evidence", "SONDA_KESIF", "records.jsonl")
    if os.path.exists(p):
        with open(p) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                c = r.get("cell") or {}
                if c.get("tag") == "sonda_R5_fe0" and float(c.get("D_over_r0", -1)) == 8.0 and "arm" in r:
                    S[(r["arm"], c["sequence"], int(r["window"]))] = r
    return S


def basarili(r):
    return r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged")


def ate(r):
    v = r.get("m_ate_se3")
    return float(v) if (r.get("status") == "ok" and r.get("m_valid") and v is not None and math.isfinite(v) and v > 0) else None


def anees(r):
    v = r.get("n_nees_pos_mean")
    return float(v) if (r.get("status") == "ok" and r.get("n_valid") and v is not None and math.isfinite(v) and v > 0) else None


S = oku(dosya)
IKIZ = os.path.join(EV, "ikiz.jsonl")
if not ONTEST and os.path.exists(IKIZ):
    S.update(oku(IKIZ))                     # kip "Ni": N dökümü, fx ve fy × (1 + 1e-12), A0 ve A1
ANA = ana_izgara()
print(f"[analiz] {dosya}: {len(S)} kayıt; durum {dict(collections.Counter(r.get('status') for r in S.values()))}; "
      f"code_fp {dict(collections.Counter(r.get('code_fp') for r in S.values()))}; "
      f"GPU {dict(collections.Counter(r.get('fe_gpu') for r in S.values()))}")
eksik = [k for k, r in S.items() if r.get("status") == "ok" and (r.get("m_ate_se3") is None or r.get("ms_ate_se3") is None
                                                                  or r.get("n_nees_pos_mean") is None)]
print(f"  alanı eksik kayıt: {len(eksik)}")


def kiyas(baslik, es):
    """es: [(etiket, r_yeni, r_basvuru)] -> ATE / ANEES göreli farkları."""
    print(f"\n== {baslik} ==")
    da, dn = [], []
    for et, a, b in es:
        if a is None or b is None:
            print(f"  {et}: kayıt yok")
            continue
        xa, xb, na, nb = ate(a), ate(b), anees(a), anees(b)
        fa = abs(xa / xb - 1) if xa and xb else float("nan")
        fn = abs(na / nb - 1) if na and nb else float("nan")
        da.append(fa); dn.append(fn)
        print(f"  {et}: ATE {xa:.5f} / {xb:.5f} m (fark %{100 * fa:.2f})   ANEES {na:.1f} / {nb:.1f} (fark %{100 * fn:.2f})"
              f"{'   BİT-AYNI' if a.get('m_ate_se3') == b.get('m_ate_se3') else ''}")
    if da:
        print(f"  -> ATE farkı medyan %{100 * np.nanmedian(da):.2f}, en çok %{100 * np.nanmax(da):.2f}; "
              f"ANEES farkı medyan %{100 * np.nanmedian(dn):.2f}, en çok %{100 * np.nanmax(dn):.2f}")


def fe_satir(et, r):
    f = r.get("frontend_stats") or {}
    print(f"  {et:28s} kp1 {f.get('avg_kp1_count', float('nan')):6.1f}  stereo_es {f.get('avg_stereo_matches', float('nan')):6.2f}  "
          f"ucgen {f.get('avg_current_stereo_triangulated_points', float('nan')):6.2f}  zamansal {f.get('avg_temporal_matches', float('nan')):6.1f}  "
          f"aktif_derinlik {f.get('avg_active_depth_tracks', float('nan')):6.1f}  t_fe {r.get('t_frontend', float('nan')):5.0f} s  "
          f"t_be {r.get('t_backend', float('nan')):5.0f} s")


if ONTEST:
    seq = "V1_01_easy"
    ws = sorted({k[4] for k in S})
    kiyas("G3: dosya yolu (R0, D/r0 8, N) / kanca yolu (ana ızgara tohum 0)",
          [(f"{a} w{w}", S.get(("R0", "N", a, seq, w)), ANA.get((8.0, a, seq, w))) for a in ("A0", "A1", "A3") for w in ws])
    kiyas("G4: temiz dosya yolu (D/r0 0, N) / ana ızgara D/r0 0 A0",
          [(f"A0 w{w}", S.get(("T", "N", "A0", seq, w)), ANA.get((0.0, "A0", seq, w))) for w in ws])
    kiyas("G5: temiz + onarım (O) / temiz dosya yolu (N)",
          [(f"A0 w{w}", S.get(("T", "O", "A0", seq, w)), S.get(("T", "N", "A0", seq, w))) for w in ws])
    a = [ate(S[("T", "O", "A0", seq, w)]) for w in ws if ("T", "O", "A0", seq, w) in S]
    b = [ate(S[("T", "N", "A0", seq, w)]) for w in ws if ("T", "N", "A0", seq, w) in S]
    if a and b and None not in a and None not in b:
        print(f"  -> sekans ATE oranı (O / N) = {np.mean(a) / np.mean(b):.3f}")
    print("\n== ön-uç istatistikleri ==")
    for et, k in (("R0 dr8 N (dosya)", ("R0", "N", "A0", seq, ws[0])), ("temiz N (dosya)", ("T", "N", "A0", seq, ws[0])),
                  ("temiz O (dosya)", ("T", "O", "A0", seq, ws[0]))):
        if k in S:
            fe_satir(et, S[k])
    for et, k in (("R0 dr8 (kanca, ana ızgara)", (8.0, "A0", seq, ws[0])), ("temiz (kanca, ana ızgara)", (0.0, "A0", seq, ws[0]))):
        if k in ANA:
            fe_satir(et, ANA[k])
    print("\n== G7: sigma_turb_used (A1) ve A3 sigma_hat ==")
    for k in sorted(S):
        if k[2] in ("A1", "A3") and k[4] == ws[0]:
            st = S[k].get("arm_stats") or {}
            print(f"  {k[0]} {k[1]} {k[2]}: sigma_turb_used_mean {st.get('sigma_turb_used_mean')}  sigma_hat_median {st.get('sigma_hat_median')}  "
                  f"t_be {S[k].get('t_backend', float('nan')):.0f} s")
    sys.exit(0)

# ----------------------------------------------------------------------------- ızgara
SND = sonda()
seqs = sorted({k[3] for k in S})
print(f"  sekans: {len(seqs)}; koşullar: {sorted({(k[0], k[1]) for k in S})}")

es = [(f"{s} {a} w{k[4]}", S[k], SND.get((a, s, k[4]))) for k in sorted(S) for s, a in [(k[3], k[2])]
      if k[0] == "R5" and k[1] == "N" and a in ("A0", "A1") and (a, s, k[4]) in SND]
if es:
    kiyas("R5 N (dosya yolu) / sonda R5 fe0 (kanca yolu)", es)
es = [(f"{k[3]} {k[2]} w{k[4]}", S[k], ANA.get((8.0, k[2], k[3], k[4]))) for k in sorted(S) if k[0] == "R0" and k[1] == "N"]
if es:
    kiyas("R0 N (dosya yolu) / ana ızgara tohum 0 (kanca yolu)", es)


def birim(rejim, mod, kol, f):
    g = collections.defaultdict(dict)
    for (rj, m, k, s, w), r in S.items():
        if (rj, m, k) == (rejim, mod, kol):
            v = f(r)
            if v is not None:
                g[s][w] = v
    return g


def temiz_birim(f):
    g = collections.defaultdict(dict)
    for (dr, k, s, w), r in ANA.items():
        if dr == 0.0 and k == "A0":
            v = f(r)
            if v is not None:
                g[s][w] = v
    return g


def oran(A, B, tur, ad, aralik=False):
    """-> None | dict(medyan, ga, n, enk, enb, tipik); aralik=True: birim oranlarının aralığı da yazdırılır."""
    lr = []
    for s in sorted(set(A) & set(B)):
        ws = sorted(set(A[s]) & set(B[s]))
        if not ws:
            continue
        xa, xb = np.array([A[s][w] for w in ws]), np.array([B[s][w] for w in ws])
        lr.append(float(np.mean(np.log(xa)) - np.mean(np.log(xb))) if tur == "geo" else float(np.log(np.mean(xa) / np.mean(xb))))
    if not lr:
        print(f"  {ad}: veri yok")
        return None
    lr = np.array(lr)
    rng = np.random.default_rng(0)
    bs = np.median(lr[rng.integers(0, len(lr), size=(10000, len(lr)))], axis=1)
    lo, hi = np.exp(np.percentile(bs, [2.5, 97.5]))
    o = dict(medyan=math.exp(np.median(lr)), ga=[float(lo), float(hi)], n=len(lr), enk=math.exp(lr.min()),
             enb=math.exp(lr.max()), tipik=math.exp(np.median(np.abs(lr))))
    print(f"  {ad}: {o['medyan']:.3f} [{lo:.3f}, {hi:.3f}]  (n = {len(lr)}, oran < 1: {int(np.sum(lr < 0))})"
          + (f"   birim aralığı [{o['enk']:.3f}, {o['enb']:.3f}], tipik çarpan ×/÷ {o['tipik']:.3f}" if aralik else ""))
    return o


print("\n== Öngörüler ==")
hd1 = oran(birim("R0", "O", "A0", ate), birim("R0", "N", "A0", ate), "ari", "HD1  R0 ATE  O/N (A0)        beklenen [0.8, 1.25]")
hd2a = oran(birim("R5", "O", "A0", ate), birim("R5", "N", "A0", ate), "ari", "HD2a R5 ATE  O/N (A0)        beklenen >= 0.7    ")
oran(birim("R5", "O", "A0", ate), temiz_birim(ate), "ari", "HD2b R5 ATE  O/temiz (A0)    beklenen >= 2      ")
oran(birim("R0", "O", "A1", anees), birim("R0", "O", "A0", anees), "geo", "HD3  R0 ANEES A1/A0 (O)      beklenen <= 0.33   ")
oran(birim("R5", "O", "A1", anees), birim("R5", "O", "A0", anees), "geo", "HD3  R5 ANEES A1/A0 (O)      beklenen <= 0.33   ")
oran(birim("R5", "C", "A1", anees), birim("R5", "C", "A0", anees), "geo", "HD3  R5 ANEES A1/A0 (C)      beklenen <= 0.33   ")
hd4 = oran(birim("R5", "C", "A0", ate), birim("R5", "N", "A0", ate), "ari", "HD4  R5 ATE  C/N (A0)        beklenen >= 0.8    ")

print("\n== Sayısal ikiz denetimi N′ (Sapma 3 eki b): aynı N dökümü, fx ve fy × (1 + 1e-12) ==")
IK = {}
for rj in ("R0", "R5"):
    for kol in ("A0", "A1"):
        IK[(rj, kol, "ate")] = oran(birim(rj, "Ni", kol, ate), birim(rj, "N", kol, ate), "ari",
                                    f"{rj} ATE   N′/N ({kol})", aralik=True)
        IK[(rj, kol, "anees")] = oran(birim(rj, "Ni", kol, anees), birim(rj, "N", kol, anees), "geo",
                                      f"{rj} ANEES N′/N ({kol})", aralik=True)
    oran(birim(rj, "N", "A1", anees), birim(rj, "N", "A0", anees), "geo", f"{rj} ANEES A1/A0 (N) ")
    oran(birim(rj, "Ni", "A1", anees), birim(rj, "Ni", "A0", anees), "geo", f"{rj} ANEES A1/A0 (N′)")
print("  karar notu (ön kayıt, Sapma 3 eki b): O/N ya da C/N medyanı, N′/N birim oranlarının [en küçük, en büyük]"
      " aralığının içindeyse karar 'sayısal gürültüden ayırt edilemiyor' notuyla raporlanır; eşikler değişmez")
for ad, o, ik in (("HD1  R0 ATE O/N (A0)", hd1, IK[("R0", "A0", "ate")]), ("HD2a R5 ATE O/N (A0)", hd2a, IK[("R5", "A0", "ate")]),
                  ("HD4  R5 ATE C/N (A0)", hd4, IK[("R5", "A0", "ate")])):
    if o and ik:
        ic = ik["enk"] <= o["medyan"] <= ik["enb"]
        print(f"  {ad}: medyan {o['medyan']:.3f}; N′/N birim aralığı [{ik['enk']:.3f}, {ik['enb']:.3f}] (n = {ik['n']}) -> "
              f"{'sayısal gürültüden AYIRT EDİLEMİYOR' if ic else 'aralığın dışında'}")
    else:
        print(f"  {ad}: veri yok")

print("\n== Betimsel oranlar ==")
for rj, mods in (("R0", ("O",)), ("R5", ("O", "C"))):
    for m in mods:
        for kol in ("A0", "A1", "A3"):
            oran(birim(rj, m, kol, ate), birim(rj, "N", kol, ate), "ari", f"{rj} ATE   {m}/N ({kol})")
            oran(birim(rj, m, kol, anees), birim(rj, "N", kol, anees), "geo", f"{rj} ANEES {m}/N ({kol})")
    for m in ("N",) + mods:
        oran(birim(rj, m, "A1", ate), birim(rj, m, "A0", ate), "ari", f"{rj} ATE   A1/A0 ({m})")
        oran(birim(rj, m, "A3", ate), birim(rj, m, "A0", ate), "ari", f"{rj} ATE   A3/A0 ({m})")
        oran(birim(rj, m, "A1", anees), birim(rj, m, "A0", anees), "geo", f"{rj} ANEES A1/A0 ({m})")
        oran(birim(rj, m, "A3", anees), birim(rj, m, "A0", anees), "geo", f"{rj} ANEES A3/A0 ({m})")
    oran(birim(rj, "N", "A0", ate), temiz_birim(ate), "ari", f"{rj} ATE   N/temiz (A0)")
oran(birim("R5", "N", "A0", ate), birim("R0", "N", "A0", ate), "ari", "ATE   R5/R0 (N, A0)")

print("\n== düzeyler: sekans medyanı ==")
print("  rejim mod kol   n_pencere  ATE_med[m]  ANEES_med  başarısız  sigma_hat_med(A3)")
G = collections.defaultdict(list)
for (rj, m, k, s, w), r in S.items():
    G[(rj, m, k)].append(r)
for k in sorted(G):
    v = G[k]
    a = [ate(r) for r in v if ate(r) is not None]
    n = [anees(r) for r in v if anees(r) is not None]
    sh = [(r.get("arm_stats") or {}).get("sigma_hat_median") for r in v]
    sh = [x for x in sh if x is not None and math.isfinite(x)]
    print(f"  {k[0]}  {k[1]}  {k[2]:3s}  {len(v):4d}  {np.median(a) if a else float('nan'):10.4f}  "
          f"{np.median(n) if n else float('nan'):9.1f}  {sum(not basarili(r) for r in v):4d}  "
          f"{np.median(sh) if sh else float('nan'):8.3f}")

print("\n== ön-uç istatistikleri (sekans ortalaması) ==")
F = collections.defaultdict(dict)
for (rj, m, k, s, w), r in S.items():
    if k == "A0" and r.get("frontend_stats"):
        F[(rj, m)][s] = r
for k in sorted(F):
    v = list(F[k].values())
    g = lambda n: float(np.mean([(x["frontend_stats"]).get(n, np.nan) for x in v]))      # noqa: E731
    print(f"  {k[0]} {k[1]}  n {len(v)}  kp1 {g('avg_kp1_count'):6.1f}  stereo_es {g('avg_stereo_matches'):6.2f}  "
          f"ucgen {g('avg_current_stereo_triangulated_points'):6.2f}  zamansal {g('avg_temporal_matches'):6.1f}  "
          f"aktif_derinlik {g('avg_active_depth_tracks'):6.1f}  t_fe {np.mean([x.get('t_frontend', np.nan) for x in v]):5.0f} s")

p = os.path.join(EV, "goruntu_olcu.jsonl")
if os.path.exists(p):
    print("\n== görüntü düzeyi (temize göre; cam0, 20 karede bir) ==")
    R = collections.defaultdict(list)
    with open(p) as f:
        for line in f:
            r = json.loads(line)
            R[(r["rejim"], r["mod"])].append(r)
    for k in sorted(R):
        v = R[k]
        print(f"  {k[0]} {k[1]}  n {len(v)}  PSNR {np.mean([x['psnr'] for x in v]):.2f} dB  SSIM {np.mean([x['ssim'] for x in v]):.3f}  "
              f"RMS akış {np.mean([x['rms_akis'] for x in v]):.3f} px  dOrt {np.mean([x['d_ort'] for x in v]):+.2f} DN")
