"""Gerçek türbülans verisiyle (Saha ve ark. TSD) çevrimiçi DIMM kestiricisinin dış denetimi.

Ön bildirim: evidence/GERCEK_VERI/ON_BILDIRIM.md (yöntem ve ölçütler orada dondurulmuştur).
Kestirici kodu (turb/dimm.py) DEĞİŞTİRİLMEZ; burada yalnız çağrılır.

Kullanım (veri güvenilmez olduğundan -I ile):
  python -I gercek_veri_dogrula.py kos   --veri <saha_tsd/ham> --cikti <jsonl> [--lk 15] [--sinir N]
  python -I gercek_veri_dogrula.py ozet  --girdi <jsonl> --cikti <json>
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cv2  # noqa: E402
import numpy as np  # noqa: E402

from turb import optics  # noqa: E402
from turb.covariance import EMP_A, EMP_P  # noqa: E402
from turb.dimm import DifferentialMotionMonitor, r0_from_sigma_px  # noqa: E402

PFOV = 5.95e-6          # rad/piksel (Saha ve ark.)
LAM = 0.55e-6
D_ANA, D_DUY = 0.06, 0.03
KUME = {"2021-07": ("temmuz/July Dataset", 724.0), "2021-10": ("ekim/October Dataset", 1450.0)}
MAX_SEP = 180.0
N_BINS = 12


def dakikalar(veri):
    import pandas as pd
    d = pd.read_csv(os.path.join(veri, "Cn2 Values", "Both_Combined.csv"), index_col=0)
    d["rel"] = d["Nikon Crop"].str.extract(r"/(\d+NIKON/DSCN\d+\.JPG)")[0]
    d["num"] = d["rel"].str.extract(r"DSCN(\d+)")[0].astype(int)
    for t, g in d.groupby("Nikon Time", sort=True):
        yer, L = KUME[t[:7]]
        g = g.sort_values(["rel"])
        yield dict(zaman=t, ay=t[:7], L=L, cn2=float(g["Cn2"].iloc[0]),
                   cn2sig=float(g["Cn2Sig"].iloc[0]), ruzgar=float(g["WindSpeed"].iloc[0]),
                   yollar=[os.path.join(veri, yer, r) for r in g["rel"]])


def oku(yollar):
    ims = [cv2.imread(p, cv2.IMREAD_GRAYSCALE) for p in yollar]
    return [i for i in ims if i is not None and i.shape == (256, 256)]


def izle(ims, lk, ref):
    """ref = "ort": ön bildirimdeki ortalama görüntü; "kare": Değişiklik 1 (ortadaki kare)."""
    if len(ims) < 60:
        return None, dict(n_kare=len(ims), n_ozn=0)
    if ref == "ort":
        ort = np.mean(np.stack(ims).astype(np.float32), axis=0)
        ort8 = np.clip(ort + 0.5, 0, 255).astype(np.uint8)
        duzey = 2
    else:
        ort8 = ims[len(ims) // 2]
        duzey = 3
    p0 = cv2.goodFeaturesToTrack(ort8, maxCorners=300, qualityLevel=0.02, minDistance=5)
    if p0 is None or len(p0) < 30:
        return None, dict(n_kare=len(ims), n_ozn=0 if p0 is None else len(p0))
    kw = dict(winSize=(lk, lk), maxLevel=duzey,
              criteria=(cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 30, 0.01))
    P = np.full((len(ims), len(p0), 2), np.nan)
    for k, im in enumerate(ims):
        p1, st, _ = cv2.calcOpticalFlowPyrLK(ort8, im, p0, None, **kw)
        pb, sb, _ = cv2.calcOpticalFlowPyrLK(im, ort8, p1, None, **kw)
        fb = np.linalg.norm((pb - p0).reshape(-1, 2), axis=1)
        ok = (st.ravel() == 1) & (sb.ravel() == 1) & (fb < 0.5)
        P[k, ok] = p1.reshape(-1, 2)[ok]
    tut = np.mean(np.isfinite(P[:, :, 0]), axis=0) >= 0.9
    P = P[:, tut]
    return P, dict(n_kare=len(ims), n_ozn=int(tut.sum()))


def model(sep, s2, tc):
    return 4.0 * s2 * (1.0 - EMP_A * np.exp(-((np.maximum(sep, 0) / max(tc, 1e-6)) ** EMP_P)))


def isle(dk, ims, lk, ref):
    P, bilgi = izle(ims, lk, ref)
    kay = dict(zaman=dk["zaman"], ay=dk["ay"], L=dk["L"], cn2=dk["cn2"], cn2sig=dk["cn2sig"],
               ruzgar=dk["ruzgar"], lk=lk, ref=ref, **bilgi)
    if P is None or bilgi["n_ozn"] < 30:
        kay["durum"] = "dislandi"
        return kay
    ort = np.nanmean(P, axis=0)
    R = P - ort[None]
    dimm = DifferentialMotionMonitor(px_per_rad=1.0 / PFOV, aperture_D=D_ANA, wavelength=LAM,
                                     n_bins=N_BINS, max_sep_px=MAX_SEP, window_frames=len(P))
    for k in range(len(P)):
        ok = np.isfinite(R[k, :, 0])
        dimm.add_frame(ort[ok], R[k, ok])
    s, tc, r0 = dimm.estimate()
    sep, Dv, C = dimm.structure_function()
    kay.update(sf_sep=[float(x) for x in sep], sf_D=[None if not np.isfinite(x) else float(x) for x in Dv],
               sf_n=[float(x) for x in C], rms_px=float(np.sqrt(np.nanmean(R ** 2))))
    if not np.isfinite(s):
        kay["durum"] = "dislandi"
        return kay
    ok = np.isfinite(Dv) & (C > 10)
    w = C[ok]
    y, yh = Dv[ok], model(sep[ok], s * s, tc)
    ybar = np.sum(w * y) / np.sum(w)
    r2 = 1.0 - np.sum(w * (y - yh) ** 2) / max(np.sum(w * (y - ybar) ** 2), 1e-30)
    r0d = r0_from_sigma_px(s, D_DUY, LAM, 1.0 / PFOV)
    r0s = optics.r0_from_cn2(dk["cn2"], dk["L"], LAM)
    kay.update(durum="ok", sigma_px=float(s), tc_px=float(tc), yapisik=bool(dimm.last_pinned),
               r0_m=float(r0), r2=float(r2),
               cn2_kest=float(optics.cn2_from_r0(r0, dk["L"], LAM)),
               cn2_kest_D003=float(optics.cn2_from_r0(r0d, dk["L"], LAM)),
               sigma_kuram_duzlem_px=float(optics.sigma_tilt_px(D_ANA / r0s, LAM, D_ANA, 1.0 / PFOV)))
    return kay


def spearman(a, b):
    from scipy.stats import rankdata
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def blok_ozyukleme(n, rng, blok=10):
    bas = rng.integers(0, max(n - blok + 1, 1), size=int(np.ceil(n / blok)))
    return np.concatenate([np.arange(b, min(b + blok, n)) for b in bas])[:n]


def ozetle(kayitlar, etiket):
    tum = [k for k in kayitlar if k["ay"] == etiket]
    ok = sorted([k for k in tum if k["durum"] == "ok"], key=lambda k: k["zaman"])
    cik = dict(n_dakika=len(tum), n_ok=len(ok), n_dislanan=len(tum) - len(ok))
    if len(ok) < 20:
        cik["not"] = "yetersiz dakika"
        return cik
    s2 = np.array([k["sigma_px"] ** 2 for k in ok])
    c = np.array([k["cn2"] for k in ok])
    ck = np.array([k["cn2_kest"] for k in ok])
    ck3 = np.array([k["cn2_kest_D003"] for k in ok])
    yap = np.array([k["yapisik"] for k in ok])

    def olc(ix):
        egim = float(np.polyfit(np.log(c[ix]), np.log(s2[ix]), 1)[0])
        return spearman(s2[ix], c[ix]), egim, float(np.median(ck[ix] / c[ix]))

    rng = np.random.default_rng(0)
    B = np.array([olc(blok_ozyukleme(len(ok), rng)) for _ in range(2000)])
    ga = lambda j: [float(x) for x in np.percentile(B[:, j], [2.5, 97.5])]  # noqa: E731
    g1, g2, g3 = olc(np.arange(len(ok)))
    r2 = float(np.median([k["r2"] for k in ok]))
    cik.update(
        G1_spearman=g1, G1_ga=ga(0), G1_gecer=bool(g1 >= 0.5 and ga(0)[0] > 0),
        G2_egim=g2, G2_ga=ga(1), G2_gecer=bool(0.5 <= g2 <= 1.3),
        G3_oran=g3, G3_ga=ga(2), G3_gecer=bool(0.12 <= g3 <= 1.2),
        G4_r2_ortanca=r2, G4_yapisik_oran=float(yap.mean()),
        G4_tc_ortanca_px=float(np.median([k["tc_px"] for k in ok])),
        G4_gecer=bool(r2 >= 0.8 and yap.mean() <= 0.5),
        sigma_px_ceyrek=[float(x) for x in np.percentile(np.sqrt(s2), [5, 25, 50, 75, 95])],
        sigma_kuram_kuresel_ortanca_px=float(np.median([k["sigma_kuram_duzlem_px"] for k in ok]) * np.sqrt(3 / 8)),
        n_ozn_ortanca=float(np.median([k["n_ozn"] for k in ok])),
        duy_D003_oran=float(np.median(ck3 / c)),
        duy_yapisiksiz=dict(n=int((~yap).sum()), **(dict(zip(("spearman", "egim", "oran"), olc(np.where(~yap)[0])))
                                                     if (~yap).sum() >= 20 else {})),
        duy_rms_spearman=spearman(np.array([k["rms_px"] for k in ok]), c),
    )
    return cik


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kip", choices=["kos", "ozet"])
    ap.add_argument("--veri")
    ap.add_argument("--girdi")
    ap.add_argument("--cikti", required=True)
    ap.add_argument("--ay", default="")
    ap.add_argument("--sinir", type=int, default=0, help="duman testi: ay başına eşit aralıklı N dakika")
    a = ap.parse_args()
    if a.kip == "kos":
        dks = [d for d in dakikalar(a.veri) if not a.ay or d["ay"] == a.ay]
        if a.sinir:
            sec = []
            for ay in KUME:
                x = [d for d in dks if d["ay"] == ay]
                sec += [x[i] for i in np.linspace(0, len(x) - 1, a.sinir).astype(int)]
            dks = sec
        with open(a.cikti, "w") as f:
            for i, dk in enumerate(dks):
                ims = oku(dk["yollar"])
                for ref, lk in (("ort", 15), ("kare", 15), ("ort", 21)):
                    k = isle(dk, ims, lk, ref)
                    f.write(json.dumps(k, ensure_ascii=False) + "\n")
                    f.flush()
                    print(f"{i + 1}/{len(dks)} {k['zaman']} {ref}{lk} {k['durum']} kare={k['n_kare']} ozn={k['n_ozn']}", flush=True)
    else:
        kay = [json.loads(s) for g in a.girdi.split(",") for s in open(g)]
        cik = {}
        for ref, lk in (("ort", 15), ("kare", 15), ("ort", 21)):
            alt = [k for k in kay if k["ref"] == ref and k["lk"] == lk]
            cik[f"{ref}{lk}"] = {ay: ozetle(alt, ay) for ay in KUME}
        json.dump(cik, open(a.cikti, "w"), indent=1, ensure_ascii=False)
        print(json.dumps(cik, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
