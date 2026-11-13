"""ON_KAYIT_FAZ2 v0.3 §4 analizi (HA2, HA3, HA5, keşif). Önceden kayıtlı testlerin birebir uygulaması.

    python analiz_v03.py                                  # ana ızgara: cikti_izgara/ana_v03/w*.jsonl
    python analiz_v03.py --girdi a.jsonl b.jsonl --etiket pilot   # sınama (pilot verisi)

Birim: (sekans, pencere, tohum) çifti; etki = ln(RPE1s_C0 / RPE1s_temiz).
- TOST ±ln(1.15), α=0.05: log-oran ortalamasının %90 t-GA'sı sınır içinde.
- Bozulma: tek yönlü Wilcoxon (log-oran > 0), Holm aile içinde (aile = sekans × C0>0).
- Medyan oran + 10 000 bootstrap %95 GA. Doz-yanıt: Spearman(log-oran, log10 C0), tek yönlü.
- HA5: bütün hücrelerde |Sim3 ölçek − 1| < 0.02.
"""
import argparse, json, math, sys
from pathlib import Path
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from rpe_ozet import ozet

HERE = Path(__file__).parent
SINIR = math.log(1.15)


def oku(dosyalar):
    satir, uyari = [], []
    for f in dosyalar:
        for l in open(f):
            r = json.loads(l)
            if r.get("durum") != "tamam":
                uyari.append(f"tamam değil: {r.get('turb_sekans')} {r.get('durum')}"); continue
            od = HERE / r["cikti_dizin"] if r.get("cikti_dizin") else HERE / "cikti_izgara" / r["turb_sekans"]
            traj = od / f"gtsam_traj_{r['turb_sekans']}.csv"
            pen = r["pencere"].replace(":", "-") if r.get("pencere") else None
            o = ozet(traj, r["sekans"], pen)
            if abs(o["ate"] - r["ate_se3"]) > 1e-6 * max(1.0, r["ate_se3"]):
                uyari.append(f"ezilmiş yörünge, dışlandı: {r['turb_sekans']}"); continue
            satir.append(dict(sekans=r["sekans"], pencere=r["pencere"], C0=float(r["C0"]), tohum=int(r["tohum"]),
                              gpu=str(r.get("gpu", "?")), olcek=r.get("sim3_olcek"), **o))
    return satir, uyari


def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); adj = np.empty(m); run = 0.0
    for i, j in enumerate(o):
        run = max(run, (m - i) * p[j]); adj[j] = min(1.0, run)
    return adj


def boot_med(x, B=10000, seed=0):
    rng = np.random.default_rng(seed)
    m = np.median(x[rng.integers(0, len(x), (B, len(x)))], 1)
    return np.percentile(m, [2.5, 97.5])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--girdi", nargs="*")
    ap.add_argument("--etiket", default="ana_v03")
    a = ap.parse_args()
    dosyalar = a.girdi or sorted((HERE / "cikti_izgara/ana_v03").glob("w*.jsonl"))
    satir, uyari = oku(dosyalar)
    for u in uyari:
        print("UYARI", u)
    # aile içinde tek GPU denetimi
    for sq in sorted({s["sekans"] for s in satir}):
        g = {s["gpu"] for s in satir if s["sekans"] == sq}
        print(f"GPU denetimi {sq}: {sorted(g)}", "✅" if len(g) == 1 else "❌ AİLE BİRDEN ÇOK GPU")
    idx = {(s["sekans"], s["pencere"], s["tohum"], s["C0"]): s for s in satir}
    sonuc = {"etiket": a.etiket, "testler": [], "uyari": uyari}
    for sq in sorted({s["sekans"] for s in satir}):
        c0s = sorted({s["C0"] for s in satir if s["sekans"] == sq and s["C0"] > 0})
        satirlar, plist = [], []
        tum_lr, tum_c = [], []
        for c in c0s:
            lr, lra = [], []
            for (q, pen, th, cc), t in idx.items():
                if q == sq and cc == c and (q, pen, th, 0.0) in idx:
                    k = idx[(q, pen, th, 0.0)]
                    lr.append(math.log(t["rpe1s"] / k["rpe1s"])); lra.append(math.log(t["ate"] / k["ate"]))
            lr = np.array(lr); n = len(lr)
            if n < 2:
                continue
            m, se = lr.mean(), lr.std(ddof=1) / math.sqrt(n)
            tq = stats.t.ppf(0.95, n - 1)
            ga90 = (m - tq * se, m + tq * se)
            p_tost = max(stats.t.sf((m + SINIR) / se, n - 1), stats.t.cdf((m - SINIR) / se, n - 1))
            p_w = stats.wilcoxon(lr, alternative="greater").pvalue if np.any(lr != 0) else 1.0
            bm = boot_med(np.exp(lr))
            satirlar.append(dict(sekans=sq, C0=c, n=n, medyan_oran=float(np.median(np.exp(lr))),
                                 boot95=[float(bm[0]), float(bm[1])], ort_logoran=float(m),
                                 ga90_oran=[math.exp(ga90[0]), math.exp(ga90[1])], p_tost=float(p_tost),
                                 esdeger=bool(ga90[0] > -SINIR and ga90[1] < SINIR), p_wilcoxon=float(p_w),
                                 oran_gt1=int((lr > 0).sum()), medyan_ate_oran=float(np.median(np.exp(lra)))))
            plist.append(p_w); tum_lr += list(lr); tum_c += [math.log10(c)] * n
        for s, ph in zip(satirlar, holm(plist) if plist else []):
            s["p_holm"] = float(ph); s["bozulma"] = bool(ph < 0.05)
        rho, p_rho = stats.spearmanr(tum_c, tum_lr) if len(set(tum_c)) > 1 else (float("nan"), float("nan"))
        p_rho1 = p_rho / 2 if rho > 0 else 1 - p_rho / 2
        sonuc["testler"] += satirlar
        sonuc[f"spearman_{sq}"] = dict(rho=float(rho), p_tek=float(p_rho1), n=len(tum_lr))
        print(f"\n{sq}  (Spearman ρ={rho:.3f}, tek yönlü p={p_rho1:.4f}, n={len(tum_lr)})")
        print(f"{'C0':>7} {'n':>3} {'med.oran':>8} {'boot95':>13} {'GA90 (TOST)':>13} {'eşd':>4} {'p_W':>7} {'p_Holm':>7} {'>1':>4} {'ATEmed':>6}")
        for s in satirlar:
            print(f"{s['C0']:7.0e} {s['n']:3d} {s['medyan_oran']:8.3f} {s['boot95'][0]:6.2f}–{s['boot95'][1]:<6.2f}"
                  f" {s['ga90_oran'][0]:6.3f}–{s['ga90_oran'][1]:<6.3f} {'E' if s['esdeger'] else '-':>4}"
                  f" {s['p_wilcoxon']:7.4f} {s['p_holm']:7.4f} {s['oran_gt1']:2d}/{s['n']:<2d} {s['medyan_ate_oran']:6.2f}")
    olc = [abs(s["olcek"] - 1) for s in satir if s["olcek"] is not None]
    sonuc["HA5"] = dict(n=len(olc), max_hata=float(max(olc)) if olc else None, gecti=bool(olc and max(olc) < 0.02))
    print(f"\nHA5 ölçek: n={len(olc)} en büyük |s−1| = {max(olc) if olc else float('nan'):.4f}", "✅" if sonuc["HA5"]["gecti"] else "❌")
    json.dump(sonuc, open(HERE / f"cikti_izgara/analiz_v03_{a.etiket}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
