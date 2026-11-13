"""H2 (ON_KAYIT_FAZ2 §3) analizi; işlemsel tanım ON_TEST_FAZ2_EUROC.md §9 (sonuçlara bakılmadan yazıldı).

    python -m turb.h2_analiz                 # karar: ana ızgarada D/r0 0/8/16 A0 tamamsa
    python -m turb.h2_analiz --durum         # yalnız eksik sayıları (sayı basmaz)
    python -m turb.h2_analiz --oz_test       # sentetik veriyle hesap denetimi

Koşullar: tilt/blur = FAZ2_FAKT; ikisi = ana ızgara A0 (speckle_global + tilt); yok = ana ızgara D/r0=0 A0.
Birim (sekans, tohum); birim değeri = geçerli pencere ATE'lerinin (m_ate_se3) ortalaması.
pay_c = medyan_u(ATE_c - ATE_yok) / medyan_u(ATE_ikisi - ATE_yok); %95 GA birim bootstrap (10 000, tohum 0).
Karar: pay_blur >= .75 ve pay_tilt <= .25 -> tuttu; pay_tilt > .50 -> çürüdü; aksi kısmen/kararsız.
"""
import argparse, json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
KOSUL = ("yok", "tilt", "blur", "ikisi")


def kosul(c):
    dr, b, t = float(c["D_over_r0"]), c.get("blur"), bool(c.get("tilt"))
    if dr == 0.0:
        return 0.0, "yok"
    if b == "speckle_global" and t:
        return dr, "ikisi"
    if b == "none" and t:
        return dr, "tilt"
    if b == "speckle_global" and not t:
        return dr, "blur"
    return dr, None


def oku(ana, fakt):
    """-> P[(dr, kosul, sekans, tohum)] = {pencere: ate|None}, beklenen pencereler W[(sekans, tohum)]"""
    P, W = defaultdict(dict), {}
    for f, kaynak in ((ana, "ana"), (fakt, "fakt")):
        for l in open(f):
            r = json.loads(l)
            if r.get("arm") != "A0":
                continue
            c = r["cell"]
            dr, k = kosul(c)
            if k is None or (kaynak == "ana" and k in ("tilt", "blur")) or (kaynak == "fakt" and k not in ("tilt", "blur")):
                continue
            if k != "yok" and dr not in (8.0, 16.0):
                continue
            u = (c["sequence"], int(c["seed"]))
            W[u] = tuple(c["windows"])
            v = r.get("m_ate_se3")
            ok = r.get("status") == "ok" and r.get("m_valid") and v is not None and math.isfinite(v)
            P[(dr, k, *u)][int(r["window"])] = float(v) if ok else None
    return P, W


def birimler(P, W, dr):
    deg, eksik, gecersiz = {}, [], 0
    for u, wins in sorted(W.items()):
        d = {}
        for k in KOSUL:
            p = P.get((0.0 if k == "yok" else dr, k, *u), {})
            if any(w not in p for w in wins):
                eksik.append((u, k)); break
            vals = [p[w] for w in wins if p[w] is not None]
            gecersiz += sum(p[w] is None for w in wins)
            if not vals:
                eksik.append((u, k)); break
            d[k] = float(np.mean(vals))
        else:
            deg[u] = d
    return deg, eksik, gecersiz


def paylar(deg, B=10000, seed=0):
    U = sorted(deg)
    D = {k: np.array([deg[u][k] - deg[u]["yok"] for u in U]) for k in ("tilt", "blur", "ikisi")}

    def pay(idx):
        m = np.median(D["ikisi"][idx])
        return np.median(D["blur"][idx]) / m, np.median(D["tilt"][idx]) / m

    tum = np.arange(len(U))
    pb, pt = pay(tum)
    rng = np.random.default_rng(seed)
    bs = np.array([pay(rng.integers(0, len(U), len(U))) for _ in range(B)])
    ga = lambda x: [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]
    if pb >= 0.75 and pt <= 0.25:
        karar = "TUTTU"
    elif pt > 0.50:
        karar = "ÇÜRÜDÜ (yalnız-tilt payı > %50)"
    else:
        karar = "KISMEN / KARARSIZ"
    lr = {k: np.log(np.array([deg[u][k] / deg[u]["yok"] for u in U])) for k in ("tilt", "blur", "ikisi")}
    dlt = lr["tilt"] - lr["blur"]
    p_tb = float(stats.wilcoxon(dlt).pvalue) if np.any(dlt != 0) else 1.0
    return dict(n_birim=len(U), pay_blur=float(pb), pay_blur_ga95=ga(bs[:, 0]), pay_tilt=float(pt),
                pay_tilt_ga95=ga(bs[:, 1]), etkilesim=float(pb + pt - 1),
                medyan_artis={k: float(np.median(v)) for k, v in D.items()},
                medyan_oran={k: float(np.exp(np.median(v))) for k, v in lr.items()},
                wilcoxon_tilt_vs_blur_p=p_tb, karar=karar)


def oz_test():
    rng = np.random.default_rng(1)
    deg = {}
    for i in range(45):  # blur ikisi'nin %80'i, tilt %10'u; beklenen pay_blur≈0.8, pay_tilt≈0.1, TUTTU
        y = 0.02 + 0.005 * rng.random(); a = 0.10 * (1 + 0.1 * rng.standard_normal())
        deg[("S", i)] = dict(yok=y, ikisi=y + a, blur=y + 0.8 * a, tilt=y + 0.1 * a)
    s = paylar(deg, B=500)
    assert abs(s["pay_blur"] - 0.8) < 0.05 and abs(s["pay_tilt"] - 0.1) < 0.05 and s["karar"] == "TUTTU", s
    for k in deg:
        deg[k]["tilt"] = deg[k]["yok"] + 0.6 * (deg[k]["ikisi"] - deg[k]["yok"])
    s = paylar(deg, B=500)
    assert s["karar"].startswith("ÇÜRÜDÜ") and abs(s["pay_tilt"] - 0.6) < 0.01, s
    # okuma/eşleme: geçici kayıtlarla
    import tempfile
    td = Path(tempfile.mkdtemp())
    ana, fakt = td / "ana.jsonl", td / "fakt.jsonl"
    with open(ana, "w") as fa, open(fakt, "w") as ff:
        for (dr, b, t, f, ate) in [(0, "speckle_global", True, fa, 1.0), (16, "speckle_global", True, fa, 3.0),
                                   (16, "none", True, ff, 1.5), (16, "speckle_global", False, ff, 2.5),
                                   (16, "speckle_global", True, fa, 9.0)]:
            for w in (10, 20):
                arm = "A0" if ate != 9.0 else "A3"  # A3 kaydı yok sayılmalı
                f.write(json.dumps(dict(arm=arm, status="ok", m_valid=True, m_ate_se3=ate + (w == 20) * 0.0, window=w,
                                        cell=dict(D_over_r0=dr, blur=b, tilt=t, sequence="X", seed=0, windows=[10, 20]))) + "\n")
    P, W = oku(ana, fakt)
    deg, eksik, gec = birimler(P, W, 16.0)
    assert deg == {("X", 0): dict(yok=1.0, tilt=1.5, blur=2.5, ikisi=3.0)} and not eksik and gec == 0, (deg, eksik)
    print("öz-test GEÇTİ")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ana", default=str(ROOT / "evidence/FAZ2/records.jsonl"))
    ap.add_argument("--fakt", default=str(ROOT / "evidence/FAZ2_FAKT/records.jsonl"))
    ap.add_argument("--cikti", default=str(ROOT / "evidence/FAZ2_FAKT/h2_analiz.json"))
    ap.add_argument("--durum", action="store_true")
    ap.add_argument("--oz_test", action="store_true")
    a = ap.parse_args()
    if a.oz_test:
        oz_test(); return 0
    P, W = oku(a.ana, a.fakt)
    out, tam = {"beklenen_birim": 45}, True
    for dr in (16.0, 8.0):
        deg, eksik, gec = birimler(P, W, dr)
        eksik_kosul = defaultdict(int)
        for _, k in eksik:
            eksik_kosul[k] += 1
        print(f"D/r0={dr:.0f}: tam birim {len(deg)}/{len(W)} (beklenen 45); eksik ilk koşul {dict(eksik_kosul)}; geçersiz pencere {gec}")
        if len(deg) < 45 or len(W) < 45:
            tam = False
        out[f"dr{dr:.0f}"] = dict(n_tam=len(deg), eksik=[[*u, k] for u, k in eksik], gecersiz_pencere=gec, deg=deg)
    if a.durum or not tam:
        if not tam:
            print("KARAR YOK: birimler eksik (ana ızgara D/r0 0/8/16 A0 bekleniyor)")
        return 0 if a.durum else 2
    for dr in (16.0, 8.0):
        s = paylar(out[f"dr{dr:.0f}"]["deg"])
        sek = {}
        for q in sorted({u[0] for u in out[f"dr{dr:.0f}"]["deg"]}):
            sd = {u: v for u, v in out[f"dr{dr:.0f}"]["deg"].items() if u[0] == q}
            sek[q] = {k: s2 for k, s2 in paylar(sd, B=1000).items() if k in ("pay_blur", "pay_tilt", "n_birim")}
        s["sekans"] = sek
        out[f"dr{dr:.0f}"].update(s)
        out[f"dr{dr:.0f}"]["deg"] = {f"{u[0]}/s{u[1]}": v for u, v in out[f"dr{dr:.0f}"]["deg"].items()}
        print(f"\nD/r0={dr:.0f} ({'birincil' if dr == 16 else 'ikincil'}): n={s['n_birim']}  pay_blur={s['pay_blur']:.3f} "
              f"[{s['pay_blur_ga95'][0]:.2f}, {s['pay_blur_ga95'][1]:.2f}]  pay_tilt={s['pay_tilt']:.3f} "
              f"[{s['pay_tilt_ga95'][0]:.2f}, {s['pay_tilt_ga95'][1]:.2f}]  etkileşim={s['etkilesim']:+.3f}")
        print(f"  medyan oran (c/yok): {', '.join(f'{k} {v:.2f}' for k, v in s['medyan_oran'].items())}; "
              f"Wilcoxon tilt vs blur p={s['wilcoxon_tilt_vs_blur_p']:.2g}")
        print(f"  KARAR ({'H2' if dr == 16 else 'bilgi'}): {s['karar']}")
    json.dump(out, open(a.cikti, "w"), indent=1, ensure_ascii=False)
    print(f"-> {a.cikti}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
