"""
Faz 2 rakip sonuçlarının eşli özeti (ON_KAYIT_FAZ2 §2: D/r0 etkisi = türbülanslı / temiz oranı, aynı
(sekans, tohum, pencere) üzerinden eşli; istatistik birimi (sekans, tohum), pencereler birim içinde ortalanır).

    python -m turb.rakip_analiz [--kayit evidence/FAZ2_rakip/rakip_kayitlari.jsonl] [--out ozet.json] \
        [--protokol onek] [--metrik m_ate_se3]

Kurallar:
  * Altyapı hatası taşıyan kayıtlar (GPU/TensorRT, OOM) dışarıda.
  * Tekrarlı koşumlar (ör. ORB-SLAM3 tohum başına 3): pencere değeri = başarılı tekrarların medyanı; pencere
    başarısız ⇔ tekrarların yarısından fazlası başarısız (ON_TEST_RAKIPLER.md §1 R5 kuralı).
  * Eşli oran (pencere düzeyi) r = ATE(D/r0) / ATE(0), ikisi de başarılıysa. Başarısızlık ayrıca oranlanır.
  * "başarısızlık dahil" medyan: temiz başarılı, türbülanslı başarısız pencereye r = +inf verilir (sıra tabanlı medyan).
  * Birim (sekans, tohum): pencere log-oranlarının ortalaması -> birim oranı exp(ort). Medyan, 10 000 örnekli
    bootstrap %95 GA, tek yönlü Wilcoxon (log-oran > 0) — keşifsel; ön kayıttaki Holm ailesi PyEco kolları içindir.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
VARSAYILAN = ROOT / "evidence" / "FAZ2_rakip" / "rakip_kayitlari.jsonl"


def oku(path, protokol):
    R = []
    for ln in open(path):
        try:
            r = json.loads(ln)
        except json.JSONDecodeError:
            continue
        if protokol and r.get("protokol") != protokol:
            continue
        if (r.get("kosum") or {}).get("altyapi_hatasi"):
            continue
        R.append(r)
    return R


def pencere_degerleri(R, metrik):
    """(rakip, sekans, dr, tohum, pencere) -> {'deger': medyan ya da None, 'basarisiz': bool, 'n_tekrar': n}"""
    g = defaultdict(list)
    for r in R:
        g[(r["rakip"], r["sequence"], float(r["D_over_r0"]), int(r["seed"]), int(r["window"]))].append(r)
    out = {}
    for k, rr in g.items():
        ok = [x[metrik] for x in rr if x["durum_p"] == "ok" and x.get(metrik) is not None]
        n_bas = sum(x["durum_p"] == "basarisiz" for x in rr)
        bas = n_bas > len(rr) / 2
        out[k] = {"deger": (float(np.median(ok)) if ok and not bas else None), "basarisiz": bas,
                  "n_tekrar": len(rr), "n_tekrar_basarisiz": n_bas,
                  "tekrar_cv": (float(np.std(ok, ddof=1) / np.mean(ok)) if len(ok) >= 2 else None)}
    return out


def bootstrap_ga(x, n=10000, seed=0):
    x = np.asarray(x, float)
    if len(x) < 2:
        return [None, None]
    rng = np.random.default_rng(seed)
    med = np.median(x[rng.integers(0, len(x), (n, len(x)))], axis=1)
    return [float(np.percentile(med, 2.5)), float(np.percentile(med, 97.5))]


def wilcoxon_tek(logr):
    try:
        from scipy.stats import wilcoxon
        if len(logr) < 3 or np.allclose(logr, 0):
            return None
        return float(wilcoxon(logr, alternative="greater").pvalue)
    except Exception:  # noqa: BLE001
        return None


def _sonlu(x):
    """JSON için: +inf -> "inf" (başarısız pencereler çoğunlukta demektir)."""
    return "inf" if math.isinf(x) else x


def ozet(R, metrik="m_ate_se3"):
    P = pencere_degerleri(R, metrik)
    rakipler = sorted({k[0] for k in P})
    drs = sorted({k[2] for k in P})
    res = {"metrik": metrik, "rakipler": {}}
    for rk in rakipler:
        d_rk = {"temiz": {}, "D_over_r0": {}}
        temiz = {k: v for k, v in P.items() if k[0] == rk and k[2] == 0}
        tv = [v["deger"] for v in temiz.values() if v["deger"] is not None]
        d_rk["temiz"] = {"n_pencere": len(temiz), "n_basarisiz": sum(v["basarisiz"] for v in temiz.values()),
                         "medyan": float(np.median(tv)) if tv else None,
                         "p10_p90": [float(np.percentile(tv, 10)), float(np.percentile(tv, 90))] if tv else None,
                         "tekrar_cv_medyan": (float(np.median([v["tekrar_cv"] for v in temiz.values()
                                                               if v["tekrar_cv"] is not None]))
                                              if any(v["tekrar_cv"] is not None for v in temiz.values()) else None)}
        for dr in drs:
            cells = {k: v for k, v in P.items() if k[0] == rk and k[2] == dr}
            ciftler, ciftler_inf, birim = [], [], defaultdict(list)
            for (_, s, _, sd, w), v in cells.items():
                t = P.get((rk, s, 0.0, sd, w))
                if t is None or t["deger"] is None:
                    continue
                if v["deger"] is not None:
                    r = v["deger"] / t["deger"]
                    ciftler.append({"sekans": s, "tohum": sd, "pencere": w, "oran": r})
                    ciftler_inf.append(r)
                    birim[(s, sd)].append(math.log(r))
                elif v["basarisiz"]:
                    ciftler_inf.append(math.inf)
            oran = [c["oran"] for c in ciftler]
            b = {k: math.exp(np.mean(v)) for k, v in birim.items()}
            bl = [math.log(x) for x in b.values()]
            n_bas = sum(v["basarisiz"] for v in cells.values())
            dd = {"n_pencere": len(cells), "n_basarisiz": n_bas,
                  "basarisizlik_orani": n_bas / len(cells) if cells else None,
                  "n_cift": len(oran),
                  "medyan_oran": float(np.median(oran)) if oran else None,
                  "medyan_oran_ga95": bootstrap_ga(oran) if oran else [None, None],
                  "medyan_oran_basarisiz_dahil": (_sonlu(float(np.median(ciftler_inf))) if ciftler_inf else None),
                  "medyan_ate": (float(np.median([v["deger"] for v in cells.values() if v["deger"] is not None]))
                                 if any(v["deger"] is not None for v in cells.values()) else None),
                  "birim_sekans_tohum": {"n": len(b), "medyan_oran": float(np.median(list(b.values()))) if b else None,
                                         "ga95": bootstrap_ga(list(b.values())) if b else [None, None],
                                         "wilcoxon_tek_p": wilcoxon_tek(bl) if dr > 0 else None,
                                         "oranlar": {f"{k[0]}/s{k[1]}": v for k, v in sorted(b.items())}},
                  "sekans_basina": {}}
            for s in sorted({c["sekans"] for c in ciftler} | {k[1] for k in cells}):
                os_ = [c["oran"] for c in ciftler if c["sekans"] == s]
                cs = {k: v for k, v in cells.items() if k[1] == s}
                dd["sekans_basina"][s] = {"n_cift": len(os_), "medyan_oran": float(np.median(os_)) if os_ else None,
                                          "n_basarisiz": sum(v["basarisiz"] for v in cs.values()),
                                          "n_pencere": len(cs)}
            dd["ciftler"] = ciftler
            d_rk["D_over_r0"][f"{dr:g}"] = dd
        res["rakipler"][rk] = d_rk
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kayit", default=str(VARSAYILAN))
    ap.add_argument("--out", default=None)
    ap.add_argument("--protokol", default="onek")
    ap.add_argument("--metrik", default="m_ate_se3", help="m_ate_se3 (nedensel, birincil) | ms_ate_se3 | m_rpe_trans")
    a = ap.parse_args()
    R = oku(a.kayit, a.protokol)
    res = {"kayit": a.kayit, "protokol": a.protokol, "n_kayit": len(R),
           "birincil": ozet(R, a.metrik), "ikincil_son_yorunge": ozet(R, "ms_ate_se3")}
    out = a.out or str(Path(a.kayit).with_name(Path(a.kayit).stem + "_ozet.json"))
    Path(out).write_text(json.dumps(res, indent=1, ensure_ascii=False))
    for rk, d in res["birincil"]["rakipler"].items():
        print(f"{rk}: temiz medyan {d['temiz']['medyan']}, başarısız {d['temiz']['n_basarisiz']}/{d['temiz']['n_pencere']}")
        for dr, dd in d["D_over_r0"].items():
            if dr != "0":
                print(f"   D/r0={dr}: medyan oran {dd['medyan_oran']} (GA {dd['medyan_oran_ga95']}), "
                      f"başarısız {dd['n_basarisiz']}/{dd['n_pencere']}, başarısız dahil {dd['medyan_oran_basarisiz_dahil']}")
    print("->", out)


if __name__ == "__main__":
    main()
