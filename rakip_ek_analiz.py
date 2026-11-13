#!/usr/bin/env python
"""Sapma 2 rakip analizleri (ON_KAYIT_FAZ2 §5): HR1 (rakipler D/r0 = 4'te) ve HI1 (ORB-SLAM3 stereo-ataletsel).
Tanımlar ön kayıttaki gibi; pencere/birim kuralları turb.rakip_analiz ile aynı (işlevleri oradan alınır).
  HR1: eşli medyan ATE oranı (4 / 0; aynı sekans, tohum, pencere; ikisinde de başarılı) >= 2, üç rakipten en az ikisinde.
       Metrik: nedensel yörünge (m_ate_se3; ana rakip analizinin birincil metriği).
  HI1: son yörünge (ms_ate_se3; başarısızlık ms_durum_p — Sapma 2 eki): D/r0 = 16'da SI başarısız pencere oranı <
       stereo ORB-SLAM3'ünki; SI eşli medyan ATE oranı (16 / 0) >= 2.
Betimsel (Holm ailesi yok). Çıktı: evidence/FAZ2_rakip_ek/rakip_ek_analiz.json"""
import json, sys, collections
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from turb.rakip_analiz import oku, ozet, pencere_degerleri  # noqa: E402
from turb.faz2_analiz import oran_testi  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

EV = ROOT / "evidence"
ANA, EK = EV / "FAZ2_rakip" / "rakip_kayitlari.jsonl", EV / "FAZ2_rakip_ek" / "rakip_kayitlari.jsonl"
SI = [EV / f"FAZ2_rakip_si_t{t}" / "rakip_kayitlari.jsonl" for t in (0, 1, 2)]


def son(R):
    """Son yörünge görünümü: başarısızlık ölçütü ms_durum_p."""
    out = []
    for r in R:
        q = dict(r); q["durum_p"] = r.get("ms_durum_p", "basarisiz"); out.append(q)
    return out


def kisa(d):
    b = d["birim_sekans_tohum"]
    return dict(n_pencere=d["n_pencere"], n_basarisiz=d["n_basarisiz"], basarisizlik=d["basarisizlik_orani"],
                n_cift=d["n_cift"], medyan_oran=d["medyan_oran"], ga95=d["medyan_oran_ga95"],
                medyan_oran_basarisiz_dahil=d["medyan_oran_basarisiz_dahil"], medyan_ate=d["medyan_ate"],
                birim_n=b["n"], birim_medyan_oran=b["medyan_oran"], birim_ga95=b["ga95"], birim_wilcoxon_tek_p=b["wilcoxon_tek_p"])


def birim_aritmetik(R, rk, dr, metrik="m_ate_se3"):
    """PyEco ana analizindeki birim tanımı (turb.faz2_analiz.eslesik): birim (sekans, tohum); birim değeri = iki düzeyde
    de başarılı pencerelerin ATE'sinin aritmetik ortalaması; oran = birim(dr) / birim(0); medyan + bootstrap + Wilcoxon."""
    P = pencere_degerleri(R, metrik)
    g = collections.defaultdict(lambda: ([], []))
    for (r_, s, d, sd, w), v in P.items():
        if r_ != rk or d != dr or v["deger"] is None:
            continue
        t = P.get((rk, s, 0.0, sd, w))
        if t is None or t["deger"] is None:
            continue
        g[(s, sd)][0].append(v["deger"]); g[(s, sd)][1].append(t["deger"])
    a = np.array([np.mean(v[0]) for v in g.values()]); b = np.array([np.mean(v[1]) for v in g.values()])
    if len(a) < 3:
        return None
    t = oran_testi(a, b)
    return {k: t[k] for k in ("n", "medyan_oran", "ga95")} | {"p": t.get("p", t.get("wilcoxon_p"))}


def sayim(R):
    c = collections.Counter((r["rakip"], float(r["D_over_r0"])) for r in R)
    return {f"{k[0]}@{k[1]:g}": v for k, v in sorted(c.items())}


def main():
    Ra, Re = oku(ANA, "onek"), oku(EK, "onek")
    Rs = [r for p in SI for r in oku(p, "onek")]
    res = {"kayit_sayimi": {"ana": sayim(Ra), "ek": sayim(Re), "si": sayim(Rs)},
           "kimlik_fp": {"ana": sorted({r.get("kimlik_fp") for r in Ra}), "ek": sorted({r.get("kimlik_fp") for r in Re}),
                         "si": sorted({r.get("kimlik_fp") for r in Rs})}}
    # ---- HR1
    def hr1(R, ad):
        o = ozet(R, "m_ate_se3")["rakipler"]
        t = {rk: {dr: kisa(o[rk]["D_over_r0"][dr]) for dr in ("4", "8", "16", "24") if dr in o[rk]["D_over_r0"]}
             | {"temiz": o[rk]["temiz"]} for rk in ("orbslam3", "airslam", "droid") if rk in o}
        for rk in t:
            for dr in ("4", "8", "16", "24"):
                if dr in t[rk]:
                    t[rk][dr]["birim_aritmetik"] = birim_aritmetik(R, rk, float(dr))
        n2 = [rk for rk in t if t[rk]["4"]["medyan_oran"] is not None and t[rk]["4"]["medyan_oran"] >= 2]
        n2b = [rk for rk in t if (t[rk]["4"]["birim_medyan_oran"] or 0) >= 2]
        n2c = [rk for rk in t if ((t[rk]["4"]["birim_aritmetik"] or {}).get("medyan_oran") or 0) >= 2]
        return {"ad": ad, "rakipler": t, "oran>=2_olanlar": n2, "karar": "TUTTU" if len(n2) >= 2 else "ÇÜRÜDÜ",
                "toplama_duzeyine_gore": {"pencere çiftleri": n2, "birim (pencere log-oran ortalaması)": n2b,
                                          "birim (ATE aritmetik ortalaması; PyEco ana analiziyle aynı)": n2c}}
    R1 = [r for r in Ra + Re if r["rakip"] in ("orbslam3", "airslam", "droid")]
    res["HR1"] = hr1(R1, "ön kayıtlı")
    res["HR1_duyarlilik_airslam_V2_01_w100_haric"] = hr1(
        [r for r in R1 if not (r["rakip"] == "airslam" and r["sequence"].startswith("V2_01") and int(r["window"]) == 100)],
        "AirSLAM V2_01 / pencere 100 hariç (sonradan duyarlılık)")
    res["HR1_son_yorunge"] = {rk: kisa(v["D_over_r0"]["4"]) for rk, v in ozet(son(R1), "ms_ate_se3")["rakipler"].items()}
    # ---- HI1
    St = [r for r in Ra if r["rakip"] == "orbslam3"] + [r for r in Re if r["rakip"] == "orbslam3"]
    o_si, o_st = ozet(son(Rs), "ms_ate_se3")["rakipler"]["orbslam3_si"], ozet(son(St), "ms_ate_se3")["rakipler"]["orbslam3"]
    hi = {"si": {dr: kisa(d) for dr, d in o_si["D_over_r0"].items()}, "stereo": {dr: kisa(d) for dr, d in o_st["D_over_r0"].items()},
          "si_temiz": o_si["temiz"], "stereo_temiz": o_st["temiz"]}
    f_si, f_st, oran16 = hi["si"]["16"]["basarisizlik"], hi["stereo"]["16"]["basarisizlik"], hi["si"]["16"]["medyan_oran"]
    hi["karar"] = {"basarisizlik_si_16": f_si, "basarisizlik_stereo_16": f_st, "ilk_yari (SI < stereo)": bool(f_si < f_st),
                   "si_oran_16_0": oran16, "si_oran_ga95": hi["si"]["16"]["ga95"], "ikinci_yari (oran >= 2)": bool(oran16 >= 2),
                   "HI1": "TUTTU" if (f_si < f_st and oran16 >= 2) else "KISMEN / ÇÜRÜDÜ"}
    # aynı pencerelerde eşli başarısızlık (son yörünge), düzey başına; temizde başarısız SI pencereleri (yapısal)
    Psi, Pst = pencere_degerleri(son(Rs), "ms_ate_se3"), pencere_degerleri(son(St), "ms_ate_se3")
    yap = sorted({(k[1], k[4]) for k, v in Psi.items() if k[2] == 0.0 and v["basarisiz"]})
    hi["si_temizde_basarisiz_pencereler"] = [f"{s}/w{w}" for s, w in yap]
    es = {}
    for dr in (0.0, 4.0, 8.0, 16.0, 24.0):
        a = {(k[1], k[3], k[4]): v["basarisiz"] for k, v in Psi.items() if k[2] == dr}
        b = {(k[1], k[3], k[4]): v["basarisiz"] for k, v in Pst.items() if k[2] == dr}
        ort = sorted(set(a) & set(b))
        ty = [k for k in ort if (k[0], k[2]) not in yap]
        es[f"{dr:g}"] = dict(n_ortak=len(ort), ikisi_de=sum(a[k] and b[k] for k in ort), yalniz_si=sum(a[k] and not b[k] for k in ort),
                             yalniz_stereo=sum(b[k] and not a[k] for k in ort), si=sum(a[k] for k in ort), stereo=sum(b[k] for k in ort),
                             yapisal_haric_n=len(ty), yapisal_haric_si=sum(a[k] for k in ty), yapisal_haric_stereo=sum(b[k] for k in ty))
        # sonradan eklenen betimsel sınama (ön kayıtta yok): eşli başarısızlık farkı, kesin McNemar (uyumsuz çiftler üzerinde binom)
        for ad_, kume in (("mcnemar_kesin_p (sonradan, betimsel)", ort), ("mcnemar_kesin_p_yapisal_haric (sonradan, betimsel)", ty)):
            y_si = sum(a[k] and not b[k] for k in kume); y_st = sum(b[k] and not a[k] for k in kume)
            es[f"{dr:g}"][ad_] = float(binomtest(y_si, y_si + y_st, 0.5).pvalue) if y_si + y_st else None
            es[f"{dr:g}"][ad_.replace("mcnemar_kesin_p", "uyumsuz_si_stereo")] = [y_si, y_st]
    hi["esli_basarisizlik_son_yorunge"] = es
    o_n = ozet(Rs, "m_ate_se3")["rakipler"]["orbslam3_si"]
    hi["si_nedensel (ikincil)"] = {dr: kisa(d) for dr, d in o_n["D_over_r0"].items()}
    res["HI1"] = hi
    out = EV / "FAZ2_rakip_ek" / "rakip_ek_analiz.json"
    out.write_text(json.dumps(res, indent=1, ensure_ascii=False))
    f = lambda x: "—" if x is None else (x if isinstance(x, str) else f"{x:.3g}")  # noqa: E731
    print("kayıt:", res["kayit_sayimi"]); print("kimlik_fp:", res["kimlik_fp"])
    for ad in ("HR1", "HR1_duyarlilik_airslam_V2_01_w100_haric"):
        h = res[ad]; print(f"\n{ad} ({h['ad']}): pencere düzeyi {h['karar']}; oran ≥ 2 olan rakipler: {h['toplama_duzeyine_gore']}")
        for rk, t in h["rakipler"].items():
            for dr in ("4", "8", "16", "24"):
                if dr in t:
                    d = t[dr]
                    print(f"  {rk:9s} D/r0 {dr:>2}: oran {f(d['medyan_oran'])} [{f(d['ga95'][0])}, {f(d['ga95'][1])}] n_çift {d['n_cift']}; "
                          f"birim {f(d['birim_medyan_oran'])} [{f(d['birim_ga95'][0])}, {f(d['birim_ga95'][1])}] n {d['birim_n']}; "
                          f"birim-aritmetik {f((d['birim_aritmetik'] or {}).get('medyan_oran'))} {[round(x, 2) for x in (d['birim_aritmetik'] or {}).get('ga95', [])]}; "
                          f"başarısız {d['n_basarisiz']}/{d['n_pencere']}; ATE medyan {f(d['medyan_ate'])} m (temiz {f(t['temiz']['medyan'])})")
    print("\nHI1 (son yörünge):", res["HI1"]["karar"])
    for kim in ("si", "stereo"):
        for dr, d in res["HI1"][kim].items():
            print(f"  {kim:6s} D/r0 {dr:>2}: başarısız {d['n_basarisiz']}/{d['n_pencere']} ({100 * d['basarisizlik']:.0f}%); oran {f(d['medyan_oran'])} "
                  f"[{f(d['ga95'][0])}, {f(d['ga95'][1])}] n_çift {d['n_cift']}; birim {f(d['birim_medyan_oran'])}; ATE medyan {f(d['medyan_ate'])} m")
    print("  SI temizde başarısız pencereler:", res["HI1"]["si_temizde_basarisiz_pencereler"])
    for dr, e in es.items():
        print(f"  eşli D/r0 {dr:>2}: SI {e['si']}/{e['n_ortak']}, stereo {e['stereo']}/{e['n_ortak']}; yalnız SI {e['yalniz_si']}, yalnız stereo "
              f"{e['yalniz_stereo']}, ikisi {e['ikisi_de']}; yapısal hariç SI {e['yapisal_haric_si']}/{e['yapisal_haric_n']}, stereo {e['yapisal_haric_stereo']}/{e['yapisal_haric_n']}; "
              f"McNemar p {e['mcnemar_kesin_p (sonradan, betimsel)']}, yapısal hariç p {e['mcnemar_kesin_p_yapisal_haric (sonradan, betimsel)']}")
    print("->", out)


if __name__ == "__main__":
    main()
