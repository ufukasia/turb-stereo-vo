"""
AirSim kolları analizi: H3 (temiz A0 ANEES), HA4 (derinlik dilimine göre ek artık), kollar (betimsel), k kalibrasyonu.
Tanımlar SAPMA_TASLAK_KOLLAR.md / ON_TEST_KOLLAR.md §2'dedir (türbülanslı sonuçlardan ÖNCE yazıldı).

    python analiz_kollar.py hepsi   [--etiket izgara] [--tohumlar 10 11 12 13 14]
    python analiz_kollar.py h3|ha4|kollar [--etiket duman --tohumlar 0]      # betik sınaması
    python analiz_kollar.py kalibre                                           # k_airsim.json yazar

H3  : birim (pencere, tohum); A0, F5, C0=0. ANEES_poz = pencere NEES_poz ortalaması (n örnek), bant χ²(3n)/n %95.
      Karar: medyan ANEES_poz > 10 -> H3 öngörüsü tuttu; medyan ANEES_poz bant içinde (medyan n'nin bandı) -> H3
      çürüdü (EuRoC aşırı güveni GT/kalibrasyon kaynaklı); arada -> "bandın üstünde, ≤ 10" (kararsız) raporlanır.
HA4 : çift (pencere, tohum), F5, C0 = 1e-11 ile C0 = 0. Dilim (esit8) başına Δσ² = σ_rob,türb² − σ_rob,temiz²
      (GT-pozlu artık, artik_dok.py). Dilim bir çiftte geçerli: iki hücrede de n ≥ 200. Dilim başına çiftler üzerinden
      medyan Δσ²; dilim, çiftlerin en az yarısında geçerliyse girer. ρ = Spearman(dilim sırası, medyan Δσ²).
      Karar: ρ > 0.8 -> HA4 öngörüsü tuttu; aksi halde tutmadı.
"""
import argparse
import datetime
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
CIKTI = HERE / "cikti_kollar"
sys.path.insert(0, str(HERE))
CHI3_MED = 2.366                # χ²₃ medyanı (EuRoC k ölçütüyle aynı; 2.37)
HA4_MIN_N = 200
HA4_ESIK = 0.8
H3_ESIK = 10.0


def oku(etiket):
    f = CIKTI / etiket / "kayit.jsonl"
    R = [json.loads(l) for l in open(f)] if f.exists() else []
    return [r for r in R if r.get("status") == "ok"], [r for r in R if r.get("status") != "ok"]


def boot_med(x, B=10000, seed=0):
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    m = np.median(x[rng.integers(0, len(x), (B, len(x)))], 1)
    return [float(v) for v in np.percentile(m, [2.5, 97.5])]


def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); adj = np.empty(m); run = 0.0
    for i, j in enumerate(o):
        run = max(run, (m - i) * p[j]); adj[j] = min(1.0, run)
    return adj


def bant(n, dof=3):
    return [float(stats.chi2.ppf(0.025, dof * n) / n), float(stats.chi2.ppf(0.975, dof * n) / n)]


# --------------------------------------------------------------------------------------------- H3
def h3(R, tohumlar, sekans="AF_F5_h5"):
    U = [r for r in R if r.get("tur") == "kol" and r["airsim"]["kol"] == "A0" and r["airsim"]["sekans"] == sekans
         and float(r["airsim"]["C0"]) == 0.0 and int(r["airsim"]["tohum"]) in tohumlar and r.get("n_valid")]
    if not U:
        return {"n_birim": 0}
    rows = sorted([(r["airsim"]["pencere"], int(r["airsim"]["tohum"]), float(r["n_nees_pos_mean"]),
                    float(r["n_nees_rot_mean"]), int(r["n_n"])) for r in U])
    ap = np.array([x[2] for x in rows]); ar = np.array([x[3] for x in rows]); nn = np.array([x[4] for x in rows])
    med, n_med = float(np.median(ap)), int(np.median(nn))
    b = bant(n_med)
    ust = sum(1 for a, n in zip(ap, nn) if a > bant(int(n))[1])
    ici = sum(1 for a, n in zip(ap, nn) if bant(int(n))[0] <= a <= bant(int(n))[1])
    # Bandın ALTI (aşırı güven yok, eksik güven) tanımda öngörülmemişti; 2026-11-08 BI-1'de (F15 temiz, ANEES 1.09)
    # görüldükten sonra AYRI etiketle eklendi (karar kuralı değişmedi: öngörü "> 10" bu durumda da tutmaz).
    karar = ("TUTTU (medyan ANEES_poz > 10)" if med > H3_ESIK else
             "ÇÜRÜDÜ (medyan ANEES_poz bant içinde)" if b[0] <= med <= b[1] else
             "KARARSIZ (bandın üstünde, ≤ 10)" if med > b[1] else
             "TUTMADI — bandın ALTINDA (eksik güven; tanımda öngörülmemiş durum, koordinatör kararı)")
    return {"n_birim": len(rows), "birimler": [dict(pencere=p, tohum=t, anees_poz=a, anees_don=d, n=n,
                                                    bant=bant(n)) for p, t, a, d, n in rows],
            "anees_poz_medyan": med, "anees_poz_medyan_ga95": boot_med(ap), "anees_don_medyan": float(np.median(ar)),
            "n_medyan": n_med, "bant_medyan_n": b, "bant_ustu": ust, "bant_ici": ici, "karar": karar}


# --------------------------------------------------------------------------------------------- HA4
def _dilim_tablo(r, ad, alan="artik_gt"):
    return [(d["n"], d["sigma_rob"]) for d in r[alan]["dilim"][ad]]


def ha4(R, tohumlar, sekans="AF_F5_h5", c0=1e-11, dilim="esit8", alan="artik_gt"):
    src = [r for r in R if (r.get("tur") == "artik" if alan == "artik_gt" else
                           (r.get("tur") == "kol" and r["airsim"]["kol"] == "A0" and "dilim" in r.get(alan, {})))
           and r["airsim"]["sekans"] == sekans and int(r["airsim"]["tohum"]) in tohumlar]
    idx = {(r["airsim"]["pencere"], int(r["airsim"]["tohum"]), float(r["airsim"]["C0"])): r for r in src}
    ciftler = sorted({(p, t) for (p, t, c) in idx if (p, t, 0.0) in idx and (p, t, float(c0)) in idx})
    if not ciftler:
        return {"n_cift": 0}
    import artik_dok as A
    kenar = A.DILIMLER[dilim]
    S = len(kenar) - 1
    D = np.full((len(ciftler), S), np.nan)
    sig_c = np.full_like(D, np.nan); sig_t = np.full_like(D, np.nan)
    for i, (p, t) in enumerate(ciftler):
        tc = _dilim_tablo(idx[(p, t, 0.0)], dilim, alan); tt = _dilim_tablo(idx[(p, t, float(c0))], dilim, alan)
        for s in range(S):
            (nc, sc), (nt, st) = tc[s], tt[s]
            if nc >= HA4_MIN_N and nt >= HA4_MIN_N:
                D[i, s] = st ** 2 - sc ** 2; sig_c[i, s] = sc; sig_t[i, s] = st
    gecerli = np.sum(np.isfinite(D), 0)
    giren = [s for s in range(S) if gecerli[s] >= math.ceil(len(ciftler) / 2)]
    med = np.array([np.nanmedian(D[:, s]) for s in giren]) if giren else np.array([])
    rho = float(stats.spearmanr(np.arange(len(giren)), med)[0]) if len(giren) >= 3 else float("nan")
    per = []
    for i in range(len(ciftler)):
        ok = [s for s in giren if np.isfinite(D[i, s])]
        if len(ok) >= 3:
            per.append(float(stats.spearmanr(ok, D[i, ok])[0]))
    # model σ_tilt (ana ızgara kaydı, uret_derinlik model dilimleri) — yalnız bilgi
    model = {}
    for (p, t, c), r in idx.items():
        if c == float(c0) and r["airsim"].get("turb_dilim"):
            for dd, v in r["airsim"]["turb_dilim"].items():
                model.setdefault(dd, []).append(v.get("sigma_tot_px"))
    return {"dilim": dilim, "alan": alan, "n_cift": len(ciftler), "ciftler": ciftler,
            "kenar": kenar, "giren_dilim": [f"{kenar[s]:g}-{kenar[s + 1]:g}" for s in giren],
            "gecerli_cift": [int(gecerli[s]) for s in range(S)],
            "medyan_dsig2": [float(x) for x in med],
            "medyan_ek_artik_px": [float(math.sqrt(max(x, 0.0))) for x in med],
            "medyan_sigma_temiz": [float(np.nanmedian(sig_c[:, s])) for s in giren],
            "medyan_sigma_turb": [float(np.nanmedian(sig_t[:, s])) for s in giren],
            "rho": rho, "cift_rho": per, "cift_rho_medyan": float(np.median(per)) if per else float("nan"),
            "cift_rho_gt08": int(sum(1 for x in per if x > HA4_ESIK)),
            "model_sigma_tot_px": {k: float(np.median([x for x in v if x is not None])) for k, v in model.items()},
            "karar": ("TUTTU (ρ > 0.8)" if rho > HA4_ESIK else "TUTMADI (ρ ≤ 0.8)") if np.isfinite(rho) else "YETERSİZ"}


# ------------------------------------------------------------------------------------- kollar (betimsel)
def kollar(R, tohumlar, sekans="AF_F5_h5"):
    K = [r for r in R if r.get("tur") == "kol" and r["airsim"]["sekans"] == sekans
         and int(r["airsim"]["tohum"]) in tohumlar and r.get("n_valid")]
    idx = {(r["airsim"]["pencere"], int(r["airsim"]["tohum"]), float(r["airsim"]["C0"]), r["airsim"]["kol"]): r
           for r in K}
    out, pl = [], []
    c0s = sorted({k[2] for k in idx}); arms = sorted({k[3] for k in idx} - {"A0"})
    for c in c0s:
        for kol in arms:
            u = [(p, t) for (p, t, cc, kk) in idx if cc == c and kk == kol and (p, t, c, "A0") in idx]
            if not u:
                continue
            ln = np.array([math.log(idx[(p, t, c, kol)]["n_nees_pos_mean"] / idx[(p, t, c, "A0")]["n_nees_pos_mean"])
                           for p, t in u])
            lr = np.array([math.log(idx[(p, t, c, kol)]["a_rpe1s"] / idx[(p, t, c, "A0")]["a_rpe1s"]) for p, t in u])
            pv = float(stats.wilcoxon(ln).pvalue) if len(ln) >= 5 and np.any(ln != 0) else float("nan")
            out.append(dict(C0=c, kol=kol, n=len(u), anees_oran_medyan=float(np.exp(np.median(ln))),
                            anees_oran_ga95=[float(np.exp(x)) for x in boot_med(ln)] if len(ln) > 1 else None,
                            anees_kol_medyan=float(np.median([idx[(p, t, c, kol)]["n_nees_pos_mean"] for p, t in u])),
                            rpe1s_oran_medyan=float(np.exp(np.median(lr))), p_wilcoxon=pv))
            pl.append(pv)
    ok = [i for i, p in enumerate(pl) if np.isfinite(p)]
    if ok:
        adj = holm([pl[i] for i in ok])
        for j, i in enumerate(ok):
            out[i]["p_holm"] = float(adj[j])
    # A3 vs AG (H4 benzeri, betimsel)
    a3ag = []
    for c in c0s:
        u = [(p, t) for (p, t, cc, kk) in idx if cc == c and kk == "A3" and (p, t, c, "AG") in idx]
        if u:
            ln = [math.log(idx[(p, t, c, "A3")]["n_nees_pos_mean"] / idx[(p, t, c, "AG")]["n_nees_pos_mean"]) for p, t in u]
            a3ag.append(dict(C0=c, n=len(u), anees_A3_AG_oran_medyan=float(np.exp(np.median(ln)))))
    return {"kol_A0": out, "A3_AG": a3ag}


# --------------------------------------------------------------------------------------------- k
def kalibre():
    R, F = oku("kalibre")
    tab = defaultdict(list)
    a0 = {}
    for r in R:
        a = r["airsim"]
        if a["kol"] == "A0C" and r.get("n_valid"):
            tab[float(a["k"])].append((a["turb_sekans"], float(r["n_nees_pos_median"]), float(r["n_nees_rot_median"])))
        if a["kol"] == "A0":
            a0[a["turb_sekans"]] = r
    pencereler = sorted({x[0] for v in tab.values() for x in v})
    satir = []
    for k in sorted(tab):
        v = tab[k]
        satir.append(dict(k=k, n_pencere=len(v), nees_poz_medyan_pencere_medyani=float(np.median([x[1] for x in v])),
                          nees_don_medyan_pencere_medyani=float(np.median([x[2] for x in v]))))
    tam = [s for s in satir if s["n_pencere"] == len(pencereler)]
    if not tam:
        print("kalibrasyon kaydı yok/eksik"); return None
    sec = min(tam, key=lambda s: abs(s["nees_poz_medyan_pencere_medyani"] - CHI3_MED))
    # plasebo: k=1 A0C == A0 (bayt)
    plasebo = []
    for r in R:
        a = r["airsim"]
        if a["kol"] == "A0C" and float(a["k"]) == 1.0 and a["turb_sekans"] in a0:
            plasebo.append(r.get("traj_sha1") == a0[a["turb_sekans"]].get("traj_sha1"))
    a0n = [float(r["n_nees_pos_median"]) for r in a0.values() if r.get("n_valid")]
    out = {"tarih": f"{datetime.datetime.now():%Y-%m-%d %H:%M}", "k": sec["k"], "olcut": "medyan(pencere medyan NEES_poz) χ²₃ medyanına (2.37) en yakın",
           "tablo": satir, "pencereler": pencereler, "A0_pencere_medyan_nees_poz_medyani": float(np.median(a0n)) if a0n else None,
           "plasebo_k1_A0_bayt_ayni": f"{sum(plasebo)}/{len(plasebo)}", "basarisiz_kayit": len(F),
           "sinirda": sec["k"] in (min(tab), max(tab)),
           "veri": "F5 temiz (C0=0), pilot tohumları 0-2, W1+W2, σ=1 DN, RTX 3090; türbülanslı veri kullanılmadı"}
    for s in satir:
        print(f"k={s['k']:>4g}  n={s['n_pencere']}  NEES_poz {s['nees_poz_medyan_pencere_medyani']:.3f}  "
              f"NEES_don {s['nees_don_medyan_pencere_medyani']:.3f}")
    print(f"SEÇİLEN k = {sec['k']:g}  (sınırda: {out['sinirda']}; plasebo {out['plasebo_k1_A0_bayt_ayni']})")
    if out["sinirda"]:
        print("!! k ızgara sınırında — ızgara başlatılmamalı, koordinatöre danış")
    json.dump(out, open(CIKTI / "k_airsim.json", "w"), indent=1, ensure_ascii=False)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kip", choices=["hepsi", "h3", "ha4", "kollar", "kalibre"])
    ap.add_argument("--etiket", default="izgara")
    ap.add_argument("--tohumlar", nargs="+", type=int, default=[10, 11, 12, 13, 14])
    ap.add_argument("--cikti", default=None)
    a = ap.parse_args()
    if a.kip == "kalibre":
        kalibre(); return
    R, F = oku(a.etiket)
    th = set(a.tohumlar)
    out = {"etiket": a.etiket, "tohumlar": sorted(th), "n_kayit": len(R), "basarisiz": len(F),
           "code_fp": sorted({r.get("code_fp") for r in R if r.get("tur") == "kol"}),
           "a0_bit_ayni": [r.get("a0_bit_ayni") for r in R if r.get("tur") == "kol" and r["airsim"]["kol"] == "A0"]}
    if a.kip in ("hepsi", "h3"):
        out["H3"] = h3(R, th)
    if a.kip in ("hepsi", "ha4"):
        out["HA4"] = ha4(R, th)
        out["HA4_ikincil_model_dilim"] = ha4(R, th, dilim="model")
        out["HA4_ikincil_kestirim_poz"] = ha4(R, th, alan="artik_est")
    if a.kip in ("hepsi", "kollar"):
        out["kollar"] = kollar(R, th)
    s = json.dumps(out, indent=1, ensure_ascii=False, default=str)
    print(s)
    if a.cikti:
        open(a.cikti, "w").write(s)


if __name__ == "__main__":
    main()
