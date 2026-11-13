#!/usr/bin/env python
"""A3S ayırma kolu analizi (ON_KAYIT_FAZ2.md §5 Sapma 4; işlemsel tanım evidence/A3S/ON_TEST.md §4).

    python a3s_analiz.py              # evidence/A3S/records.jsonl -> evidence/A3S/analiz.json + özet
    python a3s_analiz.py ontest       # evidence/A3S/on_test.jsonl ön-test tablosu
    python a3s_analiz.py oz_test      # hesap denetimi: ana ızgaranın A3 kayıtları "A3S" diye verilir

Birim, metrik ve testler ana analizle aynıdır (turb/faz2_analiz.py'den içe aktarılır): birim (sekans, tohum);
ATE pencerelerde aritmetik, ANEES geometrik ortalama; medyan oran, bootstrap %95 GA (10 000, tohum 0), iki yönlü
Wilcoxon; Holm ailesi = beş D/r0 seviyesi (karşılaştırma ve metrik başına). Eşdeğerlik: t-TOST, sınır ±ln 1.15,
alfa = 0.05 (birim log-oranlarının ortalamasının %90 t-GA'sı sınır içinde); bootstrap %90 GA da raporlanır.
"""
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
from turb.faz2_analiz import oku, eslesik, oran_testi, holm  # noqa: E402

EV = os.path.join(ROOT, "evidence", "A3S")
ANA = os.path.join(ROOT, "evidence", "FAZ2", "records.jsonl")
FP = "08e28b83407c"
KOL = "A0C+A3"
SINIR = math.log(1.15)
DUZ = (8.0, 12.0, 16.0, 20.0, 24.0)       # doğrulayıcı seviyeler
ON3 = (8.0, 12.0, 16.0)                   # HS1 / HS3 karar seviyeleri
EK = (0.0, 4.0)                           # betimsel (koşulduysa)
FE_SURE = "total_execution_time_seconds"
BIT = ("m_ate_se3", "n_nees_pos_mean", "n_nees_rot_mean", "m_rpe_trans")


def tost(a, b, seed=0, B=10000):
    lr = np.log(np.asarray(a, float) / np.asarray(b, float))
    n = len(lr)
    if n < 2:
        return dict(esdeger=False)
    m, se = float(lr.mean()), float(lr.std(ddof=1) / math.sqrt(n))
    rng = np.random.default_rng(seed)
    bs = lr[rng.integers(0, n, (B, n))].mean(axis=1)
    b5, b95 = float(np.percentile(bs, 5)), float(np.percentile(bs, 95))
    if se == 0.0:
        ga, p = (m, m), (0.0 if abs(m) < SINIR else 1.0)
    else:
        tq = float(stats.t.ppf(0.95, n - 1))
        ga = (m - tq * se, m + tq * se)
        p = float(max(stats.t.sf((m + SINIR) / se, n - 1), stats.t.cdf((m - SINIR) / se, n - 1)))
    return dict(ort_oran=math.exp(m), ga90_t=[math.exp(ga[0]), math.exp(ga[1])], p_tost=p,
                esdeger=bool(ga[0] > -SINIR and ga[1] < SINIR),
                ga90_boot=[math.exp(b5), math.exp(b95)], esdeger_boot=bool(b5 > -SINIR and b95 < SINIR))


def a3s_oku(path):
    """-> ({(etiket, kol, dr, sekans, tohum, pencere): kayıt} (son kayıt geçerli), hücre hatası listesi)"""
    R, hata = {}, []
    with open(path) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            c = r.get("cell") or {}
            if "arm" not in r:
                hata.append((r.get("key"), r.get("error")))
                continue
            R[(c.get("tag", ""), r["arm"], float(c["D_over_r0"]), c["sequence"], int(c["seed"]), int(r["window"]))] = r
    return R, hata


def fe_ayni(a, b):
    if not a or not b:
        return False
    return ({k: v for k, v in a.items() if k != FE_SURE} == {k: v for k, v in b.items() if k != FE_SURE})


def bit_ayni(r, m):
    return m is not None and all(r.get(k) == m.get(k) for k in BIT)


def hucre_denetimi(R, D, etiket="a3s"):
    """Sapma 4 denetim (a), (b) + kol kimliği. -> (geçerli hücre kümesi {(dr, sekans, tohum)}, rapor)"""
    H = defaultdict(list)
    for (tag, kol, dr, s, t, w), r in R.items():
        if tag == etiket:
            H[(dr, s, t)].append(r)
    gecerli, rapor = set(), dict(n_hucre=len(H), fe_farkli=[], gpu_farkli=[], fp_farkli=[], kol_kimligi=[],
                                 a0_bit=[0, 0], a0_farkli=[], eksik_pencere=[])
    for (dr, s, t), recs in sorted(H.items()):
        ad = f"{s}/dr{dr:g}/s{t}"
        ok = True
        ref = None                                   # ana ızgaranın aynı hücredeki herhangi bir kaydı (ön-uç ortak)
        for kol in ("A0", "A3", "AG", "A1"):
            for r in recs:
                ref = D.get((kol, dr, s, t, int(r["window"])))
                if ref is not None:
                    break
            if ref is not None:
                break
        if ref is None or not fe_ayni(recs[0].get("frontend_stats"), ref.get("frontend_stats")):
            rapor["fe_farkli"].append(ad)
            ok = False
        if "3090" not in str(recs[0].get("fe_gpu")):
            rapor["gpu_farkli"].append(ad)
            ok = False
        if any(r.get("code_fp") != FP for r in recs):
            rapor["fp_farkli"].append(ad)
            ok = False
        a3s = [r for r in recs if r["arm"] == KOL]
        for r in a3s:
            st = r.get("arm_stats") or {}
            if not (st.get("cal_factor") == 1.0 and st.get("n_patch_factors") == 0
                    and st.get("dimm_sigma_scale") == 1.0):
                rapor["kol_kimligi"].append(f"{ad}/w{r['window']}")
                ok = False
        if len(a3s) != len(recs[0]["cell"]["windows"]):
            rapor["eksik_pencere"].append(ad)
            ok = False
        for r in recs:
            if r["arm"] == "A0":
                rapor["a0_bit"][1] += 1
                if bit_ayni(r, D.get(("A0", dr, s, t, int(r["window"])))):
                    rapor["a0_bit"][0] += 1
                else:
                    rapor["a0_farkli"].append(f"{ad}/w{r['window']}")
                    ok = False
        if ok:
            gecerli.add((dr, s, t))
    rapor["n_gecerli"] = len(gecerli)
    return gecerli, rapor


def karsilastir(D, ref, metrik, duzeyler=DUZ):
    T = {}
    for dr in duzeyler:
        a, b, U, dus = eslesik(D, "A3S", ref, dr, metrik)
        if len(a) < 2:
            continue
        T[f"{dr:g}"] = dict(oran_testi(a, b), **tost(a, b), dusen_birim=dus)
    ph = holm({k: v["p"] for k, v in T.items()})
    for k in T:
        T[k]["p_holm"] = ph[k]
    return T


def betimsel(D, duzeyler):
    out = {}
    for metrik in ("ANEES_poz", "ATE", "ANEES_don"):
        for kol in ("A3S", "A3", "AG", "A1", KOL):
            for dr in duzeyler:
                a, b, U, dus = eslesik(D, kol, "A0", dr, metrik)
                if len(a):
                    out[f"{metrik}:{kol}/A0@{dr:g}"] = dict(n=len(a), medyan_oran=float(np.exp(np.median(np.log(a / b)))))
    return out


def sigma_ve_sure(D, duzeyler):
    """Hücre başına sigma_turb_used_mean ve t_backend (kol koşumu başına bir değer): A3S / A3."""
    out = {}
    for dr in duzeyler:
        s_or, t_a, t_b, bas = [], [], [], {"A3S": [0, 0], "A3": [0, 0]}
        gor = set()
        for (kol, d, s, t, w), r in D.items():
            if d != dr or kol not in ("A3S", "A3"):
                continue
            bas[kol][1] += 1
            bas[kol][0] += int(not (r.get("status") == "ok" and r.get("m_valid") and r.get("n_valid")
                                    and not r.get("m_diverged")))
            if kol != "A3S" or (s, t) in gor:
                continue
            r3 = next((D[("A3", d, s, t, w2)] for w2 in (w,) if ("A3", d, s, t, w2) in D), None)
            if r3 is None:
                continue
            gor.add((s, t))
            x, y = (r.get("arm_stats") or {}).get("sigma_turb_used_mean"), (r3.get("arm_stats") or {}).get("sigma_turb_used_mean")
            if x and y and math.isfinite(x) and math.isfinite(y):
                s_or.append(x / y)
            if r.get("t_backend") and r3.get("t_backend"):
                t_a.append(float(r["t_backend"]))
                t_b.append(float(r3["t_backend"]))
        if not s_or:
            continue
        s_or = np.array(s_or)
        out[f"{dr:g}"] = dict(n=len(s_or), sigma_oran_medyan=float(np.median(s_or)),
                              sigma_oran_min=float(s_or.min()), sigma_oran_max=float(s_or.max()),
                              sigma_5_icinde=int(np.sum(np.abs(s_or - 1) <= 0.05)),
                              t_be_medyan_A3S=float(np.median(t_a)), t_be_medyan_A3_ana=float(np.median(t_b)),
                              t_be_oran=float(np.median(t_a) / np.median(t_b)),
                              gecersiz_pencere_A3S=bas["A3S"], gecersiz_pencere_A3=bas["A3"])
    return out


def kararlar(S):
    hs1, hs2, hs3 = S["HS1_ANEES_A3S_A3"], S["HS2_ATE_A3S_A3"], S["HS3_ANEES_A3S_AG"]
    k3 = [f"{d:g}" for d in ON3]
    k5 = [f"{d:g}" for d in DUZ]
    tam = all(k in hs1 for k in k5)
    K = dict(veri_tam=tam)
    K["HS1"] = {k: bool(hs1[k]["esdeger"]) for k in k5 if k in hs1}
    K["HS1_karar"] = ("TUTTU" if all(K["HS1"].get(k) for k in k3) else "TUTMADI") if all(k in hs1 for k in k3) else "veri eksik"
    K["HS2"] = {k: bool(hs2[k]["esdeger"]) for k in k5 if k in hs2}
    K["HS2_karar"] = ("TUTTU" if all(K["HS2"].get(k) for k in k5) else "TUTMADI") if tam else "veri eksik"
    K["HS3"] = {k: dict(gecer=bool(hs3[k]["medyan_oran"] < 1 and hs3[k]["p_holm"] < 0.05),
                        medyan_beklenen_aralikta=bool(0.78 <= hs3[k]["medyan_oran"] <= 0.92)) for k in k5 if k in hs3}
    K["HS3_karar"] = ("TUTTU" if all(K["HS3"][k]["gecer"] for k in k3) else "TUTMADI") if all(k in hs3 for k in k3) else "veri eksik"
    buyuk = [k for k in k3 if k in hs1 and hs1[k]["medyan_oran"] > 1.15 and hs1[k]["p_holm"] < 0.05]
    kucuk = [k for k in k3 if k in hs1 and hs1[k]["medyan_oran"] < 0.87 and hs1[k]["p_holm"] < 0.05]
    if K["HS1_karar"] == "veri eksik":
        y = "veri eksik"
    elif K["HS1_karar"] == "TUTTU":
        y = ("mekanizma: A3'ün AG'ye üstünlüğü yama korelasyon modelinden değil, DIMM büyüklüğünden gelir; "
             "yarı-statik rejim için önerilen kol A3S")
    elif len(buyuk) >= 2:
        y = f"korelasyon yapısı tutarlılığa katkı veriyor (D/r0 {', '.join(buyuk)}); A3 önerilen kol olarak kalır"
    elif len(kucuk) >= 2:
        y = f"skaler kol daha tutarlı (D/r0 {', '.join(kucuk)})"
    else:
        y = "kararsız: seviye seviye GA'larıyla raporlanır"
    K["yorum"] = y
    return K


def calistir(D, duzeyler=DUZ):
    S = dict(HS1_ANEES_A3S_A3=karsilastir(D, "A3", "ANEES_poz", duzeyler),
             HS2_ATE_A3S_A3=karsilastir(D, "A3", "ATE", duzeyler),
             HS3_ANEES_A3S_AG=karsilastir(D, "AG", "ANEES_poz", duzeyler),
             ikincil_ANEES_don_A3S_A3=karsilastir(D, "A3", "ANEES_don", duzeyler),
             ikincil_ATE_A3S_AG=karsilastir(D, "AG", "ATE", duzeyler))
    S["kararlar"] = kararlar(S)
    S["betimsel_A0a_gore"] = betimsel(D, duzeyler)
    S["sigma_ve_sure"] = sigma_ve_sure(D, duzeyler)
    return S


def yazdir(S):
    for ad in ("HS1_ANEES_A3S_A3", "HS2_ATE_A3S_A3", "HS3_ANEES_A3S_AG", "ikincil_ANEES_don_A3S_A3",
               "ikincil_ATE_A3S_AG"):
        print(f"\n== {ad} ==")
        print("  D/r0   n  medyan  boot95            p_W      p_Holm   <1    ort.oran  GA90 t (TOST)     p_TOST   eşd  GA90 boot        eşd_b")
        for k, v in S[ad].items():
            print(f"  {k:>4} {v['n']:3d}  {v['medyan_oran']:6.3f}  [{v['ga95'][0]:.3f}, {v['ga95'][1]:.3f}]  "
                  f"{v['p']:8.2e} {v['p_holm']:8.2e} {v['yon_azalma']:2d}/{v['n']:<2d}  {v['ort_oran']:6.3f}   "
                  f"[{v['ga90_t'][0]:.3f}, {v['ga90_t'][1]:.3f}]  {v['p_tost']:8.2e}  {'E' if v['esdeger'] else '-'}   "
                  f"[{v['ga90_boot'][0]:.3f}, {v['ga90_boot'][1]:.3f}]  {'E' if v['esdeger_boot'] else '-'}")
    print("\n== kararlar ==")
    for k, v in S["kararlar"].items():
        print(f"  {k}: {v}")
    print("\n== betimsel: medyan oran (kol / A0) ==")
    B = S["betimsel_A0a_gore"]
    for metrik in ("ANEES_poz", "ATE", "ANEES_don"):
        dz = sorted({float(k.split("@")[1]) for k in B if k.startswith(metrik + ":")})
        print(f"  {metrik:10s} " + "".join(f"{d:>8g}" for d in dz))
        for kol in ("A3S", "A3", "AG", "A1", KOL):
            print(f"    {kol + '/A0':9s} " + "".join(
                f"{B[f'{metrik}:{kol}/A0@{d:g}']['medyan_oran']:8.3f}" if f"{metrik}:{kol}/A0@{d:g}" in B else "       -"
                for d in dz))
    print("\n== sigma (A3S / A3) ve arka-uç süresi ==")
    print("  D/r0   n  sigma oranı medyan [min, max]  ±%5 içinde   t_be A3S / A3(ana) [s]   oran   geçersiz pencere A3S / A3")
    for k, v in S["sigma_ve_sure"].items():
        print(f"  {k:>4} {v['n']:3d}  {v['sigma_oran_medyan']:.3f} [{v['sigma_oran_min']:.3f}, {v['sigma_oran_max']:.3f}]"
              f"       {v['sigma_5_icinde']:2d}/{v['n']:<2d}        {v['t_be_medyan_A3S']:6.0f} / {v['t_be_medyan_A3_ana']:6.0f}"
              f"        {v['t_be_oran']:.2f}   {v['gecersiz_pencere_A3S']} / {v['gecersiz_pencere_A3']}")


def ana_komut():
    D = oku(ANA)
    R, hata = a3s_oku(os.path.join(EV, "records.jsonl"))
    gecerli, rapor = hucre_denetimi(R, D)
    print(f"[a3s] {len(R)} kayıt, {rapor['n_hucre']} hücre, geçerli {rapor['n_gecerli']}; hücre hatası {len(hata)}")
    for k in ("fe_farkli", "gpu_farkli", "fp_farkli", "kol_kimligi", "eksik_pencere", "a0_farkli"):
        if rapor[k]:
            print(f"  !! {k}: {len(rapor[k])}  {rapor[k][:6]}")
    print(f"  denetim (b): A0 bit-aynı {rapor['a0_bit'][0]}/{rapor['a0_bit'][1]}")
    for (tag, kol, dr, s, t, w), r in R.items():
        if tag == "a3s" and kol == KOL and (dr, s, t) in gecerli:
            D[("A3S", dr, s, t, w)] = r
    duz = tuple(d for d in DUZ + EK if any(k[0] == "A3S" and k[1] == d for k in D))
    S = calistir(D, tuple(d for d in duz if d in DUZ))
    ek = tuple(d for d in duz if d in EK)
    if ek:
        S["ek_betimsel"] = dict(ANEES_A3S_A3=karsilastir(D, "A3", "ANEES_poz", ek), ATE_A3S_A3=karsilastir(D, "A3", "ATE", ek),
                                A0a_gore=betimsel(D, ek), sigma_ve_sure=sigma_ve_sure(D, ek))
    S["denetim"] = rapor
    S["hucre_hatasi"] = hata
    S["kaynak"] = dict(a3s="evidence/A3S/records.jsonl", ana="evidence/FAZ2/records.jsonl", n_kayit=len(R))
    yazdir(S)
    if ek:
        print("\n== ek betimsel seviyeler (D/r0 0, 4) ==")
        print(json.dumps(S["ek_betimsel"], indent=1)[:3000])
    with open(os.path.join(EV, "analiz.json"), "w") as f:
        json.dump(S, f, indent=1)
    print(f"\n[a3s] yazıldı: {os.path.join(EV, 'analiz.json')}")


def oz_test():
    """Ana ızgaranın A3 kayıtları 'A3S' diye verilir: A3S/A3 tam 1 ve eşdeğer; A3S/AG = analiz.json H4 ii."""
    D = oku(ANA)
    for k in [k for k in D if k[0] == "A3"]:
        D[("A3S",) + k[1:]] = D[k]
    S = calistir(D)
    yazdir(S)
    H = json.load(open(os.path.join(ROOT, "evidence", "FAZ2", "analiz.json")))["H4"]["ii"]
    ok = True
    for k, v in S["HS1_ANEES_A3S_A3"].items():
        ok &= v["medyan_oran"] == 1.0 and v["esdeger"] and v["n"] == 45
    for k, v in S["HS3_ANEES_A3S_AG"].items():
        h = H[k]
        ok &= (v["medyan_oran"] == h["medyan_oran"] and v["ga95"] == h["ga95"] and v["p"] == h["p"]
               and v["p_holm"] == h["p_holm"] and v["n"] == h["n"])
    # sentetik TOST: n = 45, sd = 0.20; ortalama 0 -> eşdeğer, ortalama ln 1.15 -> eşdeğer değil
    rng = np.random.default_rng(1)
    z = rng.normal(0, 0.20, 45)
    z = (z - z.mean()) / z.std(ddof=1) * 0.20
    t0, t1 = tost(np.exp(z), np.ones(45)), tost(np.exp(z + SINIR), np.ones(45))
    se, tq = 0.20 / math.sqrt(45), stats.t.ppf(0.95, 44)
    ok &= t0["esdeger"] and not t1["esdeger"] and abs(t1["p_tost"] - 0.5) < 1e-9
    ok &= abs(t0["ga90_t"][1] - math.exp(tq * se)) < 1e-9 and abs(t0["ga90_t"][0] - math.exp(-tq * se)) < 1e-9
    print(f"\nsentetik TOST: ortalama 0 -> GA90 [{t0['ga90_t'][0]:.4f}, {t0['ga90_t'][1]:.4f}] (beklenen ±{math.exp(tq * se):.4f}), "
          f"p {t0['p_tost']:.2e}, eşdeğer {t0['esdeger']}; ortalama ln 1.15 -> p {t1['p_tost']:.3f}, eşdeğer {t1['esdeger']}")
    print("\nÖZ TEST:", "GEÇTİ" if ok else "KALDI")
    sys.exit(0 if ok else 1)


def ontest():
    D = oku(ANA)
    R, hata = a3s_oku(os.path.join(EV, "on_test.jsonl"))
    print(f"[a3s ön-test] {len(R)} kayıt; hücre hatası {len(hata)} {hata[:3]}")
    print("\n== denetim hücresi (k = 12) ve A0: ana ızgarayla bit-aynılık ==")
    n = nb = 0
    for (tag, kol, dr, s, t, w), r in sorted(R.items()):
        if kol == "A0" or tag == "a3s_ontest_k12":
            m = D.get((kol, dr, s, t, w))
            n += 1
            b = bit_ayni(r, m)
            nb += b
            print(f"  {tag:15s} {kol:7s} dr{dr:g} w{w:4d}  ATE {r.get('m_ate_se3')} / {None if m is None else m.get('m_ate_se3')}  "
                  f"ANEES {r.get('n_nees_pos_mean'):.4f} / {m.get('n_nees_pos_mean'):.4f}  {'bit-aynı' if b else 'FARKLI'}")
    print(f"  toplam bit-aynı {nb}/{n}")
    print("\n== ön-uç istatistikleri (süre alanı hariç) ana ızgarayla ==")
    gor = set()
    for (tag, kol, dr, s, t, w), r in sorted(R.items()):
        if (tag, dr) in gor:
            continue
        gor.add((tag, dr))
        m = D.get(("A0", dr, s, t, w))
        fs = r.get("frontend_stats") or {}
        print(f"  {tag:15s} dr{dr:g}: {'aynı' if fe_ayni(fs, m.get('frontend_stats')) else 'FARKLI'}  "
              f"stereo {fs.get('avg_stereo_matches'):.3f} üçgen {fs.get('avg_current_stereo_triangulated_points'):.3f} "
              f"zamansal {fs.get('avg_temporal_matches'):.3f} aktif {fs.get('avg_active_depth_tracks'):.3f}  "
              f"t_fe {r.get('t_frontend'):.0f} s (ana {m.get('t_frontend'):.0f} s)  gpu {r.get('fe_gpu')}  fp {r.get('code_fp')}")
    print("\n== A3S kayıtları (k = 1) ==")
    print("  dr   w    ATE      ANEESp   ANEESr   n_n  durum  cal  n_yama  sigma_hat  sigma_used  /A3(ana)  t_be   | ana A3: ATE ANEESp | ana A0: ATE ANEESp")
    G = defaultdict(lambda: defaultdict(list))
    for (tag, kol, dr, s, t, w), r in sorted(R.items()):
        if tag != "a3s_ontest" or kol != KOL:
            continue
        st = r.get("arm_stats") or {}
        m3, m0 = D[("A3", dr, s, t, w)], D[("A0", dr, s, t, w)]
        su3 = m3["arm_stats"]["sigma_turb_used_mean"]
        print(f"  {dr:<4g} {w:4d} {r['m_ate_se3']:.5f} {r['n_nees_pos_mean']:8.2f} {r['n_nees_rot_mean']:8.2f} {r['n_n']:4d}  "
              f"{r['status']}/{int(bool(r['m_valid']))}{int(bool(r['n_valid']))}  {st.get('cal_factor')}  {st.get('n_patch_factors')}      "
              f"{st.get('sigma_hat_mean'):.4f}    {st.get('sigma_turb_used_mean'):.4f}     {st.get('sigma_turb_used_mean') / su3:.3f}   "
              f"{r['t_backend']:5.0f}  | {m3['m_ate_se3']:.5f} {m3['n_nees_pos_mean']:8.2f} | {m0['m_ate_se3']:.5f} {m0['n_nees_pos_mean']:8.2f}")
        G[dr]["A3S"].append(r)
        for kk in ("A3", "A0", "A1", "AG"):
            G[dr][kk].append(D[(kk, dr, s, t, w)])
    print("\n== sekans düzeyi (3 pencere): ANEES geometrik, ATE aritmetik ortalama ==")
    for dr, g in sorted(G.items()):
        ge = lambda L: float(np.exp(np.mean([math.log(x["n_nees_pos_mean"]) for x in L])))      # noqa: E731
        ar = lambda L: float(np.mean([x["m_ate_se3"] for x in L]))                               # noqa: E731
        print(f"  dr{dr:g}: ANEES " + "  ".join(f"{k} {ge(v):.1f}" for k, v in g.items())
              + f"  | A3S/A3 {ge(g['A3S']) / ge(g['A3']):.3f}  A3S/A0 {ge(g['A3S']) / ge(g['A0']):.3f}")
        print(f"         ATE   " + "  ".join(f"{k} {ar(v):.4f}" for k, v in g.items())
              + f"  | A3S/A3 {ar(g['A3S']) / ar(g['A3']):.3f}  A3S/A0 {ar(g['A3S']) / ar(g['A0']):.3f}")
    alan = ("m_ate_se3", "m_rpe_trans", "n_nees_pos_mean", "n_nees_rot_mean", "n_n", "t_backend", "t_frontend",
            "frontend_stats", "fe_gpu", "code_fp", "status", "m_valid", "n_valid")
    bos = [(k[0], k[1], k[2], k[5], a) for k, r in R.items() for a in alan if r.get(a) is None]
    print(f"\n== alan doluluğu: {len(R)} kayıt, boş alan {len(bos)} {bos[:5]}")


if __name__ == "__main__":
    komut = sys.argv[1] if len(sys.argv) > 1 else ""
    if komut == "ontest":
        ontest()
    elif komut == "oz_test":
        oz_test()
    else:
        ana_komut()
