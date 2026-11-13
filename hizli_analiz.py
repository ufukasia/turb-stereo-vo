#!/usr/bin/env python
"""Hızlı rejim doğrulayıcı ızgarası analizi (ON_KAYIT_FAZ2.md §5 Sapma 5).

Bu betik ızgaranın hiçbir sonucu görülmeden yazıldı; istatistik ana ızgaranınkiyle aynıdır (turb/faz2_analiz.py):
birim (sekans, tohum); birim değeri = iki tarafta da geçerli pencerelerin ortalaması (ATE aritmetik, ANEES
geometrik); etki = eşli medyan oran, birim bootstrap 10 000 ile %95 GA; iki yönlü Wilcoxon işaretli sıra; Holm.

    python hizli_analiz.py                 # evidence/HIZLI/records.jsonl + ana ızgara -> evidence/HIZLI/analiz.json
    python hizli_analiz.py --butunluk      # yalnız bütünlük sayımları (koşum sürerken; test tablosu basmaz)
    python hizli_analiz.py --oz_test       # sentetik veri + ana ızgara sayılarını yeniden üretme denetimi

Rejim R0 = ana ızgara (tau0 = dt; evidence/FAZ2/records.jsonl), hızlı rejim = evidence/HIZLI/records.jsonl
(tau0/dt = 0.314 (D / (D/r0)) / (v dt)). Eşleme (sekans, tohum, D/r0, pencere) üzerinden; iki rejimde ön-uç koşumu
farklıdır, kollar rejim içinde aynı dökümü paylaşır.

Aileler (Holm aile içinde, metrik başına):
  F1  rejim içi kol / A0        : kollar {A1, A3, AG, A0C, A0C+A3} × düzeyler (ATE, ANEES_poz, ANEES_don)
  F2  A3 / AG (ANEES_poz)       : düzeyler
  F3  rejim etkisi hızlı / R0   : kollar {A0, A1, A3, AG, A0C, A0C+A3} × düzeyler (ATE, ANEES_poz)
Betimsel: A0 hızlı@düzey / A0 R0@temiz ATE oranı ve kırılma ölçütü (ana ızgaradaki tanım), başarısız pencere oranı,
kol × düzey birim medyanları, ANEES'in adlandırılmış banda göre konumu.
Kararlar (HZ1–HZ5) eşikleri ON_KAYIT Sapma 5'te yazılı olan sabitlerle verilir (aşağıda ESIK).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
from turb.faz2_analiz import (METRIK, _ort, ate, anees, bant, eslesik, holm, kirilma_tablosu,  # noqa: E402
                              oran_testi)

EV = os.path.join(ROOT, "evidence", "HIZLI")
ANA = os.path.join(ROOT, "evidence", "FAZ2", "records.jsonl")
ANA_ANALIZ = os.path.join(ROOT, "evidence", "FAZ2", "analiz.json")
FP = "08e28b83407c"
KOLLAR = ("A1", "A3", "AG", "A0C", "A0C+A3")
HEPSI = ("A0",) + KOLLAR
D_AP, DT = 0.10, 0.05
REJIM_HIZ = {"R5": 5.0, "R1": 1.0}
# Sapma 5'te yazılı eşikler (ana ızgaranın H1, H4, H5 eşikleriyle aynı)
ESIK = dict(bozulma_oran=1.5,        # HZ1: A0 ATE hızlı@8 / temiz medyan oranı > 1.5 -> bozulma sürüyor
            rejim_bandi=1.25,        # HZ2: rejim etkisi medyan oranı [1/1.25, 1.25] içindeyse "rejimden bağımsız"
            ate_iyilesme=0.90)       # HZ4: kol / A0 ATE < 0.90 ve p_holm < 0.05 -> anlamlı iyileşme


def tau0_over_dt(rejim, dr):
    return 0.314 * (D_AP / float(dr)) / (REJIM_HIZ[rejim] * DT)


# ----------------------------------------------------------------------------- okuma
def oku_etiketli(path, etiket, D=None, sayim=None, suz=None):
    """-> D[(etiket|kol, dr, sekans, tohum, pencere)] = kayıt. Bozuk satır ve hücre hatası sayılır, atlanır."""
    D = {} if D is None else D
    sayim = Counter() if sayim is None else sayim
    with open(path) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                sayim["bozuk_satir"] += 1
                continue
            if "arm" not in r:
                sayim["hucre_hatasi"] += 1
                continue
            c = r["cell"]
            if suz is not None and not suz(r):
                sayim["suzulen"] += 1
                continue
            k = (f"{etiket}|{r['arm']}", float(c["D_over_r0"]), c["sequence"], int(c["seed"]), int(r["window"]))
            if k in D:
                sayim["yinelenen"] += 1
            D[k] = r
            sayim["kayit"] += 1
    return D, sayim


def yukle(hizli, ana, rejim, duzeyler, tohumlar, etiket=None):
    """Hızlı rejim kayıtları (yalnız istenen rejim, düzey, tohum, etiket) + ana ızgaranın aynı tohumlardaki kayıtları."""
    etiket = etiket or f"hizli_{rejim}"

    def suz_h(r):
        c = r["cell"]
        dr = float(c["D_over_r0"])
        return (dr in duzeyler and int(c["seed"]) in tohumlar
                and abs(float(c.get("tau0_over_dt", 1.0)) - tau0_over_dt(rejim, dr)) < 1e-9
                and str(c.get("tag", "")) == etiket)

    def suz_a(r):
        c = r["cell"]
        return float(c["D_over_r0"]) in tuple(duzeyler) + (0.0,) and int(c["seed"]) in tohumlar

    D, s_h = oku_etiketli(hizli, rejim, suz=suz_h)
    D, s_a = oku_etiketli(ana, "R0", D=D, suz=suz_a)
    return D, dict(hizli=dict(s_h), ana=dict(s_a))


def butunluk(D, rejim, duzeyler, tohumlar, seqs, kollar=HEPSI):
    """Beklenen (kol, düzey, sekans, tohum) hücrelerinden kaçının bütün pencereleri var; durum sayımı; boş alan."""
    pencere = defaultdict(set)
    for (k, d, s, t, w) in D:
        if k.startswith("R0|") and d in duzeyler:
            pencere[s].add(w)
    eksik, durum, bos, fp = [], Counter(), 0, Counter()
    for kol in kollar:
        for d in duzeyler:
            for s in seqs:
                for t in tohumlar:
                    for w in sorted(pencere[s]):
                        r = D.get((f"{rejim}|{kol}", d, s, t, w))
                        if r is None:
                            eksik.append((kol, d, s, t, w))
                            continue
                        durum[r.get("status")] += 1
                        fp[r.get("code_fp")] += 1
                        if r.get("status") == "ok" and any(r.get(x) is None for x in
                                                            ("m_ate_se3", "m_rpe_trans", "n_nees_pos_mean",
                                                             "n_nees_rot_mean")):
                            bos += 1
    beklenen = len(kollar) * len(duzeyler) * len(tohumlar) * sum(len(pencere[s]) for s in seqs)
    return dict(beklenen_kayit=beklenen, var=beklenen - len(eksik), eksik=len(eksik), eksik_ornek=eksik[:5],
                durum=dict(durum), bos_alanli_ok=bos, code_fp=dict(fp),
                tam=bool(beklenen > 0 and not eksik and bos == 0 and set(fp) <= {FP}))


# ----------------------------------------------------------------------------- aileler
def _holmla(out):
    ph = holm({k: v["p"] for k, v in out.items()})
    for k in out:
        out[k]["p_holm"] = ph[k]
    return out


def aile_ici(D, rejim, metrik, duzeyler, kollar=KOLLAR):
    out = {}
    for kol in kollar:
        for dr in duzeyler:
            a, b, U, dus = eslesik(D, f"{rejim}|{kol}", f"{rejim}|A0", dr, metrik)
            if len(a):
                out[f"{kol}@{dr:g}"] = dict(oran_testi(a, b), dusen_birim=dus)
    return _holmla(out)


def a3_ag(D, rejim, duzeyler):
    out = {}
    for dr in duzeyler:
        a, b, U, dus = eslesik(D, f"{rejim}|A3", f"{rejim}|AG", dr, "ANEES_poz")
        if len(a):
            out[f"{dr:g}"] = dict(oran_testi(a, b), dusen_birim=dus)
    out = _holmla(out)
    for v in out.values():
        v["A3_gecer"] = bool(v["medyan_oran"] < 1 and v["p_holm"] < 0.05)
        v["fark_25_icinde"] = bool(0.75 <= v["medyan_oran"] <= 1.25)
    return out


def rejim_etkisi(D, rejim, metrik, duzeyler, kollar=HEPSI):
    out = {}
    for kol in kollar:
        for dr in duzeyler:
            a, b, U, dus = eslesik(D, f"{rejim}|{kol}", f"R0|{kol}", dr, metrik)
            if len(a):
                out[f"{kol}@{dr:g}"] = dict(oran_testi(a, b), dusen_birim=dus)
    return _holmla(out)


def bozulma(D, rejim, duzeyler):
    """A0: hızlı@düzey / R0@temiz ve R0@düzey / R0@temiz (aynı birimler); kırılma tablosu ana ızgara tanımıyla."""
    out = {"oran": {}}
    for dr in duzeyler:
        for et in (rejim, "R0"):
            a, b, U, dus = eslesik(D, f"{et}|A0", "R0|A0", dr, "ATE", dr_ref=0.0)
            if len(a):
                out["oran"][f"{et}@{dr:g}"] = dict(oran_testi(a, b), dusen_birim=dus)
    for et in (rejim, "R0"):
        kay = [(0.0, r) for (k, d, s, t, w), r in D.items() if k == "R0|A0" and d == 0.0]
        kay += [(d, r) for (k, d, s, t, w), r in D.items() if k == f"{et}|A0" and d in duzeyler]
        if any(d == 0.0 for d, _ in kay) and any(d > 0 for d, _ in kay):
            tab, K = kirilma_tablosu(kay)
            out[f"kirilma_{et}"] = dict(K=K, tablo=tab)
    return out


def basarisizlik(D, rejim, duzeyler, kollar=HEPSI):
    """(rejim, kol, düzey): geçersiz / ıraksamış pencere sayısı ve oranı."""
    out = {}
    for et in (rejim, "R0"):
        for kol in kollar:
            for dr in duzeyler:
                R = [r for (k, d, s, t, w), r in D.items() if k == f"{et}|{kol}" and d == dr]
                if not R:
                    continue
                kotu = sum(0 if (r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged")) else 1
                           for r in R)
                out[f"{et}|{kol}@{dr:g}"] = dict(n=len(R), basarisiz=kotu, oran=kotu / len(R))
    return out


def betimsel(D, rejim, duzeyler, kollar=HEPSI):
    out = {}
    for et in (rejim, "R0"):
        for kol in kollar:
            for dr in duzeyler:
                for metrik in ("ATE", "ANEES_poz", "ANEES_don"):
                    f, tur = METRIK[metrik]
                    g = defaultdict(list)
                    for (k, d, s, t, w), r in D.items():
                        if k == f"{et}|{kol}" and d == dr and f(r) is not None:
                            g[(s, t)].append(f(r))
                    if g:
                        v = np.array([_ort(x, tur) for x in g.values()])
                        out.setdefault(f"{et}|{kol}", {}).setdefault(f"{dr:g}", {})[metrik] = dict(
                            medyan=float(np.median(v)), q25=float(np.percentile(v, 25)),
                            q75=float(np.percentile(v, 75)), n_birim=len(g))
    return out


def anees_bandi(D, rejim, duzeyler, kollar=HEPSI):
    """Pencere düzeyi: ANEES_poz'un %95 ki-kare bandına (3 serbestlik × n) göre konumu."""
    out = {}
    for et in (rejim, "R0"):
        for kol in kollar:
            for dr in duzeyler:
                ici = ust = alt = 0
                for (k, d, s, t, w), r in D.items():
                    if k != f"{et}|{kol}" or d != dr:
                        continue
                    x = anees(r, "pos")
                    if x is None or not r.get("n_n"):
                        continue
                    lo, hi = bant(int(r["n_n"]))
                    ici += lo <= x <= hi
                    ust += x > hi
                    alt += x < lo
                if ici + ust + alt:
                    out[f"{et}|{kol}@{dr:g}"] = dict(bant_ici=ici, bant_ustu=ust, bant_alti=alt)
    return out


# ----------------------------------------------------------------------------- kararlar
def kararlar(S, rejim, duzeyler):
    K = {}
    ust = max(duzeyler)
    o = S["bozulma"]["oran"].get(f"{rejim}@{ust:g}")
    kk = S["bozulma"].get(f"kirilma_{rejim}", {})
    if o:
        K["HZ1"] = dict(oran=o["medyan_oran"], ga95=o["ga95"], kirilma=kk.get("K"),
                        karar=("bozulma sürüyor (oran > %.2f)" % ESIK["bozulma_oran"]
                               if o["medyan_oran"] > ESIK["bozulma_oran"] else "bozulma eşiğin altında"))
    hz2 = {}
    b = ESIK["rejim_bandi"]
    for dr in duzeyler:
        v = S["rejim_ATE"].get(f"A0@{dr:g}")
        if v:
            anl = v["p_holm"] < 0.05
            hz2[f"{dr:g}"] = dict(oran=v["medyan_oran"], ga95=v["ga95"], p_holm=v["p_holm"],
                                  karar=("hızlı rejim daha kötü" if anl and v["medyan_oran"] > 1 else
                                         "hızlı rejim daha iyi" if anl and v["medyan_oran"] < 1 else
                                         "fark anlamlı değil" + (" (GA bant içinde)" if 1 / b <= v["ga95"][0]
                                                                 and v["ga95"][1] <= b else " (GA bandı aşıyor)")))
    K["HZ2"] = hz2
    hz3 = {"i": {}, "ii": S["a3_ag"]}
    for kol in ("A3", "A1", "AG"):
        for dr in duzeyler:
            v = S["aile_ANEES_poz"].get(f"{kol}@{dr:g}")
            if v:
                hz3["i"][f"{kol}@{dr:g}"] = dict(medyan_oran=v["medyan_oran"], p_holm=v["p_holm"],
                                                 anlamli_azalma=bool(v["medyan_oran"] < 1 and v["p_holm"] < 0.05))
    gecer = [k for k, v in S["a3_ag"].items() if v["A3_gecer"]]
    hz3["karar"] = ("A3 AG'yi Holm sonrası anlamlı geçiyor (D/r0: " + ", ".join(gecer) + ")" if gecer
                    else "A3 AG'yi anlamlı geçemiyor")
    K["HZ3"] = hz3
    iyi = sorted(k for k, v in S["aile_ATE"].items()
                 if v["medyan_oran"] < ESIK["ate_iyilesme"] and v["p_holm"] < 0.05)
    K["HZ4"] = dict(anlamli_10_iyilesme=iyi, karar="ÇÜRÜDÜ" if iyi else "TUTTU")
    hz5 = {}
    for kol in ("A1", "A3", "AG"):
        for dr in duzeyler:
            v = S["rejim_ANEES_poz"].get(f"{kol}@{dr:g}")
            if v:
                hz5[f"{kol}@{dr:g}"] = dict(oran=v["medyan_oran"], ga95=v["ga95"], p_holm=v["p_holm"],
                                            anlamli_azalma=bool(v["medyan_oran"] < 1 and v["p_holm"] < 0.05))
    K["HZ5"] = hz5
    return K


def calistir(D, rejim, duzeyler):
    S = dict(aile_ATE=aile_ici(D, rejim, "ATE", duzeyler),
             aile_ANEES_poz=aile_ici(D, rejim, "ANEES_poz", duzeyler),
             aile_ANEES_don=aile_ici(D, rejim, "ANEES_don", duzeyler),
             a3_ag=a3_ag(D, rejim, duzeyler),
             rejim_ATE=rejim_etkisi(D, rejim, "ATE", duzeyler),
             rejim_ANEES_poz=rejim_etkisi(D, rejim, "ANEES_poz", duzeyler),
             bozulma=bozulma(D, rejim, duzeyler),
             basarisizlik=basarisizlik(D, rejim, duzeyler),
             betimsel=betimsel(D, rejim, duzeyler),
             anees_bandi=anees_bandi(D, rejim, duzeyler))
    S["kararlar"] = kararlar(S, rejim, duzeyler)
    return S


def _satir(k, v):
    return (f"   {k:12s} oran {v['medyan_oran']:.3f} [{v['ga95'][0]:.3f}, {v['ga95'][1]:.3f}] n={v['n']} "
            f"p={v['p']:.3g} p_holm={v['p_holm']:.3g} düşen={v.get('dusen_birim', 0)}")


def yazdir(S, rejim):
    for ad, bas in (("aile_ATE", f"F1 ATE kol / A0 ({rejim})"), ("aile_ANEES_poz", f"F1 ANEES_poz kol / A0 ({rejim})"),
                    ("rejim_ATE", f"F3 ATE {rejim} / R0"), ("rejim_ANEES_poz", f"F3 ANEES_poz {rejim} / R0")):
        print(bas)
        for k, v in S[ad].items():
            print(_satir(k, v))
    print(f"F2 ANEES_poz A3 / AG ({rejim})")
    for k, v in S["a3_ag"].items():
        print(_satir("D/r0=" + k, v), "geçer" if v["A3_gecer"] else "-")
    print("A0 ATE düzey / temiz (R0, D/r0 = 0)")
    for k, v in S["bozulma"]["oran"].items():
        print(f"   {k:12s} oran {v['medyan_oran']:.2f} [{v['ga95'][0]:.2f}, {v['ga95'][1]:.2f}] n={v['n']}")
    for et in (rejim, "R0"):
        kk = S["bozulma"].get(f"kirilma_{et}")
        if kk:
            print(f"   kırılma ({et}): K = {kk['K']}; " + "; ".join(
                f"{t['D_over_r0']:g}: {t['oran']:.2f}x, ıraksama {t['irak_oran']:.2f}" for t in kk["tablo"]))
    print("başarısız pencere oranı")
    for k, v in S["basarisizlik"].items():
        if v["basarisiz"]:
            print(f"   {k:16s} {v['basarisiz']}/{v['n']}")
    print("KARARLAR")
    K = S["kararlar"]
    if "HZ1" in K:
        print(f"   HZ1: {K['HZ1']['oran']:.2f} [{K['HZ1']['ga95'][0]:.2f}, {K['HZ1']['ga95'][1]:.2f}], "
              f"K = {K['HZ1']['kirilma']} -> {K['HZ1']['karar']}")
    for k, v in K["HZ2"].items():
        print(f"   HZ2 D/r0={k}: {v['oran']:.3f} [{v['ga95'][0]:.3f}, {v['ga95'][1]:.3f}] p_holm={v['p_holm']:.3g} "
              f"-> {v['karar']}")
    print("   HZ3 (i):", {k: v["anlamli_azalma"] for k, v in K["HZ3"]["i"].items()})
    print(f"   HZ3 (ii): {K['HZ3']['karar']}")
    print(f"   HZ4: anlamlı >%10 ATE iyileşmesi: {K['HZ4']['anlamli_10_iyilesme'] or 'yok'} -> {K['HZ4']['karar']}")
    print("   HZ5:", {k: (round(v["oran"], 3), v["anlamli_azalma"]) for k, v in K["HZ5"].items()})


# ----------------------------------------------------------------------------- öz-test
def oz_test():
    rng = np.random.default_rng(5)
    D, seqs, toh, duz = {}, ("S1", "S2", "S3"), (1, 2, 3, 4), (4.0, 8.0)
    kol_n = {"A0": 1.0, "A1": 0.02, "A3": 0.2, "AG": 0.3, "A0C": 1 / 12, "A0C+A3": 0.02}
    for s in seqs:
        for t in toh:
            for w in (10, 20):
                for dr in (0.0,) + duz:
                    ba = 0.02 * (1 + dr / 4) * (1 + 0.05 * rng.standard_normal())
                    bn = 150 * (1 + dr) * (1 + 0.05 * rng.standard_normal())
                    for et, fa, fn in (("R0", 1.0, 1.0), ("R5", 1.4, 0.5)):
                        if et == "R5" and dr == 0.0:
                            continue
                        for kol, f in kol_n.items():
                            D[(f"{et}|{kol}", dr, s, t, w)] = dict(
                                status="ok", m_valid=True, n_valid=True, n_n=50, m_diverged=False, code_fp=FP,
                                m_ate_se3=ba * fa, m_rpe_trans=0.01, n_nees_pos_mean=bn * fn * f,
                                n_nees_rot_mean=bn * fn * f)
    b = butunluk(D, "R5", duz, toh, seqs)
    assert b["tam"] and b["beklenen_kayit"] == 6 * 2 * 4 * 3 * 2 == b["var"], b
    S = calistir(D, "R5", duz)
    assert abs(S["rejim_ATE"]["A0@8"]["medyan_oran"] - 1.4) < 1e-9 and S["rejim_ATE"]["A0@8"]["n"] == 12
    assert abs(S["rejim_ANEES_poz"]["A3@4"]["medyan_oran"] - 0.5) < 1e-9
    assert abs(S["aile_ANEES_poz"]["A3@8"]["medyan_oran"] - 0.2) < 1e-9
    assert abs(S["a3_ag"]["8"]["medyan_oran"] - 2 / 3) < 1e-9 and S["a3_ag"]["8"]["A3_gecer"]
    assert len(S["aile_ATE"]) == 10 and len(S["rejim_ATE"]) == 12 and len(S["a3_ag"]) == 2
    assert all(abs(v["medyan_oran"] - 1) < 1e-12 and v["p"] == 1.0 for v in S["aile_ATE"].values())
    assert S["kararlar"]["HZ4"]["karar"] == "TUTTU" and S["kararlar"]["HZ2"]["8"]["karar"] == "hızlı rejim daha kötü"
    o = S["bozulma"]["oran"]["R5@8"]["medyan_oran"]
    assert 3.7 < o < 4.7 and S["bozulma"]["kirilma_R5"]["K"] == 4.0, (o, S["bozulma"])     # 1.4 × (1 + 8/4) = 4.2
    # eksik kayıt bütünlükte yakalanmalı
    D2 = dict(D)
    del D2[("R5|A3", 8.0, "S1", 1, 10)]
    assert not butunluk(D2, "R5", duz, toh, seqs)["tam"]
    # tau0/dt tanımı sonda ile aynı
    assert abs(tau0_over_dt("R5", 8) - 0.0157) < 1e-12 and abs(tau0_over_dt("R5", 4) - 0.0314) < 1e-12
    print("öz-test (sentetik) GEÇTİ")
    # ana ızgara sayılarını yeniden üretme: R0 içinde kol / A0, D/r0 4 ve 8, 45 birim -> analiz.json ile aynı ham değer
    if os.path.exists(ANA) and os.path.exists(ANA_ANALIZ):
        Dm, _ = oku_etiketli(ANA, "R0")
        ref = json.load(open(ANA_ANALIZ))
        n = 0
        for metrik in ("ATE", "ANEES_poz"):
            F = aile_ici(Dm, "R0", metrik, (4.0, 8.0))
            for k, v in F.items():
                r = ref[f"aile_{metrik}"][k]
                assert v["n"] == r["n"] == 45 and v["medyan_oran"] == r["medyan_oran"] and v["p"] == r["p"] \
                    and v["ga95"] == r["ga95"], (metrik, k, v, r)
                n += 1
        G = a3_ag(Dm, "R0", (8.0,))
        r = ref["H4"]["ii"]["8"]
        assert G["8"]["medyan_oran"] == r["medyan_oran"] and G["8"]["p"] == r["p"], (G, r)
        E = rejim_etkisi(Dm, "R0", "ATE", (8.0,), kollar=("A0",))
        assert E["A0@8"]["medyan_oran"] == 1.0 and E["A0@8"]["p"] == 1.0
        print(f"öz-test (ana ızgara): {n} aile girdisi + H4 (ii) D/r0 8 analiz.json ile sayı sayı aynı; "
              f"kendisiyle rejim etkisi 1.000 GEÇTİ")
    else:
        print("öz-test (ana ızgara): kayıt ya da analiz.json yok, atlandı")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hizli", default=os.path.join(EV, "records.jsonl"))
    ap.add_argument("--ana", default=ANA)
    ap.add_argument("--rejim", default="R5")
    ap.add_argument("--duzeyler", default="8,4")
    ap.add_argument("--tohumlar", default="1,2,3,4")
    ap.add_argument("--cikti", default=os.path.join(EV, "analiz.json"))
    ap.add_argument("--etiket", default="", help="hücre etiketi (varsayılan hizli_<rejim>; ön test: hizli_ontest_R5)")
    ap.add_argument("--sekanslar", default="", help="yalnız bu sekanslar (ön test çıktısını okumak için)")
    ap.add_argument("--butunluk", action="store_true", help="yalnız bütünlük sayımları; test tablosu basılmaz")
    ap.add_argument("--oz_test", action="store_true")
    a = ap.parse_args()
    if a.oz_test:
        oz_test()
        return 0
    from hizli_kos import SEQS
    duz = tuple(sorted(float(x) for x in a.duzeyler.split(",") if x))
    toh = tuple(int(x) for x in a.tohumlar.split(",") if x)
    seqs = [x for x in a.sekanslar.split(",") if x] or SEQS
    D, sayim = yukle(a.hizli, a.ana, a.rejim, duz, toh, a.etiket or None)
    b = butunluk(D, a.rejim, duz, toh, seqs)
    print(f"okuma: {sayim}")
    print(f"bütünlük: {b['var']}/{b['beklenen_kayit']} kayıt; eksik {b['eksik']} {b['eksik_ornek']}; durum {b['durum']}; "
          f"boş alanlı 'ok' {b['bos_alanli_ok']}; code_fp {b['code_fp']}; tam = {b['tam']}")
    if a.butunluk:
        return 0 if b["tam"] else 3
    if not b["tam"]:
        print("DUR: veri tam değil; karar tabloları yalnız tam veriyle üretilir (--butunluk ile sayımlara bakın)")
        return 3
    if len(seqs) * len(toh) < 6:
        print(f"DUR: {len(seqs) * len(toh)} birim; test tablosu en az 6 birimle basılır (ön test çıktısı için --butunluk)")
        return 3
    S = calistir(D, a.rejim, duz)
    out = dict(kaynak=dict(hizli=a.hizli, ana=a.ana, rejim=a.rejim, duzeyler=list(duz), tohumlar=list(toh),
                           hiz_m_s=REJIM_HIZ[a.rejim], tau0_over_dt={f"{d:g}": tau0_over_dt(a.rejim, d) for d in duz},
                           sayim=sayim), butunluk=b, esikler=ESIK, **S)
    json.dump(out, open(a.cikti, "w"), indent=1, ensure_ascii=False, default=float)
    yazdir(S, a.rejim)
    print(f"-> {a.cikti}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
