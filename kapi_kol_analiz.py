#!/usr/bin/env python
"""Kapı kolu geliştirme koşumunun analizi (evidence/KAPI/ON_BEYAN_3.md §3). KEŞİF: tohum 0, birim = sekans, test yok.
ATE pencerelerde aritmetik, konum ANEES'i geometrik ortalama; etki = sekanslarda eşli medyan oran.

    python kapi_kol_analiz.py        # veri tam değilse yalnız bütünlük dökümü (oran tablosu basılmaz)
"""
import collections
import json
import math
import os

import numpy as np

import kapi_kol as kk

ADAYLAR = list(kk.GELISTIRME)


def ate(r):
    v = r.get("m_ate_se3")
    return float(v) if (r.get("status") == "ok" and r.get("m_valid") and v is not None and math.isfinite(v) and v > 0) else None


def anees(r):
    v = r.get("n_nees_pos_mean")
    return float(v) if (r.get("status") == "ok" and r.get("n_valid") and v is not None and math.isfinite(v) and v > 0) else None


def basarili(r):
    return r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged")


def main():
    C = kk.oku()
    S = dict(kk.basvuru())                      # (ad, sekans, pencere) -> kayıt  [eski kod: A0, A1, A3S, A0g, ...]
    eksik, n_kayit = [], 0
    for v in ADAYLAR:
        for s in kk.SEQS:
            recs = C.get((v, s), [])
            if not kk.bitti(recs):
                eksik.append((v, s, len(recs)))
            for r in recs:
                if "arm" in r:
                    n_kayit += 1
                    S[(v, s, int(r["window"]))] = r
    fps = collections.Counter(r.get("code_fp") for (v, s), recs in C.items() if v in ADAYLAR for r in recs)
    print(f"[analiz] {n_kayit} kayıt (beklenen {3 * len(ADAYLAR) * len(kk.SEQS)}), eksik koşum {len(eksik)}; code_fp {dict(fps)}")
    for v in ADAYLAR:
        var = [r for (k, s, w), r in S.items() if k == v]
        print(f"  {v:3s} kayıt {len(var):2d}; status ok {sum(r.get('status') == 'ok' for r in var)}; m_valid "
              f"{sum(bool(r.get('m_valid')) for r in var)}; n_valid {sum(bool(r.get('n_valid')) for r in var)}; "
              f"ıraksayan {sum(bool(r.get('m_diverged')) for r in var)}")
    if eksik:
        print(f"[analiz] VERİ TAM DEĞİL, oran tablosu basılmıyor. Eksik: {eksik}")
        return

    def birim(kol, f):
        g = collections.defaultdict(dict)
        for (k, s, w), r in S.items():
            if k == kol:
                x = f(r)
                if x is not None:
                    g[s][w] = x
        return g

    def oran(pay, payda, f, tur):
        A, B = birim(pay, f), birim(payda, f)
        lr, ad = [], []
        for s in kk.SEQS:
            ws = sorted(set(A.get(s, {})) & set(B.get(s, {})))
            if not ws:
                continue
            xa, xb = np.array([A[s][w] for w in ws]), np.array([B[s][w] for w in ws])
            lr.append(float(np.mean(np.log(xa)) - np.mean(np.log(xb))) if tur == "geo"
                      else float(np.log(np.mean(xa) / np.mean(xb))))
            ad.append(s)
        if not lr:
            return None
        lr = np.array(lr)
        return dict(med=math.exp(float(np.median(lr))), lt1=int(np.sum(lr < 0)), n=len(lr),
                    mn=math.exp(float(lr.min())), mx=math.exp(float(lr.max())), sek=dict(zip(ad, np.exp(lr).round(3))))

    def yaz(o):
        return "      -" if o is None else f"{o['med']:6.3f} [{o['mn']:.2f}, {o['mx']:.2f}] ({o['lt1']}/{o['n']})"

    ciftler = [(v, p) for p in ("A0", "A3S", "A1", "A0g") for v in ADAYLAR] + [("A0g", "A0"), ("A1", "A0"), ("A3S", "A0"),
                                                                              ("D2", "D1")]
    for baslik, f, tur in (("ATE (m_ate_se3)", ate, "ari"), ("konum ANEES'i (n_nees_pos_mean)", anees, "geo")):
        print(f"\n== {baslik}: medyan eşli oran [en küçük, en büyük sekans] (oranı < 1 olan sekans / sekans) ==")
        for pay, payda in ciftler:
            print(f"  {pay:3s} / {payda:3s}: {yaz(oran(pay, payda, f, tur))}")
    print("\n== sekans bazında ATE oranı (kol / A0) ==")
    for v in ADAYLAR + ["A0g", "A3S"]:
        o = oran(v, "A0", ate, "ari")
        print(f"  {v:3s}: {o['sek'] if o else None}")

    print("\n== düzeyler: medyan ATE [m] / medyan konum ANEES'i / başarısız pencere / n ==")
    for kol in ["A0", "A1", "A3S", "A0g"] + ADAYLAR:
        v = [r for (k, s, w), r in S.items() if k == kol]
        a = [ate(r) for r in v if ate(r) is not None]
        n = [anees(r) for r in v if anees(r) is not None]
        print(f"  {kol:3s} {np.median(a) if a else float('nan'):8.4f} / {np.median(n) if n else float('nan'):8.1f} / "
              f"{sum(not basarili(r) for r in v)} / {len(v)}")

    print("\n== kapı sigması ve sigma_hat (sekans başına; arm_stats) ==")
    for v in ADAYLAR:
        med, mx, sh = [], [], []
        for s in kk.SEQS:
            st = (S[(v, s, min(w for (k, ss, w) in S if k == v and ss == s))].get("arm_stats") or {})
            med.append(st.get("gate_sigma_median", float("nan")))
            mx.append(st.get("gate_sigma_max", float("nan")))
            sh.append(st.get("sigma_hat_median", float("nan")))
        ref_sh = []
        for s in kk.SEQS:
            st = (S[("A3S", s, min(w for (k, ss, w) in S if k == "A3S" and ss == s))].get("arm_stats") or {})
            ref_sh.append(st.get("sigma_hat_median", float("nan")))
        print(f"  {v}: kapı sigması medyanı (sekanslar) {np.round(med, 3).tolist()}; sekans medyanı {np.median(med):.3f}; "
              f"en büyük maks {np.max(mx):.3f}")
        print(f"      sigma_hat medyanı {np.round(sh, 3).tolist()} (kapısız A3S: {np.round(ref_sh, 3).tolist()}; "
              f"oracle sigma_tilt {kk.sigma_tilt(kk.DR):.3f}, toplam {kk.CARPAN * kk.sigma_tilt(kk.DR):.3f})")

    print("\n== ön beyan: beklenen değerler ve seçim kuralı ==")
    def bas(v):
        var = [r for (k, s, w), r in S.items() if k == v]
        return sum(not basarili(r) for r in var), len(var)

    gecer = {}
    for v, alt, ust in (("D2", 0.50, 0.85), ("D1", 0.60, 0.95)):
        o = oran(v, "A0", ate, "ari")
        b, n = bas(v)
        gecer[v] = bool(o and o["med"] <= 0.85 and o["lt1"] >= 7 and b <= 1)
        print(f"  {'G1' if v == 'D2' else 'G2'} {v} / A0 ATE: {yaz(o)}; beklenen {alt:.2f}-{ust:.2f} "
              f"{'tuttu' if o and alt <= o['med'] <= ust else 'tutmadı'}; başarısız {b}/{n}; GEÇER mi: "
              f"{'EVET' if gecer[v] else 'HAYIR'}")
    for v in ADAYLAR:
        print(f"  G3 {v} / A3S ATE: {yaz(oran(v, 'A3S', ate, 'ari'))}; beklenen < 1 (D2: 0.5-0.9)")
        print(f"  G4 {v} / A0g ATE: {yaz(oran(v, 'A0g', ate, 'ari'))}; beklenen 1.0-1.6")
        print(f"  G5 {v} / A0 konum ANEES'i: {yaz(oran(v, 'A0', anees, 'geo'))}; beklenen < 0.5")
    o1, o2 = oran("D1", "A0", ate, "ari"), oran("D2", "A0", ate, "ari")
    if gecer["D1"] and gecer["D2"]:
        q = o2["med"] / o1["med"]
        sec = "D1" if 0.90 <= q <= 1.10 else ("D2" if o2["med"] < o1["med"] else "D1")
        print(f"  seçim: ikisi de geçer; medyan oranların oranı D2/D1 = {q:.3f} -> {sec}")
    elif gecer["D1"] or gecer["D2"]:
        sec = "D1" if gecer["D1"] else "D2"
        print(f"  seçim: yalnız {sec} geçer -> {sec}")
    else:
        print("  seçim: hiçbiri geçmedi -> çevrimiçi kol bu hâliyle doğrulayıcı ızgaraya götürülmez (kural 3)")


if __name__ == "__main__":
    main()
