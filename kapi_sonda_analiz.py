#!/usr/bin/env python
"""Keşif sondası 2 analizi (evidence/SONDA_KESIF/ON_BEYAN_2.md §2-§3). Birim: sekans; ATE pencerelerde aritmetik, konum
ANEES'i geometrik ortalama; etki = sekanslarda eşli medyan oran (sonda_analiz.py ile aynı kural). Test yok.

    python kapi_sonda_analiz.py            # veri tam değilse yalnız bütünlük dökümü (oran tablosu basılmaz)
"""
import collections
import json
import math
import os

import numpy as np

import kapi_sonda as ks
import sonda_kesif as sk

BEKLENEN = {"K": list(ks.KONTROL_SEQ), "Gd": list(sk.SEQS), "Gn": list(sk.SEQS), "Gp": list(sk.SEQS)}


def ate(r):
    v = r.get("m_ate_se3")
    return float(v) if (r.get("status") == "ok" and r.get("m_valid") and v is not None and math.isfinite(v) and v > 0) else None


def anees(r):
    v = r.get("n_nees_pos_mean")
    return float(v) if (r.get("status") == "ok" and r.get("n_valid") and v is not None and math.isfinite(v) and v > 0) else None


def basarili(r):
    return r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged")


def main():
    C = ks.oku_cikti()
    ref = ks.sonda_basvuru()
    S = {}                                   # (kol, sekans, pencere) -> kayıt
    for (kol, s, w), r in ref.items():
        S[(kol, s, w)] = r
    n_kayit = n_hata = 0
    eksik = []
    for v, seqs in BEKLENEN.items():
        for s in seqs:
            recs = C.get((v, s), [])
            if not ks.bitti(recs):
                eksik.append((v, s, len(recs)))
            for r in recs:
                if "arm" not in r:
                    n_hata += 1
                    continue
                n_kayit += 1
                S[(v, s, int(r["window"]))] = r
    fps = collections.Counter(r.get("code_fp") for recs in C.values() for r in recs)
    print(f"[analiz2] {ks.OUT}: {n_kayit} kayıt (beklenen {3 * sum(len(x) for x in BEKLENEN.values())}), "
          f"hücre hatası {n_hata}, eksik koşum {len(eksik)}; code_fp {dict(fps)}")
    for v in BEKLENEN:
        var = [r for (k, s, w), r in S.items() if k == v]
        print(f"  {v:3s} kayıt {len(var):2d}; status ok {sum(r.get('status') == 'ok' for r in var)}; m_valid "
              f"{sum(bool(r.get('m_valid')) for r in var)}; n_valid {sum(bool(r.get('n_valid')) for r in var)}; "
              f"ıraksayan {sum(bool(r.get('m_diverged')) for r in var)}")
    if eksik:
        print(f"[analiz2] VERİ TAM DEĞİL, oran tablosu basılmıyor. Eksik: {eksik}")
        return

    # --- B1 ---
    n_ayni = n_k = 0
    for (k, s, w), r in S.items():
        if k == "K":
            n_k += 1
            m = S.get(("A0g", s, w)) or {}
            n_ayni += bool(m) and all(r.get(f) == m.get(f) for f in ks.BIT)
    print(f"\n== B1: K kayıtları ilk sondanın A0g kayıtlarıyla bit-aynı: {n_ayni}/{n_k} ==")

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
        for s in sk.SEQS:
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

    for baslik, f, tur in (("ATE (m_ate_se3)", ate, "ari"), ("konum ANEES'i (n_nees_pos_mean)", anees, "geo")):
        print(f"\n== {baslik}: medyan eşli oran [en küçük, en büyük sekans] (oranı < 1 olan sekans / sekans) ==")
        for pay, payda in (("A0g", "A0"), ("A1", "A0"), ("Gd", "A0"), ("Gn", "A0"), ("Gp", "A0"),
                           ("Gd", "A0g"), ("Gn", "A0g"), ("Gp", "A0g"), ("Gd", "A1"), ("Gn", "A1")):
            print(f"  {pay:3s} / {payda:3s}: {yaz(oran(pay, payda, f, tur))}")
    print("\n== sekans bazında ATE oranı (varyant / A0) ==")
    for v in ("A0g", "Gd", "Gn", "Gp"):
        o = oran(v, "A0", ate, "ari")
        print(f"  {v:3s}: {o['sek'] if o else None}")

    print("\n== düzeyler: medyan ATE [m] / medyan konum ANEES'i / başarısız pencere / n ==")
    for kol in ("A0", "A1", "A0g", "Gd", "Gn", "Gp"):
        v = [r for (k, s, w), r in S.items() if k == kol]
        a = [ate(r) for r in v if ate(r) is not None]
        n = [anees(r) for r in v if anees(r) is not None]
        print(f"  {kol:3s} {np.median(a) if a else float('nan'):8.4f} / {np.median(n) if n else float('nan'):8.1f} / "
              f"{sum(not basarili(r) for r in v)} / {len(v)}")

    print("\n== ön beyan kuralları (B3-B6) ==")
    kural = lambda o: o is not None and o["med"] <= 0.85 and o["lt1"] >= 7      # noqa: E731
    gd, gn, gp = (oran(v, "A0", ate, "ari") for v in ("Gd", "Gn", "Gp"))
    print(f"  B3 (Q1) Gd / A0 ATE: {yaz(gd)}; beklenen 0.50-0.80; kural (<= 0.85 ve >= 7/9): {'TUTTU' if kural(gd) else 'TUTMADI'}")
    print(f"          Gd / A0g ATE: {yaz(oran('Gd', 'A0g', ate, 'ari'))}; beklenen 1.0-1.5")
    print(f"  B4 (Q2) Gn / A0 ATE: {yaz(gn)}; beklenen 0.55-0.90; taşıyor mu: {'EVET' if kural(gn) else 'HAYIR'}")
    print(f"          Gp / A0 ATE: {yaz(gp)}; beklenen 0.85-1.15; taşıyor mu: {'EVET' if kural(gp) else 'HAYIR'}")
    for v, bek in (("Gn", "< 0.3"), ("Gd", "< 0.3"), ("Gp", "0.5-1.5")):
        print(f"  B5 {v} / A0 konum ANEES'i: {yaz(oran(v, 'A0', anees, 'geo'))}; beklenen {bek}")
    for v in ("Gd", "Gn", "Gp"):
        var = [r for (k, s, w), r in S.items() if k == v]
        print(f"  B6 {v}: ıraksayan / başarısız pencere {sum(not basarili(r) for r in var)}/{len(var)} (beklenen <= 1/27)")


if __name__ == "__main__":
    main()
