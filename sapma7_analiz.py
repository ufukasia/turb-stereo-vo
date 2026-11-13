#!/usr/bin/env python
"""Kapılı genel uyarlamalı taban çizgisi (AGg) ızgarasının analizi (ON_KAYIT_FAZ2.md §5 Sapma 7).

Izgaranın hiçbir sonucu görülmeden yazıldı. İstatistik Sapma 6 ile aynıdır (turb/faz2_analiz.py): birim (sekans,
tohum); eşli medyan oran, birim bootstrap 10 000 ile %95 GA; iki yönlü Wilcoxon; Holm aile içinde; TOST ±%15.

    python sapma7_analiz.py --oz_test
    python sapma7_analiz.py --butunluk             # koşum sürerken: yalnız sayımlar
    python sapma7_analiz.py --ontest --butunluk    # ön test çıktısının okunması
    python sapma7_analiz.py                        # yalnız veri tam ve doğrulamalar tutuyorken tablo üretir

Etiketler: A0, AG (kapısız; ana ızgara R0 / hızlı ızgara R5, eski kod), KAPI, A1g (Sapma 6 kayıtları), AGg ve A0y
(bu ızgara). Aileler:
  HG1  AGg / A0, ATE; dört türbülans koşulu (m = 4)
  HG2  KAPI / AGg, ATE; dört koşul (m = 4), iki yönlü + TOST
  HG3  AGg / A0, ATE; temiz: %95 GA üst sınırı < 1.10
  S1   AGg / AG (kapısız), ATE (m = 4)          S2  AGg / A1g, ATE (m = 4)
  S3   ANEES_poz: AGg / A0, KAPI / AGg, AGg / AG (m = 12)
  T    temiz, betimsel: KAPI / AGg, AGg / AG (ATE); AGg / A0 (ANEES_poz)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
from turb.faz2_analiz import METRIK, _ort, eslesik, holm, oran_testi  # noqa: E402
import kapi_dogrula_analiz as ka  # noqa: E402
from a3s_analiz import tost  # noqa: E402

EV7 = os.path.join(ROOT, "evidence", "SAPMA7")
KOSUL, TURB, BIT, FP_ANA = ka.KOSUL, ka.TURB, ka.BIT, ka.FP_ANA
TOHUMLAR = (1, 2, 3, 4)
YENI = {"I": ("AGg",), "II": ("A0y", "AGg"), "III": ("A0y", "AGg")}
ESIK = dict(alfa=0.05, hg3_ust_sinir=1.10)
ARGV = {"A0y": ["--turb_cal_factor", "12.0"], "AGg": ["--turb_cal_factor", "12.0", "--turb_gate", "3"],
        "KAPI": ["--turb_cal_factor", "1.0", "--turb_dimm_sigma_scale", "1.0", "--turb_gate", "3"]}


def oku_ag(ana, tohumlar, D):
    """Kapısız AG başvuruları (eski kod): ana ızgara (R0) ve hızlı ızgara (R5)."""
    S = Counter()
    for r in ka._satirlar(os.path.join(ana, "FAZ2", "records.jsonl"), S):
        c = r.get("cell") or {}
        if ("arm" in r and c.get("tag") == "faz2" and r["arm"] == "AG" and int(c["seed"]) in tohumlar
                and float(c["D_over_r0"]) in (0.0, 8.0, 16.0) and abs(float(c.get("tau0_over_dt", 1.0)) - 1.0) < 1e-12
                and r.get("code_fp") == FP_ANA):
            ka._koy(D, S, "R0|AG", c["D_over_r0"], c["sequence"], c["seed"], r["window"], r)
    for r in ka._satirlar(os.path.join(ana, "HIZLI", "records.jsonl"), S):
        c = r.get("cell") or {}
        if ("arm" in r and c.get("tag") == "hizli_R5" and r["arm"] == "AG" and int(c["seed"]) in tohumlar
                and float(c["D_over_r0"]) in (8.0, 4.0) and r.get("code_fp") == FP_ANA):
            ka._koy(D, S, f"R5|AG", c["D_over_r0"], c["sequence"], c["seed"], r["window"], r)
    return dict(S)


def oku_s7(yollar, tohumlar, D):
    D7, S = ka.oku_dogrula(yollar, tohumlar, D={})
    for (lab, dr, s, t, w), r in D7.items():
        rej, kol = lab.split("|")
        D[(f"{rej}|{'KAPIy' if kol == 'KAPI' else kol}", dr, s, t, w)] = r        # ön testteki KAPI yinelemesi: KAPIy
    return S


def beklenen_kollar(parca, ontest):
    """(etiket, tür) listesi; tür: yeni (bu ızgara, yeni kod), s6 (Sapma 6, yeni kod), eski (başvuru)."""
    out = [(k, "yeni") for k in YENI[parca]] + [("A0", "eski")]
    if not ontest:
        out += [("AG", "eski"), ("KAPI", "s6")] + ([("A1g", "s6")] if parca != "III" else [])
    return out


def butunluk(D, kosullar, tohumlar, seqs, pencere, fp_yeni, ontest=False):
    tablo, eksik = {}, []
    for ad in kosullar:
        rej, dr, parca = KOSUL[ad]
        kollar = beklenen_kollar(parca, ontest)
        if ontest and "A0y" not in [k for k, _ in kollar]:
            kollar.append(("A0y", "yeni"))                      # ön testte R5 hücresinde de A0 koşar
        for kol, tur in kollar:
            durum, fp, bos, nb, nv = Counter(), Counter(), 0, 0, 0
            for s in seqs:
                for t in tohumlar:
                    for w in pencere[s]:
                        nb += 1
                        r = D.get((f"{rej}|{kol}", dr, s, t, int(w)))
                        if r is None:
                            eksik.append((ad, kol, s, t, int(w)))
                            continue
                        nv += 1
                        durum[r.get("status")] += 1
                        fp[r.get("code_fp")] += 1
                        bos += r.get("status") == "ok" and any(r.get(x) is None for x in BIT)
            tablo[f"{ad}|{kol}"] = dict(tur=tur, beklenen=nb, var=nv, durum=dict(durum), bos_alanli_ok=int(bos),
                                        code_fp=dict(fp), fp_tamam=bool(set(fp) <= {FP_ANA if tur == "eski" else fp_yeni}))
    nb = sum(v["beklenen"] for v in tablo.values())
    tam = bool(nb > 0 and not eksik and all(v["bos_alanli_ok"] == 0 and v["fp_tamam"] for v in tablo.values()))
    return dict(tablo=tablo, beklenen_kayit=nb, var=nb - len(eksik), eksik=len(eksik), eksik_ornek=eksik[:6], tam=tam)


def dogrulamalar(D, kosullar):
    """V1 A0y ≡ A0 (dört alan); V2 yeniden üretilen ön-uç istatistikleri; V3 AGg kapı istatistiği (gate = 3, en küçük
    sigma = 1.0, en büyük > 1.0); V5 bayraklar; V6 (yalnız ön test) KAPIy ≡ Sapma 6 KAPI kaydı varsa bit-aynı."""
    V = {k: Counter() for k in ("V1", "V2", "V3", "V5", "V6")}
    for (lab, dr, s, t, w), r in D.items():
        rej, kol = lab.split("|")
        if kol not in ("A0y", "AGg", "KAPIy") or not any(KOSUL[a][0] == rej and KOSUL[a][1] == dr for a in kosullar):
            continue
        k = r.get("kapid") or {}
        ref = D.get((f"{rej}|A0", dr, s, t, w))
        if kol == "A0y":
            V["V1"]["basvuru_yok" if ref is None else ("ayni" if all(r.get(x) == ref.get(x) for x in BIT) else "FARKLI")] += 1
            if k.get("kaynak") == "fe":
                a = ka._fe_ayni(r.get("frontend_stats"), (ref or {}).get("frontend_stats"))
                V["V2"]["basvuru_yok" if a is None else ("ayni" if a else "FARKLI")] += 1
        st = r.get("arm_stats") or {}
        if kol == "AGg":
            mn, mx = st.get("gate_sigma_min"), st.get("gate_sigma_max")
            iyi = st.get("gate") == 3 and mn is not None and abs(mn - 1.0) < 1e-9 and mx is not None and mx > 1.0
            V["V3"]["tamam" if iyi else "HATALI"] += 1
        if kol == "KAPIy":
            r6 = D.get((f"{rej}|KAPI", dr, s, t, w))
            if r6 is not None:
                V["V6"]["ayni" if all(r.get(x) == r6.get(x) for x in BIT) else "FARKLI"] += 1
        V["V5"]["tamam" if list(k.get("argv") or []) == ARGV["KAPI" if kol == "KAPIy" else kol] else "HATALI"] += 1
    kotu = sum(v.get("FARKLI", 0) + v.get("HATALI", 0) + v.get("basvuru_yok", 0) for v in V.values())
    return dict({k: dict(v) for k, v in V.items()}, tamam=bool(kotu == 0))


def test(D, pay, payda, kosul, metrik, es=False):
    rej, dr, _ = KOSUL[kosul]
    a, b, U, dus = eslesik(D, f"{rej}|{pay}", f"{rej}|{payda}", dr, metrik)
    if not len(a):
        return None
    out = dict(oran_testi(a, b), dusen_birim=dus)
    if es:
        out.update(tost(a, b))
    return out


def aile(D, ciftler, kosullar, metrik, es=False):
    out = {}
    for pay, payda in ciftler:
        for ad in kosullar:
            t = test(D, pay, payda, ad, metrik, es)
            if t is not None:
                out[f"{pay}/{payda}@{ad}"] = t
    ph = holm({k: v["p"] for k, v in out.items()})
    for k in out:
        out[k]["p_holm"] = ph[k]
    return out


def betimsel(D, kosullar):
    out = {}
    for ad in kosullar:
        rej, dr, _ = KOSUL[ad]
        for kol in ("A0", "AG", "AGg", "KAPI", "A1g"):
            R = {k: r for k, r in D.items() if k[0] == f"{rej}|{kol}" and k[1] == dr}
            if not R:
                continue
            e = dict(n=len(R), basarisiz=sum(0 if (r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged"))
                                             else 1 for r in R.values()))
            for metrik in ("ATE", "ANEES_poz"):
                f, tur = METRIK[metrik]
                g = defaultdict(list)
                for (k, d, s, t, w), r in R.items():
                    if f(r) is not None:
                        g[(s, t)].append(f(r))
                if g:
                    e[metrik] = float(np.median([_ort(x, tur) for x in g.values()]))
            if kol in ("AGg", "KAPI"):
                g = {}
                for (k, d, s, t, w), r in sorted(R.items()):
                    g.setdefault((s, t), r.get("arm_stats") or {})
                med = np.array([x.get("gate_sigma_median", np.nan) for x in g.values()], float)
                mx = np.array([x.get("gate_sigma_max", np.nan) for x in g.values()], float)
                e.update(kapi_sigma_medyan=float(np.nanmedian(med)), kapi_sigma_aralik=[float(np.nanmin(med)),
                         float(np.nanmax(med))], kapi_sigma_maks=float(np.nanmax(mx)))
            out.setdefault(ad, {})[kol] = e
    return out


def kararlar(S):
    K, a = {}, ESIK["alfa"]
    for ad in TURB:
        v = S["HG1"].get(f"AGg/A0@{ad}")
        if v:
            anl = v["p_holm"] < a
            K[f"HG1@{ad}"] = ("TUTTU" if anl and v["medyan_oran"] < 1 else "TERS YÖN" if anl and v["medyan_oran"] > 1
                              else "TUTMADI")
        v = S["HG2"].get(f"KAPI/AGg@{ad}")
        if v:
            anl = v["p_holm"] < a
            K[f"HG2@{ad}"] = (("KAPI daha doğru" if v["medyan_oran"] < 1 else "AGg daha doğru") if anl else
                              ("ayırt edilemedi; TOST eşdeğer" if v.get("esdeger") else "ayırt edilemedi; eşdeğerlik gösterilemedi"))
    v = S["HG3"].get("AGg/A0@R0-0")
    if v:
        K["HG3"] = "GEÇTİ" if v["ga95"][1] < ESIK["hg3_ust_sinir"] else "KALDI"
    return K


def calistir(D):
    S = dict(HG1=aile(D, [("AGg", "A0")], TURB, "ATE"), HG2=aile(D, [("KAPI", "AGg")], TURB, "ATE", es=True),
             HG3=aile(D, [("AGg", "A0")], ("R0-0",), "ATE"), S1=aile(D, [("AGg", "AG")], TURB, "ATE"),
             S2=aile(D, [("AGg", "A1g")], TURB, "ATE"),
             S3=aile(D, [("AGg", "A0"), ("KAPI", "AGg"), ("AGg", "AG")], TURB, "ANEES_poz"),
             T_ATE=aile(D, [("KAPI", "AGg"), ("AGg", "AG")], ("R0-0",), "ATE", es=True),
             T_ANEES=aile(D, [("AGg", "A0"), ("KAPI", "AGg")], ("R0-0",), "ANEES_poz"),
             betimsel=betimsel(D, tuple(KOSUL)))
    S["kararlar"] = kararlar(S)
    return S


def yazdir(S):
    for ad in ("HG1", "HG2", "HG3", "S1", "S2", "S3", "T_ATE", "T_ANEES"):
        print(ad)
        for k, v in S[ad].items():
            ek = (f" | TOST ort {v['ort_oran']:.3f} GA90 [{v['ga90_t'][0]:.3f}, {v['ga90_t'][1]:.3f}] p {v['p_tost']:.3g} "
                  f"{'EŞDEĞER' if v['esdeger'] else '-'}") if "p_tost" in v else ""
            print(ka._satir(k, v) + ek)
    print("betimsel: kol n/başarısız, birim medyanı ATE [m], ANEES_poz, kapı sigması medyanı [aralık] maks")
    for ad, kol in S["betimsel"].items():
        for k, v in kol.items():
            g = (f" kapı σ {v['kapi_sigma_medyan']:.3f} [{v['kapi_sigma_aralik'][0]:.3f}, {v['kapi_sigma_aralik'][1]:.3f}] "
                 f"maks {v['kapi_sigma_maks']:.2f}") if "kapi_sigma_medyan" in v else ""
            print(f"   {ad:6s} {k:5s} {v['basarisiz']}/{v['n']} ATE {v.get('ATE', float('nan')):.4f} ANEES {v.get('ANEES_poz', float('nan')):.1f}{g}")
    print("KARARLAR")
    for k, v in S["kararlar"].items():
        print(f"   {k}: {v}")


def oz_test():
    rng = np.random.default_rng(7)
    seqs, toh, pencere = ("S1", "S2", "S3"), (1, 2, 3, 4), {s: (10, 20) for s in ("S1", "S2", "S3")}
    carp = dict(AG=1.0, AGg=0.6, KAPI=0.63, A1g=0.5)
    D = {}
    for ad, (rej, dr, parca) in KOSUL.items():
        for s in seqs:
            for t in toh:
                for w in pencere[s]:
                    ba, bn = 0.02 * math.exp(0.3 * rng.standard_normal()), 150 * math.exp(0.3 * rng.standard_normal())
                    for kol, tur in beklenen_kollar(parca, False):
                        f = 1.0 if kol in ("A0", "A0y") else carp[kol] * (math.exp(0.02 * rng.standard_normal()) if kol == "KAPI" else 1)
                        r = dict(status="ok", m_valid=True, n_valid=True, m_diverged=False, m_rpe_trans=0.01,
                                 code_fp=FP_ANA if tur == "eski" else "yeni", m_ate_se3=ba * f, n_nees_pos_mean=bn * f * f,
                                 n_nees_rot_mean=bn * f, frontend_stats={"a": 1})
                        if tur == "yeni":
                            r["kapid"] = dict(argv=ARGV[kol], kaynak="fe")
                            r["arm_stats"] = dict(gate=3, gate_sigma_min=1.0, gate_sigma_max=2.5, gate_sigma_median=2.0) if kol == "AGg" else {}
                        if kol == "KAPI":
                            r["arm_stats"] = dict(gate=3, gate_sigma_median=2.0, gate_sigma_max=3.0)
                        D[(f"{rej}|{kol}", dr, s, t, w)] = r
    b = butunluk(D, tuple(KOSUL), toh, seqs, pencere, "yeni")
    assert b["tam"] and b["beklenen_kayit"] == 24 * (2 * 5 + 2 * 6 + 5), b["beklenen_kayit"]
    V = dogrulamalar(D, tuple(KOSUL))
    assert V["tamam"] and V["V1"] == {"ayni": 72} and V["V3"] == {"tamam": 120}, V
    S = calistir(D)
    assert (len(S["HG1"]), len(S["HG2"]), len(S["S1"]), len(S["S2"]), len(S["S3"])) == (4, 4, 4, 4, 12)
    for ad in TURB:
        assert abs(S["HG1"][f"AGg/A0@{ad}"]["medyan_oran"] - 0.6) < 1e-12 and S["kararlar"][f"HG1@{ad}"] == "TUTTU"
        assert abs(S["S2"][f"AGg/A1g@{ad}"]["medyan_oran"] - 1.2) < 1e-12
        assert abs(S["S3"][f"AGg/A0@{ad}"]["medyan_oran"] - 0.36) < 1e-12
        assert S["HG2"][f"KAPI/AGg@{ad}"]["esdeger"] and "p_tost" in S["HG2"][f"KAPI/AGg@{ad}"]
    assert S["kararlar"]["HG3"] == "GEÇTİ" and abs(S["HG3"]["AGg/A0@R0-0"]["medyan_oran"] - 0.6) < 1e-12
    D2 = {k: (dict(r, m_ate_se3=r["m_ate_se3"] * 2.0) if k[0].endswith("|AGg") else r) for k, r in D.items()}
    K2 = calistir(D2)["kararlar"]
    assert K2["HG1@R0-8"] == "TERS YÖN" and K2["HG3"] == "KALDI" and K2["HG2@R5-8"] == "KAPI daha doğru", K2
    D3 = dict(D)
    del D3[("R5|AGg", 4.0, "S2", 3, 20)]
    assert not butunluk(D3, tuple(KOSUL), toh, seqs, pencere, "yeni")["tam"]
    k0 = ("R0|A0y", 8.0, "S1", 1, 10)
    assert not dogrulamalar({**D, k0: dict(D[k0], m_rpe_trans=0.011)}, tuple(KOSUL))["tamam"]
    k1 = ("R5|AGg", 8.0, "S1", 2, 10)
    assert not dogrulamalar({**D, k1: dict(D[k1], arm_stats=dict(gate=3, gate_sigma_min=1.2, gate_sigma_max=2))}, tuple(KOSUL))["tamam"]
    print("öz-test (sentetik) GEÇTİ")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--butunluk", action="store_true")
    ap.add_argument("--ontest", action="store_true")
    ap.add_argument("--oz_test", action="store_true")
    ap.add_argument("--cikti", default=os.path.join(EV7, "s7_analiz.json"))
    a = ap.parse_args()
    if a.oz_test:
        oz_test()
        return 0
    import kapi_dogrula as kd
    import sapma7_kos as s7
    fp = kd.code_fingerprint()
    pencere = {s: tuple(int(w) for w in kd.windows_for(s, kd.FAZ2_WINDOWS_FILE)) for s in kd.SEQS}
    if a.ontest:
        kosullar, toh, seqs = ("R0-8", "R5-8", "R0-16", "R0-0"), (0,), [kd.ONTEST_SEQ]
        y7, y6 = [s7.cikti_yolu(True)], [os.path.join(kd.EV, "gelistirme.jsonl")]
    else:
        kosullar, toh, seqs = tuple(KOSUL), TOHUMLAR, list(kd.SEQS)
        y7, y6 = [s7.cikti_yolu(False, p) for p in ("I", "II", "III")], [kd.cikti_yolu(False, p) for p in ("I", "II", "III")]
    D, s_b = ka.oku_basvuru(kd.ANA, toh, sonda=a.ontest)
    s_ag = oku_ag(kd.ANA, toh, D)
    s_6 = {}
    if not a.ontest:
        D, s_6 = ka.oku_dogrula(y6, toh, D=D)
        D = {k: r for k, r in D.items() if not k[0].endswith("|A0y")}       # Sapma 6'nın A0y kayıtları kullanılmaz
    s_7 = oku_s7(y7, toh, D)
    b = butunluk(D, kosullar, toh, seqs, pencere, fp, ontest=a.ontest)
    V = dogrulamalar(D, kosullar)
    print(f"kod parmak izi {fp}; okuma: başvuru {s_b}; AG {s_ag}; Sapma 6 {s_6}; Sapma 7 {s_7}")
    print(f"bütünlük: {b['var']}/{b['beklenen_kayit']}; eksik {b['eksik']} {b['eksik_ornek']}; tam = {b['tam']}")
    for k, v in b["tablo"].items():
        print(f"   {k:12s} {v['tur']:5s} {v['var']:4d}/{v['beklenen']:4d} durum {v['durum']} boş {v['bos_alanli_ok']} "
              f"fp {'tamam' if v['fp_tamam'] else v['code_fp']}")
    print(f"doğrulamalar: V1 A0y ≡ A0 {V['V1']}; V2 ön-uç {V['V2']}; V3 AGg kapı {V['V3']}; V5 bayraklar {V['V5']}; "
          f"V6 {V['V6']}; tamam = {V['tamam']}")
    if a.butunluk:
        return 0 if (b["tam"] and V["tamam"]) else 3
    if a.ontest or not (b["tam"] and V["tamam"]):
        print("DUR: tablo yalnız tam ve doğrulanmış ızgara verisiyle üretilir")
        return 3
    S = calistir(D)
    ozet = {os.path.basename(p): hashlib.sha256(open(p, "rb").read()).hexdigest() for p in y7 + y6}
    json.dump(dict(kaynak=dict(code_fp=fp, kayit_sha256=ozet, tohumlar=list(toh), okuma=dict(basvuru=s_b, ag=s_ag, s6=s_6, s7=s_7)),
                   butunluk=b, dogrulamalar=V, esikler=ESIK, **S), open(a.cikti, "w"), indent=1, ensure_ascii=False, default=float)
    yazdir(S)
    print(f"-> {a.cikti}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
