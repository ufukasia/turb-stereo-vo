#!/usr/bin/env python
"""Uyarlanır izleme kapısı kolunun doğrulayıcı ızgarasının analizi (ON_KAYIT_FAZ2.md §5 Sapma 6).

Bu betik ızgaranın hiçbir sonucu görülmeden yazıldı. İstatistik ana ızgaranınkiyle aynıdır (turb/faz2_analiz.py):
birim (sekans, tohum); birim değeri = iki tarafta da geçerli pencerelerin ortalaması (ATE aritmetik, ANEES geometrik);
etki = eşli medyan oran, birim bootstrap 10 000 ile %95 GA; iki yönlü Wilcoxon işaretli sıra (log oranlar); Holm.

    python kapi_dogrula_analiz.py                 # veri tam değilse yalnız bütünlük ve doğrulama sayımları
    python kapi_dogrula_analiz.py --butunluk      # koşum sürerken: yalnız sayımlar (metrik tablosu basılmaz)
    python kapi_dogrula_analiz.py --oz_test       # sentetik veri + ana ızgara sayılarını yeniden üretme denetimi
    python kapi_dogrula_analiz.py --ontest --butunluk     # ön test çıktısı (tohum 0, V1_01_easy): okuma denetimi

Koşullar (rejim, D/r0): R0-8, R0-16 (tau0 = dt), R5-8, R5-4 (v = 5 m/s), R0-0 (temiz).
Etiketler (rejim|kol) ve kaynakları:
  A0, A1   eski kod: ana ızgara (R0; evidence/FAZ2) ve hızlı rejim ızgarası (R5; evidence/HIZLI, Sapma 5)
  A3S      çevrimiçi skaler kovaryans, kapı kapalı: R0-8 ve R0-16'da A3S ızgarası (eski kod; Sapma 4), R5 ve temizde
           doğrulayıcı ızgara (yeni kod, kapı kapalı: eski kodla bit-aynılığı ön testte gösterildi)
  KAPI     A3S + uyarlanır izleme kapısı (dondurulmuş aday; evidence/KAPI/SECIM.json)
  A1g      oracle skaler kovaryans + kapı (üst sınır)
  A0y      yeni kodla A0 (doğrulama: A0 kaydıyla bit-aynı olmalı; teste girmez)

Aileler (Holm aile içinde):
  B   (birincil, HK1-HK3)  KAPI / A0, ATE; R0-8, R5-8, R5-4                       m = 3
  T   (HK4, temiz)         KAPI / A0, ATE; R0-0: %95 GA üst sınırı < 1.10 (eşdeğer-altı olmama)
  I1  KAPI / A0, ATE; R0-16                                                        m = 1
  I2  kapı etkisi: KAPI / A3S (çevrimiçi) ve A1g / A1 (oracle), ATE; dört koşul    m = 8
  I3  çevrimiçi ↔ oracle: KAPI / A1g, ATE; dört koşul                              m = 4
  I4  A1g / A0, A1 / A0, A3S / A0, ATE; dört koşul                                 m = 12
  I5  ANEES_poz (ve ayrı aile olarak ANEES_don): KAPI / A0, A1g / A0, A3S / A0, A1 / A0, KAPI / A3S; dört koşul   m = 20
  T2  temiz, betimsel: A3S / A0 ve KAPI / A3S (ATE), KAPI / A0 ve A3S / A0 (ANEES_poz)
Betimsel: kol × koşul birim medyanları, başarısız pencere sayıları, kapı sigması ve sigma_hat dağılımı.
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

EV = os.path.join(ROOT, "evidence", "KAPI")
FP_ANA = "08e28b83407c"
BIT = ("m_ate_se3", "n_nees_pos_mean", "n_nees_rot_mean", "m_rpe_trans")
FE_SURE_ALANI = "total_execution_time_seconds"
TOHUMLAR = (1, 2, 3, 4)
KOSUL = {"R0-8": ("R0", 8.0, "II"), "R5-8": ("R5", 8.0, "I"), "R5-4": ("R5", 4.0, "I"),
         "R0-16": ("R0", 16.0, "II"), "R0-0": ("R0", 0.0, "III")}
BIRINCIL = ("R0-8", "R5-8", "R5-4")
TURB = ("R0-8", "R5-8", "R5-4", "R0-16")
YENI = {"I": ("KAPI", "A3S", "A1g"), "II": ("A0y", "KAPI", "A1g"), "III": ("A0y", "A3S", "KAPI")}
BASVURU = {"I": ("A0", "A1"), "II": ("A0", "A1", "A3S"), "III": ("A0",)}
A0Y_TOHUM_I = (1,)
KAPILI = ("KAPI", "A1g")
# Sapma 6'da yazılı sabitler
ESIK = dict(alfa=0.05, beklenen_oran=0.85, hk4_ust_sinir=1.10)
D_AP, DT, V_R5 = 0.10, 0.05, 5.0


def tau0_over_dt(rejim, dr):
    return 1.0 if rejim == "R0" else 0.314 * (D_AP / float(dr)) / (V_R5 * DT)


# ----------------------------------------------------------------------------- okuma
def _satirlar(p, sayim):
    if not os.path.exists(p):
        sayim["dosya_yok"] += 1
        return
    with open(p) as f:
        for line in f:
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                sayim["bozuk_satir"] += 1


def _koy(D, sayim, lab, dr, s, t, w, r):
    k = (lab, float(dr), s, int(t), int(w))
    if k in D:
        sayim["yinelenen"] += 1
    D[k] = r
    sayim["kayit"] += 1


def oku_basvuru(ana, tohumlar, D=None, sonda=False):
    """Eski kodun kayıtları -> D[(rejim|kol, dr, sekans, tohum, pencere)]."""
    D = {} if D is None else D
    S = {k: Counter() for k in ("ana", "a3s", "hizli", "sonda")}
    for r in _satirlar(os.path.join(ana, "FAZ2", "records.jsonl"), S["ana"]):
        c = r.get("cell") or {}
        if ("arm" in r and c.get("tag") == "faz2" and r["arm"] in ("A0", "A1") and int(c["seed"]) in tohumlar
                and float(c["D_over_r0"]) in (0.0, 8.0, 16.0) and abs(float(c.get("tau0_over_dt", 1.0)) - 1.0) < 1e-12
                and r.get("code_fp") == FP_ANA):
            _koy(D, S["ana"], f"R0|{r['arm']}", c["D_over_r0"], c["sequence"], c["seed"], r["window"], r)
    for r in _satirlar(os.path.join(ana, "A3S", "records.jsonl"), S["a3s"]):
        c = r.get("cell") or {}
        if ("arm" in r and c.get("tag") == "a3s" and r["arm"] == "A0C+A3" and int(c["seed"]) in tohumlar
                and tuple(c.get("backend_argv") or ()) == ("--turb_cal_factor", "1.0")
                and float(c["D_over_r0"]) in (8.0, 16.0) and r.get("code_fp") == FP_ANA):
            _koy(D, S["a3s"], "R0|A3S", c["D_over_r0"], c["sequence"], c["seed"], r["window"], r)
    for r in _satirlar(os.path.join(ana, "HIZLI", "records.jsonl"), S["hizli"]):
        c = r.get("cell") or {}
        if ("arm" in r and c.get("tag") == "hizli_R5" and r["arm"] in ("A0", "A1") and int(c["seed"]) in tohumlar
                and float(c["D_over_r0"]) in (8.0, 4.0) and r.get("code_fp") == FP_ANA
                and abs(float(c.get("tau0_over_dt", 1.0)) - tau0_over_dt("R5", c["D_over_r0"])) < 1e-9):
            _koy(D, S["hizli"], f"R5|{r['arm']}", c["D_over_r0"], c["sequence"], c["seed"], r["window"], r)
    if sonda:                                        # yalnız ön test (tohum 0): R5 başvurusu keşif sondasından
        for r in _satirlar(os.path.join(ana, "SONDA_KESIF", "records.jsonl"), S["sonda"]):
            c = r.get("cell") or {}
            if ("arm" in r and c.get("tag") == "sonda_R5_fe0" and r["arm"] in ("A0", "A1")
                    and int(c["seed"]) in tohumlar and r.get("code_fp") == FP_ANA):
                _koy(D, S["sonda"], f"R5|{r['arm']}", c["D_over_r0"], c["sequence"], c["seed"], r["window"], r)
    return D, {k: dict(v) for k, v in S.items()}


def oku_dogrula(yollar, tohumlar, D=None):
    """Doğrulayıcı ızgaranın kayıtları (yeni kod). Aynı (kol, hücre, pencere) yinelenirse son kayıt geçerlidir."""
    D = {} if D is None else D
    S = Counter()
    for p in yollar:
        for r in _satirlar(p, S):
            k = r.get("kapid") or {}
            if "varyant" not in k:
                S["on_uc_hatasi" if "kapid_fe" in r else "yabanci"] += 1
                continue
            if "arm" not in r:
                S["hucre_hatasi"] += 1
                continue
            if int(k["tohum"]) not in tohumlar:
                S["suzulen"] += 1
                continue
            v = "A0y" if k["varyant"] == "A0" else k["varyant"]
            _koy(D, S, f"{k['rejim']}|{v}", k["dr"], k["sekans"], k["tohum"], r["window"], r)
    return D, dict(S)


# ----------------------------------------------------------------------------- bütünlük ve doğrulama
def _kollar(parca, tohumlar, a0y_tohum_I):
    out = [(k, "yeni", tuple(tohumlar)) for k in YENI[parca]] + [(k, "basvuru", tuple(tohumlar)) for k in BASVURU[parca]]
    if parca == "I":
        out.append(("A0y", "yeni", tuple(t for t in tohumlar if t in a0y_tohum_I)))
    return out


def butunluk(D, kosullar, tohumlar, seqs, pencere, fp_yeni, a0y_tohum_I=A0Y_TOHUM_I):
    """Beklenen (koşul, kol, sekans, tohum, pencere) kayıtlarının varlığı; durum, boş alan ve kod parmak izi sayımı."""
    tablo, eksik = {}, []
    for ad in kosullar:
        rej, dr, parca = KOSUL[ad]
        for kol, tur, toh in _kollar(parca, tohumlar, a0y_tohum_I):
            durum, fp, bos, n_bek, n_var = Counter(), Counter(), 0, 0, 0
            for s in seqs:
                for t in toh:
                    for w in pencere[s]:
                        n_bek += 1
                        r = D.get((f"{rej}|{kol}", dr, s, t, int(w)))
                        if r is None:
                            eksik.append((ad, kol, s, t, int(w)))
                            continue
                        n_var += 1
                        durum[r.get("status")] += 1
                        fp[r.get("code_fp")] += 1
                        bos += r.get("status") == "ok" and any(r.get(x) is None for x in BIT)
            fp_ok = set(fp) <= {fp_yeni if tur == "yeni" else FP_ANA}
            tablo[f"{ad}|{kol}"] = dict(tur=tur, beklenen=n_bek, var=n_var, durum=dict(durum), bos_alanli_ok=int(bos),
                                        code_fp=dict(fp), fp_tamam=bool(fp_ok))
    n_bek = sum(v["beklenen"] for v in tablo.values())
    tam = bool(n_bek > 0 and not eksik and all(v["bos_alanli_ok"] == 0 and v["fp_tamam"] for v in tablo.values()))
    return dict(tablo=tablo, beklenen_kayit=n_bek, var=n_bek - len(eksik), eksik=len(eksik), eksik_ornek=eksik[:6],
                tam=tam)


def _fe_ayni(a, b):
    if not a or not b:
        return None
    return ({k: v for k, v in a.items() if k != FE_SURE_ALANI} == {k: v for k, v in b.items() if k != FE_SURE_ALANI})


def dogrulamalar(D, kosullar, aday, c_deger, s_oracle, taban):
    """Koşumun kendi denetimlerinin yeniden hesabı (metrik düzeyi raporlanmaz):
      V1  A0y kayıtları A0 (eski kod) kayıtlarıyla dört alanda bit-aynı
      V2  yeniden üretilen ön-ucun istatistikleri başvuru A0 kaydınınkiyle aynı (süre alanı hariç)
      V3  kapılı kollarda gate = 3; A1g kapı sigması = oracle bireysel sigma; KAPI en küçük kapı sigması = taban
      V4  kapısız yeni-kod kayıtlarında (A0y, A3S) kapı kapalı
      V5  kol bayrakları ve seçilen aday kayıtlarda beklenen gibi"""
    V = {k: Counter() for k in ("V1", "V2", "V3", "V4", "V5")}
    for (lab, dr, s, t, w), r in D.items():
        rej, kol = lab.split("|")
        if not any(KOSUL[a][0] == rej and KOSUL[a][1] == dr for a in kosullar):
            continue
        k = r.get("kapid")
        if k is None:
            continue                                   # başvuru kaydı
        ref = D.get((f"{rej}|A0", dr, s, t, w))
        if kol == "A0y":
            if ref is None:
                V["V1"]["basvuru_yok"] += 1
            else:
                V["V1"]["ayni" if all(r.get(x) == ref.get(x) for x in BIT) else "FARKLI"] += 1
        if k.get("kaynak") == "fe" and kol == "A0y":
            a = _fe_ayni(r.get("frontend_stats"), (ref or {}).get("frontend_stats"))
            V["V2"]["basvuru_yok" if a is None else ("ayni" if a else "FARKLI")] += 1
        st = r.get("arm_stats") or {}
        if kol in KAPILI:
            iyi = st.get("gate") == 3
            mn, mx = st.get("gate_sigma_min"), st.get("gate_sigma_max")
            if kol == "A1g":
                iyi = iyi and mn is not None and abs(mn - s_oracle(dr)) < 1e-9 and abs(mx - s_oracle(dr)) < 1e-9
            else:
                iyi = iyi and mn is not None and abs(mn - taban) < 1e-9
            V["V3"]["tamam" if iyi else "HATALI"] += 1
        else:
            V["V4"]["tamam" if not st.get("gate") else "HATALI"] += 1
        argv = list(k.get("argv") or [])
        bek = {"A0y": ["--turb_cal_factor", "12.0"], "A3S": ["--turb_cal_factor", "1.0"],
               "KAPI": ["--turb_cal_factor", "1.0", "--turb_dimm_sigma_scale", repr(float(c_deger)), "--turb_gate", "3"],
               "A1g": ["--turb_cal_factor", "12.0", "--turb_gate", "3"]}[kol]
        V["V5"]["tamam" if (argv == bek and k.get("secim") == aday) else "HATALI"] += 1
    kotu = sum(v.get("FARKLI", 0) + v.get("HATALI", 0) for v in V.values())
    return dict({k: dict(v) for k, v in V.items()}, tamam=bool(kotu == 0 and V["V1"].get("basvuru_yok", 0) == 0))


# ----------------------------------------------------------------------------- testler
def test(D, pay, payda, kosul, metrik):
    rej, dr, _ = KOSUL[kosul]
    a, b, U, dus = eslesik(D, f"{rej}|{pay}", f"{rej}|{payda}", dr, metrik)
    return dict(oran_testi(a, b), dusen_birim=dus) if len(a) else None


def aile(D, ciftler, kosullar, metrik):
    out = {}
    for pay, payda in ciftler:
        for ad in kosullar:
            t = test(D, pay, payda, ad, metrik)
            if t is not None:
                out[f"{pay}/{payda}@{ad}"] = t
    ph = holm({k: v["p"] for k, v in out.items()})
    for k in out:
        out[k]["p_holm"] = ph[k]
    return out


def betimsel(D, kosullar):
    out = {}
    for ad in kosullar:
        rej, dr, parca = KOSUL[ad]
        for kol in sorted(set(YENI[parca]) | set(BASVURU[parca])):
            if kol == "A0y":
                continue
            for metrik in ("ATE", "ANEES_poz", "ANEES_don"):
                f, tur = METRIK[metrik]
                g = defaultdict(list)
                for (k, d, s, t, w), r in D.items():
                    if k == f"{rej}|{kol}" and d == dr and f(r) is not None:
                        g[(s, t)].append(f(r))
                if g:
                    v = np.array([_ort(x, tur) for x in g.values()])
                    out.setdefault(ad, {}).setdefault(kol, {})[metrik] = dict(
                        medyan=float(np.median(v)), q25=float(np.percentile(v, 25)), q75=float(np.percentile(v, 75)),
                        n_birim=len(g))
    return out


def basarisizlik(D, kosullar):
    out = {}
    for ad in kosullar:
        rej, dr, parca = KOSUL[ad]
        for kol in sorted(set(YENI[parca]) | set(BASVURU[parca])):
            R = [r for (k, d, s, t, w), r in D.items() if k == f"{rej}|{kol}" and d == dr]
            if R and kol != "A0y":
                kotu = sum(0 if (r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged")) else 1
                           for r in R)
                out[f"{kol}@{ad}"] = dict(n=len(R), basarisiz=kotu)
    return out


def kapi_sigmasi(D, kosullar, s_oracle):
    """Birim başına (ilk pencere kaydının arm_stats'ı; koşum başına tek değer) kapı sigması ve sigma_hat medyanı."""
    out = {}
    for ad in kosullar:
        rej, dr, _ = KOSUL[ad]
        g = {}
        for (k, d, s, t, w), r in sorted(D.items(), key=lambda x: x[0]):
            if k == f"{rej}|KAPI" and d == dr and (s, t) not in g:
                g[(s, t)] = r.get("arm_stats") or {}
        if g:
            med = np.array([x.get("gate_sigma_median", np.nan) for x in g.values()], float)
            mx = np.array([x.get("gate_sigma_max", np.nan) for x in g.values()], float)
            sh = np.array([x.get("sigma_hat_median", np.nan) for x in g.values()], float)
            out[ad] = dict(n_birim=len(g), kapi_sigma_medyan=float(np.nanmedian(med)),
                           kapi_sigma_medyan_aralik=[float(np.nanmin(med)), float(np.nanmax(med))],
                           kapi_sigma_maks=float(np.nanmax(mx)), sigma_hat_medyan=float(np.nanmedian(sh)),
                           oracle_sigma=float(s_oracle(dr)))
    return out


def kararlar(S):
    K, a = {}, ESIK["alfa"]
    for ad, hk in zip(BIRINCIL, ("HK1", "HK2", "HK3")):
        v = S["B"].get(f"KAPI/A0@{ad}")
        if v:
            anl = v["p_holm"] < a
            K[hk] = dict(kosul=ad, oran=v["medyan_oran"], ga95=v["ga95"], p=v["p"], p_holm=v["p_holm"], n=v["n"],
                         azalan_birim=v["yon_azalma"], beklenen_buyukluk=bool(v["medyan_oran"] <= ESIK["beklenen_oran"]),
                         karar=("TUTTU (ATE anlamlı azalıyor)" if anl and v["medyan_oran"] < 1 else
                                "TERS YÖN (ATE anlamlı artıyor)" if anl and v["medyan_oran"] > 1 else
                                "TUTMADI (fark anlamlı değil)"))
    v = S["T"].get("KAPI/A0@R0-0")
    if v:
        K["HK4"] = dict(oran=v["medyan_oran"], ga95=v["ga95"], p=v["p"], n=v["n"], ust_sinir=ESIK["hk4_ust_sinir"],
                        karar=("GEÇTİ (GA üst sınırı < %.2f)" % ESIK["hk4_ust_sinir"]
                               if v["ga95"][1] < ESIK["hk4_ust_sinir"] else "KALDI (GA üst sınırı >= %.2f)"
                               % ESIK["hk4_ust_sinir"]))
    return K


def calistir(D, s_oracle):
    S = dict(B=aile(D, [("KAPI", "A0")], BIRINCIL, "ATE"),
             T=aile(D, [("KAPI", "A0")], ("R0-0",), "ATE"),
             I1=aile(D, [("KAPI", "A0")], ("R0-16",), "ATE"),
             I2=aile(D, [("KAPI", "A3S"), ("A1g", "A1")], TURB, "ATE"),
             I3=aile(D, [("KAPI", "A1g")], TURB, "ATE"),
             I4=aile(D, [("A1g", "A0"), ("A1", "A0"), ("A3S", "A0")], TURB, "ATE"),
             I5_poz=aile(D, [("KAPI", "A0"), ("A1g", "A0"), ("A3S", "A0"), ("A1", "A0"), ("KAPI", "A3S")], TURB,
                         "ANEES_poz"),
             I5_don=aile(D, [("KAPI", "A0"), ("A1g", "A0"), ("A3S", "A0"), ("A1", "A0"), ("KAPI", "A3S")], TURB,
                         "ANEES_don"),
             T2_ATE=aile(D, [("A3S", "A0"), ("KAPI", "A3S")], ("R0-0",), "ATE"),
             T2_ANEES_poz=aile(D, [("KAPI", "A0"), ("A3S", "A0")], ("R0-0",), "ANEES_poz"),
             betimsel=betimsel(D, tuple(KOSUL)), basarisizlik=basarisizlik(D, tuple(KOSUL)),
             kapi_sigmasi=kapi_sigmasi(D, tuple(KOSUL), s_oracle))
    S["kararlar"] = kararlar(S)
    return S


def _satir(k, v):
    return (f"   {k:16s} oran {v['medyan_oran']:.3f} [{v['ga95'][0]:.3f}, {v['ga95'][1]:.3f}] n={v['n']} "
            f"azalan {v['yon_azalma']}/{v['n']} p={v['p']:.3g} p_holm={v['p_holm']:.3g} düşen={v.get('dusen_birim', 0)}")


def yazdir(S):
    for ad, bas in (("B", "B (birincil; HK1-HK3): KAPI / A0, ATE"), ("T", "T (HK4; temiz): KAPI / A0, ATE"),
                    ("I1", "I1: KAPI / A0, ATE, R0-16"), ("I2", "I2 kapı etkisi: KAPI / A3S, A1g / A1, ATE"),
                    ("I3", "I3 çevrimiçi / oracle: KAPI / A1g, ATE"), ("I4", "I4: A1g, A1, A3S / A0, ATE"),
                    ("I5_poz", "I5: ANEES_poz"), ("I5_don", "I5: ANEES_don"), ("T2_ATE", "T2 temiz, ATE"),
                    ("T2_ANEES_poz", "T2 temiz, ANEES_poz")):
        print(bas)
        for k, v in S[ad].items():
            print(_satir(k, v))
    print("birim medyanları: ATE [m] / ANEES_poz")
    for ad, kol in S["betimsel"].items():
        print("   " + ad + ": " + "; ".join(f"{k} {v['ATE']['medyan']:.4f} / {v.get('ANEES_poz', {}).get('medyan', float('nan')):.1f}"
                                            for k, v in kol.items() if "ATE" in v))
    print("başarısız pencere (geçersiz ya da ıraksamış)")
    print("   " + "; ".join(f"{k} {v['basarisiz']}/{v['n']}" for k, v in S["basarisizlik"].items()))
    print("kapı sigması (KAPI): medyan [birimlerde en küçük, en büyük] / en büyük maks / sigma_hat medyanı / oracle")
    for ad, v in S["kapi_sigmasi"].items():
        print(f"   {ad}: {v['kapi_sigma_medyan']:.3f} [{v['kapi_sigma_medyan_aralik'][0]:.3f}, "
              f"{v['kapi_sigma_medyan_aralik'][1]:.3f}] / {v['kapi_sigma_maks']:.3f} / {v['sigma_hat_medyan']:.3f} / "
              f"{v['oracle_sigma']:.3f}")
    print("KARARLAR")
    for hk, v in S["kararlar"].items():
        if hk == "HK4":
            print(f"   HK4 (temiz): {v['oran']:.3f} [{v['ga95'][0]:.3f}, {v['ga95'][1]:.3f}] n={v['n']} -> {v['karar']}")
        else:
            print(f"   {hk} ({v['kosul']}): {v['oran']:.3f} [{v['ga95'][0]:.3f}, {v['ga95'][1]:.3f}] n={v['n']} "
                  f"p_holm={v['p_holm']:.3g} -> {v['karar']}; beklenen büyüklük (<= {ESIK['beklenen_oran']}): "
                  f"{'evet' if v['beklenen_buyukluk'] else 'hayır'}")


# ----------------------------------------------------------------------------- öz-test
def oz_test(ana):
    rng = np.random.default_rng(6)
    seqs, toh, pencere = ("S1", "S2", "S3"), (1, 2, 3, 4), {"S1": (10, 20), "S2": (10, 20), "S3": (10, 20)}
    FPY = "yeni"
    carp = {"R0-8": dict(KAPI=0.7, A1g=0.5, A3S=1.0, A1=1.1), "R5-8": dict(KAPI=0.8, A1g=0.6, A3S=1.05, A1=1.2),
            "R5-4": dict(KAPI=0.9, A1g=0.85, A3S=1.0, A1=1.0), "R0-16": dict(KAPI=0.6, A1g=0.4, A3S=0.9, A1=1.0),
            "R0-0": dict(KAPI=1.02, A3S=1.01)}
    so = lambda dr: math.sqrt(1.0 + (0.5 * dr) ** 2)             # noqa: E731
    D = {}
    for ad, (rej, dr, parca) in KOSUL.items():
        for s in seqs:
            for t in toh:
                for w in pencere[s]:
                    ba = 0.02 * (1 + dr / 4) * math.exp(0.3 * rng.standard_normal())
                    bn = 150 * (1 + dr) * math.exp(0.3 * rng.standard_normal())
                    for kol, tur, tt in _kollar(parca, toh, A0Y_TOHUM_I):
                        if t not in tt:
                            continue
                        f = 1.0 if kol in ("A0", "A0y") else carp[ad][kol]
                        r = dict(status="ok", m_valid=True, n_valid=True, m_diverged=False,
                                 code_fp=FPY if tur == "yeni" else FP_ANA, m_ate_se3=ba * f, m_rpe_trans=0.01,
                                 n_nees_pos_mean=bn * f * f, n_nees_rot_mean=bn * f, frontend_stats={"a": 1, FE_SURE_ALANI: 5})
                        if tur == "yeni":
                            argv = {"A0y": ["--turb_cal_factor", "12.0"], "A3S": ["--turb_cal_factor", "1.0"],
                                    "KAPI": ["--turb_cal_factor", "1.0", "--turb_dimm_sigma_scale", "1.0", "--turb_gate",
                                             "3"], "A1g": ["--turb_cal_factor", "12.0", "--turb_gate", "3"]}[kol]
                            r["kapid"] = dict(secim="D1", argv=argv, kaynak="fe")
                            r["frontend_stats"] = {"a": 1, FE_SURE_ALANI: 7}
                            if kol in KAPILI:
                                g = so(dr) if kol == "A1g" else 1.5
                                r["arm_stats"] = dict(gate=3, gate_sigma_min=so(dr) if kol == "A1g" else 1.25,
                                                      gate_sigma_max=g, gate_sigma_median=g, sigma_hat_median=1.4)
                        D[(f"{rej}|{kol}", dr, s, t, w)] = r
    b = butunluk(D, tuple(KOSUL), toh, seqs, pencere, FPY)
    n_bek = 24 * (3 * 2 + 2 * 2) + 6 * 2 + 24 * 6 * 2 + 24 * 4          # I: (3 yeni + 2 başvuru) × 2 düzey + A0y tohum 1
    assert b["tam"] and b["beklenen_kayit"] == n_bek == b["var"], (b["beklenen_kayit"], n_bek, b["eksik_ornek"])
    V = dogrulamalar(D, tuple(KOSUL), "D1", 1.0, so, 1.25)
    assert V["tamam"] and V["V1"] == {"ayni": 24 * 3 + 6 * 2} and V["V2"] == {"ayni": 24 * 3 + 6 * 2}, V
    assert V["V3"] == {"tamam": 24 * 2 * 4 + 24} and "HATALI" not in V["V5"], V
    S = calistir(D, so)
    for ad in TURB:
        assert abs(S["B" if ad in BIRINCIL else "I1"][f"KAPI/A0@{ad}"]["medyan_oran"] - carp[ad]["KAPI"]) < 1e-12
        assert abs(S["I2"][f"KAPI/A3S@{ad}"]["medyan_oran"] - carp[ad]["KAPI"] / carp[ad]["A3S"]) < 1e-12
        assert abs(S["I2"][f"A1g/A1@{ad}"]["medyan_oran"] - carp[ad]["A1g"] / carp[ad]["A1"]) < 1e-12
        assert abs(S["I3"][f"KAPI/A1g@{ad}"]["medyan_oran"] - carp[ad]["KAPI"] / carp[ad]["A1g"]) < 1e-12
        assert abs(S["I5_poz"][f"KAPI/A0@{ad}"]["medyan_oran"] - carp[ad]["KAPI"] ** 2) < 1e-12
        assert abs(S["I5_don"][f"A1g/A0@{ad}"]["medyan_oran"] - carp[ad]["A1g"]) < 1e-12
    assert (len(S["B"]), len(S["I1"]), len(S["I2"]), len(S["I3"]), len(S["I4"]), len(S["I5_poz"])) == (3, 1, 8, 4, 12, 20)
    assert all(v["n"] == 12 and v["dusen_birim"] == 0 for v in S["B"].values())
    K = S["kararlar"]
    assert [K[h]["karar"][:5] for h in ("HK1", "HK2", "HK3")] == ["TUTTU"] * 3, K
    assert K["HK1"]["beklenen_buyukluk"] and K["HK2"]["beklenen_buyukluk"] and not K["HK3"]["beklenen_buyukluk"]
    assert K["HK4"]["karar"].startswith("GEÇTİ") and abs(K["HK4"]["oran"] - 1.02) < 1e-12
    # Holm: aynı p üç kez -> p_holm = 3p (m = 3); tek üyeli ailede p_holm = p
    p = S["B"]["KAPI/A0@R0-8"]["p"]
    assert abs(S["B"]["KAPI/A0@R0-8"]["p_holm"] - min(1.0, 3 * p)) < 1e-15 and S["I1"]["KAPI/A0@R0-16"]["p_holm"] == \
        S["I1"]["KAPI/A0@R0-16"]["p"]
    assert abs(S["kapi_sigmasi"]["R0-8"]["kapi_sigma_medyan"] - 1.5) < 1e-12 and S["kapi_sigmasi"]["R0-8"]["n_birim"] == 12
    # HK4 başarısızlığı, etkisizlik ve ters yön
    D2 = {k: (dict(r, m_ate_se3=r["m_ate_se3"] * 1.2) if k[0] == "R0|KAPI" and k[1] == 0.0 else r) for k, r in D.items()}
    assert calistir(D2, so)["kararlar"]["HK4"]["karar"].startswith("KALDI")
    D3 = {k: (dict(D[("R0|A0",) + k[1:]]) if k[0] == "R0|KAPI" and k[1] == 8.0 else r) for k, r in D.items()}
    k3 = calistir(D3, so)["kararlar"]["HK1"]
    assert k3["karar"].startswith("TUTMADI") and k3["oran"] == 1.0 and k3["p"] == 1.0, k3
    D4 = {k: (dict(r, m_ate_se3=r["m_ate_se3"] * 2.0) if k[0] == "R5|KAPI" and k[1] == 8.0 else r) for k, r in D.items()}
    assert calistir(D4, so)["kararlar"]["HK2"]["karar"].startswith("TERS")
    # eksik kayıt, yanlış parmak izi, A0y farkı, kapı istatistiği ve bayrak hatası yakalanmalı
    D5 = dict(D)
    del D5[("R5|KAPI", 4.0, "S2", 3, 20)]
    b5 = butunluk(D5, tuple(KOSUL), toh, seqs, pencere, FPY)
    assert not b5["tam"] and b5["eksik"] == 1
    assert not butunluk(D, tuple(KOSUL), toh, seqs, pencere, "baska")["tam"]
    k0 = ("R0|A0y", 8.0, "S1", 1, 10)
    assert not dogrulamalar({**D, k0: dict(D[k0], m_rpe_trans=0.011)}, tuple(KOSUL), "D1", 1.0, so, 1.25)["tamam"]
    k1 = ("R5|A1g", 8.0, "S1", 2, 10)
    assert not dogrulamalar({**D, k1: dict(D[k1], arm_stats=dict(D[k1]["arm_stats"], gate_sigma_max=9.0))},
                            tuple(KOSUL), "D1", 1.0, so, 1.25)["tamam"]
    assert not dogrulamalar(D, tuple(KOSUL), "D2", 1.0, so, 1.25)["tamam"]
    k2 = ("R0|A0y", 0.0, "S3", 4, 20)
    assert not dogrulamalar({**D, k2: dict(D[k2], frontend_stats={"a": 2})}, tuple(KOSUL), "D1", 1.0, so,
                            1.25)["tamam"]
    assert abs(tau0_over_dt("R5", 8) - 0.0157) < 1e-12 and abs(tau0_over_dt("R5", 4) - 0.0314) < 1e-12
    print("öz-test (sentetik) GEÇTİ")
    # ana ızgara sayılarını yeniden üretme: bu betiğin okuyucusu + eşleme ile A1 / A0 (ATE, ANEES_poz), D/r0 8 ve 16,
    # 45 birim -> evidence/FAZ2/analiz.json ile sayı sayı aynı
    ap = os.path.join(ana, "FAZ2", "analiz.json")
    if os.path.exists(ap) and os.path.exists(os.path.join(ana, "FAZ2", "records.jsonl")):
        Dm, _ = oku_basvuru(ana, (0, 1, 2, 3, 4))
        ref = json.load(open(ap))
        n = 0
        for metrik in ("ATE", "ANEES_poz"):
            for ad in ("R0-8", "R0-16"):
                v, r = test(Dm, "A1", "A0", ad, metrik), ref[f"aile_{metrik}"][f"A1@{KOSUL[ad][1]:g}"]
                assert v["n"] == r["n"] == 45 and v["medyan_oran"] == r["medyan_oran"] and v["p"] == r["p"] \
                    and v["ga95"] == r["ga95"], (metrik, ad, v, r)
                n += 1
        v = test(Dm, "A0", "A0", "R0-8", "ATE")
        assert v["medyan_oran"] == 1.0 and v["p"] == 1.0
        print(f"öz-test (ana ızgara): A1 / A0 için {n} girdi analiz.json ile sayı sayı aynı; kendisiyle oran 1.000 GEÇTİ")
    else:
        print("öz-test (ana ızgara): kayıt ya da analiz.json yok, atlandı")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--butunluk", action="store_true", help="yalnız bütünlük ve doğrulama sayımları")
    ap.add_argument("--ontest", action="store_true", help="ön test çıktısı (tohum 0, V1_01_easy; R5-4 yok)")
    ap.add_argument("--oz_test", action="store_true")
    ap.add_argument("--cikti", default=os.path.join(EV, "dogrula_analiz.json"))
    a = ap.parse_args()
    import kapi_dogrula as kd
    if a.oz_test:
        oz_test(kd.ANA)
        return 0
    aday, gecici = kd.secim_oku()
    if aday is None:
        print("DUR: evidence/KAPI/SECIM.json yok")
        return 3
    fp_yeni = kd.code_fingerprint()
    c = kd.kk.C_DEGER[aday]
    taban = math.sqrt(1.0 + (c * kd.kk.SIGMA_TABAN) ** 2)
    pencere = {s: tuple(int(w) for w in kd.windows_for(s, kd.FAZ2_WINDOWS_FILE)) for s in kd.SEQS}
    if a.ontest:
        kosullar, toh, seqs, a0y = ("R0-8", "R5-8", "R0-16", "R0-0"), (0,), [kd.ONTEST_SEQ], (0,)
        yollar = [kd.cikti_yolu(True)]
    else:
        kosullar, toh, seqs, a0y = tuple(KOSUL), TOHUMLAR, list(kd.SEQS), A0Y_TOHUM_I
        yollar = [kd.cikti_yolu(False, p) for p in ("I", "II", "III")]
    D, s_b = oku_basvuru(kd.ANA, toh, sonda=a.ontest)
    D, s_d = oku_dogrula(yollar, toh, D=D)
    b = butunluk(D, kosullar, toh, seqs, pencere, fp_yeni, a0y_tohum_I=a0y)
    V = dogrulamalar(D, kosullar, aday, c, kd.s_oracle, taban)
    print(f"seçilen aday {aday} (c = {c!r}); yeni kod parmak izi {fp_yeni}; okuma: başvuru {s_b}; doğrulayıcı {s_d}")
    print(f"bütünlük: {b['var']}/{b['beklenen_kayit']} kayıt; eksik {b['eksik']} {b['eksik_ornek']}; tam = {b['tam']}")
    for k, v in b["tablo"].items():
        print(f"   {k:14s} {v['tur']:8s} {v['var']:4d}/{v['beklenen']:4d} durum {v['durum']} boş {v['bos_alanli_ok']} "
              f"fp {'tamam' if v['fp_tamam'] else v['code_fp']}")
    print(f"doğrulamalar: V1 A0y ≡ A0 {V['V1']}; V2 ön-uç istatistikleri {V['V2']}; V3 kapı istatistikleri {V['V3']}; "
          f"V4 kapısız kollar {V['V4']}; V5 bayraklar {V['V5']}; tamam = {V['tamam']}")
    if a.butunluk:
        return 0 if (b["tam"] and V["tamam"]) else 3
    if a.ontest:
        print("DUR: ön test çıktısından test tablosu basılmaz (1 birim); --butunluk kullanın")
        return 3
    if not (b["tam"] and V["tamam"]):
        print("DUR: veri tam değil ya da doğrulama tutmadı; karar tabloları yalnız tam ve doğrulanmış veriyle üretilir")
        return 3
    S = calistir(D, kd.s_oracle)
    ozet = {os.path.basename(p): hashlib.sha256(open(p, "rb").read()).hexdigest() for p in yollar}
    out = dict(kaynak=dict(aday=aday, c=c, secim_sha256=hashlib.sha256(open(kd.SECIM, "rb").read()).hexdigest(),
                           code_fp_yeni=fp_yeni, code_fp_ana=FP_ANA, kayit_sha256=ozet, tohumlar=list(toh),
                           kosullar={k: dict(rejim=v[0], D_over_r0=v[1], parca=v[2], tau0_over_dt=(
                               tau0_over_dt(v[0], v[1]) if v[1] > 0 else None)) for k, v in KOSUL.items()},
                           okuma=dict(basvuru=s_b, dogrulayici=s_d)),
               butunluk=b, dogrulamalar=V, esikler=ESIK, **S)
    json.dump(out, open(a.cikti, "w"), indent=1, ensure_ascii=False, default=float)
    yazdir(S)
    print(f"-> {a.cikti}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
