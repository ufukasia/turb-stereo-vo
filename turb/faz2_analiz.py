"""Faz 2 ana ızgara (H1, H3-EuRoC, H4, H5) + H6 analizi.
İşlemsel tanım: ON_TEST_FAZ2_EUROC.md §10 (bu betikten önce yazıldı); ön kayıt ON_KAYIT_FAZ2.md §2–3.

    python -m turb.faz2_analiz              # -> evidence/FAZ2/analiz.json + özet
    python -m turb.faz2_analiz --oz_test    # sentetik veriyle hesap denetimi

Birim (sekans, tohum); birim değeri = eşli karşılaştırmada iki tarafta da geçerli pencerelerin ortalaması
(ATE aritmetik, ANEES geometrik). Etki = medyan oran exp(medyan log(kol/ref)); %95 GA birim bootstrap 10 000
(tohum 0); iki yönlü Wilcoxon; Holm ailesi = kollar × D/r0, metrik başına.
"""
import argparse, json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DRS = (0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0)
KOLLAR = ("A1", "A3", "AG", "A0C", "A0C+A3")
B_BOOT = 10000


# ----------------------------------------------------------------------------- okuma ve metrikler
def oku(path):
    """-> D[(kol, dr, sekans, tohum, pencere)] = kayıt"""
    D = {}
    for l in open(path):
        r = json.loads(l)
        c = r["cell"]
        D[(r["arm"], float(c["D_over_r0"]), c["sequence"], int(c["seed"]), int(r["window"]))] = r
    return D


def ate(r):
    v = r.get("m_ate_se3")
    ok = r.get("status") == "ok" and r.get("m_valid") and v is not None and math.isfinite(v) and v > 0
    return float(v) if ok else None


def anees(r, k="pos"):
    v = r.get(f"n_nees_{k}_mean")
    ok = r.get("status") == "ok" and r.get("n_valid") and v is not None and math.isfinite(v) and v > 0
    return float(v) if ok else None


METRIK = {"ATE": (ate, "aritmetik"), "ANEES_poz": (lambda r: anees(r, "pos"), "geometrik"),
          "ANEES_don": (lambda r: anees(r, "rot"), "geometrik")}


def _ort(x, tur):
    x = np.asarray(x, float)
    return float(np.exp(np.mean(np.log(x)))) if tur == "geometrik" else float(np.mean(x))


def eslesik(D, kol, ref, dr, metrik, dr_ref=None):
    """Birim başına (kol, ref) değer çiftleri; yalnız iki tarafta da geçerli pencereler. -> (a, b, birimler, düşen)"""
    f, tur = METRIK[metrik]
    dr_ref = dr if dr_ref is None else dr_ref
    grup = defaultdict(lambda: ([], []))
    birim_hepsi = set()
    for (k, d, s, t, w), r in D.items():
        if k != kol or d != dr:
            continue
        birim_hepsi.add((s, t))
        r2 = D.get((ref, dr_ref, s, t, w))
        if r2 is None:
            continue
        va, vb = f(r), f(r2)
        if va is None or vb is None:
            continue
        grup[(s, t)][0].append(va)
        grup[(s, t)][1].append(vb)
    U = sorted(grup)
    a = np.array([_ort(grup[u][0], tur) for u in U])
    b = np.array([_ort(grup[u][1], tur) for u in U])
    return a, b, U, len(birim_hepsi - set(U))


def oran_testi(a, b, seed=0, B=B_BOOT):
    lr = np.log(a / b)
    n = len(lr)
    if n == 0:
        return dict(n=0)
    p = float(stats.wilcoxon(lr).pvalue) if np.any(lr != 0) else 1.0
    rng = np.random.default_rng(seed)
    bs = np.exp(np.median(lr[rng.integers(0, n, (B, n))], axis=1))
    return dict(n=n, medyan_oran=float(np.exp(np.median(lr))), ga95=[float(np.percentile(bs, 2.5)),
                float(np.percentile(bs, 97.5))], p=p, yon_azalma=int(np.sum(lr < 0)))


def holm(pler):
    """pler: {anahtar: p} -> {anahtar: p_holm}"""
    k = sorted(pler, key=lambda x: pler[x])
    m, out, onceki = len(k), {}, 0.0
    for i, key in enumerate(k):
        v = min(1.0, max(onceki, (m - i) * pler[key]))
        out[key] = v
        onceki = v
    return out


def bant(n):
    return float(stats.chi2.ppf(0.025, 3 * n) / n), float(stats.chi2.ppf(0.975, 3 * n) / n)


# ----------------------------------------------------------------------------- hipotezler
def aile(D, metrik):
    out = {}
    for kol in KOLLAR:
        for dr in DRS:
            a, b, U, dus = eslesik(D, kol, "A0", dr, metrik)
            if len(a):
                out[f"{kol}@{dr:g}"] = dict(oran_testi(a, b), dusen_birim=dus)
    ph = holm({k: v["p"] for k, v in out.items()})
    for k in out:
        out[k]["p_holm"] = ph[k]
    return out


def h1(D, kirilma_json):
    # (i) (sekans, tohum, pencere) eşli ATE_8 / ATE_0, A0; birim bootstrap
    g = defaultdict(list)
    for (k, d, s, t, w), r in D.items():
        if k != "A0" or d != 8.0:
            continue
        r0 = D.get(("A0", 0.0, s, t, w))
        a, b = ate(r), (ate(r0) if r0 else None)
        g[(s, t)].append(a / b if (a is not None and b is not None) else math.inf)
    U = sorted(g)
    tum = np.concatenate([g[u] for u in U])
    rng = np.random.default_rng(0)
    bs = [np.median(np.concatenate([g[U[i]] for i in rng.integers(0, len(U), len(U))])) for _ in range(B_BOOT)]
    oran8 = float(np.median(tum))
    K = json.load(open(kirilma_json))
    karar_i = "TUTTU" if oran8 <= 1.2 else ("ÇÜRÜDÜ (> 1.5)" if oran8 > 1.5 else "KISMEN (1.2–1.5)")
    k = K["kirilma"]
    karar_ii = "kırılma yok (test edilen aralıkta eşik yok)" if k is None else (
        "TUTTU" if 12 <= k <= 20 else f"TUTMADI (K = {k:g}, öngörü 12–20)")
    return dict(oran8_medyan=oran8, oran8_ga95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                n_pencere=int(len(tum)), n_birim=len(U), karar_ilk_yari=karar_i, kirilma=k, karar_kirilma=karar_ii,
                kirilma_tablo=K["tablo"])


def h3_euroc(D):
    v, ici, ust, alt = [], 0, 0, 0
    for (k, d, s, t, w), r in D.items():
        if k != "A0" or d != 0.0:
            continue
        x = anees(r, "pos")
        if x is None or not r.get("n_n"):
            continue
        lo, hi = bant(int(r["n_n"]))
        v.append(x)
        ici += lo <= x <= hi
        ust += x > hi
        alt += x < lo
    return dict(n_pencere=len(v), medyan_anees_poz=float(np.median(v)), bant_ici=ici, bant_ustu=ust, bant_alti=alt,
                not_="karar AirSim'de verildi (SONUC_AIRSIM_KOLLAR.md); burada betimsel")


def h4(D, aile_anees):
    sonuc = {"i": {}, "ii": {}, "ii_ikincil": {}}
    for kol in ("A3", "A1", "AG"):
        for dr in DRS:
            if dr < 8:
                continue
            x = aile_anees[f"{kol}@{dr:g}"]
            sonuc["i"][f"{kol}@{dr:g}"] = dict(medyan_oran=x["medyan_oran"], p_holm=x["p_holm"],
                                               anlamli_azalma=bool(x["medyan_oran"] < 1 and x["p_holm"] < 0.05))
    pler, ikincil = {}, {}
    alt_var = False
    for dr in DRS:
        if dr < 8:
            continue
        a, b, U, dus = eslesik(D, "A3", "AG", dr, "ANEES_poz")
        t = oran_testi(a, b)
        t["dusen_birim"] = dus
        sonuc["ii"][f"{dr:g}"] = t
        pler[f"{dr:g}"] = t["p"]
        if np.any(a < 3) or np.any(b < 3):
            alt_var = True
        da, db = np.abs(np.log(a / 3)), np.abs(np.log(b / 3))
        dd = da - db
        ikincil[f"{dr:g}"] = dict(medyan_fark=float(np.median(dd)),
                                  p=float(stats.wilcoxon(dd).pvalue) if np.any(dd != 0) else 1.0)
    ph = holm(pler)
    for k in sonuc["ii"]:
        sonuc["ii"][k]["p_holm"] = ph[k]
        sonuc["ii"][k]["A3_gecer"] = bool(sonuc["ii"][k]["medyan_oran"] < 1 and ph[k] < 0.05)
        sonuc["ii"][k]["fark_25_icinde"] = bool(0.75 <= sonuc["ii"][k]["medyan_oran"] <= 1.25)
    sonuc["ii_ikincil"] = dict(bant_alti_birim_var=alt_var, testler=ikincil if alt_var else "gerekmedi")
    gecer = [k for k, v in sonuc["ii"].items() if v["A3_gecer"]]
    sonuc["karar"] = ("A3 AG'yi Holm sonrası anlamlı geçiyor (D/r0: " + ", ".join(gecer) + ") -> fizik-tabanlı katkı korunur"
                      if gecer else "A3 AG'yi anlamlı geçemiyor -> 'türbülans fiziği gerekli' iddiası düşer, yeniden çerçevele")
    sonuc["ongoru_i_A3_A1_AG"] = {kol: all(v["anlamli_azalma"] for k, v in sonuc["i"].items() if k.startswith(kol + "@"))
                                  for kol in ("A3", "A1", "AG")}
    return sonuc


def h5(aile_ate):
    iyi = {k: v for k, v in aile_ate.items() if float(k.split("@")[1]) <= 8 and v["medyan_oran"] < 0.90
           and v["p_holm"] < 0.05}
    return dict(anlamli_10_iyilesme=sorted(iyi), karar="ÇÜRÜDÜ" if iyi else "TUTTU",
                d_r0_le8={k: dict(medyan_oran=v["medyan_oran"], p_holm=v["p_holm"]) for k, v in aile_ate.items()
                          if float(k.split("@")[1]) <= 8})


def kirilma_tablosu(kayitlar):
    """kayitlar: [(dr, kayıt)] -> kirilma_hesap ile aynı tanım (geçersiz ATE +inf; ıraksama ölçütü)"""
    A, I = defaultdict(list), defaultdict(list)
    for dr, r in kayitlar:
        v = r.get("m_ate_se3")
        A[dr].append(float(v) if (r.get("m_valid") and v is not None and math.isfinite(v)) else math.inf)
        ok = r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged")
        I[dr].append(0 if ok else 1)
    taban = float(np.median(A[0.0]))
    tab, K = [], None
    for dr in sorted(A):
        med, irak = float(np.median(A[dr])), float(np.mean(I[dr]))
        neden = ""
        if dr > 0 and K is None and (med > 2 * taban or irak > 0.20):
            K, neden = dr, ("ATE>2x" if med > 2 * taban else "irak>%20")
        tab.append(dict(D_over_r0=dr, n=len(A[dr]), medyan_ate=med, oran=med / taban, irak_oran=irak, neden=neden))
    return tab, K


def h6(D_ana, h6_yollar):
    yap = defaultdict(list)
    for (k, d, s, t, w), r in D_ana.items():
        if k == "A0" and d in (0.0, 4.0, 8.0):
            yap["r0.2_kp256"].append((d, r))
    for p in h6_yollar:
        for l in open(p):
            r = json.loads(l)
            tag = r["key"].split("/")[-1].replace("faz2_h6_", "")
            yap[tag].append((float(r["cell"]["D_over_r0"]), r))
    out = {}
    for tag in sorted(yap):
        tab, K = kirilma_tablosu(yap[tag])
        # betimsel: D/r0 = 4'te birim (sekans, tohum) oranı, kendi temizine göre
        g = defaultdict(lambda: {0.0: [], 4.0: []})
        for d, r in yap[tag]:
            if d in (0.0, 4.0) and ate(r) is not None:
                g[(r["cell"]["sequence"], int(r["cell"]["seed"]))][d].append(ate(r))
        U = [u for u in sorted(g) if g[u][0.0] and g[u][4.0]]
        a = np.array([np.mean(g[u][4.0]) for u in U])
        b = np.array([np.mean(g[u][0.0]) for u in U])
        out[tag] = dict(tablo=tab, kirilma=K, kayma=(None if K is None else K - 4.0), n_kayit=len(yap[tag]),
                        oran4=oran_testi(a, b))
    birincil = out.get("r0.4_kp512")
    if birincil is None:
        karar = "veri yok"
    elif birincil["kirilma"] is None:
        karar = "TUTTU (8'e kadar kırılma yok: kayma > 4)"
    else:
        kay = birincil["kayma"]
        karar = "TUTTU" if kay >= 4 else ("ÇÜRÜDÜ (kayma < 2: kırılma eşleştirme kalitesi kaynaklı)" if kay < 2
                                         else "KISMEN (2 ≤ kayma < 4)")
    return dict(yapilar=out, birincil="r0.4_kp512", karar=karar)


def betimsel_tablo(D):
    """kol × D/r0: birim değerlerinin medyanı (ATE, ANEES_poz)"""
    out = {}
    for kol in ("A0",) + KOLLAR + ("A5",):
        for dr in DRS:
            for metrik in ("ATE", "ANEES_poz"):
                f, tur = METRIK[metrik]
                g = defaultdict(list)
                for (k, d, s, t, w), r in D.items():
                    if k == kol and d == dr and f(r) is not None:
                        g[(s, t)].append(f(r))
                if g:
                    out.setdefault(kol, {}).setdefault(f"{dr:g}", {})[metrik] = dict(
                        medyan=float(np.median([_ort(v, tur) for v in g.values()])), n_birim=len(g))
    return out


# ----------------------------------------------------------------------------- öz-test
def oz_test():
    rng = np.random.default_rng(3)
    D = {}
    for s in ("S1", "S2", "S3"):
        for t in range(5):
            for w in (10, 20):
                for dr in DRS:
                    base_ate = 0.02 * (1 + dr / 4) * (1 + 0.05 * rng.standard_normal())
                    base_n = 150 * (1 + dr) * (1 + 0.05 * rng.standard_normal())
                    for kol, fa, fn in (("A0", 1, 1), ("A1", 1, 0.02), ("A3", 1, 0.2), ("AG", 1, 0.3),
                                        ("A0C", 1, 1 / 12), ("A0C+A3", 1, 0.02)):
                        D[(kol, dr, s, t, w)] = dict(status="ok", m_valid=True, n_valid=True, n_n=50,
                                                     m_ate_se3=base_ate * fa, n_nees_pos_mean=base_n * fn,
                                                     n_nees_rot_mean=base_n * fn, m_diverged=False,
                                                     cell=dict(sequence=s, seed=t, D_over_r0=dr))
    a, b, U, dus = eslesik(D, "A3", "A0", 8.0, "ANEES_poz")
    t = oran_testi(a, b, B=500)
    assert len(U) == 15 and dus == 0 and abs(t["medyan_oran"] - 0.2) < 1e-9 and t["p"] < 0.001, t
    a, b, U, _ = eslesik(D, "A3", "AG", 16.0, "ANEES_poz")
    assert abs(oran_testi(a, b, B=500)["medyan_oran"] - 2 / 3) < 1e-9
    ph = holm({"x": 0.01, "y": 0.04, "z": 0.03})
    assert abs(ph["x"] - 0.03) < 1e-12 and abs(ph["z"] - 0.06) < 1e-12 and abs(ph["y"] - 0.06) < 1e-12, ph
    lo, hi = bant(400)
    assert 2.7 < lo < 2.8 and 3.2 < hi < 3.3, (lo, hi)
    tab, K = kirilma_tablosu([(d, D[("A0", d, s, t, w)]) for (k, d, s, t, w) in D if k == "A0"])
    assert K == 4.0, tab        # 0.02*(1+4/4) = 2x -> 4'te tam 2x civarı; gürültüyle > 2x beklenir
    print("öz-test GEÇTİ")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ana", default=str(ROOT / "evidence/FAZ2/records.jsonl"))
    ap.add_argument("--kirilma", default=str(ROOT / "evidence/FAZ2/kirilma.json"))
    ap.add_argument("--h6", nargs="+", default=[str(ROOT / f"evidence/FAZ2_H6/records_p{i}.jsonl") for i in (0, 1, 3)])
    ap.add_argument("--cikti", default=str(ROOT / "evidence/FAZ2/analiz.json"))
    ap.add_argument("--oz_test", action="store_true")
    a = ap.parse_args()
    if a.oz_test:
        oz_test()
        return 0
    D = oku(a.ana)
    birimler = {(s, t) for (k, d, s, t, w) in D}
    assert len(birimler) == 45, f"birim {len(birimler)} != 45"
    out = {"kaynak": {"ana": a.ana, "h6": a.h6, "n_kayit": len(D)}}
    out["aile_ATE"] = aile(D, "ATE")
    out["aile_ANEES_poz"] = aile(D, "ANEES_poz")
    out["aile_ANEES_don"] = aile(D, "ANEES_don")
    out["H1"] = h1(D, a.kirilma)
    out["H3_euroc"] = h3_euroc(D)
    out["H4"] = h4(D, out["aile_ANEES_poz"])
    out["H5"] = h5(out["aile_ATE"])
    out["H6"] = h6(D, a.h6)
    out["betimsel"] = betimsel_tablo(D)
    json.dump(out, open(a.cikti, "w"), indent=1, ensure_ascii=False)
    h = out["H1"]
    print(f"H1: ATE oranı 8/0 = {h['oran8_medyan']:.2f} [{h['oran8_ga95'][0]:.2f}, {h['oran8_ga95'][1]:.2f}] "
          f"-> {h['karar_ilk_yari']}; kırılma K = {h['kirilma']} -> {h['karar_kirilma']}")
    h = out["H3_euroc"]
    print(f"H3 (EuRoC, betimsel): A0 temiz ANEES_poz medyan {h['medyan_anees_poz']:.1f}; bant içi {h['bant_ici']}, "
          f"üstü {h['bant_ustu']}, altı {h['bant_alti']} / {h['n_pencere']}")
    h = out["H4"]
    print("H4 (i) ANEES_poz kol/A0 (medyan oran, p_holm):")
    for k, v in h["i"].items():
        print(f"   {k:10s} {v['medyan_oran']:.3f}  p_holm={v['p_holm']:.2g}  {'anlamlı azalma' if v['anlamli_azalma'] else '-'}")
    print("H4 (ii) A3/AG:")
    for k, v in h["ii"].items():
        print(f"   D/r0={k:>2s}: oran {v['medyan_oran']:.3f} [{v['ga95'][0]:.3f}, {v['ga95'][1]:.3f}] n={v['n']} "
              f"p_holm={v['p_holm']:.2g} geçer={v['A3_gecer']} fark<%25={v['fark_25_icinde']}")
    print(f"   bant altı birim: {h['ii_ikincil']['bant_alti_birim_var']}")
    print(f"   KARAR: {h['karar']}")
    h = out["H5"]
    print(f"H5: anlamlı >%10 ATE iyileşmesi (D/r0 ≤ 8): {h['anlamli_10_iyilesme'] or 'yok'} -> {h['karar']}")
    h = out["H6"]
    for tag, v in h["yapilar"].items():
        print(f"H6 {tag:12s}: K={v['kirilma']} kayma={v['kayma']} oran@4={v['oran4']['medyan_oran']:.2f} "
              f"[{v['oran4']['ga95'][0]:.2f}, {v['oran4']['ga95'][1]:.2f}] n={v['n_kayit']}")
    print(f"H6 KARAR (birincil {h['birincil']}): {h['karar']}")
    print(f"-> {a.cikti}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
