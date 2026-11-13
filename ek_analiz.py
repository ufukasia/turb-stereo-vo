"""Ana ızgara kayıtlarından BETİMSEL ek tablolar (ön kayıt dışı; yeni koşum yok, karar üretmez).

    python ek_analiz.py        # -> evidence/FAZ2/ek_analiz.json + özet

T1 kestirici doğruluğu : DIMM σ̂ (A3) / oracle σ_tilt, AG artık-RMS'i, kollarda kullanılan σ (D/r0 başına)
T2 ön-uç verimi        : stereo eşleşme, üçgenlenen, verim = üçgenlenen / eşleşme, zamansal eşleşme, PSNR, SSIM
T3 skaler açıklama     : A3 / AG ANEES oranı yalnız kullanılan skaler σ'dan öngörülebiliyor mu?
                         öngörü = (σ0² + σ_AG²) / (σ0² + σ_A3²), σ0 = 1 px; ölçülen = ANEES(A3) / ANEES(AG)

Birim (sekans, tohum); birim değeri pencerelerin ortalaması (σ aritmetik, oranlar geometrik); tablo değeri
birimlerin medyanı, [çeyrekler], n. Bu dosya parmak izi kümesinin dışındadır.
"""
import json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from turb.faz2_analiz import oku, eslesik, anees  # noqa: E402

KAYIT = ROOT / "evidence" / "FAZ2" / "records.jsonl"
CIKTI = ROOT / "evidence" / "FAZ2" / "ek_analiz.json"
DRS = (0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0)
S_TOT = math.hypot(1.0, 1.75)  # A1'in kullandığı toplam ölçek: s_tot = hypot(1, 1.75)·σ_tilt


def ozet(x):
    x = np.asarray([v for v in x if v is not None and math.isfinite(v)], float)
    if len(x) == 0:
        return dict(n=0)
    return dict(n=int(len(x)), medyan=float(np.median(x)), c25=float(np.percentile(x, 25)),
                c75=float(np.percentile(x, 75)))


def birim_ort(D, kol, dr, f, geo=False):
    """-> {(sekans, tohum): pencere ortalaması}; f(kayıt) -> sayı | None"""
    g = defaultdict(list)
    for (k, d, s, t, w), r in D.items():
        if k != kol or d != dr or r.get("status") != "ok":
            continue
        v = f(r)
        if v is not None and math.isfinite(v):
            g[(s, t)].append(v)
    return {u: (float(np.exp(np.mean(np.log(v)))) if geo else float(np.mean(v))) for u, v in g.items() if v}


def ks(alan):
    return lambda r: (r.get("arm_stats") or {}).get(alan)


def t1_kestirici(D):
    out = {}
    for dr in DRS:
        orc = birim_ort(D, "A3", dr, ks("sigma_px_oracle"))
        sh = birim_ort(D, "A3", dr, ks("sigma_hat_mean"))
        ag = birim_ort(D, "AG", dr, ks("ag_rms_mean"))
        kul = {k: birim_ort(D, k, dr, ks("sigma_turb_used_mean")) for k in ("A1", "A3", "AG", "A0C+A3")}
        U = sorted(set(orc) & set(sh))
        satir = dict(
            oracle_sigma_tilt_px=ozet(orc.values()),
            dimm_sigma_hat_px=ozet(sh.values()),
            ag_rms_px=ozet(ag.values()),
            kullanilan_sigma_px={k: ozet(v.values()) for k, v in kul.items()},
        )
        if dr > 0:
            satir["dimm_bolu_sigma_tilt"] = ozet([sh[u] / orc[u] for u in U])
            satir["dimm_bolu_s_tot"] = ozet([sh[u] / (S_TOT * orc[u]) for u in U])
            Ua = sorted(set(kul["A3"]) & set(kul["A1"]))
            satir["A3_kullanilan_bolu_A1_kullanilan"] = ozet([kul["A3"][u] / kul["A1"][u] for u in Ua])
        out[f"{dr:g}"] = satir
    return out


def t2_on_uc(D):
    out = {}
    for dr in DRS:
        hucre = {}
        for (k, d, s, t, w), r in D.items():
            if k == "A0" and d == dr and (s, t) not in hucre:
                hucre[(s, t)] = r
        fs = lambda a: [(r.get("frontend_stats") or {}).get(a) for r in hucre.values()]
        esl, ucg = fs("avg_stereo_matches"), fs("avg_current_stereo_triangulated_points")
        out[f"{dr:g}"] = dict(
            stereo_eslesme=ozet(esl), ucgenlenen=ozet(ucg),
            verim=ozet([u / e for u, e in zip(ucg, esl) if e]),
            zamansal_eslesme=ozet(fs("avg_temporal_matches")),
            etkin_derinlik_izi=ozet(fs("avg_active_depth_tracks")),
            anahtar_nokta=ozet(fs("avg_kp1_count")),
            psnr_db=ozet([(r.get("quality") or {}).get("psnr") for r in hucre.values()]),
            ssim=ozet([(r.get("quality") or {}).get("ssim") for r in hucre.values()]),
        )
    return out


def t3_skaler(D):
    out = {}
    for dr in DRS[2:]:
        g = defaultdict(lambda: ([], []))
        for (k, d, s, t, w), r in D.items():
            if k != "A3" or d != dr:
                continue
            r2 = D.get(("AG", dr, s, t, w))
            if r2 is None:
                continue
            a3, ag = anees(r), anees(r2)
            s3, sg = ks("sigma_turb_used_mean")(r), ks("sigma_turb_used_mean")(r2)
            if None in (a3, ag, s3, sg) or not (math.isfinite(s3) and math.isfinite(sg)):
                continue
            g[(s, t)][0].append(math.log((1.0 + sg * sg) / (1.0 + s3 * s3)))
            g[(s, t)][1].append(math.log(a3 / ag))
        U = sorted(g)
        ong = np.array([np.mean(g[u][0]) for u in U])
        olc = np.array([np.mean(g[u][1]) for u in U])
        rho, p = stats.spearmanr(ong, olc)
        a, b, U2, _ = eslesik(D, "A3", "AG", dr, "ANEES_poz")
        out[f"{dr:g}"] = dict(
            n=len(U), ongoru_medyan=float(np.exp(np.median(ong))), olculen_medyan=float(np.exp(np.median(olc))),
            olculen_medyan_eslesik=float(np.exp(np.median(np.log(a / b)))), n_eslesik=len(U2),
            artik_medyan=float(np.exp(np.median(olc - ong))),
            artik_c25_c75=[float(np.exp(np.percentile(olc - ong, 25))), float(np.exp(np.percentile(olc - ong, 75)))],
            spearman_rho=float(rho), spearman_p=float(p),
        )
    return out


def main():
    D = oku(KAYIT)
    S = dict(kaynak=str(KAYIT), n_kayit=len(D), statu="betimsel (ön kayıt dışı)",
             T1_kestirici=t1_kestirici(D), T2_on_uc=t2_on_uc(D), T3_skaler_aciklama=t3_skaler(D))
    CIKTI.write_text(json.dumps(S, indent=1, ensure_ascii=False))
    m = lambda o, b="%.3f": (b % o["medyan"]) if o.get("n") else "—"
    print("T1 kestirici (birim medyanı, px)")
    print(" D/r0  σ_tilt  σ̂_DIMM  σ̂/σ_tilt  σ̂/s_tot  AG_rms  kullanılan: A1     A3     AG   A0C+A3   n")
    for dr, s in S["T1_kestirici"].items():
        k = s["kullanilan_sigma_px"]
        print(f" {dr:>4}  {m(s['oracle_sigma_tilt_px'])}   {m(s['dimm_sigma_hat_px'])}    "
              f"{m(s.get('dimm_bolu_sigma_tilt', {}))}     {m(s.get('dimm_bolu_s_tot', {}))}   {m(s['ag_rms_px'])}"
              f"        {m(k['A1'])}  {m(k['A3'])}  {m(k['AG'])}  {m(k['A0C+A3'])}   {s['dimm_sigma_hat_px'].get('n')}")
    print("\nT2 ön-uç (hücre medyanı)")
    print(" D/r0  eşleşme  üçgen.  verim  zamansal  derinlik izi  PSNR    SSIM    n")
    for dr, s in S["T2_on_uc"].items():
        print(f" {dr:>4}  {m(s['stereo_eslesme'], '%6.2f')}  {m(s['ucgenlenen'], '%6.2f')}  {m(s['verim'])}  "
              f"{m(s['zamansal_eslesme'], '%7.2f')}   {m(s['etkin_derinlik_izi'], '%7.2f')}     "
              f"{m(s['psnr_db'], '%5.2f')}  {m(s['ssim'])}  {s['verim'].get('n')}")
    print("\nT3 skaler açıklama: ANEES(A3)/ANEES(AG), öngörü (1+σ_AG²)/(1+σ_A3²)")
    print(" D/r0  öngörü  ölçülen  ölçülen/öngörü [ç25, ç75]   Spearman ρ (p)   n")
    for dr, s in S["T3_skaler_aciklama"].items():
        print(f" {dr:>4}  {s['ongoru_medyan']:.3f}   {s['olculen_medyan']:.3f}    {s['artik_medyan']:.3f} "
              f"[{s['artik_c25_c75'][0]:.3f}, {s['artik_c25_c75'][1]:.3f}]      {s['spearman_rho']:+.2f} "
              f"({s['spearman_p']:.2g})   {s['n']}")
    print(f"\n-> {CIKTI}")


if __name__ == "__main__":
    main()
