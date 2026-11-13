#!/usr/bin/env python
"""Hızlı rejim ızgarası ön testi: ölçülen / beklenen tablosu (evidence/HIZLI/ON_TEST.md §2'deki kapılar).
Bu betik ön test koşulmadan önce yazıldı; beklenenler ON_TEST.md'de. Çıkış kodu 0 = bütün DUR kapıları geçti."""
import json, math, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
import hizli_kos as hk  # noqa: E402

BIT = hk.BIT_ALANLARI
SEQ = "V1_01_easy"


def oku(p, kosul):
    R = []
    if os.path.exists(p):
        for l in open(p):
            try:
                r = json.loads(l)
            except json.JSONDecodeError:
                R.append({"_bozuk": True}); continue
            if r.get("_bozuk") or "arm" not in r or kosul(r):
                R.append(r)
    return R


def gort(v):
    return float(np.exp(np.mean(np.log(v))))


def main():
    T = oku(os.path.join(hk.EV, "on_test.jsonl"), lambda r: True)
    bozuk = sum(1 for r in T if r.get("_bozuk"))
    hata = [r for r in T if not r.get("_bozuk") and "arm" not in r]
    K = [r for r in T if "arm" in r]
    A = {}
    for l in open(os.path.join(ROOT, "evidence", "FAZ2", "records.jsonl")):
        r = json.loads(l); c = r["cell"]
        if c["sequence"] == SEQ and c["seed"] == 1 and c["D_over_r0"] in (0.0, 4.0, 8.0):
            A[(r["arm"], c["D_over_r0"], int(r["window"]))] = r
    sat, dur = [], []

    def ekle(no, ad, olc, bek, ok, kapi=True):
        sat.append((no, ad, olc, bek, "geçti" if ok else ("TUTMADI" if kapi else "tutmadı (not)")))
        if kapi and not ok:
            dur.append(no)

    # G2: tohum 0 hücresi sondayla bit-aynı
    ref = hk.sonda_kayitlari().get((SEQ, 0, 8.0, round(hk.tau0_over_dt("R5", 8), 6)), {})
    s0 = [r for r in K if r["cell"]["seed"] == 0]
    nb = sum(1 for r in s0 if (m := ref.get((r["arm"], int(r["window"])))) and all(r.get(k) == m.get(k) for k in BIT))
    fe = bool(s0) and hk.fe_ayni(s0[0].get("frontend_stats"), ref.get("fe"))
    ekle("G2", "tohum 0, R5, D/r0 8, A0 + A1 ↔ sonda", f"{nb}/{len(s0)} kayıt bit-aynı; ön-uç istatistikleri aynı: {fe}",
         "6/6 bit-aynı; ön-uç aynı", nb == 6 and len(s0) == 6 and fe)
    # G3: duman
    s1 = [r for r in K if r["cell"]["seed"] == 1]
    dolu = all(all(r.get(k) is not None for k in hk.DOLU_ALANLAR) for r in s1 if r.get("status") == "ok")
    durum = {}
    for r in s1:
        durum[r.get("status")] = durum.get(r.get("status"), 0) + 1
    ekle("G3", "tohum 1, D/r0 8 ve 4, altı kol", f"{len(s1)} kayıt, durum {durum}, hücre hatası {len(hata)}, bozuk satır "
         f"{bozuk}, 'ok' kayıtlarda alanlar dolu: {dolu}", "36 kayıt; hücre hatası 0; bozuk 0; alanlar dolu; ok ≥ 30",
         len(s1) == 36 and not hata and not bozuk and dolu and durum.get("ok", 0) >= 30)
    fp = {r.get("code_fp") for r in K}
    ekle("G1", "kod parmak izi (kayıtlarda)", str(fp), "{'08e28b83407c'}", fp == {hk.BEKLENEN_FP})
    for dr, t0 in ((8.0, 0.0157), (4.0, 0.0314)):
        C = [r for r in s1 if r["cell"]["D_over_r0"] == dr]
        if not C:
            ekle("G4", f"D/r0 {dr:g} kayıtları", "yok", "var", False); continue
        g = lambda kol, alan: [r.get(alan) for r in C if r["arm"] == kol and r.get("status") == "ok" and r.get(alan)]  # noqa: E731
        ist = lambda kol: next((r["arm_stats"] for r in C if r["arm"] == kol and r.get("arm_stats")), {})  # noqa: E731
        ana_ist = lambda kol: next((A[k]["arm_stats"] for k in sorted(A) if k[0] == kol and k[1] == dr), {})  # noqa: E731
        t = {round(float(r["cell"]["tau0_over_dt"]), 6) for r in C}
        ekle("G4a", f"D/r0 {dr:g}: tau0/dt", str(t), f"{{{t0}}} (v = 5 m/s)", t == {t0})
        fa, fb = C[0]["frontend_stats"], A[("A0", dr, min(k[2] for k in A if k[1] == dr))]["frontend_stats"]
        ekle("G4b", f"D/r0 {dr:g}: rejim uygulandı (ön-uç istatistikleri R0'dan farklı; zamansal eşleşme R5 / R0)",
             f"aynı: {hk.fe_ayni(fa, fb)}; avg_temporal_matches {fa.get('avg_temporal_matches'):.1f} / "
             f"{fb.get('avg_temporal_matches'):.1f}; kp1 {fa.get('avg_kp1_count'):.0f} / {fb.get('avg_kp1_count'):.0f}",
             "farklı; anahtar nokta sayısı ±%10 içinde",
             (not hk.fe_ayni(fa, fb)) and abs(fa["avg_kp1_count"] / fb["avg_kp1_count"] - 1) <= 0.10)
        s_a1, s_a1m = ist("A1").get("sigma_turb_used_median"), ana_ist("A1").get("sigma_turb_used_median")
        ekle("G4c", f"D/r0 {dr:g}: A1 oracle σ (kullanılan)", f"{s_a1}", f"{s_a1m} (ana ızgarayla aynı)", s_a1 == s_a1m)
        ekle("G4c", f"D/r0 {dr:g}: A0C cal_factor; A3 yama faktörü", f"{ist('A0C').get('cal_factor')}; "
             f"{ist('A3').get('n_patch_factors')}", "12.0; > 0", ist("A0C").get("cal_factor") == 12.0
             and (ist("A3").get("n_patch_factors") or 0) > 0)
        sh, shm = ist("A3").get("sigma_hat_median"), ana_ist("A3").get("sigma_hat_median")
        ekle("G4d", f"D/r0 {dr:g}: A3 σ̂ medyanı R5 / R0 (aynı hücre)", f"{sh} / {shm} = {sh / shm:.2f}" if sh and shm
             else f"{sh} / {shm}", "0.75–1.33", bool(sh and shm and 0.75 <= sh / shm <= 1.33))
        ag, agm = ist("AG").get("ag_rms_median"), ana_ist("AG").get("ag_rms_median")
        ekle("G4e", f"D/r0 {dr:g}: AG artık RMS medyanı R5 / R0 (öngörü: büyür)", f"{ag} / {agm} = {ag / agm:.2f}"
             if ag and agm else f"{ag} / {agm}", "> 1 (öngörü; kapı değil)", bool(ag and agm and ag > agm), kapi=False)
        an = {kol: (gort(g(kol, "n_nees_pos_mean")) if g(kol, "n_nees_pos_mean") else math.nan) for kol in hk.KOLLAR}
        ekle("G4f", f"D/r0 {dr:g}: ANEES_poz (pencere geometrik ort.) A1, A0C < A0", "A1/A0 %.3f; A0C/A0 %.3f" %
             (an["A1"] / an["A0"], an["A0C"] / an["A0"]), "ikisi de < 1 (etki uygulandı)",
             an["A1"] < an["A0"] and an["A0C"] < an["A0"])
        ekle("G4f", f"D/r0 {dr:g}: ANEES_poz A3, AG < A0", "A3/A0 %.3f; AG/A0 %.3f" % (an["A3"] / an["A0"],
             an["AG"] / an["A0"]), "ikisi de < 1 (beklenen; kapı değil)", an["A3"] < an["A0"] and an["AG"] < an["A0"],
             kapi=False)
        at = g("A0", "m_ate_se3")
        tem = [A[k]["m_ate_se3"] for k in A if k[0] == "A0" and k[1] == 0.0]
        ekle("G4g", f"D/r0 {dr:g}: A0 ATE [m] (pencereler); temiz ana ızgara medyanı", f"{[round(x, 4) for x in at]}; "
             f"{np.median(tem):.4f}", "sonlu, 0.005–3 m aralığında", bool(at) and all(0.005 <= x <= 3 for x in at))
        nc = [(r["m_n_compare"], A[(r["arm"], dr, int(r["window"]))]["m_n_compare"]) for r in C if r["arm"] == "A0"
              and r.get("status") == "ok" and (r["arm"], dr, int(r["window"])) in A]
        ekle("G4h", f"D/r0 {dr:g}: pencere tanımı ana ızgarayla aynı (karşılaştırılan poz sayısı)", str(nc),
             "çiftler ±%5 içinde", bool(nc) and all(abs(a / b - 1) <= 0.05 for a, b in nc))
        tb = {kol: next((r.get("t_backend") for r in C if r["arm"] == kol), None) for kol in hk.KOLLAR}
        sat.append(("G6", f"D/r0 {dr:g}: süre [s] ön-uç; arka-uç kol başına", f"{C[0].get('t_frontend'):.0f}; "
                    + ", ".join(f"{k} {v:.0f}" for k, v in tb.items() if v), "bilgi (toplam süre kestirimi için)", "—"))
    print("| no | denetim | ölçülen | beklenen | sonuç |\n|---|---|---|---|---|")
    for s in sat:
        print("| " + " | ".join(str(x) for x in s) + " |")
    print("\nÖN TEST:", "GEÇTİ" if not dur else f"TUTMADI ({', '.join(dur)}) → ızgara başlatılmaz")
    return 1 if dur else 0


if __name__ == "__main__":
    sys.exit(main())
