import json, sys, numpy as np
from collections import defaultdict
F = "/tmp/faz2_on_test"
R = [json.loads(l) for l in open(f"{F}/temiz_dagilim.jsonl")]
print("kayıt", len(R), "ok", sum(r.get("status") == "ok" for r in R), "fe_gpu", {r.get("fe_gpu") for r in R},
      "code_fp", {r.get("code_fp") for r in R}, "/n1 anahtar", all("/n1" in r["key"] for r in R))
# gürültüsüz referans (duman testi, tohum 0)
D = [json.loads(l) for l in open(f"{F}/duman.jsonl")]
ref = {(r["cell"]["sequence"], r["cell"]["D_over_r0"], r["window"]): r for r in D if r["kol_etiket"] == "A0"}
A = defaultdict(dict)          # (seq, dr, w) -> seed -> ate
for r in R:
    if r.get("status") != "ok":
        print("FAIL", r.get("key"), r.get("error")); continue
    c = r["cell"]
    A[(c["sequence"], c["D_over_r0"], r["window"])][c["seed"]] = r
print("\n| sekans/pencere | temiz ATE tohum 0–4 [mm] | temiz CV | gürültüsüz temiz [mm] | medyan/gürültüsüz | D/r₀=8 ATE 0–4 [mm] | D/r₀=8 CV | eşli oran 8/0 (tohum 0–4) |")
print("|---|---|---|---|---|---|---|---|")
allr, cvs0, cvs8 = [], [], []
for s in ("MH02", "V1_01_easy"):
    for w in sorted({k[2] for k in A if k[0] == s}):
        c0 = A.get((s, 0.0, w), {}); c8 = A.get((s, 8.0, w), {})
        a0 = np.array([c0[i]["m_ate_se3"] for i in sorted(c0)]); a8 = np.array([c8[i]["m_ate_se3"] for i in sorted(c8)])
        cv0 = a0.std(ddof=1) / a0.mean() if len(a0) > 1 else np.nan
        cv8 = a8.std(ddof=1) / a8.mean() if len(a8) > 1 else np.nan
        cvs0.append(cv0); cvs8.append(cv8)
        rr = [c8[i]["m_ate_se3"] / c0[i]["m_ate_se3"] for i in sorted(set(c0) & set(c8))]
        allr += rr
        nf = ref.get((s, 0.0, w), {}).get("m_ate_se3", np.nan)
        print(f"| {s}/{w} | {', '.join(f'{x*1000:.1f}' for x in a0)} | {cv0:.0%} | {nf*1000:.1f} | {np.median(a0)/nf:.2f} | "
              f"{', '.join(f'{x*1000:.0f}' for x in a8)} | {cv8:.0%} | {', '.join(f'{x:.1f}' for x in rr)} (medyan {np.median(rr) if rr else np.nan:.1f}) |")
allr = np.array(allr)
print(f"\nçift sayısı {len(allr)}; oran > 1.5: {np.sum(allr > 1.5)}/{len(allr)} ({np.mean(allr > 1.5):.0%}); "
      f"medyan {np.median(allr):.2f}; min {allr.min():.2f}; p10 {np.percentile(allr, 10):.2f}; maks {allr.max():.1f}")
print(f"temiz CV medyanı {np.nanmedian(cvs0):.0%} (aralık {np.nanmin(cvs0):.0%}–{np.nanmax(cvs0):.0%}); "
      f"D/r0=8 CV medyanı {np.nanmedian(cvs8):.0%} (aralık {np.nanmin(cvs8):.0%}–{np.nanmax(cvs8):.0%})")
# alternatif: oran temiz medyanına göre (tek temiz yerine dağılım medyanı)
alt = []
for s in ("MH02", "V1_01_easy"):
    for w in sorted({k[2] for k in A if k[0] == s}):
        c0 = A[(s, 0.0, w)]; c8 = A[(s, 8.0, w)]
        m0 = np.median([c0[i]["m_ate_se3"] for i in c0])
        alt.append(np.median([c8[i]["m_ate_se3"] for i in c8]) / m0)
print("pencere: medyan(8)/medyan(0):", [round(x, 1) for x in alt])
# NEES (bilgi), kalite, süre
q = {}
tf = defaultdict(list); tb = defaultdict(list)
for r in R:
    c = r["cell"]; k = (c["sequence"], c["D_over_r0"])
    q.setdefault(k, set()).add(None if not r.get("quality") else round(r["quality"]["psnr"], 2))
    if r["window"] == c["windows"][0]:
        tf[k].append(r["t_frontend"]); tb[k].append(r["t_backend"])
print("PSNR:", {k: sorted(v, key=lambda x: (x is None, x)) for k, v in q.items()})
print("süre [s] (ön-uç / arka-uç A0, medyan ve aralık):")
for k in sorted(tf):
    print(f"  {k}: ön-uç {np.median(tf[k]):.0f} ({min(tf[k]):.0f}–{max(tf[k]):.0f}), arka-uç {np.median(tb[k]):.0f} ({min(tb[k]):.0f}–{max(tb[k]):.0f})")
print("NEES_pos medyanı, D/r0=0 tohumlar:", {f"{k[0]}/{k[2]}": [round(v[i]['n_nees_pos_median']) for i in sorted(v)] for k, v in A.items() if k[1] == 0})
