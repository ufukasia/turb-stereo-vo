import json, sys, numpy as np
from collections import defaultdict


def yukle(p):
    return [json.loads(l) for l in open(p) if l.strip()]


def kalib(p):
    R = [r for r in yukle(p) if r.get("status") == "ok"]
    by = defaultdict(dict)       # tag -> (seq, w) -> nees
    for r in R:
        by[r["kol_etiket"]][(r["cell"]["sequence"], r["window"])] = r
    a0 = by["A0"]
    print(f"A0: {len(a0)} pencere; NEES_pos medyanları:",
          {k: round(v["n_nees_pos_median"], 1) for k, v in sorted(a0.items())})
    print("n_n:", {k: v["n_n"] for k, v in sorted(a0.items())})
    med_a0 = float(np.median([v["n_nees_pos_median"] for v in a0.values()]))
    print(f"A0 pencere-medyanlarının medyanı: {med_a0:.2f}")
    rows = []
    for k in (1, 2, 4, 8, 12, 16, 24):
        t = by.get(f"A0C_k{k}", {})
        common = sorted(set(t) & set(a0))
        act = float(np.median([t[c]["n_nees_pos_median"] for c in common])) if common else np.nan
        rot = float(np.median([t[c]["n_nees_rot_median"] for c in common])) if common else np.nan
        pred = med_a0 / k ** 2
        # pencere bazında gerçek / ölçekli oranı
        rr = [t[c]["n_nees_pos_median"] / (a0[c]["n_nees_pos_median"] / k ** 2) for c in common]
        ate = [t[c]["m_ate_se3"] / a0[c]["m_ate_se3"] for c in common]
        rows.append((k, len(common), act, pred, rot, np.median(rr) if rr else np.nan,
                     (min(rr), max(rr)) if rr else None, np.median(ate) if ate else np.nan))
    print("| k | n pencere | gerçek A0C medyan NEES_pos | A0/k² | NEES_rot | gerçek/ölçekli (medyan, min–maks) | ATE A0C/A0 |")
    for k, n, act, pred, rot, rr, mm, ate in rows:
        print(f"| {k} | {n} | {act:.2f} | {pred:.2f} | {rot:.2f} | {rr:.2f} ({mm[0]:.2f}–{mm[1]:.2f}) | {ate:.3f} |")
    best = min(rows, key=lambda x: abs(x[2] - 2.37))
    bestp = min(rows, key=lambda x: abs(x[3] - 2.37))
    print(f"SEÇİLEN k (gerçek koşum) = {best[0]}  (medyan {best[2]:.2f});  ölçekli tahminle = {bestp[0]}")
    tb = defaultdict(list)
    for r in R:
        tb[r["kol_etiket"]].append(r["t_backend"])
    print("t_backend medyan [s]:", {k: round(float(np.median(v)), 0) for k, v in tb.items()})


def duman(p):
    R = yukle(p)
    fail = [r for r in R if r.get("status") != "ok"]
    print(f"kayıt {len(R)}, ok değil {len(fail)}")
    for r in fail:
        print("  FAIL", r.get("key"), r.get("kol_etiket"), r.get("window"), r.get("error"))
    need = ["m_ate_se3", "ms_ate_se3", "m_rpe_trans", "n_nees_pos_median", "n_nees_rot_median", "n_n",
            "t_backend", "arm_stats"]
    miss = []
    for r in R:
        for k in need:
            if k == "arm_stats" and r.get("arm") == "A0":
                continue          # A0: kol nesnesi yok (PyEco yolu birebir) -> arm_stats None tasarım gereği
            v = r.get(k)
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                miss.append((r.get("key"), r.get("kol_etiket"), r.get("window"), k))
    print(f"eksik/NaN alan: {len(miss)}", miss[:10])
    T = defaultdict(dict)
    for r in R:
        if r.get("status") != "ok":
            continue
        c = r["cell"]
        T[(c["sequence"], c["D_over_r0"], r["window"])][r["kol_etiket"]] = r
    arms = sorted({r["kol_etiket"] for r in R})
    print("\n| sekans | D/r0 | pencere | kol | ATE [m] | ATE düzl. | NEES_pos | NEES_rot | n | t_be [s] | max_err |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for key in sorted(T):
        for a in arms:
            r = T[key].get(a)
            if r is None:
                continue
            print(f"| {key[0]} | {key[1]:g} | {key[2]} | {a} | {r['m_ate_se3']:.4f} | {r['ms_ate_se3']:.4f} | "
                  f"{r.get('n_nees_pos_median', np.nan):.1f} | {r.get('n_nees_rot_median', np.nan):.1f} | "
                  f"{r.get('n_n')} | {r['t_backend']:.0f} | {r['m_max_err']:.3f} |")
    # T8: A1 == A0 at dr0
    print("\nT8 (D/r0=0, A1 vs A0):")
    for key in sorted(T):
        if key[1] != 0:
            continue
        a, b = T[key].get("A0"), T[key].get("A1")
        if a and b:
            same = all(a.get(k) == b.get(k) for k in ("m_ate_se3", "ms_ate_se3", "n_nees_pos_median",
                                                      "n_nees_rot_median", "m_rpe_trans"))
            print(f"  {key}: bit-aynı={same}  ATE {a['m_ate_se3']!r} / {b['m_ate_se3']!r}")
    print("T9 (D/r0=0, A5 vs A0):")
    for key in sorted(T):
        if key[1] != 0:
            continue
        a, b = T[key].get("A0"), T[key].get("A5")
        if a and b:
            print(f"  {key}: ATE bağıl {abs(b['m_ate_se3']-a['m_ate_se3'])/a['m_ate_se3']:.3%}  "
                  f"NEES bağıl {abs(b['n_nees_pos_median']-a['n_nees_pos_median'])/a['n_nees_pos_median']:.1%}")
    print("T10 kalite / oracle sigma:")
    seen = set()
    for r in R:
        c = r["cell"]
        k = (c["sequence"], c["D_over_r0"])
        if k in seen:
            continue
        seen.add(k)
        q = r.get("quality")
        st = r.get("arm_stats") or {}
        print(f"  {k}: quality={json.dumps(q)[:300] if q else None}  fe_stats_stereo="
              f"{(r.get('frontend_stats') or {}).get('avg_kp2_count')}")
    for r in R:
        if r.get("kol_etiket") in ("A1", "A5", "A3") and r.get("window") == r["cell"]["windows"][0]:
            st = r.get("arm_stats") or {}
            print(f"  {r['cell']['sequence']} dr{r['cell']['D_over_r0']:g} {r['kol_etiket']}: "
                  f"sigma_px_oracle={st.get('sigma_px_oracle')}  sigma_hat_median={st.get('sigma_hat_median')}"
                  f"  sigma_turb_used_med={st.get('sigma_turb_used_median', st.get('sigma_turb_used'))}")
    tb = defaultdict(list)
    tf = {}
    for r in R:
        tb[r["kol_etiket"].split("_k")[0]].append(r["t_backend"])
        tf[(r["cell"]["sequence"], r["cell"]["D_over_r0"])] = r["t_frontend"]
    print("t_backend medyan/maks [s]:", {k: (round(float(np.median(v))), round(float(max(v)))) for k, v in tb.items()})
    print("t_frontend [s]:", {k: round(v) for k, v in tf.items()})


if __name__ == "__main__":
    {"kalib": kalib, "duman": duman}[sys.argv[1]](sys.argv[2])
