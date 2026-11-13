"""Rakip ön-test kayıtlarının (on_test.jsonl + r3a.jsonl) tablo özeti: R2, R3b, R4, R5. Çıktı JSON + ekrana tablo."""
from __future__ import annotations

import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
EV = ROOT / "evidence/FAZ2_rakip/on_test"


def oku(p):
    return [json.loads(l) for l in open(p)] if Path(p).exists() else []


def main():
    R = oku(EV / "on_test.jsonl")
    R3a = oku(EV / "r3a.jsonl")
    out = {}
    # ---- R2: (a) pencere ↔ (b) önek, D/r0=0 σ=1
    r2 = []
    for r in R:
        if r["protokol"] == "pencere":
            b = [x for x in R if x["protokol"] == "onek" and x["rakip"] == r["rakip"] and x["sequence"] == r["sequence"]
                 and x["D_over_r0"] == 0 and x["window"] == r["window"] and x["tekrar"] == 0]
            b = b[0] if b else None
            ilk = None
            try:
                from turb.rakip_kos import yorunge_yollari, _kestirim_damgalari, cam0_damgalari
                ned, _ = yorunge_yollari(r["rakip"], Path(r["kosum"]["out"]))
                et = _kestirim_damgalari(ned)
                ts = cam0_damgalari(Path("/data/sdata") / r["sequence"])
                ilk = (int(np.searchsorted(ts, et[0])) - r["window"]) if len(et) else None
            except Exception:  # noqa: BLE001
                pass
            r2.append({"rakip": r["rakip"], "sekans": r["sequence"], "pencere": r["window"],
                       "a_ate": r.get("m_ate_se3"), "a_kapsam": r["m_kapsam_p"], "a_durum": r["durum_p"],
                       "a_ilk_poz_kare_ofseti": ilk, "a_sure_s": r["kosum"]["sure_s"],
                       "b_ate": b and b.get("m_ate_se3"), "b_kapsam": b and b["m_kapsam_p"],
                       "b_durum": b and b["durum_p"], "b_sure_s": b and b["kosum"]["sure_s"],
                       "oran_a_b": (r["m_ate_se3"] / b["m_ate_se3"]) if (b and r.get("m_ate_se3") and b.get("m_ate_se3"))
                       else None})
    out["R2"] = r2
    # ---- R3b: rektifiye σ=0 ↔ ham bütün sekansın aynı pencereleri
    r3b = []
    for r in R:
        if r["protokol"] == "onek_sigma0":
            h = [x for x in R3a if x["rakip"] == r["rakip"] and x["window"] == r["window"]]
            h = h[0] if h else None
            r3b.append({"rakip": r["rakip"], "pencere": r["window"], "rekt_ate": r.get("m_ate_se3"),
                        "rekt_olcek": r.get("m_sim3_scale"), "ham_ate": h and h.get("m_ate_se3"),
                        "ham_olcek": h and h.get("m_sim3_scale"),
                        "oran": (r["m_ate_se3"] / h["m_ate_se3"]) if (h and r.get("m_ate_se3") and h.get("m_ate_se3")) else None})
    out["R3b"] = r3b
    # ---- R4: önek, tekrar 0
    r4 = defaultdict(dict)
    for r in R:
        if r["protokol"] == "onek" and r["tekrar"] == 0:
            r4[(r["rakip"], r["sequence"], r["window"])][r["D_over_r0"]] = r
    tab = []
    for (rk, s, w), d in sorted(r4.items()):
        row = {"rakip": rk, "sekans": s, "pencere": w}
        for dr in (0.0, 8.0, 24.0):
            x = d.get(dr)
            row[f"dr{dr:g}"] = None if x is None else {
                "ate": x.get("m_ate_se3"), "ate_son": x.get("ms_ate_se3"), "kapsam": x["m_kapsam_p"],
                "durum_p": x["durum_p"], "neden": x["neden"], "kosum_durum": x["kosum"]["durum"],
                "sure_s": x["kosum"]["sure_s"], "yuk": x["kosum"]["yuk_once"][0], "n_kare": x["kosum"]["n_kare_girdi"],
                "gpu": x["kosum"].get("gpu_adi"), "olcek": x.get("m_sim3_scale")}
            if dr > 0 and x is not None and d.get(0.0) and x.get("m_ate_se3") and d[0.0].get("m_ate_se3") \
                    and x["durum_p"] == "ok" and d[0.0]["durum_p"] == "ok":
                row[f"oran{dr:g}"] = x["m_ate_se3"] / d[0.0]["m_ate_se3"]
        tab.append(row)
    out["R4"] = tab
    # ---- R5: V1_01 D/r0=8, tekrarlar
    r5 = defaultdict(list)
    for r in R:
        if r["protokol"] == "onek" and r["sequence"] == "V1_01_easy" and r["D_over_r0"] == 8:
            r5[(r["rakip"], r["window"])].append(r)
    res5 = {}
    for (rk, w), rr in sorted(r5.items()):
        rr = sorted(rr, key=lambda x: x["tekrar"])
        a = [x.get("m_ate_se3") for x in rr]
        hs = []
        for x in rr:
            from turb.rakip_kos import yorunge_yollari
            ned, _ = yorunge_yollari(rk, Path(x["kosum"]["out"]))
            hs.append(hashlib.sha256(Path(ned).read_bytes()).hexdigest()[:16] if Path(ned).exists() else None)
        ok = [v for v in a if v is not None]
        res5.setdefault(rk, {})[str(w)] = {"ate": a, "durum_p": [x["durum_p"] for x in rr], "sha256_nedensel": hs,
                                           "bayt_ayni": len(set(hs)) == 1 and len(hs) > 1,
                                           "cv": float(np.std(ok, ddof=1) / np.mean(ok)) if len(ok) >= 2 else None}
    for rk, d in res5.items():
        cvs = [v["cv"] for v in d.values() if v["cv"] is not None]
        d["cv_medyan"] = float(np.median(cvs)) if cvs else None
    out["R5"] = res5
    (EV / "on_test_ozet.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    f = lambda v, n=4: "NA" if v is None else (f"{v:.{n}f}" if isinstance(v, float) else str(v))
    print("R2 (a)/(b):")
    for x in r2:
        print(f"  {x['rakip']:9s} {x['sekans']:11s} w{x['pencere']:<5d} a={f(x['a_ate'])} ({x['a_durum']}, kapsam {f(x['a_kapsam'], 3)}, "
              f"ilk poz +{x['a_ilk_poz_kare_ofseti']}, {x['a_sure_s']:.0f}s)  b={f(x['b_ate'])} ({x['b_durum']})  a/b={f(x['oran_a_b'], 2)}")
    print("R3b:")
    for x in r3b:
        print(f"  {x['rakip']:9s} w{x['pencere']:<5d} rekt={f(x['rekt_ate'])} (s {f(x['rekt_olcek'], 3)}) ham={f(x['ham_ate'])} oran={f(x['oran'], 2)}")
    print("R4:")
    for x in tab:
        cells = []
        for dr in ("0", "8", "24"):
            c = x.get(f"dr{dr}")
            cells.append("—" if c is None else f"{f(c['ate'])}/{c['durum_p'][:3]}/k{f(c['kapsam'], 2)}")
        print(f"  {x['rakip']:9s} {x['sekans']:11s} w{x['pencere']:<5d} " + "  ".join(cells)
              + f"  oran8={f(x.get('oran8'), 2)} oran24={f(x.get('oran24'), 2)}")
    print("R5:", json.dumps(res5, ensure_ascii=False, default=str)[:1500])


if __name__ == "__main__":
    sys.path.insert(0, str(ROOT))
    main()
