#!/usr/bin/env python
"""Tables of README.md, generated from the record and analysis files in results/.

    python scripts/readme_tables.py            print every table as GitHub markdown
    python scripts/readme_tables.py --check    verify that every generated table row is present, verbatim, in README.md

No number in the README tables is typed by hand: the README holds the lines printed here.
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COND = [("R0-8", "Quasi-static, D/r0 = 8"), ("R0-16", "Quasi-static, D/r0 = 16"),
        ("R5-8", "Fast, D/r0 = 8"), ("R5-4", "Fast, D/r0 = 4")]


def J(*p):
    return json.load(open(os.path.join(ROOT, *p)))


def L(*p):
    return [json.loads(x) for x in open(os.path.join(ROOT, *p)) if x.strip()]


def ci(d, nd=3):
    return f"{d['medyan_oran']:.{nd}f} [{d['ga95'][0]:.{nd}f}, {d['ga95'][1]:.{nd}f}]"


def bold_min(vals, fmt="{:.3f}"):
    """Format a row of numbers, the smallest in bold."""
    m = min(v for v in vals if v is not None)
    return [("–" if v is None else (f"**{fmt.format(v)}**" if v == m else fmt.format(v))) for v in vals]


def hero_cell():
    """ATE of every window of the hero cell (MH05, seed 3, fast regime, D/r0 = 8)."""
    sel = lambda r: (r["cell"]["sequence"] == "MH05" and r["cell"]["seed"] == 3 and r["cell"]["D_over_r0"] == 8.0)
    ate = {}
    for r in L("results", "HIZLI", "records.jsonl"):
        if sel(r) and r["arm"] == "A0" and r["cell"]["tag"] == "hizli_R5":
            ate[("A0", r["window"])] = r["m_ate_se3"]
    for r in L("results", "KAPI", "dogrula_I.jsonl"):
        if sel(r):
            ate[(r["cell"]["tag"].replace("kapid_", ""), r["window"])] = r["m_ate_se3"]
    wins = sorted({w for _, w in ate})
    out = ["| Window (frames) | Baseline | Scalar Variant | **Adaptive Gate** | Oracle Gate | Gate / Baseline |",
           "|:---|---:|---:|---:|---:|---:|"]
    for w in wins:
        v = [ate[(a, w)] for a in ("A0", "A3S", "KAPI")]
        cells = bold_min(v)
        out.append(f"| {w}–{w + 399} | {cells[0]} | {cells[1]} | {cells[2]} | {ate[('A1g', w)]:.3f} | "
                   f"{ate[('KAPI', w)] / ate[('A0', w)]:.3f} |")
    return out


def highlights():
    K = J("results", "KAPI", "dogrula_analiz.json")
    G = J("results", "SAPMA7", "s7_analiz.json")
    out = ["| Condition | **Adaptive Gate** / Baseline | Gated Generic / Baseline | Adaptive Gate / Gated Generic | "
           "Scalar Variant / Baseline | Oracle Gate / Baseline | Position ANEES, Gate / Baseline |",
           "|:---|---:|---:|---:|---:|---:|---:|"]
    for c, name in COND:
        fam = "B" if f"KAPI/A0@{c}" in K["B"] else "I1"
        out.append(f"| {name} | **{ci(K[fam][f'KAPI/A0@{c}'])}** | {G['HG1'][f'AGg/A0@{c}']['medyan_oran']:.3f} | "
                   f"{ci(G['HG2'][f'KAPI/AGg@{c}'])} | {K['I4'][f'A3S/A0@{c}']['medyan_oran']:.3f} | "
                   f"{K['I4'][f'A1g/A0@{c}']['medyan_oran']:.3f} | {K['I5_poz'][f'KAPI/A0@{c}']['medyan_oran']:.3f} |")
    out.append(f"| No turbulence | {ci(K['T']['KAPI/A0@R0-0'])} | {G['HG3']['AGg/A0@R0-0']['medyan_oran']:.3f} | "
               f"{ci(G['T_ATE']['KAPI/AGg@R0-0'])} | {K['T2_ATE']['A3S/A0@R0-0']['medyan_oran']:.3f} | – | "
               f"{K['T2_ANEES_poz']['KAPI/A0@R0-0']['medyan_oran']:.3f} |")
    return out


def degradation():
    F = J("results", "FAZ2", "analiz.json")["H1"]["kirilma_tablo"]
    R = J("results", "FAZ2_rakip_ek", "rakip_ek_analiz.json")["HR1"]["rakipler"]
    lv = [4, 8, 16, 24]
    out = ["| System | D/r0 = 4 | D/r0 = 8 | D/r0 = 16 | D/r0 = 24 | Failed Windows at 16 [%] | Failed Windows at 24 [%] |",
           "|:---|---:|---:|---:|---:|---:|---:|"]
    f = {int(r["D_over_r0"]): r for r in F}
    out.append("| PyEco-SLAM (studied baseline) | " + " | ".join(f"{f[d]['oran']:.2f}" for d in lv) +
               f" | {100 * f[16]['irak_oran']:.1f} | {100 * f[24]['irak_oran']:.1f} |")
    for key, name in (("orbslam3", "ORB-SLAM3 (stereo)"), ("airslam", "AirSLAM"), ("droid", "DROID-SLAM")):
        cells = []
        for d in lv:
            v = R[key][str(d)]["medyan_oran"]
            cells.append("no valid window" if v is None else f"{v:.2f}")
        out.append(f"| {name} | " + " | ".join(cells) +
                   f" | {100 * R[key]['16']['basarisizlik']:.1f} | {100 * R[key]['24']['basarisizlik']:.1f} |")
    return out


def real():
    S = J("results", "GERCEK_VERI", "saha_ozet.json")["ort15"]["2021-07"]
    rows = [("Rank correlation with scintillometer", "G1_spearman", "G1_ga", "at least 0.6", "G1_gecer"),
            ("Log-log slope", "G2_egim", "G2_ga", "0.6 to 1.4", "G2_gecer"),
            ("Median ratio, estimated / measured Cn2", "G3_oran", "G3_ga", "0.2 to 0.8", "G3_gecer")]
    out = ["| Criterion (declared before processing) | Value [95% CI] | Outcome |", "|:---|---:|:---|"]
    for name, k, g, _, ok in rows:
        out.append(f"| {name} | {S[k]:.3f} [{S[g][0]:.3f}, {S[g][1]:.3f}] | {'passed' if S[ok] else 'failed'} |")
    out.append(f"| Median R2 of the structure-function fit | {S['G4_r2_ortanca']:.3f} | "
               f"{'passed' if S['G4_gecer'] else 'failed'} |")
    return out


def factorial():
    H = J("results", "FAZ2_FAKT", "h2_analiz.json")
    out = ["| Level | Share of Tilt Alone | Share of Blur Alone |", "|:---|---:|---:|"]
    for k, name in (("dr8", "D/r0 = 8"), ("dr16", "D/r0 = 16")):
        out.append(f"| {name} | {H[k]['pay_tilt']:.3f} | {H[k]['pay_blur']:.3f} |")
    return out


def at_a_glance():
    """Headline strip at the top of the README: four numbers, each read from an analysis file."""
    F = J("results", "FAZ2", "analiz.json")["H1"]
    K = J("results", "KAPI", "dogrula_analiz.json")
    S = J("results", "GERCEK_VERI", "saha_ozet.json")["ort15"]["2021-07"]
    best = min(K["B"][f"KAPI/A0@{c}"]["medyan_oran"] for c in ("R0-8", "R5-8", "R5-4"))
    clean = K["T"]["KAPI/A0@R0-0"]
    return [f"| **{F['oran8_medyan']:.2f}×** | **−{100 * (1 - best):.0f}%** | **{clean['medyan_oran']:.3f}** | "
            f"**{S['G1_spearman']:.3f}** |",
            "|:---:|:---:|:---:|:---:|",
            "| ATE of the baseline at D/r0 = 8 relative to clear air (median paired ratio) | median ATE of the adaptive gate against the "
            "baseline (fast regime, D/r0 = 8) | ATE ratio in clear air: the gate costs nothing when it is not needed "
            f"(95% CI up to {clean['ga95'][1]:.3f}) | rank correlation of the online estimate with a scintillometer "
            "on real imagery |"]


TABLES = [("at a glance", at_a_glance), ("hero cell", hero_cell), ("highlights", highlights), ("degradation", degradation),
          ("real imagery", real), ("tilt/blur factorial", factorial)]

if __name__ == "__main__":
    if "--check" in sys.argv:
        readme = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
        bad = n = 0
        for name, fn in TABLES:
            for line in fn():
                n += 1
                if line not in readme:
                    bad += 1
                    print(f"MISSING in README ({name}): {line}")
        print(f"{n} table lines checked, {bad} missing")
        sys.exit(1 if bad else 0)
    for name, fn in TABLES:
        print(f"<!-- {name} -->")
        print("\n".join(fn()))
        print()
