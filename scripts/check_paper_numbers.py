#!/usr/bin/env python
"""Re-read the headline numbers of the paper from the analysis files in results/ and compare.

Each entry: (label, value printed in the paper, file, path inside the JSON).  A value passes when the
stored number rounds to the printed one at the printed precision.  Nothing is typed in from the paper
except the printed value itself.   Usage:  python scripts/check_paper_numbers.py
"""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
K = "results/KAPI/dogrula_analiz.json"
F = "results/FAZ2/analiz.json"
R = "results/FAZ2_rakip_ek/rakip_ek_analiz.json"
CHECKS = [
    ("ATE growth of the studied pipeline at D/r0 = 8", "9.03", F, ["H1", "oran8_medyan"]),
    ("ORB-SLAM3 paired ratio at D/r0 = 8", "5.35", R, ["HR1", "rakipler", "orbslam3", "8", "medyan_oran"]),
    ("ORB-SLAM3 paired ratio at D/r0 = 16", "20.18", R, ["HR1", "rakipler", "orbslam3", "16", "medyan_oran"]),
    ("ORB-SLAM3 paired ratio at D/r0 = 24", "45.50", R, ["HR1", "rakipler", "orbslam3", "24", "medyan_oran"]),
    ("AirSLAM paired ratio at D/r0 = 8", "2.62", R, ["HR1", "rakipler", "airslam", "8", "medyan_oran"]),
    ("Blur share of the ATE increase at D/r0 = 16", "0.362", "results/FAZ2_FAKT/h2_analiz.json", ["dr16", "pay_blur"]),
    ("AirSim turbulence-free median position ANEES", "1.490", "results/AIRSIM/kollar/analiz.json", ["H3", "anees_poz_medyan"]),
    ("Real imagery: rank correlation with scintillometer", "0.863", "results/GERCEK_VERI/saha_ozet.json", ["ort15", "2021-07", "G1_spearman"]),
]
COND = ["R0-8", "R5-8", "R5-4", "R0-16"]
GATE = {
    ("B", "KAPI/A0"): ["0.549", "0.362", "0.613", None],
    ("I1", "KAPI/A0"): [None, None, None, "0.794"],
    ("I4", "A1g/A0"): ["0.412", "0.343", "0.662", "0.582"],
    ("I4", "A3S/A0"): ["1.026", "1.006", "0.937", "1.010"],
    ("I4", "A1/A0"): ["1.027", "1.030", "0.898", "1.039"],
    ("I2", "KAPI/A3S"): ["0.450", "0.368", "0.681", "0.813"],
    ("I3", "KAPI/A1g"): ["1.244", "1.236", "1.046", "1.423"],
    ("I5_poz", "KAPI/A0"): ["0.156", "0.075", "0.401", "0.239"],
}
for (fam, cmp_), vals in GATE.items():
    for c, v in zip(COND, vals):
        if v is not None:
            CHECKS.append((f"Gate grid {fam}: {cmp_} at {c}", v, K, [fam, f"{cmp_}@{c}", "medyan_oran"]))
CHECKS += [
    ("Gate grid, clean: KAPI/A0", "1.020", K, ["T", "KAPI/A0@R0-0", "medyan_oran"]),
    ("Gate grid, clean: A3S/A0", "0.991", K, ["T2_ATE", "A3S/A0@R0-0", "medyan_oran"]),
    ("Gate grid, clean: KAPI/A3S", "1.029", K, ["T2_ATE", "KAPI/A3S@R0-0", "medyan_oran"]),
    ("Gate grid, clean ANEES: KAPI/A0", "0.683", K, ["T2_ANEES_poz", "KAPI/A0@R0-0", "medyan_oran"]),
]
CI = [  # (label, printed lower, printed upper, file, path)
    ("CI KAPI/A0 at R0-8", "0.415", "0.638", K, ["B", "KAPI/A0@R0-8", "ga95"]),
    ("CI KAPI/A0 at R5-8", "0.347", "0.529", K, ["B", "KAPI/A0@R5-8", "ga95"]),
    ("CI KAPI/A0 at R5-4", "0.564", "0.712", K, ["B", "KAPI/A0@R5-4", "ga95"]),
    ("CI KAPI/A0 at R0-16", "0.752", "0.849", K, ["I1", "KAPI/A0@R0-16", "ga95"]),
    ("CI KAPI/A0 clean", "0.976", "1.047", K, ["T", "KAPI/A0@R0-0", "ga95"]),
]
S7 = "results/SAPMA7/s7_analiz.json"
for fam, cmp_, vals in (("HG1", "AGg/A0", ["0.542", "0.440", "0.672", "0.849"]),
                        ("HG2", "KAPI/AGg", ["0.890", "0.895", "0.952", "0.953"]),
                        ("S3", "KAPI/AGg", ["0.819", "0.743", "0.843", "0.794"])):
    for c, v in zip(COND, vals):
        CHECKS.append((f"Deviation 7 {fam}: {cmp_} at {c}", v, S7, [fam, f"{cmp_}@{c}", "medyan_oran"]))
CHECKS += [
    ("Deviation 7, clean: AGg/A0", "0.989", S7, ["HG3", "AGg/A0@R0-0", "medyan_oran"]),
    ("Deviation 7, clean: KAPI/AGg", "1.019", S7, ["T_ATE", "KAPI/AGg@R0-0", "medyan_oran"]),
    ("Deviation 7 HG2 Holm p at R0-8", "0.031", S7, ["HG2", "KAPI/AGg@R0-8", "p_holm"]),
    ("Deviation 7 HG2 Holm p at R5-8", "0.031", S7, ["HG2", "KAPI/AGg@R5-8", "p_holm"]),
    ("Deviation 7 HG2 Holm p at R5-4", "0.108", S7, ["HG2", "KAPI/AGg@R5-4", "p_holm"]),
    ("Deviation 7 HG2 Holm p at R0-16", "0.011", S7, ["HG2", "KAPI/AGg@R0-16", "p_holm"]),
]
CI += [
    ("CI AGg/A0 clean", "0.956", "1.023", S7, ["HG3", "AGg/A0@R0-0", "ga95"]),
    ("CI KAPI/AGg at R0-8", "0.842", "1.001", S7, ["HG2", "KAPI/AGg@R0-8", "ga95"]),
]
_cache = {}


def get(path_file, keys):
    J = _cache.setdefault(path_file, json.load(open(os.path.join(ROOT, path_file))))
    for k in keys:
        J = J[k]
    return J


def same(printed, stored):
    nd = len(printed.split(".")[1]) if "." in printed else 0
    return f"{stored:.{nd}f}" == printed


bad = 0
for label, printed, pf, keys in CHECKS:
    try:
        v = get(pf, keys)
        ok = same(printed, float(v))
    except (KeyError, TypeError, OSError) as e:
        v, ok = f"NOT FOUND ({e!r})", False
    bad += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {label:55s} paper {printed:>7s}  file {v}")
for label, lo, hi, pf, keys in CI:
    a, b = get(pf, keys)
    ok = same(lo, a) and same(hi, b)
    bad += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {label:55s} paper [{lo}, {hi}]  file [{a}, {b}]")
print(f"\n{len(CHECKS) + len(CI)} numbers checked, {bad} mismatches")
sys.exit(1 if bad else 0)
