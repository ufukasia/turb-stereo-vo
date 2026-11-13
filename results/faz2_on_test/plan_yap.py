import json, sys
F = "/tmp/faz2_on_test"
W = json.load(open("/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb/evidence/pencereler_faz2.json"))["sequences"]
BASE = "/media/masa24/2tb ana/sdata"
def cell(seq, dr, seed=0):
    return {"dataset": "euroc", "base_path": BASE, "sequence": seq, "seed": seed, "D_over_r0": float(dr),
            "windows": W[seq]["starts"], "window_len": 400, "nees": True, "cov_every": 10, "tag": "faz2_ontest"}
def fe(seq, dr, seed=0):
    return f"{F}/fe/{seq}_dr{dr:g}_s{seed}"
mode = sys.argv[1]
plan = []
if mode == "kalib":
    for s in ["V1_01_easy", "V1_02_medium", "V2_01_easy", "V2_02_medium"]:
        arms = [["A0", [], "A0"]] + [["A0C", ["--turb_cal_factor", str(k)], f"A0C_k{k}"] for k in (1, 2, 4, 8, 12, 16, 24)]
        plan.append({"cell": cell(s, 0), "arms": arms, "fe_dir": fe(s, 0)})
elif mode == "duman":
    k = sys.argv[2]
    for s in ["V1_01_easy", "MH02"]:
        for dr in (0, 8, 16):
            arms = [["A0", [], "A0"], ["A0C", ["--turb_cal_factor", k], f"A0C_k{k}"], ["A1", [], "A1"],
                    ["A3", [], "A3"], ["A5", [], "A5"], ["AG", [], "AG"],
                    ["A0C+A3", ["--turb_cal_factor", k], f"A0C+A3_k{k}"]]
            plan.append({"cell": cell(s, dr), "arms": arms, "fe_dir": fe(s, dr)})
json.dump(plan, open(sys.argv[-1], "w"), indent=1)
print(len(plan), "hücre", sum(len(p["arms"]) for p in plan), "kol işi")
