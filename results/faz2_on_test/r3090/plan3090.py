import json, sys
F = "/tmp/faz2_on_test/r3090"
W = json.load(open("/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb/evidence/pencereler_faz2_kapsam1m.json"))["sequences"]
BASE = "/media/masa24/2tb ana/sdata"
def cell(seq, dr, seed=0, ce=2):
    return {"dataset": "euroc", "base_path": BASE, "sequence": seq, "seed": seed, "D_over_r0": float(dr),
            "windows": W[seq]["starts"], "window_len": 400, "nees": True, "cov_every": ce, "tag": "faz2_3090",
            "sensor_noise_dn": 1.0, "psf_legacy_scale": False}
fe = lambda s, dr: f"{F}/fe/{s}_dr{dr:g}_s0_n1"
plan = []
for s in ["V1_01_easy", "V1_02_medium", "V2_01_easy", "V2_02_medium"]:
    arms = [["A0", [], "A0"]] + [["A0C", ["--turb_cal_factor", str(k)], f"A0C_k{k}"] for k in (1, 2, 4, 8, 12, 16, 24)]
    plan.append({"cell": cell(s, 0), "arms": arms, "fe_dir": fe(s, 0)})
    plan.append({"cell": cell(s, 0, ce=10), "arms": [["A0", [], "A0_cov10"]], "fe_dir": fe(s, 0)})
k = "12"
full = [["A0", [], "A0"], ["A0C", ["--turb_cal_factor", k], f"A0C_k{k}"], ["A1", [], "A1"], ["A3", [], "A3"],
        ["A5", [], "A5"], ["AG", [], "AG"], ["A0C+A3", ["--turb_cal_factor", k], f"A0C+A3_k{k}"]]
# yavaş kollar (A3, A5) önce kuyruğa girsin diye MH02 hücreleri başa alınır
plan = [{"cell": cell("MH02", dr), "arms": full, "fe_dir": fe("MH02", dr)} for dr in (8, 0)] + plan
json.dump(plan, open(f"{F}/plan3090.json", "w"), indent=1)
print(len(plan), "iş", sum(len(p["arms"]) for p in plan), "kol işi")
