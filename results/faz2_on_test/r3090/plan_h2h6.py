import json, sys, os, dataclasses
sys.path.insert(0, "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb"); os.chdir("/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb")
from turb.experiments import cells_FAZ2_faktoriyel, cells_FAZ2_h6
F = "/tmp/faz2_on_test/r3090"
cs = [c for c in cells_FAZ2_faktoriyel(seeds=1, seqs=["MH02"])] + \
     [c for c in cells_FAZ2_h6(12, seeds=1, seqs=["MH02"], varsayilan_dahil=True) if c.D_over_r0 == 12.0]
plan = []
for c in cs:
    d = dataclasses.asdict(c)
    d = {k: (list(v) if isinstance(v, tuple) else v) for k, v in d.items()}
    d["cfg_overrides"] = [list(x) for x in c.cfg_overrides]
    nm = c.key().split("/s0/")[0].replace("euroc/", "").replace("/", "_") + "_" + c.tag
    plan.append({"cell": d, "arms": [] if "--bos" in sys.argv else [["A0", [], "A0"]], "fe_dir": f"{F}/fe_h2h6/{nm}"})
json.dump(plan, open(sys.argv[1], "w"), indent=1)
print(len(plan), [p["fe_dir"].split("/")[-1][-60:] for p in plan])
