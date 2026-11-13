"""Regresyon kapısı: eski varsayılan pencerelerle bir hücre run_cell ile (runner'ın kendi yolu)."""
import json, os, sys, time
ROOT = "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from turb.runner import Cell, run_cell, WINDOW_STARTS, _build_cfg, _arm_argv
import torch
c = Cell(dataset="euroc", base_path=os.environ["SLAM_DATA_ROOT"], sequence="MH01",
         D_over_r0=0.0, seed=0, arms=("A0",), nees=True, cov_every=10)
snap = {"WINDOW_STARTS": list(WINDOW_STARTS), "cell_windows": list(c.windows),
        "argv": _arm_argv(c, "A0", "/x/o.csv", "/x/c.npz", "/x/s.json"),
        "gpu": torch.cuda.get_device_name(0)}
cfg = _build_cfg(c, "/x"); snap["cfg"] = cfg
t0 = time.time()
recs = run_cell(c, verbose=False)
snap["t_total"] = time.time() - t0
out = sys.argv[1]
with open(out, "w") as f:
    json.dump({"snap": snap, "recs": recs}, f, default=lambda o: o.item() if hasattr(o, "item") else str(o), indent=0)
print("BITTI", snap["t_total"], [r.get("status") for r in recs])
