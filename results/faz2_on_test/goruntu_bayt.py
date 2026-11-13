"""Kanca çıktısının bayt özeti: injector_from_cfg ile 30 gerçek EuRoC karesi (cam0+cam1), D/r0 ∈ {0, 8}.
   python goruntu_bayt.py out.json [noise_dn] [noise_seed]"""
import sys, os, json, hashlib, glob
ROOT = "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb"
sys.path.insert(0, ROOT); os.chdir(ROOT)
import cv2, numpy as np
import config
from turb.runner import Cell, _build_cfg
out = sys.argv[1]
dn = float(sys.argv[2]) if len(sys.argv) > 2 else None
ns = int(sys.argv[3]) if len(sys.argv) > 3 else None
base = os.environ["SLAM_DATA_ROOT"]
seq = "MH02"
f0 = sorted(glob.glob(f"{base}/{seq}/mav0/cam0/data/*.png"))[600:630]
res = {}
for dr in (0.0, 8.0):
    kw = dict(dataset="euroc", base_path=base, sequence=seq, D_over_r0=dr, seed=0, windows=(600,))
    if dn is not None:
        kw["sensor_noise_dn"] = dn
        kw["noise_seed"] = ns
    c = Cell(**kw)
    cfg = _build_cfg(c, "/tmp")
    class M: turb = None
    m = M()
    _t = cfg.get("turbulence") or {}
    # config.apply_frontend'in türbülans kapısıyla AYNI mantık: kapı config.py'den okunur
    import inspect
    src = inspect.getsource(config.apply_frontend)
    from turb import hooks
    gate = ("sensor_noise_dn" in src and float(_t.get("sensor_noise_dn", 0) or 0) > 0) or \
        bool(_t.get("enabled", False)) or str(_t.get("restoration", "none")) != "none"
    inj = hooks.injector_from_cfg(cfg) if gate else None
    h = hashlib.sha1()
    for i, p in enumerate(f0):
        a = cv2.imread(p, cv2.IMREAD_GRAYSCALE)
        b = cv2.imread(p.replace("/cam0/", "/cam1/"), cv2.IMREAD_GRAYSCALE)
        if inj is not None:
            a, b = inj.apply_stereo(a, b, 600 + i)
        h.update(a.tobytes()); h.update(b.tobytes())
        if i == 0:
            res[f"dr{dr:g}_kare0_dtype"] = str(a.dtype)
    res[f"dr{dr:g}"] = h.hexdigest()
    res[f"dr{dr:g}_kanca"] = inj is not None
json.dump(res, open(out, "w"), indent=1); print(res)
