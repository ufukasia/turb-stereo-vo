#!/usr/bin/env python
"""Makale görseli: MH05 (tohum 3, R5, D/r0 8) hücresinin bir karesini ön-uçla aynı rektifikasyon ve aynı enjektörle
temiz ve bozulmuş olarak yazar.   python hero_goruntu.py CIKTI_DIZINI DOKUM KARE [KARE ...]"""
import os, pickle, sys, tempfile
from pathlib import Path
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
import cv2, numpy as np
import kapi_dogrula as kd
from turb.runner import _build_cfg
from turb.hooks import injector_from_cfg
from stereo_front_end import StereoMatcher

cikti, dump, kareler = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3:]]
os.makedirs(cikti, exist_ok=True)
cell = kd.hucre("A0", "MH05", 3, 8.0, "R5", "D1")
cfg = _build_cfg(cell, tempfile.mkdtemp())
print("turbulence:", cfg["turbulence"])
inj = injector_from_cfg(cfg)
m = StereoMatcher.__new__(StereoMatcher)
m.base_path, m.sequence = Path(cfg["dataset"]["base_path"]), cfg["dataset"]["sequence"]
m.STEREO_EXTRINSICS = str(cfg["dataset"].get("stereo_extrinsics", "gomulu"))
m.load_calibration(); m.compute_rectification(); m.load_timestamps()
with open(dump, "rb") as f:
    f_ts = [str(int(fr["timestamp"])) for fr in pickle.load(f)["frames"]]
for i in kareler:
    ts = f_ts[i]
    fn = m.cam0_df.loc[m.cam0_df["timestamp"] == ts, "filename"].values[0]
    fn1 = m.cam1_df.loc[m.cam1_df["timestamp"] == ts, "filename"].values[0]
    d = m.base_path / m.sequence / "mav0"
    r1, r2 = m.rectify_images(cv2.imread(str(d / "cam0/data" / fn), 0), cv2.imread(str(d / "cam1/data" / fn1), 0))
    t1, _ = inj.apply_stereo(r1.copy(), r2.copy(), i)
    u, v = np.meshgrid(np.arange(0, r1.shape[1], 8.0), np.arange(0, r1.shape[0], 8.0))
    dx = inj.gt_displacement(i, 0, u.ravel(), v.ravel())
    dx = np.asarray(dx).reshape(-1, 2) if np.asarray(dx).ndim == 2 else np.stack(dx, 1)
    print(i, ts, r1.shape, "tilt RMS/eksen %.3f px, en büyük |d| %.2f px" % (float(np.sqrt((dx ** 2).mean())), float(np.hypot(dx[:, 0], dx[:, 1]).max())),
          "fark RMS %.2f DN" % float(np.sqrt(((t1.astype(float) - r1) ** 2).mean())))
    cv2.imwrite(os.path.join(cikti, f"kare{i}_temiz.png"), r1)
    cv2.imwrite(os.path.join(cikti, f"kare{i}_turb.png"), t1)
    np.save(os.path.join(cikti, f"kare{i}_d.npy"), dx)
