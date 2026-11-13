"""H2 koşullarının kanca çıktısı: MH02 30 kare, tohum 0, sigma=1; 'yok'@D/r0=8 ile D/r0=0 bayt karşılaştırması
ve tilt/blur kanallarının gerçekten ayrı ayrı uygulandığı (fark > 0)."""
import sys, os, glob, hashlib, json
ROOT = "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb"
sys.path.insert(0, ROOT); os.chdir(ROOT)
import cv2, numpy as np
from turb.runner import Cell, _build_cfg
from turb import hooks
base = os.environ["SLAM_DATA_ROOT"]
fr = sorted(glob.glob(f"{base}/MH02/mav0/cam0/data/*.png"))[600:630]
def run(dr, tilt, blur):
    c = Cell(dataset="euroc", base_path=base, sequence="MH02", D_over_r0=dr, seed=0, tilt=tilt, blur=blur,
             windows=(600,), sensor_noise_dn=1.0)
    inj = hooks.injector_from_cfg(_build_cfg(c, "/tmp"))
    h = hashlib.sha1(); imgs = []
    for i, p in enumerate(fr):
        a = cv2.imread(p, 0); b = cv2.imread(p.replace("/cam0/", "/cam1/"), 0)
        a, b = inj.apply_stereo(a, b, 600 + i)
        h.update(a.tobytes()); h.update(b.tobytes()); imgs.append(a)
    return h.hexdigest(), np.array(imgs, dtype=np.float32)
res = {}
h0, I0 = run(0.0, True, "speckle_global")
for name, (t, b) in {"yok": (False, "none"), "tilt": (True, "none"), "blur": (False, "speckle_global"),
                     "ikisi": (True, "speckle_global")}.items():
    for dr in (8.0, 16.0):
        h, I = run(dr, t, b)
        res[f"{name}_dr{dr:g}"] = {"sha1": h, "D/r0=0 ile bayt-aynı": h == h0,
                                   "ortalama |fark| DN": float(np.abs(I - I0).mean())}
res["dr0_sha1"] = h0
print(json.dumps(res, indent=1, ensure_ascii=False))
