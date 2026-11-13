import sys, os, glob, numpy as np, cv2
sys.path.insert(0, "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb")
from turb.hooks import sensor_noise, TurbulenceInjector, injector_from_cfg
from turb.uret_derinlik import gurultu_ekle
base = "/media/masa24/2tb ana/sdata"
p = sorted(glob.glob(f"{base}/MH02/mav0/cam0/data/*.png"))[600]
a = cv2.imread(p, cv2.IMREAD_GRAYSCALE); b = cv2.imread(p.replace("/cam0/", "/cam1/"), cv2.IMREAD_GRAYSCALE)
r = {}
r["uret_derinlik_ile_bayt_ayni"] = all(np.array_equal(sensor_noise(a, 1.0, s, i, c), gurultu_ekle(a, 1.0, s, i, c))
                                       for s in (0, 3) for i in (0, 77) for c in (0, 1))
d0 = sensor_noise(a, 1.0, 0, 600, 0).astype(float) - a
d1 = sensor_noise(a, 1.0, 1, 600, 0).astype(float) - a
dc = sensor_noise(b, 1.0, 0, 600, 1).astype(float) - b
dk = sensor_noise(a, 1.0, 0, 601, 0).astype(float) - a
m = (a > 5) & (a < 250)
r["std_DN"] = float(d0[m].std()); r["ortalama_DN"] = float(d0[m].mean())
r["tohum_korelasyonu"] = float(np.corrcoef(d0.ravel(), d1.ravel())[0, 1])
r["cam0_cam1_korelasyonu"] = float(np.corrcoef(d0.ravel(), dc.ravel())[0, 1])
r["kare_korelasyonu"] = float(np.corrcoef(d0.ravel(), dk.ravel())[0, 1])
r["sigma0_ayni_nesne"] = sensor_noise(a, 0.0, 0, 0, 0) is a
print(r)
