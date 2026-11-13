"""
Kuşbakışı yükseklik haritası (ComputerVision modu, "harita" kamerası dik aşağı, FOV 90).

Her karo için kamera H yükseklikte konumlanır; DepthPlanar = dikey uzaklık olduğundan
yüzey NED z = kamera_z + derinlik. Çıktı:
  DSM  : yüzey (ağaç tepeleri dahil)  — en yüksek nokta (en küçük NED z)
  DTM  : çıplak arazi ≈ DSM'nin morfolojik açılışı (pencere ~20 m) + yumuşatma
  agac : DSM - DTM > 1 m
  ortomozaik (renkli) — rota çizimi için

    python harita.py --x -100 400 --y -300 200 --H 60 --cikti harita_africa.npz
"""
import argparse
import json
import math
import os
import sys
import time

import cv2
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kayit import temizle, _quat_to_R  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", type=float, nargs=2, default=[-100, 400])
    ap.add_argument("--y", type=float, nargs=2, default=[-300, 200])
    ap.add_argument("--H", type=float, default=60.0)
    ap.add_argument("--res", type=float, default=0.5)
    ap.add_argument("--adim", type=float, default=90.0)
    ap.add_argument("--temizle", default="Water_Plane|Poacher|Crocodile|Hippopotamus|Rhinoceros")
    ap.add_argument("--cikti", required=True)
    a = ap.parse_args()
    sys.path.insert(0, os.path.join(HERE, "lib_ms"))
    import airsim
    c = airsim.VehicleClient(); c.confirmConnection()
    silinen = temizle(c, a.temizle) if a.temizle else []
    nx = int((a.x[1] - a.x[0]) / a.res) + 1; ny = int((a.y[1] - a.y[0]) / a.res) + 1
    dsm = np.full((nx, ny), np.nan, np.float32)
    rgb = np.zeros((nx, ny, 3), np.uint8)
    req = [airsim.ImageRequest("harita", airsim.ImageType.DepthPlanar, True, False),
           airsim.ImageRequest("harita", airsim.ImageType.Scene, False, False)]
    xs = np.arange(a.x[0] + a.adim / 2, a.x[1], a.adim); ys = np.arange(a.y[0] + a.adim / 2, a.y[1], a.adim)
    t0 = time.time()
    for X in xs:
        for Y in ys:
            c.simSetVehiclePose(airsim.Pose(airsim.Vector3r(float(X), float(Y), -a.H),
                                            airsim.to_quaternion(0, 0, 0)), True)
            time.sleep(1.0)                                # doku akışı otursun
            r = c.simGetImages(req)
            D = np.array(r[0].image_data_float, np.float32).reshape(r[0].height, r[0].width)
            im = np.frombuffer(r[1].image_data_uint8, np.uint8).reshape(r[1].height, r[1].width, -1)[:, :, :3]
            W, Hh = D.shape[1], D.shape[0]
            f = W / 2 / math.tan(math.radians(45))
            cp, cq = r[0].camera_position, r[0].camera_orientation
            R = _quat_to_R([cq.w_val, cq.x_val, cq.y_val, cq.z_val])   # dünya <- AirSim kamera (x optik)
            u, v = np.meshgrid(np.arange(W) + 0.5, np.arange(Hh) + 0.5)
            ok = D < 1e3
            # AirSim kamera ekseni: x optik, y sağ, z aşağı(görüntü)
            Pc = np.stack([D, (u - W / 2) / f * D, (v - Hh / 2) / f * D], -1)[ok]
            Pw = Pc @ R.T + np.array([cp.x_val, cp.y_val, cp.z_val])
            gi = np.round((Pw[:, 0] - a.x[0]) / a.res).astype(int)
            gj = np.round((Pw[:, 1] - a.y[0]) / a.res).astype(int)
            m = (gi >= 0) & (gi < nx) & (gj >= 0) & (gj < ny)
            gi, gj, z = gi[m], gj[m], Pw[m, 2]
            col = im[ok][m]
            cur = dsm[gi, gj]
            yaz = np.isnan(cur) | (z < cur)
            # aynı hücreye düşen çoklu noktada en yüksek (en küçük NED z) kalır
            order = np.argsort(-z[yaz])
            gi2, gj2, z2, c2 = gi[yaz][order], gj[yaz][order], z[yaz][order], col[yaz][order]
            dsm[gi2, gj2] = z2; rgb[gi2, gj2] = c2
        print(f"  x={X:.0f} tamam  {time.time() - t0:.0f} s", flush=True)
    # boşlukları doldur
    bos = np.isnan(dsm)
    if bos.any():
        dsm_f = cv2.inpaint(np.nan_to_num(dsm, nan=0).astype(np.float32), bos.astype(np.uint8), 3, cv2.INPAINT_NS)
        dsm = np.where(bos, dsm_f, dsm)
    k = int(20 / a.res) | 1
    ker = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
    # NED'de yükseklik = -z; açılış (önce erozyon) yüksekliğe uygulanır
    h = -dsm
    dtm_h = cv2.morphologyEx(h, cv2.MORPH_OPEN, ker)
    dtm_h = cv2.GaussianBlur(dtm_h, (0, 0), 5 / a.res)
    agac = (h - dtm_h) > 1.0
    np.savez_compressed(a.cikti, dsm=dsm, dtm=-dtm_h, agac=agac, rgb=rgb, x0=a.x[0], y0=a.y[0],
                        res=a.res, H=a.H, silinen=np.array(silinen))
    print(json.dumps({"boyut": [nx, ny], "bos_oran": float(bos.mean()),
                      "arazi_yukseklik_aralik_m": [float(dtm_h.min()), float(dtm_h.max())],
                      "agac_orani": float(agac.mean()),
                      "agac_yukseklik_p90_m": float(np.percentile((h - dtm_h)[agac], 90)) if agac.any() else 0}))


if __name__ == "__main__":
    main()
