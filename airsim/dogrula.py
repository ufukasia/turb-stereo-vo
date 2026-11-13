"""
Kayıt doğrulaması — TartanAir/KITTI'de yaşanan kusurların AirSim'de tekrarlanmaması için.

D1 zaman: damgalar 19 hane, sözlüksel sıra == sayısal sıra, sabit dt
D2 stereo geometri: sol-sağ SIFT eşleşmelerinde |dy| < 1 px (rektifiye) ve
   disparite > 0 (cam1 gerçekten +x'te)
D3 ölçek: eşleşme disparitesi ile render derinliğinden beklenen fx*B/Z karşılaştırması
   (mesafe doğruluğunun doğrudan testi)
D4 poz: ardışık kareler arası GT göreli pozu ile görüntüden kestirilen esas matris
   (dönme hatası, öteleme yönü kosinüsü) — GT çerçeve sözleşmesinin testi

    python dogrula.py <sekans_dizini>
"""
import json
import os
import sys

import cv2
import numpy as np
import yaml


def oku(seq):
    m = os.path.join(seq, "mav0")
    c0 = yaml.safe_load(open(os.path.join(m, "cam0/sensor.yaml")))
    c1 = yaml.safe_load(open(os.path.join(m, "cam1/sensor.yaml")))
    fx, fy, cx, cy = c0["intrinsics"]
    B = -np.array(c1["T_cn_cnm1"]).reshape(4, 4)[0, 3]
    ts = [l.split(",")[0] for l in open(os.path.join(m, "cam0/data.csv")) if not l.startswith("#")]
    gt = np.loadtxt(os.path.join(m, "state_groundtruth_estimate0/data.csv"), delimiter=",", comments="#")
    return m, fx, cx, cy, B, ts, gt


def R_q(w, x, y, z):
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def main(seq, d4_aralik=5):
    m, fx, cx, cy, B, ts, gt = oku(seq)
    sonuc = {}
    # D1
    ti = [int(t) for t in ts]
    d = np.diff(ti)
    sonuc["D1_zaman"] = bool(all(len(t) == 19 for t in ts) and sorted(ts) == ts
                             and sorted(ti) == ti and d.min() == d.max())
    sift = cv2.SIFT_create(2000)
    bf = cv2.BFMatcher()
    K = np.array([[fx, 0, cx], [0, fx, cy], [0, 0, 1]])
    idx = np.linspace(20, len(ts) - 30, 8).astype(int)
    dy_all, disp_err, rot_err, cos_t = [], [], [], []
    for i in idx:
        L = cv2.imread(os.path.join(m, "cam0/data", ts[i] + ".png"), 0)
        R = cv2.imread(os.path.join(m, "cam1/data", ts[i] + ".png"), 0)
        Z = np.load(os.path.join(m, "depth0/data", ts[i] + ".npy")).astype(np.float32)
        kl, dl = sift.detectAndCompute(L, None)
        kr, dr = sift.detectAndCompute(R, None)
        if dl is None or dr is None:
            continue
        mm = [a for a, b in bf.knnMatch(dl, dr, k=2) if a.distance < 0.7 * b.distance]
        pl = np.array([kl[a.queryIdx].pt for a in mm]); pr = np.array([kr[a.trainIdx].pt for a in mm])
        if len(pl) < 20:
            continue
        dy = pr[:, 1] - pl[:, 1]; disp = pl[:, 0] - pr[:, 0]
        ok = np.abs(dy) < 2
        dy_all += list(dy[ok])
        sx = Z.shape[1] / L.shape[1]; sy = Z.shape[0] / L.shape[0]   # derinlik düşük çözünürlüklü olabilir
        u = np.clip((pl[ok, 0] * sx).astype(int), 0, Z.shape[1] - 1)
        v = np.clip((pl[ok, 1] * sy).astype(int), 0, Z.shape[0] - 1)
        z = Z[v, u]; good = (z > 0.5) & (z < 1e3)
        disp_err += list(disp[ok][good] - fx * B / z[good])
        # D4: kare i ile i+d4_aralik (varsayılan 5; tele/ileri harekette paralaks için büyütülür)
        j = min(i + d4_aralik, len(ts) - 1)
        L2 = cv2.imread(os.path.join(m, "cam0/data", ts[j] + ".png"), 0)
        k2, d2 = sift.detectAndCompute(L2, None)
        m2 = [a for a, b in bf.knnMatch(dl, d2, k=2) if a.distance < 0.7 * b.distance]
        if len(m2) < 30:
            continue
        p1 = np.array([kl[a.queryIdx].pt for a in m2]); p2 = np.array([k2[a.trainIdx].pt for a in m2])
        E, msk = cv2.findEssentialMat(p1, p2, K, cv2.RANSAC, 0.999, 1.0)
        _, Re, te, _ = cv2.recoverPose(E, p1, p2, K, mask=msk)
        g1, g2 = gt[i], gt[j]
        R1, R2 = R_q(*g1[4:8]), R_q(*g2[4:8])
        R12 = R2.T @ R1                          # kamera1 -> kamera2 (x2 = R12 x1 + t)
        t12 = R2.T @ (g1[1:4] - g2[1:4])
        ang = np.degrees(np.arccos(np.clip((np.trace(Re.T @ R12) - 1) / 2, -1, 1)))
        rot_err.append(ang)
        if np.linalg.norm(t12) > 0.05:
            cos_t.append(float(te.ravel() @ t12 / np.linalg.norm(t12)))
    dy_all = np.array(dy_all); disp_err = np.array(disp_err)
    sonuc["D2_dy_medyan_px"] = float(np.median(np.abs(dy_all))) if dy_all.size else None
    sonuc["D2_disparite_pozitif_oran"] = None
    sonuc["D3_disparite_hatasi_medyan_px"] = float(np.median(disp_err)) if disp_err.size else None
    sonuc["D3_disparite_hatasi_mad_px"] = float(np.median(np.abs(disp_err - np.median(disp_err)))) if disp_err.size else None
    sonuc["D4_donme_hatasi_medyan_deg"] = float(np.median(rot_err)) if rot_err else None
    sonuc["D4_oteleme_kosinus_medyan"] = float(np.median(cos_t)) if cos_t else None
    sonuc["GECTI"] = bool(sonuc["D1_zaman"]
                          and (sonuc["D2_dy_medyan_px"] or 9) < 0.5
                          and abs(sonuc["D3_disparite_hatasi_medyan_px"] or 9) < 0.3
                          and (sonuc["D4_donme_hatasi_medyan_deg"] or 9) < 0.3
                          and (sonuc["D4_oteleme_kosinus_medyan"] or 0) > 0.99)
    print(json.dumps(sonuc, indent=2))
    sonuc["D4_aralik"] = d4_aralik
    json.dump(sonuc, open(os.path.join(seq, f"dogrulama{'' if d4_aralik == 5 else '_d4_' + str(d4_aralik)}.json"), "w"), indent=2)
    return sonuc


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 5)
