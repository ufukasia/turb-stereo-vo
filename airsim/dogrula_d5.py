"""
D5: derinlik + GT göreli poz ile yeniden izdüşüm testi (esas-matris belirsizliğinden bağımsız).
Kare i'deki SIFT noktası p1, derinlik Z(p1) ile 3B'ye kaldırılır, GT ile kare j'ye taşınır ve
izdüşürülür; eşleşen p2 ile fark (px). Dar FOV / düzleme yakın sahnede D4'ün yerine geçer.

    python dogrula_d5.py <sekans_dizini> [aralik=5]
Geçme: medyan hata < 1.0 px ve %90'lık < 3 px.
"""
import json, os, sys
import cv2
import numpy as np
from dogrula import oku, R_q


def main(seq, gap=5):
    m, fx, cx, cy, B, ts, gt = oku(seq)
    sift = cv2.SIFT_create(2000); bf = cv2.BFMatcher()
    err_all = []; per = []
    for i in np.linspace(20, len(ts) - gap - 5, 12).astype(int):
        j = i + gap
        L = cv2.imread(os.path.join(m, "cam0/data", ts[i] + ".png"), 0)
        L2 = cv2.imread(os.path.join(m, "cam0/data", ts[j] + ".png"), 0)
        Z = np.load(os.path.join(m, "depth0/data", ts[i] + ".npy")).astype(np.float32)
        kl, dl = sift.detectAndCompute(L, None); k2, d2 = sift.detectAndCompute(L2, None)
        if dl is None or d2 is None:
            continue
        mm = [a for a, b in bf.knnMatch(dl, d2, k=2) if a.distance < 0.7 * b.distance]
        if len(mm) < 20:
            continue
        p1 = np.array([kl[a.queryIdx].pt for a in mm]); p2 = np.array([k2[a.trainIdx].pt for a in mm])
        sx = Z.shape[1] / L.shape[1]; sy = Z.shape[0] / L.shape[0]
        u = np.clip((p1[:, 0] * sx).astype(int), 0, Z.shape[1] - 1)
        v = np.clip((p1[:, 1] * sy).astype(int), 0, Z.shape[0] - 1)
        z = Z[v, u]; ok = (z > 0.5) & (z < 1e3)
        # AirSim DepthPerspective: optik eksen boyunca değil ışın boyunca olabilir -> ikisi de denenir
        X = np.stack([(p1[:, 0] - cx) / fx, (p1[:, 1] - cy) / fx, np.ones(len(p1))], 1)
        R1, R2 = R_q(*gt[i][4:8]), R_q(*gt[j][4:8])
        R12 = R2.T @ R1; t12 = R2.T @ (gt[i][1:4] - gt[j][1:4])
        res = {}
        for ad, P in (("planar", X * z[:, None]), ("isin", X / np.linalg.norm(X, axis=1, keepdims=True) * z[:, None])):
            Q = P @ R12.T + t12
            q = np.stack([fx * Q[:, 0] / Q[:, 2] + cx, fx * Q[:, 1] / Q[:, 2] + cy], 1)
            res[ad] = np.linalg.norm(q - p2, axis=1)[ok]
        per.append({"i": int(i), "n": int(ok.sum()), "akis_px": float(np.median(np.linalg.norm(p2 - p1, axis=1))),
                    "med_planar": float(np.median(res["planar"])), "med_isin": float(np.median(res["isin"]))})
        err_all.append(res)
    e = {k: np.concatenate([r[k] for r in err_all]) for k in ("planar", "isin")}
    best = min(e, key=lambda k: np.median(e[k]))
    s = {"D5_aralik": gap, "D5_derinlik_kipi": best, "D5_med_px": float(np.median(e[best])),
         "D5_p90_px": float(np.percentile(e[best], 90)),
         "D5_med_px_diger": float(np.median(e["isin" if best == "planar" else "planar"])), "kareler": per}
    s["D5_GECTI"] = bool(s["D5_med_px"] < 1.0 and s["D5_p90_px"] < 3.0)
    json.dump(s, open(os.path.join(seq, f"dogrulama_d5_{gap}.json"), "w"), indent=2)
    print(json.dumps({k: v for k, v in s.items() if k != "kareler"}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 5)
