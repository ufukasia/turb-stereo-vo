"""
Nokta başına yeniden-izdüşüm artığı dökümü (HA4: derinlik dilimine göre ek artık) — AirSim.

    python artik_dok.py --dump <landmarks_*.pkl> --sekans AF_F5_h5 --pencere 20:857 --cikti artik.npz
    python artik_dok.py ... --poz <traj.csv>     # ikincil: kestirilen yörüngeyle (A0 düzleştirilmiş)

YÖNTEM (2026-11-08; türbülanslı veriye bakılmadan sabitlendi, ON_TEST_KOLLAR.md §2):
- Girdi yalnız ön-uç dökümü (gözlemler: [lid, u, v, has_depth, fresh, X, Y, Z], rektifiye cam0) ve GT.
  Arka-uca, kollara ve PyEco kod parmak izi dosyalarına DOKUNMAZ.
- Birincil poz kaynağı GT (AirSim'de kusursuz): kamera pozu T_WC = T_WB · T_B_C (T_B_C = T_BS · R1ᵀ).
  Böylece artık, poz kestirim hatasını İÇERMEZ; yalnız ölçüm gürültüsü (anahtar nokta konumu + türbülans
  eğimi) + üçgenleme kalıntısı kalır.
- Her landmark (≥ MIN_GOZLEM gözlem) GT pozlarıyla üçgenlenir: ilk değer stereo noktaların ortanca dünya
  konumu (yoksa ışın-orta-nokta en küçük kareleri), sonra Huber (k = 1 px) ağırlıklı Levenberg–Marquardt,
  8 yineleme. Artık r = (u, v)_gözlem − π(T_WC⁻¹ X).
- Serbestlik düzeltmesi: r · sqrt(2n / (2n − 3)) (n gözlem, 3 parametre) — iz uzunluğu dağılımı temiz ve
  türbülanslı arasında farklı olabildiği için.
- Gözlem derinliği = üçgenlenmiş noktanın gözlem kamerasındaki z'si.
- Dilim özeti: σ_rob = 1.4826 · medyan(|r|) (u ve v birlikte), n, medyan derinlik.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
VERI = Path("/data/sdata_airsim")
MIN_GOZLEM = 5
HUBER_PX = 1.0
LM_YINELEME = 8
# Dilim tanımları. "model": uret_derinlik.DERIN_DILIM (modelin kendi dilimleri, ikincil).
# "esit8": BİRİNCİL. Temiz pilot verisinin (F5, tohum 0-2, W1+W2, C0=0) GT-üçgenlenmiş gözlem derinliklerinin
# eşit sayılı 8 dilimi, 5 m'ye yuvarlanmış kenarlar (geçerli aralık 25-400 m; havuz 302 764 gözlem; nicelikler
# 64.6/79.8/91.2/102.9/114.9/128.9/147.2 m; türbülanslı veriye bakılmadan, 2026-11-08, ON_TEST_KOLLAR.md §2).
# 25 m altı ve 400 m üstü (toplam < %0.3; yozlaşmış üçgenleme) analiz dışıdır.
DILIMLER = {
    "model": [0, 25, 50, 100, 200, 400, 1e9],
    "esit8": None,      # aşağıda ESIT8 ile doldurulur
}
ESIT8 = [25, 65, 80, 90, 105, 115, 130, 145, 400]
DILIMLER["esit8"] = ESIT8


def _R(q):
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def gt_kesit(sekans: str, pencere: str | None) -> np.ndarray:
    """uret_derinlik --kareler a:b ile AYNI kesit: kaynak GT satırları [a, b)."""
    g = np.loadtxt(VERI / sekans / "mav0/state_groundtruth_estimate0/data.csv", delimiter=",", comments="#")
    if pencere:
        a, b = (int(x) for x in pencere.replace("-", ":").split(":"))
        g = g[a:b]
    return g


def pozlar_gt(g: np.ndarray, T_B_C: np.ndarray) -> dict:
    """GT satırlarından {ts(int): T_WC}. GT sütunları: ts, p(3), q(w,x,y,z)."""
    out = {}
    for row in g:
        T = np.eye(4)
        T[:3, :3] = _R(row[4:8])
        T[:3, 3] = row[1:4]
        out[int(round(row[0]))] = T @ T_B_C
    return out


def pozlar_traj(traj_csv: str, T_B_C: np.ndarray) -> dict:
    """PyEco yörüngesi (ts, p, q_x, q_y, q_z, q_w) -> {ts: T_WC}."""
    e = np.loadtxt(traj_csv, delimiter=",", comments="#")
    out = {}
    for row in e:
        T = np.eye(4)
        x, y, z, w = row[4:8]
        T[:3, :3] = _R((w, x, y, z))
        T[:3, 3] = row[1:4]
        out[int(round(row[0]))] = T @ T_B_C
    return out


def _ucgenle(Ts, uv, K, X0):
    """Huber-IRLS Levenberg–Marquardt; Ts: (n,4,4) T_CW; uv: (n,2). Dönüş X (3,) ya da None."""
    fx, fy, cx, cy = K
    X = X0.astype(np.float64).copy()
    lam = 1e-3
    R, t = Ts[:, :3, :3], Ts[:, :3, 3]

    def art(Xv):
        Pc = np.einsum("nij,j->ni", R, Xv) + t
        z = Pc[:, 2]
        if np.any(z <= 1e-3):
            return None, None, None
        pr = np.stack([fx * Pc[:, 0] / z + cx, fy * Pc[:, 1] / z + cy], 1)
        return uv - pr, Pc, z

    r, Pc, z = art(X)
    if r is None:
        return None
    for _ in range(LM_YINELEME):
        a = np.abs(r).max(1)
        w = np.where(a <= HUBER_PX, 1.0, HUBER_PX / np.maximum(a, 1e-12))
        # d(pr)/d(Pc)·d(Pc)/dX = J_proj · R
        Jp = np.zeros((len(z), 2, 3))
        Jp[:, 0, 0] = fx / z; Jp[:, 0, 2] = -fx * Pc[:, 0] / z ** 2
        Jp[:, 1, 1] = fy / z; Jp[:, 1, 2] = -fy * Pc[:, 1] / z ** 2
        J = np.einsum("nij,njk->nik", Jp, R).reshape(-1, 3)
        ww = np.repeat(w, 2)
        rr = r.reshape(-1)
        H = (J * ww[:, None]).T @ J
        gvec = (J * ww[:, None]).T @ rr
        cost0 = float(np.sum(ww * rr ** 2))
        ok = False
        for _k in range(6):
            try:
                dx = np.linalg.solve(H + lam * np.diag(np.diag(H) + 1e-12), gvec)
            except np.linalg.LinAlgError:
                return None
            r2, Pc2, z2 = art(X + dx)
            if r2 is not None and float(np.sum(ww * r2.reshape(-1) ** 2)) < cost0:
                X, r, Pc, z = X + dx, r2, Pc2, z2
                lam = max(lam / 10, 1e-9)
                ok = True
                break
            lam *= 10
        if not ok or np.linalg.norm(dx) < 1e-9 * max(1.0, np.linalg.norm(X)):
            break
    return X


def artik_hesapla(dump_pkl: str, pozlar_fn, min_gozlem: int = MIN_GOZLEM) -> dict:
    """Bütün gözlemler için artık. pozlar_fn(T_B_C) -> {ts: T_WC}."""
    d = pickle.load(open(dump_pkl, "rb"))
    K = [float(x) for x in d["rect_K"]]
    T_B_C = np.asarray(d["T_BS"], float) @ np.block([[np.asarray(d["R1"], float).T, np.zeros((3, 1))],
                                                    [np.zeros((1, 3)), np.ones((1, 1))]])
    P = pozlar_fn(T_B_C)
    frames = d["frames"]
    rows = []
    for fi, fr in enumerate(frames):
        ts = int(fr["timestamp"])
        if ts not in P:
            continue
        o = fr["obs"]
        if len(o):
            rows.append(np.column_stack([np.full(len(o), fi), o]))
    O = np.concatenate(rows)                  # [fi, lid, u, v, hd, fresh, X, Y, Z]
    TWC = {fi: P[int(frames[fi]["timestamp"])] for fi in np.unique(O[:, 0]).astype(int)}
    TCW = {k: np.linalg.inv(v) for k, v in TWC.items()}
    order = np.lexsort((O[:, 0], O[:, 1]))
    O = O[order]
    lids, start, cnt = np.unique(O[:, 1], return_index=True, return_counts=True)
    fx, fy, cx, cy = K
    out = {k: [] for k in ("fi", "lid", "u", "v", "z", "ru", "rv", "n")}
    n_tri = n_fail = 0
    for lid, s0, c in zip(lids, start, cnt):
        if c < min_gozlem:
            continue
        B = O[s0:s0 + c]
        fis = B[:, 0].astype(int)
        Ts = np.stack([TCW[f] for f in fis])
        uv = B[:, 2:4]
        st = (B[:, 4] > 0.5) & np.all(np.isfinite(B[:, 6:9]), 1) & (B[:, 8] > 0)
        if st.any():
            Xw = np.stack([TWC[f][:3, :3] @ p + TWC[f][:3, 3] for f, p in zip(fis[st], B[st, 6:9])])
            X0 = np.median(Xw, 0)
        else:                                  # ışın orta-nokta en küçük kareleri
            A = np.zeros((3, 3)); bb = np.zeros(3)
            for f, (u, v) in zip(fis, uv):
                Rw, c_ = TWC[f][:3, :3], TWC[f][:3, 3]
                dv = Rw @ np.array([(u - cx) / fx, (v - cy) / fy, 1.0]); dv /= np.linalg.norm(dv)
                M = np.eye(3) - np.outer(dv, dv); A += M; bb += M @ c_
            try:
                X0 = np.linalg.solve(A, bb)
            except np.linalg.LinAlgError:
                n_fail += 1
                continue
        X = _ucgenle(Ts, uv, K, X0)
        if X is None:
            n_fail += 1
            continue
        Pc = np.einsum("nij,j->ni", Ts[:, :3, :3], X) + Ts[:, :3, 3]
        z = Pc[:, 2]
        pr = np.stack([fx * Pc[:, 0] / z + cx, fy * Pc[:, 1] / z + cy], 1)
        r = (uv - pr) * np.sqrt(2.0 * c / (2.0 * c - 3.0))
        n_tri += 1
        out["fi"].append(fis); out["lid"].append(np.full(c, lid)); out["u"].append(uv[:, 0]); out["v"].append(uv[:, 1])
        out["z"].append(z); out["ru"].append(r[:, 0]); out["rv"].append(r[:, 1]); out["n"].append(np.full(c, c))
    res = {k: np.concatenate(v).astype(np.float32 if k not in ("fi", "lid", "n") else np.int32)
           for k, v in out.items()}
    res["meta"] = dict(K=K, n_landmark_ucgen=n_tri, n_landmark_hata=n_fail, n_gozlem=int(len(res["z"])),
                       n_kare=len(frames), min_gozlem=min_gozlem, huber_px=HUBER_PX)
    return res


def dilim_ozet(z, ru, rv, kenar) -> list:
    """Dilim başına {a, b, n, z_med, sigma_rob, rms}. σ_rob = 1.4826·medyan(|r|), u ve v birlikte."""
    out = []
    for a, b in zip(kenar[:-1], kenar[1:]):
        m = np.isfinite(z) & (z >= a) & (z < b)
        n = int(m.sum())
        if n == 0:
            out.append(dict(a=a, b=b, n=0, z_med=None, sigma_rob=None, rms=None)); continue
        r = np.concatenate([ru[m], rv[m]]).astype(np.float64)
        out.append(dict(a=a, b=b, n=n, z_med=float(np.median(z[m])),
                        sigma_rob=float(1.4826 * np.median(np.abs(r))), rms=float(np.sqrt(np.mean(r ** 2)))))
    return out


def ozetle(res: dict) -> dict:
    return {"meta": res["meta"], "dilim": {ad: dilim_ozet(res["z"], res["ru"], res["rv"], k)
                                           for ad, k in DILIMLER.items()}}


def kaydet(res: dict, yol: str):
    np.savez_compressed(yol, **{k: v for k, v in res.items() if k != "meta"},
                        meta=np.array(json.dumps(res["meta"])))


def sha1(p, n=None):
    h = hashlib.sha1()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True)
    ap.add_argument("--sekans", required=True)
    ap.add_argument("--pencere", default=None, help="a:b (uret_derinlik --kareler ile aynı)")
    ap.add_argument("--poz", default="gt", help="gt | <yörünge.csv>")
    ap.add_argument("--cikti", default=None, help="npz (gözlem başına)")
    ap.add_argument("--kenar_oner", type=int, default=0, help="eşit sayılı N dilim kenarı öner (yalnız temiz veri!)")
    a = ap.parse_args()
    g = gt_kesit(a.sekans, a.pencere)
    fn = (lambda T: pozlar_gt(g, T)) if a.poz == "gt" else (lambda T: pozlar_traj(a.poz, T))
    res = artik_hesapla(a.dump, fn)
    if a.cikti:
        kaydet(res, a.cikti)
    oz = ozetle(res)
    if a.kenar_oner:
        q = np.percentile(res["z"], np.linspace(0, 100, a.kenar_oner + 1)[1:-1])
        oz["kenar_oner"] = [float(x) for x in q]
    print(json.dumps(oz))


if __name__ == "__main__":
    main()
