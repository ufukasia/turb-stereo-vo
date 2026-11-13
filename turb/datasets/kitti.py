"""KITTI odometri -> ASL/EuRoC (mav0) dönüştürücü.

NEDEN: EuRoC dışı üç veri kümesi (TartanAir, TUM VI, 4Seasons) hattın track
ömrünü 1.8-3.5'e düşürdü ve ATE %2.9-16'ya çıktı (bkz.
`evidence/VERI_KUMESI_DURUM_RAPORU.md`). KITTI'nin STEREO GEOMETRİSİ bu üçünden
belirgin biçimde güçlü:

    fx = 718.86 (seq00) / 707.09,  taban çizgisi = 0.5372 m  ->  fx*B = 386
    EuRoC 47.2 · TartanAir 80.0 · TUM VI 31.2 · 4Seasons 150.7 · KITTI 386

20 m'de disparite 19.3 px (4Seasons'ta 5.0 px) -> derinlik %2.6 belirsiz.

EKSEN KONVANSİYONU — ÖLÇÜLEREK DOĞRULANDI (varsayılmadı):
poses/XX.txt satırları 3x4 T_w_cam0 (satır-major) olup KITTI dokümantasyonuna
göre SOL GRİ KAMERA çerçevesindedir. Ölçüm: ardışık pozların yerel çerçevedeki
hareket yönü ort=[-0.001, -0.017, +0.995], z>0 oranı **%99.9**. Yani poz zaten
GÖRÜNTÜ çerçevesindedir (x sağ, y aşağı, z ileri) ve TartanAir'de gereken
NED->görüntü permütasyonu BURADA GEREKMEZ. T_BS = I yazılır.

GÖRÜNTÜLER: `image_0` (sol gri) ve `image_1` (sağ gri) ZATEN rektifiyedir
(P0/P1 aynı fx, fy, cy; P1 yalnız -fx*B kadar x ötelemesi taşır). Hattın kendi
cv2.stereoRectify çağrısı D=0, R=I, T=[-B,0,0] ile özdeşlik döndürür.

ZAMAN: times.txt saniye cinsindendir; nanosaniyeye çevrilir. Damgalar SAYISAL
sıralanır (sözlüksel sıralama karışık basamakta sırayı bozar; bkz.
`evidence/TARTANAIR_SIRALAMA_KUSURU.md`).

    python -m turb.datasets.kitti --src /home/.../kitti_raw/dataset \
        --dst /home/.../sdata_kitti --seqs 07
"""
from __future__ import annotations

import argparse
import os
import shutil

import numpy as np


def _read_calib(path):
    d = {}
    for line in open(path):
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        vals = [float(x) for x in v.split()]
        if len(vals) == 12:
            d[k.strip()] = np.array(vals).reshape(3, 4)
    P0, P1 = d["P0"], d["P1"]
    fx, fy = float(P0[0, 0]), float(P0[1, 1])
    cx, cy = float(P0[0, 2]), float(P0[1, 2])
    baseline = float(-P1[0, 3] / fx)
    return fx, fy, cx, cy, baseline


def _sensor_yaml(path, name, T_BS, fx, fy, cx, cy, W, H, rate, T_cn_cnm1=None):
    lines = [
        "sensor_type: camera",
        f"comment: KITTI odometri -> ASL donusturuldu ({name}); image_0/1 zaten rektifiye, bozulma SIFIR.",
        "T_BS:",
        "  cols: 4",
        "  rows: 4",
        "  data: [" + ", ".join(f"{v:.9f}" for v in np.asarray(T_BS).ravel()) + "]",
        f"rate_hz: {rate}",
        f"resolution: [{W}, {H}]",
        "camera_model: pinhole",
        f"intrinsics: [{fx}, {fy}, {cx}, {cy}]",
        "distortion_model: radial-tangential",
        "distortion_coefficients: [0.0, 0.0, 0.0, 0.0]",
    ]
    if T_cn_cnm1 is not None:
        lines.append("T_cn_cnm1: [" + ", ".join(f"{v:.9f}" for v in np.asarray(T_cn_cnm1).ravel()) + "]")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def convert(src_seq, poses_txt, dst_seq, max_frames=None, link=True):
    from scipy.spatial.transform import Rotation as Rot
    import cv2

    l_dir = os.path.join(src_seq, "image_0")
    r_dir = os.path.join(src_seq, "image_1")
    if not (os.path.isdir(l_dir) and os.path.isdir(r_dir)):
        return None
    fx, fy, cx, cy, B = _read_calib(os.path.join(src_seq, "calib.txt"))
    times = np.loadtxt(os.path.join(src_seq, "times.txt"))       # saniye

    names = sorted(f for f in os.listdir(l_dir) if f.endswith(".png"))
    have_r = {f for f in os.listdir(r_dir) if f.endswith(".png")}
    names = [f for f in names if f in have_r]
    n = min(len(names), len(times))
    if max_frames:
        n = min(n, max_frames)
    if n < 50:
        return None
    ts = [int(round(times[i] * 1e9)) + 1_000_000_000 for i in range(n)]
    # sayısal artan olmalı
    assert all(ts[i] < ts[i + 1] for i in range(n - 1)), "zaman damgaları artan değil"
    rate = int(round(1e9 / np.median(np.diff(ts))))

    im0 = cv2.imread(os.path.join(l_dir, names[0]), cv2.IMREAD_GRAYSCALE)
    H, W = im0.shape

    for sub in ("mav0/cam0/data", "mav0/cam1/data", "mav0/state_groundtruth_estimate0"):
        os.makedirs(os.path.join(dst_seq, sub), exist_ok=True)

    for cam, srcdir in ((0, l_dir), (1, r_dir)):
        ddir = os.path.join(dst_seq, f"mav0/cam{cam}/data")
        rows = []
        for i in range(n):
            out = f"{ts[i]}.png"
            dstp = os.path.join(ddir, out)
            if not os.path.exists(dstp):
                srcp = os.path.abspath(os.path.join(srcdir, names[i]))
                if link:
                    os.symlink(srcp, dstp)
                else:
                    shutil.copy(srcp, dstp)
            rows.append(f"{ts[i]},{out}")
        with open(os.path.join(dst_seq, f"mav0/cam{cam}/data.csv"), "w") as f:
            f.write("#timestamp [ns],filename\n" + "\n".join(rows) + "\n")

    T_BS = np.eye(4)                    # gövde = rektifiye cam0 (GT zaten bu çerçevede)
    T_cn = np.eye(4); T_cn[0, 3] = -B
    _sensor_yaml(os.path.join(dst_seq, "mav0/cam0/sensor.yaml"), "cam0", T_BS,
                 fx, fy, cx, cy, W, H, rate)
    T_BS1 = np.eye(4); T_BS1[0, 3] = B
    _sensor_yaml(os.path.join(dst_seq, "mav0/cam1/sensor.yaml"), "cam1", T_BS1,
                 fx, fy, cx, cy, W, H, rate, T_cn)

    P = np.loadtxt(poses_txt).reshape(-1, 3, 4)
    lines = ["#timestamp,p_RS_R_x [m],p_RS_R_y [m],p_RS_R_z [m],"
             "q_RS_w [],q_RS_x [],q_RS_y [],q_RS_z []"]
    for i in range(min(n, len(P))):
        R = P[i, :3, :3]; t = P[i, :3, 3]
        q = Rot.from_matrix(R).as_quat()           # xyzw
        lines.append(f"{ts[i]},{t[0]:.9f},{t[1]:.9f},{t[2]:.9f},"
                     f"{q[3]:.9f},{q[0]:.9f},{q[1]:.9f},{q[2]:.9f}")
    with open(os.path.join(dst_seq, "mav0/state_groundtruth_estimate0/data.csv"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return n, rate, fx, B, W, H


def main():
    ap = argparse.ArgumentParser(description="KITTI odometri -> ASL donusturucu")
    ap.add_argument("--src", default="/data/kitti_raw/dataset")
    ap.add_argument("--poses", default=None, help="poses dizini (varsayilan: <src>/../poses)")
    ap.add_argument("--dst", default="/data/sdata_kitti")
    ap.add_argument("--seqs", nargs="+", default=["07"])
    ap.add_argument("--max_frames", type=int, default=None)
    ap.add_argument("--copy", action="store_true")
    a = ap.parse_args()

    poses_dir = a.poses or os.path.join(os.path.dirname(a.src.rstrip("/")), "poses")
    os.makedirs(a.dst, exist_ok=True)
    for s in a.seqs:
        src = os.path.join(a.src, "sequences", s)
        pz = os.path.join(poses_dir, f"{s}.txt")
        if not os.path.isdir(src) or not os.path.isfile(pz):
            print(f"- {s}: eksik ({src} / {pz})")
            continue
        out = convert(src, pz, os.path.join(a.dst, f"kitti{s}"), a.max_frames, link=not a.copy)
        if out:
            n, rate, fx, B, W, H = out
            print(f"+ kitti{s}: {n} kare @{rate}Hz {W}x{H} fx={fx:.1f} B={B:.4f} fx*B={fx*B:.1f}")


if __name__ == "__main__":
    main()
