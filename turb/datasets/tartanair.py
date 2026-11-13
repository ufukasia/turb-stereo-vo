"""
TartanAir -> ASL/EuRoC (mav0) dönüştürücü.

NEDEN: rapor §5.1, TartanAir'i "fiziksel olarak makul rejim" (dış mekân, uzun menzil,
açık gökyüzü) veri kümesi olarak konumlandırır. Hat ise ASL/EuRoC dizin düzenini okur.
Bu modül TartanAir yörüngelerini o düzene çevirir; hattın kendisinde HİÇBİR değişiklik
gerekmez (T_cn_cnm1 + resolution anahtarları TUM VI için zaten destekleniyordu).

TartanAir kalibrasyonu (sabit, tüm ortamlar):
    640x480, fx = fy = 320.0, cx = 320.0, cy = 240.0, bozulma YOK
    stereo taban çizgisi 0.25 m (sol -> sağ, +x)

POZ KONVANSİYONU: pose_left.txt satırları  tx ty tz qx qy qz qw  olup kameranın
NED dünya çerçevesindeki pozudur ve kamera eksenleri de NED hizalıdır
(x ileri, y sağ, z aşağı). Görüntü çerçevesi (x sağ, y aşağı, z ileri) ile ilişki:

    R_ned_img = [[0, 0, 1],
                 [1, 0, 0],
                 [0, 1, 0]]        (img z -> ned x,  img x -> ned y,  img y -> ned z)

Dönüştürücü GT'yi GÖRÜNTÜ çerçevesinde yazar ve sensor.yaml'da T_BS = I verir;
böylece hattın "gövde" çerçevesi doğrudan rektifiye cam0'dır ve GT ile tutarlıdır.

ZAMAN: TartanAir zaman damgası içermez. Nominal 10 Hz atanır (dt = 0.1 s); türbülans
zaman ekseni (tau0/dt) bu dt ile tutarlı olacak şekilde Cell.frame_dt'ye geçirilir.

    python -m turb.datasets.tartanair --src /home/.../tartanair_raw \
        --dst /home/.../sdata_tartanair --envs westerndesert seasonsforest --min_frames 300
"""
from __future__ import annotations

import argparse
import os
import shutil

import numpy as np

FX = FY = 320.0
CX, CY = 320.0, 240.0
WIDTH, HEIGHT = 640, 480
BASELINE = 0.25
DT_NS = 100_000_000            # 10 Hz

R_NED_IMG = np.array([[0.0, 0.0, 1.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0]])


def _sensor_yaml(path, name, T_BS, T_cn_cnm1=None):
    lines = [
        f"sensor_type: camera",
        f"comment: TartanAir -> ASL dönüştürülmüş ({name})",
        "T_BS:",
        "  cols: 4",
        "  rows: 4",
        "  data: [" + ", ".join(f"{v:.9f}" for v in np.asarray(T_BS).ravel()) + "]",
        "rate_hz: 10",
        f"resolution: [{WIDTH}, {HEIGHT}]",
        "camera_model: pinhole",
        f"intrinsics: [{FX}, {FY}, {CX}, {CY}]",
        "distortion_model: radial-tangential",
        "distortion_coefficients: [0.0, 0.0, 0.0, 0.0]",
    ]
    if T_cn_cnm1 is not None:
        flat = ", ".join(f"{v:.9f}" for v in np.asarray(T_cn_cnm1).ravel())
        lines.append(f"T_cn_cnm1: [{flat}]")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def convert_trajectory(src_traj: str, dst_seq: str, link: bool = True):
    """Tek bir TartanAir yörüngesini (P00x) mav0 düzenine çevir."""
    left = os.path.join(src_traj, "image_left")
    right = os.path.join(src_traj, "image_right")
    pose_l = os.path.join(src_traj, "pose_left.txt")
    if not (os.path.isdir(left) and os.path.isdir(right) and os.path.isfile(pose_l)):
        return None
    imgs_l = sorted(f for f in os.listdir(left) if f.endswith(".png"))
    imgs_r = sorted(f for f in os.listdir(right) if f.endswith(".png"))
    poses = np.loadtxt(pose_l)
    n = min(len(imgs_l), len(imgs_r), len(poses))
    if n < 50:
        return None

    for sub in ("mav0/cam0/data", "mav0/cam1/data", "mav0/state_groundtruth_estimate0"):
        os.makedirs(os.path.join(dst_seq, sub), exist_ok=True)

    ts = [1_000_000_000 + i * DT_NS for i in range(n)]

    for cam, imgs, srcdir in ((0, imgs_l, left), (1, imgs_r, right)):
        ddir = os.path.join(dst_seq, f"mav0/cam{cam}/data")
        rows = []
        for i in range(n):
            name = f"{ts[i]}.png"
            dstp = os.path.join(ddir, name)
            if not os.path.exists(dstp):
                srcp = os.path.abspath(os.path.join(srcdir, imgs[i]))
                if link:
                    os.symlink(srcp, dstp)
                else:
                    shutil.copy(srcp, dstp)
            rows.append(f"{ts[i]},{name}")
        with open(os.path.join(dst_seq, f"mav0/cam{cam}/data.csv"), "w") as f:
            f.write("#timestamp [ns],filename\n" + "\n".join(rows) + "\n")

    # --- kalibrasyon ---
    T_BS = np.eye(4)                       # gövde = cam0 görüntü çerçevesi
    T_cn_cnm1 = np.eye(4)                  # cam1 <- cam0 : +x yönünde taban çizgisi
    T_cn_cnm1[0, 3] = -BASELINE
    _sensor_yaml(os.path.join(dst_seq, "mav0/cam0/sensor.yaml"), "cam0", T_BS)
    T_BS1 = np.eye(4); T_BS1[0, 3] = BASELINE
    _sensor_yaml(os.path.join(dst_seq, "mav0/cam1/sensor.yaml"), "cam1", T_BS1, T_cn_cnm1)

    # --- ground truth (görüntü çerçevesinde) ---
    from scipy.spatial.transform import Rotation as Rot
    lines = ["#timestamp,p_RS_R_x [m],p_RS_R_y [m],p_RS_R_z [m],"
             "q_RS_w [],q_RS_x [],q_RS_y [],q_RS_z []"]
    for i in range(n):
        t_ned = poses[i, :3]
        q_ned = poses[i, 3:7]                       # qx qy qz qw
        R_w_ned = Rot.from_quat(q_ned).as_matrix()
        R_w_img = R_w_ned @ R_NED_IMG
        q = Rot.from_matrix(R_w_img).as_quat()      # xyzw
        lines.append(f"{ts[i]},{t_ned[0]:.9f},{t_ned[1]:.9f},{t_ned[2]:.9f},"
                     f"{q[3]:.9f},{q[0]:.9f},{q[1]:.9f},{q[2]:.9f}")
    with open(os.path.join(dst_seq, "mav0/state_groundtruth_estimate0/data.csv"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return n


def main():
    ap = argparse.ArgumentParser(description="TartanAir -> ASL/EuRoC dönüştürücü")
    ap.add_argument("--src", default="/data/tartanair_raw")
    ap.add_argument("--dst", default="/data/sdata_tartanair")
    ap.add_argument("--envs", nargs="+", default=None)
    ap.add_argument("--difficulty", default="Easy")
    ap.add_argument("--min_frames", type=int, default=300)
    ap.add_argument("--max_trajs", type=int, default=6)
    ap.add_argument("--copy", action="store_true", help="sembolik bağ yerine kopyala")
    a = ap.parse_args()

    envs = a.envs or sorted(d for d in os.listdir(a.src)
                            if os.path.isdir(os.path.join(a.src, d)))
    os.makedirs(a.dst, exist_ok=True)
    total = 0
    for env in envs:
        base = os.path.join(a.src, env, a.difficulty)
        if not os.path.isdir(base):
            continue
        trajs = sorted(d for d in os.listdir(base) if d.startswith("P"))
        kept = 0
        for tr in trajs:
            if kept >= a.max_trajs:
                break
            src_traj = os.path.join(base, tr)
            nimg = len([f for f in os.listdir(os.path.join(src_traj, "image_left"))
                        if f.endswith(".png")]) if os.path.isdir(
                            os.path.join(src_traj, "image_left")) else 0
            if nimg < a.min_frames:
                continue
            name = f"{env}_{a.difficulty}_{tr}"
            n = convert_trajectory(src_traj, os.path.join(a.dst, name), link=not a.copy)
            if n:
                print(f"✓ {name}: {n} kare")
                kept += 1
                total += 1
    print(f"\nToplam {total} sekans -> {a.dst}")


if __name__ == "__main__":
    main()
