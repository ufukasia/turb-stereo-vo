"""4Seasons -> ASL/EuRoC (mav0) dönüştürücü.

NEDEN: TartanAir'de taban çizgisi kurulamadı (öğrenilmiş dedektör sentetik dokuda
tekrarlanabilir nokta üretmiyor; bkz. evidence/TARTANAIR_NEDEN_CALISMIYOR.md) ve
TUM VI outdoors'ta stereo taban çizgisi çok kısa (fx*B=31, EuRoC 50.5) + el
kamerası 3.26 deg/kare dönüyor. 4Seasons ikisini de çözer:

    fx = 501.48, taban çizgisi = 0.3005 m  ->  fx*B = 150.7  (EuRoC'un 3 KATI)
    araç montajlı  ->  el titremesi yok
    RTK-GNSS yer gerçeği, TÜM yörünge boyunca (TUM VI'da yalnız %9-50)

ÖLÇÜLEN takip kalitesi (9-kare zincir, ALIKED+LightGlue):
    EuRoC MH01 %45.7 · 4Seasons %29.9 (30 Hz) · TUM VI outdoors1 %15.9 · TartanAir %0.3
30 Hz'in TAMAMI kullanılmalı: 10 Hz'e seyreltilince %10.5'e düşüyor.

GÖRÜNTÜLER: `undistorted_images/` zaten rektifiye ve bozulmasız (kalibrasyon
`undistorted_calib_*.txt`). Hattın kendi cv2.stereoRectify çağrısı D=0, R=I,
T=[-B,0,0] ile özdeşlik döndürür -> hat DEĞİŞMEDEN çalışır.

POZ: GNSSPoses.txt satırları  ts_ns, tx,ty,tz, qx,qy,qz,qw, scale, quality
olup ÖLÇÜLDÜ ki hareket yönü gövde çerçevesinde +z (%98.8) — yani poz zaten
GÖRÜNTÜ çerçevesindedir (x sağ, y aşağı, z ileri). Bu yüzden T_BS = I yazılır
ve ek bir eksen permütasyonu UYGULANMAZ.

ÖLÇEK: GNSSPoses'un 9. sütunu (`scale`, bu kayıtta 0.940013) VO çözümünü metrik
dünyaya taşıyan çarpandır; ham öteleme onunla ÇARPILARAK metrik hale getirilir.
Uygulanmazsa yörünge %6 kısa çıkar ve ATE sistematik olarak şişer.

    python -m turb.datasets.fourseasons --src /home/.../4seasons \
        --dst /home/.../sdata_4seasons --recs recording_2020-06-12_10-10-57
"""
from __future__ import annotations

import argparse
import os

import numpy as np


def _read_undistorted_calib(cal_dir):
    """undistorted_calib_0.txt + undistorted_calib_stereo.txt -> (fx,fy,cx,cy,W,H,B)."""
    with open(os.path.join(cal_dir, "undistorted_calib_0.txt")) as f:
        parts = f.readline().split()
    assert parts[0] == "Pinhole", parts[0]
    fx, fy, cx, cy = (float(x) for x in parts[1:5])
    with open(os.path.join(cal_dir, "undistorted_calib_0.txt")) as f:
        f.readline()
        W, H = (int(x) for x in f.readline().split())
    T = np.loadtxt(os.path.join(cal_dir, "undistorted_calib_stereo.txt"))
    baseline = abs(float(T[0, 3]))
    return fx, fy, cx, cy, W, H, baseline


def _sensor_yaml(path, name, T_BS, fx, fy, cx, cy, W, H, rate, T_cn_cnm1=None):
    lines = [
        "sensor_type: camera",
        f"comment: 4Seasons -> ASL donusturuldu ({name}); undistorted_images kullanildi, bozulma SIFIR.",
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


def convert(src_rec, cal_dir, dst_seq, max_frames=None, link=True):
    from scipy.spatial.transform import Rotation as Rot
    import shutil

    fx, fy, cx, cy, W, H, B = _read_undistorted_calib(cal_dir)
    c0 = os.path.join(src_rec, "undistorted_images", "cam0")
    c1 = os.path.join(src_rec, "undistorted_images", "cam1")
    if not (os.path.isdir(c0) and os.path.isdir(c1)):
        return None

    # Damgalar SAYISAL sıralanır (sözlüksel sıralama karışık basamakta sırayı bozar;
    # bkz. evidence/TARTANAIR_SIRALAMA_KUSURU.md).
    ts0 = sorted(int(f[:-4]) for f in os.listdir(c0) if f.endswith(".png"))
    have1 = {int(f[:-4]) for f in os.listdir(c1) if f.endswith(".png")}
    ts = [t for t in ts0 if t in have1]
    if max_frames:
        ts = ts[:max_frames]
    if len(ts) < 50:
        return None
    rate = int(round(1e9 / np.median(np.diff(ts))))

    for sub in ("mav0/cam0/data", "mav0/cam1/data", "mav0/state_groundtruth_estimate0"):
        os.makedirs(os.path.join(dst_seq, sub), exist_ok=True)

    for cam, srcdir in ((0, c0), (1, c1)):
        ddir = os.path.join(dst_seq, f"mav0/cam{cam}/data")
        rows = []
        for t in ts:
            name = f"{t}.png"
            dstp = os.path.join(ddir, name)
            if not os.path.exists(dstp):
                srcp = os.path.abspath(os.path.join(srcdir, name))
                if link:
                    os.symlink(srcp, dstp)
                else:
                    shutil.copy(srcp, dstp)
            rows.append(f"{t},{name}")
        with open(os.path.join(dst_seq, f"mav0/cam{cam}/data.csv"), "w") as f:
            f.write("#timestamp [ns],filename\n" + "\n".join(rows) + "\n")

    T_BS = np.eye(4)                       # gövde = rektifiye cam0 (GT zaten bu çerçevede)
    T_cn = np.eye(4); T_cn[0, 3] = -B      # cam1 <- cam0
    _sensor_yaml(os.path.join(dst_seq, "mav0/cam0/sensor.yaml"), "cam0", T_BS,
                 fx, fy, cx, cy, W, H, rate)
    T_BS1 = np.eye(4); T_BS1[0, 3] = B
    _sensor_yaml(os.path.join(dst_seq, "mav0/cam1/sensor.yaml"), "cam1", T_BS1,
                 fx, fy, cx, cy, W, H, rate, T_cn)

    # --- yer gerçeği ---
    g = np.loadtxt(os.path.join(src_rec, "GNSSPoses.txt"), delimiter=",", skiprows=1)
    gts = g[:, 0].astype(np.int64)
    gxyz = g[:, 1:4] * float(g[0, 8])      # metrik ölçek (bkz. modül başlığı)
    gq = g[:, 4:8]                          # qx qy qz qw
    lines = ["#timestamp,p_RS_R_x [m],p_RS_R_y [m],p_RS_R_z [m],"
             "q_RS_w [],q_RS_x [],q_RS_y [],q_RS_z []"]
    for i in range(len(gts)):
        p = gxyz[i]; q = gq[i]
        lines.append(f"{gts[i]},{p[0]:.9f},{p[1]:.9f},{p[2]:.9f},"
                     f"{q[3]:.9f},{q[0]:.9f},{q[1]:.9f},{q[2]:.9f}")
    with open(os.path.join(dst_seq, "mav0/state_groundtruth_estimate0/data.csv"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # kaç kare GT kapsamında?
    lo, hi = gts[0], gts[-1]
    cov = sum(1 for t in ts if lo <= t <= hi)
    return len(ts), rate, cov / len(ts), fx, B


def main():
    ap = argparse.ArgumentParser(description="4Seasons -> ASL/EuRoC donusturucu")
    ap.add_argument("--src", default="/data/4seasons")
    ap.add_argument("--dst", default="/data/sdata_4seasons")
    ap.add_argument("--recs", nargs="+", default=None)
    ap.add_argument("--max_frames", type=int, default=None)
    ap.add_argument("--copy", action="store_true")
    a = ap.parse_args()

    cal = os.path.join(a.src, "calibration")
    recs = a.recs or sorted(d for d in os.listdir(a.src) if d.startswith("recording_"))
    os.makedirs(a.dst, exist_ok=True)
    for r in recs:
        srcr = os.path.join(a.src, r)
        if not os.path.isdir(os.path.join(srcr, "undistorted_images")):
            print(f"- {r}: goruntu yok, atlandi")
            continue
        out = convert(srcr, cal, os.path.join(a.dst, r), a.max_frames, link=not a.copy)
        if out:
            n, rate, cov, fx, B = out
            print(f"+ {r}: {n} kare @{rate}Hz, GT kapsam %{cov*100:.1f}, "
                  f"fx={fx:.1f} B={B:.4f} fx*B={fx*B:.1f}")


if __name__ == "__main__":
    main()
