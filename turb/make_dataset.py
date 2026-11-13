"""
Turb-VIO türetilmiş veri kümesi üreteci (D2 / C1).

Bozulmuş sekansları DİSKE yazar; böylece
  (i)  benchmark açık kaynak olarak yayımlanabilir (rapor D2/D8),
  (ii) BAŞKA sistemler (ORB-SLAM3, VINS-Fusion, OpenVINS, pyslam) hattımıza
       dokunmadan aynı bozulmayla koşulabilir (E8-a genellenebilirlik).

Yazılan sekans REKTİFİYE düzlemdedir: türbülans alanı rektifiye görüntü düzleminde
tanımlıdır ve dışa aktarılan kalibrasyon (bozulmasız pinhole + saf taban çizgisi)
bunu birebir yansıtır. Böylece tüketen sistemin kendi rektifikasyonu bozulmayı
YENİDEN örneklemez ve GT eğim alanı geçerli kalır.

    python -m turb.make_dataset --sequence MH03 --D_over_r0 5 --seed 0 \
        --out /home/.../sdata_turb --n_frames 2000
"""
from __future__ import annotations

import argparse
import json
import os

import cv2
import numpy as np


def _sensor_yaml(path, T_BS, fx, fy, cx, cy, w, h, T_cn_cnm1=None, comment=""):
    lines = ["sensor_type: camera",
             f"comment: {comment}",
             "T_BS:", "  cols: 4", "  rows: 4",
             "  data: [" + ", ".join(f"{v:.9f}" for v in np.asarray(T_BS).ravel()) + "]",
             "rate_hz: 20",
             f"resolution: [{w}, {h}]",
             "camera_model: pinhole",
             f"intrinsics: [{fx:.9f}, {fy:.9f}, {cx:.9f}, {cy:.9f}]",
             "distortion_model: radial-tangential",
             "distortion_coefficients: [0.0, 0.0, 0.0, 0.0]"]
    if T_cn_cnm1 is not None:
        lines.append("T_cn_cnm1: [" + ", ".join(f"{v:.9f}" for v in
                                                np.asarray(T_cn_cnm1).ravel()) + "]")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def make(sequence, base_path, out_root, D_over_r0, seed, n_frames,
         theta_c_over_fov=0.15, tau0_over_dt=1.0, blur="speckle_global",
         scint=0.0, tag=None, save_gt_field=True, gt_field_every=25):
    import config
    from stereo_front_end import StereoMatcher
    from .simulator import TurbParams, TurbulenceSimulator

    cfg = config.load()
    f = cfg["frontend"]
    m = StereoMatcher(base_path=base_path, sequence=sequence,
                      matcher_type=f["matcher"], filter_method=f["filter"],
                      y_axis_threshold=f["y_axis_threshold"], extractor=f["extractor"],
                      max_track_length=0, max_keypoints=8,
                      detection_threshold=f["detection_threshold"])
    name = tag or (f"{sequence}_dr{D_over_r0:g}_tc{theta_c_over_fov:g}"
                   f"_t0{tau0_over_dt:g}_s{seed}")
    dst = os.path.join(out_root, name)
    for sub in ("mav0/cam0/data", "mav0/cam1/data",
                "mav0/state_groundtruth_estimate0", "turb"):
        os.makedirs(os.path.join(dst, sub), exist_ok=True)

    p = TurbParams(enabled=D_over_r0 > 0, seed=seed, D_over_r0=D_over_r0,
                   theta_c_over_fov=theta_c_over_fov, tau0_over_dt=tau0_over_dt,
                   blur=blur, scint_sigma_I2=scint)
    # SAYISAL sırala: sözlüksel sıralama, damga basamak sayısı karışık olan veri
    # kümelerinde (TartanAir 10/11 hane) kare sırasını bozar. Sabit-hane damgalarda
    # (EuRoC/TUM-VI, 19 hane) sonuç birebir aynıdır.
    ts_all = [str(t) for t in sorted({int(t) for t in m.cam0_df["timestamp"]})][:n_frames]
    sim = None
    rows0, rows1 = [], []
    fields = {}
    for i, ts in enumerate(ts_all):
        n0 = m.cam0_df.loc[m.cam0_df["timestamp"] == ts, "filename"].values
        n1 = m.cam1_df.loc[m.cam1_df["timestamp"] == ts, "filename"].values
        if len(n0) == 0 or len(n1) == 0:
            continue
        i0 = cv2.imread(str(m.base_path / sequence / "mav0/cam0/data" / n0[0]),
                        cv2.IMREAD_GRAYSCALE)
        i1 = cv2.imread(str(m.base_path / sequence / "mav0/cam1/data" / n1[0]),
                        cv2.IMREAD_GRAYSCALE)
        if i0 is None or i1 is None:
            continue
        r0, r1 = m.rectify_images(i0, i1)
        if sim is None:
            sim = TurbulenceSimulator(p, image_size=(r0.shape[1], r0.shape[0]))
        d0 = sim.apply(r0, i, cam=0)
        d1 = sim.apply(r1, i, cam=1)
        cv2.imwrite(os.path.join(dst, "mav0/cam0/data", f"{ts}.png"), d0)
        cv2.imwrite(os.path.join(dst, "mav0/cam1/data", f"{ts}.png"), d1)
        rows0.append(f"{ts},{ts}.png")
        rows1.append(f"{ts},{ts}.png")
        if save_gt_field and sim.active and (i % gt_field_every == 0):
            du, dv = sim.displacement(i, 0)
            fields[f"f{i}"] = np.stack([du[::8, ::8], dv[::8, ::8]]).astype(np.float32)

    for cam, rows in ((0, rows0), (1, rows1)):
        with open(os.path.join(dst, f"mav0/cam{cam}/data.csv"), "w") as f_:
            f_.write("#timestamp [ns],filename\n" + "\n".join(rows) + "\n")

    # rektifiye kalibrasyon
    fx, fy, cx, cy = m.P1[0, 0], m.P1[1, 1], m.P1[0, 2], m.P1[1, 2]
    h, w = m.map1x.shape
    T_BS0 = np.asarray(m.cam0_params["T_BS"]["data"], dtype=np.float64).reshape(4, 4)
    R1 = np.eye(4); R1[:3, :3] = m.R1
    T_BS_rect = T_BS0 @ np.linalg.inv(R1)
    base = float(-m.P2[0, 3] / m.P2[0, 0])
    T_cn = np.eye(4); T_cn[0, 3] = -base
    T_BS1 = T_BS_rect.copy(); T_BS1[:3, 3] += T_BS_rect[:3, :3] @ np.array([base, 0, 0])
    _sensor_yaml(os.path.join(dst, "mav0/cam0/sensor.yaml"), T_BS_rect, fx, fy, cx, cy,
                 w, h, comment=f"Turb-VIO rektifiye cam0 ({name})")
    _sensor_yaml(os.path.join(dst, "mav0/cam1/sensor.yaml"), T_BS1, fx, fy, cx, cy,
                 w, h, T_cn, comment=f"Turb-VIO rektifiye cam1 ({name})")

    # GT'yi olduğu gibi kopyala
    import shutil
    shutil.copy(os.path.join(base_path, sequence,
                             "mav0/state_groundtruth_estimate0/data.csv"),
                os.path.join(dst, "mav0/state_groundtruth_estimate0/data.csv"))

    meta = {"source": sequence, "n_frames": len(rows0), "params": p.__dict__,
            "rect_K": [float(fx), float(fy), float(cx), float(cy)],
            "baseline_m": base, "info": sim.info() if sim else {}}
    with open(os.path.join(dst, "turb", "meta.json"), "w") as f_:
        json.dump(meta, f_, indent=1, default=float)
    if fields:
        np.savez_compressed(os.path.join(dst, "turb", "gt_tilt_fields.npz"), **fields)
    print(f"✓ {dst}  ({len(rows0)} kare)")
    return dst


def main():
    ap = argparse.ArgumentParser(description="Turb-VIO veri kümesi üreteci")
    ap.add_argument("--sequence", default="MH03")
    ap.add_argument("--base_path", default="/data/sdata/")
    ap.add_argument("--out", default="/data/sdata_turb")
    ap.add_argument("--D_over_r0", type=float, nargs="+", default=[0.0, 3.0, 5.0, 8.0])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n_frames", type=int, default=2000)
    ap.add_argument("--theta_c", type=float, default=0.15)
    ap.add_argument("--tau0", type=float, default=1.0)
    ap.add_argument("--blur", default="speckle_global")
    a = ap.parse_args()
    for dr in a.D_over_r0:
        make(a.sequence, a.base_path, a.out, dr, a.seed, a.n_frames,
             a.theta_c, a.tau0, a.blur)


if __name__ == "__main__":
    main()
