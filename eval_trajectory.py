"""
Backend trajectory'sini EuRoC ground-truth ile karşılaştırır: GT zaman damgalarına
interpolasyonla hizalama, ATE (RMSE) + eksen bazlı hata, X-Y / X-Z / 3D plot PNG.

Kullanım:
    python eval_trajectory.py --traj outputs/gtsam_traj_MH_05_difficult.csv \
        --gt MH_05_difficult/mav0/state_groundtruth_estimate0/data.csv
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # başsız ortamda da PNG üretebilsin
import matplotlib.pyplot as plt
from pathlib import Path


def load_gt(gt_csv):
    df = pd.read_csv(gt_csv)
    df.columns = [c.strip() for c in df.columns]
    t = df['#timestamp'].values.astype(np.float64)
    xyz = df[['p_RS_R_x [m]', 'p_RS_R_y [m]', 'p_RS_R_z [m]']].values.astype(np.float64)
    return t, xyz


def load_traj(traj_csv):
    df = pd.read_csv(traj_csv)
    df.columns = [c.strip() for c in df.columns]
    t = df['#timestamp [ns]'].values.astype(np.float64)
    xyz = df[['p_RS_R_x [m]', 'p_RS_R_y [m]', 'p_RS_R_z [m]']].values.astype(np.float64)
    return t, xyz


def umeyama(src, dst, with_scale):
    """dst ~ s R src + t (Umeyama 1991)."""
    mu_s, mu_d = src.mean(0), dst.mean(0)
    xs, xd = src - mu_s, dst - mu_d
    cov = xd.T @ xs / len(src)
    U, D, Vt = np.linalg.svd(cov)
    S = np.eye(3)
    if np.linalg.det(U) * np.linalg.det(Vt) < 0:
        S[2, 2] = -1
    R = U @ S @ Vt
    var_s = (xs ** 2).sum() / len(src)
    s = float(np.trace(np.diag(D) @ S) / var_s) if with_scale else 1.0
    t = mu_d - s * R @ mu_s
    return s, R, t


def align_and_ate(gt_t, gt_xyz, est_t, est_xyz):
    """Kestirim GT damgalarina enterpole edilir, SE(3) Umeyama ile hizalanir (olcek sabit); ATE RMSE.
    Kestirim ilk karenin govde cercevesinde oldugundan hizalama sarttir."""
    t_start = max(gt_t[0], est_t[0]); t_end = min(gt_t[-1], est_t[-1])
    valid = (gt_t >= t_start) & (gt_t <= t_end)
    gtt, gxyz = gt_t[valid], gt_xyz[valid]
    interp = np.stack([np.interp(gtt, est_t, est_xyz[:, k]) for k in range(3)], axis=1)
    _, R, t = umeyama(interp, gxyz, with_scale=False)
    interp = interp @ R.T + t
    err = gxyz - interp
    ate = float(np.sqrt(np.mean(np.sum(err**2, axis=1))))
    axis_rmse = np.sqrt(np.mean(err**2, axis=0))
    return ate, axis_rmse, gtt, gxyz, interp


def evaluate(traj_csv, gt_csv, out_png, dataset_name="", subtitle="Pure-Visual Stereo VIO vs Ground Truth"):
    """GT ile kıyasla, ATE + eksen RMSE hesapla, başlıklı plot PNG üret. Dönüş: (ate, axis_rmse, n_compare)."""
    gt_t, gt_xyz = load_gt(gt_csv)
    est_t, est_xyz = load_traj(traj_csv)
    ate, axis_rmse, gtt, gxyz, interp = align_and_ate(gt_t, gt_xyz, est_t, est_xyz)

    print(f"[{dataset_name}] ATE (RMSE, SE(3) hizali): {ate:.4f} m  | eksen "
          f"x={axis_rmse[0]:.3f} y={axis_rmse[1]:.3f} z={axis_rmse[2]:.3f} m | "
          f"{len(gtt)} kıyas noktası")

    seg = (gt_t >= max(gt_t[0], est_t[0])) & (gt_t <= min(gt_t[-1], est_t[-1]))
    gseg = gt_xyz[seg]

    fig = plt.figure(figsize=(16, 5.4))
    # Dataset adı + ATE başlıkta.
    fig.suptitle(f"{dataset_name}     |     ATE (RMSE) = {ate:.3f} m",
                 fontsize=16, fontweight='bold')
    fig.text(0.5, 0.925, f"{subtitle}   "
             f"(axis RMSE  x={axis_rmse[0]:.3f}  y={axis_rmse[1]:.3f}  z={axis_rmse[2]:.3f} m)",
             ha='center', fontsize=10, color='0.35')

    ax1 = fig.add_subplot(1, 3, 1)
    ax1.plot(gseg[:, 0], gseg[:, 1], 'r-', lw=2, label='Ground Truth')
    ax1.plot(est_xyz[:, 0], est_xyz[:, 1], 'b-', lw=1, label='Estimate')
    ax1.scatter([est_xyz[0, 0]], [est_xyz[0, 1]], c='g', s=40, zorder=5, label='start')
    ax1.set_xlabel('X [m]'); ax1.set_ylabel('Y [m]'); ax1.set_title('Top view (X-Y)')
    ax1.axis('equal'); ax1.grid(True); ax1.legend()
    # ATE'yi panelin içine de yaz (görselden bağımsız okunsun)
    ax1.text(0.03, 0.97, f"{dataset_name}\nATE = {ate:.3f} m", transform=ax1.transAxes,
             va='top', ha='left', fontsize=11, fontweight='bold',
             bbox=dict(boxstyle='round', fc='white', ec='0.5', alpha=0.85))

    ax2 = fig.add_subplot(1, 3, 2)
    ax2.plot(gseg[:, 0], gseg[:, 2], 'r-', lw=2, label='Ground Truth')
    ax2.plot(est_xyz[:, 0], est_xyz[:, 2], 'b-', lw=1, label='Estimate')
    ax2.set_xlabel('X [m]'); ax2.set_ylabel('Z [m]'); ax2.set_title('Side view (X-Z)')
    ax2.axis('equal'); ax2.grid(True); ax2.legend()

    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    ax3.plot(gseg[:, 0], gseg[:, 1], gseg[:, 2], 'r-', lw=2, label='Ground Truth')
    ax3.plot(est_xyz[:, 0], est_xyz[:, 1], est_xyz[:, 2], 'b-', lw=1, label='Estimate')
    ax3.set_xlabel('X'); ax3.set_ylabel('Y'); ax3.set_zlabel('Z'); ax3.set_title('3D')
    ax3.legend()

    plt.tight_layout(rect=[0, 0, 1, 0.9])
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f"[{dataset_name}] Plot kaydedildi: {out_png}")
    return ate, axis_rmse, len(gtt)


def _default_gt(sequence='MH01'):
    """GT yolunu $SLAM_DATA_ROOT'tan türet (config.py ile aynı değişken; yoksa 'data')."""
    root = os.environ.get('SLAM_DATA_ROOT', 'data')
    return str(Path(root) / sequence / 'mav0/state_groundtruth_estimate0/data.csv')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--traj', default='outputs/gtsam_traj_MH01.csv')
    ap.add_argument('--gt', default=_default_gt('MH01'),
                    help='Ground-truth CSV. Varsayılan $SLAM_DATA_ROOT/MH01/mav0/... yolundan türetilir.')
    ap.add_argument('--out', default='outputs/trajectory_eval.png')
    ap.add_argument('--dataset', default='MH01', help='Plot başlığına yazılacak dataset adı')
    args = ap.parse_args()
    evaluate(args.traj, args.gt, args.out, dataset_name=args.dataset)


if __name__ == '__main__':
    main()
