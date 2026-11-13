"""
Saf görsel GTSAM backend'i (IMU yok).

Front-end'in kare başına landmark dump'ını okur; PnP ile poz init, motion-only BA
ve iSAM2 / fixed-lag smoother ile trajectory kestirir. Optimizasyon VO world
frame'inde (X(0)=I) çalışır, sonunda GT'ye tek SE(3) ile oturtulup gövde
frame'inde trajectory CSV'si yazılır.
"""
import pickle
import os
import time
import threading
import math
import numpy as np
import cv2
import gtsam
import gtsam_unstable as gu  # BatchFixedLagSmoother için
from gtsam import (Pose3, Rot3, Point3, Cal3_S2, ISAM2, ISAM2Params,
                   NonlinearFactorGraph, Values)
from gtsam.symbol_shorthand import X, L
from scipy.spatial.transform import Rotation as Rot, Slerp
from scipy.spatial import cKDTree
import pandas as pd
from pathlib import Path
from collections import Counter

# Yöntem sayaçları; yalnız --audit verildiğinde basılır.
AUDIT = Counter()

# ----------------------------- [TÜRBÜLANS/E5] pencere-göreli NEES kovaryansı -----------------------------
# Yalnız --cov_out verilince kullanılır; kestirime dokunmaz (ayrı bir grafın doğrusallaştırılması).
_COND_LOG = []


class _CovSkip(Exception):
    """Kovaryansı tanımsız örnek (ör. gözlemsiz poz): hata değil, atlanır."""


def _pose_marginal_schur(gfg, keys, dims, n_pose_keys, jitter=1e-9):
    """keys[0] pozunun 6x6 marjinal kovaryansı; landmark'lar Schur tümleyeniyle elenir.

    Eski hattaki `_pose_marginal_schur` ile aynı matematik (GTSAM Marginals/CHOLESKY koşul sayısı
    yüzünden kullanılmaz; landmark bloğunun blok-köşegenliği VARSAYILMAZ, yama faktörleri
    landmark-landmark blokları yaratır). Fark yalnız uygulamada: PyEco'nun anahtar kare penceresinde
    binlerce landmark olduğu için yoğun Hessian (gfg.hessian) yerine doğrusallaştırılmış (beyazlatılmış,
    robust ağırlıklı) Jakobiyenlerden SEYREK H = A^T A kurulur ve Hll seyrek LU ile çözülür.
    Sıra: keys[:n_pose_keys] pozlar (ilki hedef poz), kalanlar landmark. Dönüş (6,6) ya da None.
    """
    import scipy.sparse as _sp
    import scipy.sparse.linalg as _spl
    col, off = {}, 0
    for k, d in zip(keys, dims):
        col[k] = (off, int(d))
        off += int(d)
    N = off
    np_dim = int(sum(dims[:n_pose_keys]))
    rows_l, cols_l, vals_l = [], [], []
    r0 = 0
    for i in range(gfg.size()):
        fac = gfg.at(i)
        if fac is None:
            continue
        A, _b = fac.jacobian()
        A = np.asarray(A, dtype=np.float64)
        m = A.shape[0]
        c = 0
        for k in fac.keys():
            o, d = col[k]
            blk = A[:, c:c + d]
            rr, cc = np.nonzero(blk)
            rows_l.append(rr + r0)
            cols_l.append(cc + o)
            vals_l.append(blk[rr, cc])
            c += d
        r0 += m
    if not rows_l or np_dim < 6:
        return None
    Am = _sp.csr_matrix((np.concatenate(vals_l), (np.concatenate(rows_l), np.concatenate(cols_l))),
                        shape=(r0, N))
    H = (Am.T @ Am).tocsc()
    Hpp = H[:np_dim, :np_dim].toarray()
    S = Hpp.copy()
    if N > np_dim:
        Hpl = H[:np_dim, np_dim:].toarray()
        Hll = (H[np_dim:, np_dim:] + _sp.identity(N - np_dim, format='csc') * jitter).tocsc()
        try:
            X_ = _spl.splu(Hll).solve(np.ascontiguousarray(Hpl.T))
        except (RuntimeError, ValueError):
            return None
        S -= Hpl @ X_
    S[np.diag_indices_from(S)] += jitter
    try:
        cov_p = np.linalg.inv(S)
    except np.linalg.LinAlgError:
        return None
    C6 = cov_p[:6, :6]
    if not np.all(np.isfinite(C6)):
        return None
    C6 = 0.5 * (C6 + C6.T)
    ev = np.linalg.eigvalsh(C6)
    if not np.all(ev > 0):
        return None            # kovaryans tanım gereği PD; belirsiz Schur sonucu NEES'e girmez
    _COND_LOG.append(float(ev[-1] / max(ev[0], 1e-300)))
    return C6


# ----------------------------- Ground truth -----------------------------
def load_groundtruth(gt_csv):
    """EuRoC state_groundtruth_estimate0/data.csv -> (t_ns, pos Nx3, quat_wxyz Nx4)."""
    df = pd.read_csv(gt_csv)
    df.columns = [c.strip() for c in df.columns]
    t = df['#timestamp'].values.astype(np.int64)
    pos = df[['p_RS_R_x [m]', 'p_RS_R_y [m]', 'p_RS_R_z [m]']].values.astype(np.float64)
    quat_wxyz = df[['q_RS_w []', 'q_RS_x []', 'q_RS_y []', 'q_RS_z []']].values.astype(np.float64)
    return t, pos, quat_wxyz


def gt_pose_at(t_query, gt_t, gt_pos, gt_quat_wxyz):
    """Sorgu zaman damgasında GT gövde pozunu (T_WB, 4x4) interpolasyonla döndür."""
    t_query = int(t_query)
    if t_query <= gt_t[0]:
        idx = 0
    elif t_query >= gt_t[-1]:
        idx = len(gt_t) - 1
    else:
        idx = np.searchsorted(gt_t, t_query)
    # pozisyon lineer, rotasyon slerp
    i0 = max(0, min(idx - 1, len(gt_t) - 2))
    i1 = i0 + 1
    t0, t1 = gt_t[i0], gt_t[i1]
    alpha = 0.0 if t1 == t0 else float(np.clip((t_query - t0) / (t1 - t0), 0.0, 1.0))
    pos = (1 - alpha) * gt_pos[i0] + alpha * gt_pos[i1]
    key_rots = Rot.from_quat(np.stack([
        gt_quat_wxyz[i0][[1, 2, 3, 0]],  # wxyz -> xyzw
        gt_quat_wxyz[i1][[1, 2, 3, 0]],
    ]))
    slerp = Slerp([0.0, 1.0], key_rots)
    R = slerp([alpha])[0].as_matrix()
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = pos
    return T


# ----------------------------- Yardımcılar -----------------------------
def T_body_rectcam0(T_BS, R1):
    """Gövde <- rektifiye cam0 statik dönüşümü (T_BS: gövde <- cam0, R1: rektifikasyon rotasyonu)."""
    T_orig_rect = np.eye(4)
    T_orig_rect[:3, :3] = np.asarray(R1).T
    return T_BS @ T_orig_rect


#: numpy motion-only BA ayarlari (main doldurur; None ise GTSAM yolu)
_MOBA_NP = None


def _so3_exp(w):
    th = float(np.linalg.norm(w))
    K = np.array([[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]])
    if th < 1e-10:
        return np.eye(3) + K, np.eye(3) + 0.5 * K
    a, b_ = np.sin(th) / th, (1 - np.cos(th)) / th ** 2
    c_ = (th - np.sin(th)) / th ** 3
    return np.eye(3) + a * K + b_ * K @ K, np.eye(3) + b_ * K + c_ * K @ K


def _moba_np(T_init, obs, lm_world, max_init_depth, skip_lids):
    """Tek pozlu motion-only BA, noktalar sabit (GTSAM yolundaki 1e-4 m onselin limiti): GTSAM ile ayni izdusum,
    Huber agirligi ve LM kurali (lambda 1e-5, x10, goreli/mutlak tolerans 1e-5, en cok 100 yineleme), Pose3 sag
    pertürbasyonu T*Exp(xi), xi=[w, v]. Tek degiskenli 6x6 sistem; ~200 noktalik GTSAM grafini kurmaktan ~10x hizli."""
    c = _MOBA_NP
    fx, fy, cx, cy, bl = c['fx'], c['fy'], c['cx'], c['cy'], c['b']
    P, M, S = [], [], []
    for r in obs:
        lid = int(r[0])
        if lid not in lm_world or lid in skip_lids:
            continue
        u, vv, Zd = float(r[1]), float(r[2]), r[7]
        st = r[4] > 0.5 and r[3] > 0.5 and np.isfinite(Zd) and 0 < Zd < max_init_depth
        P.append(lm_world[lid])
        M.append((u, (u - fx * bl / Zd) if st else 0.0, vv))
        S.append(st)
    used = len(P)
    if used < 6:
        return T_init, used
    P = np.asarray(P, dtype=np.float64); M = np.asarray(M); S = np.asarray(S, dtype=bool)
    sig = np.where(S, c['ss'], c['sm'])
    kh = np.where(S, c['ks'], c['km'])

    def _hesap(R, t, jac):
        pc = (P - t) @ R
        z = pc[:, 2]
        on = z > 0
        iz = np.where(on, 1.0 / np.where(on, z, 1.0), 0.0)
        u = fx * pc[:, 0] * iz + cx
        v = fy * pc[:, 1] * iz + cy
        ur = u - fx * bl * iz
        e = np.stack([u - M[:, 0], np.where(S, ur - M[:, 1], 0.0), v - M[:, 2]], axis=1)
        e[~on] = 2.0 * fx          # GTSAM: kamera arkasi (throwCheirality=False) sabit hata, Jakobiyen sifir
        e[~on & ~S, 1] = 0.0
        e = e / sig[:, None]
        n = np.linalg.norm(e, axis=1)
        w = np.where(n <= kh, 1.0, kh / np.maximum(n, 1e-12))
        w = np.where(kh > 0, w, 1.0)
        rho = np.where((kh > 0) & (n > kh), kh * (n - 0.5 * kh), 0.5 * n * n)
        if not jac:
            return float(rho.sum()), None, None
        x, y = pc[:, 0], pc[:, 1]
        Dp = np.zeros((len(P), 3, 3))
        Dp[:, 0, 0] = fx * iz; Dp[:, 0, 2] = -fx * x * iz * iz
        Dp[:, 1, 0] = fx * iz; Dp[:, 1, 2] = -fx * (x - bl) * iz * iz
        Dp[:, 2, 1] = fy * iz; Dp[:, 2, 2] = -fy * y * iz * iz
        Dp[~S, 1, :] = 0.0
        Dp[~on] = 0.0
        Jpc = np.zeros((len(P), 3, 6))
        Jpc[:, 0, 1] = -pc[:, 2]; Jpc[:, 0, 2] = pc[:, 1]
        Jpc[:, 1, 0] = pc[:, 2]; Jpc[:, 1, 2] = -pc[:, 0]
        Jpc[:, 2, 0] = -pc[:, 1]; Jpc[:, 2, 1] = pc[:, 0]
        Jpc[:, 0, 3] = Jpc[:, 1, 4] = Jpc[:, 2, 5] = -1.0
        J = np.einsum('nij,njk->nik', Dp, Jpc) / sig[:, None, None]
        H = np.einsum('n,nik,nil->kl', w, J, J)
        g = np.einsum('n,nik,ni->k', w, J, e)
        return float(rho.sum()), H, g

    R, t = T_init[:3, :3].copy(), T_init[:3, 3].copy()
    lam = 1e-5
    err, H, g = _hesap(R, t, True)
    for _ in range(100):
        while True:
            try:
                d = np.linalg.solve(H + lam * np.eye(6), -g)
            except np.linalg.LinAlgError:
                d = None
            if d is not None:
                Rd, V = _so3_exp(d[:3])
                R2, t2 = R @ Rd, t + R @ (V @ d[3:])
                e2 = _hesap(R2, t2, False)[0]
                if e2 <= err:
                    break
            lam *= 10.0
            if lam > 1e5:
                d = None
                break
        if d is None:
            break
        dE = err - e2
        R, t = R2, t2
        lam = max(lam / 10.0, 1e-30)
        err_eski, (err, H, g) = err, _hesap(R, t, True)
        if dE < 1e-5 or dE / max(err_eski, 1e-300) < 1e-5:
            break
    T = np.eye(4); T[:3, :3] = R; T[:3, 3] = t
    return T, used


def _solve_motion_only(T_init, obs, lm_world, K_gtsam, K_stereo, fx, baseline,
                       meas_noise, stereo_noise, max_init_depth, skip_lids=None):
    """Tek pozlu motion-only BA (landmark'lar sıkı prior ile sabit); skip_lids dışlanır.
    (T_opt 4x4, kullanılan gözlem sayısı) döndürür."""
    skip_lids = skip_lids or set()
    if _MOBA_NP is not None:
        return _moba_np(T_init, obs, lm_world, max_init_depth, skip_lids)
    g = NonlinearFactorGraph(); v = Values()
    v.insert(X(0), Pose3(T_init))
    fixed = gtsam.noiseModel.Isotropic.Sigma(3, 1e-4)
    used = 0
    for r in obs:
        lid = int(r[0])
        if lid not in lm_world or lid in skip_lids:
            continue
        Lk = L(lid)
        if not v.exists(Lk):
            p = lm_world[lid]
            v.insert(Lk, Point3(*p))
            g.add(gtsam.PriorFactorPoint3(Lk, Point3(*p), fixed))
        u, vv = float(r[1]), float(r[2]); Z = r[7]
        if r[4] > 0.5 and r[3] > 0.5 and np.isfinite(Z) and 0 < Z < max_init_depth:
            g.add(gtsam.GenericStereoFactor3D(
                gtsam.StereoPoint2(u, u - fx * baseline / Z, vv), stereo_noise, X(0), Lk, K_stereo))
        else:
            g.add(gtsam.GenericProjectionFactorCal3_S2(
                np.array([u, vv]), meas_noise, X(0), Lk, K_gtsam, False, False))
        used += 1
    if used < 6:
        return T_init, used
    try:
        return gtsam.LevenbergMarquardtOptimizer(g, v).optimize().atPose3(X(0)).matrix(), used
    except Exception:
        return T_init, used


def _pencere_disi(T_WC, obs, lm_world, fx, fy, cx, cy, r_px):
    """Baslangic pozunda izdusumu olcumden r_px'ten uzak (ya da kamera arkasi) gozlemlerin landmark ID'leri."""
    T_CW = np.linalg.inv(T_WC)
    out = set()
    for r in obs:
        lid = int(r[0])
        if lid not in lm_world:
            continue
        pw = lm_world[lid]
        pc = T_CW[:3, :3] @ np.asarray(pw) + T_CW[:3, 3]
        if pc[2] <= 1e-6:
            out.add(lid); continue
        du = fx * pc[0] / pc[2] + cx - float(r[1]); dv = fy * pc[1] / pc[2] + cy - float(r[2])
        if du * du + dv * dv > r_px * r_px:
            out.add(lid)
    return out


def _chi2_outliers(T, obs, lm_world, fx, fy, cx, cy, baseline, max_init_depth,
                   sigma_px, sigma_stereo, th_mono, th_stereo):
    """Poz T ile χ² reprojection hatası eşiği aşan landmark ID'lerini döndürür."""
    if _MOBA_NP is not None:
        # vektorlestirilmis esdegeri (numpy hizli yol)
        sec = [i for i in range(len(obs)) if int(obs[i, 0]) in lm_world]
        if not sec:
            return set()
        O = obs[sec]
        ids = O[:, 0].astype(np.int64)
        Pw = np.array([lm_world[int(l)] for l in ids], dtype=np.float64)
        Tc = np.linalg.inv(T)
        pc = Pw @ Tc[:3, :3].T + Tc[:3, 3]
        z = pc[:, 2]
        on = z > 1e-6
        zs = np.where(on, z, 1.0)
        up = fx * pc[:, 0] / zs + cx
        vp = fy * pc[:, 1] / zs + cy
        du = up - O[:, 1]; dv = vp - O[:, 2]
        Zm = O[:, 7]
        with np.errstate(divide='ignore', invalid='ignore'):
            st = (O[:, 4] > 0.5) & (O[:, 3] > 0.5) & np.isfinite(Zm) & (Zm > 0) & (Zm < max_init_depth)
            dr = (up - fx * baseline / zs) - (O[:, 1] - fx * baseline / np.where(st, Zm, 1.0))
        c2 = np.where(st, (du * du + dv * dv + dr * dr) / (sigma_stereo * sigma_stereo),
                      (du * du + dv * dv) / (sigma_px * sigma_px))
        kotu = ~on | np.where(st, c2 > th_stereo, c2 > th_mono)
        return set(int(l) for l in ids[kotu])
    T_CW = np.linalg.inv(T)
    out = set()
    for r in obs:
        lid = int(r[0])
        if lid not in lm_world:
            continue
        pw = lm_world[lid]
        pc = T_CW @ np.array([pw[0], pw[1], pw[2], 1.0])
        if pc[2] <= 1e-6:
            out.add(lid); continue  # kamera arkası
        u_pred = fx * pc[0] / pc[2] + cx
        v_pred = fy * pc[1] / pc[2] + cy
        du = u_pred - float(r[1]); dv = v_pred - float(r[2])
        Z = r[7]
        stereo_ok = r[4] > 0.5 and r[3] > 0.5 and np.isfinite(Z) and 0 < Z < max_init_depth
        if stereo_ok:
            uR_pred = u_pred - fx * baseline / pc[2]
            uR_meas = float(r[1]) - fx * baseline / Z
            dr = uR_pred - uR_meas
            chi2 = (du * du + dv * dv + dr * dr) / (sigma_stereo * sigma_stereo)
            if chi2 > th_stereo:
                out.add(lid)
        else:
            chi2 = (du * du + dv * dv) / (sigma_px * sigma_px)
            if chi2 > th_mono:
                out.add(lid)
    return out


def refine_pose_motion_only(T_init, obs, lm_world, K_gtsam, K_stereo, fx, baseline,
                            meas_noise, stereo_noise, max_init_depth, args=None,
                            fy=None, cx=None, cy=None, return_outliers=False, init_skip=None):
    """Motion-only BA ile poz rafinasyonu; chi2_reject açıksa iteratif χ² outlier elemesi.
    return_outliers=True ise (T_opt, gate) döner; gate = gevşek eşikle (chi2_gate_mult) grafa sokulmayacak landmark'lar."""
    # init_skip: baslangic pozuyla tutarsiz gozlemler (PnP-RANSAC disi ya da ongoru penceresi disi)
    # en bastan disarida; ORB-SLAM3 Tracking.cc:3700-3714 PoseOptimization'a yalniz RANSAC inlier'larini verir.
    skip = set(init_skip) if init_skip else set()
    AUDIT['moba_calls'] += 1
    T_opt, used = _solve_motion_only(T_init, obs, lm_world, K_gtsam, K_stereo, fx, baseline,
                                     meas_noise, stereo_noise, max_init_depth, skip_lids=skip)
    AUDIT['moba_obs_used'] += used
    if used < 6:
        AUDIT['moba_abort_lt6'] += 1
        return (T_init, set()) if return_outliers else T_init
    AUDIT['moba_solved'] += 1
    if args is not None and int(getattr(args, 'chi2_reject', 0)) and fy is not None:
        AUDIT['chi2_reject_frames'] += 1
        for _ in range(int(getattr(args, 'chi2_iters', 4))):
            AUDIT['chi2_reject_rounds'] += 1
            new_out = _chi2_outliers(
                T_opt, obs, lm_world, fx, fy, cx, cy, baseline, max_init_depth,
                float(args.px_noise), float(args.stereo_sigma),
                float(args.chi2_thresh_mono), float(args.chi2_thresh_stereo))
            if new_out <= skip:
                break
            AUDIT['chi2_reject_removed'] += len(new_out - skip)
            skip |= new_out
            T_ref, used = _solve_motion_only(T_init, obs, lm_world, K_gtsam, K_stereo, fx, baseline,
                                             meas_noise, stereo_noise, max_init_depth, skip_lids=skip)
            if used < 6:
                break  # fazla eleme; son geçerli çözümde kal
            T_opt = T_ref
    if not return_outliers:
        return T_opt
    # graf kapısı: gevşek eşikte yalnız gross outlier'lar
    gate = set()
    if args is not None and fy is not None:
        mult = float(getattr(args, 'chi2_gate_mult', 9.0))
        gate = _chi2_outliers(
            T_opt, obs, lm_world, fx, fy, cx, cy, baseline, max_init_depth,
            float(args.px_noise), float(args.stereo_sigma),
            float(args.chi2_thresh_mono) * mult, float(args.chi2_thresh_stereo) * mult)
        AUDIT['gate_calls'] += 1
        AUDIT['gate_lids'] += len(gate)
        if gate:
            AUDIT['gate_frames_nonempty'] += 1
    return T_opt, gate


def _triangulate_multiview(poses_T_WC, uvs, fx, fy, cx, cy, T_CW=None):
    """N-görüşlü DLT üçgenleme. (nokta, S[-2]/S[-1] koşullama oranı) döner; dejenere ise (None, 0).
    T_CW verilirse önceden hesaplanmış tersler kullanılır."""
    K = np.array([[fx, 0.0, cx], [0.0, fy, cy], [0.0, 0.0, 1.0]])
    Tc = np.linalg.inv(np.asarray(poses_T_WC)) if T_CW is None else np.asarray(T_CW)
    Kp = K @ Tc[:, :3, :4]                      # (n,3,4); gorus basina K @ T_CW[:3,:4] ile ayni
    uv = np.asarray(uvs, dtype=np.float64)
    A = np.empty((2 * len(uvs), 4))
    A[0::2] = uv[:, :1] * Kp[:, 2] - Kp[:, 0]
    A[1::2] = uv[:, 1:2] * Kp[:, 2] - Kp[:, 1]
    _, S, Vt = np.linalg.svd(A)
    Xh = Vt[-1]
    if abs(Xh[3]) < 1e-9:
        return None, 0.0
    cond = float(S[-2] / S[-1]) if S[-1] > 1e-12 else 1e12
    return Xh[:3] / Xh[3], cond


def _mv_parallax_deg(poses_T_WC, p_world):
    """Kamera merkezleri arasında p_world'e bakan en büyük paralaks açısı (derece).
    max arccos(cos_ij) = arccos(min cos_ij); Gram matrisi üç dış çarpımla kurulur."""
    centers = np.array([T[:3, 3] for T in poses_T_WC])
    rays = centers - np.asarray(p_world)
    norms = np.linalg.norm(rays, axis=1, keepdims=True)
    if np.any(norms < 1e-6):
        return 0.0
    rays_n = rays / norms
    if len(rays_n) < 2:
        return 0.0
    x, y, z = rays_n[:, 0], rays_n[:, 1], rays_n[:, 2]
    G = np.outer(x, x) + np.outer(y, y) + np.outer(z, z)
    np.fill_diagonal(G, 2.0)   # köşegen min'e girmesin
    c = float(np.clip(G.min(), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


_NEDENSEL = None


def main(args=None):
    global _NEDENSEL
    _NEDENSEL = {} if os.environ.get('PYECO_NEDENSEL_OUT') else None
    # parametreler config.yaml'dan (CLI override mümkün)
    if args is None:
        import config
        args = config.backend_args()
    _t0_backend = time.time()
    AUDIT.clear()
    _audit = int(getattr(args, 'audit', 0) or 0)
    _lm_ever = set()   # şimdiye kadar init edilmiş landmark'lar (re-init sayımı)

    # boru hatti kipinde kareler on uctan kuyrukla gelir (boru_hatti.py); aksi halde dokumden
    data = getattr(args, '_veri', None) or pickle.load(open(args.dump, 'rb'))
    fx, fy, cx, cy = data['rect_K']
    frames = data['frames']
    # değerlendirme penceresi frames[start_frame:max_frames]
    start_frame = int(getattr(args, 'start_frame', 0) or 0)
    frames = frames[start_frame:args.max_frames]
    if start_frame > 0:
        print(f"Değerlendirme penceresi: kare [{start_frame}:{args.max_frames}] "
              f"({len(frames)} kare; dump'tan dilimlendi)")
    T_BS = np.asarray(data['T_BS'], dtype=np.float64)
    R1 = np.asarray(data['R1'], dtype=np.float64)
    T_B_C = T_body_rectcam0(T_BS, R1)          # gövde <- rektifiye cam0
    T_C_B = np.linalg.inv(T_B_C)

    baseline = float(data.get('baseline', 0.11))
    # max_init_depth verilmezse ThDepth*b (ORB-SLAM2 makalesi III-A; ORB-SLAM3 TUM-VI.yaml:61 ThDepth=40)
    if getattr(args, 'max_init_depth', None) is None or float(args.max_init_depth) <= 0:
        args.max_init_depth = float(getattr(args, 'th_depth_baselines', 40.0)) * baseline
    print(f"max_init_depth = {args.max_init_depth:.3f} m (b={baseline:.5f} m)")

    # Disparite-bias düzeltmesi: d'=d+δ, 3D nokta ışın boyunca Z'/Z ile ölçeklenir.
    disp_bias = float(getattr(args, 'disparity_bias', 0.0) or 0.0)
    if disp_bias != 0.0:
        fxb = fx * baseline; n_corr = 0
        for _fr in frames:
            _obs = _fr['obs']
            for _r in _obs:
                _Z = _r[7]
                if _r[3] > 0.5 and np.isfinite(_Z) and _Z > 0:
                    _dp = fxb / _Z + disp_bias
                    if _dp > 1e-6:
                        _k = (fxb / _dp) / _Z
                        _r[5] *= _k; _r[6] *= _k; _r[7] *= _k; n_corr += 1
        print(f"Disparite-bias ölçek düzeltmesi: δ={disp_bias:+.3f}px → {n_corr} stereo gözlem düzeltildi")

    K_gtsam = Cal3_S2(fx, fy, 0.0, cx, cy)
    K_stereo = gtsam.Cal3_S2Stereo(fx, fy, 0.0, cx, cy, baseline)
    K_cv = np.array([[fx, 0, cx], [0, fy, cy], [0, 0, 1]], dtype=np.float64)

    # Gürültü modelleri
    meas_base = gtsam.noiseModel.Isotropic.Sigma(2, args.px_noise)
    # Huber esigi (beyazlatilmis birimde): ORB-SLAM3 Optimizer.cc:852-853 deltaMono=sqrt(5.991),
    # deltaStereo=sqrt(7.815). Eski ayar: mono 1.345, stereo robust degil (huber_stereo_k=0).
    meas_noise = gtsam.noiseModel.Robust.Create(
        gtsam.noiseModel.mEstimator.Huber.Create(float(getattr(args, 'huber_mono_k', 1.345))), meas_base)
    stereo_noise = gtsam.noiseModel.Isotropic.Sigma(3, args.stereo_sigma)
    if float(getattr(args, 'huber_stereo_k', 0.0) or 0.0) > 0:
        stereo_noise = gtsam.noiseModel.Robust.Create(
            gtsam.noiseModel.mEstimator.Huber.Create(float(args.huber_stereo_k)), stereo_noise)
    _gurultu_onbellek = {}
    global _MOBA_NP
    _MOBA_NP = None
    if int(getattr(args, 'moba_numpy', 0)) and not int(getattr(args, 'karanlik_sigma', 0)):
        _MOBA_NP = {'fx': fx, 'fy': fy, 'cx': cx, 'cy': cy, 'b': baseline,
                    'sm': float(args.px_noise), 'km': float(getattr(args, 'huber_mono_k', 1.345)),
                    'ss': float(args.stereo_sigma), 'ks': float(getattr(args, 'huber_stereo_k', 0.0) or 0.0)}

    def _kare_gurultu(fr):
        """Aday 3: S_d > tau karede piksel gurultusu m = 1 + (S_d - tau)/(1 - tau) katina cikar; aksi halde taban."""
        if not int(getattr(args, 'karanlik_sigma', 0)) or 's_d' not in fr:
            return meas_noise, stereo_noise
        tau = float(args.karanlik_tau)
        s_d = float(fr['s_d'])
        if s_d <= tau:
            return meas_noise, stereo_noise
        m = round(1.0 + (s_d - tau) / (1.0 - tau), 3)
        if m not in _gurultu_onbellek:
            mn = gtsam.noiseModel.Robust.Create(
                gtsam.noiseModel.mEstimator.Huber.Create(float(getattr(args, 'huber_mono_k', 1.345))),
                gtsam.noiseModel.Isotropic.Sigma(2, args.px_noise * m))
            sn = gtsam.noiseModel.Isotropic.Sigma(3, args.stereo_sigma * m)
            if float(getattr(args, 'huber_stereo_k', 0.0) or 0.0) > 0:
                sn = gtsam.noiseModel.Robust.Create(gtsam.noiseModel.mEstimator.Huber.Create(float(args.huber_stereo_k)), sn)
            _gurultu_onbellek[m] = (mn, sn)
        return _gurultu_onbellek[m]

    pose0_prior = gtsam.noiseModel.Diagonal.Sigmas(np.array([1e-4]*6))  # gauge anchor
    # gevşek landmark prior'ı: yalnız indeterminacy emniyeti, ölçeği kilitlemez
    lm_prior = gtsam.noiseModel.Isotropic.Sigma(3, 10.0)
    # mv-üçgenleme prior'ı daha sıkı; DLT çözümü gevşek prior'la indeterminant kalıyor
    mv_prior = gtsam.noiseModel.Isotropic.Sigma(3, float(getattr(args, 'mv_prior_sigma', 0.3)))
    # sabit-hız modeli gürültüsü [rot x3, trans x3]; gevşek, yalnız regularizasyon
    motion_noise = gtsam.noiseModel.Diagonal.Sigmas(
        np.array([args.motion_rot_sigma]*3 + [args.motion_trans_sigma]*3))

    # Hız-orantılı hareket prior'ı (motion_sigma_mode='adaptive'):
    #   σ_trans(f) = clip(motion_sigma_rel · ||t(prev_rel)||, σ_min, motion_trans_sigma)
    motion_sigma_mode = str(getattr(args, 'motion_sigma_mode', 'fixed') or 'fixed')
    motion_sigma_rel = float(getattr(args, 'motion_sigma_rel', 0.5) or 0.5)
    motion_sigma_min = float(getattr(args, 'motion_sigma_min', 0.02) or 0.02)

    def _motion_noise_for(rel_pose):
        """Karenin CV factor gürültü modeli (adaptive modda hıza bağlı)."""
        if motion_sigma_mode != 'adaptive':
            return motion_noise
        speed_m = float(np.linalg.norm(np.asarray(rel_pose)[:3, 3]))
        sig_t = min(max(motion_sigma_rel * speed_m, motion_sigma_min),
                    float(args.motion_trans_sigma))
        sig_r = float(args.motion_rot_sigma)
        if int(getattr(args, 'motion_rot_uyarlamali', 0)):
            # donme gurultusu de hizla olceklenir (otelemeyle ayni katsayi); hizli donuste sabit-hiz onseli gevser
            ang = float(np.linalg.norm(Rot.from_matrix(np.asarray(rel_pose)[:3, :3]).as_rotvec()))
            sig_r = max(sig_r, motion_sigma_rel * ang)
        return gtsam.noiseModel.Diagonal.Sigmas(
            np.array([sig_r]*3 + [sig_t]*3))

    params = ISAM2Params()
    params.setRelinearizeThreshold(float(getattr(args, 'relinearize_threshold', 0.01)))
    params.relinearizeSkip = args.relinearize_skip
    # keyframe modu: yalnız keyframe'ler graf node'u olur, ara kareler PnP ile pozlanır
    kf_mode = int(getattr(args, 'keyframe_mode', 0)) > 0
    # Fixed-lag smoother: fixed_lag>0 (kare penceresi) veya keyframe + lba_window>0 (KF penceresi).
    # Zaman damgası kf_seq; keyframe kapalıyken kf_seq==f.
    use_fl = (int(getattr(args, 'fixed_lag', 0)) > 0) or (kf_mode and int(getattr(args, 'lba_window', 0)) > 0)
    if use_fl:
        win = int(args.fixed_lag) if int(getattr(args, 'fixed_lag', 0)) > 0 else int(args.lba_window)
        # Batch smoother + LM erken durma (sıcak başlangıç). IncrementalFixedLagSmoother bu GTSAM sürümünde kararsız.
        lm_params = gtsam.LevenbergMarquardtParams()
        lm_params.setMaxIterations(int(getattr(args, 'fl_lm_max_iters', 5)))
        lm_params.setRelativeErrorTol(float(getattr(args, 'fl_lm_rel_err_tol', 1e-3)))
        # GTSAM varsayilani 1e-5: zayif kisitli noktalar yuzunden kucuk sonumde sistem tekil, LM sonumu buyutene kadar
        # her guncellemede birkac eleme bosa gider
        lm_params.setlambdaInitial(float(getattr(args, 'fl_lm_lambda0', 1e-5)))
        _kayan = str(getattr(args, 'smoother_tipi', 'batch')) == 'kayan'

        def _kayan_yeni():
            from kayan_pencere import KayanPencere
            return KayanPencere(win, lm_params, sabit_ek=float(getattr(args, 'kayan_sabit_ek', 0)),
                                sabit_n=int(getattr(args, 'kayan_sabit_n', 1)),
                                kf_aralik=int(getattr(args, 'kayan_kf_aralik', 0)),
                                yakin=float(getattr(args, 'kayan_yakin', 10)),
                                dusen_oncel=bool(int(getattr(args, 'kayan_dusen_oncel', 1))))
        if str(getattr(args, 'smoother_tipi', 'batch')) == 'isam2':
            # ISAM2 tabanli artimli sabit gecikmeli duzlestirici: ayni pencere, yalniz degisen degiskenler yeniden
            # dogrusallastirilir (GTSAM IncrementalFixedLagSmoother)
            _ip = ISAM2Params()
            _ip.setRelinearizeThreshold(float(getattr(args, 'relinearize_threshold', 0.01)))
            _ip.relinearizeSkip = int(getattr(args, 'relinearize_skip', 1))
            smoother = gu.IncrementalFixedLagSmoother(float(win), _ip)
            print("Duzlestirici: ISAM2 tabanli artimli (IncrementalFixedLagSmoother)")
        elif _kayan:
            smoother = _kayan_yeni()
            print(f"Duzlestirici: marjinalizasyonsuz kayan pencere (sabit ek {args.kayan_sabit_ek:g}, "
                  f"sabit n {args.kayan_sabit_n})")
        else:
            smoother = gu.BatchFixedLagSmoother(float(win), lm_params)
        unit = "kare" if not kf_mode else "keyframe"
        print(f"FIXED-LAG/LOCAL-BA smoother AÇIK: pencere={win} {unit} "
              f"(batch, LM erken-durma: maxIter={int(getattr(args, 'fl_lm_max_iters', 5))}, "
              f"relTol={float(getattr(args, 'fl_lm_rel_err_tol', 1e-3)):g})")
    else:
        isam = ISAM2(params)
    if kf_mode:
        print(f"KEYFRAME modu AÇIK: min_trans={args.kf_min_trans}m, min_rot={args.kf_min_rot}rad, "
              f"tracked_ratio<{args.kf_max_tracked_ratio}, aralık[{args.kf_min_interval},{args.kf_max_interval}], "
              f"lba_window={args.lba_window}")
    REFRESH_EVERY = args.refresh_every  # PnP landmark cache'i kaç karede bir tazelenir
    if kf_mode and float(getattr(args, 'kf_max_interval_s', 0) or 0) > 0 and len(frames) > 1:
        # en uzun anahtar kare araligi saniyeden, sensorun kare hizindan (ORB3 Tracking.cc:3156 c1a:
        # mMaxFrames = fps, yani en az saniyede bir anahtar kare)
        _ts = np.array([int(fr_['timestamp']) for fr_ in frames], dtype=np.int64)
        _dt = float(np.median(np.diff(_ts))) * 1e-9
        args.kf_max_interval = max(1, int(round(float(args.kf_max_interval_s) / max(_dt, 1e-6))))
        print(f"kf_max_interval = {args.kf_max_interval} kare ({args.kf_max_interval_s} s, dt={_dt*1e3:.1f} ms)")

    # Loop closure (varsayılan kapalı); frontend.export_keypoints gerektirir.
    loop_closer = None
    loop_warned = [False]
    loop_stride = 1 if kf_mode else 3   # tam-graf modunda her kareyi DB'ye atma
    if int(getattr(args, 'loop_enable', 0)):
        from loop_closing import LoopCloser
        loop_closer = LoopCloser(args, fx, fy, cx, cy, baseline)
        print("LOOP CLOSURE AÇIK: " + ("PGO" if int(args.loop_pgo) else "PGO-yok") +
              (" + GBA" if int(args.loop_gba) else "") +
              f" | appearance_sim>{args.loop_appearance_sim}, min_inliers={args.loop_min_inliers}, "
              f"min_gap={args.loop_min_keyframe_gap} KF")

    initialized_lms = set()   # grafa eklenmiş landmark'lar (tam-graf modu)
    live_lms = set()          # fixed-lag: şu an pencerede olan landmark'lar
    live_poses = set()        # pencerede duran poz düğümleri (marjinalize olana factor eklenmez)
    lm_world = {}             # id -> world 3D (PnP cache'i, periyodik tazelenir)
    lm_desc = {}              # id -> temsili descriptor (search-by-projection)
    pending_mv = {}           # çok-görüşlü üçgenleme tamponu: lid -> [(f, T_WC, u, v), ...]
    n_mv_init = 0
    # landmark culling; bayraklar kapalıyken culler.active=False ve kancalar no-op
    from landmark_culling import LandmarkCuller
    culler = LandmarkCuller(args, fx, fy, cx, cy)
    resident_lms = set()       # fixed-lag: smoother'da fiilen duran landmark'lar
    n_culled_total = 0
    if culler.active:
        print(f"LANDMARK CULLING AÇIK: {culler.describe()}")
    pose_estimates = {}       # frame_idx -> T_W0_C (4x4)
    frame_ts = []             # frame_idx -> timestamp (int)
    prev_pose = np.eye(4)
    prev_rel = np.eye(4)      # son kare-kare bağıl hareket (sabit-hız modeli)
    # keyframe defteri: kf_seq commit edilen node sayacı, prev_graph_* son commit edilen node
    kf_seq = 0
    prev_graph_f = 0
    prev_graph_pose = np.eye(4)
    last_kf_pose = np.eye(4)
    ref_tracked = 0            # son KF'de izlenen harita noktası sayısı
    kf_flags = {}              # frame_idx -> keyframe mi
    # ara-kare pozları referans KF'ye göreli saklanır, çıktıda rafine KF'den kurulur
    inter_ref = {}             # ara_frame -> referans keyframe
    inter_rel = {}             # ara_frame -> inv(ham_KF_pozu) @ T_WC
    keyframes = []             # commit edilen keyframe'ler (sıralı)
    # Iraksama koruması: ardışık ya da yaygın update hatasında RuntimeError ile temiz çık.
    consec_err = 0
    total_err = 0
    DIVERGE_CONSEC = 20
    MAX_ERR_RATIO = 0.05
    MIN_FRAMES_FOR_RATIO = 200
    err_print_cap = 5
    # Sağlamlık kapıları (0 = kapalı): pnp_min_inliers, jump_max_trans (öngörüden sıçrama sınırı),
    # chi2_gate_graph (χ² outlier'larını grafa sokma).
    min_inl = int(getattr(args, 'pnp_min_inliers', 0) or 0)
    jump_gate = float(getattr(args, 'jump_max_trans', 0.0) or 0.0)
    jump_pnp_max = float(getattr(args, 'jump_pnp_max', 0.0) or 0.0)
    chi2_gate = int(getattr(args, 'chi2_gate_graph', 0) or 0) > 0
    # Sıçrama kapısı kilitlenmesin: güçlü PnP en fazla 1 kare, her durumda en fazla
    # JUMP_MAX_CONSEC kare reddedilir; sonra ölçüme geri dönülür (re-sync).
    JUMP_MAX_CONSEC = 5
    PNP_STRONG_RATIO = 0.6      # inlier oranı bunun üstünde = güçlü konsensüs
    jump_rejects = 0
    jump_gate_ref = str(getattr(args, 'jump_gate_ref', 'prediction') or 'prediction')
    resync_velocity = str(getattr(args, 'resync_velocity', 'zero') or 'zero')
    last_good_rel = np.eye(4)   # kapıdan geçmiş son bağıl hareket
    n_recover = 0               # smoother yeniden kurulma sayısı
    debug_frames = int(getattr(args, 'debug_frames', 0) or 0)
    # çok-görüşlü init kipi: 'prior' | 'factors'
    mv_init_mode = str(getattr(args, 'mv_init_mode', 'prior') or 'prior')
    bootstrap_depth = int(getattr(args, 'bootstrap_from_depth', 0) or 0) > 0
    pnp_min_frac = float(getattr(args, 'pnp_min_inlier_frac', 0.0) or 0.0)
    if min_inl > 0 or jump_gate > 0 or chi2_gate:
        print(f"SAĞLAMLIK KAPILARI: pnp_min_inliers={min_inl}, jump_max_trans={jump_gate}m, "
              f"chi2_gate_graph={'AÇIK' if chi2_gate else 'kapalı'}")

    # [Asama B] kalici yerel harita + yerel harita takibi (ORB-SLAM3 Tracking.cc TrackLocalMap / SearchLocalPoints).
    # Pencereden cikan nokta unutulmaz: son kestirim konumu, son tanimlayicisi ve gozlem sayisi saklanir. Her karede
    # takip pozundan sonra bu karede gozlenmeyen harita noktalari izdusurulur, tanimlayiciyla eslenir ve poz bu
    # eslesmelerle yeniden rafine edilir; pencere disindaki eslesen noktalar haritadaki konumunda onselle grafa doner.
    yh_on = int(getattr(args, 'yerel_harita', 0)) > 0
    harita_p, harita_d, harita_n = {}, {}, {}
    harita_son = {}            # lid -> son goruldugu kare (yerel harita = son fixed_lag karede gorulenler)
    yh_prior = gtsam.noiseModel.Isotropic.Sigma(3, float(getattr(args, 'yh_prior_sigma', 0.1)))
    # Zayif landmark onseli: kisitsiz kalan nokta dogrusal sistemi tekil yapmasin (g2o/LM sonumlemesinin karsiligi);
    # sigma buyuk oldugundan kestirime katkisi ihmal edilir. 0 = kapali.
    _zo = float(getattr(args, 'lm_zayif_onsel_sigma', 0.0) or 0.0)
    zayif_onsel = gtsam.noiseModel.Isotropic.Sigma(3, _zo) if _zo > 0 else None
    if yh_on:
        print(f"YEREL HARITA TAKIBI ACIK: yaricap {args.yh_yaricap_px} px (kayipta {args.yh_yaricap_kayip_px}), "
              f"tanimlayici {args.yh_desc_dist}, en az gozlem {args.yh_min_obs}, onsel {args.yh_prior_sigma} m")

    def _yh_kare(fr_, obs_):
        """Karenin anahtar noktalari (birim tanimlayicilarla) + gozlenen izlerin tanimlayicilarini haritada tazeler."""
        if fr_.get('kps') is None or fr_.get('desc') is None or len(fr_['kps']) == 0:
            return None
        kps_ = np.asarray(fr_['kps'], dtype=np.float64).reshape(-1, 2)
        desc_ = np.asarray(fr_['desc'], dtype=np.float64)
        desc_ = desc_ / np.maximum(np.linalg.norm(desc_, axis=1, keepdims=True), 1e-12)
        tree_ = cKDTree(kps_)
        used_ = set()
        if len(obs_):
            dd, jj = tree_.query(obs_[:, 1:3], k=1, distance_upper_bound=2.0)
            for r_, d_, j_ in zip(obs_, dd, jj):
                if np.isfinite(d_):
                    harita_d[int(r_[0])] = desc_[int(j_)]
                    used_.add(int(j_))
        return kps_, desc_, tree_, used_

    def _yh_esle(T_WC_, kare_, obs_ids_, r_px):
        """Gozlenmeyen harita noktalarini T_WC_ pozuna izdusur, yaricap icindeki en benzer serbest anahtar noktayla esle.
        Her anahtar nokta en fazla bir harita noktasina (en yuksek benzerlik). Donus: {lid: (u, v)}."""
        kps_, desc_, tree_, used_ = kare_
        lids_ = [l for l, n_ in harita_n.items()
                 if n_ >= int(args.yh_min_obs) and l not in obs_ids_ and l in harita_d and l in harita_p]
        if not lids_:
            return {}
        P_ = np.array([harita_p[l] for l in lids_])
        T_CW_ = np.linalg.inv(T_WC_)
        pc_ = P_ @ T_CW_[:3, :3].T + T_CW_[:3, 3]
        z_ = pc_[:, 2]
        ok_ = z_ > float(baseline)                       # ORB3 minZ = taban (Frame.cc:841)
        u_ = np.where(ok_, fx * pc_[:, 0] / np.where(ok_, z_, 1.0) + cx, -1.0)
        v_ = np.where(ok_, fy * pc_[:, 1] / np.where(ok_, z_, 1.0) + cy, -1.0)
        ok_ &= (u_ >= 0) & (u_ < 2 * cx) & (v_ >= 0) & (v_ < 2 * cy)
        idx_ = np.where(ok_)[0]
        AUDIT['yh_izdusen'] += len(idx_)
        if len(idx_) == 0:
            return {}
        komsu = tree_.query_ball_point(np.stack([u_[idx_], v_[idx_]], axis=1), r=float(r_px))
        en_iyi = {}
        thr = 1.0 - float(args.yh_desc_dist)
        nn_oran = float(getattr(args, 'yh_nn_orani', 0.0))
        for i_, cand_ in zip(idx_, komsu):
            cand_ = [c for c in cand_ if c not in used_]
            if not cand_:
                continue
            sims_ = desc_[cand_] @ harita_d[lids_[i_]]
            b_ = int(np.argmax(sims_))
            if float(sims_[b_]) < thr:
                continue
            # ORB3 TrackLocalMap: ORBmatcher(0.8) en iyi / ikinci en iyi mesafe orani (ORBmatcher.cc SearchByProjection)
            if nn_oran > 0 and len(cand_) > 1:
                d2_ = 1.0 - float(np.partition(sims_, -2)[-2])
                if (1.0 - float(sims_[b_])) > nn_oran * d2_:
                    AUDIT['yh_oran_red'] += 1
                    continue
            k_ = int(cand_[b_])
            if k_ not in en_iyi or sims_[b_] > en_iyi[k_][0]:
                en_iyi[k_] = (float(sims_[b_]), lids_[i_])
        return {lid_: (float(kps_[k_, 0]), float(kps_[k_, 1])) for k_, (_, lid_) in en_iyi.items()}

    def _ts_yeni():
        if use_fl and str(getattr(args, 'smoother_tipi', 'batch')) == 'kayan':
            from kayan_pencere import ZamanHaritasi
            return ZamanHaritasi()
        return gu.FixedLagSmootherKeyTimestampMap()

    def _tekil_sabitle(hata, values_):
        """Belirsiz (tekil) dogrusal sistem: hata mesajindaki degiskeni o anki kestiriminde sabitleyip guncellemeyi
        tekrarla (en cok 5 kez). Pencerenin geri kalani korunur; rakiplerde aykiri harita noktasinin silinmesine denk.
        Basariliysa True."""
        import re as _re
        for _ in range(5):
            m_ = _re.search(r"Symbol: ([lx])(\d+)", str(hata))
            if not m_:
                return False
            ch_, i_ = m_.group(1), int(m_.group(2))
            g_ = NonlinearFactorGraph()
            if ch_ == 'l':
                try:
                    p_ = smoother.calculateEstimatePoint3(L(i_))
                    p_ = np.array([p_[0], p_[1], p_[2]])
                except Exception:
                    p_ = (lm_world.get(i_) if i_ in lm_world else
                          (np.array(values_.atPoint3(L(i_))) if values_.exists(L(i_)) else harita_p.get(i_)))
                if p_ is None:
                    return False
                g_.add(gtsam.PriorFactorPoint3(L(i_), Point3(*p_), yh_prior))
            else:
                try:
                    T_ = smoother.calculateEstimatePose3(X(i_)).matrix()
                except Exception:
                    T_ = pose_estimates.get(i_)
                if T_ is None:
                    return False
                g_.add(gtsam.PriorFactorPose3(X(i_), Pose3(T_), pose0_prior))
            AUDIT['tekil_sabitleme'] += 1
            try:
                smoother.update(g_, Values(), _ts_yeni())
                return True
            except Exception as e_:
                hata = e_
        return False

    def _yh_tanimlayici_konumla(f_, kare_):
        """ORB3 TrackReferenceKeyFrame karsiligi: poz ongorusu olmadan, karenin butun anahtar noktalari son pencerede
        gorulen harita noktalarina yalniz tanimlayiciyla eslenir (en iyi / ikinci en iyi orani 0.7, ORB3
        Tracking.cc TrackReferenceKeyFrame ORBmatcher(0.7); karsilikli en iyi), PnP-RANSAC ile poz cozulur.
        Donus: (T_WC, {lid: (u, v)}) ya da None."""
        kps_, desc_, _, _ = kare_
        pen_ = int(getattr(args, 'fixed_lag', 80) or 80)
        lids_ = [l for l, s_ in harita_son.items() if s_ >= f_ - pen_ and l in harita_d and l in harita_p
                 and harita_n.get(l, 0) >= int(args.yh_min_obs)]
        if len(lids_) < 10 or len(kps_) < 10:
            return None
        H_ = np.stack([harita_d[l] for l in lids_])
        S_ = desc_ @ H_.T                                    # (K, M) kosinüs benzerligi
        k_en = np.argmax(S_, axis=0)                         # her harita noktasi icin en iyi anahtar nokta
        Ss_ = np.sort(S_, axis=0)
        d1_, d2_ = 1.0 - Ss_[-1], 1.0 - Ss_[-2]
        m_en = np.argmax(S_, axis=1)                         # karsilikli kontrol
        sec_ = [(j_, int(k_en[j_])) for j_ in range(len(lids_))
                if d1_[j_] <= float(args.yh_desc_dist) and d1_[j_] < 0.7 * d2_[j_] and m_en[k_en[j_]] == j_]
        AUDIT['yh_konumla_aday'] += len(sec_)
        if len(sec_) < 10:
            return None
        obj_ = np.array([harita_p[lids_[j_]] for j_, _ in sec_], dtype=np.float64)
        img_ = np.array([kps_[k_] for _, k_ in sec_], dtype=np.float64)
        ok_, rv_, tv_, inl_ = cv2.solvePnPRansac(obj_.reshape(-1, 1, 3), img_.reshape(-1, 1, 2), K_cv, None,
                                                 reprojectionError=args.pnp_reproj,
                                                 iterationsCount=int(getattr(args, 'pnp_iters', 200)),
                                                 confidence=float(getattr(args, 'pnp_confidence', 0.99)),
                                                 flags=cv2.SOLVEPNP_ITERATIVE)
        if not ok_ or inl_ is None or len(inl_) < 10:     # ORB3 Relocalization: PnP en az 10 inlier
            return None
        R_, _ = cv2.Rodrigues(rv_)
        T_CW_ = np.eye(4); T_CW_[:3, :3] = R_; T_CW_[:3, 3] = tv_.ravel()
        AUDIT['yh_konumla_kabul'] += 1
        es_ = {lids_[sec_[i_][0]]: (float(img_[i_][0]), float(img_[i_][1])) for i_ in np.asarray(inl_).ravel()}
        return np.linalg.inv(T_CW_), es_

    _isam2_fl = use_fl and str(getattr(args, 'smoother_tipi', 'batch')) == 'isam2'
    n_frames = len(frames)
    print(f"{n_frames} kare işleniyor (BA, VO frame)...")
    _son_sure = [time.time()]
    _zaman = {'update': 0.0, 'estimate': 0.0}
    be_aralik = int(getattr(args, 'be_aralik', 1))
    be_paralel = int(getattr(args, 'be_paralel', 0)) > 0 and use_fl
    # 2: hizli hareketle tetiklenen guncelleme esz zamanli (gecikme hizli harekette izi bozuyor), sakin harekette arka plan
    be_paralel_hizli_esz = int(getattr(args, 'be_paralel', 0)) >= 2
    _par = {'t': None, 'hata': None, 'sure': 0.0}
    _ts_ayna = {}       # arka plan kipinde duzlestiricideki anahtarlarin zaman damgalari (is bittiginde)

    def _par_is(g_, v_, ts_):
        """Arka plan pencere guncellemesi (GIL'i birakan derlemede izlemeyle eszamanli calisir)."""
        t0_ = time.time()
        try:
            try:
                smoother.update(g_, v_, ts_)
            except Exception as e_ilk:
                if not (int(getattr(args, 'tekil_sabitle', 0)) and _tekil_sabitle(e_ilk, v_)):
                    raise
        except Exception as e_:
            _par['hata'] = e_
        _par['sure'] += time.time() - t0_

    def _par_katil():
        """Suren arka plan guncellemesini bekle (bekleme suresi 'update'e yazilir); hata verdiyse burada yukselt.
        Bekleyen is vardiysa True."""
        t_ = _par['t']
        if t_ is None:
            return False
        t0_ = time.time()
        t_.join()
        _zaman['update'] += time.time() - t0_
        _par['t'] = None
        e_, _par['hata'] = _par['hata'], None
        if e_ is not None:
            raise e_
        return True
    seyrek_k = int(getattr(args, 'gozlem_seyrek', 1))
    seyrek_min = int(getattr(args, 'gozlem_seyrek_min', 60))
    be_ac_deg = float(getattr(args, 'be_aralik_deg', 0.0))
    be_ac_m = float(getattr(args, 'be_aralik_m', 0.0))
    _son_upd_T = None
    _bek_g, _bek_v, _bek_ts = NonlinearFactorGraph(), Values(), {}   # aralikli guncellemede bekleyen faktorler

    # [TÜRBÜLANS] ölçüm modeli kolu (turb/backend.py). turb_arm=A0 -> None ve aşağıdaki bütün türbülans dalları
    # atlanır: PyEco kod yolu birebir. Kollar yalnız pencere faktörlerinin gürültü modelini değiştirir.
    turb_builder = None
    if str(getattr(args, 'turb_arm', 'A0') or 'A0').upper() != 'A0':
        from turb.backend import build_from_args as _turb_build
        turb_builder = _turb_build(args, fx, fy, cx, cy, baseline, K_gtsam, K_stereo)
        print(f"TÜRBÜLANS KOLU: {turb_builder.c.arm} (korelasyon={turb_builder.c.corr_model}, "
              f"sigma_oracle={turb_builder.c.sigma_px_oracle:.3f}px, robust: Huber mono "
              f"{turb_builder.c.huber_mono_k:g} / stereo {turb_builder.c.huber_stereo_k:g}, gözlem başına)")
    _tb_pf = turb_builder.model_for if (turb_builder is not None and turb_builder.c.arm == 'A2') else None
    # [TÜRBÜLANS/E5] pencere-göreli NEES: graf karelerinin faktörlerinin yuvarlanan geçmişi (yalnız --cov_out)
    cov_out = str(getattr(args, 'cov_out', '') or '')
    cov_every = max(1, int(getattr(args, 'cov_every', 1) or 1))
    cov_records = []
    _cov_n = 0
    _COND_LOG.clear()
    from collections import deque as _deque
    cov_hist = _deque(maxlen=(int(win) if use_fl else 100) + 4)

    for f, fr in enumerate(frames):
        ts = int(fr['timestamp'])
        frame_ts.append(ts)
        obs = fr['obs']  # (N,8): id,u,v,has_depth,fresh_depth,X,Y,Z
        f_meas, f_stereo = _kare_gurultu(fr)   # aday 3; kapaliyken meas_noise, stereo_noise
        AUDIT['frames'] += 1
        AUDIT['obs_total'] += len(obs)
        if obs.shape[0] and bool(np.any(obs[:, 4] > 0.5)):
            AUDIT['frames_fresh_stereo'] += 1

        graph = NonlinearFactorGraph()
        values = Values()
        chi2_skip = set()   # bu karede χ² ile elenen gözlemler
        basarisiz = False   # izleme basarisiz kare (motion-only BA inlier < izleme_basari_min)
        dondur = False      # basarisiz karede poz dondurulup factor'ler atlandi mi (kip 1)
        new_lids = []       # bu karede ilk kez eklenen landmark'lar (rollback için)
        yh_es = {}          # [Asama B] bu karede haritadan eslenen noktalar: lid -> (u, v)

        # 1) poz başlangıç tahmini
        if f == 0:
            T_WC = np.eye(4)
            graph.add(gtsam.PriorFactorPose3(X(0), Pose3(T_WC), pose0_prior))
        else:
            # PnP: lm_world cache'indeki landmark'larla (kare başına calculateEstimate yok)
            obj_pts, img_pts, ids_used = [], [], []
            inliers = None
            for row in obs:
                lid = int(row[0])
                if lid in lm_world:
                    p = lm_world[lid]
                    obj_pts.append([p[0], p[1], p[2]])
                    img_pts.append([row[1], row[2]])
                    ids_used.append(lid)
            # sabit-hız öngörüsü; PnP başarısızsa fallback
            T_pred = (prev_pose @ prev_rel).copy()
            if os.environ.get('PYECO_DBG_KAYDET') and f == int(os.environ['PYECO_DBG_KAYDET']):
                # teshis: bu karedeki PnP girdisini kaydet (yalniz ortam degiskeniyle; sonucu degistirmez)
                pickle.dump({'lm_world': dict(lm_world), 'obs': obs.copy(), 'T_pred': T_pred, 'prev_pose': prev_pose,
                             'K': (fx, fy, cx, cy), 'pose_estimates': {k: v.copy() for k, v in pose_estimates.items()}},
                            open(os.environ['PYECO_DBG_KAYDET'] + '_durum.pkl', 'wb'))
            T_WC = T_pred
            pnp_strong = False
            resync = False
            pnp_used_inliers = None      # kabul edilen PnP'nin inlier/outlier landmark kumeleri
            pnp_used_outliers = None
            if len(obj_pts) >= 6:
                obj_pts = np.asarray(obj_pts, dtype=np.float64)
                img_pts = np.asarray(img_pts, dtype=np.float64)
                ok, rvec, tvec, inliers = cv2.solvePnPRansac(
                    obj_pts.reshape(-1, 1, 3), img_pts.reshape(-1, 1, 2), K_cv, None,
                    reprojectionError=args.pnp_reproj, iterationsCount=int(getattr(args, 'pnp_iters', 200)),
                    confidence=float(getattr(args, 'pnp_confidence', 0.99)),
                    flags=cv2.SOLVEPNP_ITERATIVE)
                # inlier kapısı; pnp_min_inlier_frac>0 ise adaylarla oransal (6 ile min_inl arası)
                _need = min_inl
                if pnp_min_frac > 0 and min_inl > 0:
                    _need = max(6, min(min_inl, int(math.ceil(pnp_min_frac * len(obj_pts)))))
                AUDIT['pnp_called'] += 1
                AUDIT['pnp_cand'] += len(obj_pts)
                if inliers is not None:
                    AUDIT['pnp_inliers'] += len(inliers)
                if _need != min_inl:
                    AUDIT['pnp_frac_gate_active'] += 1
                # RANSAC dışı kalan landmark'lar strike alır
                if culler.pnp_strikes > 0 and ok and inliers is not None:
                    _inl = set(int(i) for i in np.asarray(inliers).ravel())
                    culler.on_pnp_outliers([ids_used[i] for i in range(len(ids_used)) if i not in _inl])
                if ok and rvec is not None and (_need <= 0 or
                                                (inliers is not None and len(inliers) >= _need)):
                    AUDIT['pnp_accepted'] += 1
                    R_cw, _ = cv2.Rodrigues(rvec)
                    T_CW = np.eye(4); T_CW[:3, :3] = R_cw; T_CW[:3, 3] = tvec.ravel()
                    T_WC = np.linalg.inv(T_CW)
                    # AirSLAM map_builder.cc:319-321: PnP sonucu son kabul edilen pozdan > jump_pnp_max
                    # metre uzaksa reddedilir; baslangic sabit-hiz ongorusune (ORB-SLAM3 TrackWithMotionModel) doner.
                    if inliers is not None:
                        _inl = set(int(i) for i in np.asarray(inliers).ravel())
                        pnp_used_inliers = {ids_used[i] for i in _inl}
                        pnp_used_outliers = {ids_used[i] for i in range(len(ids_used)) if i not in _inl}
                    if jump_pnp_max > 0 and float(np.linalg.norm(T_WC[:3, 3] - prev_pose[:3, 3])) > jump_pnp_max:
                        AUDIT['pnp_jump_rejected'] += 1
                        T_WC = T_pred
                        pnp_used_inliers = pnp_used_outliers = None
                    pnp_strong = (inliers is not None and
                                  len(inliers) >= max(min_inl, 6, int(PNP_STRONG_RATIO * len(obj_pts))))
                    if pnp_strong:
                        AUDIT['pnp_strong'] += 1
                else:
                    AUDIT['pnp_rejected'] += 1
            else:
                AUDIT['pnp_skipped_lt6'] += 1
            # --debug_frames N: ilk N karede poz kestirim zincirini bas
            if debug_frames and f < debug_frames:
                _ni = 0 if inliers is None else len(inliers)
                _acc = 'PnP-KABUL' if not np.allclose(T_WC, T_pred) else 'ONGORU(dead-reckon)'
                _why = ''
                if len(obj_pts) < 6:
                    _why = f'  <- lm_world yetersiz ({len(obj_pts)}<6)'
                elif _need > 0 and _ni < _need:
                    _why = f'  <- inlier {_ni} < gereken {_need}'
                print(f"  [dbg f={f:>3}] gozlem={len(obs):>3} lm_world={len(lm_world):>4} "
                      f"PnP-aday={len(obj_pts):>3} inlier={_ni:>3} {_acc}{_why}")

            # motion-only BA ile pozu rafine et
            if args.motion_only_ba:
                init_skip = None
                if int(getattr(args, 'moba_ilk_kapi', 0)):
                    if pnp_used_inliers is not None and int(args.moba_ilk_kapi) == 1:
                        init_skip = pnp_used_outliers
                        AUDIT['moba_kapi_pnp'] += 1
                    else:
                        # ORB-SLAM3 TrackWithMotionModel (Tracking.cc:2879-2894): stereo icin th=7 px pencere,
                        # 20'den az eslesmede 2*th. Ongoru pozunda bu pencere disina dusen gozlem disarida.
                        init_skip = _pencere_disi(T_WC, obs, lm_world, fx, fy, cx, cy, args.moba_pencere_px)
                        if sum(1 for r in obs if int(r[0]) in lm_world) - len(init_skip) < 20:
                            init_skip = _pencere_disi(T_WC, obs, lm_world, fx, fy, cx, cy, 2 * args.moba_pencere_px)
                        AUDIT['moba_kapi_ongoru'] += 1
                    AUDIT['moba_kapi_disari'] += len(init_skip)
                T_WC, chi2_skip = refine_pose_motion_only(
                    T_WC, obs, lm_world, K_gtsam, K_stereo, fx,
                    baseline, meas_noise, stereo_noise, args.max_init_depth,
                    args=args, fy=fy, cx=cx, cy=cy, return_outliers=True, init_skip=init_skip)
                culler.on_chi2_outliers(chi2_skip)
                if int(getattr(args, 'izleme_basari', 0)) and f > 0:
                    # Izleme basari olcutu: poz en az izleme_basari_min inlier gozlemle desteklenmeli
                    # (6 = ORB3 PnP cozucusunun en kucuk kumesi, Tracking.cc:3657 minSet; 10'luk mutlak esikler
                    # 1200/400 ozniteligin butcesine gore). Basarisiz karede
                    # poz son poza esitlenir (AIR map_builder.cc:321), hiz modeli gecersiz sayilir (ORB3
                    # Tracking.cc:2205-2215 hiz yalniz basarili izlemede guncellenir) ve var olan landmark'lara
                    # factor eklenmez; yalniz taze stereo yeni landmark kurabilir (yeniden ilklendirme).
                    _cand = {int(r[0]) for r in obs if int(r[0]) in lm_world}
                    _n_inl = len(_cand - chi2_skip)
                    if _n_inl < max(1, int(getattr(args, 'izleme_basari_min', 0) or min_inl)):
                        basarisiz = True
                        AUDIT['izleme_basarisiz'] += 1
                        if int(args.izleme_basari) == 1:
                            # kip 1: poz son poza esitlenir, var olan landmark'lara factor yok
                            T_WC = prev_pose.copy()
                            chi2_skip = set(_cand)
                            dondur = True
                        # kip 2: yalniz hiz modeli gecersiz (ORB3 Tracking.cc:2205-2215), gozlemler korunur
            # [Asama B] yerel harita takibi (ORB3 TrackLocalMap): takip pozundan sonra harita noktalarini esle, pozu
            # yeniden rafine et. Takip basarisizsa ongoruden genis yaricapla yeniden yakalama denenir.
            yh_es = {}
            if yh_on and args.motion_only_ba and f > 0:
                yh_kare = _yh_kare(fr, obs)
                if yh_kare is not None:
                    _kayip = basarisiz
                    _T0 = T_pred if _kayip else T_WC
                    _r = float(args.yh_yaricap_kayip_px) if _kayip else float(args.yh_yaricap_px)
                    _konum = None
                    if _kayip and int(getattr(args, 'yh_konumla', 0)):
                        _konum = _yh_tanimlayici_konumla(f, yh_kare)
                    if _konum is not None:
                        # tanimlayici konumundan izdusumle zenginlestir (ORB3: TrackReferenceKeyFrame sonrasi TrackLocalMap)
                        _T0 = _konum[0]
                        yh_es = _yh_esle(_T0, yh_kare, {int(r[0]) for r in obs}, float(args.yh_yaricap_px))
                        yh_es.update(_konum[1])
                    else:
                        yh_es = _yh_esle(_T0, yh_kare, {int(r[0]) for r in obs}, _r)
                    AUDIT['yh_aday_esles'] += len(yh_es)
                    if yh_es:
                        _rows = np.array([[l, uv[0], uv[1], 0.0, 0.0, np.nan, np.nan, np.nan]
                                          for l, uv in yh_es.items()], dtype=np.float64)
                        obs_aug = np.vstack([obs, _rows]) if len(obs) else _rows
                        lmw_aug = dict(lm_world)
                        lmw_aug.update({l: harita_p[l] for l in yh_es})
                        _skip0 = set() if _kayip else set(chi2_skip)
                        T_yh, _gate = refine_pose_motion_only(
                            _T0, obs_aug, lmw_aug, K_gtsam, K_stereo, fx, baseline, meas_noise, stereo_noise,
                            args.max_init_depth, args=args, fy=fy, cx=cx, cy=cy, return_outliers=True,
                            init_skip=_skip0)
                        # ORB3 PoseOptimization sonrasi chi2 > 5.991 eslesme aykiri: haritadan gelen eslesmeler suzulur
                        _TCW = np.linalg.inv(T_yh)
                        _th = float(args.chi2_thresh_mono) * float(args.px_noise) ** 2
                        for l in list(yh_es):
                            _pc = _TCW[:3, :3] @ harita_p[l] + _TCW[:3, 3]
                            if _pc[2] <= 1e-3:
                                yh_es.pop(l); continue
                            _du = fx * _pc[0] / _pc[2] + cx - yh_es[l][0]
                            _dv = fy * _pc[1] / _pc[2] + cy - yh_es[l][1]
                            if _du * _du + _dv * _dv > _th:
                                yh_es.pop(l)
                        AUDIT['yh_esles'] += len(yh_es)
                        _cand = {int(r[0]) for r in obs if int(r[0]) in lm_world}
                        _n_inl = len(_cand - _gate) + len(yh_es)
                        if _n_inl >= max(6, int(getattr(args, 'izleme_basari_min', 0) or 6)):
                            T_WC = T_yh
                            chi2_skip = set(_gate)
                            if _kayip:
                                basarisiz = False
                                dondur = False
                                AUDIT['yh_yeniden_yakalama'] += 1
                        else:
                            yh_es = {}
            if os.environ.get('PYECO_IZ'):
                _a, _b = (int(x) for x in os.environ['PYECO_IZ'].split(':'))
                if _a <= f < _b:
                    _ni = 0 if inliers is None else len(inliers)
                    _step = float(np.linalg.norm(T_WC[:3, 3] - prev_pose[:3, 3]))
                    _pstep = float(np.linalg.norm(T_pred[:3, 3] - prev_pose[:3, 3]))
                    print(f"IZ f={f} obs={len(obs)} lmw_aday={sum(1 for r in obs if int(r[0]) in lm_world)} "
                          f"pnp_inl={_ni} pnp_guclu={int(pnp_strong)} chi2_disi={len(chi2_skip)} yh={len(yh_es)} "
                          f"basarisiz={int(basarisiz)} adim={_step:.3f} ongoru_adim={_pstep:.3f} taze={int(np.any(obs[:, 4] > 0.5)) if len(obs) else 0}")
            # sıçrama kapısı: tekil spike reddedilir, ısrarcı güçlü ölçüm kazanır (re-sync)
            if jump_gate > 0:
                # referans 'prediction' (öngörü) ya da 'last_pose' (son kabul edilen poz; kilitlenmez)
                AUDIT['jump_gate_evals'] += 1
                _ref = T_pred if jump_gate_ref == 'prediction' else prev_pose
                if float(np.linalg.norm(T_WC[:3, 3] - _ref[:3, 3])) > jump_gate:
                    AUDIT['jump_gate_trips'] += 1
                    if (pnp_strong and jump_rejects >= 1) or jump_rejects >= JUMP_MAX_CONSEC:
                        resync = True
                        jump_rejects = 0
                        AUDIT['jump_resync'] += 1
                    else:
                        T_WC = T_pred
                        jump_rejects += 1
                        AUDIT['jump_rejected_pose'] += 1
                else:
                    jump_rejects = 0
            # keyframe kararı; taze stereo taşıyan kare daima keyframe
            if kf_mode:
                rel_kf = np.linalg.inv(last_kf_pose) @ T_WC
                d_trans = float(np.linalg.norm(rel_kf[:3, 3]))
                d_rot = float(np.linalg.norm(Rot.from_matrix(rel_kf[:3, :3]).as_rotvec()))
                n_tracked = sum(1 for r in obs if int(r[0]) in lm_world)
                if int(getattr(args, 'kf_harita_say', 0)):
                    n_tracked += len(yh_es)      # ORB3 mnMatchesInliers: yerel haritadan izlenenler de sayilir
                has_fresh = bool(np.any(obs[:, 4] > 0.5)) if obs.shape[0] else False
                gap = f - prev_graph_f
                tracked_ratio = n_tracked / max(ref_tracked, 1)
                # taze stereo karesi anahtar kare olur; kf_taze_aralik > 1 ise en erken bu kadar karede bir
                _taze_kf = has_fresh and gap >= max(1, int(getattr(args, 'kf_taze_aralik', 1)))
                is_kf = (_taze_kf or gap >= args.kf_max_interval or
                         (gap >= args.kf_min_interval and
                          (d_trans >= args.kf_min_trans or d_rot >= args.kf_min_rot or
                           tracked_ratio < args.kf_max_tracked_ratio)))
                if not is_kf:
                    pose_estimates[f] = T_WC                  # ara kare; çıktıda KF'den yeniden kurulur
                    inter_ref[f] = prev_graph_f
                    inter_rel[f] = np.linalg.inv(prev_graph_pose) @ T_WC
                    kf_flags[f] = False
                    # ara karenin hizi: izleme basarisizsa hiz modeli gecersiz (ORB3 Tracking.cc:2205-2215)
                    prev_rel = np.eye(4) if basarisiz else np.linalg.inv(prev_pose) @ T_WC
                    prev_pose = T_WC
                    # gercek zamanli poz = motion-only BA sonucu (ORB3/AirSLAM izleme pozu); duzeltme oncesi kod bu karelerde
                    # nedensel kaydi atliyordu
                    if _NEDENSEL is not None:
                        _NEDENSEL[f] = T_WC.copy()
                    AUDIT['kf_ara_kare'] += 1
                    continue
                AUDIT['kf_anahtar_kare'] += 1
            values.insert(X(f), Pose3(T_WC))
            # sabit-hız factor'u; keyframe modunda son commit edilen KF'ye odometri factor'u
            if args.cv_motion and not (kf_mode and not int(getattr(args, 'kf_odometri', 1))):
                if kf_mode:
                    rel_graph = np.linalg.inv(prev_graph_pose) @ T_WC
                    graph.add(gtsam.BetweenFactorPose3(X(prev_graph_f), X(f), Pose3(rel_graph),
                                                       _motion_noise_for(rel_graph)))
                else:
                    # factor son başarılı commit edilen poza bağlanır; boşluk/re-sync varsa odometri ölçümü
                    if basarisiz:
                        rel_cv = np.eye(4)            # hiz bilinmiyor: sifir hareket
                    elif prev_graph_f == f - 1 and not resync:
                        rel_cv = prev_rel
                    else:
                        rel_cv = np.linalg.inv(prev_graph_pose) @ T_WC
                    # adaptive σ önceki hızdan (prev_rel) türetilir, ölçümün kendisinden değil
                    graph.add(gtsam.BetweenFactorPose3(X(prev_graph_f), X(f), Pose3(rel_cv),
                                                       _motion_noise_for(prev_rel)))
                    AUDIT['cv_factors'] += 1
                    if prev_graph_f != f - 1:
                        AUDIT['cv_factor_gap'] += 1

        if f == 0:
            values.insert(X(0), Pose3(T_WC))

        # 2) gözlemler: taze stereo'da stereo factor (metrik ölçek), aksi halde mono projeksiyon
        n_proj, n_stereo, n_new, n_gated = 0, 0, 0, 0
        # bootstrap kararı kare başında bir kez verilir
        bootstrapping = bootstrap_depth and len(lm_world) < max(6, min_inl)
        if bootstrapping:
            AUDIT['bootstrap_frames'] += 1
        if len(lm_world) < max(6, min_inl):
            AUDIT['frames_map_starved'] += 1
        frame_lids = []  # bu karede factor alan landmark'lar (fixed-lag timestamp)
        if turb_builder is not None:
            # graf karesinin kol modelleri (f'ten önceki graf karelerinin artıklarından; nedensel)
            f_meas, f_stereo = turb_builder.begin_frame(X(f))
            _tb_TCW = np.linalg.inv(T_WC)

            def _turb_pred(p_w):
                """Kare başındaki harita konumunun bu karenin izleme pozuna izdüşümü (artık kestiricileri için)."""
                if p_w is None or not turb_builder._needs_residuals:
                    return None
                pc_ = _tb_TCW @ np.array([p_w[0], p_w[1], p_w[2], 1.0])
                if pc_[2] <= 1e-3:
                    return None
                return (fx * pc_[0] / pc_[2] + cx, fy * pc_[1] / pc_[2] + cy)
        # ORB-SLAM3 uzak nokta kurali: ilklendirmede (Tracking.cc:2335 StereoInitialization) z>0 olan TUM
        # stereo noktalar; sonraki stereo karelerde (CreateNewKeyFrame, Tracking.cc:3251-3327) derinlige gore
        # sirali, yakinlarin hepsi + toplam stereo nokta sayisi maxPoint'i (100, butceyle olcekli) asana dek uzaklar.
        uzak_ok = set()
        if int(getattr(args, 'stereo_uzak_min', 0)) > 0 or (int(getattr(args, 'stereo_ilk_hepsi', 0)) and f == 0):
            _fr = sorted((float(r[7]), int(r[0])) for r in obs
                         if r[4] > 0.5 and r[3] > 0.5 and np.isfinite(r[7]) and r[7] > 0)
            if int(getattr(args, 'stereo_ilk_hepsi', 0)) and f == 0:
                uzak_ok = {l for _, l in _fr}
            else:
                _n = 0
                for _z, _l in _fr:
                    _n += 1
                    if _z >= args.max_init_depth:
                        uzak_ok.add(_l)
                    if _z > args.max_init_depth and _n > int(args.stereo_uzak_min):
                        break
        T_CW_cur = np.linalg.inv(T_WC)
        active_set = live_lms if use_fl else initialized_lms
        # kare başına mv-init sayısı sınırlı; çok sayıda yeni değişken smoother'ı indeterminant yapıyor
        n_mv_this_frame = 0
        for row in obs:
            lid = int(row[0])
            u, v = float(row[1]), float(row[2])
            uv = np.array([u, v], dtype=np.float64)
            if turb_builder is not None:
                _tb_pw = lm_world.get(lid)      # bu karede ilklenen nokta artığa girmez (artık ~0 olurdu)
            has_depth = row[3] > 0.5
            fresh = row[4] > 0.5
            Z = row[7]
            stereo_ok = fresh and has_depth and np.isfinite(Z) and 0 < Z < args.max_init_depth
            if (not stereo_ok) and lid in uzak_ok:
                stereo_ok = True
                AUDIT['stereo_uzak_kabul'] += 1
            # bootstrap: harita boşken 'fresh' şartı gevşetilir, bayat derinlik de kullanılır
            if (not stereo_ok) and bootstrapping:
                stereo_ok = has_depth and np.isfinite(Z) and 0 < Z < args.max_init_depth
                if stereo_ok:
                    AUDIT['bootstrap_obs_rescued'] += 1

            # χ² outlier'ı grafa sokma (init dalından önce olmalı)
            if chi2_gate and lid in chi2_skip:
                n_gated += 1
                AUDIT['gated_obs_chi2'] += 1
                continue

            # kara-listedeki landmark'a factor eklenmez; yalnız allow_reinit + taze stereo + grafta yoksa yeniden yaratılır
            if culler.active and culler.is_bad(lid):
                _resident = (lid in resident_lms) if use_fl else (lid in initialized_lms)
                if not (culler.allow_reinit and stereo_ok and not _resident):
                    AUDIT['cull_skipped_obs'] += 1
                    continue

            active_set = live_lms if use_fl else initialized_lms
            if lid not in active_set:
                # yalnız taze stereo derinliğiyle init; smoother'da hâlâ duran değişken yeniden insert edilemez
                if stereo_ok and np.all(np.isfinite(row[5:8])) and not (use_fl and lid in resident_lms):
                    p_world = (T_WC @ np.array([row[5], row[6], row[7], 1.0]))[:3]
                    values.insert(L(lid), Point3(*p_world))
                    culler.on_init(lid, f)
                    # fixed-lag'de prior eklenmez; marjinalizasyonda smoother'ı bozuyor
                    if not use_fl:
                        graph.add(gtsam.PriorFactorPoint3(L(lid), Point3(*p_world), lm_prior))
                    elif zayif_onsel is not None:
                        graph.add(gtsam.PriorFactorPoint3(L(lid), Point3(*p_world), zayif_onsel))
                    active_set.add(lid)
                    if not use_fl:
                        initialized_lms.add(lid)
                    lm_world[lid] = p_world
                    n_new += 1
                    new_lids.append(lid)
                    if lid in _lm_ever:
                        AUDIT['lm_reinit'] += 1     # marjinalize olup yeniden init edilen
                    else:
                        AUDIT['lm_first_init'] += 1
                        _lm_ever.add(lid)
                    pending_mv.pop(lid, None)  # stereo geldi, mv tamponu gereksiz
                elif int(getattr(args, 'mv_triangulation', 0)) and not dondur:
                    # mono-only track: poz+piksel biriktir, yeterli paralaksta DLT ile init et
                    buf = pending_mv.setdefault(lid, [])
                    buf.append((f, T_WC.copy(), u, v, T_CW_cur.copy()))   # ters bir kez hesaplanır
                    AUDIT['mv_buffered'] += 1
                    if len(buf) > args.mv_max_track_age:
                        buf.pop(0)
                    if len(buf) < args.mv_min_obs:
                        AUDIT['mv_rej_min_obs'] += 1
                        continue
                    if n_mv_this_frame >= args.mv_max_per_frame:
                        AUDIT['mv_rej_per_frame_cap'] += 1
                        continue
                    AUDIT['mv_tri_attempts'] += 1
                    poses_buf = [b[1] for b in buf]
                    tcw_buf = [b[4] for b in buf]
                    p_cand, mv_cond = _triangulate_multiview(poses_buf, [(b[2], b[3]) for b in buf],
                                                             fx, fy, cx, cy, T_CW=tcw_buf)
                    if p_cand is None:
                        AUDIT['mv_rej_degenerate'] += 1
                        continue
                    if args.mv_min_condition > 0 and mv_cond < args.mv_min_condition:
                        AUDIT['mv_rej_condition'] += 1
                        continue
                    if _mv_parallax_deg(poses_buf, p_cand) < args.mv_min_parallax_deg:
                        AUDIT['mv_rej_parallax'] += 1
                        continue
                    ok = True
                    _pc_h = np.array([p_cand[0], p_cand[1], p_cand[2], 1.0])
                    for (_f_b, T_WC_b, u_b, v_b, _TCW_b) in buf:
                        pc = _TCW_b @ _pc_h
                        if pc[2] <= 1e-3 or pc[2] > args.max_init_depth:
                            ok = False; break
                        du_b = (fx * pc[0] / pc[2] + cx) - u_b
                        dv_b = (fy * pc[1] / pc[2] + cy) - v_b
                        if du_b * du_b + dv_b * dv_b > args.mv_max_reproj_err ** 2:
                            ok = False; break
                    if not ok:
                        AUDIT['mv_rej_reproj_or_depth'] += 1
                        continue
                    values.insert(L(lid), Point3(*p_cand))
                    # 'factors' kipi: biriken gözlemleri pencerede duran pozlara mono factor olarak ekle;
                    # ≥2 canlı poz yoksa 'prior' kipine düşülür (DLT çözümü prior olarak).
                    mv_added = 0
                    if mv_init_mode == 'factors':
                        for (_fb, _T, _ub, _vb, _TCWb) in buf:
                            if _fb == f or _fb not in live_poses:
                                continue      # bu karenin factor'ü aşağıda ekleniyor
                            if _isam2_fl and _fb <= kf_seq - win + 2:
                                continue      # ISAM2 sabit gecikme: bu guncellemede marjinalize olacak poza faktor eklenmez
                            graph.add(gtsam.GenericProjectionFactorCal3_S2(
                                np.array([_ub, _vb]),
                                meas_noise if turb_builder is None else
                                (f_meas if _tb_pf is None else _tb_pf(lid, False)),
                                X(_fb), L(lid), K_gtsam, False, False))
                            mv_added += 1
                        AUDIT['mv_factors_added'] += mv_added
                    if mv_added >= 2:
                        AUDIT['mv_init_by_factors'] += 1
                        if zayif_onsel is not None:
                            graph.add(gtsam.PriorFactorPoint3(L(lid), Point3(*p_cand), zayif_onsel))
                    else:
                        graph.add(gtsam.PriorFactorPoint3(L(lid), Point3(*p_cand), mv_prior))
                        AUDIT['mv_init_by_prior'] += 1
                    active_set.add(lid)
                    if not use_fl:
                        initialized_lms.add(lid)
                    lm_world[lid] = p_cand
                    n_new += 1
                    new_lids.append(lid)
                    n_mv_init += 1
                    n_mv_this_frame += 1
                    if int(getattr(args, 'mv_debug', 0)):
                        par = _mv_parallax_deg(poses_buf, p_cand)
                        print(f"  [kare {f}] MV-init lid={lid} nobs={len(buf)} paralaks={par:.2f}° "
                              f"cond={mv_cond:.2f} p={p_cand}")
                    del pending_mv[lid]
                else:
                    AUDIT['obs_dropped_no_init'] += 1
                    continue

            # piksel outlier kapısı: geri-yansıtma hatası eşiği aşarsa factor ekleme
            if args.outlier_gate_px > 0 and lid in lm_world:
                AUDIT['outlier_gate_evals'] += 1
                pw = lm_world[lid]
                pc = T_CW_cur @ np.array([pw[0], pw[1], pw[2], 1.0])
                if pc[2] > 1e-3:
                    du = (fx * pc[0] / pc[2] + cx) - u
                    dv = (fy * pc[1] / pc[2] + cy) - v
                    if du * du + dv * dv > args.outlier_gate_px ** 2:
                        n_gated += 1
                        AUDIT['gated_obs_px'] += 1
                        continue

            # gozlem seyreltme: uzun izin her karedeki gozlemi Schur maliyetini (6n)^2 buyutur; mono izler k karede bir
            # (kimlige gore kaydirmali, her karede ~1/k'si) grafa girer. Taze stereo ve bu karede ilklenen nokta her zaman.
            if (seyrek_k > 1 and not stereo_ok and len(obs) > seyrek_min and (f + lid) % seyrek_k != 0
                    and not (new_lids and new_lids[-1] == lid)):
                AUDIT['seyrek_atlanan'] += 1
                continue

            if turb_builder is not None and not turb_builder.inline:
                # [TÜRBÜLANS] yama kolu: gözlem kare sonunda yama faktörüne girer
                turb_builder.add(lid, u, v, Z if stereo_ok else None, _turb_pred(_tb_pw),
                                 must_keep=bool(new_lids) and new_lids[-1] == lid)
                if stereo_ok:
                    n_stereo += 1
                    AUDIT['fac_stereo'] += 1
                else:
                    AUDIT['fac_mono'] += 1
            elif stereo_ok:
                uR = u - (fx * baseline / Z)   # uR derinlikten
                graph.add(gtsam.GenericStereoFactor3D(
                    gtsam.StereoPoint2(u, uR, v), f_stereo if _tb_pf is None else _tb_pf(lid, True),
                    X(f), L(lid), K_stereo))
                n_stereo += 1
                AUDIT['fac_stereo'] += 1
            else:
                graph.add(gtsam.GenericProjectionFactorCal3_S2(
                    uv, f_meas if _tb_pf is None else _tb_pf(lid, False), X(f), L(lid), K_gtsam, False, False))
                AUDIT['fac_mono'] += 1
            if turb_builder is not None and turb_builder.inline:
                turb_builder.observe(lid, u, v, _turb_pred(_tb_pw))
            frame_lids.append(lid)
            culler.on_observation(lid)
            n_proj += 1

        # 2a) [Asama B] haritadan eslenen noktalar: pencerede ise dogrudan, degilse harita konumunda onselle geri
        for lid, (u, v) in yh_es.items():
            if chi2_gate and lid in chi2_skip:
                continue
            if lid not in active_set:
                if not use_fl or lid in resident_lms:
                    continue
                p_h = harita_p[lid]
                values.insert(L(lid), Point3(*p_h))
                graph.add(gtsam.PriorFactorPoint3(L(lid), Point3(*p_h), yh_prior))
                active_set.add(lid)
                lm_world[lid] = p_h
                new_lids.append(lid)
                AUDIT['yh_geri_donen'] += 1
            if turb_builder is not None and not turb_builder.inline:
                turb_builder.add(lid, u, v, None, _turb_pred(harita_p.get(lid)),
                                 must_keep=bool(new_lids) and new_lids[-1] == lid)
            else:
                graph.add(gtsam.GenericProjectionFactorCal3_S2(
                    np.array([u, v]), f_meas if _tb_pf is None else _tb_pf(lid, False),
                    X(f), L(lid), K_gtsam, False, False))
                if turb_builder is not None:
                    turb_builder.observe(lid, u, v, _turb_pred(harita_p.get(lid)))
            frame_lids.append(lid)
            n_proj += 1
            AUDIT['yh_faktor'] += 1

        # 2b) search-by-projection: takipte kopmuş harita noktalarını descriptor ile geri kazan
        n_recovered = 0
        if (int(getattr(args, 'search_by_projection', 0)) and f > 0 and not dondur
                and fr.get('kps') is not None and fr.get('desc') is not None
                and len(fr['kps']) > 0):
            AUDIT['sbp_frames_active'] += 1
            kps = np.asarray(fr['kps'], dtype=np.float64).reshape(-1, 2)
            desc = np.asarray(fr['desc'], dtype=np.float64)
            desc = desc / np.maximum(np.linalg.norm(desc, axis=1, keepdims=True), 1e-12)
            obs_ids = set(int(r[0]) for r in obs)
            used = set()  # atanmış kp indeksleri
            # izlenen landmark'ların descriptor'ını tazele
            for r in obs:
                d2 = (kps[:, 0] - r[1]) ** 2 + (kps[:, 1] - r[2]) ** 2
                j = int(np.argmin(d2)) if len(d2) else -1
                if j >= 0 and d2[j] <= 2.0 ** 2:
                    lm_desc[int(r[0])] = desc[j]
                    used.add(j)
            # görülmeyen ama grafta olan landmark'ları projekte et ve eşleştir
            for lid, pw in list(lm_world.items()):
                if lid in obs_ids or lid in yh_es:
                    AUDIT['sbp_skip_seen'] += 1
                    continue
                if lid not in lm_desc:
                    AUDIT['sbp_skip_nodesc'] += 1
                    continue
                if lid not in active_set:
                    AUDIT['sbp_skip_inactive'] += 1
                    continue
                AUDIT['sbp_candidates'] += 1
                pc = T_CW_cur @ np.array([pw[0], pw[1], pw[2], 1.0])
                if pc[2] <= 1e-3:
                    AUDIT['sbp_skip_behind'] += 1
                    continue
                up = fx * pc[0] / pc[2] + cx
                vp = fy * pc[1] / pc[2] + cy
                if not (0 <= up <= 2 * cx and 0 <= vp <= 2 * cy):
                    AUDIT['sbp_skip_offimage'] += 1
                    continue
                d2 = (kps[:, 0] - up) ** 2 + (kps[:, 1] - vp) ** 2
                cand = np.where(d2 <= float(args.sbp_reproj_px) ** 2)[0]
                cand = np.array([c for c in cand if c not in used])
                if cand.size == 0:
                    AUDIT['sbp_skip_nocand'] += 1
                    continue
                sims = desc[cand] @ lm_desc[lid]
                best = int(np.argmax(sims))
                if (1.0 - float(sims[best])) > float(args.sbp_desc_dist):
                    AUDIT['sbp_skip_desc'] += 1
                    continue
                k = int(cand[best])
                used.add(k)
                if turb_builder is not None and not turb_builder.inline:
                    turb_builder.add(lid, float(kps[k, 0]), float(kps[k, 1]), None, _turb_pred(pw))
                else:
                    graph.add(gtsam.GenericProjectionFactorCal3_S2(
                        np.array([kps[k, 0], kps[k, 1]]), f_meas if _tb_pf is None else _tb_pf(lid, False),
                        X(f), L(lid), K_gtsam, False, False))
                    if turb_builder is not None:
                        turb_builder.observe(lid, float(kps[k, 0]), float(kps[k, 1]), _turb_pred(pw))
                frame_lids.append(lid)
                n_recovered += 1
                AUDIT['sbp_recovered'] += 1
            n_proj += n_recovered
            AUDIT['sbp_desc_cached'] = len(lm_desc)

        # [TÜRBÜLANS] kare sonu: yama kolunun faktörleri grafa yazılır, artıklar kestiricilere itilir
        if turb_builder is not None:
            turb_builder.flush(graph)
        if cov_out:
            cov_hist.append((f, graph))

        # 2c) landmark culling (update'ten önce; kara-listeye alınan sonraki kareden factor almaz)
        n_culled_frame = 0
        if culler.active:
            _obs_ids = set(int(r[0]) for r in obs)
            culler.visibility_pass(T_WC, lm_world, active_set, _obs_ids)
            _culled = culler.evaluate(f, active_set)
            if _culled:
                # tam-grafta değişken grafta kalır, yalnız defterlerden düşer
                culler.apply_cull(_culled, lm_world, live_lms, dicts=(lm_desc, pending_mv))
                n_culled_frame = len(_culled)
                n_culled_total += n_culled_frame
                AUDIT['cull_landmarks'] += n_culled_frame

        # 3) iSAM2 / fixed-lag update
        update_ok = True
        _atlandi = False    # aralikli guncelleme: bu karenin faktorleri bekletildi, poz izlemeden
        _paralel_bu = False  # bu karede guncelleme arka planda baslatildi (cikti pozu izlemeden)
        _hareket_tetik = False   # guncelleme hizli hareketle tetiklendi (arka plana atilmaz)
        _est_hazir = None
        try:
            if use_fl:
                # görülen değişkenlerin timestamp'i tazelenir; görülmeyenler marjinalize olur
                _bek_ts[X(f)] = float(kf_seq)
                for lid in frame_lids:
                    _bek_ts[L(lid)] = float(kf_seq)
                _bek_g.push_back(graph)
                _bek_v.insert(values)
                # rakiplerdeki izleme/haritalama ayrimi: her kare izlenir ve grafa girer, pencere optimizasyonu
                # be_aralik karede bir (ara karenin cikti pozu harita uzerinden motion-only BA'dan)
                _atlandi = be_aralik > 1 and f > 0 and f % be_aralik != 0 and f != n_frames - 1
                if _atlandi and _son_upd_T is not None and (be_ac_deg > 0 or be_ac_m > 0):
                    # uyarlamali aralik: son optimizasyondan beri birikmis hareket buyukse hemen optimize et
                    _rel = np.linalg.inv(_son_upd_T) @ T_WC
                    _ang = float(np.degrees(np.linalg.norm(Rot.from_matrix(_rel[:3, :3]).as_rotvec())))
                    if (be_ac_deg > 0 and _ang > be_ac_deg) or (be_ac_m > 0 and float(np.linalg.norm(_rel[:3, 3])) > be_ac_m):
                        _atlandi = False
                        _hareket_tetik = True
                        AUDIT['be_aralik_hareket'] += 1
                if not _atlandi:
                    _son_upd_T = T_WC.copy()
                if not _atlandi:
                    ts_map = _ts_yeni()
                    for _k, _t in _bek_ts.items():
                        ts_map.insert((_k, _t))
                    _g_u, _v_u, _ts_gonder = _bek_g, _bek_v, _bek_ts
                    _bek_g, _bek_v, _bek_ts = NonlinearFactorGraph(), Values(), {}
                    if be_paralel:
                        # zaman damgasi aynasi: is bitince duzlestiricide kalacak anahtarlar (BatchFixedLagSmoother:
                        # en yeni damga - pencereden eski olanlar marjinalize)
                        for _k, _t in _ts_gonder.items():
                            _ts_ayna[_k] = max(_t, _ts_ayna.get(_k, _t))
                        _kes = max(_ts_ayna.values()) - win
                        for _k in [_k for _k, _t in _ts_ayna.items() if _t < _kes]:
                            del _ts_ayna[_k]
                    if be_paralel and f != n_frames - 1 and not (be_paralel_hizli_esz and _hareket_tetik):
                        # yerel haritalama is parcacigi karsiligi: onceki guncelleme bu karede beklenir ve kestirimi
                        # alinir, yenisi arka planda baslar; birlesme noktasi sabit oldugu icin sonuc zamanlamadan
                        # bagimsizdir (belirlenimci)
                        _vardi = _par_katil()
                        _est_hazir = smoother.calculateEstimate() if _vardi else None
                        _par['t'] = threading.Thread(target=_par_is, args=(_g_u, _v_u, ts_map), daemon=True)
                        _par['t'].start()
                        _paralel_bu = True
                    else:
                        _par_katil()
                        try:
                            _t0u = time.time()
                            smoother.update(_g_u, _v_u, ts_map)
                            _zaman['son'] = time.time() - _t0u
                            _zaman['update'] += _zaman['son']
                        except Exception as e_ilk:
                            if not (int(getattr(args, 'tekil_sabitle', 0)) and _tekil_sabitle(e_ilk, _v_u)):
                                raise
            else:
                _t0u = time.time()
                isam.update(graph, values)
                for _ in range(args.extra_updates):
                    isam.update()
                _zaman['update'] += time.time() - _t0u
            consec_err = 0
        except Exception as e:
            update_ok = False
            AUDIT['update_err'] += 1
            AUDIT['update_err_' + type(e).__name__] += 1
            # rollback: commit edilemeyen landmark'ları defterlerden sil, pozu ham init ile sakla
            for lid in new_lids:
                culler.forget(lid)
                initialized_lms.discard(lid)
                live_lms.discard(lid)
                lm_world.pop(lid, None)
            pose_estimates[f] = T_WC
            consec_err += 1
            total_err += 1
            # smoother_recover: hata sonrası smoother son sağlıklı poz üzerinde sıfırdan kurulur
            if int(getattr(args, 'smoother_recover', 0)) and use_fl:
                try:
                    anchor_f = prev_graph_f
                    T_anchor = pose_estimates.get(anchor_f, T_WC)
                    if str(getattr(args, 'smoother_tipi', 'batch')) == 'isam2':
                        _ip = ISAM2Params()
                        _ip.setRelinearizeThreshold(float(getattr(args, 'relinearize_threshold', 0.01)))
                        _ip.relinearizeSkip = int(getattr(args, 'relinearize_skip', 1))
                        smoother = gu.IncrementalFixedLagSmoother(float(win), _ip)
                    elif _kayan:
                        smoother = _kayan_yeni()
                    else:
                        smoother = gu.BatchFixedLagSmoother(float(win), lm_params)
                    g_r = gtsam.NonlinearFactorGraph()
                    v_r = gtsam.Values()
                    v_r.insert(X(anchor_f), Pose3(T_anchor))
                    g_r.add(gtsam.PriorFactorPose3(X(anchor_f), Pose3(T_anchor),
                                                   pose0_prior))
                    ts_r = _ts_yeni()
                    ts_r.insert((X(anchor_f), float(kf_seq)))
                    smoother.update(g_r, v_r, ts_r)
                    live_lms = set()
                    _ts_ayna = {X(anchor_f): float(kf_seq)}
                    if int(getattr(args, 'kurtarma_duzelt', 0)):
                        # HATA DUZELTMESI: yeni smoother'da yalniz capa pozu var; duzeltme oncesi kod eski pencerenin poz/landmark
                        # defterlerini tutuyordu -> sonraki karelerde olmayan anahtarlara factor ("x2064 does not exist")
                        # ve zincirleme yeniden kurulum.
                        live_poses = {anchor_f}
                        resident_lms = set()
                    consec_err = 0
                    n_recover += 1
                    AUDIT['smoother_recover'] += 1
                    if n_recover <= err_print_cap:
                        print(f"  [kare {f}] smoother YENIDEN KURULDU "
                              f"(cipa=X({anchor_f}), toplam kurtarma={n_recover})")
                except Exception as e2:
                    print(f"  [kare {f}] smoother kurtarma BASARISIZ: "
                          f"{type(e2).__name__}: {e2}")
            if total_err <= err_print_cap:
                msg = " | ".join(l.strip() for l in str(e).strip().splitlines()[:4] if l.strip()) or type(e).__name__
                print(f"  [kare {f}] update hatası: {msg}")
            elif total_err == err_print_cap + 1:
                print(f"  [kare {f}] ... (sonraki update hataları susturuldu)")
            if consec_err >= DIVERGE_CONSEC:
                raise RuntimeError(
                    f"IRAKSAK config: {consec_err} ardışık update hatası (kare {f}). "
                    f"Erken çıkılıyor (toplam {total_err} hata).")
            if f + 1 >= MIN_FRAMES_FOR_RATIO and total_err > MAX_ERR_RATIO * (f + 1):
                raise RuntimeError(
                    f"IRAKSAK config: {total_err}/{f+1} karede update hatası "
                    f"(oran={total_err/(f+1):.1%} > {MAX_ERR_RATIO:.0%}, ardışık olmasa da yaygın "
                    f"kararsızlık). Erken çıkılıyor.")

        # bayat pending_mv tamponlarını temizle
        if pending_mv and f % 200 == 0:
            stale = [lid for lid, buf in pending_mv.items() if f - buf[-1][0] > args.mv_max_track_age]
            for lid in stale:
                del pending_mv[lid]

        # sabit-hız durumunu güncelle; sıçrama kapısını aşan hareket CV durumuna alınmaz
        if f > 0:
            rel_new = np.linalg.inv(prev_pose) @ T_WC
            if resync:
                # resync_velocity: 'zero' identity'ye sıfırlar, 'last_good' son sağlıklı hızı korur
                prev_rel = np.eye(4) if resync_velocity == 'zero' else last_good_rel
            else:
                if basarisiz:
                    prev_rel = np.eye(4)              # ORB3: basarisiz izlemede hiz modeli gecersiz
                elif not (jump_gate > 0 and float(np.linalg.norm(rel_new[:3, 3])) > jump_gate):
                    prev_rel = rel_new
                    if float(np.linalg.norm(rel_new[:3, 3])) > 1e-4:   # donuk kare (rel≈I) sayılmaz
                        last_good_rel = rel_new
        prev_pose = T_WC

        # commit defteri yalnız update başarılıysa güncellenir
        if update_ok:
            prev_graph_f = f
            prev_graph_pose = T_WC
            kf_seq += 1
            if kf_mode:
                kf_flags[f] = True
                keyframes.append(f)
                last_kf_pose = T_WC
                ref_tracked = len(frame_lids)

        # loop closure: commit edilen node'u appearance DB'sine gönder
        if update_ok and loop_closer is not None and (kf_mode or (f % loop_stride == 0)):
            kps_f = fr.get('kps'); desc_f = fr.get('desc')
            if kps_f is None or desc_f is None:
                if not loop_warned[0]:
                    print("  [LOOP] UYARI: dump'ta keypoint/descriptor yok "
                          "(frontend.export_keypoints AÇIK + dump'ı yeniden üret) → loop closure devre dışı.")
                    loop_warned[0] = True
            else:
                loop_closer.add_keyframe(f, T_WC, kps_f, desc_f, obs, lm_world)

        if use_fl and (_atlandi or _paralel_bu):
            pose_estimates[f] = T_WC
        _lm_onceki = lm_world
        if use_fl and update_ok and not _atlandi and not (_paralel_bu and _est_hazir is None):
            # pencere-içi çözümden pozları sakla, live_lms / lm_world'ü yeniden kur
            _t0e = time.time()
            try:
                est_fl = _est_hazir if _paralel_bu else smoother.calculateEstimate()
            except IndexError:
                # GTSAM IncrementalFixedLagSmoother hatasi: marjinalize edilen anahtar dogrusallastirma noktasinda
                # kaliyor ama delta'da yok. Kestirim yalniz delta'daki (penceredeki) degiskenlerden kurulur.
                AUDIT['isam2_tahmin_yolu'] += 1
                _is = smoother.getISAM2()
                _th, _dl = _is.getLinearizationPoint(), _is.getDelta()
                est_fl = Values()
                for _k in _dl.keys():
                    _c = gtsam.Symbol(_k).chr()
                    if _c == ord('x'):
                        est_fl.insert(_k, _th.atPose3(_k).retract(_dl.at(_k)))
                    elif _c == ord('l'):
                        est_fl.insert(_k, _th.atPoint3(_k) + _dl.at(_k))
            _zaman['estimate'] += time.time() - _t0e
            if str(f) in os.environ.get('PYECO_ANALIZ', '').split(','):
                # tani: penceredeki grafin yapisi ve dogrusal cozum maliyeti (yalniz gelistirme)
                _fg = smoother.getFactors(); _tip = Counter(); _mk = 0; _mkd = 0
                _temiz = gtsam.NonlinearFactorGraph()
                for _i in range(_fg.size()):
                    _fa = _fg.at(_i)
                    if _fa is None:
                        continue
                    _temiz.add(_fa); _ad = type(_fa).__name__; _tip[_ad] += 1
                    if 'Linear' in _ad:
                        _mk = max(_mk, len(_fa.keys())); _mkd += len(_fa.keys())
                _t1 = time.time(); _lin = _temiz.linearize(est_fl); _t2 = time.time()
                _bn = _lin.eliminateMultifrontal(); _t3 = time.time()
                _lp = gtsam.LevenbergMarquardtParams(); _lp.setMaxIterations(5); _lp.setRelativeErrorTol(1e-3)
                _opt = gtsam.LevenbergMarquardtOptimizer(_temiz, est_fl, _lp); _opt.optimize(); _t4 = time.time()
                print(f"  [analiz] kare {f}: faktor {dict(_tip)}, marjinal faktor anahtar max {_mk} toplam {_mkd}, "
                      f"degisken {est_fl.size()}, dogrusallastir {(_t2-_t1)*1e3:.0f} ms, eleme {(_t3-_t2)*1e3:.0f} ms, "
                      f"LM {(_t4-_t3)*1e3:.0f} ms ({_opt.iterations()} it), smoother.update bu kare {_zaman.get('son', 0)*1e3:.0f} ms")
            live_lms = set()
            live_poses = set()
            lm_world = {}
            resident_lms = set()
            for key in est_fl.keys():
                s = gtsam.Symbol(key); c = s.chr()
                try:
                    if c == ord('x'):
                        pose_estimates[s.index()] = est_fl.atPose3(key).matrix()
                        live_poses.add(s.index())
                    elif c == ord('l'):
                        lid = s.index()
                        resident_lms.add(lid)
                        # kara-listedeki landmark smoother'da dursa da canlı sayılmaz
                        if culler.active and culler.is_bad(lid):
                            continue
                        p = est_fl.atPoint3(key)
                        live_lms.add(lid)
                        lm_world[lid] = np.array([p[0], p[1], p[2]])
                except Exception:
                    pass
            if culler.active:
                culler.on_window(resident_lms)
            if yh_on:
                for lid in live_lms:
                    harita_p[lid] = lm_world[lid]
                for lid in set(frame_lids):
                    harita_n[lid] = harita_n.get(lid, 0) + 1
                    harita_son[lid] = f
            AUDIT['win_lms_sum'] += len(live_lms)
            AUDIT['win_samples'] += 1
            # [TÜRBÜLANS/E5] PENCERE-GÖRELİ marjinal kovaryans (NEES), yalnız --cov_out ile; kestirime dokunmaz.
            # Örnek = bu kestirimdeki EN YENİ poz X(f_n) (anahtar kare kipinde: kendi güncellemesini içeren ilk
            # pencere kestirimindeki anahtar kare; kare kipinde: güncellemeyi kapatan kare). Referans = aynı
            # kestirimdeki en eski poz X(f0) (sıkı önselle gauge). Faktörler: pencereye yazılmış ve bütün
            # anahtarları bu kestirimde duran kare grafları (düzleştiricinin marjinal önselleri HARİÇ).
            # Göreli poz da AYNI kestirimden alınır (kovaryansın doğrusallaştırma noktası).
            if cov_out:
                _cov_n += 1
            if cov_out and (_cov_n % cov_every == 0):
                AUDIT['cov_try'] += 1
                try:
                    keys_w = set(est_fl.keys())
                    _xs = [gtsam.Symbol(k_).index() for k_ in keys_w if gtsam.Symbol(k_).chr() == ord('x')]
                    f_n = max(_xs) if _xs else None
                    f0 = next((ff for (ff, _gg) in cov_hist if X(ff) in keys_w), None)
                    if f_n is None or f0 is None or f0 == f_n:
                        AUDIT['cov_skip_window'] += 1
                        raise _CovSkip()
                    g_c = NonlinearFactorGraph()
                    used = set()
                    for (_ff, gg) in cov_hist:
                        for i_ in range(gg.size()):
                            fac = gg.at(i_)
                            if fac is None:
                                continue
                            fk = list(fac.keys())
                            if all(k_ in keys_w for k_ in fk):
                                g_c.add(fac)
                                used.update(fk)
                    if X(f_n) not in used:
                        AUDIT['cov_skip_no_factor'] += 1
                        raise _CovSkip()
                    g_c.add(gtsam.PriorFactorPose3(X(f0), est_fl.atPose3(X(f0)),
                                                   gtsam.noiseModel.Isotropic.Sigma(6, 1e-6)))
                    used.add(X(f0))
                    vals_c = Values()
                    for k_ in used:
                        if gtsam.Symbol(k_).chr() == ord('x'):
                            vals_c.insert(k_, est_fl.atPose3(k_))
                        else:
                            vals_c.insert(k_, est_fl.atPoint3(k_))
                    # zayıf landmark önseli (10 m): tek gözlemli noktayı iyi koşullu tutar, poz kovaryansını etkilemez
                    for k_ in used:
                        if gtsam.Symbol(k_).chr() != ord('x'):
                            g_c.add(gtsam.PriorFactorPoint3(k_, est_fl.atPoint3(k_), lm_prior))
                    pk = sorted((k_ for k_ in used if gtsam.Symbol(k_).chr() == ord('x')),
                                key=lambda k_: (k_ != X(f_n), k_))
                    lk = sorted(k_ for k_ in used if gtsam.Symbol(k_).chr() != ord('x'))
                    C6 = _pose_marginal_schur(g_c.linearize(vals_c), pk + lk,
                                              [6] * len(pk) + [3] * len(lk), len(pk))
                    if C6 is None:
                        AUDIT['cov_fail_schur'] += 1
                        raise _CovSkip()
                    T_rel = est_fl.atPose3(X(f0)).inverse().compose(est_fl.atPose3(X(f_n))).matrix()
                    cov_records.append((f_n, f0, np.asarray(C6, dtype=np.float64), T_rel))
                except _CovSkip:
                    pass
                except Exception as _e:
                    AUDIT['cov_fail'] += 1
                    AUDIT['cov_fail_' + type(_e).__name__] += 1
                    if AUDIT['cov_fail'] == 1:
                        print(f"  [kare {f}] marjinal kovaryans başarısız (ilk): {type(_e).__name__}: {str(_e)[:300]!r}")
        if _paralel_bu:
            # pencere defteri: uygulanan kestirim bir onceki isten; ucustaki isin ekledigi ve marjinalize edecegi
            # anahtarlar aynadan. Ucustaki yeni noktalarin konumu ekleme anindaki ilk kestirim.
            live_lms, live_poses = set(), set()
            for _k in _ts_ayna:
                _s = gtsam.Symbol(_k)
                if _s.chr() == ord('x'):
                    live_poses.add(_s.index())
                else:
                    live_lms.add(_s.index())
            resident_lms = set(live_lms)
            lm_world = {l_: (lm_world[l_] if l_ in lm_world else _lm_onceki[l_]) for l_ in live_lms
                        if l_ in lm_world or l_ in _lm_onceki}
        # Tam ISAM2 (marjinalizasyon yok): karenin pozu ve bu karede gorulen noktalar her karede okunur; graf hicbir
        # seyi unutmadigi icin kalici haritanin kendisidir. Tam kestirim yine periyodik (asagidaki dal).
        if (not use_fl) and update_ok and int(getattr(args, 'isam_kare_tahmin', 0)):
            _t0e = time.time()
            try:
                pose_estimates[f] = isam.calculateEstimatePose3(X(f)).matrix()
                for lid in set(frame_lids):
                    _p = isam.calculateEstimatePoint3(L(lid))
                    lm_world[lid] = np.array([_p[0], _p[1], _p[2]])
                    if yh_on:
                        harita_p[lid] = lm_world[lid]
                        harita_n[lid] = harita_n.get(lid, 0) + 1
                        harita_son[lid] = f
            except Exception as e_:
                AUDIT['isam_kare_tahmin_hata'] += 1
            _zaman['estimate'] += time.time() - _t0e
        # PnP cache'ini periyodik olarak tam çözümden tazele (keyframe modunda her KF'de)
        if use_fl and update_ok:
            pass
        elif update_ok and (kf_mode or (f % REFRESH_EVERY == 0 and f > 0) or f == n_frames - 1):
            try:
                est_refresh = isam.calculateEstimate()
                for lid in initialized_lms:
                    if culler.active and culler.is_bad(lid):
                        continue
                    if est_refresh.exists(L(lid)):
                        p = est_refresh.atPoint3(L(lid))
                        lm_world[lid] = np.array([p[0], p[1], p[2]])
            except Exception as e:
                print(f"  [kare {f}] cache tazeleme hatası: {e}")

        if f % 100 == 0 or f == n_frames - 1:
            _sim = time.time()
            print(f"  [sure] kare {f}: son 100 kare {(_sim - _son_sure[0]) * 10:.0f} ms/kare "
                  f"(update {_zaman['update'] * 10:.0f}, estimate {_zaman['estimate'] * 10:.0f}), pencere lm {len(live_lms)}, "
                  f"poz {len(live_poses)}, harita {len(harita_p)}, mv tampon {len(pending_mv)}")
            _son_sure[0] = _sim
            if hasattr(smoother if use_fl else None, 'sure') and smoother.sure['n']:
                _s = smoother.sure
                print(f"  [kayan] kur {_s['kur']/_s['n']*1e3:.0f} ms, LM {_s['lm']/_s['n']*1e3:.0f} ms, "
                      f"it {_s['it']/_s['n']:.1f}, faktor {_s['fak']/_s['n']:.0f}, degisken {_s['deg']/_s['n']:.0f}")
                for _k in _s:
                    _s[_k] = 0
            _zaman['update'] = _zaman['estimate'] = 0.0
            print(f"  kare {f}: stereo={n_stereo}, mono={n_proj-n_stereo}, yeni_lm={n_new}, "
                  f"gated={n_gated}, sbp={n_recovered}, toplam_lm={len(initialized_lms)}, "
                  + (f"culled={n_culled_total}, " if culler.active else "")
                  + f"mv_init(toplam)={n_mv_init}")

        # [hakem kıyası] gerçek zamanlı (nedensel) poz: kare f işlendiği andaki kestirim, sonraki
        # pencere güncellemelerinden önce. Yalnız PYECO_NEDENSEL_OUT verilince tutulur; çıktıyı değiştirmez.
        if _NEDENSEL is not None and f in pose_estimates:
            _NEDENSEL[f] = pose_estimates[f].copy()
        if debug_frames and f < debug_frames and f in pose_estimates:
            # teshis: on-tahmin (PnP + motion-only BA) ile pencere sonrasi nedensel poz farki
            print(f"  [dbg f={f:>3}] init_t={np.round(T_WC[:3, 3], 3)} ned_t={np.round(pose_estimates[f][:3, 3], 3)} "
                  f"fark={np.linalg.norm(T_WC[:3, 3] - pose_estimates[f][:3, 3]):.3f} upd_ok={update_ok} "
                  f"stereo={n_stereo} mono={n_proj - n_stereo} yeni={n_new} gated={n_gated}")

    # sayaç raporu (yalnız --audit)
    if _audit:
        AUDIT['mv_init'] = n_mv_init
        AUDIT['n_recover_final'] = n_recover
        AUDIT['keyframes'] = len(keyframes)
        AUDIT['inter_frames'] = len(inter_rel)
        AUDIT['lm_unique_ever'] = len(_lm_ever)
        for _k, _v in culler.report().items():
            AUDIT['cull_' + _k] = _v
        print("\n================ YÖNTEM SAYAÇLARI (audit) ================")
        for k in sorted(AUDIT):
            print(f"  {k:28s} = {AUDIT[k]}")
        if AUDIT.get('win_samples'):
            print(f"  {'win_lms_mean':28s} = {AUDIT['win_lms_sum']/AUDIT['win_samples']:.1f}")
        print("==========================================================\n")

    # son toplu poz okuması yalnız tam-graf iSAM2 için (fixed-lag'de kare kare toplandı)
    if not use_fl:
        est = isam.calculateEstimate()
        for f in range(n_frames):
            if est.exists(X(f)):
                pose_estimates[f] = est.atPose3(X(f)).matrix()

    # Ara kareleri kuşatan iki rafine KF arasında SE3 interpolasyonuyla yeniden kur;
    # son KF'den sonrakiler referans KF'ye göreli odometriyle.
    if inter_rel:
        graph_frames = sorted(fr for fr in pose_estimates if fr not in inter_rel)
        gf_arr = np.asarray(graph_frames)
        for f in sorted(inter_rel.keys()):
            pos = int(np.searchsorted(gf_arr, f))
            a = graph_frames[pos - 1] if pos > 0 else (graph_frames[0] if graph_frames else None)
            b = graph_frames[pos] if pos < len(graph_frames) else None
            if a is None:
                continue
            if b is None or b == a:
                rk = inter_ref.get(f, a)
                pose_estimates[f] = pose_estimates.get(rk, pose_estimates[a]) @ inter_rel[f]
                continue
            t = (f - a) / float(b - a)
            Ta, Tb = pose_estimates[a], pose_estimates[b]
            Ti = np.eye(4)
            key = Rot.from_matrix(np.stack([Ta[:3, :3], Tb[:3, :3]]))
            Ti[:3, :3] = Slerp([0.0, 1.0], key)([t])[0].as_matrix()
            Ti[:3, 3] = (1.0 - t) * Ta[:3, 3] + t * Tb[:3, 3]
            pose_estimates[f] = Ti

    # [TÜRBÜLANS/E5] NEES için pozların loop-closure DÜZELTMESİ ÖNCESİ hâli (kovaryans düzeltmeyi görmez).
    pose_estimates_precorr = {k: v.copy() for k, v in pose_estimates.items()} if cov_out else None

    # Loop closure düzeltmesi: KF başına ΔT = corrected @ inv(old), önceki en yakın KF'nin ΔT'si ara karelere yayılır.
    if loop_closer is not None and loop_closer.n_detected > 0:
        kf_frames, corrected = loop_closer.finalize()
        kf_sorted = sorted(kf_frames)
        old = {fk: pose_estimates[fk].copy() for fk in kf_sorted if fk in pose_estimates}
        deltas = {fk: corrected[fk] @ np.linalg.inv(old[fk])
                  for fk in kf_sorted if fk in old}
        if deltas:
            dk = sorted(deltas.keys())
            dk_arr = np.array(dk)
            for f in list(pose_estimates.keys()):
                pos = np.searchsorted(dk_arr, f, side='right') - 1
                fk = dk[pos] if pos >= 0 else dk[0]
                pose_estimates[f] = deltas[fk] @ pose_estimates[f]
        print(f"[LOOP] Düzeltme uygulandı ({len(deltas)} keyframe, tüm karelere yayıldı).")
    elif loop_closer is not None:
        loop_closer.finalize()

    # Kapsam doğrulaması: eksik trajectory'de ATE yalnız kısa segmentte hesaplanır ve yapay iyi çıkar.
    coverage = len(pose_estimates) / n_frames if n_frames else 1.0
    MIN_COVERAGE = 0.95
    if coverage < MIN_COVERAGE:
        raise RuntimeError(
            f"EKSİK TRAJECTORY: yalnız {len(pose_estimates)}/{n_frames} kare pozlanabildi "
            f"(kapsam={coverage:.0%} < {MIN_COVERAGE:.0%}) → yarım/güvenilmez sonuç, ıraksak sayılıyor")

    frame_ts = np.array(frame_ts, dtype=np.int64)
    if int(getattr(args, 'gt_capa', 1)):
        # gt_capa=1 (eski, yalniz kiyas icin): GT'ye anchor + gövde frame çıktısı (tek katı dönüşüm; SE3/Sim3 hizalamalı ATE buna duyarsız)
        gt_t, gt_pos, gt_quat = load_groundtruth(args.gt)
        covered = np.where(frame_ts >= gt_t[0])[0]
        if len(covered) == 0:
            raise RuntimeError("Hiçbir kare GT zaman aralığında değil!")
        f_a = int(covered[0])
        if f_a not in pose_estimates:
            f_a = min(pose_estimates.keys(), key=lambda k: abs(k - f_a))
        T_WB_gt_fa = gt_pose_at(frame_ts[f_a], gt_t, gt_pos, gt_quat)
        T_W0_C_fa = pose_estimates[f_a]
        # VO-world -> GT-world rigid dönüşümü
        T_W_W0 = T_WB_gt_fa @ T_B_C @ np.linalg.inv(T_W0_C_fa)
    else:
        # GT HIC okunmaz: dunya = ilk pozlanan karedeki govde cercevesi.
        f_a = min(pose_estimates.keys())
        T_W_W0 = T_B_C @ np.linalg.inv(pose_estimates[f_a])

    rows = []
    for f in range(n_frames):
        if f not in pose_estimates:
            continue
        T_WC = T_W_W0 @ pose_estimates[f]
        T_WB = T_WC @ T_C_B
        p = T_WB[:3, 3]
        q_xyzw = Rot.from_matrix(T_WB[:3, :3]).as_quat()
        rows.append([frame_ts[f], p[0], p[1], p[2],
                     q_xyzw[0], q_xyzw[1], q_xyzw[2], q_xyzw[3]])

    out_df = pd.DataFrame(rows, columns=[
        '#timestamp [ns]', ' p_RS_R_x [m]', ' p_RS_R_y [m]', ' p_RS_R_z [m]',
        ' q_RS_x []', ' q_RS_y []', ' q_RS_z []', ' q_RS_w []'])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False)
    if _NEDENSEL:
        rows_n = []
        for f in range(n_frames):
            if f not in _NEDENSEL:
                continue
            T_WB = T_W_W0 @ _NEDENSEL[f] @ T_C_B
            q_xyzw = Rot.from_matrix(T_WB[:3, :3]).as_quat()
            rows_n.append([frame_ts[f], *T_WB[:3, 3], *q_xyzw])
        pd.DataFrame(rows_n, columns=out_df.columns).to_csv(os.environ['PYECO_NEDENSEL_OUT'], index=False)
        print(f"Nedensel trajectory yazıldı: {os.environ['PYECO_NEDENSEL_OUT']} ({len(rows_n)} poz)")
    if kf_mode:
        print(f"ANAHTAR KARE: {AUDIT['kf_anahtar_kare']} anahtar, {AUDIT['kf_ara_kare']} ara kare")
    # [TÜRBÜLANS/E5] NEES kovaryans dökümü. C6: X(f)'in X(f0) sabitken YEREL teğet kovaryansı ([rot, trans],
    # GTSAM retract T*Exp(xi)). rel: aynı pencere kestiriminden T_cam(f0)<-T_cam(f) (NEES bunu kullanır);
    # rel_final: loop-closure düzeltmesi öncesi SON pozlardan (eski hattın tanımı; karşılaştırma için).
    if cov_out and cov_records:
        PE = pose_estimates_precorr if pose_estimates_precorr is not None else pose_estimates
        ts_l, ts0_l, rel_l, relf_l, C_l, f_l = [], [], [], [], [], []
        for f_n, f0, C6, T_rel in cov_records:
            ts_l.append(frame_ts[f_n]); ts0_l.append(frame_ts[f0]); f_l.append((f_n, f0))
            rel_l.append(T_rel); C_l.append(C6)
            relf_l.append(np.linalg.inv(PE[f0]) @ PE[f_n] if (f_n in PE and f0 in PE) else np.full((4, 4), np.nan))
        Path(cov_out).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cov_out, ts=np.array(ts_l, dtype=np.int64), ts0=np.array(ts0_l, dtype=np.int64),
                            frames=np.array(f_l, dtype=np.int64), rel=np.array(rel_l), rel_final=np.array(relf_l),
                            cov=np.array(C_l), T_B_C=np.asarray(T_B_C),
                            mode=np.array('keyframe' if kf_mode else 'frame'))
        _c = np.array(_COND_LOG) if _COND_LOG else np.array([np.nan])
        print(f"NEES kovaryans dökümü: {cov_out} ({len(ts_l)} örnek; deneme {AUDIT.get('cov_try', 0)}, "
              f"hata {AUDIT.get('cov_fail', 0)}, Schur-red {AUDIT.get('cov_fail_schur', 0)}, "
              f"koşul medyan {np.nanmedian(_c):.2e})")
    if turb_builder is not None and str(getattr(args, 'turb_stats_out', '') or ''):
        import json as _json
        Path(args.turb_stats_out).parent.mkdir(parents=True, exist_ok=True)
        _st = turb_builder.summary()
        _st['cov'] = {k_: int(v_) for k_, v_ in AUDIT.items() if k_.startswith('cov_')}
        with open(args.turb_stats_out, 'w') as _fj:
            _json.dump(_st, _fj, indent=1)
    if yh_on:
        print(f"YEREL HARITA: izdusen {AUDIT['yh_izdusen']}, aday {AUDIT['yh_aday_esles']}, kabul {AUDIT['yh_esles']}, "
              f"faktor {AUDIT['yh_faktor']}, grafa geri donen {AUDIT['yh_geri_donen']}, "
              f"yeniden yakalama {AUDIT['yh_yeniden_yakalama']}, oran red {AUDIT['yh_oran_red']}, harita {len(harita_p)} nokta, "
              f"tekil sabitleme {AUDIT['tekil_sabitleme']}, tanimlayici konumlama {AUDIT['yh_konumla_kabul']} "
              f"(aday {AUDIT['yh_konumla_aday']})")
    print(f"\nTrajectory yazıldı: {args.out} ({len(rows)} poz)")
    print(f"Anchor karesi: index {f_a}, ts {frame_ts[f_a]}")

    if not int(getattr(args, 'gt_capa', 1)):
        args.n_frames_ = n_frames
        args.backend_elapsed_s_ = time.time() - _t0_backend
        # GT varsa yalniz rapor: SE(3) hizali ATE + PNG (yorunge dosyasi yazildiktan sonra; kestirimi etkilemez)
        if Path(args.out).stem.startswith('gtsam_traj_') and args.gt and Path(args.gt).exists():
            try:
                from eval_trajectory import evaluate as _plot
                seq = Path(args.out).stem.replace('gtsam_traj_', '')
                return _plot(args.out, args.gt, str(Path(args.out).with_name(f"trajectory_{seq}.png")),
                             dataset_name=seq)[0]
            except Exception as e:
                print(f"[DEGERLENDIRME] atlandi: {e}")
        return None
    # ATE, GT ile örtüşen segmentte
    est_t = out_df['#timestamp [ns]'].values.astype(np.float64)
    est_xyz = out_df[[' p_RS_R_x [m]', ' p_RS_R_y [m]', ' p_RS_R_z [m]']].values
    t_start = max(gt_t[0], est_t[0]); t_end = min(gt_t[-1], est_t[-1])
    valid = (gt_t >= t_start) & (gt_t <= t_end)
    gtt = gt_t[valid].astype(np.float64)
    gxyz = gt_pos[valid]
    interp = np.stack([np.interp(gtt, est_t, est_xyz[:, k]) for k in range(3)], axis=1)
    ate = float(np.sqrt(np.mean(np.sum((gxyz - interp)**2, axis=1))))
    print(f"ATE (RMSE, GT ile örtüşen segment, {valid.sum()} nokta): {ate:.4f} m")

    # FPS raporu için backend süresi ve kare sayısı args'a iliştirilir (PNG çizimi hariç)
    args.n_frames_ = n_frames
    args.backend_elapsed_s_ = time.time() - _t0_backend

    # trajectory PNG'si yalnız gtsam_traj_<seq>.csv çıktısında (tavlama koşumlarında atlanır)
    stem = Path(args.out).stem
    if stem.startswith('gtsam_traj_'):
        try:
            from eval_trajectory import evaluate as _plot
            seq = stem.replace('gtsam_traj_', '')
            png = str(Path(args.out).with_name(f"trajectory_{seq}.png"))
            _plot(args.out, args.gt, png, dataset_name=seq)
        except Exception as e:
            print(f"[PLOT] görsel atlandı: {e}")
    return ate


if __name__ == '__main__':
    main()
