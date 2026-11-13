"""
Yörünge, tutarlılık ve ön-uç metrikleri (rapor §7).

Yörünge metrikleri üç hizalama düzeyinde raporlanır:
    ATE_noalign  hattın kendi GT-anchor'lı çıktısı (iç tutarlılık)
    ATE_se3      Umeyama SE(3) (ölçek sabit)  -> ANA raporlama metriği
    ATE_sim3     Umeyama Sim(3) (ölçek serbest)
RPE 1 s pencerede hesaplanır: ardışık ATE örnekleri güçlü otokorelasyonlu olduğu
için istatistiksel testler RPE üzerinden de raporlanır.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as Rot, Slerp


# ============================================================== yükleme
def load_gt(gt_csv):
    df = pd.read_csv(gt_csv)
    df.columns = [c.strip() for c in df.columns]
    t = df["#timestamp"].values.astype(np.float64)
    xyz = df[["p_RS_R_x [m]", "p_RS_R_y [m]", "p_RS_R_z [m]"]].values.astype(np.float64)
    q_wxyz = df[["q_RS_w []", "q_RS_x []", "q_RS_y []", "q_RS_z []"]].values.astype(np.float64)
    return t, xyz, q_wxyz


def load_traj(traj_csv):
    df = pd.read_csv(traj_csv)
    df.columns = [c.strip() for c in df.columns]
    t = df["#timestamp [ns]"].values.astype(np.float64)
    xyz = df[["p_RS_R_x [m]", "p_RS_R_y [m]", "p_RS_R_z [m]"]].values.astype(np.float64)
    q_xyzw = df[["q_RS_x []", "q_RS_y []", "q_RS_z []", "q_RS_w []"]].values.astype(np.float64)
    return t, xyz, q_xyzw


def associate(est_t, est_xyz, est_q_xyzw, gt_t, gt_xyz, gt_q_wxyz, max_dt_ns=2e7):
    """Kestirim zaman damgalarında GT'yi interpolasyonla eşleştir (poz + yönelim)."""
    lo = max(est_t[0], gt_t[0])
    hi = min(est_t[-1], gt_t[-1])
    m = (est_t >= lo) & (est_t <= hi)
    et = est_t[m]
    if et.size < 2:
        return None
    gp = np.stack([np.interp(et, gt_t, gt_xyz[:, k]) for k in range(3)], axis=1)
    R_gt = Rot.from_quat(gt_q_wxyz[:, [1, 2, 3, 0]])       # wxyz -> xyzw
    slerp = Slerp(gt_t, R_gt)
    gq = slerp(et)
    return et, est_xyz[m], Rot.from_quat(est_q_xyzw[m]), gp, gq


# ============================================================== hizalama
def umeyama(src: np.ndarray, dst: np.ndarray, with_scale: bool):
    """dst ≈ s R src + t  (Umeyama 1991). Dönüş: (R, t, s)."""
    mu_s = src.mean(axis=0)
    mu_d = dst.mean(axis=0)
    S = src - mu_s
    Dm = dst - mu_d
    C = (Dm.T @ S) / src.shape[0]
    U, sig, Vt = np.linalg.svd(C)
    W = np.eye(3)
    if np.linalg.det(U) * np.linalg.det(Vt) < 0:
        W[2, 2] = -1.0
    R = U @ W @ Vt
    if with_scale:
        var_s = np.mean(np.sum(S ** 2, axis=1))
        s = float(np.trace(np.diag(sig) @ W) / max(var_s, 1e-30))
    else:
        s = 1.0
    t = mu_d - s * R @ mu_s
    return R, t, s


def ate(est_xyz, gt_xyz, mode="se3"):
    """ATE RMSE [m]. mode: none | se3 | sim3."""
    if mode == "none":
        err = gt_xyz - est_xyz
    else:
        R, t, s = umeyama(est_xyz, gt_xyz, with_scale=(mode == "sim3"))
        aligned = (s * (R @ est_xyz.T).T) + t
        err = gt_xyz - aligned
    return float(np.sqrt(np.mean(np.sum(err ** 2, axis=1))))


def rpe(est_t, est_xyz, est_R, gt_xyz, gt_R, delta_ns=1e9):
    """Relative Pose Error: delta (varsayılan 1 s) aralıklı bağıl hareketler.

    Dönüş: (rpe_trans_rmse [m], rpe_rot_rmse [deg], n).
    """
    n = len(est_t)
    j = np.searchsorted(est_t, est_t + delta_ns)
    idx = [(i, int(j[i])) for i in range(n) if j[i] < n]
    if not idx:
        return np.nan, np.nan, 0
    te, re_ = [], []
    for i, k in idx:
        d_est_t = est_R[i].inv().apply(est_xyz[k] - est_xyz[i])
        d_gt_t = gt_R[i].inv().apply(gt_xyz[k] - gt_xyz[i])
        te.append(np.linalg.norm(d_est_t - d_gt_t))
        dR = (est_R[i].inv() * est_R[k]).inv() * (gt_R[i].inv() * gt_R[k])
        re_.append(np.rad2deg(np.linalg.norm(dR.as_rotvec())))
    return (float(np.sqrt(np.mean(np.square(te)))),
            float(np.sqrt(np.mean(np.square(re_)))), len(idx))


# ============================================================== birleşik
def trajectory_metrics(traj_csv, gt_csv, divergence_m=5.0,
                       t_min: float | None = None, t_max: float | None = None):
    """Bir yörünge dosyası için tüm metrikler (sözlük).

    t_min/t_max verilirse metrikler yalnız o ZAMAN PENCERESİNDE hesaplanır. Bu,
    rapor §6.0'daki "ortak ısınma geçmişi + 400 karelik değerlendirme penceresi"
    protokolünün uygulanışıdır: arka-uç sekansın başından koşar (adaptif ön-uçta
    taze stereo kareleri seyrek olduğu için pencere ortasından başlatmak mümkün
    değildir — kestirici hiç landmark init edemez), değerlendirme ise yalnız
    pencerede yapılır.
    """
    est_t, est_xyz, est_q = load_traj(traj_csv)
    if t_min is not None or t_max is not None:
        m = np.ones(len(est_t), dtype=bool)
        if t_min is not None:
            m &= est_t >= t_min
        if t_max is not None:
            m &= est_t <= t_max
        est_t, est_xyz, est_q = est_t[m], est_xyz[m], est_q[m]
        if est_t.size < 5:
            return {"valid": False}
    gt_t, gt_xyz, gt_q = load_gt(gt_csv)
    a = associate(est_t, est_xyz, est_q, gt_t, gt_xyz, gt_q)
    if a is None:
        return {"valid": False}
    et, exyz, eR, gxyz, gR = a
    out = {
        "valid": True,
        "n_poses": int(len(est_t)),
        "n_compare": int(len(et)),
        "ate_noalign": ate(exyz, gxyz, "none"),
        "ate_se3": ate(exyz, gxyz, "se3"),
        "ate_sim3": ate(exyz, gxyz, "sim3"),
        "path_length_gt": float(np.sum(np.linalg.norm(np.diff(gxyz, axis=0), axis=1))),
    }
    rt, rr, nr = rpe(et, exyz, eR, gxyz, gR)
    out.update({"rpe_trans": rt, "rpe_rot_deg": rr, "rpe_n": nr})
    # Sim(3) ölçek katsayısı — ölçek kaybı göstergesi
    _, _, s = umeyama(exyz, gxyz, with_scale=True)
    out["sim3_scale"] = float(s)
    # ıraksama: hizalanmış konum hatasının maksimumu
    R, t, _ = umeyama(exyz, gxyz, with_scale=False)
    al = (R @ exyz.T).T + t
    e = np.linalg.norm(gxyz - al, axis=1)
    out["max_err"] = float(e.max())
    out["diverged"] = bool(e.max() > divergence_m)
    return out


# ============================================================== NEES
# TEK KAYNAK: NEES tutarlilik bandi. Bar-Shalom 2001 §5.4, TEK ORNEK icin iki
# yanli %95 olasilik bolgesi: [chi2.ppf(0.025, dof), chi2.ppf(0.975, dof)].
# (2026-08-21: metrics %95 kullanirken analyze ve sanity %90 [0.35, 7.81]
# kullaniyordu; ayni tabloda iki farkli bant tanimi vardi. Hicbir kapinin
# sonucunu degistirmiyordu -- olculen degerler 11-661 arasi, iki bandin da cok
# disinda -- ama makalede beyan edilen bant TEK olmali.)
def nees_band(dof: int = 3) -> tuple[float, float]:
    from scipy.stats import chi2
    return float(chi2.ppf(0.025, dof)), float(chi2.ppf(0.975, dof))


NEES_BAND = (0.21579528262389785, 9.348403604496148)   # dof=3


def _gt_logmap(E: np.ndarray) -> np.ndarray:
    """SE(3) logaritması, GTSAM'in Pose3 teğet konvansiyonunda: [omega(3), v(3)].

    Kovaryansı üreten `linearize` ile BİREBİR aynı konvansiyon olsun diye
    kasıtlı olarak GTSAM'e delege edilir (bkz. nees_from_dump belgesi)."""
    import gtsam
    return np.asarray(gtsam.Pose3.Logmap(gtsam.Pose3(np.asarray(E, dtype=np.float64))),
                      dtype=np.float64)


def _T_from(p, q_wxyz):
    T = np.eye(4)
    T[:3, :3] = Rot.from_quat(q_wxyz[[1, 2, 3, 0]]).as_matrix()
    T[:3, 3] = p
    return T


def nees_from_dump(nees_npz, gt_csv, dof=3,
                   t_min: float | None = None, t_max: float | None = None,
                   rel_key: str = "rel"):
    """PENCERE-GÖRELİ NEES (konum ve yönelim ayrı) — rapor E5.

    gtsam_backend `--cov_out` dosyası:
        ts, ts0   pencerenin son ve ilk karesinin zaman damgaları
        rel       (N,4,4) kestirilen göreli KAMERA pozu  T_cam(f0) <- T_cam(f)
        cov       (N,6,6) X(f)'in X(f0)'a göre YEREL teğet kovaryansı [rot, trans]
        T_B_C     gövde <- rektifiye kamera dönüşümü

    Hata tanımı (GTSAM konvansiyonu, retract T <- T*Exp(xi)):
        xi_err = Log( T_rel_est^{-1} * T_rel_gt )
    NEES = xi_err^T C^{-1} xi_err; konum ve yönelim blokları ayrı raporlanır.
    Tutarlılık bandı: 3 serbestlik derecesi için chi2 %95 aralığı.

    TEĞET UZAYI (2026-08-21 düzeltmesi — KRİTİK):
    Kovaryans, arka-uçta `g_c.linearize(...)` ile üretilir; yani GTSAM'in
    Pose3::retract konvansiyonundaki teğet koordinatlardadır. Bu build'de
    ÖLÇÜLDÜ: `T.retract(xi) == T.compose(Pose3.Expmap(xi))` — yani teğet TAM
    SE(3) üsteli, ötelenme bileşeni v ile t = V(omega)·v ilişkisiyle bağlı.
    Bu kod eskiden hata olarak HAM ötelenmeyi (E[:3,3] = t) alıyordu; o vektör
    kovaryansın yaşadığı uzayda DEĞİLDİR. Sapma ölçüldü:
        |omega| = 0.1 rad -> %3.9 · 0.3 -> %11.9 · 0.6 -> %24.1 · 1.0 -> %40.2
    Pencere 80 kare (4 s, EuRoC MH) olduğundan |omega| rutin olarak 0.3-1.0 rad
    aralığındadır, yani sapma ihmal edilebilir DEĞİLDİR ve dönme miktarına
    bağlı olduğu için pencereden pencereye ve KOLDAN KOLA farklı davranır
    (kollar farklı yörünge kestirir) — yani kol karşılaştırmasını kirletir.
    Çözüm: hatayı `gtsam.Pose3.Logmap` ile al. Kendi kapalı-form V^{-1}
    uygulamamızı yazmak yerine GTSAM'i çağırmak, konvansiyonun kovaryansla
    TANIMI GEREĞİ aynı kalmasını garanti eder (build bayrağı değişse bile).

    PYECO (pyeco_turb, 2026-11-07): `rel` artık kovaryansla AYNI pencere kestiriminden
    (doğrusallaştırma noktası) alınır; eski tanım (son, düzeltme öncesi pozlar) dökümde
    `rel_final` olarak durur ve `rel_key="rel_final"` ile seçilebilir. Döküm `mode`
    alanı örneklerin anahtar kare mi kare mi olduğunu söyler.
    """
    d = np.load(nees_npz)
    ts, ts0, rel, C = d["ts"], d["ts0"], d[rel_key], d["cov"]
    T_B_C = d["T_B_C"] if "T_B_C" in d else np.eye(4)
    T_C_B = np.linalg.inv(T_B_C)
    gt_t, gt_xyz, gt_q = load_gt(gt_csv)
    m = (ts >= gt_t[0]) & (ts <= gt_t[-1]) & (ts0 >= gt_t[0]) & (ts0 <= gt_t[-1])
    if t_min is not None:
        m &= ts >= t_min
    if t_max is not None:
        m &= ts <= t_max
    m &= np.all(np.isfinite(rel.reshape(len(rel), -1)), axis=1)
    ts, ts0, rel, C = ts[m], ts0[m], rel[m], C[m]
    if ts.size == 0:
        return {"valid": False}

    R_gt = Rot.from_quat(gt_q[:, [1, 2, 3, 0]])
    slerp = Slerp(gt_t, R_gt)

    def gt_pose(t):
        p = np.array([np.interp(t, gt_t, gt_xyz[:, k]) for k in range(3)])
        T = np.eye(4)
        T[:3, :3] = slerp([t])[0].as_matrix()
        T[:3, 3] = p
        return T

    from scipy.stats import chi2
    lo, hi = chi2.ppf(0.025, dof), chi2.ppf(0.975, dof)
    nees_p, nees_r = [], []
    for i in range(len(ts)):
        T_WB0 = gt_pose(float(ts0[i]))
        T_WB1 = gt_pose(float(ts[i]))
        # GT göreli KAMERA pozu
        T_rel_gt = T_C_B @ np.linalg.inv(T_WB0) @ T_WB1 @ T_B_C
        E = np.linalg.inv(rel[i]) @ T_rel_gt
        # GTSAM teğetinde tam SE(3) logaritması -> [omega(3), v(3)]
        xi = _gt_logmap(E)
        xi_r, xi_t = xi[:3], xi[3:]
        Cr, Cp = C[i][:3, :3], C[i][3:, 3:]
        try:
            nees_p.append(float(xi_t @ np.linalg.solve(Cp, xi_t)))
            nees_r.append(float(xi_r @ np.linalg.solve(Cr, xi_r)))
        except np.linalg.LinAlgError:
            continue
    nees_p = np.array(nees_p)
    nees_r = np.array(nees_r)
    return {
        "valid": bool(nees_p.size),
        "n": int(nees_p.size),
        "nees_pos_mean": float(np.mean(nees_p)) if nees_p.size else np.nan,
        "nees_pos_median": float(np.median(nees_p)) if nees_p.size else np.nan,
        "nees_rot_mean": float(np.mean(nees_r)) if nees_r.size else np.nan,
        "nees_rot_median": float(np.median(nees_r)) if nees_r.size else np.nan,
        "coverage_pos": float(np.mean((nees_p >= lo) & (nees_p <= hi))) if nees_p.size else np.nan,
        "coverage_rot": float(np.mean((nees_r >= lo) & (nees_r <= hi))) if nees_r.size else np.nan,
        "band": [float(lo), float(hi)],
    }


# ============================================================== görüntü kalitesi
def psnr(a: np.ndarray, b: np.ndarray) -> float:
    a = a.astype(np.float64); b = b.astype(np.float64)
    mse = np.mean((a - b) ** 2)
    return float(10 * np.log10(255.0 ** 2 / mse)) if mse > 0 else float("inf")


def ssim(a: np.ndarray, b: np.ndarray) -> float:
    """Gauss pencereli SSIM (skimage bağımlılığı olmadan)."""
    import cv2
    a = a.astype(np.float64); b = b.astype(np.float64)
    C1, C2 = (0.01 * 255) ** 2, (0.03 * 255) ** 2
    mu_a = cv2.GaussianBlur(a, (11, 11), 1.5)
    mu_b = cv2.GaussianBlur(b, (11, 11), 1.5)
    saa = cv2.GaussianBlur(a * a, (11, 11), 1.5) - mu_a ** 2
    sbb = cv2.GaussianBlur(b * b, (11, 11), 1.5) - mu_b ** 2
    sab = cv2.GaussianBlur(a * b, (11, 11), 1.5) - mu_a * mu_b
    num = (2 * mu_a * mu_b + C1) * (2 * sab + C2)
    den = (mu_a ** 2 + mu_b ** 2 + C1) * (saa + sbb + C2)
    return float(np.mean(num / den))
