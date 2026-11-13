"""
Loop closure + pose-graph optimization (PGO) + opsiyonel global BA. Varsayılan kapalı (loop.enable).

Harici sözlük yok: keyframe descriptor ortalaması global imza olur, eski KF'lere kosinüs
benzerliğiyle bakılır; adaylar karşılıklı-NN eşleştirme + PnP RANSAC ile doğrulanır. Sekans
sonunda KF pozları üzerinde PGO (odometri + Huber loop kenarları), istenirse GBA.
Gereksinim: frontend.export_keypoints açık (dump'ta kare başına 'kps'/'desc').
"""
import numpy as np
import cv2
import gtsam
from gtsam import Pose3, Point3, Cal3_S2
from gtsam.symbol_shorthand import X, L


def _normalize_rows(D):
    D = np.asarray(D, dtype=np.float64)
    if D.ndim != 2 or D.shape[0] == 0:
        return D.reshape(0, 0) if D.ndim != 2 else D
    return D / np.maximum(np.linalg.norm(D, axis=1, keepdims=True), 1e-12)


def _mutual_nn(desc_a, desc_b, ratio=0.85):
    """Normalize descriptor kümeleri arası karşılıklı en-yakın-komşu eşleşmeleri; (idx_a, idx_b) listesi."""
    if desc_a.shape[0] == 0 or desc_b.shape[0] == 0:
        return []
    S = desc_a @ desc_b.T                      # (Na, Nb) kosinüs benzerliği
    nn_ab = np.argmax(S, axis=1)               # a -> en iyi b
    nn_ba = np.argmax(S, axis=0)               # b -> en iyi a
    out = []
    for ia in range(S.shape[0]):
        ib = int(nn_ab[ia])
        if int(nn_ba[ib]) != ia:               # karşılıklı değil
            continue
        # ratio test: ikinci en iyi b'ye göre
        row = S[ia]
        best = row[ib]
        if row.shape[0] >= 2:
            second = np.partition(row, -2)[-2]
            if best < 1e-9 or (1.0 - best) > ratio * (1.0 - second) + 1e-9:
                # Lowe testi yaklaşıklığı; şu an eleme yapmıyor
                pass
        out.append((ia, ib))
    return out


class LoopCloser:
    """Keyframe tabanlı loop closure + PGO + opsiyonel GBA.
    Backend her KF'de add_keyframe çağırır; sekans sonunda finalize düzeltilmiş KF pozlarını verir."""

    def __init__(self, args, fx, fy, cx, cy, baseline):
        self.args = args
        self.fx, self.fy, self.cx, self.cy = float(fx), float(fy), float(cx), float(cy)
        self.baseline = float(baseline)
        self.K_cv = np.array([[fx, 0, cx], [0, fy, cy], [0, 0, 1]], dtype=np.float64)
        # Keyframe deposu (paralel listeler; index = kf sırası)
        self.kf_f = []          # frame_idx
        self.kf_pose = []       # T_WC (4x4) — commit anındaki tahmin
        self.kf_sig = []        # global appearance imzası (D,)
        self.kf_desc = []       # (M,D) normalize descriptor'lar
        self.kf_kps = []        # (M,2) keypoint xy
        self.kf_lmworld = []    # {lid: world_xyz} snapshot
        self.kf_obs = []        # (N,8) obs satırları (GBA için)
        self.loops = []         # (i, j, T_ij 4x4) — i<j, T_ij = inv(pose_i) @ pose_j (loop ölçümü)
        self.n_detected = 0

    # ----------------------------------------------------------------- ekleme
    def add_keyframe(self, f, T_WC, kps, desc, obs, lm_world):
        """Bir keyframe'i depoya ekle ve eski KF'lerle loop ara (varsa biriktir)."""
        desc = _normalize_rows(desc) if desc is not None else np.zeros((0, 0))
        kps = np.asarray(kps, dtype=np.float64).reshape(-1, 2) if kps is not None else np.zeros((0, 2))
        sig = desc.mean(axis=0) if desc.shape[0] > 0 else None
        if sig is not None:
            sig = sig / max(np.linalg.norm(sig), 1e-12)
        # Landmark dünya konumu snapshot'ı (sonraki drift/marjinalizasyondan etkilenmesin)
        lmw = {}
        for r in obs:
            lid = int(r[0])
            if lid in lm_world:
                lmw[lid] = np.asarray(lm_world[lid], dtype=np.float64)
        idx = len(self.kf_f)
        self.kf_f.append(int(f))
        self.kf_pose.append(np.asarray(T_WC, dtype=np.float64).copy())
        self.kf_sig.append(sig)
        self.kf_desc.append(desc)
        self.kf_kps.append(kps)
        self.kf_lmworld.append(lmw)
        self.kf_obs.append(np.asarray(obs, dtype=np.float64).copy())

        if sig is None:
            return None
        # ---------------- appearance sorgusu ----------------
        gap = int(self.args.loop_min_keyframe_gap)
        if idx < gap:
            return None
        cand_range = range(0, idx - gap + 1)
        sims = []
        for c in cand_range:
            sc = self.kf_sig[c]
            if sc is None:
                continue
            sims.append((float(np.dot(sig, sc)), c))
        if not sims:
            return None
        sims.sort(reverse=True)
        thr = float(self.args.loop_appearance_sim)
        topk = int(self.args.loop_max_candidates)
        for s, c in sims[:topk]:
            if s < thr:
                break
            res = self._verify(idx, c)
            if res is not None:
                T_ij = res
                self.loops.append((c, idx, T_ij))
                self.n_detected += 1
                print(f"  [LOOP] KF#{idx} (kare {self.kf_f[idx]}) <-> KF#{c} (kare {self.kf_f[c]}) "
                      f"appearance={s:.3f} DOĞRULANDI")
                return (c, idx, T_ij)
        return None

    # --------------------------------------------------------- geometrik doğrulama
    def _verify(self, cur, cand):
        """cur KF noktalarını cand KF'nin 3D landmark'larına PnP ile eşle; yeterli inlier varsa
        T_cand_cur (4x4), yoksa None."""
        desc_cur, kps_cur = self.kf_desc[cur], self.kf_kps[cur]
        desc_cand, kps_cand = self.kf_desc[cand], self.kf_kps[cand]
        lmw_cand = self.kf_lmworld[cand]
        if desc_cur.shape[0] == 0 or desc_cand.shape[0] == 0 or len(lmw_cand) == 0:
            return None
        # Aday KF'de keypoint -> landmark eşlemesi (obs(u,v) en yakın kp'ye)
        obs_cand = self.kf_obs[cand]
        kp_to_lid = {}
        if kps_cand.shape[0] > 0:
            for r in obs_cand:
                lid = int(r[0])
                if lid not in lmw_cand:
                    continue
                d2 = (kps_cand[:, 0] - r[1]) ** 2 + (kps_cand[:, 1] - r[2]) ** 2
                j = int(np.argmin(d2))
                if d2[j] <= 2.0 ** 2:
                    kp_to_lid[j] = lid
        if len(kp_to_lid) < self.args.loop_min_inliers:
            return None
        # Descriptor eşleştir (cur <-> cand)
        matches = _mutual_nn(desc_cur, desc_cand)
        obj_pts, img_pts = [], []
        for ia, ib in matches:
            lid = kp_to_lid.get(ib)
            if lid is None:
                continue
            obj_pts.append(lmw_cand[lid])
            img_pts.append(kps_cur[ia])
        if len(obj_pts) < int(self.args.loop_min_inliers):
            return None
        obj_pts = np.asarray(obj_pts, dtype=np.float64).reshape(-1, 1, 3)
        img_pts = np.asarray(img_pts, dtype=np.float64).reshape(-1, 1, 2)
        try:
            ok, rvec, tvec, inliers = cv2.solvePnPRansac(
                obj_pts, img_pts, self.K_cv, None,
                reprojectionError=4.0, iterationsCount=300, flags=cv2.SOLVEPNP_ITERATIVE)
        except cv2.error:
            return None
        if not ok or inliers is None or len(inliers) < int(self.args.loop_min_inliers):
            return None
        R_cw, _ = cv2.Rodrigues(rvec)
        T_CW = np.eye(4); T_CW[:3, :3] = R_cw; T_CW[:3, 3] = tvec.ravel()
        T_pnp = np.linalg.inv(T_CW)                 # cur kamerasının (aday haritasına göre) dünya pozu
        T_cand_cur = np.linalg.inv(self.kf_pose[cand]) @ T_pnp
        return T_cand_cur

    # ----------------------------------------------------------------- PGO
    def run_pgo(self):
        """KF pozları üzerinde poz-graf optimizasyonu; {kf_index: T_WC} döndürür, loop yoksa değişmez."""
        n = len(self.kf_pose)
        if n == 0 or not self.loops:
            return {i: self.kf_pose[i] for i in range(n)}
        graph = gtsam.NonlinearFactorGraph()
        values = gtsam.Values()
        for i in range(n):
            values.insert(X(i), Pose3(self.kf_pose[i]))
        # Gauge: ilk KF'yi sabitle
        graph.add(gtsam.PriorFactorPose3(
            X(0), Pose3(self.kf_pose[0]),
            gtsam.noiseModel.Diagonal.Sigmas(np.array([1e-4] * 6))))
        # Odometri kenarları (ardışık KF'ler arası VO bağıl pozu)
        odo_noise = gtsam.noiseModel.Diagonal.Sigmas(np.array([0.05] * 3 + [0.10] * 3))
        for i in range(1, n):
            rel = np.linalg.inv(self.kf_pose[i - 1]) @ self.kf_pose[i]
            graph.add(gtsam.BetweenFactorPose3(X(i - 1), X(i), Pose3(rel), odo_noise))
        # Loop kenarları (robust Huber — yanlış loop'lara dayanıklılık)
        ls = float(self.args.loop_pgo_sigma)
        loop_noise = gtsam.noiseModel.Robust.Create(
            gtsam.noiseModel.mEstimator.Huber.Create(1.345),
            gtsam.noiseModel.Diagonal.Sigmas(np.array([ls] * 3 + [ls] * 3)))
        for (i, j, T_ij) in self.loops:
            graph.add(gtsam.BetweenFactorPose3(X(i), X(j), Pose3(T_ij), loop_noise))
        try:
            params = gtsam.LevenbergMarquardtParams()
            result = gtsam.LevenbergMarquardtOptimizer(graph, values, params).optimize()
        except Exception as e:
            print(f"  [LOOP] PGO başarısız: {e}")
            return {i: self.kf_pose[i] for i in range(n)}
        return {i: result.atPose3(X(i)).matrix() for i in range(n)}

    # ----------------------------------------------------------------- GBA
    def run_gba(self, pgo_poses):
        """PGO pozlarından başlayarak tüm KF poz + landmark'lar üzerinde global BA; {kf_index: T_WC}."""
        n = len(self.kf_pose)
        K_gtsam = Cal3_S2(self.fx, self.fy, 0.0, self.cx, self.cy)
        K_stereo = gtsam.Cal3_S2Stereo(self.fx, self.fy, 0.0, self.cx, self.cy, self.baseline)
        meas_noise = gtsam.noiseModel.Robust.Create(
            gtsam.noiseModel.mEstimator.Huber.Create(1.345),
            gtsam.noiseModel.Isotropic.Sigma(2, float(self.args.px_noise)))
        stereo_noise = gtsam.noiseModel.Isotropic.Sigma(3, float(self.args.stereo_sigma))
        graph = gtsam.NonlinearFactorGraph()
        values = gtsam.Values()
        # Landmark dünya konumlarını PGO düzeltmesiyle taşı (her KF'nin kendi düzeltmesiyle)
        lm_init = {}
        for i in range(n):
            dT = pgo_poses[i] @ np.linalg.inv(self.kf_pose[i])   # bu KF'nin rijit düzeltmesi
            values.insert(X(i), Pose3(pgo_poses[i]))
            for lid, pw in self.kf_lmworld[i].items():
                if lid not in lm_init:
                    lm_init[lid] = (dT @ np.array([pw[0], pw[1], pw[2], 1.0]))[:3]
        for lid, pw in lm_init.items():
            values.insert(L(lid), Point3(*pw))
        graph.add(gtsam.PriorFactorPose3(
            X(0), Pose3(pgo_poses[0]), gtsam.noiseModel.Diagonal.Sigmas(np.array([1e-4] * 6))))
        max_d = float(self.args.max_init_depth)
        for i in range(n):
            for r in self.kf_obs[i]:
                lid = int(r[0])
                if lid not in lm_init:
                    continue
                u, v = float(r[1]), float(r[2]); Z = r[7]
                stereo_ok = r[4] > 0.5 and r[3] > 0.5 and np.isfinite(Z) and 0 < Z < max_d
                if stereo_ok:
                    uR = u - self.fx * self.baseline / Z
                    graph.add(gtsam.GenericStereoFactor3D(
                        gtsam.StereoPoint2(u, uR, v), stereo_noise, X(i), L(lid), K_stereo))
                else:
                    graph.add(gtsam.GenericProjectionFactorCal3_S2(
                        np.array([u, v]), meas_noise, X(i), L(lid), K_gtsam, False, False))
        try:
            params = gtsam.LevenbergMarquardtParams()
            params.setMaxIterations(int(self.args.loop_gba_iters))
            result = gtsam.LevenbergMarquardtOptimizer(graph, values, params).optimize()
        except Exception as e:
            print(f"  [LOOP] GBA başarısız: {e}")
            return pgo_poses
        return {i: result.atPose3(X(i)).matrix() for i in range(n)}

    # ----------------------------------------------------------------- finalize
    def finalize(self):
        """Sekans sonunda çağrılır; (kf_frame_indices, {frame_index: T_WC_düzeltilmiş}) döndürür."""
        n = len(self.kf_pose)
        if n == 0:
            return [], {}
        print(f"\n[LOOP] Toplam {n} keyframe, {self.n_detected} loop doğrulandı.")
        if not self.loops:
            print("[LOOP] Loop bulunamadı — düzeltme uygulanmıyor.")
            return self.kf_f, {self.kf_f[i]: self.kf_pose[i] for i in range(n)}
        poses = self.kf_pose
        if int(self.args.loop_pgo):
            poses = self.run_pgo()
            print("[LOOP] PGO tamamlandı.")
        else:
            poses = {i: self.kf_pose[i] for i in range(n)}
        if int(self.args.loop_gba):
            poses = self.run_gba(poses)
            print("[LOOP] GBA tamamlandı.")
        # kf_index -> frame_index anahtarına çevir
        return self.kf_f, {self.kf_f[i]: poses[i] for i in range(n)}
