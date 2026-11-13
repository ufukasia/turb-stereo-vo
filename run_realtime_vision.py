"""
Komple gerçek zamanlı VIO: ön uç ve arka uç tek döngüde, canlı görüntü + Pangolin yörünge.
Dump kullanmaz; her kare canlı işlenir ve anında iSAM2 arka ucuna verilir, ölçülen hız gerçek sistem hızıdır.
Çıkış: Live Tracking penceresinde 'q' ya da Ctrl-C. Sonunda trajectory CSV, ATE ve hız raporu yazılır.

    python run_realtime_vision.py                 # config.yaml'daki sekans
    python run_realtime_vision.py --max_frames 300
"""
import _pango_preload
_pango_preload.ensure()  # kütüphane çakışmasını gerekirse çöz; torch'tan önce çağrılmalı

import csv
import time
import argparse
from pathlib import Path
import numpy as np
import cv2
import gtsam
import gtsam_unstable as gu  # BatchFixedLagSmoother için
from gtsam import (Pose3, Point3, Cal3_S2, ISAM2, ISAM2Params,
                   NonlinearFactorGraph, Values)
from gtsam.symbol_shorthand import X, L
from scipy.spatial.transform import Rotation as Rot

import config
from stereo_front_end import StereoMatcher
from gtsam_backend import (load_groundtruth, gt_pose_at, T_body_rectcam0,
                           refine_pose_motion_only)

# Pangolin isteğe bağlı; yoksa yörünge penceresi atlanır
try:
    import pypangolin as pango
    from OpenGL.GL import (glClear, glClearColor, glEnable, glLineWidth, glColor3f,
                           glBegin, glEnd, glVertex3f, glPointSize,
                           GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, GL_DEPTH_TEST,
                           GL_LINE_STRIP, GL_LINES, GL_POINTS)
    _PANGO = True
except Exception as _e:
    _PANGO = False
    print(f"Pangolin yüklenemedi (yörünge penceresi atlanır): {_e}")


class Cv2TrajView:
    """Pangolin yoksa devreye giren yedek yörünge penceresi (OpenCV, tepeden X-Y görünümü).
    OpenCV zaten zorunlu bağımlılık olduğu için her kurulumda çalışır."""

    def __init__(self, gt_xyz, w=760, h=760, title="Trajectory (yedek görünüm) — GT kırmızı / kestirim mavi"):
        self.w, self.h, self.title = int(w), int(h), title
        g = np.asarray(gt_xyz, dtype=np.float64)[:, :2]
        lo, hi = g.min(0), g.max(0)
        pad = 0.08 * max(float((hi - lo).max()), 1e-6)
        self.lo, self.hi = lo - pad, hi + pad
        span = np.maximum(self.hi - self.lo, 1e-6)
        self.s = min((self.w - 20) / span[0], (self.h - 20) / span[1])   # eşit ölçek (şekil bozulmasın)
        self.gt_px = self._to_px(g[::5])
        cv2.namedWindow(self.title, cv2.WINDOW_NORMAL)
        cv2.resizeWindow(self.title, self.w, self.h)

    def _to_px(self, xy):
        p = (np.asarray(xy, dtype=np.float64)[:, :2] - self.lo) * self.s
        px = np.empty((len(p), 2), np.int32)
        px[:, 0] = np.clip(p[:, 0] + 10, 0, self.w - 1)
        px[:, 1] = np.clip(self.h - 10 - p[:, 1], 0, self.h - 1)      # y yukarı
        return px

    def render(self, traj):
        img = np.full((self.h, self.w, 3), 255, np.uint8)
        cv2.polylines(img, [self.gt_px], False, (40, 40, 210), 2, cv2.LINE_AA)     # GT (kırmızı)
        if traj is not None and len(traj) >= 2:
            tp = self._to_px(np.asarray(traj))
            cv2.polylines(img, [tp], False, (210, 90, 30), 2, cv2.LINE_AA)         # kestirim (mavi)
            cv2.circle(img, tuple(tp[0]), 6, (40, 170, 40), -1, cv2.LINE_AA)       # başlangıç
            cv2.circle(img, tuple(tp[-1]), 7, (0, 200, 240), -1, cv2.LINE_AA)      # anlık uç
            cv2.putText(img, f"{len(traj)} poz", (12, self.h - 14),
                        cv2.FONT_HERSHEY_SIMPLEX, .55, (90, 90, 90), 1, cv2.LINE_AA)
        cv2.putText(img, "GT (kirmizi)  vs  kestirim (mavi)", (12, 24),
                    cv2.FONT_HERSHEY_SIMPLEX, .6, (60, 60, 60), 1, cv2.LINE_AA)
        cv2.imshow(self.title, img)


def _draw_poly(pts, color, width):
    if len(pts) < 2:
        return
    glLineWidth(width); glColor3f(*color); glBegin(GL_LINE_STRIP)
    for p in pts:
        glVertex3f(float(p[0]), float(p[1]), float(p[2]))
    glEnd()


def _draw_pt(p, color, size):
    glPointSize(size); glColor3f(*color); glBegin(GL_POINTS)
    glVertex3f(float(p[0]), float(p[1]), float(p[2])); glEnd()


class RealtimeBackend:
    """Kare başına artımlı BA (iSAM2 ya da fixed-lag); gtsam_backend ile aynı mantık, canlı gözlemle beslenir.
    step() bir kareyi işler, finalize() GT'ye ankorlar ve ATE döndürür."""

    def __init__(self, calib, a):
        fx, fy, cx, cy = calib['rect_K']
        self.fx, self.fy, self.cx, self.cy = fx, fy, cx, cy
        self.baseline = float(calib['baseline'])
        self.T_B_C = T_body_rectcam0(calib['T_BS'], calib['R1'])
        self.T_C_B = np.linalg.inv(self.T_B_C)
        self.K_gtsam = Cal3_S2(fx, fy, 0.0, cx, cy)
        self.K_stereo = gtsam.Cal3_S2Stereo(fx, fy, 0.0, cx, cy, self.baseline)
        self.K_cv = np.array([[fx, 0, cx], [0, fy, cy], [0, 0, 1]], dtype=np.float64)
        self.a = a

        self.meas_noise = gtsam.noiseModel.Robust.Create(
            gtsam.noiseModel.mEstimator.Huber.Create(1.345),
            gtsam.noiseModel.Isotropic.Sigma(2, a.px_noise))
        self.stereo_noise = gtsam.noiseModel.Isotropic.Sigma(3, a.stereo_sigma)
        self.pose0_prior = gtsam.noiseModel.Diagonal.Sigmas(np.array([1e-4] * 6))
        self.lm_prior = gtsam.noiseModel.Isotropic.Sigma(3, a.landmark_prior_sigma)
        self.motion_noise = gtsam.noiseModel.Diagonal.Sigmas(
            np.array([a.motion_rot_sigma] * 3 + [a.motion_trans_sigma] * 3))

        # fixed_lag>0 ise gtsam_backend.py ile aynı BatchFixedLagSmoother; aksi hâlde tam-graf ISAM2.
        self.use_fl = int(getattr(a, 'fixed_lag', 0)) > 0
        if self.use_fl:
            # LM erken durma: sıcak başlangıçlı pencerede az iterasyon yeter (gtsam_backend.py ile aynı).
            lm_params = gtsam.LevenbergMarquardtParams()
            lm_params.setMaxIterations(int(getattr(a, 'fl_lm_max_iters', 5)))
            lm_params.setRelativeErrorTol(float(getattr(a, 'fl_lm_rel_err_tol', 1e-3)))
            self.smoother = gu.BatchFixedLagSmoother(float(a.fixed_lag), lm_params)
            print(f"[REALTIME] FIXED-LAG smoother AÇIK: pencere={int(a.fixed_lag)} kare "
                  f"(batch, LM erken-durma — gtsam_backend.py ile aynı; config.yaml'daki fixed_lag'e sadık).")
        else:
            params = ISAM2Params()
            params.setRelinearizeThreshold(float(getattr(a, 'relinearize_threshold', 0.01)))
            params.relinearizeSkip = a.relinearize_skip
            self.isam = ISAM2(params)
        self.refresh_every = int(getattr(a, 'refresh_every', 20))

        self.initialized = set()   # tam-graf modu: grafa eklenmiş landmark ID'leri
        self.live_lms = set()      # fixed-lag modu: şu an pencerede olan landmark ID'leri
        self.lm_world = {}
        # Landmark ayıklama; gtsam_backend.py ile aynı sınıf ve kancalar.
        from landmark_culling import LandmarkCuller
        self.culler = LandmarkCuller(a, fx, fy, cx, cy)
        self.resident_lms = set()  # fixed-lag: smoother'da fiilen duran landmark değişkenleri
        self.n_culled_total = 0
        if self.culler.active:
            print(f"[REALTIME] LANDMARK CULLING AÇIK: {self.culler.describe()}")
        self.prev_pose = np.eye(4)
        self.prev_rel = np.eye(4)
        self.pose_estimates = {}
        self.frame_ts = []
        self.f = 0
        # CV factor son başarıyla commit edilen düğüme bağlanır; başarısız karenin X(f)'i grafta yoktur.
        self.prev_graph_f = 0
        self.prev_graph_pose = np.eye(4)
        self.consec_err = 0
        self.total_err = 0
        # Sağlamlık kapıları (gtsam_backend.py ile aynı; 0 = kapalı):
        self.min_inl = int(getattr(a, 'pnp_min_inliers', 0) or 0)
        self.jump_gate = float(getattr(a, 'jump_max_trans', 0.0) or 0.0)
        self.chi2_gate = int(getattr(a, 'chi2_gate_graph', 0) or 0) > 0
        # Sıçrama kapısı koşulsuz uygulanırsa doğru PnP sonsuza dek reddedilebilir; en fazla
        # JUMP_MAX_CONSEC kare (güçlü konsensüste 1 kare) reddedilir, sonra ölçüme dönülür.
        self.JUMP_MAX_CONSEC = 5
        self.PNP_STRONG_RATIO = 0.6
        self.jump_rejects = 0
        if self.min_inl > 0 or self.jump_gate > 0 or self.chi2_gate:
            print(f"[REALTIME] SAĞLAMLIK KAPILARI: pnp_min_inliers={self.min_inl}, "
                  f"jump_max_trans={self.jump_gate}m, chi2_gate_graph={'AÇIK' if self.chi2_gate else 'kapalı'}")

        # Loop closure + PGO (config loop.enable). Tespitler konsola yazılır; düzeltme finalize'da global uygulanır.
        self.loop_closer = None
        self._loop_warned = False
        self._last_obs = None
        self.loop_stride = 3   # her kareyi appearance DB'sine atma (maliyet)
        if int(getattr(a, 'loop_enable', 0)):
            from loop_closing import LoopCloser
            self.loop_closer = LoopCloser(a, fx, fy, cx, cy, self.baseline)
            print("LOOP CLOSURE AÇIK (realtime): " + ("PGO" if int(a.loop_pgo) else "PGO-yok") +
                  (" + GBA" if int(a.loop_gba) else "") +
                  f" | appearance_sim>{a.loop_appearance_sim}, min_inliers={a.loop_min_inliers}, "
                  f"min_gap={a.loop_min_keyframe_gap} KF (düzeltme finalize'da uygulanır)")

    def loop_add(self, kps, desc):
        """Bu karenin keypoint+descriptor'larını loop DB'sine gönder (stride'lı); veri yoksa bir kez uyarır."""
        if self.loop_closer is None:
            return
        f = self.f - 1   # step() içinde işlenmiş son kare
        if f < 0 or f % self.loop_stride != 0:
            return
        if kps is None or desc is None or len(kps) == 0:
            if not self._loop_warned:
                print("  [LOOP] UYARI: front-end keypoint/descriptor vermiyor "
                      "(frontend.export_keypoints AÇIK olmalı) → loop closure boş kalır.")
                self._loop_warned = True
            return
        try:
            self.loop_closer.add_keyframe(f, self.prev_pose, kps, desc, self._last_obs, self.lm_world)
        except Exception as e:
            if not self._loop_warned:
                print(f"  [LOOP] add_keyframe hatası: {e}")
                self._loop_warned = True

    def step(self, timestamp, obs):
        f = self.f
        self.frame_ts.append(int(timestamp))
        a = self.a
        # Disparite-bias ölçek düzeltmesi (δ px): d'=d+δ, 3D nokta ışın boyunca ölçeklenir. 0.0 = kapalı.
        _db = float(getattr(a, 'disparity_bias', 0.0) or 0.0)
        if _db != 0.0:
            _fxb = self.fx * self.baseline
            for _r in obs:
                _Z = _r[7]
                if _r[3] > 0.5 and np.isfinite(_Z) and _Z > 0:
                    _dp = _fxb / _Z + _db
                    if _dp > 1e-6:
                        _k = (_fxb / _dp) / _Z
                        _r[5] *= _k; _r[6] *= _k; _r[7] *= _k
        graph = NonlinearFactorGraph(); values = Values()
        chi2_skip = set()   # bu karede χ² ile elenen gözlemler (chi2_gate_graph)
        new_lids = []       # bu karede grafa ilk kez eklenen landmark'lar (hatada geri alınır)
        resync = False      # sıçrama kapısı bu karede ölçüme geri döndü mü

        if f == 0:
            T_WC = np.eye(4)
            graph.add(gtsam.PriorFactorPose3(X(0), Pose3(T_WC), self.pose0_prior))
            values.insert(X(0), Pose3(T_WC))
        else:
            obj, img, ids_used = [], [], []
            for r in obs:
                lid = int(r[0])
                if lid in self.lm_world:
                    p = self.lm_world[lid]; obj.append([p[0], p[1], p[2]]); img.append([r[1], r[2]])
                    ids_used.append(lid)
            T_pred = (self.prev_pose @ self.prev_rel).copy()   # sabit-hız önkestirim
            T_WC = T_pred
            pnp_strong = False
            if len(obj) >= 6:
                ok, rvec, tvec, inliers = cv2.solvePnPRansac(
                    np.asarray(obj, np.float64).reshape(-1, 1, 3),
                    np.asarray(img, np.float64).reshape(-1, 1, 2),
                    self.K_cv, None, reprojectionError=a.pnp_reproj, iterationsCount=200,
                    flags=cv2.SOLVEPNP_ITERATIVE)
                # PnP RANSAC dışı kalan landmark'lar strike alır (cull_pnp_strikes)
                if self.culler.pnp_strikes > 0 and ok and inliers is not None:
                    _inl = set(int(i) for i in np.asarray(inliers).ravel())
                    self.culler.on_pnp_outliers([ids_used[i] for i in range(len(ids_used)) if i not in _inl])
                # az inlier'lı PnP pozu güvenilmez; öngörüde kal (pnp_min_inliers)
                if ok and rvec is not None and (self.min_inl <= 0 or
                                                (inliers is not None and len(inliers) >= self.min_inl)):
                    Rcw, _ = cv2.Rodrigues(rvec)
                    Tcw = np.eye(4); Tcw[:3, :3] = Rcw; Tcw[:3, 3] = tvec.ravel()
                    T_WC = np.linalg.inv(Tcw)
                    pnp_strong = (inliers is not None and
                                  len(inliers) >= max(self.min_inl, 6,
                                                      int(self.PNP_STRONG_RATIO * len(obj))))
            if a.motion_only_ba:
                # χ² outlier turları için args/fy/cx/cy geçilir (gtsam_backend.py ile aynı çağrı).
                T_WC, chi2_skip = refine_pose_motion_only(
                    T_WC, obs, self.lm_world, self.K_gtsam,
                    self.K_stereo, self.fx, self.baseline,
                    self.meas_noise, self.stereo_noise, a.max_init_depth,
                    args=a, fy=self.fy, cx=self.cx, cy=self.cy, return_outliers=True)
                self.culler.on_chi2_outliers(chi2_skip)   # χ² outlier'lar strike alır
            # Sıçrama kapısı: tekil spike reddedilir; ısrarcı ve güçlü konsensüslü ölçüm re-sync ile kazanır.
            if self.jump_gate > 0:
                if float(np.linalg.norm(T_WC[:3, 3] - T_pred[:3, 3])) > self.jump_gate:
                    if (pnp_strong and self.jump_rejects >= 1) or self.jump_rejects >= self.JUMP_MAX_CONSEC:
                        resync = True          # ölçüme geri dön
                        self.jump_rejects = 0
                    else:
                        T_WC = T_pred
                        self.jump_rejects += 1
                else:
                    self.jump_rejects = 0
            values.insert(X(f), Pose3(T_WC))
            if a.cv_motion:
                # Factor son başarılı commit'e bağlanır; re-sync'te bozuk prev_rel yerine gerçek bağıl hareket kullanılır.
                if self.prev_graph_f == f - 1 and not resync:
                    rel_cv = self.prev_rel
                else:
                    rel_cv = np.linalg.inv(self.prev_graph_pose) @ T_WC
                graph.add(gtsam.BetweenFactorPose3(X(self.prev_graph_f), X(f), Pose3(rel_cv), self.motion_noise))

        T_CW_cur = np.linalg.inv(T_WC)
        active_set = self.live_lms if self.use_fl else self.initialized
        frame_lids = []  # fixed-lag zaman-damgası tazelemesi için bu karede factor alan lid'ler
        for r in obs:
            lid = int(r[0]); u, v = float(r[1]), float(r[2]); Z = r[7]
            has_depth = r[3] > 0.5; fresh = r[4] > 0.5
            stereo_ok = fresh and has_depth and np.isfinite(Z) and 0 < Z < a.max_init_depth
            # χ² outlier grafa girmez (chi2_gate_graph)
            if self.chi2_gate and lid in chi2_skip:
                continue
            # Kara listedeki landmark factor almaz; yalnız allow_reinit + taze stereo + grafta yoksa yeniden yaratılır.
            if self.culler.active and self.culler.is_bad(lid):
                _resident = (lid in self.resident_lms) if self.use_fl else (lid in self.initialized)
                if not (self.culler.allow_reinit and stereo_ok and not _resident):
                    continue
            if lid not in active_set:
                if stereo_ok and np.all(np.isfinite(r[5:8])) and not (self.use_fl and lid in self.resident_lms):
                    pw = (T_WC @ np.array([r[5], r[6], r[7], 1.0]))[:3]
                    values.insert(L(lid), Point3(*pw))
                    self.culler.on_init(lid, f)
                    # Fixed-lag'de prior eklenmez: stereo factor tek gözlemde tam kısıtlar, prior marjinalizasyonda sarkar.
                    if not self.use_fl:
                        graph.add(gtsam.PriorFactorPoint3(L(lid), Point3(*pw), self.lm_prior))
                    active_set.add(lid); self.lm_world[lid] = pw
                    new_lids.append(lid)
                else:
                    continue
            if a.outlier_gate_px > 0 and lid in self.lm_world:
                pw = self.lm_world[lid]; pc = T_CW_cur @ np.array([pw[0], pw[1], pw[2], 1.0])
                if pc[2] > 1e-3:
                    du = (self.fx * pc[0] / pc[2] + self.cx) - u
                    dv = (self.fy * pc[1] / pc[2] + self.cy) - v
                    if du * du + dv * dv > a.outlier_gate_px ** 2:
                        continue
            if stereo_ok:
                uR = u - (self.fx * self.baseline / Z)
                graph.add(gtsam.GenericStereoFactor3D(
                    gtsam.StereoPoint2(u, uR, v), self.stereo_noise, X(f), L(lid), self.K_stereo))
            else:
                graph.add(gtsam.GenericProjectionFactorCal3_S2(
                    np.array([u, v]), self.meas_noise, X(f), L(lid), self.K_gtsam, False, False))
            frame_lids.append(lid)
            self.culler.on_observation(lid)

        # Görünürlük sayımı + ayıklama kararı; update öncesi (gtsam_backend.py ile aynı sıra).
        if self.culler.active:
            self.culler.visibility_pass(T_WC, self.lm_world, active_set, set(int(r[0]) for r in obs))
            _culled = self.culler.evaluate(f, active_set)
            if _culled:
                self.culler.apply_cull(_culled, self.lm_world, self.live_lms)
                self.n_culled_total += len(_culled)

        update_ok = True
        try:
            if self.use_fl:
                ts_map = gu.FixedLagSmootherKeyTimestampMap()
                ts_map.insert((X(f), float(f)))
                for lid in frame_lids:
                    ts_map.insert((L(lid), float(f)))
                self.smoother.update(graph, values, ts_map)
            else:
                self.isam.update(graph, values)
                for _ in range(a.extra_updates):
                    self.isam.update()
            self.consec_err = 0
        except Exception as e:
            update_ok = False
            # Commit edilemeyen yeni landmark'ları geri al; yoksa sonraki kareler var olmayan L(lid)'e factor ekler.
            for lid in new_lids:
                self.culler.forget(lid)
                self.initialized.discard(lid)
                self.live_lms.discard(lid)
                self.lm_world.pop(lid, None)
            self.consec_err += 1
            self.total_err += 1
            if self.total_err <= 5:
                print(f"  [kare {f}] iSAM2/smoother: {e}")
            elif self.total_err == 6:
                print("  [kare {}] ... (sonraki update hataları susturuldu)".format(f))
            if self.consec_err == 20:
                print(f"  ⚠ [kare {f}] 20 ARDIŞIK update hatası — backend toparlanamıyor, "
                      f"trajectory ölü-hesap (PnP/CV) modunda ilerliyor (toplam {self.total_err} hata).")
        # Bozuk bağıl hareket CV durumuna taşınmaz; re-sync'te hız durumu sıfırlanır.
        if f > 0:
            if resync:
                self.prev_rel = np.eye(4)
            else:
                rel_new = np.linalg.inv(self.prev_pose) @ T_WC
                if not (self.jump_gate > 0 and float(np.linalg.norm(rel_new[:3, 3])) > self.jump_gate):
                    self.prev_rel = rel_new
        self.prev_pose = T_WC
        self.pose_estimates[f] = T_WC  # geçici: ham init (refresh'te optimize ile değiştirilir)
        self._last_obs = obs           # loop_add için bu karenin gözlemleri
        self.f += 1
        if update_ok:                  # yalnız gerçekten grafa giren kare CV-factor hedefi olabilir
            self.prev_graph_f = f
            self.prev_graph_pose = T_WC

        if self.use_fl and update_ok:
            # Fixed-lag: pencere küçük, her karede calculateEstimate ucuz. Marjinalize olanların son değeri
            # pose_estimates/lm_world'de kalır; live_lms yeniden kurulur.
            try:
                est_fl = self.smoother.calculateEstimate()
                self.live_lms = set()
                self.resident_lms = set()
                for key in est_fl.keys():
                    s = gtsam.Symbol(key); c = s.chr()
                    try:
                        if c == ord('x'):
                            self.pose_estimates[s.index()] = est_fl.atPose3(key).matrix()
                        elif c == ord('l'):
                            lidx = s.index()
                            self.resident_lms.add(lidx)
                            if self.culler.active and self.culler.is_bad(lidx):
                                continue   # kara liste canlı sayılmaz
                            p = est_fl.atPoint3(key)
                            self.live_lms.add(lidx)
                            self.lm_world[lidx] = np.array([p[0], p[1], p[2]])
                    except Exception:
                        pass
                if self.culler.active:
                    self.culler.on_window(self.resident_lms)
            except Exception as e:
                print(f"  [kare {f}] fixed-lag refresh: {e}")
        # Tam-graf: PnP önbelleğini ve pozları periyodik olarak optimize çözümden tazele; yoksa PnP
        # ham üçgenlemeyle sapar ve çizim titrer. FL modunda ya da update başarısızsa buraya girilmez.
        elif (not self.use_fl) and update_ok and (self.f % self.refresh_every == 0 or self.f == 1):
            try:
                est = self.isam.calculateEstimate()
                for lid in self.initialized:
                    if self.culler.active and self.culler.is_bad(lid):
                        continue
                    if est.exists(L(lid)):
                        p = est.atPoint3(L(lid))
                        self.lm_world[lid] = np.array([p[0], p[1], p[2]])
                for k in list(self.pose_estimates.keys()):
                    if est.exists(X(k)):
                        self.pose_estimates[k] = est.atPose3(X(k)).matrix()
            except Exception as e:
                print(f"  [kare {f}] refresh: {e}")

    def display_trajectory(self, gt_t, gt_pos, gt_quat):
        """Kestirilen pozları GT'ye ankorlayıp gövde-frame pozisyon listesi döndür (canlı çizim için).
        Anchor henüz yoksa boş liste."""
        if not self.pose_estimates:
            return []
        fts = np.asarray(self.frame_ts, dtype=np.int64)
        covered = np.where(fts >= gt_t[0])[0]
        if len(covered) == 0:
            return []
        f_a = int(covered[0])
        if f_a not in self.pose_estimates:
            ks = [k for k in self.pose_estimates if k >= f_a]
            if not ks:
                return []
            f_a = min(ks)
        T_W_W0 = gt_pose_at(int(fts[f_a]), gt_t, gt_pos, gt_quat) @ self.T_B_C @ \
            np.linalg.inv(self.pose_estimates[f_a])
        return [((T_W_W0 @ self.pose_estimates[k]) @ self.T_C_B)[:3, 3]
                for k in sorted(self.pose_estimates)]

    def finalize(self, gt_path, out_csv):
        # Fixed-lag'de pozlar her karede toplandı; toplu okuma yalnız tam-graf modunda.
        if not self.use_fl:
            est = self.isam.calculateEstimate()
            for f in range(self.f):
                if est.exists(X(f)):
                    self.pose_estimates[f] = est.atPose3(X(f)).matrix()

        # Loop PGO düzeltmesi: keyframe pozlarını düzelt, düzeltmeyi tüm karelere yay. GT ankorundan önce çalışır.
        if self.loop_closer is not None and self.loop_closer.n_detected > 0:
            kf_frames, corrected = self.loop_closer.finalize()
            kf_sorted = sorted(kf_frames)
            old = {fk: self.pose_estimates[fk].copy() for fk in kf_sorted if fk in self.pose_estimates}
            deltas = {fk: corrected[fk] @ np.linalg.inv(old[fk]) for fk in kf_sorted if fk in old}
            if deltas:
                dk = sorted(deltas.keys()); dk_arr = np.array(dk)
                for f in list(self.pose_estimates.keys()):
                    pos = np.searchsorted(dk_arr, f, side='right') - 1
                    fk = dk[pos] if pos >= 0 else dk[0]
                    self.pose_estimates[f] = deltas[fk] @ self.pose_estimates[f]
            print(f"[LOOP] Düzeltme uygulandı ({len(deltas)} keyframe, tüm karelere yayıldı).")
        elif self.loop_closer is not None:
            self.loop_closer.finalize()  # özet yazdır (loop yok)

        gt_t, gt_pos, gt_quat = load_groundtruth(gt_path)
        frame_ts = np.array(self.frame_ts, dtype=np.int64)
        covered = np.where(frame_ts >= gt_t[0])[0]
        if len(covered) == 0:
            print("UYARI: hiçbir kare GT aralığında değil, ATE atlandı."); return None
        f_a = int(covered[0])
        if f_a not in self.pose_estimates:
            f_a = min(self.pose_estimates.keys(), key=lambda k: abs(k - f_a))
        T_W_W0 = gt_pose_at(frame_ts[f_a], gt_t, gt_pos, gt_quat) @ self.T_B_C @ \
            np.linalg.inv(self.pose_estimates[f_a])

        rows = []
        for f in range(self.f):
            if f not in self.pose_estimates:
                continue
            T_WB = (T_W_W0 @ self.pose_estimates[f]) @ self.T_C_B
            p = T_WB[:3, 3]; q = Rot.from_matrix(T_WB[:3, :3]).as_quat()
            rows.append([frame_ts[f], p[0], p[1], p[2], q[0], q[1], q[2], q[3]])
        import pandas as pd
        from pathlib import Path
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows, columns=['#timestamp [ns]', ' p_RS_R_x [m]', ' p_RS_R_y [m]', ' p_RS_R_z [m]',
                                    ' q_RS_x []', ' q_RS_y []', ' q_RS_z []', ' q_RS_w []']).to_csv(out_csv, index=False)
        et = np.array([r[0] for r in rows], float); ex = np.array([[r[1], r[2], r[3]] for r in rows])
        ts, te = max(gt_t[0], et[0]), min(gt_t[-1], et[-1]); m = (gt_t >= ts) & (gt_t <= te)
        ip = np.stack([np.interp(gt_t[m], et, ex[:, i]) for i in range(3)], 1)
        ate = float(np.sqrt(np.mean(np.sum((gt_pos[m] - ip) ** 2, 1))))
        print(f"Trajectory yazıldı: {out_csv} ({len(rows)} poz)")
        return ate


REALTIME_SUMMARY_FIELDS = ['sequence', 'ate_rmse_m', 'rmse_x', 'rmse_y', 'rmse_z',
                          'compare_points', 'n_frames', 'wall_s', 'fps',
                          'fe_ms_per_frame', 'fe_fps', 'be_ms_per_frame', 'be_fps',
                          'extractor', 'matcher']


def _write_realtime_summary(cfg, row):
    """outputs/realtime_summary.csv'ye bu sekansın satırını yaz; sekans başına tek satır tutulur, varsa üzerine yazılır."""
    path = Path(cfg['paths']['outputs']) / 'realtime_summary.csv'
    rows = {}
    if path.exists():
        with open(path, newline='') as fh:
            for r in csv.DictReader(fh):
                rows[r['sequence']] = r
    rows[row['sequence']] = row
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=REALTIME_SUMMARY_FIELDS)
        w.writeheader()
        w.writerows(rows.values())
    print(f"Realtime özet: {path}")


def main():
    ap = argparse.ArgumentParser(description="Komple realtime VIO + FPS (config.yaml temelli)")
    ap.add_argument('--max_frames', type=int, default=None)
    # Görselleştirme bayrakları; varsayılanları config.yaml belirler. Yörüngeyi etkilemez.
    ap.add_argument('--viz_interval', type=int, default=None,
                    help='her N karede bir görselleştir (config.yaml frontend.viz_interval)')
    ap.add_argument('--no_save_viz', action='store_true',
                    help='PNG kaydetmeyi kapat (canlı pencere açık kalır)')
    # Arka ucu ayrı süreçte koşturur; ATE değişmez (gerekçe: realtime_pipeline.py). Varsayılan kapalı.
    ap.add_argument('--pipeline', action='store_true',
                    help='arka ucu ayrı süreçte koştur (ön uçla örtüşür; ATE aynı kalır)')
    cli, _ = ap.parse_known_args()

    cfg = config.load()
    ds, fe = cfg['dataset'], cfg['frontend']
    if cli.viz_interval is not None:
        fe['viz_interval'] = int(cli.viz_interval)
    if cli.no_save_viz:
        fe['save_visualizations'] = False
    bargs = config.backend_args(cfg, argv=[])
    max_frames = cli.max_frames if cli.max_frames is not None else ds['max_frames']

    print("Front-end kuruluyor (canlı işleme)...")
    matcher = StereoMatcher(
        base_path=ds['base_path'], sequence=ds['sequence'], matcher_type=fe['matcher'],
        filter_method=fe['filter'], y_axis_threshold=fe['y_axis_threshold'],
        extractor=fe['extractor'], max_track_length=fe['max_track_length'],
        max_keypoints=fe['max_keypoints'], detection_threshold=fe['detection_threshold'])
    config.apply_frontend(matcher, cfg)

    # Boru hattında bu nesne ebeveynde yalnız yapılandırmayı okumak için; gerçek kestirim
    # çocuk süreçte aynı bargs ile kurulur.
    backend = RealtimeBackend(matcher.get_calibration_export(), bargs)
    pipe = {'q_in': None, 'q_out': None, 'proc': None} if cli.pipeline else None
    backend_has_loop = backend.loop_closer is not None

    # Loop closure descriptor ister; export_keypoints otomatik açılır.
    if backend.loop_closer is not None and not matcher.EXPORT_KEYPOINTS:
        matcher.EXPORT_KEYPOINTS = True
        print("  [LOOP] frontend.export_keypoints otomatik AÇILDI (loop closure descriptor gerektirir).")

    # Ground-truth (Pangolin kırmızı referans + anchor)
    gt_t, gt_pos, gt_quat = load_groundtruth(config.gt_path(cfg))
    gt_draw = gt_pos[::5]  # çizim için seyreltilmiş GT (hız)

    # Pangolin yörünge penceresi kurulumu
    seq = ds['sequence']
    pview = {'ok': False, 'cv': None, 'traj': []}
    if _PANGO:
        try:
            w, h = cfg['live']['window_w'], cfg['live']['window_h']
            cen = gt_pos.mean(0); ext = float(np.linalg.norm(gt_pos.max(0) - gt_pos.min(0))) + 1e-6
            title = (f"REALTIME — {seq}  |  extractor={fe['extractor']} matcher={fe['matcher']}"
                     f"  |  GT(kirmizi) vs kestirim(mavi)")
            pango.CreateWindowAndBind(title, w, h)
            glEnable(GL_DEPTH_TEST)
            scam = pango.OpenGlRenderState(
                pango.ProjectionMatrix(w, h, 600, 600, w // 2, h // 2, 0.1, 5000),
                pango.ModelViewLookAt(cen[0], cen[1] - ext, cen[2] + ext * 0.6,
                                      cen[0], cen[1], cen[2], 0.0, 0.0, 1.0))
            handler = pango.Handler3D(scam)
            dcam = pango.CreateDisplay()
            dcam.SetBounds(pango.Attach(0.0), pango.Attach(1.0), pango.Attach(0.0), pango.Attach(1.0), -w / h)
            dcam.SetHandler(handler)
            pview['ok'] = True
            print("Pangolin yörünge penceresi açık (Live Tracking ile beraber).")
        except Exception as e:
            print(f"Pangolin penceresi açılamadı: {e}")
    # Pangolin yoksa ya da açılamadıysa OpenCV yedek pencere.
    if not pview['ok']:
        try:
            pview['cv'] = Cv2TrajView(gt_pos, cfg['live']['window_w'], cfg['live']['window_h'])
            print("Pangolin yok → YEDEK yörünge penceresi (OpenCV, tepeden X-Y) açıldı.")
            print("  Pangolin istiyorsanız: pypangolin, SLAM'in koştuğu Python sürümü için "
                  "derlenmiş olmalı (sürüme-özgü uzantı).")
        except Exception as e:
            print(f"Yedek yörünge penceresi de açılamadı (görselsiz devam): {e}")

    refresh_every = max(1, int(cfg['live']['refresh_every']))

    def render_pango():
        glClearColor(1, 1, 1, 1)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        dcam.Activate(scam)
        _draw_poly(gt_draw, (0.85, 0.1, 0.1), 2.5)          # GT (kırmızı)
        _draw_poly(pview['traj'], (0.1, 0.3, 0.9), 2.0)     # kestirim (mavi)
        if pview['traj']:
            _draw_pt(pview['traj'][0], (0.1, 0.8, 0.1), 14)  # başlangıç (yeşil)
            _draw_pt(pview['traj'][-1], (1.0, 0.85, 0.0), 16)  # anlık uç (sarı)
        pango.FinishFrame()

    # Zamanlama: ön uç ve arka uç ayrı sayılır.
    timing = {'fe': 0.0, 'be': 0.0, 'n': 0, 'last': None}

    def on_frame(timestamp, obs, frame_idx):
        now = time.time()
        if timing['last'] is not None:
            timing['fe'] += now - timing['last']      # bu kareye kadarki front-end süresi
        t0 = time.time()
        if pipe is not None:
            # Boru hattı: kareyi arka uç sürecine yolla ve hemen dön; kuyruk sınırlı, arka uç geride kalırsa burada bloklanır.
            want = (timing['n'] % refresh_every == 0)
            _kps = _desc = None
            if backend_has_loop:
                _kps = getattr(matcher, '_cur_all_kps', None)
                _desc = getattr(matcher, '_cur_all_descs', None)
            pipe['q_in'].put(('f', timestamp, obs, _kps, _desc, want))
            timing['be'] += time.time() - t0          # yalnız kuyruğa koyma maliyeti
            # Arka uçtan gelen yörünge anlık görüntülerini topla (bloklamadan)
            while True:
                try:
                    kind, *rest = pipe['q_out'].get_nowait()
                except Exception:                     # noqa: BLE001 — kuyruk boş
                    break
                if kind == 'traj':
                    pview['traj'] = rest[0]
                elif kind == 'err':
                    raise RuntimeError(f'arka uç süreci düştü: {rest[0]}')
        else:
            backend.step(timestamp, obs)              # back-end iSAM2 adımı
            timing['be'] += time.time() - t0
            if backend.loop_closer is not None:       # bu karenin keypoint/desc'lerini loop DB'sine ver
                backend.loop_add(getattr(matcher, '_cur_all_kps', None),
                                 getattr(matcher, '_cur_all_descs', None))
        timing['n'] += 1
        # Yörüngeyi periyodik tazele, her kare çiz. Boru hattında yörünge q_out ile çocuktan gelir.
        if pipe is None and (pview['ok'] or pview['cv'] is not None):
            if timing['n'] % refresh_every == 0:
                pview['traj'] = backend.display_trajectory(gt_t, gt_pos, gt_quat)
        if pview['cv'] is not None:
            pview['cv'].render(pview['traj'])
        if pview['ok']:
            render_pango()
            # Pangolin penceresi kapatılırsa (X / ESC) KeyboardInterrupt ile temiz durdur; finalize ve rapor çalışır.
            if pango.ShouldQuit():
                print("Pangolin penceresi kapatıldı → durduruluyor...")
                raise KeyboardInterrupt
        timing['last'] = time.time()

    out_csv = config.traj_path(cfg, '_realtime')

    # Arka uç sürecini başlat; spawn ile temiz yorumlayıcı, CUDA bağlamı miras alınmaz.
    if pipe is not None:
        import multiprocessing as _mp
        import realtime_pipeline as _rp
        _ctx = _mp.get_context('spawn')
        pipe['q_in'] = _ctx.Queue(maxsize=_rp.QUEUE_MAX)
        pipe['q_out'] = _ctx.Queue()
        pipe['proc'] = _ctx.Process(
            target=_rp.backend_worker,
            args=(matcher.get_calibration_export(), bargs, (gt_t, gt_pos, gt_quat),
                  config.gt_path(cfg), out_csv, pipe['q_in'], pipe['q_out']),
            daemon=True)
        pipe['proc'].start()
        print(f"  [BORU HATTI] arka uç ayrı süreçte (pid={pipe['proc'].pid}), "
              f"kuyruk={_rp.QUEUE_MAX} kare.")

    print(f"REALTIME başladı (sekans={ds['sequence']}, extractor={fe['extractor']}, "
          f"matcher={fe['matcher']}). 'q' ile çıkış.")
    timing['last'] = time.time()
    t_start = time.time()
    try:
        matcher.process_sequence(max_frames=max_frames, visualize=True,
                                 viz_interval=fe['viz_interval'],
                                 save_visualizations=fe['save_visualizations'],
                                 on_frame=on_frame)
    finally:
        # Ön uç bitti (ya da Ctrl+C): kuyruğu kapat, arka ucun birikmiş kareleri tüketip finalize etmesini bekle.
        if pipe is not None and pipe.get('proc') is not None:
            pipe['q_in'].put(None)

    ate = None
    if pipe is not None:
        t_drain = time.time()
        while True:
            kind, *rest = pipe['q_out'].get()          # 'done' veya 'err' gelene dek bekle
            if kind == 'traj':
                pview['traj'] = rest[0]
                continue
            if kind == 'err':
                pipe['proc'].join(timeout=10)
                raise RuntimeError(f'arka uç süreci düştü: {rest[0]}')
            ate, be_secs, n_be = rest                  # 'done'
            timing['be'] = be_secs                     # çocukta ölçülen gerçek arka uç süresi
            break
        pipe['proc'].join(timeout=60)
        drain = time.time() - t_drain
        print(f"  [BORU HATTI] arka uç {n_be} kare işledi; ön uç bittikten sonra "
              f"{drain:.1f} s kuyruk boşalttı.")
    else:
        ate = backend.finalize(config.gt_path(cfg), out_csv)
    wall = time.time() - t_start

    # Final görsel: tüm düzeltmeler (PGO dahil) uygulanmış son trajectory.
    axis_rmse, n_cmp = None, None
    try:
        from eval_trajectory import evaluate as _plot
        png = config.plot_path(cfg, '_realtime')
        _, axis_rmse, n_cmp = _plot(out_csv, config.gt_path(cfg), png,
                                    dataset_name=f"{seq} (realtime, son hâl)")
        print(f"Final görsel: {png}")
    except Exception as e:
        print(f"[PLOT] görsel atlandı: {e}")

    n = max(timing['n'], 1)
    fe_fps = n / max(timing['fe'], 1e-9)
    be_fps = n / max(timing['be'], 1e-9)
    print("\n" + "=" * 60)
    print("  REALTIME HIZ RAPORU (akademik değerlendirme)")
    print("=" * 60)
    print(f"  İşlenen kare           : {timing['n']}")
    print(f"  Toplam süre (wall)     : {wall:.2f} s")
    print(f"  ORTALAMA HIZ           : {timing['n'] / wall:.2f} FPS")
    print(f"  Front-end (canlı+viz)  : {1000 * timing['fe'] / n:.1f} ms/kare  ({fe_fps:.1f} FPS)")
    print(f"  Back-end (iSAM2 BA)    : {1000 * timing['be'] / n:.1f} ms/kare  ({be_fps:.1f} FPS)")
    if pipe is not None:
        # Boru hattında ön uç ve arka uç örtüşür; duvar süresi toplam değil yaklaşık büyüğüdür.
        seri = timing['fe'] + timing['be']
        print(f"  Boru hattı örtüşmesi   : seri {seri:.1f} s → gerçek {wall:.1f} s "
              f"(kazanç {seri / max(wall, 1e-9):.2f}x)")
    if ate is not None:
        print(f"  ATE (RMSE)             : {ate:.4f} m")
    print("=" * 60)

    _write_realtime_summary(cfg, {
        'sequence': seq,
        'ate_rmse_m': round(ate, 4) if ate is not None else '',
        'rmse_x': round(float(axis_rmse[0]), 4) if axis_rmse is not None else '',
        'rmse_y': round(float(axis_rmse[1]), 4) if axis_rmse is not None else '',
        'rmse_z': round(float(axis_rmse[2]), 4) if axis_rmse is not None else '',
        'compare_points': n_cmp if n_cmp is not None else '',
        'n_frames': timing['n'],
        'wall_s': round(wall, 2),
        'fps': round(timing['n'] / wall, 2) if wall > 0 else '',
        'fe_ms_per_frame': round(1000 * timing['fe'] / n, 2),
        'fe_fps': round(fe_fps, 2),
        'be_ms_per_frame': round(1000 * timing['be'] / n, 2),
        'be_fps': round(be_fps, 2),
        'extractor': fe['extractor'],
        'matcher': fe['matcher'],
    })


if __name__ == '__main__':
    main()
