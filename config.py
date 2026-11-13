"""
Tek kontrol noktası yükleyici: config.yaml'ı okur, frontend/backend/canlı argparse
varsayılanlarını oradan besler. Öncelik: CLI override > config.yaml > kod.

    cfg = config.load()             # dict
    args = config.backend_args()    # argparse.Namespace (config + CLI)
"""
import argparse
import os
import numpy as np
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parent

#: [TÜRBÜLANS] Alternatif yapılandırma dosyası ortam değişkeniyle seçilir (ör. SLAM_CONFIG=config_turb.yaml).
#: Tanımsızsa ROOT/config.yaml (PyEco davranışı birebir). Göreli yol ROOT'a göre çözülür.
CONFIG_ENV = "SLAM_CONFIG"


def _resolve_default_path():
    p = os.environ.get(CONFIG_ENV)
    if not p:
        return ROOT / "config.yaml"
    p = Path(p).expanduser()
    return p if p.is_absolute() else ROOT / p


DEFAULT_PATH = _resolve_default_path()

#: Veri kökü: $SLAM_DATA_ROOT > config.yaml:dataset.base_path.
#: config.yaml mutlak yol tuttuğundan makine değişince yollar kırılıyordu; config dosyasına dokunmadan taşımak için.
DATA_ROOT_ENV = "SLAM_DATA_ROOT"


def load(path=None):
    """config.yaml'ı dict olarak yükle; veri kökünü ortama göre çöz + sekansı doğrula."""
    with open(path or DEFAULT_PATH, "r") as f:
        cfg = yaml.safe_load(f)
    env_root = os.environ.get(DATA_ROOT_ENV)
    if env_root:
        cfg.setdefault("dataset", {})["base_path"] = env_root
    _check_sequence(cfg)
    _scale_with_budget(cfg)
    return cfg


def _scale_with_budget(cfg):
    """Budget-proportional settings: T, N_replenish and the far-stereo minimum follow max_keypoints (kp)."""
    f, b = cfg.get("frontend", {}), cfg.get("backend", {})
    kp = int(f["max_keypoints"])
    if f.get("target_ratio_kp") is not None:
        f["target_feature_count"] = round(float(f["target_ratio_kp"]) * kp)
    if f.get("replenish_ratio_kp") is not None:
        f["replenish_threshold"] = round(float(f["replenish_ratio_kp"]) * kp)
    if b.get("stereo_uzak_min_ratio_kp") is not None:
        b["stereo_uzak_min"] = round(float(b.pop("stereo_uzak_min_ratio_kp")) * kp)
    # Budget-adaptive window: per-frame window (butce_kare_ayar) at kp <= butce_kare_esik, keyframe window above
    if kp <= int(b.get("butce_kare_esik", 0) or 0):
        b.update(b.get("butce_kare_ayar") or {})


def _check_sequence(cfg):
    """`dataset.sequence` diskte var mı? Yoksa mevcut sekansları listeleyerek hata ver.
    Sekans adı klasör adıdır; kısaltma yazılınca hata çok sonra alakasız bir yerde çıkıyordu."""
    ds = cfg.get("dataset") or {}
    base, seq = ds.get("base_path"), ds.get("sequence")
    if not base or not seq:
        return
    root = Path(base)
    if (root / seq / "mav0").is_dir():
        return
    if not root.is_dir():
        raise SystemExit(
            f"\ndataset.base_path bulunamadı: {root}\n"
            f"  config.yaml'da düzeltin ya da $SLAM_DATA_ROOT ile geçersiz kılın.\n")
    have = sorted(d.name for d in root.iterdir()
                  if (d / "mav0" / "state_groundtruth_estimate0" / "data.csv").is_file())
    near = [h for h in have if h.lower().startswith(str(seq).lower()[:5])]
    raise SystemExit(
        f"\ndataset.sequence = '{seq}' → '{root / seq / 'mav0'}' yok.\n"
        + (f"  Bunu mu demek istediniz: {', '.join(near)}\n" if near else "")
        + f"  Kullanılabilir sekanslar: {', '.join(have) if have else '(hiçbiri)'}\n"
        + "  (Sekans adı KLASÖR adıdır; kısaltma değil tam adı yazın.)\n")


# ----------------------------- Yol yardımcıları -----------------------------
def gt_path(cfg):
    d = cfg["dataset"]
    return str(Path(d["base_path"]) / d["sequence"] /
               "mav0/state_groundtruth_estimate0/data.csv")


def dump_path(cfg):
    return str(Path(cfg["paths"]["outputs"]) / f"landmarks_{cfg['dataset']['sequence']}.pkl")


def traj_path(cfg, tag=""):
    return str(Path(cfg["paths"]["outputs"]) / f"gtsam_traj_{cfg['dataset']['sequence']}{tag}.csv")


def plot_path(cfg, tag=""):
    """Trajectory görseli (X-Y / X-Z / 3B), sekans adına göre PNG."""
    return str(Path(cfg["paths"]["outputs"]) / f"trajectory_{cfg['dataset']['sequence']}{tag}.png")


def _culling_section(cfg, legacy):
    """Birleşik `culling:` bölümü; anahtar orada yoksa eski konuma (frontend/backend) düş."""
    cu = dict(legacy or {})
    cu.update(cfg.get("culling") or {})
    return cu


# ----------------------------- Frontend -------------------------------------
def apply_frontend(matcher, cfg):
    """config.yaml frontend bölümünü StereoMatcher örneğine uygula (yapım sonrası nitelikler)."""
    f = cfg["frontend"]
    matcher.REPLENISH_THRESHOLD = int(f["replenish_threshold"])
    matcher.TARGET_FEATURE_COUNT = int(f["target_feature_count"])
    # Track ömrü + dağılım gridi (tavlamada aranır).
    matcher.MAX_TRACK_AGE_UNSEEN = int(f.get("max_track_age_unseen", 5))
    if str(f.get("grid_rows", 10)) != "auto" and str(f.get("grid_cols", 12)) != "auto":
        matcher.GRID_ROWS = int(f.get("grid_rows", 10))
        matcher.GRID_COLS = int(f.get("grid_cols", 12))
    matcher.MIN_STEREO_DISPARITY = float(f["min_stereo_disparity"])
    # Derinlik penceresi sensorden: [min_depth_baselines*b, th_depth_baselines*b] (ORB-SLAM3 Frame.cc:841
    # minZ=b; ORB-SLAM2 makalesi III-A ve ORB-SLAM3 Examples/Stereo/TUM-VI.yaml:61 ThDepth=40). Acik bir
    # max_landmark_depth verilirse o kullanilir.
    _b = float(np.linalg.norm(matcher.T))
    if f.get("max_landmark_depth") is None:
        matcher.MAX_LANDMARK_DEPTH = float(f["th_depth_baselines"]) * _b
    else:
        # .inf: ORB-SLAM3 Frame::ComputeStereoMatches gibi ust derinlik siniri yok (disparite > 0, minZ = b)
        matcher.MAX_LANDMARK_DEPTH = float(f["max_landmark_depth"])
    matcher.TRIGGER_CLOSE_DEPTH = (float(f["th_depth_baselines"]) * _b) if f.get("trigger_close_only") else 0.0
    matcher.MIN_LANDMARK_DEPTH = float(f.get("min_depth_baselines", 0.0)) * _b
    matcher.TEMPORAL_RANSAC_THRESHOLD = float(f["temporal_ransac_threshold"])
    matcher.TEMPORAL_RANSAC_CONF = float(f.get("temporal_ransac_confidence", 0.99))
    matcher.TEMPORAL_RANSAC_ITERS = int(f.get("temporal_ransac_iters", 2000))
    matcher.TEMPORAL_MATCHER = str(f.get("temporal_matcher", "bf_lowe"))
    matcher.TEMPORAL_LG_KAPILI = bool(f.get("temporal_lg_kapili", False))
    matcher.TEMPORAL_LG_SAGLIK = float(f.get("temporal_lg_saglik", 0.0))
    matcher.EXPORT_S_D = bool(f.get("export_s_d", False))
    matcher.TEMPORAL_LOWE_RATIO = float(f.get("temporal_lowe_ratio", 0.8))
    matcher.STEREO_DEPTH_ASSOC_PX = float(f.get("stereo_depth_assoc_px", 1.5))
    matcher.STEREO_TRACK_RATE = float(f.get("stereo_track_rate", 0.0))
    matcher.LG_IMAGE_SIZE_WH = bool(f.get("lg_image_size_wh", False))
    matcher.BOS_KARE_DUZELT = bool(f.get("bos_kare_duzelt", False))
    # Izgara: 'auto' -> hucre sayisi ~ hedef track sayisi T, en-boy orani goruntuyle ayni (ORB-SLAM
    # ORBextractor.cc DistributeOctTree: dugum sayisi istenen nokta sayisina ulasana kadar bolunur).
    if str(f.get("grid_rows")) == "auto" or str(f.get("grid_cols")) == "auto":
        _h, _w = matcher.map1x.shape[:2]
        _T = max(1, int(f["target_feature_count"]))
        _r = max(1, int(round((_T * _h / _w) ** 0.5)))
        matcher.GRID_ROWS, matcher.GRID_COLS = _r, max(1, int(round(_T / _r)))
        print(f"Izgara (auto): {matcher.GRID_ROWS}x{matcher.GRID_COLS} (T={_T}, {_w}x{_h})")
    print(f"derinlik penceresi [{matcher.MIN_LANDMARK_DEPTH:.3f}, {matcher.MAX_LANDMARK_DEPTH:.3f}] m, "
          f"tetik_yakin<{matcher.TRIGGER_CLOSE_DEPTH:.2f} m, "
          f"y_tol={matcher.y_axis_threshold}, zamansal={matcher.TEMPORAL_MATCHER}, "
          f"stereo_orani={matcher.STEREO_TRACK_RATE}, izgara={matcher.GRID_ROWS}x{matcher.GRID_COLS}")
    # Ölçek-duyarlı tetik bayrakları (ablation)
    matcher.SCALE_AWARE_TRIGGER = bool(f["scale_aware_trigger"])
    matcher.SCALE_NEAR_TRIGGER = bool(f["scale_near_trigger"])
    matcher.SCALE_NEAR_DEPTH = float(f["scale_near_depth"])
    matcher.SCALE_MIN_NEAR_ANCHORS = int(f["scale_min_near_anchors"])
    matcher.SCALE_STALE_TRIGGER = bool(f["scale_stale_trigger"])
    matcher.SCALE_MAX_DEPTH_AGE = int(f["scale_max_depth_age"])
    matcher.SCALE_MAX_STALE_FRAC = float(f["scale_max_stale_frac"])
    # Paralaks-açısı üçgenleme kapısı
    matcher.PARALLAX_CHECK = bool(f.get("parallax_check", False))
    matcher.COS_MAX_PARALLAX = float(f.get("cos_max_parallax", 0.9998))
    # Kare başına keypoint+descriptor dışa aktarımı (loop closure bunu ister)
    matcher.EXPORT_KEYPOINTS = bool(f.get("export_keypoints", False))
    # Adaptif özellik çıkarımı: entropi güdümlü CLAHE + dedektör eşiği uyarlaması (varsayılan kapalı).
    matcher.ADAPTIVE_PREPROCESSING_ENABLED = bool(f.get("adaptive_preprocessing", False))
    matcher.ADAPTIVE_CLAHE_THRESHOLD_S_D = float(f.get("adaptive_clahe_s_d_threshold", 0.40))
    matcher.ADAPTIVE_CLAHE_CLIP_LIMIT_MIN = float(f.get("adaptive_clahe_clip_min", 0.5))
    matcher.ADAPTIVE_CLAHE_CLIP_LIMIT_MAX = float(f.get("adaptive_clahe_clip_max", 2.5))
    matcher.ADAPTIVE_CLAHE_TILE = int(f.get("adaptive_clahe_tile", 8))
    # Bütçe-eşlemeli CLAHE (BE-CLAHE): varsayılan kapalı; kapalıyken davranış bayt-özdeş.
    matcher.BE_CLAHE = bool(f.get("adaptive_clahe_butce_eslesik", False))
    matcher.BE_MERDIVEN = tuple(float(x) for x in f.get("adaptive_clahe_be_merdiven", (0.5, 1.0, 1.5, 2.0, 2.5)))
    matcher.BE_HEDEF = int(f["max_keypoints"])
    matcher.BE_ARAMA = str(f.get("adaptive_clahe_be_arama", "tam"))
    matcher.BE_PROFIL = bool(f.get("adaptive_clahe_be_profil", False))
    matcher.BE_ONBELLEK = bool(f.get("adaptive_clahe_be_onbellek", False))  # seçilen basamağın haritasını yeniden kullan
    if matcher.BE_CLAHE:
        print(f"BE-CLAHE ACIK: C={matcher.BE_MERDIVEN} hedef n0>={matcher.BE_HEDEF} arama={matcher.BE_ARAMA} "
              f"profil={matcher.BE_PROFIL} onbellek={matcher.BE_ONBELLEK} tau={matcher.ADAPTIVE_CLAHE_THRESHOLD_S_D}")
    matcher.ADAPTIVE_DET_ENABLED = bool(f.get("adaptive_det_threshold", True))
    matcher.ADAPTIVE_DET_MIN_THRESHOLD = float(f.get("adaptive_det_threshold_min", 1e-5))
    # Dedektör kolunu tau kapısına bağla (varsayılan kapalı).
    matcher.ADAPTIVE_DET_GATED = bool(f.get("adaptive_det_gated", False))
    # Taban-oranlı eşik: theta = clamp(s_(T), kappa*s_ort, theta_max), karenin kendi skor dağılımından.
    # Diğer kollarla dışlayıcı; bkz. TabanOranliCikarici.
    if bool(f.get("adaptive_det_floor", False)):
        from stereo_front_end import TabanOranliCikarici
        if not isinstance(matcher.aliked, TabanOranliCikarici):
            _tmax = float(f.get("adaptive_det_threshold_max") or f["detection_threshold"])
            matcher.aliked = TabanOranliCikarici(
                matcher.aliked, int(f["target_feature_count"]), _tmax,
                kappa=float(f.get("adaptive_det_floor_kappa", 2.0)),
                gamma=float(f.get("adaptive_det_floor_gamma", 1.0)))
        matcher.ADAPTIVE_DET_ENABLED = False
        print(f"Taban-oranli esik ACIK: hedef={f['target_feature_count']} "
              f"kappa={f.get('adaptive_det_floor_kappa', 2.0)} theta_max={_tmax:g}")

    # Arz açığı kapılı kol: çıkarıcı ArzKapiliCikarici ile sarılır, S_d güdümlü kol kapatılır.
    if bool(f.get("adaptive_det_deficit", False)):
        from stereo_front_end import ArzKapiliCikarici
        _ref = str(f.get("adaptive_det_deficit_ref", "target")).lower()
        _n = (int(f["target_feature_count"]) if _ref == "target"
              else int(f["max_keypoints"]))
        _tmax = float(f.get("adaptive_det_threshold_max") or f["detection_threshold"])
        if not isinstance(matcher.aliked, ArzKapiliCikarici):
            _pay = float(f.get("adaptive_det_deficit_margin", 1.0))
            matcher.aliked = ArzKapiliCikarici(
                matcher.aliked, _n, _tmax, matcher.ADAPTIVE_DET_MIN_THRESHOLD, _pay,
                us=float(f.get("adaptive_det_deficit_us", 0.0)))
        matcher.ADAPTIVE_DET_ENABLED = False
        print(f"Arz açığı kapılı kol AÇIK: referans={_ref} N={_n} x pay={_pay:g} "
              f"= {matcher.aliked.referans_n}, "
              f"theta {_tmax:g} -> {matcher.ADAPTIVE_DET_MIN_THRESHOLD:g}")
    # Eşik tavanı verilmezse kurucudaki detection_threshold kalır.
    if f.get("adaptive_det_threshold_max", None) is not None:
        matcher.ADAPTIVE_DET_MAX_THRESHOLD = float(f["adaptive_det_threshold_max"])
    # CLAHE tile değişmiş olabilir → nesne cache'ini sıfırla.
    matcher._clahe_cache = {}
    # Ön-uç track ayıklama; anahtarlar birleşik `culling:` bölümünden, hepsi varsayılan kapalı.
    cu = _culling_section(cfg, f)
    matcher.CULL_RANSAC_STRIKES = int(cu.get("cull_ransac_strikes", 0))
    matcher.CULL_DESC_DRIFT_MAX = float(cu.get("cull_desc_drift_max", 0.0))
    matcher.CULL_FLICKER_RATIO = float(cu.get("cull_flicker_ratio", 0.0))
    matcher.CULL_FLICKER_MIN_LIFE = int(cu.get("cull_flicker_min_life", 10))
    matcher.STEREO_REPROJ_CHECK = bool(cu.get("stereo_reproj_check", False))
    matcher.STEREO_REPROJ_MAX_PX = float(cu.get("stereo_reproj_max_px", 5.0))
    matcher.MAX_DEPTH_SIGMA_REL = float(cu.get("max_depth_sigma_rel", 0.0))
    matcher.STEREO_DISP_SIGMA_PX = float(cu.get("stereo_disp_sigma_px", 1.0))
    # [TÜRBÜLANS] `turbulence:` bölümü. enabled=false (ve restoration=none) iken turb paketi hiç içe
    # aktarılmaz, matcher.turb = None kalır ve ön-uç PyEco ile bit düzeyinde aynıdır (E0).
    _t = cfg.get("turbulence") or {}
    # Faz 2 (2026-11-08): sensor_noise_dn > 0 da kancayı açar (D/r0=0'da tohumlu temiz dağılım); 0 -> eski kapı.
    if (bool(_t.get("enabled", False)) or str(_t.get("restoration", "none")) != "none"
            or float(_t.get("sensor_noise_dn", 0.0) or 0.0) > 0):
        from turb.hooks import injector_from_cfg
        matcher.turb = injector_from_cfg(cfg)


# ----------------------------- Backend / Canlı ------------------------------
def _add_backend_args(ap, cfg):
    """Backend ve canlı görselleştirmenin paylaştığı argümanlar (config-temelli)."""
    b, ds = cfg["backend"], cfg["dataset"]
    ap.add_argument("--dump", default=dump_path(cfg))
    ap.add_argument("--gt", default=gt_path(cfg))
    ap.add_argument("--start_frame", type=int, default=int(ds.get("start_frame", 0) or 0))
    ap.add_argument("--max_frames", type=int, default=ds["max_frames"])
    ap.add_argument("--px_noise", type=float, default=b["px_noise"])
    ap.add_argument("--stereo_sigma", type=float, default=b["stereo_sigma"])
    # Aday 3: karanlik karede olcum gurultusu m = 1 + (S_d - tau)/(1 - tau) ile carpilir (S_d > tau; tau = kapi esigi)
    ap.add_argument("--karanlik_sigma", type=int, default=int(b.get("karanlik_sigma", 0)))
    ap.add_argument("--karanlik_tau", type=float,
                    default=float((cfg.get("frontend") or {}).get("adaptive_clahe_s_d_threshold", 0.40)))
    ap.add_argument("--landmark_prior_sigma", type=float, default=b["landmark_prior_sigma"])
    ap.add_argument("--pnp_reproj", type=float, default=b["pnp_reproj"])
    # null -> th_depth_baselines * b (backend dump'taki tabandan)
    ap.add_argument("--max_init_depth", type=float, default=b.get("max_init_depth"))
    ap.add_argument("--th_depth_baselines", type=float, default=float(b.get("th_depth_baselines", 40.0)))
    ap.add_argument("--huber_mono_k", type=float, default=float(b.get("huber_mono_k", 1.345)))
    ap.add_argument("--huber_stereo_k", type=float, default=float(b.get("huber_stereo_k", 0.0)))
    ap.add_argument("--pnp_iters", type=int, default=int(b.get("pnp_iters", 200)))
    ap.add_argument("--pnp_confidence", type=float, default=float(b.get("pnp_confidence", 0.99)))
    ap.add_argument("--jump_pnp_max", type=float, default=float(b.get("jump_pnp_max", 0.0)))
    ap.add_argument("--gt_capa", type=int, default=int(b.get("gt_capa", 1)))
    ap.add_argument("--moba_ilk_kapi", type=int, default=int(b.get("moba_ilk_kapi", 0)))
    ap.add_argument("--moba_pencere_px", type=float, default=float(b.get("moba_pencere_px", 7.0)))
    ap.add_argument("--stereo_uzak_min", type=int, default=int(b.get("stereo_uzak_min", 0)))
    ap.add_argument("--stereo_ilk_hepsi", type=int, default=int(b.get("stereo_ilk_hepsi", 0)))
    ap.add_argument("--kurtarma_duzelt", type=int, default=int(b.get("kurtarma_duzelt", 0)))
    ap.add_argument("--izleme_basari", type=int, default=int(b.get("izleme_basari", 0)))
    ap.add_argument("--izleme_basari_min", type=int, default=int(b.get("izleme_basari_min", 0)))
    ap.add_argument("--kf_max_interval_s", type=float, default=float(b.get("kf_max_interval_s", 0.0)))
    ap.add_argument("--kf_odometri", type=int, default=int(b.get("kf_odometri", 1)))
    ap.add_argument("--disparity_bias", type=float, default=float(b.get("disparity_bias", 0.0)))
    ap.add_argument("--relinearize_skip", type=int, default=b["relinearize_skip"])
    ap.add_argument("--relinearize_threshold", type=float,
                    default=float(b.get("relinearize_threshold", 0.01)))
    ap.add_argument("--fixed_lag", type=int, default=int(b.get("fixed_lag", 0)))
    # Fixed-lag LM erken durma: pencere sıcak başladığından çoğu kare 1 iterasyonda yakınsar.
    ap.add_argument("--fl_lm_max_iters", type=int, default=int(b.get("fl_lm_max_iters", 5)))
    ap.add_argument("--fl_lm_rel_err_tol", type=float, default=float(b.get("fl_lm_rel_err_tol", 1e-3)))
    ap.add_argument("--fl_lm_lambda0", type=float, default=float(b.get("fl_lm_lambda0", 1e-5)))   # LM baslangic sonumu
    ap.add_argument("--cv_motion", type=int, default=int(b["cv_motion"]))
    ap.add_argument("--motion_trans_sigma", type=float, default=b["motion_trans_sigma"])
    ap.add_argument("--motion_rot_sigma", type=float, default=b["motion_rot_sigma"])
    # Hız-orantılı hareket prior'ı ('fixed' = sabit σ). 'adaptive':
    # σ_trans = clip(motion_sigma_rel·||prev_rel_t||, motion_sigma_min, motion_trans_sigma).
    ap.add_argument("--motion_sigma_mode", default=str(b.get("motion_sigma_mode", "fixed")),
                    choices=["fixed", "adaptive"])
    ap.add_argument("--motion_sigma_rel", type=float, default=float(b.get("motion_sigma_rel", 0.5)))
    ap.add_argument("--motion_sigma_min", type=float, default=float(b.get("motion_sigma_min", 0.02)))
    ap.add_argument("--extra_updates", type=int, default=b["extra_updates"])
    ap.add_argument("--motion_only_ba", type=int, default=int(b["motion_only_ba"]))
    ap.add_argument("--outlier_gate_px", type=float, default=b["outlier_gate_px"])
    # Sağlamlık kapıları (0 = kapalı): pnp_min_inliers PnP inlier eşiği; jump_max_trans öngörüden
    # bu kadar metreden fazla sapan pozu reddeder; chi2_gate_graph χ² outlier'larını grafa sokmaz.
    ap.add_argument("--pnp_min_inliers", type=int, default=int(b.get("pnp_min_inliers", 0)))
    ap.add_argument("--jump_max_trans", type=float, default=float(b.get("jump_max_trans", 0.0)))
    ap.add_argument("--chi2_gate_graph", type=int, default=int(b.get("chi2_gate_graph", 0)))
    # Graf kapısı eşiği = chi2_thresh_* × chi2_gate_mult; sıkı eşik sağlıklı stereo gözlemleri de atıyor.
    ap.add_argument("--chi2_gate_mult", type=float, default=float(b.get("chi2_gate_mult", 9.0)))
    # Sıçrama kapısı kilitlenme düzeltmeleri: jump_gate_ref=last_pose kapıyı son kabul edilmiş
    # poza göre uygular; resync_velocity=last_good resync'te hızı sıfırlamaz.
    # Cold-start: bootstrap_from_depth=1 harita boşken taze-stereo şartını gevşetir;
    # pnp_min_inlier_frac>0 mutlak PnP eşiğini oransal yapar (taban 6).
    ap.add_argument("--bootstrap_from_depth", type=int, default=int(b.get("bootstrap_from_depth", 0)))
    ap.add_argument("--pnp_min_inlier_frac", type=float, default=float(b.get("pnp_min_inlier_frac", 0.0)))
    ap.add_argument("--debug_frames", type=int, default=0,
                    help="Ilk N karede poz kestirim zincirini satir satir bas (cold-start teshisi).")
    # Yöntem sayaçları: her bayrağın kaç kez tetiklendiğini koşum sonunda bas (ölü-kod denetimi).
    ap.add_argument("--audit", type=int, default=0,
                    help="Kosum sonunda yontem tetiklenme sayaclarini bas (olu-kod denetimi).")
    ap.add_argument("--jump_gate_ref", default=str(b.get("jump_gate_ref", "prediction")),
                    choices=["prediction", "last_pose"])
    ap.add_argument("--resync_velocity", default=str(b.get("resync_velocity", "zero")),
                    choices=["zero", "last_good"])
    # İteratif χ² outlier reddi (motion-only BA içinde)
    ap.add_argument("--chi2_reject", type=int, default=int(b.get("chi2_reject", 0)))
    ap.add_argument("--chi2_iters", type=int, default=int(b.get("chi2_iters", 4)))
    ap.add_argument("--chi2_thresh_mono", type=float, default=float(b.get("chi2_thresh_mono", 5.991)))
    ap.add_argument("--chi2_thresh_stereo", type=float, default=float(b.get("chi2_thresh_stereo", 7.815)))
    # Keyframe + local BA
    ap.add_argument("--keyframe_mode", type=int, default=int(b.get("keyframe_mode", 0)))
    ap.add_argument("--kf_min_trans", type=float, default=float(b.get("kf_min_trans", 0.10)))
    ap.add_argument("--kf_min_rot", type=float, default=float(b.get("kf_min_rot", 0.15)))
    ap.add_argument("--kf_max_tracked_ratio", type=float, default=float(b.get("kf_max_tracked_ratio", 0.80)))
    ap.add_argument("--kf_min_interval", type=int, default=int(b.get("kf_min_interval", 2)))
    ap.add_argument("--kf_max_interval", type=int, default=int(b.get("kf_max_interval", 12)))
    ap.add_argument("--lba_window", type=int, default=int(b.get("lba_window", 0)))
    # Search-by-projection
    ap.add_argument("--search_by_projection", type=int, default=int(b.get("search_by_projection", 0)))
    # Asama B: kalici yerel harita + yerel harita takibi (ORB-SLAM3 TrackLocalMap); 0 = kapali (bit-ozdes)
    ap.add_argument("--yerel_harita", type=int, default=int(b.get("yerel_harita", 0)))
    ap.add_argument("--yh_yaricap_px", type=float, default=float(b.get("yh_yaricap_px", 4.0)))
    ap.add_argument("--yh_yaricap_kayip_px", type=float, default=float(b.get("yh_yaricap_kayip_px", 15.0)))
    ap.add_argument("--yh_desc_dist", type=float, default=float(b.get("yh_desc_dist", 0.78)))
    ap.add_argument("--yh_min_obs", type=int, default=int(b.get("yh_min_obs", 3)))
    ap.add_argument("--yh_prior_sigma", type=float, default=float(b.get("yh_prior_sigma", 0.1)))
    ap.add_argument("--yh_nn_orani", type=float, default=float(b.get("yh_nn_orani", 0.0)))   # ORB3: 0.8
    ap.add_argument("--lm_zayif_onsel_sigma", type=float, default=float(b.get("lm_zayif_onsel_sigma", 0.0)))
    ap.add_argument("--tekil_sabitle", type=int, default=int(b.get("tekil_sabitle", 0)))
    ap.add_argument("--motion_rot_uyarlamali", type=int, default=int(b.get("motion_rot_uyarlamali", 0)))
    ap.add_argument("--yh_konumla", type=int, default=int(b.get("yh_konumla", 0)))
    ap.add_argument("--kf_taze_aralik", type=int, default=int(b.get("kf_taze_aralik", 1)))
    ap.add_argument("--kf_harita_say", type=int, default=int(b.get("kf_harita_say", 0)))
    ap.add_argument("--smoother_tipi", default=str(b.get("smoother_tipi", "batch")))   # batch | isam2 | kayan
    ap.add_argument("--kayan_sabit_ek", type=float, default=float(b.get("kayan_sabit_ek", 0)))
    ap.add_argument("--kayan_sabit_n", type=int, default=int(b.get("kayan_sabit_n", 1)))
    ap.add_argument("--kayan_kf_aralik", type=int, default=int(b.get("kayan_kf_aralik", 0)))
    ap.add_argument("--kayan_yakin", type=float, default=float(b.get("kayan_yakin", 10)))
    ap.add_argument("--kayan_dusen_oncel", type=int, default=int(b.get("kayan_dusen_oncel", 1)))
    ap.add_argument("--be_aralik", type=int, default=int(b.get("be_aralik", 1)))
    ap.add_argument("--be_paralel", type=int, default=int(b.get("be_paralel", 0)))
    ap.add_argument("--moba_numpy", type=int, default=int(b.get("moba_numpy", 0)))
    ap.add_argument("--gozlem_seyrek", type=int, default=int(b.get("gozlem_seyrek", 1)))
    ap.add_argument("--gozlem_seyrek_min", type=int, default=int(b.get("gozlem_seyrek_min", 60)))
    ap.add_argument("--be_aralik_deg", type=float, default=float(b.get("be_aralik_deg", 0.0)))
    ap.add_argument("--be_aralik_m", type=float, default=float(b.get("be_aralik_m", 0.0)))
    ap.add_argument("--isam_kare_tahmin", type=int, default=int(b.get("isam_kare_tahmin", 0)))
    ap.add_argument("--sbp_reproj_px", type=float, default=float(b.get("sbp_reproj_px", 6.0)))
    ap.add_argument("--sbp_desc_dist", type=float, default=float(b.get("sbp_desc_dist", 0.7)))
    # IndeterminantLinearSystem sonrası smoother'ı yeniden kur (yumuşak sıfırlama, varsayılan kapalı).
    ap.add_argument("--smoother_recover", type=int, default=int(b.get("smoother_recover", 0)))
    # Çok-görüşlü üçgenleme: mono-only track'leri yeterli paralaks birikince DLT ile 3B'ye kavuştur.
    ap.add_argument("--mv_triangulation", type=int, default=int(b.get("mv_triangulation", 0)))
    ap.add_argument("--mv_min_obs", type=int, default=int(b.get("mv_min_obs", 3)))
    ap.add_argument("--mv_min_parallax_deg", type=float, default=float(b.get("mv_min_parallax_deg", 5.0)))
    ap.add_argument("--mv_max_reproj_err", type=float, default=float(b.get("mv_max_reproj_err", 2.0)))
    ap.add_argument("--mv_max_track_age", type=int, default=int(b.get("mv_max_track_age", 30)))
    # Aynı update'te çok yeni landmark eklenince smoother indeterminant oluyor; kare başına kabulü sınırla.
    ap.add_argument("--mv_max_per_frame", type=int, default=int(b.get("mv_max_per_frame", 3)))
    # DLT son iki tekil değer oranı (S[-2]/S[-1]) düşükse derinlik zayıf belirlenmiş; reddet.
    ap.add_argument("--mv_min_condition", type=float, default=float(b.get("mv_min_condition", 5.0)))
    # Çok-görüşlü init kipi: 'prior' DLT çözümünü 3B prior olarak çakar; 'factors' tampondaki
    # gözlemleri pencerede duran pozlara mono factor olarak ekler.
    ap.add_argument("--mv_init_mode", default=str(b.get("mv_init_mode", "prior")),
                    choices=["prior", "factors"])
    ap.add_argument("--mv_prior_sigma", type=float, default=float(b.get("mv_prior_sigma", 0.3)))
    # Landmark ayıklama (bkz. landmark_culling.py); anahtarlar birleşik `culling:` bölümünden,
    # hepsi varsayılan kapalı, sayaçları --audit 1 basar.
    cu = _culling_section(cfg, b)
    ap.add_argument("--cull_found_ratio", type=float, default=float(cu.get("cull_found_ratio", 0.0)))
    ap.add_argument("--cull_found_min_visible", type=int, default=int(cu.get("cull_found_min_visible", 8)))
    ap.add_argument("--cull_chi2_strikes", type=int, default=int(cu.get("cull_chi2_strikes", 0)))
    ap.add_argument("--cull_pnp_strikes", type=int, default=int(cu.get("cull_pnp_strikes", 0)))
    ap.add_argument("--cull_allow_reinit", type=int, default=int(cu.get("cull_allow_reinit", 1)))
    # Loop closure / PGO / GBA (varsayılan kapalı)
    lp = cfg.get("loop", {}) or {}
    ap.add_argument("--loop_enable", type=int, default=int(lp.get("enable", 0)))
    ap.add_argument("--loop_pgo", type=int, default=int(lp.get("pgo", 1)))
    ap.add_argument("--loop_gba", type=int, default=int(lp.get("gba", 0)))
    ap.add_argument("--loop_min_keyframe_gap", type=int, default=int(lp.get("min_keyframe_gap", 30)))
    ap.add_argument("--loop_appearance_sim", type=float, default=float(lp.get("appearance_sim", 0.80)))
    ap.add_argument("--loop_min_inliers", type=int, default=int(lp.get("min_inliers", 30)))
    ap.add_argument("--loop_max_candidates", type=int, default=int(lp.get("max_candidates", 3)))
    ap.add_argument("--loop_pgo_sigma", type=float, default=float(lp.get("pgo_loop_sigma", 0.05)))
    ap.add_argument("--loop_gba_iters", type=int, default=int(lp.get("gba_iters", 10)))
    # [TÜRBÜLANS] Ölçüm modeli kolları (turb/backend.py). turb_arm=A0 -> kol nesnesi yok, PyEco birebir.
    # Kollar yalnız pencere (smoother) ölçüm faktörlerinin kovaryansını değiştirir; robust çekirdek
    # (huber_mono_k / huber_stereo_k) BÜTÜN kollarda PyEco'nunkiyle aynıdır, kol başına ayarı yoktur.
    tb = cfg.get("turbulence", {}) or {}
    ap.add_argument("--turb_arm", default=str(tb.get("arm", "A0")),
                    help="A0 | A0C | A1 | A2 | A3 | A4 | A5 | AG | A0C+A3 (bkz. turb/backend.py)")
    ap.add_argument("--turb_corr_model", default=str(tb.get("corr_model", "empirical")))
    ap.add_argument("--turb_D_over_r0", type=float, default=float(tb.get("D_over_r0", 0.0)))
    ap.add_argument("--turb_theta_c_over_fov", type=float, default=float(tb.get("theta_c_over_fov", 0.15)))
    ap.add_argument("--turb_theta_c_nominal_over_fov", type=float,
                    default=float(tb.get("theta_c_nominal_over_fov", 0.15)))
    ap.add_argument("--turb_aperture_D", type=float, default=float(tb.get("aperture_D", 0.10)))
    ap.add_argument("--turb_wavelength", type=float, default=float(tb.get("wavelength", 550e-9)))
    ap.add_argument("--turb_px_per_rad", type=float, default=float(tb.get("px_per_rad", 170648.0)))
    ap.add_argument("--turb_image_width", type=int, default=int(tb.get("image_width", 752)))
    ap.add_argument("--turb_max_per_patch", type=int, default=int(tb.get("max_per_patch", 5)))
    ap.add_argument("--turb_a4_keep", type=int, default=int(tb.get("a4_keep_per_patch", 5)))
    ap.add_argument("--turb_shrinkage", type=float, default=float(tb.get("shrinkage", 0.05)))
    ap.add_argument("--turb_cal_factor", type=float, default=float(tb.get("cal_factor", 16.0)))
    ap.add_argument("--turb_ind_over_tilt", type=float, default=float(tb.get("ind_over_tilt", 1.75)))
    ap.add_argument("--turb_dimm_every", type=int, default=int(tb.get("dimm_update_every", 10)))
    ap.add_argument("--turb_dimm_window", type=int, default=int(tb.get("dimm_window", 60)))
    ap.add_argument("--turb_dimm_version", type=int, default=int(tb.get("dimm_version", 2)))
    ap.add_argument("--turb_dimm_sigma_scale", type=float, default=float(tb.get("dimm_sigma_scale", 1.0)))
    ap.add_argument("--turb_ag_window", type=int, default=int(tb.get("ag_window", 60)))
    ap.add_argument("--turb_gate", type=int, default=int(tb.get("gate", 0)),
                    help="uyarlanır izleme kapısı: 0 kapalı | 1 gürültü grubu | 2 PnP eşiği | 3 ikisi (turb/backend.py)")
    ap.add_argument("--turb_stats_out", default="")
    # [TÜRBÜLANS/E5] pencere-göreli NEES kovaryans dökümü (boş = kapalı; kestirimi değiştirmez)
    ap.add_argument("--cov_out", default="")
    ap.add_argument("--cov_every", type=int, default=1)


def backend_args(cfg=None, argv=None):
    """gtsam_backend.main için Namespace (config.yaml + CLI override)."""
    cfg = cfg or load()
    ap = argparse.ArgumentParser(description="GTSAM backend (config.yaml temelli)")
    _add_backend_args(ap, cfg)
    ap.add_argument("--out", default=traj_path(cfg))
    ap.add_argument("--refresh_every", type=int, default=cfg["backend"]["refresh_every"])
    return ap.parse_args(argv)


def live_args(cfg=None, argv=None):
    """run_vio_live için Namespace (backend + canlı görselleştirme; config.yaml + CLI)."""
    cfg = cfg or load()
    lv = cfg["live"]
    ap = argparse.ArgumentParser(description="Canlı VIO (Pangolin, config.yaml temelli)")
    _add_backend_args(ap, cfg)
    ap.add_argument("--out", default=traj_path(cfg, "_live"))
    ap.add_argument("--refresh_every", type=int, default=lv["refresh_every"])
    ap.add_argument("--delay", type=float, default=lv["delay"])
    ap.add_argument("--w", type=int, default=lv["window_w"])
    ap.add_argument("--h", type=int, default=lv["window_h"])
    return ap.parse_args(argv)
