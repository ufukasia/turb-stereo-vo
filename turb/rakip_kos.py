"""
Faz 2 rakip koşum hattı (ON_KAYIT_FAZ2 §1; ON_TEST_RAKIPLER.md): ORB-SLAM3 stereo, AirSLAM, DROID-SLAM.

Rakipler önceki makalenin koşucularıyla, DEĞİŞTİRİLMEDEN çağrılır
(`slam/slam_hakem_kosum/rakip/{orbslam3,airslam,droid}/run_*_asl.py`, sözleşme `rakip/SOZLESME.md`).
Bu modül yalnız (1) komutu kurar, (2) süre / yük ortalaması / GPU adını kaydeder, (3) pencere
değerlendirmesini PyEco ile BİREBİR aynı fonksiyonla (`turb/metrics.py::trajectory_metrics`) yapar.

Yapılandırmalar (ON_TEST_RAKIPLER.md §0):
  orbslam3  rakipler/ORB_SLAM3, stereo, nFeatures 1200 (EuRoC.yaml), döngü kapama kapalı, kamera hızında besleme
  airslam   airslam_ws_ek11c (yazar kodu d1265771 + ek yaması), YALNIZ ek.belirlenimci=2, vo_euroc.yaml, VO
  droid     rakipler/DROID-SLAM, stereo, test_euroc.py varsayılanları, çevrimiçi (nedensel) anahtar kare çıktısı
Birincil yörünge: gerçek zamanlı (nedensel); ikincil: son yörünge.

Tek koşum (ön-test ve elle kullanım):
    python -m turb.rakip_kos kos --rakip orbslam3 --seq <ASL> --kaynak <ham ASL> --sekans V1_01_easy \
        --dr 8 --tohum 0 --sigma 1 --protokol onek --out_dir <dizin> --kayit <jsonl> [--tekrar 1]
    python -m turb.rakip_kos kimlik                # rakip kimlik özeti (JSON)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

RAKIP_DIR = Path("/opt/rivals/wrappers")
RAKIPLER = Path("/opt/rivals")
PY_TURB = "python"
PY_DROID = os.path.expanduser("~/anaconda3/envs/droid/bin/python")
AIRSLAM_WS = RAKIPLER / "airslam_ws_ek11c"
ORB_ROOT = RAKIPLER / "ORB_SLAM3"
DROID_ROOT = RAKIPLER / "DROID-SLAM"
RAKIP_LISTE = ("orbslam3", "airslam", "droid")
#: kare başına poz veren sistemlerde 0.5 s; DROID çevrimiçi yalnız anahtar kare -> 2.0 s (önceki makalenin kuralı)
MAX_GAP_S = {"orbslam3": 0.5, "orbslam3_si": 0.5, "airslam": 0.5, "droid": 2.0}
#: Sapma 2 (b): ORB-SLAM3 stereo-ataletsel; ön kayıtlı listede değil, --rakipler orbslam3_si ile seçilir
RAKIP_EK = ("orbslam3_si",)
#: SAPMA 2 (ON_TEST_RAKIPLER.md §4.5, türbülanslı sonuç görülmeden): DROID çevrimiçi çıktısı yalnız anahtar karedir ve
#: DROID kamera durağanken (son anahtar kareye akış < eşik) anahtar kare EKLEMEZ. Kapsam bu yüzden sıfırıncı derece
#: tutmayla ölçülür: pencere karesi, kendinden önce bir anahtar kare varsa ve koşum tamamlandıysa (ya da kare son
#: anahtar kareden önceyse) kapsanır. Iraksama yayılım / max_err ile, çökme / ilklenememe yörünge yokluğuyla yakalanır.
KAPSAM_TUTMA = {"droid"}
#: başarısızlık ölçütü (ON_TEST_RAKIPLER.md §0; koşumdan önce sabitlendi)
KAPSAM_ESIK = 0.90
YAYILIM_BANT = (0.5, 2.0)
IRAKSAMA_M = 5.0
ZAMAN_ASIMI_S = 3600
#: GPU bellek kapısı [MB] (3090'da boş bellek bu değerin altındaysa bekle)
GPU_IHTIYAC_MB = {"airslam": 3500, "droid": 12000}
BEKLENEN_GPU = "NVIDIA GeForce RTX 3090"
WINDOW_LEN = 400


# =============================================================================
#  Yardımcılar
# =============================================================================
def sha256_dosya(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sekans_sha256(seq: Path) -> str:
    """Görüntü dosyalarının (cam0+cam1 PNG) sıralı sha256 listesinin sha256'sı + data.csv'ler."""
    h = hashlib.sha256()
    m = Path(seq) / "mav0"
    for cam in ("cam0", "cam1"):
        h.update(sha256_dosya(m / cam / "data.csv").encode())
        for f in sorted((m / cam / "data").iterdir()):
            h.update(f.name.encode())
            h.update(sha256_dosya(f).encode())
    return h.hexdigest()


def cam0_damgalari(seq: Path) -> list[int]:
    """Kanca/üretici kuralı: cam0 data.csv damgalarının SAYISAL sıralı listesi (mutlak kare indeksi)."""
    out = []
    with open(Path(seq) / "mav0/cam0/data.csv") as f:
        for ln in f:
            ln = ln.strip()
            if ln and not ln.startswith("#"):
                out.append(int(ln.split(",")[0]))
    return sorted(set(out))


def pencereler(sekans: str, path=None) -> list[int]:
    path = path or (ROOT / "evidence" / "pencereler_faz2_kapsam1m.json")
    return [int(s) for s in json.load(open(path))["sequences"][sekans]["starts"]]


def gpu_listesi() -> list[dict]:
    o = subprocess.run(["nvidia-smi", "--query-gpu=index,uuid,name,memory.free,pci.bus_id",
                        "--format=csv,noheader,nounits"], capture_output=True, text=True).stdout
    out = []
    for ln in o.strip().splitlines():
        i, u, n, fr, bus = [x.strip() for x in ln.split(",")]
        out.append({"index": int(i), "uuid": u, "name": n, "free_mb": int(fr), "bus": bus})
    return out


def gpu_3090() -> dict:
    g = [x for x in gpu_listesi() if "3090" in x["name"]]
    if len(g) != 1:
        raise SystemExit(f"RTX 3090 tek olarak bulunamadı: {gpu_listesi()}")
    return g[0]


def gpu_bekle(rakip: str, en_cok_s: float = 6 * 3600):
    ihtiyac = GPU_IHTIYAC_MB.get(rakip)
    if not ihtiyac:
        return None
    t0 = time.time()
    while True:
        g = gpu_3090()
        if g["free_mb"] >= ihtiyac:
            return g["free_mb"]
        if time.time() - t0 > en_cok_s:
            raise RuntimeError(f"3090'da {ihtiyac} MB boş bellek {en_cok_s:.0f} s içinde açılmadı")
        time.sleep(30)


def _git(d: Path) -> str | None:
    try:
        return subprocess.run(["git", "-C", str(d), "log", "-1", "--format=%H"], capture_output=True,
                              text=True, timeout=20).stdout.strip() or None
    except Exception:  # noqa: BLE001
        return None


def kimlik() -> dict:
    """Rakip kimlikleri: commit / ikili parmak izi, koşucu sha256, yapılandırma dosyaları sha256."""
    airslam_src = Path(open(AIRSLAM_WS / ".airslam_src").read().strip())
    k = {
        "orbslam3": {
            "kok": str(ORB_ROOT), "git": _git(ORB_ROOT),
            "lib_sha256": sha256_dosya(ORB_ROOT / "lib/libORB_SLAM3.so"),
            "ikili_sha256": sha256_dosya(ORB_ROOT / "Examples/Stereo/stereo_asl_headless"),
            "Tracking.cc_sha256": sha256_dosya(ORB_ROOT / "src/Tracking.cc"),
            "voc_sha256": sha256_dosya(ORB_ROOT / "Vocabulary/ORBvoc.txt"),
            "yazar_ayari_sha256": sha256_dosya(ORB_ROOT / "Examples/Stereo/EuRoC.yaml"),
            "kosucu": str(RAKIP_DIR / "orbslam3/run_orbslam3_asl.py"),
            "kosucu_sha256": sha256_dosya(RAKIP_DIR / "orbslam3/run_orbslam3_asl.py"),
            "ayar": "koşucu sensor.yaml'dan settings.yaml üretir (EuRoC.yaml ORB değerleri, nFeatures 1200, "
                    "loopClosing 0); üretilen dosya her koşumun yan dizininde (settings.yaml) ve sha256'sı kayıtta",
        },
        "airslam": {
            "ws": str(AIRSLAM_WS), "src": str(airslam_src), "src_git": _git(airslam_src),
            "yazar_commit": "d1265771026091db97ee8f7a6d80692f0edd5fef",
            "visual_odometry_sha256": sha256_dosya(AIRSLAM_WS / "devel/lib/air_slam/visual_odometry"),
            "vo_euroc.yaml_sha256": sha256_dosya(airslam_src / "configs/visual_odometry/vo_euroc.yaml"),
            "docker_imaj": "xukuanhit/air_slam:v4",
            "kosucu": str(RAKIP_DIR / "airslam/run_airslam_asl.py"),
            "kosucu_sha256": sha256_dosya(RAKIP_DIR / "airslam/run_airslam_asl.py"),
            "ek": {"belirlenimci": 2},
            "motorlar": {e: sha256_dosya(RAKIPLER / "airslam_cache/models_NVIDIA_GeForce_RTX_3090" / e)
                         for e in ("plnet_s0.engine", "plnet_s1.engine", "superpoint_lightglue.engine",
                                   "superpoint_v1_sim_int32.engine")
                         if (RAKIPLER / "airslam_cache/models_NVIDIA_GeForce_RTX_3090" / e).exists()},
        },
        "droid": {
            "kok": str(DROID_ROOT), "git": _git(DROID_ROOT),
            "agirlik_sha256": sha256_dosya(DROID_ROOT / "droid.pth"),
            "kosucu": str(RAKIP_DIR / "droid/run_droid_asl.py"),
            "kosucu_sha256": sha256_dosya(RAKIP_DIR / "droid/run_droid_asl.py"),
            "python": PY_DROID,
            "ayar": "test_euroc.py varsayilanlari (koşucu argümanları), stereo, --causal-out",
        },
        "orbslam3_si": {
            "kok": str(ORB_ROOT),
            "lib_sha256": sha256_dosya(ORB_ROOT / "lib/libORB_SLAM3.so"),
            "ikili_sha256": sha256_dosya(ORB_ROOT / "Examples/Stereo-Inertial/stereo_inertial_asl_headless"),
            "ikili_kaynak_sha256": sha256_dosya(ORB_ROOT / "Examples/Stereo-Inertial/stereo_inertial_asl_headless.cc"),
            "yazar_ayari_sha256": sha256_dosya(ORB_ROOT / "Examples/Stereo-Inertial/EuRoC.yaml"),
            "kosucu_sha256": sha256_dosya(RAKIP_DIR / "orbslam3/run_orbslam3_asl.py"),
            "ayar": "orbslam3 ile aynı + --imu (IMU_STEREO; IMU bloğu Stereo-Inertial/EuRoC.yaml, gürültü imu0/sensor.yaml)",
        },
        "bu_modul_sha256": sha256_dosya(Path(__file__)),
        "metrics.py_sha256": sha256_dosya(ROOT / "turb/metrics.py"),
    }
    try:
        k["airslam"]["docker_imaj_id"] = subprocess.run(
            ["docker", "image", "inspect", "-f", "{{.Id}}", "xukuanhit/air_slam:v4"],
            capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception:  # noqa: BLE001
        pass
    return k


def kimlik_ozeti(k: dict | None = None) -> str:
    """Kayıtlara giren kısa parmak izi (kimlik sözlüğünün sha256'sının ilk 12 hanesi)."""
    k = k or kimlik()
    return hashlib.sha256(json.dumps(k, sort_keys=True).encode()).hexdigest()[:12]


# =============================================================================
#  Komut ve koşum
# =============================================================================
def komut(rakip: str, seq: Path, out: Path, max_frames: int | None = None):
    seq = str(seq)
    if rakip in ("orbslam3", "orbslam3_si"):
        cmd = [PY_TURB, str(RAKIP_DIR / "orbslam3/run_orbslam3_asl.py"), "--timeout", str(ZAMAN_ASIMI_S),
               "--seq", seq, "--calib", seq, "--out", str(out)]
        if rakip == "orbslam3_si":       # IMU: üretilen sekanstaki özgün imu0 kopyası (uret_sekans kopyalar)
            cmd.append("--imu")
        env = {"ORBSLAM3_ROOT": str(ORB_ROOT)}
    elif rakip == "airslam":
        cmd = [PY_TURB, str(RAKIP_DIR / "airslam/run_airslam_asl.py"), "--ws", str(AIRSLAM_WS),
               "--no-scale-thresholds", "--ek", "belirlenimci=2", "--gpu", gpu_3090()["uuid"],
               "--timeout", str(ZAMAN_ASIMI_S), "--seq", seq, "--calib", seq, "--out", str(out)]
        env = {}
    elif rakip == "droid":
        glob_out = out.with_name(out.stem + "_global.csv")
        cmd = [PY_DROID, str(RAKIP_DIR / "droid/run_droid_asl.py"), "--no-viz", "--timeout", str(ZAMAN_ASIMI_S),
               "--seq", seq, "--calib", seq, "--out", str(glob_out), "--causal-out", str(out)]
        env = {"CUDA_DEVICE_ORDER": "PCI_BUS_ID", "CUDA_VISIBLE_DEVICES": "1"}
    else:
        raise ValueError(rakip)
    if max_frames:
        cmd += ["--max-frames", str(int(max_frames))]
    return cmd, env


def yorunge_yollari(rakip: str, out: Path) -> tuple[Path, Path]:
    """(nedensel = birincil, son = ikincil) yörünge dosyaları."""
    if rakip in ("orbslam3", "orbslam3_si"):
        return out.with_name(out.stem + "_nedensel.csv"), out
    if rakip == "airslam":
        c = out.with_name(out.stem + "_airslam") / "frames_raw_causal.csv"
        # koşum bitmediyse koşucu ana çıktıya akışla yazılmış nedensel pozları koyar
        return (c if c.exists() else out), out
    if rakip == "droid":
        return out, out.with_name(out.stem + "_global.csv")
    raise ValueError(rakip)


def kos(rakip: str, seq: Path, out: Path, max_frames: int | None = None) -> dict:
    """Rakibi bir kez koş; süre, yük ortalaması, GPU ve koşucunun yan dosyasını döndür."""
    out = Path(out).resolve()                # koşucular cwd=RAKIP_DIR ile çalışır: yollar MUTLAK olmalı
    seq = Path(os.path.abspath(seq))
    out.parent.mkdir(parents=True, exist_ok=True)
    for f in list(out.parent.glob(out.stem + ".*")) + list(out.parent.glob(out.stem + "_*")):  # yarım deneme
        shutil.rmtree(f) if f.is_dir() else f.unlink()
    bos = gpu_bekle(rakip)
    cmd, env = komut(rakip, seq, out, max_frames)
    full_env = {**os.environ, **env}
    for k_ in ("PYTHONPATH", "PYTHONHOME"):
        full_env.pop(k_, None)
    yuk0 = os.getloadavg()
    t0 = time.time()
    with open(out.with_suffix(".log"), "w") as lg:
        p = subprocess.run(cmd, stdout=lg, stderr=subprocess.STDOUT, env=full_env, cwd=str(RAKIP_DIR))
    sure = time.time() - t0
    yuk1 = os.getloadavg()
    side_p = out.with_suffix(".json")
    side = json.loads(side_p.read_text()) if side_p.exists() else {"durum": "cokme", "not": "yan dosya yok"}
    if rakip == "droid":
        g = out.with_name(out.stem + "_global.json")
        side_g = json.loads(g.read_text()) if g.exists() else {}
        side["ana_durum_global"] = side_g.get("durum")
        side["gpu_adi"] = side_g.get("gpu_adi")
        side["tepe_gpu_mb"] = side_g.get("tepe_gpu_mb")
    gpu_adi = None
    if rakip == "droid":
        gpu_adi = side.get("gpu_adi")
    elif rakip == "airslam":
        lg_ = out.with_name(out.stem + "_airslam") / "airslam.log"
        t = lg_.read_text(errors="ignore") if lg_.exists() else ""
        i = t.find("airslam_cache/models_")
        gpu_adi = t[i + len("airslam_cache/models_"):].split(",")[0].split()[0] if i >= 0 else None
    altyapi = None
    if rakip == "airslam" and (side.get("gpu_hatasi") or
                               (side.get("durum") == "ok" and side.get("n_islenen") and side.get("n_kare_girdi")
                                and side["n_islenen"] < 0.99 * side["n_kare_girdi"])):
        altyapi = "airslam GPU/TensorRT hatası ya da yarıda kesildi (sonuç değil)"
    if rakip == "droid" and "OOM" in str(side.get("not", "")):
        altyapi = "DROID CUDA bellek yetersiz (sonuç değil)"
    return {"rakip": rakip, "seq": str(seq), "out": str(out), "rc": p.returncode, "sure_s": round(sure, 2),
            "yuk_once": [round(x, 1) for x in yuk0], "yuk_sonra": [round(x, 1) for x in yuk1],
            "gpu_bos_mb_once": bos, "gpu_adi": gpu_adi, "durum": side.get("durum"), "not": side.get("not"),
            "n_kare_girdi": side.get("n_kare_girdi"), "n_poz": side.get("n_poz"), "fps": side.get("fps"),
            "tepe_gpu_mb": side.get("tepe_gpu_mb"), "altyapi_hatasi": altyapi, "komut": cmd}


# =============================================================================
#  Pencere değerlendirmesi (PyEco ile aynı metrik fonksiyonu)
# =============================================================================
def _kestirim_damgalari(csv: Path) -> np.ndarray:
    if not Path(csv).exists():
        return np.zeros(0, np.int64)
    t = []
    with open(csv) as f:
        for ln in f:
            if ln and ln[0].isdigit():
                t.append(int(ln.split(",")[0]))
    return np.unique(np.array(t, np.int64))


def kapsam(est_t: np.ndarray, win_t: np.ndarray, max_gap_s: float) -> float:
    """Pencere damgalarından kestirimin kapsadığı oran (tam damga ya da boşluğu <= max_gap iki kestirim arası)."""
    if len(est_t) == 0:
        return 0.0
    tam = np.isin(win_t, est_t)
    j = np.searchsorted(est_t, win_t, side="right")
    ic = (j > 0) & (j < len(est_t))
    jj = np.clip(j, 1, len(est_t) - 1)
    ara = ic & ((est_t[jj] - est_t[jj - 1]) <= max_gap_s * 1e9)
    return float(np.mean(tam | ara))


def _yayilim(est_csv: Path, gt_csv: Path, t0: float, t1: float) -> float | None:
    from turb.metrics import load_traj, load_gt, associate
    et, ex, eq = load_traj(est_csv)
    m = (et >= t0) & (et <= t1)
    if m.sum() < 5:
        return None
    gt_t, gx, gq = load_gt(gt_csv)
    a = associate(et[m], ex[m], eq[m], gt_t, gx, gq)
    if a is None:
        return None
    _, exyz, _, gxyz, _ = a
    r_e = np.sqrt(((exyz - exyz.mean(0)) ** 2).sum(1).mean())
    r_g = np.sqrt(((gxyz - gxyz.mean(0)) ** 2).sum(1).mean())
    return float(r_e / r_g) if r_g > 0 else None


def kapsam_tutma(est_t: np.ndarray, win_t: np.ndarray, kosum_tamam: bool) -> float:
    """Sıfırıncı derece tutma: ilk kestirimden sonraki ve (koşum yarıda kaldıysa) son kestirimden önceki kareler."""
    if len(est_t) == 0:
        return 0.0
    ic = win_t >= est_t[0]
    if not kosum_tamam:
        ic &= win_t <= est_t[-1]
    return float(np.mean(ic))


def pencere_degerlendir(rakip: str, est_csv: Path, gt_csv: Path, ts_all: list[int], w: int,
                        kosum_tamam: bool = True) -> dict:
    """Tek pencere: trajectory_metrics (PyEco) + kapsam + yayılım + başarısızlık kodu."""
    from turb.metrics import trajectory_metrics
    i0, i1 = int(w), int(min(w + WINDOW_LEN, len(ts_all)) - 1)
    t0, t1 = float(ts_all[i0]), float(ts_all[i1])
    win_t = np.array(ts_all[i0:i1 + 1], np.int64)
    est_t = _kestirim_damgalari(est_csv)
    r = {"t_win_start": t0, "t_win_end": t1}
    try:
        m = trajectory_metrics(str(est_csv), str(gt_csv), t_min=t0, t_max=t1) if len(est_t) else {"valid": False}
    except Exception as e:  # noqa: BLE001
        m = {"valid": False, "hata": f"{type(e).__name__}: {e}"}
    r.update(m)
    r["kapsam_p"] = (kapsam_tutma(est_t, win_t, kosum_tamam) if rakip in KAPSAM_TUTMA
                     else kapsam(est_t, win_t, MAX_GAP_S[rakip]))
    if rakip in KAPSAM_TUTMA:
        r["kapsam_p_bosluk2s"] = kapsam(est_t, win_t, MAX_GAP_S[rakip])   # eski tanım, bilgi için
    r["yayilim"] = _yayilim(est_csv, gt_csv, t0, t1) if m.get("valid") else None
    neden = []
    if not m.get("valid"):
        neden.append("pencerede < 5 kestirim (yörünge yok)")
    if r["kapsam_p"] < KAPSAM_ESIK:
        neden.append(f"kapsam {r['kapsam_p']:.3f} < {KAPSAM_ESIK}")
    if r["yayilim"] is not None and not (YAYILIM_BANT[0] <= r["yayilim"] <= YAYILIM_BANT[1]):
        neden.append(f"yayilim {r['yayilim']:.3f} bant disi")
    if m.get("valid") and m.get("max_err", 0) > IRAKSAMA_M:
        neden.append(f"max_err {m['max_err']:.2f} m > {IRAKSAMA_M}")
    r["durum_p"] = "basarisiz" if neden else "ok"
    r["neden"] = "; ".join(neden)
    return r


def degerlendir(rakip: str, out: Path, kaynak_seq: Path, wins: list[int], kosum_tamam: bool = True) -> list[dict]:
    """Her pencere için birincil (nedensel, m_*) ve ikincil (son, ms_*) metrikler."""
    gt = Path(kaynak_seq) / "mav0/state_groundtruth_estimate0/data.csv"
    ts_all = cam0_damgalari(kaynak_seq)
    ned, son = yorunge_yollari(rakip, Path(out))
    recs = []
    for w in wins:
        a = pencere_degerlendir(rakip, ned, gt, ts_all, w, kosum_tamam)
        b = pencere_degerlendir(rakip, son, gt, ts_all, w, kosum_tamam)
        rec = {"window": int(w), "durum_p": a["durum_p"], "neden": a["neden"],
               "yorunge_nedensel": str(ned), "yorunge_son": str(son)}
        rec.update({"m_" + k: v for k, v in a.items() if k not in ("durum_p", "neden")})
        rec.update({"ms_" + k: v for k, v in b.items()})
        recs.append(rec)
    return recs


# =============================================================================
#  Pencere görünümü (protokol (a): yalnız o pencerenin kareleri; görüntüler kopyalanmaz)
# =============================================================================
def pencere_gorunumu(seq: Path, hedef: Path, w: int) -> Path:
    """seq'in [w, w+400) mutlak karelerini içeren ASL görünümü (data.csv süzülür, data/ sembolik bağ)."""
    ts_kaynak = None
    meta = Path(seq) / "turb_meta.json"
    if meta.exists():
        ts_kaynak = cam0_damgalari(Path(json.loads(meta.read_text())["kaynak"]["yol"]))
    ts_all = ts_kaynak or cam0_damgalari(seq)
    sec = set(ts_all[w:w + WINDOW_LEN])
    if hedef.exists():
        shutil.rmtree(hedef)
    for cam in ("cam0", "cam1"):
        d = hedef / "mav0" / cam
        d.mkdir(parents=True)
        os.symlink(Path(seq) / "mav0" / cam / "data", d / "data")
        shutil.copyfile(Path(seq) / "mav0" / cam / "sensor.yaml", d / "sensor.yaml")
        with open(Path(seq) / "mav0" / cam / "data.csv") as f, open(d / "data.csv", "w") as g:
            for ln in f:
                if ln.startswith("#") or (ln[0].isdigit() and int(ln.split(",")[0]) in sec):
                    g.write(ln)
    return hedef


# =============================================================================
#  CLI
# =============================================================================
def kayit_yaz(path: Path, recs: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as f:
        for r in recs:
            f.write(json.dumps(r, ensure_ascii=False, default=lambda o: o.item() if hasattr(o, "item") else str(o))
                    + "\n")


def tek_kosum(rakip, seq, kaynak, sekans, dr, tohum, sigma, protokol, out_dir, kayit, tekrar=0,
              wins=None, max_frames=None, etiket="", kimlik_fp=None, ek=None) -> list[dict]:
    seq, kaynak, out_dir = Path(os.path.abspath(seq)), Path(os.path.abspath(kaynak)), Path(out_dir).resolve()
    wins = list(wins) if wins is not None else pencereler(sekans)
    out = out_dir / rakip / f"{etiket or seq.name}_r{tekrar}.csv"
    info = kos(rakip, seq, out, max_frames)
    recs = degerlendir(rakip, out, kaynak, wins, kosum_tamam=(info["durum"] == "ok"))
    base = {"rakip": rakip, "sequence": sekans, "D_over_r0": float(dr), "seed": int(tohum),
            "sensor_noise_dn": float(sigma), "protokol": protokol, "tekrar": int(tekrar),
            "seq_dizin": str(seq), "kimlik_fp": kimlik_fp or kimlik_ozeti(), "kosum": info,
            "zaman": time.strftime("%Y-%m-%d %H:%M:%S"), **(ek or {})}
    out_recs = [{**base, **r} for r in recs]
    kayit_yaz(Path(kayit), out_recs)
    for r in out_recs:
        a = r.get("m_ate_se3")
        print(f"[{rakip}] {sekans} dr{dr:g} s{tohum} {protokol} w{r['window']} r{tekrar}: "
              f"durum={info['durum']} ATE={'NA' if a is None else f'{a:.4f}'} kapsam={r['m_kapsam_p']:.3f} "
              f"{r['durum_p']} {r['neden']} ({info['sure_s']:.0f} s, yük {info['yuk_once'][0]})", flush=True)
    return out_recs


def yeniden_degerlendir(kayit: Path):
    """Koşum bilgisini koruyarak pencere metriklerini güncel kurallarla yeniden hesapla (eski dosya .eski_<zaman>)."""
    R = [json.loads(l) for l in open(kayit)]
    yedek = kayit.with_name(kayit.name + time.strftime(".eski_%Y%m%d_%H%M%S"))
    shutil.copyfile(kayit, yedek)
    out = []
    for r in R:
        k = r["kosum"]
        w = degerlendir(r["rakip"], Path(k["out"]), Path("/data/sdata") / r["sequence"],
                        [r["window"]], kosum_tamam=(k.get("durum") == "ok"))[0]
        nr = {kk: vv for kk, vv in r.items() if not (kk.startswith("m_") or kk.startswith("ms_")
                                                    or kk in ("durum_p", "neden", "window"))}
        out.append({**nr, **w})
    with open(kayit, "w") as f:
        for r in out:
            f.write(json.dumps(r, ensure_ascii=False, default=lambda o: o.item() if hasattr(o, "item") else str(o))
                    + "\n")
    print(f"{len(out)} kayıt yeniden değerlendirildi; eski: {yedek}")


def main():
    ap = argparse.ArgumentParser(description="Faz 2 rakip koşum hattı")
    sp = ap.add_subparsers(dest="komut", required=True)
    k = sp.add_parser("kos")
    k.add_argument("--rakip", required=True, choices=RAKIP_LISTE)
    k.add_argument("--seq", required=True, help="koşulacak ASL dizini (üretilmiş ya da ham)")
    k.add_argument("--kaynak", required=True, help="ham kaynak ASL (GT ve mutlak kare indeksi)")
    k.add_argument("--sekans", required=True, help="pencere listesindeki ad (ör. V1_01_easy)")
    k.add_argument("--dr", type=float, required=True)
    k.add_argument("--tohum", type=int, required=True)
    k.add_argument("--sigma", type=float, required=True)
    k.add_argument("--protokol", default="onek", help="onek | pencere | bütun | ...")
    k.add_argument("--out_dir", required=True)
    k.add_argument("--kayit", required=True)
    k.add_argument("--tekrar", type=int, default=0)
    k.add_argument("--pencereler", default=None, help="virgüllü başlangıçlar (vars.: kapsam1m listesi)")
    k.add_argument("--max_frames", type=int, default=None)
    k.add_argument("--etiket", default="")
    sp.add_parser("kimlik")
    y = sp.add_parser("yeniden", help="kayıtların pencere metriklerini yörünge dosyalarından yeniden hesapla")
    y.add_argument("--kayit", required=True)
    a = ap.parse_args()
    if a.komut == "yeniden":
        yeniden_degerlendir(Path(a.kayit))
        return
    if a.komut == "kimlik":
        kk = kimlik()
        print(json.dumps({"kimlik_fp": kimlik_ozeti(kk), **kk}, indent=1, ensure_ascii=False))
        return
    wins = [int(x) for x in a.pencereler.split(",")] if a.pencereler else None
    tek_kosum(a.rakip, a.seq, a.kaynak, a.sekans, a.dr, a.tohum, a.sigma, a.protokol, a.out_dir, a.kayit,
              a.tekrar, wins, a.max_frames, a.etiket)


if __name__ == "__main__":
    main()
