"""
AirSim kol koşucusu: PyEco kol koşucusunu (`pyeco_turb/turb/runner.py`) AirSim verisine bağlar (H3 + HA4).

    # tek hücre (duman testi; açık döküm yolu):
    python kol_kos.py hucre --sekans AF_F5_h5 --C0 1e-11 --tohum 0 --pencere 20:857 \
        --dizin cikti_izgara/gpu1/AF_F5_h5__C01e-11_s0_g1t0_k20-857 --kollar A0 AG A3 A0C+A3 --k 12 --etiket duman
    # ızgara (ana ızgaranın ön-uç dökümleri yeniden kullanılır):
    python kol_kos.py izgara --isci 3 [--bekle]
    # k kalibrasyonu (yalnız temiz pilot verisi):
    python kol_kos.py kalibre --isci 3
    # tam yeniden üretim (uret_derinlik + runner ön-ucu + kollar) — bit-özdeşlik testi:
    python kol_kos.py uret --sekans AF_F15_h5 --C0 0 --tohum 10 --pencere 20:857 --kollar A0 --etiket bi2

BAĞLANTI (dosya değişikliği YOK; runner.py/gtsam_backend.py/config.yaml ve code_fp 08e28b83407c aynen):
- Hücre = runner.Cell(dataset="airsim", base_path=<GT iskeleti>, sequence=<türbülanslı sekans adı>,
  windows=(0,), window_len=<kare sayısı>). Türbülans + sensör gürültüsü sekansa ÖNCEDEN işlenmiştir
  (turb.uret_derinlik, ana ızgarayla aynı); bu yüzden D_over_r0 = 0 ve sensor_noise_dn = 0 -> ön-uç kancası
  hiç içe aktarılmaz (PyEco ön-ucu birebir).
- Kamera iç parametreleri ve taban çizgisi arka-uca ön-uç dökümünden gelir (rect_K, baseline; PyEco'nun kendi
  yolu). Kolların türbülans parametreleri: --turb_px_per_rad = rektifiye fx, --turb_image_width = 1024,
  --turb_aperture_D = 14.4 mm (F5) / 4.8 mm (F15), --turb_cal_factor = k (A0C, A0C+A3).
- Arka-uç `runner._run_arms` ile koşar (EuRoC ızgarasıyla aynı fonksiyon): nedensel yörünge, düzleştirilmiş
  yörünge, pencere-göreli kovaryans dökümü (cov_every = 2), kol istatistikleri.
- AirSim metrikleri ana ızgaranınkiyle AYNI fonksiyonlarla eklenir: turb_izgara.metrik (ATE SE3, Sim3 ölçek) ve
  rpe_ozet.ozet (RPE1s/RPE10s) düzleştirilmiş yörünge (`traj_<kol>.csv` = ana ızgaranın gtsam_traj_*.csv'si) üzerinde.
- A0 her hücrede ana ızgaranın A0 yörüngesiyle bayt bayt karşılaştırılır (`a0_bit_ayni`).
- Nokta başına artık dökümü: artik_dok.py (ön-uç dökümü + GT; arka-uçtan bağımsız) -> <hücre>/artik_gt.npz.

ÇIKTI: cikti_kollar/<etiket>/kayit.jsonl ("tur": "kol" | "artik"), hücre dizinleri cikti_kollar/<etiket>/<sekans>/.
Aynı turb_sekans adı yalnız etiket dizini içinde tekildir; ana ızgaranın cikti_izgara/ dizinine YAZILMAZ.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PE = HERE.parent / "pyeco_turb"
PY = str(HERE.parent / "kod/turbvio/bin/python")
VERI = Path("/data/sdata_airsim")
CIKTI = HERE / "cikti_kollar"
TURB_KOK_KOLLAR = Path("/data/sdata_turb_airsim/kollar")    # ana ızgaranın gpu{G}/ kökünden AYRI
EUROC_KOK = "/data/sdata"      # yalnız config.load()'un varsayılan-sekans denetimi için
DTM = HERE / "africa/harita_africa.npz"
D_FOV = {15: 4.8e-3, 5: 14.4e-3}
BEKLENEN_CODE_FP = "08e28b83407c"
BEKLENEN_GPU = "NVIDIA GeForce RTX 3090"
KOLLAR = ("A0", "AG", "A3", "A0C+A3")
COV_EVERY = 2
GURULTU_DN = 1.0

# Izgara (SAPMA_TASLAK_KOLLAR.md): F5, C0 ∈ {0, 1e-12, 1e-11}, tohum 10-14, W1+W2 (ana ızgaranın dökümleri)
IZGARA_SEKANS = "AF_F5_h5"
IZGARA_C0 = (0.0, 1e-12, 1e-11)
IZGARA_TOHUM = (10, 11, 12, 13, 14)
PENCERELER = ("20:857", "1157:1994")
# k kalibrasyonu: temiz F5 pilot tohumları 0-2 (doğrulayıcı analize girmez), W1+W2, yalnız RTX 3090 dökümleri
K_IZGARA = (1, 2, 4, 8, 12, 16, 24)
K_TOHUM = (0, 1, 2)
K_DOSYA = CIKTI / "k_airsim.json"


# ----------------------------------------------------------------------------- yardımcılar
def sha(p, alg="sha1"):
    h = hashlib.new(alg)
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def kod_ozeti() -> dict:
    return {f: sha(HERE / f, "sha256")[:16] for f in ("kol_kos.py", "artik_dok.py", "turb_izgara.py", "rpe_ozet.py")}


def kamera(sekans: str) -> dict:
    import yaml
    c0 = yaml.safe_load(open(VERI / sekans / "mav0/cam0/sensor.yaml"))
    fov = int(round(json.load(open(VERI / sekans / "kayit_meta.json"))["arguman"]["fov"]))
    return {"W": int(c0["resolution"][0]), "H": int(c0["resolution"][1]), "fov": fov, "D": D_FOV[fov]}


def gt_iskelet(veri_kok: Path, sekans: str, ad: str, pencere: str | None) -> Path:
    """uret_derinlik ile BAYT BAYT aynı GT kesiti (başlık + kaynak satırları [a, b))."""
    out = veri_kok / ad / "mav0/state_groundtruth_estimate0/data.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    gl = list(open(VERI / sekans / "mav0/state_groundtruth_estimate0/data.csv"))
    gh = [l for l in gl if l.startswith("#")]
    gd = [l for l in gl if not l.startswith("#")]
    a, b = (0, len(gd)) if not pencere else (int(x) for x in pencere.split(":"))
    txt = "".join(gh + [gd[i] for i in range(a, b)])
    if not out.exists() or out.read_text() != txt:
        out.write_text(txt)
    return out


def anees_bandi(n: int, dof: int = 3):
    from scipy.stats import chi2
    if not n:
        return [float("nan"), float("nan")]
    return [float(chi2.ppf(0.025, dof * n) / n), float(chi2.ppf(0.975, dof * n) / n)]


def turb_adi(sekans, C0, tohum, pencere):
    """uret_derinlik'in adlandırması (gürültü 1 DN, gürültü tohumu = tohum)."""
    sys.path.insert(0, str(PE))
    from turb.derinlik_model import DerinlikParam
    p = DerinlikParam(C0=float(C0), tohum=int(tohum), D=kamera(sekans)["D"])
    a, b = pencere.split(":")
    return f"{sekans}__{p.anahtar()}_g{GURULTU_DN:g}t{tohum}_k{a}-{b}"


# ----------------------------------------------------------------------------- hücre tanımları (spec)
def spec_ana(sekans, C0, tohum, pencere):
    """Ana ızgara (cikti_izgara/ana_v03) kaydından spec; yoksa None. Ezilme denetimi: ATE yeniden hesaplanır."""
    from turb_izgara import metrik
    for f in sorted((HERE / "cikti_izgara/ana_v03").glob("w*.jsonl")):
        for l in open(f):
            r = json.loads(l)
            if (r.get("durum") == "tamam" and r["sekans"] == sekans and float(r["C0"]) == float(C0)
                    and int(r["tohum"]) == int(tohum) and r["pencere"] == pencere):
                od = HERE / r["cikti_dizin"]
                ad = r["turb_sekans"]
                traj = od / f"gtsam_traj_{ad}.csv"
                gt = gt_iskelet(CIKTI / "_gt_denetim", sekans, ad, pencere)
                ate = metrik(traj, gt)["ate_se3"]
                if abs(ate - r["ate_se3"]) > 1e-9 * max(1.0, r["ate_se3"]):
                    raise SystemExit(f"ana ızgara yörüngesi ezilmiş: {od}")
                if str(r.get("gpu")) != "1":
                    raise SystemExit(f"ana ızgara kaydı GPU1 değil: {r.get('gpu')}")
                return {"sekans": sekans, "C0": float(C0), "tohum": int(tohum), "pencere": pencere, "turb_sekans": ad,
                        "dump": str(od / f"landmarks_{ad}.pkl"), "ana_traj": str(traj), "ana_ate": r["ate_se3"],
                        "kaynak": "ana_v03", "turb_dilim": (r.get("turb") or {}).get("dilim"),
                        "uretim_kod": (r.get("turb") or {}).get("kod")}
    return None


def spec_dizin(sekans, C0, tohum, pencere, dizin, kaynak="pilot"):
    od = Path(dizin) if Path(dizin).is_absolute() else HERE / dizin
    ad = od.name
    if ad != turb_adi(sekans, C0, tohum, pencere):
        raise SystemExit(f"dizin adı hücreyle uyuşmuyor: {ad} != {turb_adi(sekans, C0, tohum, pencere)}")
    traj = od / f"gtsam_traj_{ad}.csv"
    return {"sekans": sekans, "C0": float(C0), "tohum": int(tohum), "pencere": pencere, "turb_sekans": ad,
            "dump": str(od / f"landmarks_{ad}.pkl"), "ana_traj": str(traj) if traj.exists() else None,
            "ana_ate": None, "kaynak": kaynak, "turb_dilim": None, "uretim_kod": None}


# ----------------------------------------------------------------------------- işçi: tek kol
def _ortam_hazirla():
    if os.environ.get("SLAM_CONFIG"):
        raise SystemExit("SLAM_CONFIG tanımlı olmamalı (code_fp config.yaml'ı özetler)")
    os.environ["SLAM_DATA_ROOT"] = EUROC_KOK     # config.load() varsayılan MH01 denetimi; base_path hücreden ezilir
    os.chdir(PE)
    for p in (str(PE), str(HERE)):
        if p not in sys.path:
            sys.path.insert(0, p)


def _hucre(spec, kol, k, veri_kok, n_kare, fx, B):
    from turb.runner import Cell
    cam = kamera(spec["sekans"])
    return Cell(dataset="airsim", base_path=str(veri_kok), sequence=spec["turb_sekans"], seed=int(spec["tohum"]),
                D_over_r0=0.0, arms=(kol,), windows=(0,), window_len=int(n_kare), nees=True, cov_every=COV_EVERY,
                px_per_rad=float(fx), image_width=cam["W"], stereo_baseline=float(B), frame_dt=0.05,
                sensor_noise_dn=0.0, psf_legacy_scale=False,
                backend_argv=("--turb_cal_factor", str(float(k)), "--turb_aperture_D", str(cam["D"])),
                tag=f"airsim_kollar_C0{float(spec['C0']):g}_p{spec['pencere'].replace(':', '-')}")


def _cov_ozet(cov):
    import numpy as np
    d = np.load(cov)
    C = d["cov"]
    tp = np.trace(C[:, 3:, 3:], axis1=1, axis2=2)
    tr = np.trace(C[:, :3, :3], axis1=1, axis2=2)
    return {"cov_n": int(len(C)), "cov_poz_iz_med": float(np.median(tp)), "cov_don_iz_med": float(np.median(tr))}


def kol_isi(a):
    """(spec, kol, k, etiket) -> kayıt. Arka-uç yalnız CPU."""
    spec, kol, k, etiket = a
    _ortam_hazirla()
    import pickle
    import config
    from turb import runner
    from turb_izgara import metrik
    from rpe_ozet import ozet
    t0 = time.time()
    hd = CIKTI / etiket / spec["turb_sekans"]
    rd = hd / (kol.replace("+", "p") + (f"_k{float(k):g}" if kol in ("A0C", "A0C+A3") else ""))
    rd.mkdir(parents=True, exist_ok=True)
    veri_kok = CIKTI / etiket / "veri"
    gt = gt_iskelet(veri_kok, spec["sekans"], spec["turb_sekans"], spec["pencere"])
    dl = rd / f"landmarks_{spec['turb_sekans']}.pkl"
    if dl.is_symlink() or dl.exists():
        dl.unlink()
    os.symlink(spec["dump"], dl)
    with open(spec["dump"], "rb") as f:
        d = pickle.load(f)
    f_ts = [int(fr["timestamp"]) for fr in d["frames"]]
    fx, B = float(d["rect_K"][0]), float(d["baseline"])
    del d
    # duman testi kısaltması: arka-uç yalnız ilk n kareyi koşar (--max_frames); metrik GT'si de kesilir
    n_sinir = spec.get("n_kare_sinir")
    n_kare = len(f_ts) if not n_sinir else min(len(f_ts), int(n_sinir))
    pen_m = spec["pencere"]
    gt_m = gt
    if n_kare < len(f_ts):
        a0_, _b = (int(x) for x in spec["pencere"].split(":"))
        pen_m = f"{a0_}:{a0_ + n_kare}"
        gt_m = gt_iskelet(CIKTI / etiket / "veri_kesit", spec["sekans"], spec["turb_sekans"], pen_m)
    cell = _hucre(spec, kol, k, veri_kok, n_kare, fx, B)
    cfg = runner._build_cfg(cell, str(rd))
    assert config.gt_path(cfg) == str(gt), (config.gt_path(cfg), gt)
    assert config.dump_path(cfg) == str(dl)
    base = {"t_frontend": 0.0, "fe_gpu": BEKLENEN_GPU if spec["kaynak"] != "uret" else spec.get("fe_gpu"),
            "frontend_stats": {}, "quality": None}
    rec = runner._run_arms(cell, cfg, str(rd), str(gt), f_ts, base)[0]
    rec["tur"] = "kol"
    rec["airsim"] = dict(spec, kol=kol, k=float(k), etiket=etiket, fx=fx, B=B, n_kare=n_kare, n_kare_dump=len(f_ts),
                         dump_sha1=sha(spec["dump"]), kod=kod_ozeti(), sure_s=time.time() - t0)
    rec["airsim"].pop("turb_dilim", None)
    traj = rd / f"traj_{kol}.csv"
    if traj.exists():
        rec["traj_sha1"] = sha(traj)
        try:
            rec.update({"a_" + kk: v for kk, v in metrik(traj, gt_m).items()})
            rec.update({"a_" + kk: v for kk, v in ozet(traj, spec["sekans"], pen_m.replace(":", "-")).items()})
        except Exception as e:                                       # noqa: BLE001
            rec["status"] = "fail"; rec["error"] = f"airsim metrik: {e}"
        if kol == "A0" and spec.get("ana_traj") and Path(spec["ana_traj"]).exists() and n_kare == len(f_ts):
            rec["a0_bit_ayni"] = (sha(spec["ana_traj"]) == rec["traj_sha1"])
    cov = rd / f"cov_{kol}.npz"
    if cov.exists():
        rec.update(_cov_ozet(cov))
    if rec.get("n_n"):
        rec["anees_bant"] = anees_bandi(int(rec["n_n"]))
    if kol == "A0" and traj.exists():                               # ikincil: kestirilen-poz artığı
        import artik_dok as A
        try:
            r = A.artik_hesapla(spec["dump"], lambda T: A.pozlar_traj(str(traj), T))
            rec["artik_est"] = A.ozetle(r)
        except Exception as e:                                       # noqa: BLE001
            rec["artik_est"] = {"hata": str(e)}
    return rec


def artik_isi(a):
    """(spec, etiket) -> GT-pozlu nokta başına artık (HA4 birincil); npz hücre dizinine."""
    spec, etiket = a
    _ortam_hazirla()
    import artik_dok as A
    t0 = time.time()
    hd = CIKTI / etiket / spec["turb_sekans"]
    hd.mkdir(parents=True, exist_ok=True)
    g = A.gt_kesit(spec["sekans"], spec["pencere"])
    r = A.artik_hesapla(spec["dump"], lambda T: A.pozlar_gt(g, T))
    A.kaydet(r, str(hd / "artik_gt.npz"))
    return {"tur": "artik", "key": spec["turb_sekans"], "status": "ok", "airsim": dict(spec, etiket=etiket,
            dump_sha1=sha(spec["dump"]), kod=kod_ozeti(), sure_s=time.time() - t0),
            "artik_gt": A.ozetle(r), "npz": str(hd / "artik_gt.npz")}


# ----------------------------------------------------------------------------- işçi: üretim + ön-uç (GPU1)
def uret(spec_temel, etiket, isci=6, sakla=False):
    """uret_derinlik (ana ızgarayla aynı bağımsız değişkenler) + runner ön-ucu (RTX 3090) -> spec."""
    sekans, C0, tohum, pen = spec_temel["sekans"], spec_temel["C0"], spec_temel["tohum"], spec_temel["pencere"]
    kok = TURB_KOK_KOLLAR / etiket
    kok.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    r = subprocess.run([PY, "-m", "turb.uret_derinlik", "--girdi", str(VERI / sekans), "--cikti_kok", str(kok),
                        "--C0", str(float(C0)), "--tohum", str(int(tohum)), "--dtm", str(DTM),
                        "--D", str(kamera(sekans)["D"]), "--isci", str(isci), "--kareler", pen,
                        "--gurultu", str(GURULTU_DN)], cwd=PE, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("üretim hatası:\n" + r.stderr[-2000:])
    ad = json.loads(r.stdout.strip().splitlines()[-1])["cikti"].split("/")[-1]
    assert ad == turb_adi(sekans, C0, tohum, pen), ad
    t1 = time.time()
    import multiprocessing as mp
    with mp.get_context("spawn").Pool(1) as pool:
        meta = pool.apply(onuc_isi, ((str(kok), ad, sekans, C0, tohum, pen, etiket),))
    gt_gen = kok / ad / "mav0/state_groundtruth_estimate0/data.csv"
    gt_isk = gt_iskelet(CIKTI / etiket / "veri", sekans, ad, pen)
    meta.update(uretim_s=t1 - t0, gt_bayt_ayni=(sha(gt_gen) == sha(gt_isk)),
                turb_meta=json.load(open(kok / ad / "turb_meta.json")).get("kod_ozeti"))
    if not sakla:
        import shutil
        shutil.rmtree(kok / ad, ignore_errors=True)
    spec = dict(spec_temel, turb_sekans=ad, kaynak="uret", **meta)
    return spec


def onuc_isi(a):
    kok, ad, sekans, C0, tohum, pen, etiket = a
    _ortam_hazirla()
    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    import torch
    n = torch.cuda.device_count()
    name = torch.cuda.get_device_name(0) if n else None
    if n != 1 or name != BEKLENEN_GPU:
        raise SystemExit(f"GPU DUR: {n} kart, {name} (CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 gerekli)")
    import run_frontend
    from turb.runner import Cell, _build_cfg
    hd = CIKTI / etiket / ad / "onuc"
    hd.mkdir(parents=True, exist_ok=True)
    n_kare = sum(1 for l in open(Path(kok) / ad / "mav0/cam0/data.csv") if not l.startswith("#"))
    cell = Cell(dataset="airsim", base_path=kok, sequence=ad, seed=int(tohum), D_over_r0=0.0, windows=(0,),
                window_len=n_kare, sensor_noise_dn=0.0)
    cfg = _build_cfg(cell, str(hd))
    t0 = time.time()
    st = run_frontend.main(cfg)
    return {"dump": str(hd / f"landmarks_{ad}.pkl"), "t_frontend": time.time() - t0, "fe_gpu": name,
            "frontend_stats": {k: float(v) for k, v in (st or {}).items() if isinstance(v, (int, float))}}


# ----------------------------------------------------------------------------- koşum
def _yapilan(jsonl: Path):
    d = set()
    if jsonl.exists():
        for l in open(jsonl):
            try:
                r = json.loads(l)
            except json.JSONDecodeError:
                continue
            if r.get("status") == "ok":
                a = r.get("airsim", {})
                d.add((r.get("tur"), a.get("turb_sekans"), a.get("kol"), a.get("k")))
    return d


def kos(isler, jsonl: Path, isci: int):
    """isler: [("kol", spec, kol, k) | ("artik", spec)]; kayıtlar ana süreçte yazılır, devam ettirilebilir."""
    import multiprocessing as mp
    jsonl.parent.mkdir(parents=True, exist_ok=True)
    (jsonl.parent / ".trackerignore").touch()
    etiket = jsonl.parent.name
    d = _yapilan(jsonl)
    todo = []
    for j in isler:
        if j[0] == "kol":
            if ("kol", j[1]["turb_sekans"], j[2], float(j[3])) not in d:
                todo.append(("kol", (j[1], j[2], float(j[3]), etiket)))
        elif ("artik", j[1]["turb_sekans"], None, None) not in d:
            todo.append(("artik", (j[1], etiket)))
    # yavaş kollar önce (A3), artık işleri sonda
    sira = {"A3": 0, "A0C+A3": 1, "AG": 2, "A0": 3, "A0C": 4}
    todo.sort(key=lambda t: (t[0] == "artik", sira.get(t[1][1], 9) if t[0] == "kol" else 0))
    print(f"[kollar] {len(isler)} iş, {len(todo)} koşulacak, işçi {isci} -> {jsonl}", flush=True)
    if not todo:
        return
    t0 = time.time()
    with mp.get_context("spawn").Pool(isci, maxtasksperchild=4) as pool, open(jsonl, "a") as fout:
        res = pool.imap_unordered(_disp, todo)
        for i, rec in enumerate(res, 1):
            fout.write(json.dumps(rec, default=_js) + "\n"); fout.flush()
            a = rec.get("airsim", {})
            print(f"[kollar] {i}/{len(todo)} {rec.get('tur')} {a.get('turb_sekans')} {a.get('kol', '')} "
                  f"{rec.get('status')} {a.get('sure_s', 0):.0f}s ATE={rec.get('a_ate_se3')} "
                  f"NEESp={rec.get('n_nees_pos_mean')} bit={rec.get('a0_bit_ayni')} "
                  f"geçen {(time.time() - t0) / 60:.1f} dk", flush=True)


def _disp(t):
    try:
        return kol_isi(t[1]) if t[0] == "kol" else artik_isi(t[1])
    except BaseException as e:                                       # noqa: BLE001
        import traceback
        sp = t[1][0]
        return {"tur": t[0], "status": "fail", "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()[-1500:],
                "airsim": dict(sp, kol=(t[1][1] if t[0] == "kol" else None),
                               k=(t[1][2] if t[0] == "kol" else None))}


def _js(o):
    import numpy as np
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, tuple)):
        return list(o)
    return str(o)


def denetim_ortam(gpu_gerekli=False):
    """Kapılar: SLAM_CONFIG yok, CUDA ortamı GPU1, code_fp beklenen."""
    env = {k: os.environ.get(k) for k in ("CUDA_DEVICE_ORDER", "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS",
                                          "SLAM_CONFIG", "SLAM_DATA_ROOT")}
    if env["SLAM_CONFIG"]:
        sys.exit("DUR: SLAM_CONFIG tanımlı")
    if env["CUDA_DEVICE_ORDER"] != "PCI_BUS_ID" or env["CUDA_VISIBLE_DEVICES"] != "1":
        sys.exit("DUR: CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 gerekli (yalnız RTX 3090)")
    sys.path.insert(0, str(PE))
    cwd = os.getcwd(); os.chdir(PE)
    from turb.runner import code_fingerprint
    fp = code_fingerprint(); os.chdir(cwd)
    if fp != BEKLENEN_CODE_FP:
        sys.exit(f"DUR: code_fp {fp} != {BEKLENEN_CODE_FP}")
    return {"env": env, "code_fp": fp, "kod": kod_ozeti()}


def izgara_isleri(k, bekle=False, sekans=IZGARA_SEKANS, c0s=IZGARA_C0, tohumlar=IZGARA_TOHUM, kollar=KOLLAR):
    isler, eksik = [], []
    for pen in PENCERELER:
        for th in tohumlar:
            for c in c0s:
                sp = spec_ana(sekans, c, th, pen)
                if sp is None:
                    eksik.append((sekans, c, th, pen)); continue
                isler += [("kol", sp, kol, k) for kol in kollar]
                if c in (0.0, 1e-11):
                    isler.append(("artik", sp))
    return isler, eksik


def kalibre_specleri():
    """Temiz F5 pilot tohumları 0-2: W1 dökümleri GPU1 (f5_cift ATE eşleşmesi, ON_TEST_KOLLAR §3);
    W2'nin üst-düzey dökümleri GPU0 tarafından ezilmiş -> GPU1'de yeniden üretilir (spec 'uret')."""
    out = []
    for th in K_TOHUM:
        ad = turb_adi("AF_F5_h5", 0.0, th, "20:857")
        out.append(spec_dizin("AF_F5_h5", 0.0, th, "20:857", HERE / "cikti_izgara" / ad, kaynak="pilot_gpu1"))
        out.append({"sekans": "AF_F5_h5", "C0": 0.0, "tohum": th, "pencere": "1157:1994", "uret": True})
    return out


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="kip", required=True)
    h = sub.add_parser("hucre"); u = sub.add_parser("uret")
    for p in (h, u):
        p.add_argument("--sekans", required=True); p.add_argument("--C0", type=float, required=True)
        p.add_argument("--tohum", type=int, required=True); p.add_argument("--pencere", required=True)
        p.add_argument("--kollar", nargs="+", default=list(KOLLAR)); p.add_argument("--k", type=float, default=None)
        p.add_argument("--etiket", required=True); p.add_argument("--isci", type=int, default=1)
        p.add_argument("--artik", action="store_true")
        p.add_argument("--n_kare", type=int, default=None, help="duman: arka-uç yalnız ilk n kare")
    h.add_argument("--dizin", default=None, help="ön-uç dökümünün dizini (yoksa ana ızgara kaydı aranır)")
    u.add_argument("--sakla", action="store_true"); u.add_argument("--uret_isci", type=int, default=6)
    g = sub.add_parser("izgara"); g.add_argument("--isci", type=int, default=3); g.add_argument("--bekle", action="store_true")
    g.add_argument("--kuru", action="store_true"); g.add_argument("--etiket", default="izgara")
    g.add_argument("--kollar", nargs="+", default=list(KOLLAR))
    c = sub.add_parser("kalibre"); c.add_argument("--isci", type=int, default=3); c.add_argument("--kuru", action="store_true")
    c.add_argument("--uret_isci", type=int, default=6)
    c.add_argument("--k_izgara", nargs="+", type=float, default=list(K_IZGARA),
                   help="k adayları (varsayılan EuRoC ızgarası; temiz AirSim A0 NEES'i χ²₃ altında kalırsa k<1 eklenmesi koordinatör kararıdır)")
    a = ap.parse_args()
    man = denetim_ortam()
    print(json.dumps(man, ensure_ascii=False), flush=True)

    if a.kip in ("hucre", "uret"):
        k = a.k if a.k is not None else (json.load(open(K_DOSYA))["k"] if K_DOSYA.exists() else None)
        if k is None and any(x in ("A0C", "A0C+A3") for x in a.kollar):
            sys.exit("k yok: --k ver ya da önce kalibre")
        k = 1.0 if k is None else k
        if a.kip == "uret":
            sp = uret({"sekans": a.sekans, "C0": a.C0, "tohum": a.tohum, "pencere": a.pencere}, a.etiket,
                      isci=a.uret_isci, sakla=a.sakla)
            print(json.dumps({x: sp.get(x) for x in ("turb_sekans", "dump", "t_frontend", "uretim_s", "fe_gpu",
                                                     "gt_bayt_ayni")}), flush=True)
            ana = spec_ana(a.sekans, a.C0, a.tohum, a.pencere)
            if ana:
                sp["ana_traj"], sp["ana_ate"], sp["ana_dump"] = ana["ana_traj"], ana["ana_ate"], ana["dump"]
            json.dump(sp, open(CIKTI / a.etiket / f"{sp['turb_sekans']}_spec.json", "w"), indent=1, default=_js)
        else:
            sp = (spec_dizin(a.sekans, a.C0, a.tohum, a.pencere, a.dizin) if a.dizin
                  else spec_ana(a.sekans, a.C0, a.tohum, a.pencere))
            if sp is None:
                sys.exit("hücre bulunamadı")
        if a.n_kare:
            sp["n_kare_sinir"] = int(a.n_kare)
        isler = [("kol", sp, kol, k) for kol in a.kollar] + ([("artik", sp)] if a.artik else [])
        kos(isler, CIKTI / a.etiket / "kayit.jsonl", a.isci)

    elif a.kip == "izgara":
        k_gerek = any(x in ("A0C", "A0C+A3") for x in a.kollar)
        if k_gerek and not K_DOSYA.exists():
            sys.exit(f"DUR: {K_DOSYA} yok — önce `kol_kos.py kalibre`")
        k = float(json.load(open(K_DOSYA))["k"]) if k_gerek else 1.0   # A0/AG/A3: tek-hücre kipindeki gibi k=1 (cal_factor etkisiz)
        while True:
            isler, eksik = izgara_isleri(k, kollar=tuple(a.kollar))
            print(f"[izgara] kollar={a.kollar} k={k}; hazır hücre {len({j[1]['turb_sekans'] for j in isler})}, "
                  f"ana ızgarada henüz yok {len(eksik)}", flush=True)
            if a.kuru:
                print("[kuru] eksik:", eksik); return
            kos(isler, CIKTI / a.etiket / "kayit.jsonl", a.isci)
            if not eksik or not a.bekle:
                break
            time.sleep(600)
        print("[izgara] BITTI" + ("" if not eksik else f" (eksik {len(eksik)} hücre)"), flush=True)

    elif a.kip == "kalibre":
        specs = kalibre_specleri()
        print(f"[kalibre] {len(specs)} temiz pencere; k ∈ {a.k_izgara}; GPU1 yeniden üretim: "
              f"{sum(1 for s in specs if s.get('uret'))}", flush=True)
        if a.kuru:
            return
        hazir = []
        for s in specs:
            if s.get("uret"):
                sf = CIKTI / "kalibre" / f"{turb_adi(s['sekans'], 0.0, s['tohum'], s['pencere'])}_spec.json"
                if sf.exists():
                    s = json.load(open(sf))
                else:
                    s = uret(s, "kalibre", isci=a.uret_isci)
                    (CIKTI / "kalibre").mkdir(parents=True, exist_ok=True)
                    json.dump(s, open(sf, "w"), indent=1, default=_js)
            hazir.append(s)
        isler = [("kol", s, "A0", 1.0) for s in hazir] + [("kol", s, "A0C", float(kk)) for s in hazir for kk in a.k_izgara]
        kos(isler, CIKTI / "kalibre" / "kayit.jsonl", a.isci)
        r = subprocess.run([PY, str(HERE / "analiz_kollar.py"), "kalibre"], capture_output=True, text=True)
        print(r.stdout[-3000:], r.stderr[-2000:], flush=True)


if __name__ == "__main__":
    main()
