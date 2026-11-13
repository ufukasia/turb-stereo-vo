"""
Deney orkestrasyonu — Monte-Carlo hücreleri, paralel koşum, kanıt disiplini.

HÜCRE = (veri kümesi, sekans, türbülans parametreleri, tohum). Bir hücre için:
    1. ÖN-UÇ bir kez koşar (GPU) -> geçici landmark dump'ı
    2. Aynı dump üzerinde TÜM pencereler x TÜM kollar arka-uçta çözülür (CPU)
    3. Metrikler JSONL'e yazılır, dump SİLİNİR
Bu, rapordaki R7 azaltmasının uygulamasıdır: pahalı ön-uç tek sefer, ucuz arka-uç
kolları defalarca.

PYECO-SLAM BAĞLANTISI (pyeco_turb, 2026-11-07): ön-uç `run_frontend.main(cfg)`, arka-uç
`gtsam_backend.main(config.backend_args(cfg, argv))` — PyEco'nun kendi giriş noktaları. Arka-uç
kol başına bir kez sekansın başından `max(windows)+window_len` kareye koşar; değerlendirme pencere
başına yapılır. BİRİNCİL yörünge PyEco makalesinin ölçütü olan GERÇEK ZAMANLI (nedensel) yörüngedir
(`PYECO_NEDENSEL_OUT`; anahtar kare kipinde ara karelerin izleme pozu), `m_*` alanları ondandır;
düzleştirilmiş yörünge `ms_*` alanlarındadır. NEES (`n_*`) pencere-göreli kovaryans dökümünden.
PENCERE PROTOKOLÜ: 400 kare, pencereler 500/900/1300/1600 -> [1300,1700) ile [1600,2000) 100 kare
ÖRTÜŞÜR (denetim B6: sözde-tekrar). `warmup` alanı yalnız kayıt içindir: arka-uç her zaman kare 0'dan
başlar (ortak ısınma geçmişi), ayrı bir ısınma kesimi uygulanmaz.
FAZ 2 (2026-11-07): pencereler sekans başına GT hareketiyle seçilir (`turb/pencere_sec.py`,
`evidence/pencereler_faz2.json`) ve `windows_for(seq, path)` ile hücreye verilir; `path`
verilmezse eski WINDOW_STARTS aynen kullanılır.

TEKRARLANABİLİRLİK: her kayıt hücrenin tam parametre sözlüğünü ve kod sürümünü
taşır; koşum yarıda kesilirse `--resume` ile kaldığı yerden devam eder (tamamlanmış
hücreler JSONL'den okunur).

PARALELLİK: ön-uç GPU semaforuyla sınırlıdır (VRAM), arka-uç CPU çekirdeklerine
yayılır. Her işçi kendi sürecinde çalışır; bir hücrenin çökmesi koşumu düşürmez.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
import time
import traceback
from dataclasses import dataclass, asdict, field

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Rapor §6.0'dan devralınan protokol
WINDOW_STARTS = (500, 900, 1300, 1600)
WINDOW_LEN = 400
WARMUP = 100


def windows_for(sequence: str, path: str | None = None, default=WINDOW_STARTS) -> tuple:
    """Sekans başına pencere başlangıçları (Faz 2, 2026-11-07).

    `path` verilmezse eski davranış: `default` (= WINDOW_STARTS) döner. Verilirse
    `turb/pencere_sec.py` çıktısı okunur ({"sequences": {seq: {"starts": [...]}}} ya da
    düz {seq: [...]}); sekans tabloda yoksa sessizce varsayılana DÜŞMEZ, hata verir.
    Hücre (`Cell.windows`) zaten pencere listesini taşıdığı için koşum yolu değişmez.
    """
    if not path:
        return tuple(int(x) for x in default)
    with open(path) as f:
        d = json.load(f)
    seqs = d.get("sequences", d)
    if sequence not in seqs:
        raise KeyError(f"{sequence} pencere tablosunda yok: {path}")
    v = seqs[sequence]
    starts = v["starts"] if isinstance(v, dict) else v
    if not starts:
        raise ValueError(f"{sequence} için seçilmiş pencere yok: {path}")
    return tuple(int(x) for x in starts)


@dataclass
class Cell:
    """Tek Monte-Carlo hücresi."""
    dataset: str                      # 'euroc' | 'tumvi' | 'tartanair' | 'kitti' (yalnizca etiket)
    base_path: str
    sequence: str
    seed: int = 0
    D_over_r0: float = 0.0
    theta_c_over_fov: float = 0.15
    tau0_over_dt: float = 1.0
    blur: str = "speckle_global"
    tilt: bool = True
    scint_sigma_I2: float = 0.0
    arms: tuple = ("A0",)
    windows: tuple = WINDOW_STARTS
    window_len: int = WINDOW_LEN
    warmup: int = WARMUP
    # Ön-uç tohum/şiddet dışındaki türbülans eksenleri ve kol parametreleri (config.yaml'dan farklıysa)
    psf_legacy_scale: bool = False
    nees_rel: str = "rel"          # NEES göreli pozu: rel (kovaryansla aynı kestirim) | rel_final (eski tanım)
    # E2 ile kalibre edilmis ampirik korelasyon VARSAYILANDIR. (2026-08-20: bu
    # varsayilan "physical" idi ve runner onu CLI'ya ACIKCA geciyordu, yani
    # backend'deki "empirical" varsayilani hicbir kosumda etkili olmuyordu.)
    corr_model: str = "empirical"
    frame_dt: float = 0.05
    px_per_rad: float = 170648.0
    image_width: int = 752
    stereo_baseline: float = 0.11
    nees: bool = False
    cov_every: int = 5
    # Loop closure / PGO / GBA. VARSAYILAN KAPALI — ana deneylerde kasıtlı olarak
    # kapalıdır: loop closure, ölçüm modelinden BAĞIMSIZ bir bilgi kaynağı (görünüm
    # tabanlı yer tanıma) enjekte eder ve yanlış ağırlıklandırmadan doğan sürüklenmeyi
    # kapanış anında düzeltir; bu, kollar arası farkı maskeleyip gücü düşürür. Ayrıca
    # düzeltme fixed-lag kovaryansının GÖRMEDİĞİ bir işlemdir, dolayısıyla NEES ile
    # ATE farklı kestiricilere ait olur. E9 deneyi bunu AYRI bir eksen olarak sınar.
    loop_closure: bool = False
    tag: str = ""
    restoration: str = "none"      # E3 kolu: none|temporal_median|lucky|elastic
    # Veri kümesine özgü ayarlar (ör. TartanAir: hızlı hareket, uzak sahne).
    # cfg_overrides: config.yaml sözlüğüne birleştirilir; backend_argv: ek CLI bayrakları.
    cache_dump: bool = False       # ON-UC dump'i kalici onbellege yaz/oku
    cfg_overrides: tuple = ()      # (("frontend","max_landmark_depth",150.0), ...)
    backend_argv: tuple = ()       # ("--jump_max_trans","8.0", ...)
    # Faz 2 (2026-11-08): sensör gürültüsü [DN], türbülanstan sonra, D/r0=0'da da etkin.
    # noise_seed None -> türbülans tohumu (seed). 0 -> kapalı: anahtar/önbellek/kanca eskisiyle aynı.
    sensor_noise_dn: float = 0.0
    noise_seed: int | None = None

    def key(self) -> str:
        return (f"{self.dataset}/{self.sequence}/dr{self.D_over_r0:g}"
                f"/tc{self.theta_c_over_fov:g}/t0{self.tau0_over_dt:g}"
                f"/blur{self.blur}/tilt{int(self.tilt)}/sc{self.scint_sigma_I2:g}"
                f"/res{self.restoration}"
                f"/s{self.seed}"
                + (f"/n{self.sensor_noise_dn:g}"
                   + (f"_ns{self.noise_seed}" if self.noise_seed is not None else "")
                   if self.sensor_noise_dn > 0 else "")
                + f"{('/' + self.tag) if self.tag else ''}")


# =============================================================================
#  ON-UC DUMP ONBELLEGI
# =============================================================================
# GPU tabani: her hucre 2000 karelik bir ALIKED+LightGlue kosumu ister (~40 s saf
# GPU). Deneylerin cogu AYNI on-uc kosumunu paylasir: E4 (dr 3/5/8 x tohum<6),
# E3'un "none" kolu ve E8'in ana on-uc varyanti, E1a hucrelerinin ALT KUMESIDIR.
# Bu yuzden dump'lar, ON-UCU BELIRLEYEN parametrelerin hash'i ile kalici bir
# onbellege yazilir ve sonraki deneyler GPU'ya hic ugramadan onu okur
# (hesaplama-kaynagi.md §6: "agin kaydedilmis ciktisi okunuyorsa CPU" — mesru).

def _frontend_config_fingerprint() -> str:
    """config.yaml'ın ÖN-UÇ ÇIKTISINI belirleyen bölümünün parmak izi.

    NEDEN GEREKLİ (2026-08-17'de ölçülerek eklendi): dump önbelleği yalnız `Cell`
    alanlarıyla anahtarlanıyordu. SLAM yazılımı güncellenip config.yaml'daki ön-uç
    parametreleri değiştiğinde (yeni 8-sekanslı Tabu koşumu: max_keypoints 268->256,
    detection_threshold 0.004456->0.0045, max_landmark_depth 45.6->49.0,
    temporal_ransac_threshold 2.45->2.0, ...) önbellek SESSİZCE ESKİ dump'ları
    servis ediyordu — yani "yeni SLAM ile koşuldu" denilen deneyler eski ön-uç
    çıktısını kullanmış olurdu. Parmak izi anahtara girince böyle bir değişiklik
    önbelleği ISKALAR ve dump yeniden üretilir; bu, sessiz bozulmayı imkânsız kılar.

    `visualize`/`save_visualizations`/`viz_interval` gibi yalnız görselleştirmeyi
    etkileyen anahtarlar KASITLA dışarıda: bunlar geometriyi değiştirmez ve
    dahil edilseydi gereksiz yeniden koşum üretirlerdi.
    """
    import hashlib
    import config as _cfg
    try:
        f = dict((_cfg.load() or {}).get("frontend", {}) or {})
    except BaseException:                       # config okunamıyorsa önbelleği zorlama
        return "nofe"
    for k in ("visualize", "viz_interval", "save_visualizations", "export_keypoints"):
        f.pop(k, None)
    return hashlib.sha1(repr(sorted(f.items())).encode()).hexdigest()[:8]


_FE_CFG_FP = None
_CFG_PSF_LEGACY = None


def _psf_legacy(cell: "Cell") -> bool:
    """Hucrede etkin turbulence.psf_legacy_scale (config + cfg_overrides)."""
    global _CFG_PSF_LEGACY
    if _CFG_PSF_LEGACY is None:
        import config as _cfg
        try:
            t = (_cfg.load() or {}).get("turbulence", {}) or {}
            _CFG_PSF_LEGACY = bool(t.get("psf_legacy_scale", False))
        except BaseException:
            _CFG_PSF_LEGACY = False
    val = _CFG_PSF_LEGACY
    for sec, key, v in (cell.cfg_overrides or ()):
        if sec == "turbulence" and key == "psf_legacy_scale":
            val = bool(v)
    return val


def _fe_key(cell: "Cell") -> str:
    import hashlib
    global _FE_CFG_FP
    if _FE_CFG_FP is None:
        _FE_CFG_FP = _frontend_config_fingerprint()
    parts = (cell.dataset, cell.base_path, cell.sequence, cell.seed, cell.D_over_r0,
             cell.theta_c_over_fov, cell.tau0_over_dt, cell.blur, cell.tilt,
             cell.scint_sigma_I2, cell.restoration, cell.frame_dt, cell.px_per_rad,
             cell.stereo_baseline, tuple(cell.cfg_overrides or ()),
             int(max(cell.windows) + cell.window_len), _FE_CFG_FP)
    # Piksel olcekli PSF (2026-11-07) bulanikli on-uc ciktisini degistirir; eski
    # kip anahtarlari (ve arsivdeki dump'lar) degismesin diye isaret yalniz yeni
    # kipte ve bulaniklik gercekten uygulandiginda eklenir.
    if (str(cell.blur).startswith("speckle") and float(cell.D_over_r0) > 0
            and not _psf_legacy(cell)):
        parts = parts + ("psf_px_scale",)
    # Sensör gürültüsü (Faz 2): yalnız açıkken anahtara girer (eski anahtarlar değişmez)
    if float(cell.sensor_noise_dn) > 0:
        parts = parts + ("noise", float(cell.sensor_noise_dn),
                         int(cell.seed if cell.noise_seed is None else cell.noise_seed))
    h = hashlib.sha1(repr(parts).encode()).hexdigest()[:16]
    return f"{cell.dataset}_{cell.sequence}_dr{cell.D_over_r0:g}_s{cell.seed}_{h}"


def _cache_dir() -> str:
    d = os.environ.get("TURB_DUMP_CACHE", "/home/masa24/Desktop/fable/turb_dump_cache")
    os.makedirs(d, exist_ok=True)
    return d


def _gpu_name() -> str | None:
    """Ön-ucun koştuğu GPU'nun adı (2026-11-08). PyEco ön-ucu aynı GPU'da bit-belirlenimci, farklı
    GPU modelleri arasında DEĞİL; kayıt hangi kartta üretildiğini taşımalı (aynı aile tek GPU'da)."""
    try:
        import torch as _t
        if _t.cuda.is_available():
            return str(_t.cuda.get_device_name(0))
    except Exception:
        pass
    return None


def _build_cfg(cell: Cell, out_dir: str):
    import config
    cfg = config.load(str(config.DEFAULT_PATH))
    cfg["dataset"].update(base_path=cell.base_path, sequence=cell.sequence,
                          max_frames=int(max(cell.windows) + cell.window_len))
    cfg["frontend"]["visualize"] = False
    cfg["frontend"]["save_visualizations"] = False
    cfg["paths"]["outputs"] = out_dir
    # Loop closure: kapanış tespiti + PGO açılır, GBA AÇILMAZ.
    # ÖLÇÜM (2026-08-17, evidence/gba_maliyeti_*.json, MH01+MH02, boş makine):
    #   loop kapalı -> PGO : ATE %33 / %38 iyileşiyor, süre farkı ~0 (bedava)
    #   PGO -> PGO+GBA     : süre +%16, ATE MH01'de %85 KÖTÜLEŞİYOR, MH02'de aynı
    # GBA en iyi ihtimalle nötr, en kötü ihtimalle zararlı; bu yüzden hiçbir kolda
    # açılmaz. Aday mekanizmalar ve önerilen düzeltmeler: GBA_BULGUSU.md.
    cfg["loop"]["enable"] = bool(cell.loop_closure)
    cfg["loop"]["pgo"] = bool(cell.loop_closure)
    cfg["loop"]["gba"] = False
    cfg.setdefault("turbulence", {})
    cfg["turbulence"].update(
        enabled=bool(cell.D_over_r0 > 0),
        psf_legacy_scale=bool(cell.psf_legacy_scale),
        seed=int(cell.seed), D_over_r0=float(cell.D_over_r0),
        theta_c_over_fov=float(cell.theta_c_over_fov),
        tau0_over_dt=float(cell.tau0_over_dt),
        blur=str(cell.blur), tilt=bool(cell.tilt),
        scint_sigma_I2=float(cell.scint_sigma_I2),
        frame_dt=float(cell.frame_dt), px_per_rad=float(cell.px_per_rad),
        stereo_baseline=float(cell.stereo_baseline),
        restoration=str(cell.restoration))
    if float(cell.sensor_noise_dn) > 0:      # yalnız açıkken yazılır: kapalıyken cfg sözlüğü eskisiyle aynı
        cfg["turbulence"].update(sensor_noise_dn=float(cell.sensor_noise_dn),
                                 noise_seed=(None if cell.noise_seed is None else int(cell.noise_seed)))
    # Veri kümesine özgü ön-uç/arka-uç ayarları (TartanAir: hızlı hareket, uzak sahne).
    for sec, key, val in (cell.cfg_overrides or ()):
        cfg[sec][key] = val
    return cfg



# =============================================================================
# KOD KÖKENİ (provenance) — hangi SLAM sürümüyle üretildiği her kayda yazılır.
#
# NEDEN (2026-08-17): SLAM arka-ucu köklü değiştiğinde eski kayıtlar diskte durur.
# Kayıtta köken damgası yoksa iki sürümün sonuçları AYNI tabloda birleşir ve fark,
# yöntemin etkisi sanılır — makaleyi çökertecek türden sessiz bir bozulma. Damga
# varsa `turb.sanity` karışımı KIRMIZI bulgu olarak yakalar.
# =============================================================================
_CODE_FP = None


def code_fingerprint() -> str:
    """SLAM + turb kaynak dosyalarının içerik özeti."""
    global _CODE_FP
    if _CODE_FP is not None:
        return _CODE_FP
    import hashlib
    h = hashlib.sha1()
    # PyEco-SLAM dosya kümesi (pyeco_turb): ön-uç + ALIKED/LightGlue kodu, arka-uç + kayan pencere +
    # eleme, giriş noktaları, yapılandırma; turb/ tarafında bozulma, kollar, kestiriciler, metrikler
    # ve bu koşucunun kendisi (pencere tanımı sonucu belirler).
    files = ["stereo_front_end.py", "gtsam_backend.py", "config.py", "config.yaml",
             "kayan_pencere.py", "landmark_culling.py", "loop_closing.py",
             "run_frontend.py", "run_backend.py",
             "lightglue/__init__.py", "lightglue/aliked.py", "lightglue/lightglue.py",
             "lightglue/utils.py",
             "turb/__init__.py", "turb/backend.py", "turb/covariance.py", "turb/factors.py",
             "turb/simulator.py", "turb/fields.py", "turb/screen.py", "turb/optics.py",
             "turb/correlation.py", "turb/hooks.py", "turb/restoration.py", "turb/metrics.py",
             "turb/dimm.py", "turb/dimm_v2.py", "turb/runner.py"]
    # "config.yaml" girdisi etkin yapilandirma dosyasini ozetler ($SLAM_CONFIG ile
    # secilmisse o dosyayi); etiket ayni kalir, varsayilan durumda ozet birebir eskisi.
    sys.path.insert(0, ROOT)
    import config as _cfgmod
    for fn in files:
        pth = str(_cfgmod.DEFAULT_PATH) if fn == "config.yaml" else os.path.join(ROOT, fn)
        try:
            with open(pth, "rb") as fh:
                h.update(fn.encode())
                h.update(hashlib.sha1(fh.read()).digest())
        except OSError:
            h.update(b"MISSING:" + fn.encode())
    _CODE_FP = h.hexdigest()[:12]
    return _CODE_FP



def _arm_argv(cell: "Cell", arm: str, out: str, cov: str, st: str):
    argv = ["--out", out, "--start_frame", "0",
            "--max_frames", str(int(max(cell.windows) + cell.window_len)),
            "--turb_arm", arm, "--turb_corr_model", cell.corr_model,
            "--turb_D_over_r0", str(cell.D_over_r0),
            "--turb_theta_c_over_fov", str(cell.theta_c_over_fov),
            "--turb_px_per_rad", str(cell.px_per_rad),
            "--turb_image_width", str(cell.image_width),
            "--turb_stats_out", st] + list(cell.backend_argv or ())
    if cov:
        argv += ["--cov_out", cov, "--cov_every", str(cell.cov_every)]
    return argv


def _run_arms(cell: "Cell", cfg: dict, d: str, gt: str, f_ts, base_rec: dict, verbose=False):
    """Hücrenin bütün kollarını PyEco arka-ucunda koş, pencere kayıtlarını döndür."""
    import config
    from gtsam_backend import main as backend_main
    from turb.metrics import trajectory_metrics, nees_from_dump
    recs = []
    for arm in cell.arms:
        out = os.path.join(d, f"traj_{arm}.csv")
        rt = os.path.join(d, f"traj_{arm}_realtime.csv")
        cov = os.path.join(d, f"cov_{arm}.npz") if cell.nees else ""
        st = os.path.join(d, f"stats_{arm}.json")
        argv = _arm_argv(cell, arm, out, cov, st)
        tb = time.time()
        err = None
        _old_env = os.environ.get("PYECO_NEDENSEL_OUT")
        os.environ["PYECO_NEDENSEL_OUT"] = rt        # gercek zamanli (nedensel) yorunge = birincil olcut
        try:
            backend_main(config.backend_args(cfg, argv=argv))
        except BaseException as e:
            s_ = str(e).strip().splitlines()
            err = f"{type(e).__name__}: {s_[0] if s_ else ''}"
            if verbose:
                traceback.print_exc()
        finally:
            if _old_env is None:
                os.environ.pop("PYECO_NEDENSEL_OUT", None)
            else:
                os.environ["PYECO_NEDENSEL_OUT"] = _old_env
        t_be = time.time() - tb
        arm_stats = None
        if os.path.exists(st):
            try:
                with open(st) as fjs:
                    arm_stats = json.load(fjs)
            except json.JSONDecodeError:
                pass
        for w in cell.windows:
            i0 = int(w)
            i1 = int(min(w + cell.window_len, len(f_ts)) - 1)
            rec = dict(base_rec)
            rec.update({"cell": asdict(cell), "key": cell.key(), "window": int(w), "arm": arm,
                        "code_fp": code_fingerprint(), "t_backend": t_be, "arm_stats": arm_stats})
            if err is not None:
                rec.update(status="fail", error=err)
            elif i1 <= i0:
                rec.update(status="fail", error="pencere sekans dışında")
            else:
                t0w, t1w = float(f_ts[i0]), float(f_ts[i1])
                rec["t_win_start"], rec["t_win_end"] = t0w, t1w
                try:
                    m = trajectory_metrics(rt, gt, t_min=t0w, t_max=t1w)
                    rec.update({("m_" + k): v for k, v in m.items()})
                    ms = trajectory_metrics(out, gt, t_min=t0w, t_max=t1w)
                    rec.update({("ms_" + k): v for k, v in ms.items()})
                    if cov and os.path.exists(cov):
                        n = nees_from_dump(cov, gt, t_min=t0w, t_max=t1w, rel_key=cell.nees_rel)
                        rec.update({("n_" + k): v for k, v in n.items()})
                    rec["status"] = "ok" if m.get("valid") else "fail"
                except BaseException as e:
                    rec.update(status="fail", error=f"{type(e).__name__}: {e}")
            recs.append(rec)
    return recs

def run_cell(cell: Cell, gpu_sem=None, keep_dump: str | None = None,
             verbose: bool = False):
    """Bir hücreyi koş; kayıt listesi döndür."""
    sys.path.insert(0, ROOT)
    os.chdir(ROOT)
    import config
    import run_frontend
    from gtsam_backend import main as backend_main
    from turb.metrics import trajectory_metrics, nees_from_dump

    recs = []
    tmp = tempfile.mkdtemp(prefix="turbcell_", dir=os.environ.get("TURB_TMP", None))
    try:
        cfg = _build_cfg(cell, tmp)
        t0 = time.time()
        # GPU semaforu yalnız ON-UC suresince tutulur. Ayrica on-uc biter bitmez
        # CUDA onbellegi bosaltilir: aksi halde her isci sureci kendi ~1 GB'lik
        # ALIKED+LightGlue tahsisini arka-uc boyunca da tutar ve 10 iscide 8 GB'lik
        # kart tasar (olculdu: 718/725 hucre CUDA OOM ile dustu).
        stats_fe = None
        for _try in range(3):
            if gpu_sem is not None:
                gpu_sem.acquire()
            try:
                stats_fe = run_frontend.main(cfg)
                break
            except RuntimeError as _e:
                if "out of memory" not in str(_e).lower() or _try == 2:
                    raise
                time.sleep(20 * (_try + 1))
            finally:
                if gpu_sem is not None:
                    gpu_sem.release()
                try:
                    import torch as _t
                    _t.cuda.empty_cache()
                    _t.cuda.ipc_collect()
                except Exception:
                    pass
        t_fe = time.time() - t0
        dump = config.dump_path(cfg)
        gt = config.gt_path(cfg)
        if keep_dump:
            os.makedirs(keep_dump, exist_ok=True)
            shutil.copy(dump, os.path.join(keep_dump, f"{cell.key().replace('/','_')}.pkl"))

        # [E3] görüntü kalitesi özeti (restorasyon kolu açıksa)
        qual = None
        _qp = os.path.join(tmp, "turb_quality.json")
        if os.path.exists(_qp):
            try:
                with open(_qp) as _f:
                    qual = json.load(_f)
            except json.JSONDecodeError:
                pass

        # Pencere zaman damgaları: dump'tan kare indeksine göre okunur.
        import pickle
        with open(dump, "rb") as fpk:
            _d = pickle.load(fpk)
        f_ts = [int(fr_["timestamp"]) for fr_ in _d["frames"]]
        del _d

        # ARKA-UÇ KOL BAŞINA TEK KOŞUM: sekansın BAŞINDAN pencerelerin sonuna kadar (ortak ısınma
        # geçmişi; adaptif ön-uçta pencere ortasından başlatılan arka-uç landmark ilkleyemez).
        base = {"t_frontend": t_fe, "fe_gpu": _gpu_name(),
                "frontend_stats": {k: float(v) for k, v in (stats_fe or {}).items()
                                   if isinstance(v, (int, float))},
                "quality": qual}
        recs += _run_arms(cell, cfg, tmp, gt, f_ts, base, verbose=verbose)
    except BaseException as e:
        recs.append({"cell": asdict(cell), "key": cell.key(),
                     "code_fp": code_fingerprint(), "status": "cell_fail",
                     "error": f"{type(e).__name__}: {e}"})
        if verbose:
            traceback.print_exc()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return recs


# =============================================================================


# =============================================================================
#  İKİ-FAZLI HAT (hesaplama-kaynagi.md §5: GPU işleri seri, CPU işleri paralel)
# =============================================================================
# TEŞHİS (ölçüldü): tek işçinin içinde ön-uç (GPU) ve arka-ucu (CPU) ardışık
# koşturmak, işçilerin çoğunu GPU semaforunda BEKLETİR. 12 işçi / 2 GPU yuvasıyla
# ölçülen: ön-uç 43 s -> 198 s (GPU zaman-paylaşımı + CPU çekişmesi), CPU %50
# boşta, bağlam değişimi 768k/s, verim 0.60 hücre/dk.
#
# ÇÖZÜM: fazları AYIR ve çift-tamponla örtüştür.
#   FAZ A (GPU): `gpu_jobs` süreç yalnız ön-uç koşar, dump'ı diske yazar.
#   FAZ B (CPU): `cpu_jobs` süreç yalnız arka-uç kollarını koşar, dump'ı siler.
# Bir yığının FAZ B'si koşarken sonraki yığının FAZ A'sı asenkron başlatılır;
# böylece GPU ve CPU aynı anda doludur ve hiçbir iş yanlış donanıma alınmaz.

def _frontend_worker(args):
    """FAZ A — yalnız ön-uç (GPU). Dump'ı spool dizinine yazar, meta döndürür."""
    cell, spool = args
    sys.path.insert(0, ROOT)
    os.chdir(ROOT)
    # THREAD SINIRI (olculdu): her on-uc sureci OpenCV/torch havuzlariyla 28 thread
    # aciyordu; 5 escalisan on-uc + 8 arka-uc ile load 40 / 24 cekirdek (thrash).
    # OMP_NUM_THREADS bunu KAPSAMAZ — OpenCV kendi havuzunu cekirdek sayisina gore
    # kurar. Surec-ici paralellik burada YANLIS eksendir: paralellik zaten surecler
    # arasindadir (hesaplama-kaynagi.md §5).
    try:
        import cv2 as _cv2
        _cv2.setNumThreads(2)
    except Exception:
        pass
    try:
        import torch as _tt
        _tt.set_num_threads(1)
    except Exception:
        pass
    import config
    import run_frontend
    meta = {"cell": cell, "ok": False}
    d = os.path.join(spool, cell.key().replace("/", "_"))
    os.makedirs(d, exist_ok=True)
    t0 = time.time()
    try:
        cfg = _build_cfg(cell, d)
        # --- ONBELLEK OKUMA: dump varsa GPU'ya hic dokunma ---
        cpath = os.path.join(_cache_dir(), _fe_key(cell) + ".pkl")
        if cell.cache_dump and os.path.exists(cpath):
            shutil.copy(cpath, config.dump_path(cfg))
            qsrc = cpath[:-4] + "_quality.json"
            if os.path.exists(qsrc):
                shutil.copy(qsrc, os.path.join(d, "turb_quality.json"))
            meta.update(ok=True, dir=d, dump=config.dump_path(cfg),
                        gt=config.gt_path(cfg), t_frontend=0.0, cached=True,
                        stats={})
            return meta
        for _try in range(3):
            try:
                stats = run_frontend.main(cfg)
                break
            except RuntimeError as e:
                if "out of memory" not in str(e).lower() or _try == 2:
                    raise
                time.sleep(20 * (_try + 1))
            finally:
                try:
                    import torch as _t
                    _t.cuda.empty_cache()
                except Exception:
                    pass
        meta.update(ok=True, dir=d, dump=config.dump_path(cfg), gt=config.gt_path(cfg),
                    t_frontend=time.time() - t0, cached=False, gpu=_gpu_name(),
                    stats={k: float(v) for k, v in (stats or {}).items()
                           if isinstance(v, (int, float))})
        # --- ONBELLEGE YAZMA (hata olursa kosum DUSMEZ) ---
        if cell.cache_dump:
            try:
                tmp_c = cpath + f".tmp{os.getpid()}"
                shutil.copy(config.dump_path(cfg), tmp_c)
                os.replace(tmp_c, cpath)
                qp = os.path.join(d, "turb_quality.json")
                if os.path.exists(qp):
                    shutil.copy(qp, cpath[:-4] + "_quality.json")
            except OSError:
                pass
    except BaseException as e:
        meta["error"] = f"{type(e).__name__}: {e}"
        meta["t_frontend"] = time.time() - t0
        shutil.rmtree(d, ignore_errors=True)
    return meta


def _backend_worker(meta):
    """FAZ B — yalnız arka-uç kolları (CPU; GTSAM seyrek çözücü)."""
    sys.path.insert(0, ROOT)
    os.chdir(ROOT)
    try:
        import cv2 as _cv2
        _cv2.setNumThreads(1)
    except Exception:
        pass
    import config
    import pickle
    from gtsam_backend import main as backend_main
    from turb.metrics import trajectory_metrics, nees_from_dump

    cell = meta["cell"]
    recs = []
    if not meta.get("ok"):
        return [{"cell": asdict(cell), "key": cell.key(),
                 "code_fp": code_fingerprint(), "status": "cell_fail",
                 "error": meta.get("error", "ön-uç başarısız")}]
    d, dump, gt = meta["dir"], meta["dump"], meta["gt"]
    try:
        cfg = _build_cfg(cell, d)
        qual = None
        qp = os.path.join(d, "turb_quality.json")
        if os.path.exists(qp):
            try:
                with open(qp) as f:
                    qual = json.load(f)
            except json.JSONDecodeError:
                pass
        with open(dump, "rb") as f:
            f_ts = [int(fr["timestamp"]) for fr in pickle.load(f)["frames"]]

        base = {"t_frontend": meta.get("t_frontend"), "frontend_stats": meta.get("stats"),
                "fe_gpu": meta.get("gpu"),
                "quality": qual}
        recs += _run_arms(cell, cfg, d, gt, f_ts, base)
    except BaseException as e:
        recs.append({"cell": asdict(cell), "key": cell.key(),
                     "code_fp": code_fingerprint(), "status": "cell_fail",
                     "error": f"{type(e).__name__}: {e}"})
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return recs


def run_experiment_pipelined(cells, out_jsonl: str, cpu_jobs: int = 10,
                             gpu_jobs: int = 2, chunk: int = 24,
                             resume: bool = True, spool: str | None = None):
    """İki-faz + çift-tampon koşum. GPU ve CPU aynı anda dolu tutulur."""
    os.makedirs(os.path.dirname(out_jsonl) or ".", exist_ok=True)
    spool = spool or os.environ.get("TURB_SPOOL", "/tmp/turb_spool")
    os.makedirs(spool, exist_ok=True)

    done = set()
    if resume and os.path.exists(out_jsonl):
        with open(out_jsonl) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    done.add((r.get("key"), r.get("window"), r.get("arm")))
                except json.JSONDecodeError:
                    pass
    todo = [c for c in cells
            if not all((c.key(), int(w), a) in done
                       for w in c.windows for a in c.arms)]
    print(f"[runner] {len(cells)} hücre, {len(todo)} koşulacak "
          f"({len(cells)-len(todo)} zaten tamam), cpu={cpu_jobs} gpu={gpu_jobs} "
          f"yığın={chunk}", flush=True)
    if not todo:
        return

    ctx = mp.get_context("spawn")
    chunks = [todo[i:i + chunk] for i in range(0, len(todo), chunk)]
    t0 = time.time()
    n_done = 0
    with ctx.Pool(processes=int(gpu_jobs), maxtasksperchild=8) as gpu_pool, \
            ctx.Pool(processes=int(cpu_jobs), maxtasksperchild=8) as cpu_pool, \
            open(out_jsonl, "a") as fout:
        pending = gpu_pool.map_async(_frontend_worker,
                                     [(c, spool) for c in chunks[0]])
        for i, ch in enumerate(chunks):
            metas = pending.get()                       # FAZ A biter
            if i + 1 < len(chunks):                     # sonraki yığının FAZ A'sı
                pending = gpu_pool.map_async(           # ŞİMDİ başlar (örtüşme)
                    _frontend_worker, [(c, spool) for c in chunks[i + 1]])
            for recs in cpu_pool.imap_unordered(_backend_worker, metas):
                for r in recs:
                    fout.write(json.dumps(r, default=_js) + "\n")
                fout.flush()
                n_done += 1
                el = time.time() - t0
                rate = n_done / max(el / 60.0, 1e-9)
                eta = (len(todo) - n_done) / max(rate, 1e-9)
                print(f"[runner] {n_done}/{len(todo)}  {recs[0].get('key','?')}  "
                      f"({sum(1 for x in recs if x.get('status')=='ok')}/{len(recs)} ok)"
                      f"  {rate:.2f} hücre/dk  kalan ~{eta:.0f}dk", flush=True)


def _worker(args):
    cell, gpu_sem, keep_dump = args
    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    return run_cell(cell, gpu_sem, keep_dump)


def run_experiment(cells, out_jsonl: str, jobs: int = 8, gpu_jobs: int = 3,
                   resume: bool = True, keep_dump: str | None = None):
    """Hücre listesini paralel koş, sonuçları JSONL'e ekle."""
    os.makedirs(os.path.dirname(out_jsonl) or ".", exist_ok=True)
    done = set()
    if resume and os.path.exists(out_jsonl):
        with open(out_jsonl) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    done.add((r.get("key"), r.get("window"), r.get("arm")))
                except json.JSONDecodeError:
                    pass
    todo = []
    for c in cells:
        need = [(c.key(), int(w), a) for w in c.windows for a in c.arms]
        if all(k in done for k in need):
            continue
        todo.append(c)
    print(f"[runner] {len(cells)} hücre, {len(todo)} koşulacak "
          f"({len(cells)-len(todo)} zaten tamam), jobs={jobs} gpu={gpu_jobs}")
    if not todo:
        return

    ctx = mp.get_context("spawn")
    mgr = ctx.Manager()
    sem = mgr.BoundedSemaphore(int(gpu_jobs))
    t0 = time.time()
    n_done = 0
    # maxtasksperchild: her N hucrede isci sureci yenilenir -> CUDA baglami ve
    # bellek parcalanmasi sifirlanir (uzun kosumlarda OOM emniyeti).
    with ctx.Pool(processes=int(jobs), maxtasksperchild=4) as pool, \
            open(out_jsonl, "a") as fout:
        for recs in pool.imap_unordered(_worker, [(c, sem, keep_dump) for c in todo]):
            for r in recs:
                fout.write(json.dumps(r, default=_js) + "\n")
            fout.flush()
            n_done += 1
            el = time.time() - t0
            eta = el / n_done * (len(todo) - n_done)
            ok = sum(1 for r in recs if r.get("status") == "ok")
            print(f"[runner] {n_done}/{len(todo)}  {recs[0].get('key','?')}  "
                  f"({ok}/{len(recs)} ok)  geçen {el/60:.1f}dk  kalan ~{eta/60:.1f}dk",
                  flush=True)


def _js(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (set, tuple)):
        return list(o)
    return str(o)


def load_records(path: str):
    """JSONL kayıtlarını pandas DataFrame olarak yükle."""
    import pandas as pd
    rows = []
    with open(path) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            c = r.pop("cell", {})
            for k, v in c.items():
                r.setdefault("c_" + k, v)
            rows.append(r)
    return pd.DataFrame(rows)
