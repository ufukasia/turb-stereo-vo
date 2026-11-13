"""
uret_sekans.py sınamaları.

    python turb/test_uret_sekans.py            # hepsi (veri yoksa ilgili test atlanır)
    python turb/test_uret_sekans.py derinlik   # adında 'derinlik' geçenler

Testler:
  * kare aralığı ayrıştırma
  * derinlik büyütme (sabit, düzlem, basamak, sonsuz, gerçek AirSim karesi)
  * belirleyicilik: aynı parametrelerle iki üretim (+ çok süreçli üretim) bayt bayt aynı
  * tüketici rektifikasyonu: üretilen sekansı kancanın rektifikasyonu KİMLİK görür
  * AirSim: kancanın rektifikasyonu ham AirSim görüntüsünü bayt bayt değiştirmez
  * kanca eşdeğerliği: kod/stereo_front_end.py (SLAM_CONFIG=config_euroc.yaml,
    türbülans açık) ön-uca verdiği görüntü == üreticinin yazdığı PNG
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import cv2
import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from turb import uret_sekans as U  # noqa: E402

EUROC = Path("/data/sdata")
AIRSIM = Path("~/sim_envs/sdata_airsim/NH_egik_h5_v1").expanduser()
KOD = Path(__file__).resolve().parent.parent.parent / "kod"
PY = sys.executable
SONUC = {}


class Atla(Exception):
    pass


def _gerek(p: Path):
    if not p.exists():
        raise Atla(f"{p} yok")


def _dizin_ozeti(d: Path, haric=("turb_uretim.json",)) -> dict:
    out = {}
    for f in sorted(d.rglob("*")):
        if f.is_file() and f.name not in haric:
            out[str(f.relative_to(d))] = hashlib.sha256(f.read_bytes()).hexdigest()
    return out


# ----------------------------------------------------------------------------
def test_kareleri_ayristir():
    idx, ar = U.kareleri_ayristir("400:900", 3682)
    assert idx[0] == 400 and idx[-1] == 899 and len(idx) == 500 and ar == [(400, 900)]
    idx, ar = U.kareleri_ayristir("400:900,800:1300", 3682)
    assert len(idx) == 900 and ar == [(400, 900), (800, 1300)]
    idx, _ = U.kareleri_ayristir(None, 10)
    assert idx == list(range(10))
    idx, _ = U.kareleri_ayristir("5:", 10)
    assert idx == [5, 6, 7, 8, 9]
    idx, _ = U.kareleri_ayristir(":3,7", 10)
    assert idx == [0, 1, 2, 7]
    try:
        U.kareleri_ayristir("20:30", 10)
        raise AssertionError("sekans dışı aralık kabul edildi")
    except ValueError:
        pass


# ----------------------------------------------------------------------------
def test_derinlik_buyut_sabit_ve_boyut():
    d = np.full((192, 256), 7.25, np.float32)
    o = U.derinlik_buyut(d, (1024, 768))
    assert o.shape == (768, 1024) and o.dtype == np.float32
    assert np.allclose(o, 7.25, rtol=1e-6)


def test_derinlik_buyut_duzlem():
    """Eğik düzlem: 1/Z piksel koordinatlarında afin -> iç bölgede analitik değer."""
    h, w, H, W = 192, 256, 768, 1024
    a, b, c = 0.05, 2e-4, 3e-4
    v, u = np.mgrid[0:h, 0:w].astype(np.float64)
    d = (1.0 / (a + b * u + c * v)).astype(np.float32)
    o = U.derinlik_buyut(d, (W, H), kenar_orani=1.15)
    V, Uu = np.mgrid[0:H, 0:W].astype(np.float64)
    ul = (Uu + 0.5) * w / W - 0.5          # yüksek çözünürlük -> düşük koordinat
    vl = (V + 0.5) * h / H - 0.5
    ref = 1.0 / (a + b * ul + c * vl)
    ic = (slice(4, H - 4), slice(4, W - 4))
    rel = np.abs(o[ic] - ref[ic]) / ref[ic]
    SONUC["derinlik_duzlem_maks_bagil_hata"] = float(rel.max())
    assert rel.max() < 2e-4, rel.max()     # float32 + cv2 sabit-nokta ağırlıkları


def test_derinlik_buyut_basamak_ve_sonsuz():
    d = np.full((192, 256), 5.0, np.float32)
    d[:, 128:] = 40.0
    d[:20, :] = np.inf                     # gökyüzü
    o = U.derinlik_buyut(d, (1024, 768))
    vals = set(np.unique(o[np.isfinite(o)]).round(4).tolist())
    # Süreksizliklerde ara değer ÜRETİLMEZ; düzgün bölgeler sabit kalır.
    assert vals <= {5.0, 40.0}, sorted(vals)[:10]
    assert np.isinf(o[:79, :]).all() and np.isfinite(o[81:, :]).all()
    assert (o[400, :510] == 5.0).all() and (o[400, 514:] == 40.0).all()


def test_derinlik_ic_parametre():
    k = U.derinlik_ic_parametre([731.211779, 731.211779, 512.0, 384.0], (1024, 768),
                                (256, 192))
    assert abs(k[0] - 182.80294475) < 1e-6 and abs(k[2] - 127.625) < 1e-9 and abs(k[3] - 95.625) < 1e-9


def test_derinlik_airsim_gercek():
    _gerek(AIRSIM / "mav0/depth0/data")
    ok = U.DerinlikOkuyucu(AIRSIM)
    ts = U.asl_oku(AIRSIM)["ts"][100]
    d = ok.oku(ts)
    o = U.derinlik_buyut(d, (1024, 768))
    assert o.shape == (768, 1024)
    # Geri küçültme (4x4 blok, ters derinlik ortalaması) kaynağa yakın olmalı.
    inv_o = np.where(np.isfinite(o), 1.0 / o, 0.0)
    blok = inv_o.reshape(192, 4, 256, 4).mean(axis=(1, 3))
    inv_d = np.where(np.isfinite(d), 1.0 / d, 0.0)
    m = inv_d > 0
    rel = np.abs(blok[m] - inv_d[m]) / inv_d[m]
    SONUC["derinlik_airsim_geri_kucultme_medyan_bagil"] = float(np.median(rel))
    SONUC["derinlik_airsim_geri_kucultme_p99_bagil"] = float(np.percentile(rel, 99))
    assert np.median(rel) < 1e-3, np.median(rel)
    assert np.isinf(o).sum() == 16 * np.isinf(d).sum()


# ----------------------------------------------------------------------------
def _uret(girdi, kok, kareler, isci=1, D=5.0, ek=()):
    argv = ["--girdi", str(girdi), "--cikti_kok", str(kok), "--D_over_r0", str(D),
            "--tohum", "0", "--kareler", kareler, "--isci", str(isci), *ek]
    return U.main(argv)


def test_belirleyicilik():
    _gerek(EUROC / "MH01")
    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        a = _uret(EUROC / "MH01", t / "a", "400:412,890:900")
        b = _uret(EUROC / "MH01", t / "b", "400:412,890:900")
        c = _uret(EUROC / "MH01", t / "c", "400:412,890:900", isci=3)
        oa, ob, oc = _dizin_ozeti(a), _dizin_ozeti(b), _dizin_ozeti(c)
        SONUC["belirleyicilik_dosya_sayisi"] = len(oa)
        assert oa == ob, [k for k in oa if oa[k] != ob.get(k)][:5]
        assert oa == oc, [k for k in oa if oa[k] != oc.get(k)][:5]
        # farklı tohum -> bütün PNG'ler farklı (testin kör olmadığını göster)
        png = [k for k in oa if k.endswith(".png")]
        e = U.main(["--girdi", str(EUROC / "MH01"), "--cikti_kok", str(t / "e"),
                    "--D_over_r0", "5", "--tohum", "1", "--kareler", "400:412,890:900"])
        oe = _dizin_ozeti(e)
        assert sum(oa[k] != oe[k] for k in png) == len(png)


def test_tuketici_rektifikasyonu_kimlik():
    """Üretilen rektifiye sekansı, kancanın rektifikasyonu (T_cn_cnm1 okur) kimlik görür."""
    _gerek(EUROC / "MH01")
    with tempfile.TemporaryDirectory() as t:
        s = _uret(EUROC / "MH01", Path(t), "400:402")
        cal = U.asl_oku(s)
        r = U.Rektifikator(cal, "rektifiye")      # kancanın birebir aynısı
        (m1x, m1y), (m2x, m2y) = r.maps
        Uu, V = np.meshgrid(np.arange(752, dtype=np.float32), np.arange(480, dtype=np.float32))
        dev = max(np.abs(m1x - Uu).max(), np.abs(m1y - V).max(),
                  np.abs(m2x - Uu).max(), np.abs(m2y - V).max())
        SONUC["tuketici_rekt_harita_sapma_px"] = float(dev)
        f = cal["f0"][cal["ts"][0]]
        i0 = cv2.imread(str(s / "mav0/cam0/data" / f), cv2.IMREAD_GRAYSCALE)
        i1 = cv2.imread(str(s / "mav0/cam1/data" / f), cv2.IMREAD_GRAYSCALE)
        o0, o1 = r.uygula(i0, i1)
        SONUC["tuketici_rekt_piksel_farki"] = int(max(
            np.abs(o0.astype(int) - i0).max(), np.abs(o1.astype(int) - i1).max()))
        assert dev < 1e-3
        assert U.Rektifikator(cal, "otomatik").kip == "kimlik"
        assert abs(r.baseline - 0.110074) < 1e-5


def test_airsim_kanca_rektifikasyonu_kimlik():
    _gerek(AIRSIM / "mav0")
    cal = U.asl_oku(AIRSIM)
    assert U.Rektifikator(cal, "otomatik").kip == "kimlik"
    r = U.Rektifikator(cal, "rektifiye")
    ts = cal["ts"][0]
    i0 = cv2.imread(str(AIRSIM / "mav0/cam0/data" / cal["f0"][ts]), cv2.IMREAD_GRAYSCALE)
    i1 = cv2.imread(str(AIRSIM / "mav0/cam1/data" / cal["f1"][ts]), cv2.IMREAD_GRAYSCALE)
    o0, o1 = r.uygula(i0, i1)
    SONUC["airsim_kanca_rekt_piksel_farki"] = int(max(
        np.abs(o0.astype(int) - i0).max(), np.abs(o1.astype(int) - i1).max()))
    assert np.array_equal(o0, i0) and np.array_equal(o1, i1)


# ----------------------------------------------------------------------------
_KANCA_BETIGI = r'''
import os, sys, json
from pathlib import Path
import numpy as np
KOD = sys.argv[1]; kareler = [int(x) for x in sys.argv[2].split(",")]; cikti = sys.argv[3]
sys.path.insert(0, KOD); os.chdir(KOD)
import config
cfg = config.load()                       # SLAM_CONFIG=config_euroc.yaml
cfg["turbulence"].update(enabled=True, D_over_r0=5.0, seed=0)
from turb.hooks import injector_from_cfg
import turb.hooks as H
from stereo_front_end import StereoMatcher
m = StereoMatcher.__new__(StereoMatcher)  # model yüklemeden: yalnız G/Ç+rekt+kanca
m.base_path = Path(cfg["dataset"]["base_path"]); m.sequence = cfg["dataset"]["sequence"]
m.load_calibration(); m.load_timestamps()
m.turb = injector_from_cfg(cfg)
class _Dur(Exception): pass
yak = {}
def _yakala(r1, r2):
    yak["r1"], yak["r2"] = r1.copy(), r2.copy(); raise _Dur()
m._preprocess_stereo_pair = _yakala
ts_all = [str(t) for t in sorted({int(t) for t in m.cam0_df["timestamp"]})]
out = {}
for i in kareler:
    m._seq_frame_idx = i
    yak.clear()
    try:
        m.process_frame_pair(ts_all[i])
    except _Dur:
        pass
    out[f"c0_{i}"] = yak["r1"]; out[f"c1_{i}"] = yak["r2"]; out[f"ts_{i}"] = np.int64(ts_all[i])
np.savez(cikti, **out)
print(json.dumps({"dtype": str(yak["r1"].dtype), "hooks": H.__file__,
                  "params": {k: v for k, v in vars(m.turb.p).items()}}, default=str))
'''


def test_kanca_esdegerligi(kareler=(400, 401, 650, 899)):
    _gerek(EUROC / "MH01")
    _gerek(KOD / "stereo_front_end.py")
    with tempfile.TemporaryDirectory() as t:
        t = Path(t)
        bet = t / "kanca.py"
        bet.write_text(_KANCA_BETIGI)
        env = dict(os.environ, SLAM_CONFIG="config_euroc.yaml", CUDA_VISIBLE_DEVICES="")
        r = subprocess.run([PY, str(bet), str(KOD), ",".join(map(str, kareler)),
                            str(t / "kanca.npz")], env=env, capture_output=True, text=True)
        if r.returncode != 0:
            raise AssertionError("kanca betiği düştü:\n" + r.stderr[-3000:])
        bilgi = json.loads(r.stdout.strip().splitlines()[-1])
        k = np.load(t / "kanca.npz")
        spec = ",".join(f"{i}:{i + 1}" for i in kareler)
        s = _uret(EUROC / "MH01", t / "u", spec)
        farklar = {}
        for i in kareler:
            ts = int(k[f"ts_{i}"])
            for cam in (0, 1):
                h = k[f"c{cam}_{i}"]
                g = cv2.imread(str(s / f"mav0/cam{cam}/data/{ts}.png"), cv2.IMREAD_UNCHANGED)
                assert g.dtype == np.uint8 and g.shape == h.shape
                farklar[f"{i}/cam{cam}"] = int(np.abs(h.astype(int) - g.astype(int)).max())
        SONUC["kanca_dtype"] = bilgi["dtype"]
        SONUC["kanca_hooks_modulu"] = bilgi["hooks"]
        SONUC["kanca_maks_piksel_farki"] = farklar
        # Kanca parametreleri == üreticinin meta parametreleri
        meta = json.loads((s / "turb_meta.json").read_text())["parametreler"]
        fark_p = {kk: (v, meta.get(kk)) for kk, v in bilgi["params"].items()
                  if str(v) != str(meta.get(kk)) and
                  not (isinstance(meta.get(kk), list) and str(tuple(meta[kk])) == str(v))}
        SONUC["kanca_parametre_farki"] = fark_p
        assert not fark_p, fark_p
        assert bilgi["dtype"] == "uint8"
        assert max(farklar.values()) == 0, farklar


# ----------------------------------------------------------------------------
def main(filtre=None):
    testler = [(n, f) for n, f in sorted(globals().items())
               if n.startswith("test_") and callable(f) and (not filtre or filtre in n)]
    ok = True
    for n, f in testler:
        try:
            f()
            print(f"GEÇTİ  {n}")
        except Atla as e:
            print(f"ATLANDI {n}: {e}")
        except Exception as e:  # noqa: BLE001
            ok = False
            print(f"KALDI  {n}: {type(e).__name__}: {e}")
    print(json.dumps(SONUC, indent=1, ensure_ascii=False, default=str))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1] if len(sys.argv) > 1 else None) else 1)
