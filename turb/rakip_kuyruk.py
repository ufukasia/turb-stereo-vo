"""
Faz 2 rakip ızgarası kuyruğu (ON_KAYIT_FAZ2 §1: ORB-SLAM3 stereo, AirSLAM, DROID-SLAM; D/r0 {0, 8, 16, 24};
tohum 0-2; σ = 1 DN, gürültü tohumu = hücre tohumu; 9 sekans; kapsam1m pencereleri). ON_TEST_RAKIPLER.md §4.

Protokol (ön-test R2 kararı): ÖNEK koşumu. Her hücre için `turb/uret_sekans.py` kare 0 -> son pencere sonu
(max(başlangıç) + 400) bozulmuş sekansı üretir (PyEco'nun arka-ucunun koştuğu karelerin AYNISI), üç rakip sırayla
(aynı anda tek rakip süreci) koşulur, pencereler PyEco ile aynı fonksiyonla kesilir, kayıtlar JSONL'e eklenir,
üretilen sekans silinir (--sil). Devam ettirilebilir: tamamlanmış (rakip, sekans, D/r0, tohum, tekrar) atlanır.

    export CUDA_DEVICE_ORDER=PCI_BUS_ID      # (DROID koşucusu ayrıca CUDA_VISIBLE_DEVICES=1 alır)
    python -m turb.rakip_kuyruk --kuru                   # denetim + plan, hiçbir şey koşmaz
    python -m turb.rakip_kuyruk --sil                    # koşum (tek örnek kilidi; devam ettirilebilir)
    python -m turb.rakip_kuyruk --denetle                # kayıt denetimi (koşum sırasında/sonunda)

Güvenlik kapıları (her başlatmada): RTX 3090 tek; çıktı kökü `/` dosya sisteminde değil ve >= 20 GB boş;
ön kayıt sha256 = evidence/on_kayit_sha256.txt; pencere dosyası sha256 kayıtlı; manifest_*.json yazılır.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from turb import rakip_kos as RK  # noqa: E402

VERI = Path(os.environ.get("SLAM_DATA_ROOT", "/data/sdata"))
CIKTI_KOK = Path("/data/sdata_turb_rakip/izgara")
EV = ROOT / "evidence" / "FAZ2_rakip"
KAYIT = EV / "rakip_kayitlari.jsonl"
URETIM = EV / "uretim.jsonl"
KOSUM_DIZIN = EV / "kosum"
SEKANSLAR = ("MH01", "MH02", "MH03", "MH04", "MH05", "V1_01_easy", "V1_02_medium", "V2_01_easy", "V2_02_medium")
DR = (0, 8, 16, 24)
TOHUMLAR = (0, 1, 2)
SIGMA = 1.0
#: PyEco ile aynı rektifikasyon (pyeco_turb/config.yaml frontend.stereo_extrinsics: sensor_yaml; ön-test G0f)
STEREO_DIS = "sensor_yaml"
#: "onek" (ön-test R2 kararı: kare 0 -> son pencere sonu, PyEco protokolü) | "pencere" (her pencere ayrı kısa sekans)
PROTOKOL = "onek"
#: tohum başına tekrar sayısı (ön-test R5 kararı; ON_TEST_RAKIPLER.md §4)
TEKRAR = {"orbslam3": 3, "orbslam3_si": 3, "airslam": 1, "droid": 1}
ON_KAYIT = ROOT.parent / "ON_KAYIT_FAZ2.md"


def hucreler(sekanslar=SEKANSLAR, drs=DR, tohumlar=TOHUMLAR):
    """Sıra: tohum -> sekans -> D/r0. Böylece bir (sekans, tohum)'un bütün D/r0'ları art arda biter ve eşli
    oranlar erken denetlenebilir."""
    return [(s, d, t) for t in tohumlar for s in sekanslar for d in drs]


def onek_son(sekans: str) -> int:
    return max(RK.pencereler(sekans)) + RK.WINDOW_LEN


def kareler_belirtimi(sekans: str) -> str:
    if PROTOKOL == "pencere":       # yalnız pencere kareleri (türbülans mutlak indeksle -> aynı baytlar)
        return ",".join(f"{w}:{w + RK.WINDOW_LEN}" for w in RK.pencereler(sekans))
    return f"0:{onek_son(sekans)}"


def sekans_adi(sekans: str, dr: float, tohum: int) -> str:
    from turb.uret_sekans import arg_ayristirici, parametreleri_kur, cikti_adi, kareleri_ayristir, asl_oku, Rektifikator
    g = VERI / sekans
    a = arg_ayristirici().parse_args(["--girdi", str(g), "--D_over_r0", str(dr), "--tohum", str(tohum),
                                      "--gurultu_dn", str(SIGMA), "--kareler", kareler_belirtimi(sekans),
                                      "--stereo_dis", STEREO_DIS])
    cal = asl_oku(g)
    p, _ = parametreleri_kur(a, Rektifikator(cal, "otomatik", STEREO_DIS).baseline)
    _, ar = kareleri_ayristir(a.kareler, len(cal["ts"]))
    return cikti_adi(g.name, p, ar, SIGMA, None, STEREO_DIS)


def kayitlar(path=KAYIT) -> list[dict]:
    if not Path(path).exists():
        return []
    out = []
    for ln in open(path):
        try:
            out.append(json.loads(ln))
        except json.JSONDecodeError:
            pass
    return out


def tamamlanan(R) -> set:
    """Tamamlanmış (rakip, sekans, D/r0, tohum, tekrar[, pencere]) anahtarları (altyapı hatalı kayıtlar sayılmaz)."""
    R = [r for r in R if r.get("protokol") == PROTOKOL and not (r.get("kosum") or {}).get("altyapi_hatasi")]
    if PROTOKOL == "pencere":
        return {(r["rakip"], r["sequence"], float(r["D_over_r0"]), int(r["seed"]), int(r["tekrar"]), int(r["window"]))
                for r in R}
    n = Counter((r["rakip"], r["sequence"], float(r["D_over_r0"]), int(r["seed"]), int(r["tekrar"])) for r in R)
    return {k for k, v in n.items() if v >= len(RK.pencereler(k[1]))}


def isler_hucre(h, rakipler, bitti) -> list:
    """Bir hücrenin kalan koşumları: (rakip, tekrar, pencere|None)."""
    s, dr, t = h
    if PROTOKOL == "pencere":
        return [(rk, k, w) for rk in rakipler for k in range(TEKRAR[rk]) for w in RK.pencereler(s)
                if (rk, s, float(dr), int(t), k, w) not in bitti]
    return [(rk, k, None) for rk in rakipler for k in range(TEKRAR[rk]) if (rk, s, float(dr), int(t), k) not in bitti]


def uret(sekans: str, dr: float, tohum: int, isci: int) -> tuple[Path, dict]:
    ad = sekans_adi(sekans, dr, tohum)
    hedef = CIKTI_KOK / ad
    bilgi = {"ad": ad, "yeniden_kullanildi": False}
    if hedef.exists() and (hedef / "turb_meta.json").exists():
        bilgi["yeniden_kullanildi"] = True
    else:
        cmd = [RK.PY_TURB, "-m", "turb.uret_sekans", "--girdi", str(VERI / sekans), "--cikti_kok", str(CIKTI_KOK),
               "--D_over_r0", str(dr), "--tohum", str(tohum), "--gurultu_dn", str(SIGMA),
               "--kareler", kareler_belirtimi(sekans), "--isci", str(isci), "--uzerine_yaz", "--stereo_dis", STEREO_DIS]
        y0, t0 = os.getloadavg(), time.time()
        p = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True,
                           env={**os.environ, "OMP_NUM_THREADS": "2"})
        if p.returncode != 0:
            raise RuntimeError(f"üretim başarısız: {ad}\n{p.stderr[-3000:]}")
        bilgi.update(komut=cmd, sure_s=round(time.time() - t0, 1), yuk_once=[round(x, 1) for x in y0],
                     yuk_sonra=[round(x, 1) for x in os.getloadavg()])
    u = json.loads((hedef / "turb_uretim.json").read_text())
    bilgi.update(sekans=sekans, D_over_r0=float(dr), seed=int(tohum), sensor_noise_dn=SIGMA,
                 n_kare=u.get("n_kare"), boyut_bayt=u.get("boyut_bayt"), sha256=RK.sekans_sha256(hedef),
                 meta_sha256=RK.sha256_dosya(hedef / "turb_meta.json"), zaman=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(URETIM, "a") as f:
        f.write(json.dumps(bilgi, ensure_ascii=False) + "\n")
    return hedef, bilgi


def kapilar(sert=True) -> dict:
    """Başlatma kapıları; herhangi biri tutmazsa (sert) çık."""
    sorun = []
    g = RK.gpu_3090()
    if g["name"] != RK.BEKLENEN_GPU:
        sorun.append(f"GPU adı {g['name']}")
    CIKTI_KOK.mkdir(parents=True, exist_ok=True)
    if os.stat(CIKTI_KOK).st_dev == os.stat("/").st_dev:
        sorun.append(f"çıktı kökü {CIKTI_KOK} '/' dosya sisteminde")
    bos_gb = shutil.disk_usage(CIKTI_KOK).free / 1e9
    if bos_gb < 20:
        sorun.append(f"çıktı diskinde {bos_gb:.0f} GB boş < 20 GB")
    kayitli = (ROOT / "evidence/on_kayit_sha256.txt").read_text().split()
    if RK.sha256_dosya(ON_KAYIT) not in kayitli:
        sorun.append("ON_KAYIT_FAZ2.md sha256 kayıtlı değerle aynı değil")
    for s in SEKANSLAR:
        if not (VERI / s / "mav0/cam0/data.csv").exists():
            sorun.append(f"veri yok: {VERI / s}")
    for p in (RK.PY_DROID, RK.AIRSLAM_WS / "devel/lib/air_slam/visual_odometry",
              RK.ORB_ROOT / "Examples/Stereo/stereo_asl_headless"):
        if not Path(p).exists():
            sorun.append(f"eksik: {p}")
    d = {"gpu": g, "cikti_bos_gb": round(bos_gb, 1), "sorun": sorun}
    if sorun and sert:
        raise SystemExit("KAPI KALDI: " + "; ".join(sorun))
    return d


def manifest(plan, kapi) -> Path:
    EV.mkdir(parents=True, exist_ok=True)
    k = RK.kimlik()
    m = {"zaman": time.strftime("%Y-%m-%d %H:%M:%S"), "kimlik_fp": RK.kimlik_ozeti(k), "kimlik": k,
         "on_kayit_sha256": RK.sha256_dosya(ON_KAYIT),
         "pencere_dosyasi": str(ROOT / "evidence/pencereler_faz2_kapsam1m.json"),
         "pencere_sha256": RK.sha256_dosya(ROOT / "evidence/pencereler_faz2_kapsam1m.json"),
         "uret_sekans_sha256": RK.sha256_dosya(ROOT / "turb/uret_sekans.py"),
         "stereo_dis": STEREO_DIS,
         "hooks_sha256": RK.sha256_dosya(ROOT / "turb/hooks.py"),
         "uretici_yapilandirma_sha256": RK.sha256_dosya(ROOT.parent / "kod/config_euroc.yaml"),
         "kuyruk_sha256": RK.sha256_dosya(Path(__file__)),
         "plan": {"sekanslar": SEKANSLAR, "D_over_r0": DR, "tohumlar": TOHUMLAR, "sigma_dn": SIGMA,
                  "tekrar": TEKRAR, "protokol": PROTOKOL + (" (kare 0 -> max(pencere)+400)" if PROTOKOL == "onek"
                                                            else " (her pencere ayrı 400 kare)"), "n_hucre": len(plan),
                  "n_rakip_kosumu": len(plan) * sum(TEKRAR.values())},
         "kapilar": kapi, "ortam": {k_: os.environ.get(k_) for k_ in ("CUDA_DEVICE_ORDER", "CUDA_VISIBLE_DEVICES",
                                                                     "SLAM_DATA_ROOT", "OMP_NUM_THREADS")}}
    p = EV / f"manifest_{time.strftime('%Y%m%d_%H%M%S')}.json"
    p.write_text(json.dumps(m, indent=1, ensure_ascii=False))
    return p


def denetle(path=KAYIT, sert=False) -> bool:
    import numpy as np
    R = [r for r in kayitlar(path) if r.get("protokol") == PROTOKOL]
    print(f"[denetim] {path}: {len(R)} pencere kaydı")
    if not R:
        return True
    ok = True
    fps = Counter(r.get("kimlik_fp") for r in R)
    gpus = Counter((r["rakip"], (r.get("kosum") or {}).get("gpu_adi")) for r in R if r["rakip"] not in ("orbslam3", "orbslam3_si"))
    alt = [r for r in R if (r.get("kosum") or {}).get("altyapi_hatasi")]
    print(f"  kimlik_fp: {dict(fps)}\n  GPU: {dict(gpus)}\n  altyapı hatası (sonuç değil): {len(alt)}")
    if len(fps) > 1:
        print("  !! birden çok kimlik_fp (kod/derleme koşum sırasında değişti)"); ok = False
    for (rk, g), n in gpus.items():
        if g and "3090" not in str(g):
            print(f"  !! {rk} GPU {g}"); ok = False
    R = [r for r in R if not (r.get("kosum") or {}).get("altyapi_hatasi")]
    by = defaultdict(list)
    for r in R:
        by[(r["rakip"], float(r["D_over_r0"]))].append(r)
    for (rk, dr), rr in sorted(by.items()):
        f = sum(r["durum_p"] == "basarisiz" for r in rr)
        a = [r["m_ate_se3"] for r in rr if r["durum_p"] == "ok"]
        s = [r["kosum"]["sure_s"] for r in rr]
        print(f"  {rk:9s} D/r0={dr:4.0f}: {len(rr):4d} pencere, başarısız {f:3d} ({f / len(rr):.0%}), "
              f"ATE medyan {np.median(a) if a else float('nan'):.4f} m, koşum süresi medyan {np.median(s):.0f} s")
        if dr == 0 and f:
            print(f"  !! {rk} D/r0=0'da {f} başarısız pencere: incele"); ok = False
        if dr == 0 and a and max(a) > 0.15:
            print(f"  !! {rk} temiz ATE > 0.15 m: {max(a):.3f}"); ok = False
    print(f"[denetim] {'TAMAM' if ok else 'SORUN VAR — kuyruğu durdurup incele'}")
    if sert and not ok:
        sys.exit(2)
    return ok


def main():
    ap = argparse.ArgumentParser(description="Faz 2 rakip ızgarası kuyruğu (tek rakip süreci, RTX 3090)")
    ap.add_argument("--kuru", action="store_true")
    ap.add_argument("--denetle", action="store_true")
    ap.add_argument("--sil", action="store_true", help="hücre bitince üretilen sekansı sil")
    ap.add_argument("--isci", type=int, default=4, help="üretim süreç sayısı (ön-testte 4)")
    ap.add_argument("--sekanslar", default=",".join(SEKANSLAR))
    ap.add_argument("--tohumlar", default=",".join(map(str, TOHUMLAR)))
    ap.add_argument("--dr", default=",".join(map(str, DR)))
    ap.add_argument("--rakipler", default=",".join(RK.RAKIP_LISTE))
    ap.add_argument("--en_cok_hucre", type=int, default=0, help=">0: bu kadar hücreden sonra dur (erken denetim)")
    ap.add_argument("--protokol", default="onek", choices=["onek", "pencere"],
                    help="onek (vars., ön-test R2 kararı) | pencere (her pencere ayrı kısa sekans)")
    ap.add_argument("--ev", default=None, help="kayıt dizini (vars.: evidence/FAZ2_rakip); uçtan uca test için ayrı dizin")
    a = ap.parse_args()
    global PROTOKOL
    PROTOKOL = a.protokol
    if a.ev:                            # uçtan uca test: kayıtlar, koşum dosyaları ve sekanslar ayrı yerde
        global EV, KAYIT, URETIM, KOSUM_DIZIN, CIKTI_KOK
        EV = Path(a.ev).resolve()
        KAYIT, URETIM, KOSUM_DIZIN = EV / "rakip_kayitlari.jsonl", EV / "uretim.jsonl", EV / "kosum"
        CIKTI_KOK = CIKTI_KOK.parent / ("test_" + EV.name)
        EV.mkdir(parents=True, exist_ok=True)
    if a.denetle:
        denetle(KAYIT)
        return
    plan = hucreler(a.sekanslar.split(","), [float(x) for x in a.dr.split(",")],
                    [int(x) for x in a.tohumlar.split(",")])
    rakipler = a.rakipler.split(",")
    kapi = kapilar(sert=not a.kuru)
    bitti = tamamlanan(kayitlar(KAYIT))
    isler = [(h, x) for h in plan for x in isler_hucre(h, rakipler, set())]
    kalan = [(h, x) for h in plan for x in isler_hucre(h, rakipler, bitti)]
    print(f"[plan] {len(plan)} hücre, {len(isler)} rakip koşumu; tamamlanmış {len(isler) - len(kalan)}, kalan "
          f"{len(kalan)}; kapılar: {kapi}")
    if a.kuru:
        for h in plan[:6]:
            print("  ", h, sekans_adi(*h), "kareler:", kareler_belirtimi(h[0]))
        print("  ...")
        return
    lk = open(EV / "rakip_kuyruk.kilit", "w")
    try:
        fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("başka bir rakip kuyruğu koşuyor (kilit)")
    mp = manifest(plan, kapi)
    kfp = json.loads(mp.read_text())["kimlik_fp"]
    print(f"[manifest] {mp}  kimlik_fp={kfp}", flush=True)
    n_hucre = 0
    bitti = tamamlanan(kayitlar(KAYIT))
    gerekli = [(h, isler_hucre(h, rakipler, bitti)) for h in plan]
    gerekli = [(h, g) for h, g in gerekli if g]
    # Üretim bir hücre önden, rakip koşumlarıyla eşzamanlı (diskte en çok 2 sekans; rakip süreci yine tek)
    from concurrent.futures import ThreadPoolExecutor
    ex = ThreadPoolExecutor(1)
    sonraki = ex.submit(uret, *gerekli[0][0], a.isci) if gerekli else None
    for j, ((s, dr, t), gerek) in enumerate(gerekli):
        if shutil.disk_usage(CIKTI_KOK).free / 1e9 < 20:
            raise SystemExit("çıktı diskinde < 20 GB boş: durdu")
        seq, ub = sonraki.result()
        sonraki = ex.submit(uret, *gerekli[j + 1][0], a.isci) if j + 1 < len(gerekli) else None
        print(f"[üretim] {ub['ad']} {ub.get('sure_s', 'yeniden')} s, {ub['boyut_bayt'] / 1e9:.2f} GB, "
              f"sha256 {ub['sha256'][:12]}", flush=True)
        for rk, k, w in gerek:
            if w is None:
                hedef, wins, et = seq, None, f"{s}_dr{dr:g}_s{t}"
            else:                       # protokol (a): pencere görünümü (data.csv süzülür, görüntüler bağlanır)
                hedef = RK.pencere_gorunumu(seq, seq.parent / f"{seq.name}__gorunum_w{w}", w)
                wins, et = [w], f"{s}_dr{dr:g}_s{t}_w{w}"
            for deneme in range(3):
                recs = RK.tek_kosum(rk, hedef, VERI / s, s, dr, t, SIGMA, PROTOKOL, KOSUM_DIZIN, KAYIT, tekrar=k,
                                    wins=wins, etiket=et, kimlik_fp=kfp,
                                    ek={"uretim_sha256": ub["sha256"], "uretim_ad": ub["ad"]})
                if not recs[0]["kosum"].get("altyapi_hatasi"):
                    break
                print(f"[altyapı] {rk} {s} dr{dr:g} s{t} r{k}: {recs[0]['kosum']['altyapi_hatasi']}; "
                      f"{deneme + 1}. yeniden deneme 180 s sonra", flush=True)
                time.sleep(180)
        if a.sil:
            shutil.rmtree(seq, ignore_errors=True)
            for g in seq.parent.glob(f"{seq.name}__gorunum_w*"):
                shutil.rmtree(g, ignore_errors=True)
        n_hucre += 1
        if n_hucre in (1, 4, 12) or n_hucre % 36 == 0:
            denetle(KAYIT)
        if a.en_cok_hucre and n_hucre >= a.en_cok_hucre:
            print(f"[dur] {n_hucre} hücre (--en_cok_hucre)")
            if sonraki is not None:      # önden üretilen sekans: bitmesini bekle, sil (yarım kalmasın)
                try:
                    nseq, _ = sonraki.result()
                    if a.sil:
                        shutil.rmtree(nseq, ignore_errors=True)
                except Exception as e:  # noqa: BLE001
                    print("[uyarı] önden üretim:", e)
            break
    denetle(KAYIT)
    print("BITTI", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    main()
