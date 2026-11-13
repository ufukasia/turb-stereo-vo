#!/usr/bin/env python
"""Derin onarma kolu (ON_KAYIT_FAZ2.md Sapma 3): öğrenilmiş türbülans giderme ağı VO'nun önünde ön-işlemci.

Hat: turb/uret_sekans.py (bozulmuş rektifiye ASL sekansı) -> kod/derin_onarim/betikler/onar_sekans.py (MambaTM)
-> PyEco, kanca kapalı. Parmak izli kodda değişiklik yok; bu dosya parmak izi kümesinin dışındadır.

    python onarim_kos.py plan
    python onarim_kos.py hazirla [--isci 4] [--yalniz V1_01_easy] [--kipler N,O]   # üret + onar (devam ettirilebilir)
    python onarim_kos.py kos --gpu_jobs 2 --jobs 8 [--yalniz ...]     # PyEco (devam ettirilebilir)
    python onarim_kos.py olc                                          # görüntü düzeyi ölçüler (temize göre)
    python onarim_kos.py sil                                          # kaydı tamamlanan hücrelerin sekanslarını sil

Ön-test hücreleri (--ontest): V1_01_easy R0 N (kanca yoluyla kıyas) ve V1_01_easy temiz + O (D/r0 = 0).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import threading
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from turb.runner import Cell, windows_for, _frontend_worker, _backend_worker, _js, code_fingerprint  # noqa: E402
from turb.experiments import EUROC_FAZ2, FAZ2_WINDOWS_FILE, FAZ2_COV_EVERY, FAZ2_K_A0C  # noqa: E402

BEKLENEN_FP = "08e28b83407c"
SEQS = ["V1_01_easy", "MH01", "MH02", "V1_02_medium", "MH03", "V2_01_easy", "MH04", "V2_02_medium", "MH05"]
KOK = "/data/sdata_turb_onarim"
EV = os.path.join(ROOT, "evidence", "ONARIM")
DO_KOK = os.path.join(os.path.dirname(ROOT), "kod", "derin_onarim")
PY_TURB = os.path.join(os.path.dirname(ROOT), "kod", "turbvio", "bin", "python")
V_HIZLI, D_AP, DT = 5.0, 0.10, 0.05
SIGMA_DN = 1.0
KOLLAR = ("A0", "A1", "A3")
# (rejim, D/r0, ön-işlem kipleri): N = onarım yok, O = çevrimdışı, C = nedensel L=5
KOSULLAR = [("R5", 8, ("N", "O", "C")), ("R0", 8, ("N", "O"))]
ONTEST = [("R0", 8, ("N",)), ("T", 0, ("N", "O"))]          # T = temiz (D/r0 = 0, yalnız sensör gürültüsü)
ONAR_ARG = {"O": ["--mode", "offline", "--clip", "20", "--ovl", "2"], "C": ["--mode", "causal", "--win", "5"]}


def tau0_over_dt(rejim: str, dr: float) -> float:
    return 1.0 if rejim in ("R0", "T") else 0.314 * (D_AP / dr) / (V_HIZLI * DT)


def onek_son(seq: str) -> int:
    return int(max(windows_for(seq, FAZ2_WINDOWS_FILE)) + 400)


def ad_bozuk(seq, rejim, dr) -> str:
    return f"{seq}__{rejim}_dr{dr:g}_s0"


def ad_mod(seq, rejim, dr, mod) -> str:
    return ad_bozuk(seq, rejim, dr) + ("" if mod == "N" else f"__{mod}")


def plan(kosullar, seqs):
    """[(seq, rejim, dr, mod)], sıra: sekans -> koşul -> kip (erken eşli denetim için)."""
    return [(s, r, d, m) for s in seqs for r, d, mods in kosullar for m in mods]


def hazir(seq, rejim, dr, mod) -> bool:
    d = os.path.join(KOK, ad_mod(seq, rejim, dr, mod))
    return os.path.exists(os.path.join(d, "turb_meta.json" if mod == "N" else "onarim_meta.json"))


def hucre(seq, rejim, dr, mod) -> Cell:
    kollar = ("A0",) if rejim == "T" else KOLLAR
    return Cell(dataset="euroc", base_path=KOK, sequence=ad_mod(seq, rejim, dr, mod), seed=0,
                D_over_r0=float(dr), tau0_over_dt=tau0_over_dt(rejim, dr),
                windows=windows_for(seq, FAZ2_WINDOWS_FILE), window_len=400, nees=True, cov_every=FAZ2_COV_EVERY,
                psf_legacy_scale=False, sensor_noise_dn=0.0, noise_seed=None, cache_dump=False,
                cfg_overrides=(("turbulence", "enabled", False),), arms=kollar,
                backend_argv=("--turb_cal_factor", str(float(FAZ2_K_A0C))), tag=f"onarim_{rejim}_{mod}")


# ----------------------------------------------------------------------------- hazırlık (üret + onar)
def _yaz_satir(path, d):
    with open(path, "a") as f:
        f.write(json.dumps(d, ensure_ascii=False, default=_js) + "\n")


def uret(seq, rejim, dr, isci) -> dict:
    ad = ad_bozuk(seq, rejim, dr)
    hedef = os.path.join(KOK, ad)
    if os.path.exists(os.path.join(hedef, "turb_meta.json")):
        return {"ad": ad, "yeniden": True}
    cmd = [PY_TURB, "-m", "turb.uret_sekans", "--girdi", os.path.join(EUROC_FAZ2, seq), "--cikti_kok", KOK,
           "--ad", ad, "--D_over_r0", str(float(dr)), "--tohum", "0", "--gurultu_dn", str(SIGMA_DN),
           "--tau0_over_dt", repr(tau0_over_dt(rejim, dr)), "--kareler", f"0:{onek_son(seq)}",
           "--isci", str(isci), "--uzerine_yaz", "--stereo_dis", "sensor_yaml"]
    t0 = time.time()
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, env={**os.environ, "OMP_NUM_THREADS": "2"})
    if p.returncode != 0:
        raise RuntimeError(f"üretim başarısız: {ad}\n{p.stderr[-2000:]}")
    return {"ad": ad, "yeniden": False, "sure_s": round(time.time() - t0, 1), "komut": cmd[1:],
            "son_satir": p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ""}


def onar(seq, rejim, dr, mod) -> dict:
    ad = ad_mod(seq, rejim, dr, mod)
    hedef = os.path.join(KOK, ad)
    if os.path.exists(os.path.join(hedef, "onarim_meta.json")):
        return {"ad": ad, "yeniden": True}
    env = {**os.environ, "CUDA_DEVICE_ORDER": "PCI_BUS_ID", "CUDA_VISIBLE_DEVICES": "1",
           "TRITON_PTXAS_PATH": "/usr/bin/ptxas", "TRITON_CACHE_DIR": os.path.join(DO_KOK, "triton_cache")}
    cmd = [os.path.join(DO_KOK, "venv", "bin", "python"), os.path.join(DO_KOK, "betikler", "onar_sekans.py"),
           "--girdi", os.path.join(KOK, ad_bozuk(seq, rejim, dr)), "--cikti", hedef, "--model", "mambatm"] + ONAR_ARG[mod]
    t0 = time.time()
    p = subprocess.run(cmd, cwd=DO_KOK, capture_output=True, text=True, env=env)
    if p.returncode != 0 or not os.path.exists(os.path.join(hedef, "onarim_meta.json")):
        shutil.rmtree(hedef, ignore_errors=True)
        raise RuntimeError(f"onarım başarısız: {ad}\n{p.stderr[-2000:]}")
    with open(os.path.join(hedef, "onarim_meta.json")) as f:
        m = json.load(f)
    return {"ad": ad, "yeniden": False, "sure_s": round(time.time() - t0, 1), "komut": cmd[1:],
            "kameralar": m["kameralar"], "agirlik_sha256": m["agirlik_sha256"], "depo_commit": m["depo_commit"],
            "gpu": m["gpu"]}


def hazirla(P, isci):
    os.makedirs(KOK, exist_ok=True)
    log = os.path.join(EV, "hazirlik.jsonl")
    bozuklar = list(dict.fromkeys((s, r, d) for s, r, d, _ in P))
    t_bas = time.time()
    for i, (s, r, d) in enumerate(bozuklar):
        b = uret(s, r, d, isci)
        if not b["yeniden"]:
            _yaz_satir(log, {"is": "uret", "sekans": s, "rejim": r, "dr": d, "zaman": time.strftime("%F %T"), **b})
        print(f"[hazirla] {i + 1}/{len(bozuklar)} üretildi: {b['ad']} ({b.get('sure_s', 0)} s)", flush=True)
        for m in [m for s2, r2, d2, m in P if (s2, r2, d2) == (s, r, d) and m != "N"]:
            o = onar(s, r, d, m)
            if not o["yeniden"]:
                _yaz_satir(log, {"is": "onar", "sekans": s, "rejim": r, "dr": d, "mod": m,
                                 "zaman": time.strftime("%F %T"), **o})
            sk = o.get("kameralar", {}).get("cam0", {}).get("s_kare")
            print(f"[hazirla]    onarıldı: {o['ad']} ({o.get('sure_s', 0)} s"
                  + (f", {sk:.3f} s/kare" if sk else "") + f"); toplam {(time.time() - t_bas) / 60:.1f} dk", flush=True)
    print(f"[hazirla] bitti, {(time.time() - t_bas) / 60:.1f} dk")


# ----------------------------------------------------------------------------- PyEco koşumu
def okunan(out):
    done = set()
    if os.path.exists(out):
        with open(out) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get("status") != "cell_fail" and "arm" in r:
                    done.add((r.get("key"), r.get("window"), r.get("arm")))
    return done


def kos(P, out, jobs, gpu_jobs, bekle):
    fp = code_fingerprint()
    print(f"[onarim] code_fp {fp} (beklenen {BEKLENEN_FP})")
    if fp != BEKLENEN_FP:
        sys.exit("[onarim] DUR: kod parmak izi ana ızgaradan farklı")
    spool = "/tmp/turb_spool_onarim"
    os.makedirs(spool, exist_ok=True)
    bilgi = {hucre(*p).key(): p for p in P}
    assert len(bilgi) == len(P), "anahtar çakışması"
    kilit = threading.Lock()
    t_bas = time.time()
    sayac = {"n": 0}
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=gpu_jobs, maxtasksperchild=8) as gp, \
            ctx.Pool(processes=jobs, maxtasksperchild=4) as cp, open(out, "a") as fout:
        def yaz(recs):
            with kilit:
                for r in recs:
                    s, rj, d, m = bilgi.get(r.get("key"), (None,) * 4)
                    r["onarim"] = {"sekans": s, "rejim": rj, "dr": d, "mod": m}
                    fout.write(json.dumps(r, default=_js) + "\n")
                fout.flush()
                sayac["n"] += 1
                print(f"[onarim] {sayac['n']} hücre yazıldı, {(time.time() - t_bas) / 60:.1f} dk", flush=True)

        def hata(e):
            print(f"[onarim] İŞ HATASI: {type(e).__name__}: {e}", flush=True)

        bekleyen, gonderilen = [], set()
        while True:
            done = okunan(out)
            tamam = lambda c: all((c.key(), int(w), x) in done for w in c.windows for x in c.arms)   # noqa: E731
            kalan = [p for p in P if hucre(*p).key() not in gonderilen and not tamam(hucre(*p))]
            simdi = [p for p in kalan if hazir(*p)]
            if simdi:
                print(f"[onarim] {len(simdi)} hücre gönderiliyor (hazır olmayan {len(kalan) - len(simdi)})", flush=True)
                for p in simdi:
                    gonderilen.add(hucre(*p).key())
                for meta in gp.imap(_frontend_worker, [(hucre(*p), spool) for p in simdi]):
                    c = meta["cell"]
                    print(f"[onarim] ön-uç {'tamam' if meta.get('ok') else 'BAŞARISIZ'}: {c.key()} "
                          f"({meta.get('t_frontend', 0):.0f} s) {meta.get('error', '')}", flush=True)
                    bekleyen.append(cp.apply_async(_backend_worker, (meta,), callback=yaz, error_callback=hata))
                continue
            if not kalan or not bekle:
                break
            time.sleep(60)          # hazırlık süreci sekans üretmeyi sürdürüyor
        for b in bekleyen:
            b.wait()
    print(f"[onarim] bitti, {(time.time() - t_bas) / 60:.1f} dk")


# ----------------------------------------------------------------------------- görüntü düzeyi ölçü
def olc(P, out, adim=20):
    """Temiz rektifiye kareye (D/r0 = 0, gürültüsüz) göre PSNR / SSIM / Farneback RMS; cam0, her `adim` karede bir."""
    import cv2
    import numpy as np
    sys.path.insert(0, os.path.join(DO_KOK, "betikler"))
    from analiz import psnr, ssim, rms_flow
    for s in dict.fromkeys(p[0] for p in P):
        temiz = os.path.join(KOK, f"{s}__temiz_gurultusuz")
        if not os.path.exists(os.path.join(temiz, "turb_meta.json")):
            cmd = [PY_TURB, "-m", "turb.uret_sekans", "--girdi", os.path.join(EUROC_FAZ2, s), "--cikti_kok", KOK,
                   "--ad", os.path.basename(temiz), "--D_over_r0", "0", "--tohum", "0", "--gurultu_dn", "0",
                   "--kareler", f"0:{onek_son(s)}", "--isci", "4", "--uzerine_yaz", "--stereo_dis", "sensor_yaml"]
            subprocess.run(cmd, cwd=ROOT, check=True, capture_output=True)
        adlar = sorted(os.listdir(os.path.join(temiz, "mav0", "cam0", "data")))[::adim]
        rd = lambda d, n: cv2.imread(os.path.join(d, "mav0", "cam0", "data", n), cv2.IMREAD_GRAYSCALE)   # noqa: E731
        for (s2, r, d, m) in [p for p in P if p[0] == s]:
            dz = os.path.join(KOK, ad_mod(s, r, d, m))
            if not hazir(s, r, d, m):
                continue
            v = {"psnr": [], "ssim": [], "rms_akis": [], "d_ort": []}
            for n in adlar:
                a, c = rd(dz, n), rd(temiz, n)
                v["psnr"].append(psnr(a, c)); v["ssim"].append(ssim(a, c)); v["rms_akis"].append(rms_flow(a, c))
                v["d_ort"].append(float(a.mean()) - float(c.mean()))
            rec = {"sekans": s, "rejim": r, "dr": d, "mod": m, "n_kare": len(adlar),
                   **{k: float(np.mean(x)) for k, x in v.items()}}
            _yaz_satir(out, rec)
            print(f"[olc] {s} {r} {m}: PSNR {rec['psnr']:.2f}  SSIM {rec['ssim']:.3f}  RMS akış {rec['rms_akis']:.3f} px  "
                  f"dOrt {rec['d_ort']:+.2f}", flush=True)
        shutil.rmtree(temiz, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("komut", choices=["plan", "hazirla", "kos", "olc", "sil"])
    ap.add_argument("--ontest", action="store_true")
    ap.add_argument("--yalniz", default="", help="virgülle ayrılmış sekanslar")
    ap.add_argument("--isci", type=int, default=4)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--gpu_jobs", type=int, default=2)
    ap.add_argument("--bekle", action="store_true", help="kos: hazır olmayan sekansları bekle")
    ap.add_argument("--kipler", default="", help="hazirla: yalnız bu ön-işlem kipleri (ör. N,O); boş = hepsi")
    a = ap.parse_args()
    os.makedirs(EV, exist_ok=True)
    seqs = [s for s in a.yalniz.split(",") if s] or (["V1_01_easy"] if a.ontest else SEQS)
    P = plan(ONTEST if a.ontest else KOSULLAR, seqs)
    out = os.path.join(EV, "on_test.jsonl" if a.ontest else "records.jsonl")
    if a.komut == "plan":
        done = okunan(out)
        for p in P:
            c = hucre(*p)
            n = sum((c.key(), int(w), x) in done for w in c.windows for x in c.arms)
            print(f"  {p[0]:14s} {p[1]} dr{p[2]:<2} {p[3]}  tau0/dt={c.tau0_over_dt:.4f}  hazır={int(hazir(*p))}  "
                  f"kayıt={n}/{len(c.windows) * len(c.arms)}  {c.sequence}")
        print(f"[plan] {len(P)} hücre, kollar {KOLLAR}; çıktı {out}")
    elif a.komut == "hazirla":
        kip = [k for k in a.kipler.split(",") if k]
        hazirla([p for p in P if not kip or p[3] in kip], a.isci)
    elif a.komut == "kos":
        kos(P, out, a.jobs, a.gpu_jobs, a.bekle)
    elif a.komut == "olc":
        olc(P, os.path.join(EV, "goruntu_olcu_on_test.jsonl" if a.ontest else "goruntu_olcu.jsonl"))
    elif a.komut == "sil":
        done = okunan(out)
        for p in P:
            c = hucre(*p)
            if all((c.key(), int(w), x) in done for w in c.windows for x in c.arms):
                # N dizini O/C onarımının girdisidir: yalnız aynı koşulun bütün kipleri bitince silinir
                es = [q for q in P if q[:3] == p[:3]]
                if p[3] != "N" or all(all((hucre(*q).key(), int(w), x) in done for w in hucre(*q).windows
                                          for x in hucre(*q).arms) for q in es):
                    shutil.rmtree(os.path.join(KOK, c.sequence), ignore_errors=True)
                    print("[sil]", c.sequence)


if __name__ == "__main__":
    main()
