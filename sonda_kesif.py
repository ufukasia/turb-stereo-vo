#!/usr/bin/env python
"""Keşif sondası (2026-11-11): kapı ölçekleme x zamansal rejim. Ön beyan: evidence/SONDA_KESIF/ON_BEYAN.md.

Kod değişikliği yok: yalnız Cell.cfg_overrides / backend_argv / tau0_over_dt. Bu dosya kod parmak izi
kümesinin (turb.runner.code_fingerprint) dışındadır.

    python sonda_kesif.py --kuru                      # plan + birim denetimleri, koşmaz
    python sonda_kesif.py --duman                     # tek sekans, tek pencere, D/r0=8, 4 koşul
    python sonda_kesif.py --gpu_jobs 3 --jobs 12      # sonda (devam ettirilebilir)
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import threading
import time
from dataclasses import asdict

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from turb.runner import (Cell, windows_for, _frontend_worker, _backend_worker, _fe_key,  # noqa: E402
                         _cache_dir, _js, code_fingerprint)
from turb.experiments import (EUROC_FAZ2, FAZ2_WINDOWS_FILE, FAZ2_NOISE_DN, FAZ2_COV_EVERY,  # noqa: E402
                              FAZ2_K_A0C)
from turb.backend import ArmConfig  # noqa: E402

BEKLENEN_FP = "08e28b83407c"
SEQS = ["MH01", "V1_01_easy", "MH02", "V1_02_medium", "MH03", "V2_01_easy", "MH04", "V2_02_medium", "MH05"]
V_HIZLI = 5.0                      # [m/s] R5 enine hız
D_AP, DT = 0.10, 0.05
KOSULLAR = [("R0", 8), ("R5", 8), ("R0", 16), ("R5", 4)]      # (rejim, D/r0), öncelik sırası


def sigma_tilt(dr: float) -> float:
    return float(ArmConfig(D_over_r0=float(dr)).sigma_px_oracle)


def tau0_over_dt(rejim: str, dr: float) -> float:
    return 1.0 if rejim == "R0" else 0.314 * (D_AP / dr) / (V_HIZLI * DT)


def fe_overrides(rejim: str, dr: float, fe: str) -> tuple:
    if fe == "fe0":
        return ()
    g = math.sqrt(1.0 + sigma_tilt(dr) ** 2)
    ov = [("frontend", "y_axis_threshold", round(2.0 * g, 4))]
    if rejim != "R0":
        ov.append(("frontend", "temporal_ransac_threshold", round(1.96 * g, 4)))
    return tuple(ov)


def be_genis_argv(dr: float) -> tuple:
    s = math.sqrt(1.0 + (math.hypot(1.0, 1.75) * sigma_tilt(dr)) ** 2)
    return ("--px_noise", f"{s:.6f}", "--stereo_sigma", f"{s:.6f}", "--pnp_reproj", f"{2.4477 * s:.6f}")


def gruplar(seqs, kosullar, tek_pencere=False, ek=""):
    """Her grup: (taban hücresi [A0, A1], geniş hücre [A0 + geniş arka-uç kapıları]); aynı ön-uç dökümü."""
    out = []
    for s in seqs:
        w = windows_for(s, FAZ2_WINDOWS_FILE)
        if tek_pencere:
            w = w[:1]
        for rejim, dr in kosullar:
            for fe in ("fe0", "fe1"):
                ortak = dict(dataset="euroc", base_path=EUROC_FAZ2, sequence=s, seed=0, D_over_r0=float(dr),
                             tau0_over_dt=tau0_over_dt(rejim, dr), windows=w, window_len=400, nees=True,
                             cov_every=FAZ2_COV_EVERY, psf_legacy_scale=False, sensor_noise_dn=FAZ2_NOISE_DN,
                             noise_seed=None, cache_dump=True, cfg_overrides=fe_overrides(rejim, dr, fe))
                cal = ("--turb_cal_factor", str(float(FAZ2_K_A0C)))
                taban = Cell(arms=("A0", "A1"), backend_argv=cal, tag=f"sonda{ek}_{rejim}_{fe}", **ortak)
                genis = Cell(arms=("A0",), backend_argv=cal + be_genis_argv(dr),
                             tag=f"sonda{ek}_{rejim}_{fe}_bg", **ortak)
                out.append((taban, genis))
    return out


def _genis_is(args):
    """CPU havuzunda: önbellekteki dökümü al (GPU'ya dokunmaz), geniş kapılı A0'ı koş, önbelleği sil."""
    cell, spool = args
    cpath = os.path.join(_cache_dir(), _fe_key(cell) + ".pkl")
    try:
        if not os.path.exists(cpath):
            return [{"cell": asdict(cell), "key": cell.key(), "code_fp": code_fingerprint(),
                     "status": "cell_fail", "error": "ön-uç önbelleği yok"}]
        meta = _frontend_worker((cell, spool))
        if not meta.get("cached"):
            return [{"cell": asdict(cell), "key": cell.key(), "code_fp": code_fingerprint(),
                     "status": "cell_fail", "error": "önbellek okunamadı"}]
        return _backend_worker(meta)
    finally:
        for p in (cpath, cpath[:-4] + "_quality.json"):
            try:
                os.remove(p)
            except OSError:
                pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kuru", action="store_true")
    ap.add_argument("--duman", action="store_true")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--gpu_jobs", type=int, default=3)
    a = ap.parse_args()

    ev = os.path.join(ROOT, "evidence", "SONDA_KESIF")
    os.makedirs(ev, exist_ok=True)
    if a.duman:
        G = gruplar(["V1_01_easy"], [("R0", 8), ("R5", 8)], tek_pencere=True, ek="duman")
        out = os.path.join(ev, "duman.jsonl")
        spool = "/tmp/turb_spool_sonda_duman"
    else:
        G = gruplar(SEQS, KOSULLAR)
        out = os.path.join(ev, "records.jsonl")
        spool = "/tmp/turb_spool_sonda"
    os.makedirs(spool, exist_ok=True)

    fp = code_fingerprint()
    print(f"[sonda] code_fp {fp} (beklenen {BEKLENEN_FP}); önbellek {_cache_dir()}; spool {spool}")
    if fp != BEKLENEN_FP:
        sys.exit("[sonda] DUR: kod parmak izi ana ızgaradan farklı")
    for rejim, dr in sorted(set((r, d) for r, d in (KOSULLAR if not a.duman else [("R0", 8), ("R5", 8)]))):
        t0 = tau0_over_dt(rejim, dr)
        v = 0.314 * (D_AP / dr) / (t0 * DT)
        print(f"  {rejim} D/r0={dr:2d}: tau0/dt={t0:.4f}  v={v:.3f} m/s  sigma_tilt={sigma_tilt(dr):.3f} px  "
              f"fe1={fe_overrides(rejim, dr, 'fe1')}  be_genis={be_genis_argv(dr)}")
    anahtarlar = [c.key() for g in G for c in g]
    assert len(set(anahtarlar)) == len(anahtarlar), "anahtar çakışması"
    for t, g in G:
        assert _fe_key(t) == _fe_key(g), "taban/geniş ön-uç anahtarı farklı"

    done = set()
    if os.path.exists(out):
        with open(out) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    if r.get("status") != "cell_fail":
                        done.add((r.get("key"), r.get("window"), r.get("arm")))
                except json.JSONDecodeError:
                    pass
    tamam = lambda c: all((c.key(), int(w), x) in done for w in c.windows for x in c.arms)   # noqa: E731
    todo = [g for g in G if not (tamam(g[0]) and tamam(g[1]))]
    print(f"[sonda] {len(G)} grup ({sum(len(c.arms) for g in G for c in g)} kol koşumu), kalan {len(todo)}; "
          f"gpu={a.gpu_jobs} cpu={a.jobs}; çıktı {out}")
    if a.kuru or not todo:
        print("[sonda] koşum BAŞLATILMADI" if a.kuru else "[sonda] hepsi tamam")
        return

    kilit = threading.Lock()
    t_bas = time.time()
    sayac = {"n": 0}
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=a.gpu_jobs, maxtasksperchild=8) as gp, \
            ctx.Pool(processes=a.jobs, maxtasksperchild=4) as cp, open(out, "a") as fout:
        def yaz(recs):
            with kilit:
                for r in recs:
                    fout.write(json.dumps(r, default=_js) + "\n")
                fout.flush()
                sayac["n"] += 1
                print(f"[sonda] {sayac['n']}/{2 * len(todo)} hücre yazıldı, {(time.time() - t_bas) / 60:.1f} dk",
                      flush=True)

        def hata(e):
            print(f"[sonda] İŞ HATASI: {type(e).__name__}: {e}", flush=True)

        bekleyen = []
        genis_of = {t.key(): g for t, g in todo}
        for meta in gp.imap(_frontend_worker, [(t, spool) for t, _ in todo]):
            c = meta["cell"]
            print(f"[sonda] ön-uç {'tamam' if meta.get('ok') else 'BAŞARISIZ'}: {c.key()} "
                  f"({meta.get('t_frontend', 0):.0f} s) {meta.get('error', '')}", flush=True)
            bekleyen.append(cp.apply_async(_backend_worker, (meta,), callback=yaz, error_callback=hata))
            if meta.get("ok"):
                bekleyen.append(cp.apply_async(_genis_is, ((genis_of[c.key()], spool),),
                                               callback=yaz, error_callback=hata))
        for b in bekleyen:
            b.wait()
    print(f"[sonda] bitti, {(time.time() - t_bas) / 60:.1f} dk")


if __name__ == "__main__":
    main()
