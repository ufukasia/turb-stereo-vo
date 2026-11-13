#!/usr/bin/env python
"""A3S ayırma kolu (ON_KAYIT_FAZ2.md §5 Sapma 4): A3'ün marjinal varyansı, yama korelasyonu olmadan.

Kol: mevcut "A0C+A3" kolu, `--turb_cal_factor 1.0` ile (sigma^2 = sigma0^2 + s_DIMM^2; yama faktörü yok).
Kod değişikliği yok; bu dosya kod parmak izi kümesinin (turb.runner.code_fingerprint) dışındadır.
Hücre tanımı ana ızgaranınkiyle (turb.experiments.cells_FAZ2) aynıdır; yalnız kollar, cal_factor ve etiket farklı.

    export CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=4
    python a3s_kos.py --kuru                          # plan + denetimler, koşmaz
    python a3s_kos.py --ontest --gpu_jobs 2 --jobs 4  # ön test (evidence/A3S/ON_TEST.md)
    python a3s_kos.py --gpu_jobs 3 --jobs 8           # ızgara: D/r0 8, 16, 24, 12, 20 (devam ettirilebilir)
    python a3s_kos.py --duzeyler 4,0 ...              # betimsel ek düzeyler (yalnız süre kalırsa)

Koşum sırasında her biten hücre ana ızgarayla karşılaştırılır (Sapma 4, denetim a ve b):
  (a) ön-uç istatistikleri ana ızgara kaydıyla aynı (süre alanı hariç),
  (b) A0 koşulan hücrelerde A0 kaydı ana ızgaranın A0 kaydıyla bit-aynı.
Biri tutmazsa yeni hücre gönderilmez (DUR) ve çıkış kodu 2 olur.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import threading
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from turb.runner import Cell, windows_for, _frontend_worker, _backend_worker, _js, code_fingerprint  # noqa: E402
from turb.experiments import (EUROC_FAZ2, FAZ2_WINDOWS_FILE, FAZ2_NOISE_DN, FAZ2_COV_EVERY,  # noqa: E402
                              FAZ2_K_A0C, FAZ2_SEQS)

BEKLENEN_FP = "08e28b83407c"
KOL = "A0C+A3"
K_A3S = 1.0
DUZEYLER = (8, 16, 24, 12, 20)            # Sapma 4 koşum sırası
# uzun/kısa sekanslar dönüşümlü (yük dengesi); küme FAZ2_SEQS ile aynı
SEQS = ["MH01", "V1_01_easy", "MH02", "V1_02_medium", "MH03", "V2_01_easy", "MH04", "V2_02_medium", "MH05"]
assert sorted(SEQS) == sorted(FAZ2_SEQS)
FE_SURE_ALANI = "total_execution_time_seconds"      # ölçüm (duvar saati), çıktı değil: eşitlik denetiminin dışında
BIT_ALANLARI = ("m_ate_se3", "n_nees_pos_mean", "n_nees_rot_mean", "m_rpe_trans")


def _hucre(s, sd, dr, arms, k, tag):
    return Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence=s, seed=int(sd), D_over_r0=float(dr),
                arms=tuple(arms), windows=windows_for(s, FAZ2_WINDOWS_FILE), window_len=400, nees=True,
                cov_every=FAZ2_COV_EVERY, psf_legacy_scale=False, sensor_noise_dn=FAZ2_NOISE_DN,
                noise_seed=None, backend_argv=("--turb_cal_factor", str(float(k))), tag=tag, cache_dump=False)


def hucreler(duzeyler, tohumlar, seqs):
    """Düzey sırası korunur; düzey içinde tohum-önce (D/r0 = 8, tohum 0 denetim hücreleri en başta)."""
    out = []
    for dr in duzeyler:
        for sd in tohumlar:
            for s in seqs:
                arms = ("A0", KOL) if (int(dr) == 8 and int(sd) == 0) else (KOL,)
                out.append(_hucre(s, sd, dr, arms, K_A3S, "a3s"))
    return out


def ontest_hucreleri():
    """V1_01_easy, tohum 0, D/r0 {8, 0}: denetim hücresi (k = 12: ana ızgaranın A0 ve A0C+A3 kayıtlarıyla
    bit-aynı olmalı) + A3S hücresi (ızgaradaki kol kümesiyle)."""
    out = []
    for dr in (8, 0):
        out.append(_hucre("V1_01_easy", 0, dr, ("A0", KOL), FAZ2_K_A0C, "a3s_ontest_k12"))
        out.append(_hucre("V1_01_easy", 0, dr, ("A0", KOL) if dr == 8 else (KOL,), K_A3S, "a3s_ontest"))
    return out


def ana_izgara():
    """(sekans, tohum, D/r0) -> {"fe": frontend_stats, (kol, pencere): kayıt} (yalnız karşılaştırılan kollar)."""
    ana = {}
    with open(os.path.join(ROOT, "evidence", "FAZ2", "records.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            c = r.get("cell") or {}
            if "arm" not in r or c.get("tag") != "faz2":
                continue
            d = ana.setdefault((c["sequence"], int(c["seed"]), float(c["D_over_r0"])), {})
            d.setdefault("fe", r.get("frontend_stats"))
            if r["arm"] in ("A0", KOL):
                d[(r["arm"], int(r["window"]))] = {k: r.get(k) for k in BIT_ALANLARI + ("status",)}
    return ana


def fe_ayni(a, b):
    if not a or not b:
        return False
    ka = {k: v for k, v in a.items() if k != FE_SURE_ALANI}
    kb = {k: v for k, v in b.items() if k != FE_SURE_ALANI}
    return ka == kb


def denetle(recs, ana):
    """Bir hücrenin kayıtlarını ana ızgarayla karşılaştırır -> (özet metni, sorun listesi)."""
    sorun, n_bit, n_bit_ok = [], 0, 0
    kayit = [r for r in recs if "arm" in r]
    if not kayit:
        return "hücre hatası: " + str(recs[0].get("error")), ["hücre hatası"]
    c = kayit[0]["cell"]
    ref = ana.get((c["sequence"], int(c["seed"]), float(c["D_over_r0"])))
    if ref is None:
        return "ana ızgarada karşılığı yok", ["ana ızgarada karşılığı yok"]
    fe_ok = fe_ayni(kayit[0].get("frontend_stats"), ref.get("fe"))
    if not fe_ok:
        sorun.append("ön-uç istatistikleri ana ızgaradan farklı")
    k12 = float(c["backend_argv"][1]) == float(FAZ2_K_A0C)
    for r in kayit:
        if r["arm"] == "A0" or (r["arm"] == KOL and k12):
            m = ref.get((r["arm"], int(r["window"])))
            n_bit += 1
            if m is not None and all(r.get(k) == m.get(k) for k in BIT_ALANLARI):
                n_bit_ok += 1
            else:
                sorun.append(f"{r['arm']} pencere {r['window']} ana ızgarayla bit-aynı değil "
                             f"({r.get('m_ate_se3')} / {None if m is None else m.get('m_ate_se3')})")
    n_ok = sum(1 for r in kayit if r.get("status") == "ok")
    return (f"{n_ok}/{len(kayit)} ok, ön-uç {'aynı' if fe_ok else 'FARKLI'}, "
            f"bit-aynı {n_bit_ok}/{n_bit}"), sorun


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kuru", action="store_true")
    ap.add_argument("--ontest", action="store_true")
    ap.add_argument("--duzeyler", default=",".join(str(d) for d in DUZEYLER))
    ap.add_argument("--tohumlar", default="0,1,2,3,4")
    ap.add_argument("--yalniz", default="", help="yalnız bu sekanslar (virgüllü)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--gpu_jobs", type=int, default=3)
    a = ap.parse_args()

    ev = os.path.join(ROOT, "evidence", "A3S")
    os.makedirs(ev, exist_ok=True)
    if a.ontest:
        C = ontest_hucreleri()
        out, spool = os.path.join(ev, "on_test.jsonl"), "/tmp/turb_spool_a3s_ontest"
    else:
        duz = [int(x) for x in a.duzeyler.split(",") if x != ""]
        toh = [int(x) for x in a.tohumlar.split(",") if x != ""]
        sec = [s for s in a.yalniz.split(",") if s]
        assert all(s in SEQS for s in sec), "bilinmeyen sekans"
        C = hucreler(duz, toh, [s for s in SEQS if not sec or s in sec])
        out, spool = os.path.join(ev, "records.jsonl"), "/tmp/turb_spool_a3s"
    os.makedirs(spool, exist_ok=True)

    fp = code_fingerprint()
    env = {k: os.environ.get(k) for k in ("CUDA_DEVICE_ORDER", "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS")}
    print(f"[a3s] code_fp {fp} (beklenen {BEKLENEN_FP}); ortam {env}; spool {spool}")
    if fp != BEKLENEN_FP:
        sys.exit("[a3s] DUR: kod parmak izi ana ızgaradan farklı")
    if env["CUDA_DEVICE_ORDER"] != "PCI_BUS_ID" or env["CUDA_VISIBLE_DEVICES"] != "1":
        sys.exit("[a3s] DUR: yalnız RTX 3090 (CUDA_DEVICE_ORDER=PCI_BUS_ID, CUDA_VISIBLE_DEVICES=1)")
    anahtar = [c.key() for c in C]
    assert len(set(anahtar)) == len(anahtar), "anahtar çakışması"
    ana = ana_izgara()
    eksik = [c.key() for c in C if (c.sequence, c.seed, c.D_over_r0) not in ana]
    assert not eksik, f"ana ızgarada karşılığı olmayan hücre: {eksik[:3]}"

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
    todo = [c for c in C if not all((c.key(), int(w), x) in done for w in c.windows for x in c.arms)]
    n_kol = sum(len(c.arms) for c in C)
    print(f"[a3s] {len(C)} hücre ({n_kol} kol koşumu, {sum(len(c.arms) * len(c.windows) for c in C)} kayıt), "
          f"kalan {len(todo)}; gpu={a.gpu_jobs} cpu={a.jobs}; çıktı {out}")
    if a.kuru or not todo:
        for c in C[:12] + C[-3:]:
            print(f"   {c.key()}  kollar={c.arms}  argv={c.backend_argv}  pencereler={c.windows}")
        print("[a3s] koşum BAŞLATILMADI" if a.kuru else "[a3s] hepsi tamam")
        return

    kilit = threading.Lock()
    dur = threading.Event()
    yer = threading.Semaphore(a.jobs + a.gpu_jobs + 2)      # geri basınç: spool'da bekleyen döküm sayısı sınırlı
    t_bas = time.time()
    sayac = {"n": 0, "sorun": 0}

    def uretec():
        for c in todo:
            yer.acquire()
            if dur.is_set():
                return
            yield (c, spool)

    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=a.gpu_jobs, maxtasksperchild=8) as gp, \
            ctx.Pool(processes=a.jobs, maxtasksperchild=4) as cp, open(out, "a") as fout:
        def yaz(recs):
            with kilit:
                for r in recs:
                    fout.write(json.dumps(r, default=_js) + "\n")
                fout.flush()
                sayac["n"] += 1
                ozet, sorun = denetle(recs, ana)
                el = (time.time() - t_bas) / 60.0
                print(f"[a3s] {sayac['n']}/{len(todo)}  {recs[0].get('key', '?')}  {ozet}  "
                      f"{el:.1f} dk  kalan ~{el / sayac['n'] * (len(todo) - sayac['n']):.0f} dk", flush=True)
                if sorun:
                    sayac["sorun"] += 1
                    for s in sorun[:6]:
                        print(f"[a3s] DUR: {recs[0].get('key', '?')}: {s}", flush=True)
                    dur.set()
            yer.release()

        def hata(e):
            print(f"[a3s] İŞ HATASI: {type(e).__name__}: {e}", flush=True)
            dur.set()
            yer.release()

        bekleyen = []
        for meta in gp.imap(_frontend_worker, uretec()):
            c = meta["cell"]
            print(f"[a3s] ön-uç {'tamam' if meta.get('ok') else 'BAŞARISIZ'}: {c.key()} "
                  f"({meta.get('t_frontend', 0):.0f} s) {meta.get('error', '')}", flush=True)
            bekleyen.append(cp.apply_async(_backend_worker, (meta,), callback=yaz, error_callback=hata))
        for b in bekleyen:
            b.wait()
    print(f"[a3s] bitti: {sayac['n']} hücre, sorunlu {sayac['sorun']}, {(time.time() - t_bas) / 60:.1f} dk")
    if sayac["sorun"] or dur.is_set():
        sys.exit(2)


if __name__ == "__main__":
    main()
