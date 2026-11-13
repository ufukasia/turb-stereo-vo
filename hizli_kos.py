#!/usr/bin/env python
"""Hızlı rejim doğrulayıcı ızgarası (ON_KAYIT_FAZ2.md §5 Sapma 5).

Ana ızgaranın bütün hücreleri yarı-statik rejimdeydi (tau0 = dt). Bu ızgara aynı hücre tanımını sabit enine
hız rejiminde koşar: tau0/dt = 0.314 * (D / (D/r0)) / (v * dt). Kod değişikliği yok (yalnız Cell.tau0_over_dt);
bu dosya kod parmak izi kümesinin (turb.runner.code_fingerprint) dışındadır. Hücre tanımı ana ızgaranınkiyle
(turb.experiments.cells_FAZ2) aynıdır; yalnız tau0_over_dt, kol kümesi ve etiket farklı.

    export CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=4
    python hizli_kos.py --kuru                                   # plan + birim denetimleri, koşmaz
    python hizli_kos.py --ontest --gpu_jobs 1 --jobs 2           # ön test (evidence/HIZLI/ON_TEST.md)
    python hizli_kos.py --gpu_jobs 3 --jobs 10                   # ızgara (devam ettirilebilir)

Koşum sırasında her biten hücre denetlenir:
  (a) bütün kayıtlarda metrik alanları dolu (durum "ok" olanlarda),
  (b) keşif sondasında (tohum 0, fe0; kollar A0 ve A1) ya da ön testte (evidence/HIZLI/on_test.jsonl) karşılığı olan
      hücrelerde kayıt başvuruyla bit-aynı ve ön-uç istatistikleri aynı.
(b) tutmazsa yeni hücre gönderilmez (DUR) ve çıkış kodu 2 olur.
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
                              FAZ2_K_A0C, FAZ2_SEQS, FAZ2_ARMS)

BEKLENEN_FP = "08e28b83407c"
D_AP, DT = 0.10, 0.05                      # açıklık [m], kare aralığı [s] (sonda ile aynı)
REJIM_HIZ = {"R5": 5.0, "R1": 1.0}         # rejim adı -> enine hız [m/s]
REJIMLER = ("R5",)
DUZEYLER = (8, 4)
TOHUMLAR = (1, 2, 3, 4)                    # sondanın ve onarım ızgarasının görmediği tohumlar
KOLLAR = tuple(FAZ2_ARMS)                  # ana ızgaranın kol kümesi (A5 hariç)
SEQS = ["MH01", "V1_01_easy", "MH02", "V1_02_medium", "MH03", "V2_01_easy", "MH04", "V2_02_medium", "MH05"]
assert sorted(SEQS) == sorted(FAZ2_SEQS)
FE_SURE_ALANI = "total_execution_time_seconds"
BIT_ALANLARI = ("m_ate_se3", "n_nees_pos_mean", "n_nees_rot_mean", "m_rpe_trans")
DOLU_ALANLAR = ("m_ate_se3", "m_rpe_trans", "n_nees_pos_mean", "n_nees_rot_mean")
EV = os.path.join(ROOT, "evidence", "HIZLI")
SONDA = os.path.join(ROOT, "evidence", "SONDA_KESIF", "records.jsonl")


def tau0_over_dt(rejim: str, dr: float) -> float:
    return 0.314 * (D_AP / float(dr)) / (REJIM_HIZ[rejim] * DT)


def enine_hiz(t0: float, dr: float) -> float:
    return 0.314 * (D_AP / float(dr)) / (t0 * DT)


def _hucre(s, sd, rejim, dr, arms, tag=None):
    return Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence=s, seed=int(sd), D_over_r0=float(dr),
                tau0_over_dt=tau0_over_dt(rejim, dr), arms=tuple(arms), windows=windows_for(s, FAZ2_WINDOWS_FILE),
                window_len=400, nees=True, cov_every=FAZ2_COV_EVERY, psf_legacy_scale=False,
                sensor_noise_dn=FAZ2_NOISE_DN, noise_seed=None,
                backend_argv=("--turb_cal_factor", str(float(FAZ2_K_A0C))), tag=tag or f"hizli_{rejim}",
                cache_dump=False)


def hucreler(rejimler, duzeyler, tohumlar, seqs, kollar):
    """Sıra: rejim -> düzey -> tohum -> sekans (her düzeyin ilk tohumu bütün sekanslarda en başta biter)."""
    return [_hucre(s, sd, rj, dr, kollar) for rj in rejimler for dr in duzeyler for sd in tohumlar for s in seqs]


def ontest_hucreleri(kollar):
    """V1_01_easy: (1) tohum 0, R5, D/r0 8, kollar A0 + A1 -> sonda kaydıyla bit-aynı olmalı (gerileme kapısı);
    (2) tohum 1, R5, D/r0 8 ve 4, ızgaranın kol kümesiyle (duman + süre ölçümü)."""
    return [_hucre("V1_01_easy", 0, "R5", 8, ("A0", "A1"), "hizli_ontest_R5"),
            _hucre("V1_01_easy", 1, "R5", 8, kollar, "hizli_ontest_R5"),
            _hucre("V1_01_easy", 1, "R5", 4, kollar, "hizli_ontest_R5")]


def _basvuru_oku(path, kaynak, etiket_sonu, ref):
    if not os.path.exists(path):
        return ref
    with open(path) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            c = r.get("cell") or {}
            if "arm" not in r or not str(c.get("tag", "")).endswith(etiket_sonu):
                continue
            k = (c["sequence"], int(c["seed"]), float(c["D_over_r0"]), round(float(c["tau0_over_dt"]), 6))
            if ref.get(k, {}).get("kaynak", kaynak) != kaynak:
                continue                    # aynı hücre için ilk kaynak (sonda) geçerli
            d = ref.setdefault(k, {"kaynak": kaynak})
            d.setdefault("fe", r.get("frontend_stats"))
            d[(r["arm"], int(r["window"]))] = {x: r.get(x) for x in BIT_ALANLARI + ("status",)}
    return ref


def sonda_kayitlari():
    """(sekans, tohum, D/r0, tau0/dt 6 hane) -> {"kaynak", "fe": frontend_stats, (kol, pencere): alanlar}; fe0 taban."""
    return _basvuru_oku(SONDA, "sonda", "_fe0", {})


def basvuru_kayitlari(ontest_dahil):
    """Sonda kayıtları; ızgara koşumunda ayrıca ön testin kayıtları (aynı hücre yeniden koşulunca bit-aynı olmalı)."""
    ref = sonda_kayitlari()
    if ontest_dahil:
        _basvuru_oku(os.path.join(EV, "on_test.jsonl"), "ön test", "hizli_ontest_R5", ref)
    return ref


def fe_ayni(a, b):
    if not a or not b:
        return False
    return ({k: v for k, v in a.items() if k != FE_SURE_ALANI}
            == {k: v for k, v in b.items() if k != FE_SURE_ALANI})


def denetle(recs, ref):
    """-> (özet metni, sorun listesi). Sorun = DUR gerektiren durum."""
    sorun = []
    kayit = [r for r in recs if "arm" in r]
    if not kayit:
        return "hücre hatası: " + str(recs[0].get("error")), ["hücre hatası"]
    c = kayit[0]["cell"]
    n_ok = sum(1 for r in kayit if r.get("status") == "ok")
    bos = sum(1 for r in kayit if r.get("status") == "ok" and any(r.get(k) is None for k in DOLU_ALANLAR))
    if bos:
        sorun.append(f"{bos} 'ok' kaydında boş metrik alanı")
    m = ref.get((c["sequence"], int(c["seed"]), float(c["D_over_r0"]), round(float(c["tau0_over_dt"]), 6)))
    ek = ""
    if m is not None:
        fe_ok = fe_ayni(kayit[0].get("frontend_stats"), m.get("fe"))
        if not fe_ok:
            sorun.append(f"ön-uç istatistikleri başvurudan ({m.get('kaynak', 'sonda')}) farklı")
        n_bit = n_bit_ok = 0
        for r in kayit:
            s = m.get((r["arm"], int(r["window"])))
            if s is None:
                continue
            n_bit += 1
            if all(r.get(k) == s.get(k) for k in BIT_ALANLARI):
                n_bit_ok += 1
            else:
                sorun.append(f"{r['arm']} pencere {r['window']} başvuruyla ({m.get('kaynak', 'sonda')}) bit-aynı değil "
                             f"({r.get('m_ate_se3')} / {s.get('m_ate_se3')})")
        ek = f", {m.get('kaynak', 'sonda')}: ön-uç {'aynı' if fe_ok else 'FARKLI'}, bit-aynı {n_bit_ok}/{n_bit}"
    durumlar = {}
    for r in kayit:
        durumlar[r.get("status")] = durumlar.get(r.get("status"), 0) + 1
    return f"{n_ok}/{len(kayit)} ok {durumlar if n_ok != len(kayit) else ''}{ek}", sorun


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kuru", action="store_true")
    ap.add_argument("--ontest", action="store_true")
    ap.add_argument("--rejimler", default=",".join(REJIMLER))
    ap.add_argument("--duzeyler", default=",".join(str(d) for d in DUZEYLER))
    ap.add_argument("--tohumlar", default=",".join(str(t) for t in TOHUMLAR))
    ap.add_argument("--kollar", default=",".join(KOLLAR))
    ap.add_argument("--yalniz", default="", help="yalnız bu sekanslar (virgüllü)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--gpu_jobs", type=int, default=2)
    a = ap.parse_args()

    os.makedirs(EV, exist_ok=True)
    kollar = tuple(x for x in a.kollar.split(",") if x)
    if a.ontest:
        C = ontest_hucreleri(kollar)
        out, spool = os.path.join(EV, "on_test.jsonl"), "/tmp/turb_spool_hizli_ontest"
    else:
        rej = [x for x in a.rejimler.split(",") if x]
        duz = [int(x) for x in a.duzeyler.split(",") if x != ""]
        toh = [int(x) for x in a.tohumlar.split(",") if x != ""]
        sec = [s for s in a.yalniz.split(",") if s]
        assert all(s in SEQS for s in sec), "bilinmeyen sekans"
        assert all(r in REJIM_HIZ for r in rej), "bilinmeyen rejim"
        C = hucreler(rej, duz, toh, [s for s in SEQS if not sec or s in sec], kollar)
        out, spool = os.path.join(EV, "records.jsonl"), "/tmp/turb_spool_hizli"
    os.makedirs(spool, exist_ok=True)

    fp = code_fingerprint()
    env = {k: os.environ.get(k) for k in ("CUDA_DEVICE_ORDER", "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS")}
    print(f"[hizli] code_fp {fp} (beklenen {BEKLENEN_FP}); ortam {env}; spool {spool}")
    if fp != BEKLENEN_FP:
        sys.exit("[hizli] DUR: kod parmak izi ana ızgaradan farklı")
    if env["CUDA_DEVICE_ORDER"] != "PCI_BUS_ID" or env["CUDA_VISIBLE_DEVICES"] != "1":
        sys.exit("[hizli] DUR: yalnız RTX 3090 (CUDA_DEVICE_ORDER=PCI_BUS_ID, CUDA_VISIBLE_DEVICES=1)")
    anahtar = [c.key() for c in C]
    assert len(set(anahtar)) == len(anahtar), "anahtar çakışması"
    for t0, dr in sorted(set((c.tau0_over_dt, c.D_over_r0) for c in C)):
        print(f"  D/r0={dr:4.0f}: tau0/dt={t0:.5f}  enine hız v={enine_hiz(t0, dr):.3f} m/s")
    ref = basvuru_kayitlari(ontest_dahil=not a.ontest)
    n_ref = sum(1 for c in C if (c.sequence, c.seed, c.D_over_r0, round(c.tau0_over_dt, 6)) in ref)

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
    print(f"[hizli] {len(C)} hücre ({sum(len(c.arms) for c in C)} kol koşumu, "
          f"{sum(len(c.arms) * len(c.windows) for c in C)} kayıt), kalan {len(todo)}; başvurusu (sonda / ön test) olan {n_ref}; "
          f"gpu={a.gpu_jobs} cpu={a.jobs}; çıktı {out}")
    if a.kuru or not todo:
        for c in C[:10] + C[-3:]:
            print(f"   {c.key()}  kollar={c.arms}  argv={c.backend_argv}  pencereler={c.windows}")
        print("[hizli] koşum BAŞLATILMADI" if a.kuru else "[hizli] hepsi tamam")
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
                ozet, sorun = denetle(recs, ref)
                el = (time.time() - t_bas) / 60.0
                print(f"[hizli] {sayac['n']}/{len(todo)}  {recs[0].get('key', '?')}  {ozet}  "
                      f"{el:.1f} dk  kalan ~{el / sayac['n'] * (len(todo) - sayac['n']):.0f} dk", flush=True)
                if sorun:
                    sayac["sorun"] += 1
                    for s in sorun[:6]:
                        print(f"[hizli] DUR: {recs[0].get('key', '?')}: {s}", flush=True)
                    dur.set()
            yer.release()

        def hata(e):
            print(f"[hizli] İŞ HATASI: {type(e).__name__}: {e}", flush=True)
            dur.set()
            yer.release()

        bekleyen = []
        for meta in gp.imap(_frontend_worker, uretec()):
            c = meta["cell"]
            print(f"[hizli] ön-uç {'tamam' if meta.get('ok') else 'BAŞARISIZ'}: {c.key()} "
                  f"({meta.get('t_frontend', 0):.0f} s) {meta.get('error', '')}", flush=True)
            bekleyen.append(cp.apply_async(_backend_worker, (meta,), callback=yaz, error_callback=hata))
        for b in bekleyen:
            b.wait()
    print(f"[hizli] bitti: {sayac['n']} hücre, sorunlu {sayac['sorun']}, {(time.time() - t_bas) / 60:.1f} dk")
    if sayac["sorun"] or dur.is_set():
        sys.exit(2)


if __name__ == "__main__":
    main()
