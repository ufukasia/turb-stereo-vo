#!/usr/bin/env python
"""Keşif sondası 2 (2026-11-11): kapı ayırma ve kestirimli kapı. Ön beyan: evidence/SONDA_KESIF/ON_BEYAN_2.md.

Parmak izli kodda değişiklik yok: saklı ön-uç dökümleri üstünde yalnız `backend_argv` değişir (kol A0). Bu dosya kod
parmak izi kümesinin (turb.runner.code_fingerprint) dışındadır. Ekran kartı kullanılmaz.

    python kapi_sonda.py --kuru                       # plan + bayraklar, koşmaz
    python kapi_sonda.py --dokumler D [--jobs 2]      # ön-test kümesi -> kapılar (B1, B2) -> kalanlar (sürdürülebilir)

Varyantlar: K (kontrol; ilk sondanın A0g hücresi), Gd (üç bayrak, DIMM sigma-hat düzeyi), Gn (yalnız gürültü grubu),
Gp (yalnız PnP eşiği). Çıktı: evidence/SONDA_KESIF/kapi2.jsonl (kayıt biçimi ızgarayla aynı; ek alan `kapi2`).
Çıkış kodu: 0 tamam; 2 kapı tutmadı; 3 bellek eşiği (koşum bırakıldı, yeniden başlatılabilir).
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import multiprocessing as mp
import os
import shutil
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

import sonda_kesif as sk  # noqa: E402
from turb.runner import _build_cfg, _backend_worker, _js, code_fingerprint  # noqa: E402
from turb.experiments import FAZ2_K_A0C  # noqa: E402

EV = os.path.join(ROOT, "evidence", "SONDA_KESIF")
OUT = os.path.join(EV, "kapi2.jsonl")
SONDA = os.path.join(EV, "records.jsonl")
ANA = os.path.join(ROOT, "evidence", "FAZ2", "records.jsonl")
SPOOL = "/tmp/turb_spool_kapi2"
DR = 8
BIT = ("m_ate_se3", "n_nees_pos_mean", "n_nees_rot_mean", "m_rpe_trans")
KONTROL_SEQ = ("V1_01_easy", "V2_01_easy")
ONTEST_SEQ = "V1_01_easy"
VARYANTLAR = ("K", "Gd", "Gn", "Gp")
CARPAN = math.hypot(1.0, 1.75)
PNP0 = 2.4477


def s_genis(sig: float) -> float:
    return math.sqrt(1.0 + (CARPAN * float(sig)) ** 2)


def sigma_hat_tablosu() -> dict:
    """Ana ızgara, tohum 0, R0, D/r0 8, kol A3: sekans -> sigma_hat_median [px] (sekans başına tek değer olmalı)."""
    out = {}
    with open(ANA) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            c = r.get("cell") or {}
            if (c.get("seed") == 0 and r.get("arm") == "A3" and c.get("D_over_r0") == float(DR)
                    and c.get("tau0_over_dt") == 1.0 and c.get("tag") == "faz2" and r.get("arm_stats")):
                out.setdefault(c["sequence"], set()).add(float(r["arm_stats"]["sigma_hat_median"]))
    for s, v in out.items():
        assert len(v) == 1 and all(math.isfinite(x) and x > 0 for x in v), (s, v)
    return {s: next(iter(v)) for s, v in out.items()}


def varyant_argv(v: str, seq: str, sh: dict) -> tuple:
    cal = ("--turb_cal_factor", str(float(FAZ2_K_A0C)))
    s_or = s_genis(sk.sigma_tilt(DR))
    if v == "K":
        return cal + sk.be_genis_argv(DR)
    if v == "Gd":
        s = s_genis(sh[seq])
        return cal + ("--px_noise", f"{s:.6f}", "--stereo_sigma", f"{s:.6f}", "--pnp_reproj", f"{PNP0 * s:.6f}")
    if v == "Gn":
        return cal + ("--px_noise", f"{s_or:.6f}", "--stereo_sigma", f"{s_or:.6f}")
    if v == "Gp":
        return cal + ("--pnp_reproj", f"{PNP0 * s_or:.6f}")
    raise ValueError(v)


def hucre(v: str, seq: str, argv: tuple):
    """İlk sondanın R0 / D/r0 8 / fe0 geniş hücresi; K'de birebir, ötekilerde yalnız backend_argv ve etiket değişir."""
    genis = [g for _, g in sk.gruplar([seq], [("R0", DR)]) if g.tag == "sonda_R0_fe0_bg"]
    assert len(genis) == 1, [g.tag for _, g in sk.gruplar([seq], [("R0", DR)])]
    c = genis[0]
    if v == "K":
        assert tuple(c.backend_argv) == tuple(argv), (c.backend_argv, argv)
        return c
    return dataclasses.replace(c, backend_argv=tuple(argv), tag=f"kapi2_{v}")


def dokum_yolu(D: str, seq: str) -> str:
    return os.path.join(D, f"{seq}__dr{DR}__s0.pkl")


def _isci(args):
    v, seq, argv, yol = args
    import config
    cell = hucre(v, seq, argv)
    d = os.path.join(SPOOL, cell.key().replace("/", "_") + f"__{v}")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d, exist_ok=True)
    cfg = _build_cfg(cell, d)
    dp = config.dump_path(cfg)
    shutil.copyfile(yol, dp)                # kopya: arka-uç işçisi dizini siler; saklı döküme dokunulmaz
    meta = dict(cell=cell, ok=True, dir=d, dump=dp, gt=config.gt_path(cfg), t_frontend=None, stats={}, gpu=None)
    t0 = time.time()
    recs = _backend_worker(meta)
    for r in recs:
        r["kapi2"] = {"varyant": v, "sekans": seq, "argv": list(argv), "sure_s": round(time.time() - t0, 1),
                      "dokum": os.path.basename(yol), "dokum_bayt": os.path.getsize(yol)}
    return recs


def oku_cikti() -> dict:
    """-> {(varyant, sekans): [kayıt, ...]} (bozuk satırlar atlanır; aynı pencere yinelenirse son kayıt geçerli)."""
    out = {}
    if os.path.exists(OUT):
        with open(OUT) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                k = r.get("kapi2") or {}
                if "varyant" in k:
                    out.setdefault((k["varyant"], k["sekans"]), {})[r.get("window")] = r
    return {k: list(v.values()) for k, v in out.items()}


def bitti(recs, n_pencere=3) -> bool:
    return len([r for r in recs if "arm" in r]) >= n_pencere


def sonda_basvuru() -> dict:
    """İlk sonda, R0 / D/r0 8 / fe0: {(kol, sekans, pencere): kayıt}; kol = A0, A1, A0g."""
    out = {}
    with open(SONDA) as f:
        for line in f:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            c = r.get("cell") or {}
            if "arm" not in r or c.get("D_over_r0") != float(DR) or abs(c.get("tau0_over_dt", 0) - 1.0) > 1e-9:
                continue
            tag = c.get("tag", "")
            if "_fe1" in tag:
                continue
            kol = "A0g" if tag.endswith("_bg") else r["arm"]
            out[(kol, c["sequence"], int(r["window"]))] = r
    return out


def kapilar(ontest) -> bool:
    """B1 (K bit-aynı) ve B2 (duman) — ön beyan §3. True: kalanlar koşulabilir."""
    C, ref, ok = oku_cikti(), sonda_basvuru(), True
    n_ayni = n_k = 0
    for v, seq in ontest:
        recs = C.get((v, seq), [])
        if not bitti(recs):
            print(f"[kapı] EKSİK {v} {seq}: {len(recs)} kayıt", flush=True)
            ok = False
            continue
        for r in recs:
            w = int(r["window"])
            if v == "K":
                n_k += 1
                m = ref.get(("A0g", seq, w)) or {}
                ayni = bool(m) and all(r.get(k) == m.get(k) for k in BIT)
                n_ayni += ayni
                if not ayni:
                    print(f"[kapı] B1 FARK {seq} w{w}: {[(k, r.get(k), m.get(k)) for k in BIT]}", flush=True)
                    ok = False
            dolu = (r.get("status") == "ok" and r.get("m_valid") and r.get("n_valid")
                    and abs(int(r.get("m_n_compare") or 0) - 400) <= 20)
            if not dolu:
                print(f"[kapı] B2 EKSİK ALAN {v} {seq} w{w}: status={r.get('status')} m_valid={r.get('m_valid')} "
                      f"n_valid={r.get('n_valid')} n_compare={r.get('m_n_compare')} err={r.get('error')}", flush=True)
                ok = False
            if v != "K":
                a0 = ref.get(("A0", seq, w)) or {}
                if not a0 or r.get("m_ate_se3") == a0.get("m_ate_se3"):
                    print(f"[kapı] B2 ETKİ YOK {v} {seq} w{w}: ATE {r.get('m_ate_se3')} = A0 {a0.get('m_ate_se3')}",
                          flush=True)
                    ok = False
    print(f"[kapı] B1: K kayıtları ilk sondanın A0g kayıtlarıyla bit-aynı {n_ayni}/{n_k}; B2 + B1 sonucu: "
          f"{'GEÇTİ' if ok else 'KALDI'}", flush=True)
    return ok


def bellek_mb() -> float:
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024.0
    return float("inf")


def kos(isler, jobs, bellek_alt, etiket):
    """isler: [(v, seq, argv, yol)]. Bittikçe yazar. Bellek eşiğin altına inerse havuzu bırakır (çıkış 3)."""
    if not isler:
        return
    t_bas, sayac, hatalar = time.time(), {"n": 0}, []
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=jobs, maxtasksperchild=1) as cp, open(OUT, "a") as fout:
        def yaz(recs):
            for r in recs:
                fout.write(json.dumps(r, default=_js) + "\n")
            fout.flush()
            sayac["n"] += 1
            k = (recs[0].get("kapi2") or {}) if recs else {}
            print(f"[kapi2] {etiket} {sayac['n']}/{len(isler)} {time.strftime('%T')} {(time.time() - t_bas) / 60:.1f} dk: "
                  f"{k.get('varyant')} {k.get('sekans')} {[(r.get('window'), r.get('status')) for r in recs]} "
                  f"{k.get('sure_s')} s", flush=True)

        def hata(e):
            hatalar.append(repr(e))
            print(f"[kapi2] İŞÇİ HATASI: {e!r}", flush=True)

        bek = [cp.apply_async(_isci, (i,), callback=yaz, error_callback=hata) for i in isler]
        while not all(b.ready() for b in bek):
            time.sleep(5.0)
            m = bellek_mb()
            if m < bellek_alt:
                print(f"[kapi2] BELLEK {m:.0f} MB < {bellek_alt} MB: koşum bırakılıyor (yeniden başlatılabilir)", flush=True)
                cp.terminate()
                fout.flush()
                sys.exit(3)
    if hatalar:
        sys.exit(f"[kapi2] DUR: {len(hatalar)} işçi hatası")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kuru", action="store_true")
    ap.add_argument("--dokumler", default=None)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--bellek_alt", type=float, default=3600.0, help="[MB] MemAvailable bunun altına inerse bırak")
    a = ap.parse_args()

    fp = code_fingerprint()
    print(f"[kapi2] code_fp {fp} (beklenen {sk.BEKLENEN_FP}); çıktı {OUT}", flush=True)
    if fp != sk.BEKLENEN_FP:
        sys.exit("[kapi2] DUR: kod parmak izi ana ızgaradan farklı")
    sh = sigma_hat_tablosu()
    assert sorted(sh) == sorted(sk.SEQS), sorted(sh)
    s_or = s_genis(sk.sigma_tilt(DR))
    print(f"[kapi2] oracle sigma {sk.sigma_tilt(DR):.4f} px -> s {s_or:.6f}; pnp {PNP0 * s_or:.6f}", flush=True)

    ontest = [("K", s) for s in KONTROL_SEQ] + [(v, ONTEST_SEQ) for v in ("Gd", "Gn", "Gp")]
    kalan = [(v, s) for v in ("Gd", "Gn", "Gp") for s in sk.SEQS if (v, s) not in ontest]
    for ad, P in (("ön-test", ontest), ("kalan", kalan)):
        print(f"[kapi2] {ad}: {len(P)} koşum")
        for v, s in P:
            argv = varyant_argv(v, s, sh)
            c = hucre(v, s, argv)
            ek = f" sigma_hat {sh[s]:.4f} ({sh[s] / sk.sigma_tilt(DR):.3f} x oracle)" if v == "Gd" else ""
            print(f"    {v:2s} {s:13s} {' '.join(argv)}{ek}  [{c.tag}; pencereler {tuple(c.windows)}]")
    if a.kuru:
        return
    if not a.dokumler:
        sys.exit("[kapi2] --dokumler gerekli")
    for s in sk.SEQS:
        assert os.path.getsize(dokum_yolu(a.dokumler, s)) > 50e6, dokum_yolu(a.dokumler, s)
    os.makedirs(SPOOL, exist_ok=True)

    def isler(P):
        C = oku_cikti()
        return [(v, s, varyant_argv(v, s, sh), dokum_yolu(a.dokumler, s)) for v, s in P if not bitti(C.get((v, s), []))]

    kos(isler(ontest), a.jobs, a.bellek_alt, "ön-test")
    if not kapilar(ontest):
        sys.exit(2)
    kos(isler(kalan), a.jobs, a.bellek_alt, "kalan")
    C = oku_cikti()
    eksik = [(v, s) for v, s in ontest + kalan if not bitti(C.get((v, s), []))]
    print(f"[kapi2] bitti {time.strftime('%F %T')}: {len(ontest) + len(kalan) - len(eksik)}/{len(ontest) + len(kalan)} koşum tam; "
          f"eksik {eksik}", flush=True)


if __name__ == "__main__":
    main()
