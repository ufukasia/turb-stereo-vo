#!/usr/bin/env python
"""Uyarlanır izleme kapısı kolu (`--turb_gate`): ön testler (T1-T3) ve geliştirme koşumu (KEŞİF; tohum 0).

Ön beyan: evidence/KAPI/ON_BEYAN_3.md (koşumdan önce yazılır; sha256 evidence/KAPI/kimlik.txt içinde).
Bu ağaç (pyeco_turb_kapi) geliştirme kopyasıdır; koşan ana ağaca (pyeco_turb, code_fp 08e28b83407c) dokunulmaz.
Saklı ön-uç dökümleri (tohum 0, R0, D/r0 8) üstünde yalnız arka-uç koşulur; ekran kartı kullanılmaz.

    python kapi_kol.py --kuru
    python kapi_kol.py --asama ontest --dokumler D [--jobs 5]     # T1 (A3 hariç), T2, T3 -> kapılar (çıkış 2: kaldı)
    python kapi_kol.py --asama a3 --dokumler D                     # T1'in A3 kolu (yavaş; ayrı dosya)
    python kapi_kol.py --asama gelistirme --dokumler D [--jobs 6]  # D1, D2: 9 sekans (ön test geçmişse)

Varyantlar (kol, arka-uç bayrakları):
  T1 gerileme (kapı kapalı; eski kayıtlarla bit-aynı olmalı): R_A0, R_A0C, R_A1, R_A3S, R_A3
  T2 eşdeğerlik: E_A0Cg (A0C k = s, gate 3)  ==  E_A0g (A0, --px_noise s --stereo_sigma s --pnp_reproj 2.4477 s)
                 E_bir3 (A0C k = 1, gate 3)  ==  E_bir0 (A0C k = 1, kapı kapalı)
  T3 duman + geliştirme: D1 (A3S + gate 3, c = 1), D2 (c = hypot(1, 1.75)), D1n (gate 1), D1p (gate 2),
                 A1g (A1 + gate 3; oracle üst sınır)
Çıkış kodu: 0 tamam; 2 kapı tutmadı; 3 bellek eşiği (yeniden başlatılabilir).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import shutil
import sys
import time
from dataclasses import asdict

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from turb.runner import Cell, windows_for, _build_cfg, _run_arms, _js, code_fingerprint  # noqa: E402
from turb.experiments import (EUROC_FAZ2, FAZ2_WINDOWS_FILE, FAZ2_NOISE_DN, FAZ2_COV_EVERY,  # noqa: E402
                              FAZ2_K_A0C, FAZ2_SEQS)
from turb.backend import ArmConfig  # noqa: E402

ANA = os.path.join(os.path.dirname(ROOT), "pyeco_turb", "evidence")      # başvuru kayıtları (salt okunur)
ANA_FP = "08e28b83407c"
EV = os.path.join(ROOT, "evidence", "KAPI")
CIKTI = {"ontest": os.path.join(EV, "on_test.jsonl"), "a3": os.path.join(EV, "on_test_a3.jsonl"),
         "gelistirme": os.path.join(EV, "gelistirme.jsonl")}
SPOOL = "/tmp/turb_spool_kapi_kol"
DR = 8
SEQS = ["MH01", "V1_01_easy", "MH02", "V1_02_medium", "MH03", "V2_01_easy", "MH04", "V2_02_medium", "MH05"]
assert sorted(SEQS) == sorted(FAZ2_SEQS)
ONTEST_SEQ = "V1_01_easy"
BIT = ("m_ate_se3", "n_nees_pos_mean", "n_nees_rot_mean", "m_rpe_trans")
CARPAN = math.hypot(1.0, 1.75)
PNP0, PX0, ST0 = 2.4477, 1.0, 1.0          # config.yaml tabanı (kuru aşamada doğrulanır)
SIGMA_TABAN = 0.05                          # ArmConfig.sigma_floor_px (ilk DIMM kestiriminden önce)


def sigma_tilt(dr: float) -> float:
    return float(ArmConfig(D_over_r0=float(dr)).sigma_px_oracle)


S_TAM = math.sqrt(1.0 + (CARPAN * sigma_tilt(DR)) ** 2)      # A1'in bireysel-faktör sigması (yuvarlanmamış)
S6 = f"{S_TAM:.6f}"                                           # ilk sondanın A0g bayrağındaki yuvarlanmış değer
_CAL12 = ("--turb_cal_factor", str(float(FAZ2_K_A0C)))
_CAL1 = ("--turb_cal_factor", "1.0")


def _d(c: float, gate: int) -> tuple:
    return _CAL1 + ("--turb_dimm_sigma_scale", repr(float(c)), "--turb_gate", str(int(gate)))


VARYANT = {
    # T1: kapı kapalı, eski kodun kayıtlarıyla bit-aynı
    "R_A0": ("A0", _CAL12), "R_A0C": ("A0C", _CAL12), "R_A1": ("A1", _CAL12),
    "R_A3S": ("A0C+A3", _CAL1), "R_A3": ("A3", _CAL12),
    # T2: eşdeğerlik
    "E_A0g": ("A0", _CAL12 + ("--px_noise", S6, "--stereo_sigma", S6, "--pnp_reproj", repr(PNP0 * float(S6)))),
    "E_A0Cg": ("A0C", ("--turb_cal_factor", S6, "--turb_gate", "3")),
    "E_bir0": ("A0C", _CAL1), "E_bir3": ("A0C", _CAL1 + ("--turb_gate", "3")),
    # T3 + geliştirme
    "D1": ("A0C+A3", _d(1.0, 3)), "D2": ("A0C+A3", _d(CARPAN, 3)),
    "D1n": ("A0C+A3", _d(1.0, 1)), "D1p": ("A0C+A3", _d(1.0, 2)),
    "A1g": ("A1", _CAL12 + ("--turb_gate", "3")),
}
C_DEGER = {"D1": 1.0, "D2": CARPAN, "D1n": 1.0, "D1p": 1.0}
GATE = {"E_A0Cg": 3, "E_bir3": 3, "D1": 3, "D2": 3, "D1n": 1, "D1p": 2, "A1g": 3}
ONTEST = ["E_A0g", "E_A0Cg", "R_A0", "R_A3S", "D1", "D2", "E_bir0", "E_bir3", "R_A0C", "R_A1", "A1g", "D1n", "D1p"]
GELISTIRME = ["D1", "D2"]


def hucre(v: str, seq: str, seed: int = 0, dr: float = DR, t0: float = 1.0) -> Cell:
    kol, argv = VARYANT[v]
    return Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence=seq, seed=int(seed), D_over_r0=float(dr),
                tau0_over_dt=float(t0), arms=(kol,), windows=windows_for(seq, FAZ2_WINDOWS_FILE), window_len=400,
                nees=True, cov_every=FAZ2_COV_EVERY, psf_legacy_scale=False, sensor_noise_dn=FAZ2_NOISE_DN,
                noise_seed=None, backend_argv=tuple(argv), tag=f"kapi_{v}", cache_dump=False)


def dokum_yolu(D: str, seq: str) -> str:
    return os.path.join(D, f"{seq}__dr{DR}__s0.pkl")


def _sha(p: str):
    if not os.path.exists(p):
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _isci(args):
    """Tek varyant, tek sekans: saklı döküm (sembolik bağ) üstünde arka-uç; yörünge dosyalarının sha256'sı kayda eklenir."""
    v, seq, yol = args
    os.chdir(ROOT)
    try:
        import cv2 as _cv2
        _cv2.setNumThreads(1)
    except Exception:
        pass
    import pickle
    import config
    cell = hucre(v, seq)
    kol = cell.arms[0]
    d = os.path.join(SPOOL, f"{seq}__{v}")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d, exist_ok=True)
    t0 = time.time()
    try:
        cfg = _build_cfg(cell, d)
        dp = config.dump_path(cfg)
        os.makedirs(os.path.dirname(dp), exist_ok=True)
        os.symlink(os.path.abspath(yol), dp)          # saklı döküme dokunulmaz (yalnız okunur)
        with open(dp, "rb") as f:
            f_ts = [int(fr["timestamp"]) for fr in pickle.load(f)["frames"]]
        base = {"t_frontend": None, "frontend_stats": {}, "fe_gpu": None, "quality": None}
        recs = _run_arms(cell, cfg, d, config.gt_path(cfg), f_ts, base)
        ek = {"varyant": v, "sekans": seq, "kol": kol, "argv": list(cell.backend_argv), "sure_s": round(time.time() - t0, 1),
              "dokum": os.path.basename(yol), "dokum_bayt": os.path.getsize(yol), "n_kare_dokum": len(f_ts),
              "sha_traj": _sha(os.path.join(d, f"traj_{kol}.csv")),
              "sha_rt": _sha(os.path.join(d, f"traj_{kol}_realtime.csv"))}
    except BaseException as e:                       # noqa: BLE001
        recs = [{"cell": asdict(cell), "key": cell.key(), "code_fp": code_fingerprint(), "status": "cell_fail",
                 "error": f"{type(e).__name__}: {e}"}]
        ek = {"varyant": v, "sekans": seq, "kol": kol, "argv": list(cell.backend_argv), "sure_s": round(time.time() - t0, 1)}
    finally:
        shutil.rmtree(d, ignore_errors=True)
    for r in recs:
        r["kapi"] = ek
    return recs


def oku(yollar=None) -> dict:
    """-> {(varyant, sekans): [kayıt, ...]} (bozuk satırlar atlanır; aynı pencere yinelenirse son kayıt geçerli)."""
    out = {}
    for p in (yollar or list(CIKTI.values())):
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                k = r.get("kapi") or {}
                if "varyant" in k:
                    out.setdefault((k["varyant"], k["sekans"]), {})[r.get("window")] = r
    return {k: list(v.values()) for k, v in out.items()}


def bitti(recs, n_pencere=3) -> bool:
    return len([r for r in recs if "arm" in r]) >= n_pencere


def _satirlar(p):
    with open(p) as f:
        for line in f:
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def basvuru(seqs=None) -> dict:
    """Eski kodun (ANA_FP) tohum 0, R0, D/r0 8 kayıtları: {(ad, sekans, pencere): kayıt}.
    ad: A0, A0C, A1, A3, AG, A0C+A3 (ana ızgara, k = 12); A3S (A3S ızgarası, k = 1); A0g (ilk sonda, geniş kapı)."""
    out = {}

    def uygun(c):
        return (c.get("seed") == 0 and c.get("D_over_r0") == float(DR) and abs(c.get("tau0_over_dt", 0) - 1.0) < 1e-12
                and (seqs is None or c.get("sequence") in seqs))

    for r in _satirlar(os.path.join(ANA, "FAZ2", "records.jsonl")):
        c = r.get("cell") or {}
        if "arm" in r and uygun(c) and c.get("tag") == "faz2" and r.get("code_fp") == ANA_FP:
            out[(r["arm"], c["sequence"], int(r["window"]))] = r
    for r in _satirlar(os.path.join(ANA, "A3S", "records.jsonl")):
        c = r.get("cell") or {}
        if "arm" in r and uygun(c) and c.get("tag") == "a3s" and r["arm"] == "A0C+A3" and r.get("code_fp") == ANA_FP:
            out[("A3S", c["sequence"], int(r["window"]))] = r
    for r in _satirlar(os.path.join(ANA, "SONDA_KESIF", "records.jsonl")):
        c = r.get("cell") or {}
        if "arm" in r and uygun(c) and c.get("tag") == "sonda_R0_fe0_bg" and r.get("code_fp") == ANA_FP:
            out[("A0g", c["sequence"], int(r["window"]))] = r
    return out


def _ate(r):
    v = r.get("m_ate_se3")
    return float(v) if (r.get("status") == "ok" and r.get("m_valid") and v is not None and math.isfinite(v) and v > 0) else None


def kapilar(asama: str) -> bool:
    """Ön beyan §2: T1 (gerileme), T2 (eşdeğerlik), T3 sert denetimleri. True: sonraki aşama koşulabilir.
    T3'ün yumuşak beklentileri (S1-S4) yalnız raporlanır."""
    C, ref, ok = oku(), basvuru([ONTEST_SEQ]), True
    s = ONTEST_SEQ

    def pen(v):
        return {int(r["window"]): r for r in C.get((v, s), []) if "arm" in r}

    def tam(v):
        p = pen(v)
        if len(p) < 3:
            print(f"[kapı] EKSİK {v}: {len(p)} pencere; hata {[r.get('error') for r in C.get((v, s), []) if 'arm' not in r]}",
                  flush=True)
        return len(p) >= 3

    # ---- T1 ----
    t1 = (("R_A3", "A3"),) if asama == "a3" else (("R_A0", "A0"), ("R_A0C", "A0C"), ("R_A1", "A1"), ("R_A3S", "A3S"))
    for v, ad in t1:
        if not tam(v):
            ok = False
            continue
        n = 0
        for w, r in sorted(pen(v).items()):
            m = ref.get((ad, s, w)) or {}
            ayni = bool(m) and all(r.get(k) == m.get(k) for k in BIT) and r.get("status") == "ok"
            n += ayni
            if not ayni:
                print(f"[kapı] T1 FARK {v} w{w}: {[(k, r.get(k), m.get(k)) for k in BIT]}", flush=True)
                ok = False
        print(f"[kapı] T1 {v:6s} (kapı kapalı) eski kodun {ad} kaydıyla bit-aynı: {n}/3", flush=True)
    if asama == "a3":
        print(f"[kapı] T1 (A3): {'GEÇTİ' if ok else 'KALDI'}", flush=True)
        return ok

    # ---- T2 ----
    for a, b, acik in (("E_A0Cg", "E_A0g", f"A0C(k = {S6}) + gate 3  ==  A0 + üç bayrak (s = {S6})"),
                       ("E_bir3", "E_bir0", "A0C(k = 1) + gate 3  ==  A0C(k = 1)")):
        if not (tam(a) and tam(b)):
            ok = False
            continue
        pa, pb = pen(a), pen(b)
        n = sum(all(pa[w].get(k) == pb[w].get(k) for k in BIT) and pa[w].get("status") == "ok" for w in pa)
        ka, kb = pa[min(pa)]["kapi"], pb[min(pb)]["kapi"]
        yor = ka.get("sha_traj") == kb.get("sha_traj") and ka.get("sha_traj") is not None
        rt = ka.get("sha_rt") == kb.get("sha_rt") and ka.get("sha_rt") is not None
        iyi = n == 3 and yor and rt
        ok &= iyi
        print(f"[kapı] T2 {acik}: metrik bit-aynı {n}/3; yörünge sha256 {'aynı' if yor else 'FARKLI'}; gerçek-zaman "
              f"yörüngesi sha256 {'aynı' if rt else 'FARKLI'} -> {'GEÇTİ' if iyi else 'KALDI'}", flush=True)
        if not iyi:
            for w in sorted(pa):
                print(f"         w{w}: {[(k, pa[w].get(k), pb[w].get(k)) for k in BIT]}", flush=True)
    # bilgi (kapı değil): A0C(k = 1) ile A0; yeni A0g ile ilk sondanın A0g kaydı (pnp 1e-6 yuvarlama farkı)
    for a, ad in (("E_bir0", "A0"), ("E_A0g", "A0g")):
        pa = pen(a)
        n = sum(all(pa[w].get(k) == (ref.get((ad, s, w)) or {}).get(k) for k in BIT) for w in pa)
        print(f"[kapı] bilgi: {a} ile eski {ad} kaydı bit-aynı {n}/{len(pa)}", flush=True)

    # ---- T3 sert ----
    for v in ("D1", "D2", "D1n", "D1p", "A1g"):
        if not tam(v):
            ok = False
            continue
        p = pen(v)
        for w, r in sorted(p.items()):
            dolu = (r.get("status") == "ok" and r.get("m_valid") and r.get("n_valid")
                    and abs(int(r.get("m_n_compare") or 0) - 400) <= 20
                    and all(r.get(k) is not None and math.isfinite(r.get(k)) for k in BIT))
            if not dolu:
                print(f"[kapı] T3-H1 EKSİK ALAN {v} w{w}: status={r.get('status')} m_valid={r.get('m_valid')} "
                      f"n_valid={r.get('n_valid')} n_compare={r.get('m_n_compare')} err={r.get('error')}", flush=True)
                ok = False
        r0 = p[min(p)]
        st = r0.get("arm_stats") or {}
        n_bek = int(max(r0["cell"]["windows"]) + r0["cell"]["window_len"]) - 1
        h2 = st.get("gate") == GATE[v] and st.get("gate_sigma_n") == n_bek
        if v == "A1g":
            h3 = abs(st.get("gate_sigma_min", 0) - S_TAM) < 1e-9 and abs(st.get("gate_sigma_max", 0) - S_TAM) < 1e-9
        else:
            taban = math.sqrt(1.0 + (C_DEGER[v] * SIGMA_TABAN) ** 2)
            h3 = abs(st.get("gate_sigma_min", 0) - taban) < 1e-9 and st.get("gate_sigma_max", 0) > 1.3
        ad = "A1" if v == "A1g" else "A3S"
        h4 = all(r.get("m_ate_se3") != (ref.get((ad, s, w)) or {}).get("m_ate_se3") for w, r in p.items())
        ok &= bool(h2 and h3 and h4)
        print(f"[kapı] T3 {v:3s}: H2 gate={st.get('gate')} n={st.get('gate_sigma_n')} (beklenen {GATE[v]}, {n_bek}) "
              f"{'ok' if h2 else 'KALDI'}; H3 sigma min {st.get('gate_sigma_min', float('nan')):.6f} / medyan "
              f"{st.get('gate_sigma_median', float('nan')):.4f} / maks {st.get('gate_sigma_max', float('nan')):.4f} "
              f"{'ok' if h3 else 'KALDI'}; H4 etki uygulandı (kapısız {ad} kaydından farklı) {'ok' if h4 else 'KALDI'}; "
              f"sigma_hat medyan {st.get('sigma_hat_median', float('nan')):.4f}", flush=True)

    # ---- T3 yumuşak (rapor) ----
    def ort(v_or_ad, kaynak):
        x = [_ate(r) for r in kaynak]
        return float(sum(x) / len(x)) if x and all(t is not None for t in x) else float("nan")

    a0 = ort("A0", [ref[("A0", s, w)] for w in sorted(pen("R_A0") or {100: 0}) if ("A0", s, w) in ref])
    for v, alt, ust in (("D1", 0.3, 1.0), ("D2", 0.3, 1.0), ("A1g", 0.3, 0.8), ("D1n", 0.3, 1.0), ("D1p", 0.7, 1.3)):
        if not pen(v):
            continue
        o = ort(v, list(pen(v).values())) / a0
        print(f"[kapı] T3-S3 (yumuşak) {v:3s} / A0 ATE (3 pencere ortalaması, {s}): {o:.3f}; beklenen [{alt}, {ust}] "
              f"{'tuttu' if alt <= o <= ust else 'tutmadı'}", flush=True)
    for v, alt, ust in (("D1", 1.6, 4.0), ("D2", 2.8, 9.0), ("D1n", 1.5, 4.0), ("D1p", 1.5, 4.0)):
        if not pen(v):
            continue
        m = (pen(v)[min(pen(v))].get("arm_stats") or {}).get("gate_sigma_median", float("nan"))
        print(f"[kapı] T3-S1/S2 (yumuşak) {v:3s} kapı sigması medyanı {m:.4f} px; beklenen [{alt}, {ust}] "
              f"{'tuttu' if alt <= m <= ust else 'tutmadı'}", flush=True)
    print(f"[kapı] ön test (T1 A3 hariç, T2, T3 sert): {'GEÇTİ' if ok else 'KALDI'}", flush=True)
    return ok


def bellek_mb() -> float:
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024.0
    return float("inf")


def kos(isler, jobs, bellek_alt, etiket, cikti):
    """isler: [(v, seq, yol)]. Bittikçe yazar. Bellek eşiğin altına inerse havuzu bırakır (çıkış 3)."""
    if not isler:
        return
    t_bas, sayac, hatalar = time.time(), {"n": 0}, []
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=jobs, maxtasksperchild=1) as cp, open(cikti, "a") as fout:
        def yaz(recs):
            for r in recs:
                fout.write(json.dumps(r, default=_js) + "\n")
            fout.flush()
            sayac["n"] += 1
            k = (recs[0].get("kapi") or {}) if recs else {}
            print(f"[kapi] {etiket} {sayac['n']}/{len(isler)} {time.strftime('%T')} {(time.time() - t_bas) / 60:.1f} dk: "
                  f"{k.get('varyant')} {k.get('sekans')} {[(r.get('window'), r.get('status')) for r in recs]} "
                  f"{k.get('sure_s')} s", flush=True)

        def hata(e):
            hatalar.append(repr(e))
            print(f"[kapi] İŞÇİ HATASI: {e!r}", flush=True)

        bek = [cp.apply_async(_isci, (i,), callback=yaz, error_callback=hata) for i in isler]
        while not all(b.ready() for b in bek):
            time.sleep(5.0)
            m = bellek_mb()
            if m < bellek_alt:
                print(f"[kapi] BELLEK {m:.0f} MB < {bellek_alt} MB: koşum bırakılıyor (yeniden başlatılabilir)", flush=True)
                cp.terminate()
                fout.flush()
                sys.exit(3)
    if hatalar:
        sys.exit(f"[kapi] DUR: {len(hatalar)} işçi hatası")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kuru", action="store_true")
    ap.add_argument("--asama", choices=("ontest", "a3", "gelistirme"), default="ontest")
    ap.add_argument("--dokumler", default=None)
    ap.add_argument("--jobs", type=int, default=5)
    ap.add_argument("--bellek_alt", type=float, default=3600.0, help="[MB] MemAvailable bunun altına inerse bırak")
    ap.add_argument("--yalniz_kapi", action="store_true", help="koşmadan kapıları değerlendir")
    a = ap.parse_args()

    os.makedirs(EV, exist_ok=True)
    fp = code_fingerprint()
    print(f"[kapi] code_fp {fp} (ana ağaç {ANA_FP}); kök {ROOT}; aşama {a.asama}; çıktı {CIKTI[a.asama]}", flush=True)
    if fp == ANA_FP:
        sys.exit("[kapi] DUR: bu ağaçta kapı kodu yok (parmak izi ana ağaçla aynı)")
    print(f"[kapi] oracle sigma_tilt {sigma_tilt(DR):.6f} px; s (A1 bireysel sigma) {S_TAM!r}; s6 {S6}; "
          f"pnp {PNP0 * float(S6)!r}; c(D2) {CARPAN!r}", flush=True)

    if a.asama == "ontest":
        P = [(v, ONTEST_SEQ) for v in ONTEST]
    elif a.asama == "a3":
        P = [("R_A3", ONTEST_SEQ)]
    else:
        P = [(v, s) for s in SEQS for v in GELISTIRME]
    for v, s in P:
        c = hucre(v, s)
        print(f"    {v:6s} {s:13s} kol {c.arms[0]:7s} {' '.join(c.backend_argv)}  [{c.tag}; pencereler {tuple(c.windows)}]")
    if a.kuru:
        import config
        c = hucre("R_A0", ONTEST_SEQ)
        b = config.backend_args(_build_cfg(c, "/tmp/_kapi_kuru"), argv=[])
        assert (b.pnp_reproj, b.px_noise, b.stereo_sigma) == (PNP0, PX0, ST0), (b.pnp_reproj, b.px_noise, b.stereo_sigma)
        assert int(getattr(b, "turb_gate", -1)) == 0
        ref = basvuru()
        n = {ad: sum(1 for k in ref if k[0] == ad) for ad in ("A0", "A0C", "A1", "A3", "A3S", "A0g")}
        print(f"[kapi] taban pnp/px/stereo {b.pnp_reproj}/{b.px_noise}/{b.stereo_sigma}; turb_gate varsayılanı 0; "
              f"başvuru kayıt sayıları {n} (beklenen her biri 27)")
        assert all(x == 27 for x in n.values()), n
        return
    if a.yalniz_kapi:
        sys.exit(0 if kapilar(a.asama if a.asama != "gelistirme" else "ontest") else 2)
    if not a.dokumler:
        sys.exit("[kapi] --dokumler gerekli")
    for s in sorted({s for _, s in P}):
        assert os.path.getsize(dokum_yolu(a.dokumler, s)) > 50e6, dokum_yolu(a.dokumler, s)
    os.makedirs(SPOOL, exist_ok=True)
    if a.asama == "gelistirme" and not kapilar("ontest"):
        sys.exit(2)

    C = oku()
    isler = [(v, s, dokum_yolu(a.dokumler, s)) for v, s in P if not bitti(C.get((v, s), []))]
    kos(isler, a.jobs, a.bellek_alt, a.asama, CIKTI[a.asama])
    if a.asama in ("ontest", "a3") and not kapilar(a.asama):
        sys.exit(2)
    C = oku()
    eksik = [(v, s) for v, s in P if not bitti(C.get((v, s), []))]
    print(f"[kapi] {a.asama} bitti {time.strftime('%F %T')}: {len(P) - len(eksik)}/{len(P)} koşum tam; eksik {eksik}", flush=True)


if __name__ == "__main__":
    main()
