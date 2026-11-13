#!/usr/bin/env python
"""Uyarlanır izleme kapısı kolunun doğrulayıcı ızgarası (ON_KAYIT_FAZ2.md §5 Sapma 6).

Tek, dondurulmuş yapılandırma koşulur: evidence/KAPI/SECIM.json (geliştirme koşumunun ön beyanlı seçim kuralı,
evidence/KAPI/ON_BEYAN_3.md §4). Tohumlar 1-4 (geliştirme ve sondalar yalnız tohum 0'ı gördü). Bu dosya kod parmak
izi kümesinin (turb.runner.code_fingerprint) dışındadır.

Kollar (kol, arka-uç bayrakları):
  A0    taban (doğrulama: eski kodun kaydıyla bit-aynı olmalı)
  A3S   çevrimiçi skaler kovaryans, kapı kapalı   (A0C+A3, --turb_cal_factor 1.0)
  KAPI  A3S + uyarlanır izleme kapısı             (seçilen aday: D1 ya da D2; --turb_gate 3)
  A1g   oracle skaler kovaryans + kapı            (A1, --turb_gate 3; üst sınır)
Parçalar (her biri kendi çıktı dosyasına yazar; ayrı süreçler olarak aynı anda koşabilir):
  I    R5 (v = 5 m/s), D/r0 {8, 4}: hızlı rejim ızgarasının (Sapma 5) saklı ön-uç dökümleri üstünde yalnız arka-uç.
       KAPI, A3S, A1g; tohum 1'de ayrıca A0 (hızlı ızgaranın A0 kaydıyla bit-aynı olmalı). A0 ve A1 tabanları hızlı
       ızgaranın kayıtlarından eşlenir (aynı döküm).
  II   R0 (tau0 = dt), D/r0 {8, 16}: A0 (ana ızgaranın A0 kaydıyla bit-aynı olmalı), KAPI, A1g. Döküm: A3S
       ızgarasından saklanmışsa o, yoksa ön-uç yeniden üretilir (belirlenimci; A0 denetimi her hücrede).
       A1 ana ızgaradan, A3S (kapı kapalı) A3S ızgarasından (Sapma 4) eşlenir.
  III  temiz (D/r0 0): A0 (bit-aynı denetimi), A3S, KAPI. Ön-uç yeniden üretilir.

    python kapi_dogrula.py --kuru [--aday D1]
    python kapi_dogrula.py --ontest --dokumler DIZIN[,DIZIN]            # tohum 0, V1_01_easy (ön test; kapılar)
    python kapi_dogrula.py --parca II --dokumler DIZIN[,DIZIN]          # ızgara parçası (devam ettirilebilir)
    python kapi_dogrula.py --parca I --dokumler DIZIN --bekle           # dökümler geldikçe koşar

Koşum sırasında her biten koşum denetlenir (tutmazsa yeni iş gönderilmez, çıkış kodu 2):
  (a) üç pencere kaydı var, "ok" kayıtlarda metrik alanları dolu;
  (b) A0: başvuru kaydıyla (ana ızgara / hızlı ızgara / sonda) bit-aynı; yeniden üretilen ön-ucun istatistikleri
      başvurununkiyle aynı (süre alanı hariç);
  (c) kapılı kollar: arm_stats.gate = 3, kapı sigması örnek sayısı = kare sayısı - 1; A1g'de kapı sigması = oracle
      bireysel sigma (1e-9); KAPI'de en küçük değer = taban (ilk DIMM kestiriminden önce).
Çıkış kodu: 0 tamam; 2 denetim tutmadı / hücre hatası; 3 bellek eşiği; 4 saklı döküm bekleniyor (--bekle yok).
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import multiprocessing as mp
import os
import pickle
import shutil
import sys
import threading
import time
from dataclasses import asdict

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from turb.runner import Cell, windows_for, _build_cfg, _run_arms, _frontend_worker, _js, code_fingerprint  # noqa: E402
from turb.experiments import (EUROC_FAZ2, FAZ2_WINDOWS_FILE, FAZ2_NOISE_DN, FAZ2_COV_EVERY, FAZ2_SEQS)  # noqa: E402
import hizli_kos as hk  # noqa: E402  (tau0/dt tanımı ve hücre anahtarı hızlı ızgarayla aynı olsun diye)
import kapi_kol as kk   # noqa: E402  (kol bayrakları geliştirme koşumuyla aynı olsun diye)

ANA = kk.ANA                                   # ana ağacın evidence dizini (başvuru kayıtları; salt okunur)
ANA_FP = kk.ANA_FP
EV = kk.EV
SECIM = os.path.join(EV, "SECIM.json")
SPOOL = "/tmp/turb_spool_kapi_dogrula"
SPOOL_FE = "/tmp/turb_spool_kapi_dogrula_fe"
TOHUMLAR = (1, 2, 3, 4)
SEQS = list(kk.SEQS)
assert sorted(SEQS) == sorted(FAZ2_SEQS)
ONTEST_SEQ = "V1_01_easy"
BIT = kk.BIT
KAPILI = ("KAPI", "A1g")
PARCA = {
    "I": dict(rejim="R5", duzeyler=(8, 4), kollar=("KAPI", "A3S", "A1g"), a0_tohum=(1,), kaynak="sakli"),
    "II": dict(rejim="R0", duzeyler=(8, 16), kollar=("A0", "KAPI", "A1g"), a0_tohum=(), kaynak="ikisi"),
    "III": dict(rejim="R0", duzeyler=(0,), kollar=("A0", "A3S", "KAPI"), a0_tohum=(), kaynak="fe"),
}
DOKUM_YAS_S = 180.0        # saklı döküm en az bu kadar eskiyse okunur (yazımı sürerken okunmasın)


def tau0(rejim: str, dr: float) -> float:
    return 1.0 if rejim == "R0" else float(hk.tau0_over_dt(rejim, dr))


def s_oracle(dr: float) -> float:
    """A1'in bireysel-faktör sigması: sqrt(sigma0^2 + (hypot(1, 1.75) sigma_tilt)^2)."""
    return math.sqrt(1.0 + (kk.CARPAN * kk.sigma_tilt(dr)) ** 2)


def varyant(v: str, aday: str):
    if v == "KAPI":
        return kk.VARYANT[aday]
    return {"A0": kk.VARYANT["R_A0"], "A3S": kk.VARYANT["R_A3S"], "A1g": kk.VARYANT["A1g"]}[v]


def hucre(v: str, seq: str, seed: int, dr: float, rejim: str, aday: str, tag=None) -> Cell:
    kol, argv = varyant(v, aday)
    return Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence=seq, seed=int(seed), D_over_r0=float(dr),
                tau0_over_dt=tau0(rejim, dr), arms=(kol,), windows=windows_for(seq, FAZ2_WINDOWS_FILE), window_len=400,
                nees=True, cov_every=FAZ2_COV_EVERY, psf_legacy_scale=False, sensor_noise_dn=FAZ2_NOISE_DN,
                noise_seed=None, backend_argv=tuple(argv), tag=tag or f"kapid_{v}", cache_dump=False)


def sakli_adlar(rejim: str, seq: str, seed: int, dr: float):
    """Saklayıcıların verdiği dosya adları (spool dizin adı = hücre anahtarı) + tohum 0 saklayıcısının kısa adı."""
    if rejim == "R5":
        return [hk._hucre(seq, seed, rejim, dr, hk.KOLLAR).key().replace("/", "_") + ".pkl"]
    a3s = hucre("A3S", seq, seed, dr, "R0", "D1", tag="a3s").key().replace("/", "_") + ".pkl"
    return [a3s, f"{seq}__dr{dr:g}__s{seed}.pkl"]


def sakli_dokum(rejim, seq, seed, dr, dizinler, yas=DOKUM_YAS_S):
    for d in dizinler:
        for ad in sakli_adlar(rejim, seq, seed, dr):
            p = os.path.join(d, ad)
            try:
                st = os.stat(p)
            except OSError:
                continue
            if st.st_size > 50e6 and (time.time() - st.st_mtime) >= yas:
                return p
    return None


def anahtar(e, v=None):
    k = (e["parca"], e["seq"], int(e["seed"]), float(e["dr"]))
    return k if v is None else k + (v,)


def plan(parcalar, tohumlar=TOHUMLAR, seqs=None):
    """Sıra: parça -> düzey -> tohum -> sekans; hücre içinde A0 en başta (doğrulama erken biter)."""
    out = []
    for p in parcalar:
        P = PARCA[p]
        for dr in P["duzeyler"]:
            for sd in tohumlar:
                for s in (seqs or SEQS):
                    vs = list(P["kollar"])
                    if sd in P["a0_tohum"] and "A0" not in vs:
                        vs = ["A0"] + vs
                    out.append(dict(parca=p, rejim=P["rejim"], seq=s, seed=int(sd), dr=float(dr), varyantlar=vs,
                                    kaynak=P["kaynak"]))
    return out


def ontest_plani():
    """Tohum 0, V1_01_easy: üç parça ve iki döküm yolu (yeniden üretim, saklı döküm) baştan sona."""
    s = ONTEST_SEQ
    return [dict(parca="II", rejim="R0", seq=s, seed=0, dr=8.0, varyantlar=["A0", "KAPI", "A1g"], kaynak="fe"),
            dict(parca="III", rejim="R0", seq=s, seed=0, dr=0.0, varyantlar=["A0", "A3S", "KAPI"], kaynak="fe"),
            dict(parca="I", rejim="R5", seq=s, seed=0, dr=8.0, varyantlar=["A0", "KAPI", "A3S", "A1g"], kaynak="fe"),
            dict(parca="II", rejim="R0", seq=s, seed=0, dr=16.0, varyantlar=["A0", "KAPI", "A1g"], kaynak="sakli")]


def _sha(p):
    if not os.path.exists(p):
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _isci(job):
    """Tek kol, tek hücre: döküm (sembolik bağ) üstünde arka-uç. Yörünge dosyalarının sha256'sı kayda eklenir."""
    os.chdir(ROOT)
    try:
        import cv2 as _cv2
        _cv2.setNumThreads(1)
    except Exception:
        pass
    import config
    v, seq, seed, dr, rejim, aday = job["v"], job["seq"], job["seed"], job["dr"], job["rejim"], job["aday"]
    cell = hucre(v, seq, seed, dr, rejim, aday)
    kol = cell.arms[0]
    d = os.path.join(SPOOL, f"{job['etiket']}__{job['parca']}__{seq}__dr{dr:g}__s{seed}__{v}")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d, exist_ok=True)
    t0 = time.time()
    ek = {"parca": job["parca"], "varyant": v, "sekans": seq, "tohum": int(seed), "dr": float(dr), "rejim": rejim,
          "kol": kol, "argv": list(cell.backend_argv), "secim": aday, "kaynak": job["kaynak"],
          "dokum": os.path.basename(job["dump"])}
    try:
        cfg = _build_cfg(cell, d)
        dp = config.dump_path(cfg)
        os.makedirs(os.path.dirname(dp), exist_ok=True)
        os.symlink(os.path.abspath(job["dump"]), dp)          # döküme dokunulmaz (yalnız okunur)
        with open(dp, "rb") as f:
            f_ts = [int(fr["timestamp"]) for fr in pickle.load(f)["frames"]]
        recs = _run_arms(cell, cfg, d, config.gt_path(cfg), f_ts, dict(job["base"]))
        ek.update(sure_s=round(time.time() - t0, 1), dokum_bayt=os.path.getsize(job["dump"]), n_kare_dokum=len(f_ts),
                  sha_traj=_sha(os.path.join(d, f"traj_{kol}.csv")),
                  sha_rt=_sha(os.path.join(d, f"traj_{kol}_realtime.csv")))
    except BaseException as e:                       # noqa: BLE001
        recs = [{"cell": asdict(cell), "key": cell.key(), "code_fp": code_fingerprint(), "status": "cell_fail",
                 "error": f"{type(e).__name__}: {e}"}]
        ek.update(sure_s=round(time.time() - t0, 1))
    finally:
        shutil.rmtree(d, ignore_errors=True)
    for r in recs:
        r["kapid"] = ek
    return recs


# ----------------------------------------------------------------------------- okuma / başvuru
def _satirlar(p):
    if not os.path.exists(p):
        return
    with open(p) as f:
        for line in f:
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def cikti_yolu(ontest: bool, parca: str | None = None) -> str:
    return os.path.join(EV, "dogrula_ontest.jsonl" if ontest else f"dogrula_{parca}.jsonl")


def oku(yollar) -> dict:
    """-> {(parça, sekans, tohum, dr, varyant): [kayıt, ...]} (aynı pencere yinelenirse son kayıt geçerli)."""
    out = {}
    for p in yollar:
        for r in _satirlar(p):
            k = r.get("kapid") or {}
            if "varyant" in k:
                out.setdefault((k["parca"], k["sekans"], int(k["tohum"]), float(k["dr"]), k["varyant"]),
                               {})[r.get("window")] = r
    return {k: list(v.values()) for k, v in out.items()}


def bitti(recs, n_pencere=3) -> bool:
    return len([r for r in recs if "arm" in r]) >= n_pencere


def basvuru_a0(ontest: bool) -> dict:
    """Eski kodun (ANA_FP) A0 kayıtları: {(rejim, sekans, tohum, dr): {"fe": frontend_stats, pencere: BIT alanları}}.
    R0: ana ızgara (etiket faz2). R5: hızlı ızgara (etiket hizli_R5); ön testte ayrıca sonda (sonda_R5_fe0, tohum 0)."""
    out = {}

    def ekle(rejim, r):
        c = r["cell"]
        d = out.setdefault((rejim, c["sequence"], int(c["seed"]), float(c["D_over_r0"])), {})
        d.setdefault("fe", r.get("frontend_stats"))
        d[int(r["window"])] = {x: r.get(x) for x in BIT + ("status",)}

    for r in _satirlar(os.path.join(ANA, "FAZ2", "records.jsonl")):
        c = r.get("cell") or {}
        if (r.get("arm") == "A0" and c.get("tag") == "faz2" and r.get("code_fp") == ANA_FP
                and abs(c.get("tau0_over_dt", 1.0) - 1.0) < 1e-12 and c.get("D_over_r0") in (0.0, 8.0, 16.0)):
            ekle("R0", r)
    for r in _satirlar(os.path.join(ANA, "HIZLI", "records.jsonl")):
        c = r.get("cell") or {}
        if (r.get("arm") == "A0" and c.get("tag") == "hizli_R5" and r.get("code_fp") == ANA_FP
                and abs(c.get("tau0_over_dt", 1.0) - tau0("R5", c.get("D_over_r0"))) < 1e-9):
            ekle("R5", r)
    if ontest:
        for r in _satirlar(os.path.join(ANA, "SONDA_KESIF", "records.jsonl")):
            c = r.get("cell") or {}
            if (r.get("arm") == "A0" and c.get("tag") == "sonda_R5_fe0" and r.get("code_fp") == ANA_FP
                    and c.get("seed") == 0):
                ekle("R5", r)
    return out


def secim_oku(aday_gecici=None):
    if os.path.exists(SECIM):
        with open(SECIM) as f:
            s = json.load(f)
        assert s["aday"] in ("D1", "D2") and int(s.get("gate", 3)) == 3, s
        assert abs(float(s["c"]) - kk.C_DEGER[s["aday"]]) < 1e-15, s
        return s["aday"], False
    if aday_gecici:
        return aday_gecici, True
    return None, True


# ----------------------------------------------------------------------------- denetim
def denetle(recs, e, v, aday, ref):
    """-> (özet, sorunlar, a0_bekliyor). Sorun = DUR gerektiren durum."""
    sorun = []
    kayit = [r for r in recs if "arm" in r]
    if len(kayit) < 3:
        return "hücre hatası: " + str((recs[0] if recs else {}).get("error")), ["hücre hatası"], False
    n_ok = sum(1 for r in kayit if r.get("status") == "ok")
    bos = sum(1 for r in kayit if r.get("status") == "ok" and any(r.get(k) is None for k in BIT))
    if bos:
        sorun.append(f"{bos} 'ok' kaydında boş metrik alanı")
    ek, bekliyor = "", False
    if v == "A0":
        m = ref.get((e["rejim"], e["seq"], int(e["seed"]), float(e["dr"])))
        if m is None or not all(int(r["window"]) in m for r in kayit):
            ek, bekliyor = ", A0 başvurusu henüz yok (ertelendi)", True
        else:
            n = sum(all(r.get(k) == m[int(r["window"])].get(k) for k in BIT) for r in kayit)
            ek = f", A0 başvuruyla bit-aynı {n}/{len(kayit)}"
            if n != len(kayit):
                sorun.append("A0 başvuru kaydıyla bit-aynı değil: " + str(
                    [(r["window"], r.get("m_ate_se3"), m[int(r["window"])].get("m_ate_se3")) for r in kayit]))
    if v in KAPILI:
        st = kayit[0].get("arm_stats") or {}
        n_bek = int(max(kayit[0]["cell"]["windows"]) + kayit[0]["cell"]["window_len"]) - 1
        mn, mx = st.get("gate_sigma_min", float("nan")), st.get("gate_sigma_max", float("nan"))
        if st.get("gate") != 3 or st.get("gate_sigma_n") != n_bek:
            sorun.append(f"kapı istatistiği: gate={st.get('gate')} n={st.get('gate_sigma_n')} (beklenen 3, {n_bek})")
        if v == "A1g":
            s = s_oracle(e["dr"])
            if not (abs(mn - s) < 1e-9 and abs(mx - s) < 1e-9):
                sorun.append(f"A1g kapı sigması {mn}..{mx}, beklenen {s}")
        else:
            taban = math.sqrt(1.0 + (kk.C_DEGER[aday] * kk.SIGMA_TABAN) ** 2)
            if not abs(mn - taban) < 1e-9:
                sorun.append(f"KAPI kapı sigması en küçük {mn}, beklenen taban {taban}")
        ek += f", kapı sigması medyan {st.get('gate_sigma_median', float('nan')):.3f} maks {mx:.3f}"
    return f"{n_ok}/{len(kayit)} ok{ek}", sorun, bekliyor


def a0_denetle(kayit_ozet, m):
    """Ertelenen A0 denetimi: kayit_ozet {pencere: BIT alanları}. -> None (başvuru yok) | True | False"""
    if m is None or not all(w in m for w in kayit_ozet):
        return None
    return all(kayit_ozet[w].get(k) == m[w].get(k) for w in kayit_ozet for k in BIT)


def bellek_mb() -> float:
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024.0
    return float("inf")


# ----------------------------------------------------------------------------- koşum
def kos(E, a, aday, cikti, etiket, ontest):
    """E: plan girdileri (kalan varyantlarıyla). Dökümü hazır olan hücrelerin kolları CPU havuzuna, ön-uç gereken
    hücreler GPU havuzuna gönderilir; bellek kabul eşiğinin altındayken yeni iş gönderilmez."""
    dizinler = [d for d in (a.dokumler or "").split(",") if d]
    ref = basvuru_a0(ontest)
    kilit = threading.Lock()
    dur = threading.Event()
    S = dict(be_calisan=0, fe_calisan=0, be_biten=0, sorun=0, fe_hata=0, t_ref=time.time())
    kuyruk = collections.deque()
    bekleyen_a0 = {}                       # anahtar -> {pencere: BIT alanları}
    harita = {anahtar(e): e for e in E}
    for e in E:
        e.update(durum="yok", dump=None, fe_dir=None, base=None, kalan=len(e["varyantlar"]))
    n_be = sum(e["kalan"] for e in E)
    t_bas = time.time()
    ctx = mp.get_context("spawn")
    fe_gerek = any(e["kaynak"] in ("fe", "ikisi") for e in E)
    gp = ctx.Pool(processes=max(1, a.gpu_jobs), maxtasksperchild=8) if fe_gerek else None
    cp = ctx.Pool(processes=a.jobs, maxtasksperchild=1)
    fout = open(cikti, "a")

    def kuyruga(e):
        for v in e["varyantlar"]:
            kuyruk.append((e, v))

    def fe_bitti(meta):
        with kilit:
            S["fe_calisan"] -= 1
            c = meta["cell"]
            rejim = "R0" if abs(c.tau0_over_dt - 1.0) < 1e-12 else "R5"
            e = harita[(next(p for p in ("I", "II", "III") if (p, c.sequence, int(c.seed), float(c.D_over_r0)) in harita
                             and harita[(p, c.sequence, int(c.seed), float(c.D_over_r0))]["rejim"] == rejim),
                        c.sequence, int(c.seed), float(c.D_over_r0))]
            if not meta.get("ok"):
                S["fe_hata"] += 1
                e["durum"], e["kalan"] = "hata", 0
                fout.write(json.dumps({"cell": asdict(c), "key": c.key(), "code_fp": code_fingerprint(),
                                       "status": "cell_fail", "error": meta.get("error", "ön-uç başarısız"),
                                       "kapid_fe": {"parca": e["parca"], "sekans": e["seq"], "tohum": e["seed"],
                                                    "dr": e["dr"]}}, default=_js) + "\n")
                fout.flush()
                print(f"[dogrula] ÖN-UÇ BAŞARISIZ {anahtar(e)}: {meta.get('error')}", flush=True)
                dur.set()
                return
            qual = None
            qp = os.path.join(meta["dir"], "turb_quality.json")
            if os.path.exists(qp):
                try:
                    with open(qp) as f:
                        qual = json.load(f)
                except json.JSONDecodeError:
                    pass
            e.update(durum="hazir", dump=meta["dump"], fe_dir=meta["dir"],
                     base={"t_frontend": meta.get("t_frontend"), "frontend_stats": meta.get("stats"),
                           "fe_gpu": meta.get("gpu"), "quality": qual})
            m = ref.get((e["rejim"], e["seq"], int(e["seed"]), float(e["dr"])))
            fe_ok = None if m is None else hk.fe_ayni(meta.get("stats"), m.get("fe"))
            print(f"[dogrula] ön-uç tamam {anahtar(e)} ({meta.get('t_frontend', 0):.0f} s); istatistikler başvuruyla "
                  f"{'aynı' if fe_ok else ('başvuru yok' if fe_ok is None else 'FARKLI')}", flush=True)
            if fe_ok is False:
                S["sorun"] += 1
                print(f"[dogrula] DUR: {anahtar(e)}: ön-uç istatistikleri başvurudan farklı", flush=True)
                dur.set()
            kuyruga(e)

    def be_bitti(recs):
        with kilit:
            S["be_calisan"] -= 1
            S["be_biten"] += 1
            for r in recs:
                fout.write(json.dumps(r, default=_js) + "\n")
            fout.flush()
            k = recs[0]["kapid"]
            e = harita[(k["parca"], k["sekans"], int(k["tohum"]), float(k["dr"]))]
            v = k["varyant"]
            e["kalan"] -= 1
            ozet, sorun, bekliyor = denetle(recs, e, v, aday, ref)
            if bekliyor:
                bekleyen_a0[anahtar(e)] = {int(r["window"]): {x: r.get(x) for x in BIT} for r in recs if "arm" in r}
            el = (time.time() - t_bas) / 60.0
            print(f"[dogrula] {etiket} {S['be_biten']}/{n_be} {time.strftime('%T')} {el:.1f} dk: {anahtar(e)} {v} "
                  f"{ozet}; {k.get('sure_s')} s", flush=True)
            if sorun:
                S["sorun"] += 1
                for s in sorun[:4]:
                    print(f"[dogrula] DUR: {anahtar(e)} {v}: {s}", flush=True)
                dur.set()
            if e["kalan"] <= 0 and e["fe_dir"]:
                shutil.rmtree(e["fe_dir"], ignore_errors=True)

    def hata(x):
        with kilit:
            S["sorun"] += 1
            print(f"[dogrula] İŞ HATASI: {x!r}", flush=True)
            dur.set()

    def fe_hata(x):
        with kilit:
            S["fe_calisan"] -= 1
        hata(x)

    cikis = 0
    try:
        while True:
            with kilit:
                kalan = [e for e in E if e["kalan"] > 0]
                if not kalan:
                    break
                if dur.is_set() and S["be_calisan"] == 0 and S["fe_calisan"] == 0:
                    break
                mem = bellek_mb()
                if mem < a.bellek_alt:
                    print(f"[dogrula] BELLEK {mem:.0f} MB < {a.bellek_alt} MB: koşum bırakılıyor (yeniden başlatılabilir)",
                          flush=True)
                    cikis = 3
                    break
                if not dur.is_set():
                    # saklı dökümler
                    for e in kalan:
                        if e["durum"] == "yok" and e["kaynak"] in ("sakli", "ikisi"):
                            p = sakli_dokum(e["rejim"], e["seq"], e["seed"], e["dr"], dizinler)
                            if p:
                                e.update(durum="hazir", dump=p, kaynak="sakli",
                                         base={"t_frontend": None, "frontend_stats": {}, "fe_gpu": None, "quality": None})
                                kuyruga(e)
                    # ön-uç (plan sırasıyla; açık döküm sayısı sınırlı)
                    if gp is not None:
                        acik = sum(1 for e in kalan if e["durum"] in ("kosuyor", "hazir") and e["fe_dir"] is None
                                   and e["kaynak"] in ("fe", "ikisi")) + sum(1 for e in kalan if e["fe_dir"])
                        for e in kalan:
                            if acik >= a.gpu_jobs + 2 or S["fe_calisan"] >= a.gpu_jobs or mem < a.bellek_kabul_fe:
                                break
                            if e["durum"] == "yok" and e["kaynak"] in ("fe", "ikisi"):
                                e["durum"] = "kosuyor"
                                S["fe_calisan"] += 1
                                acik += 1
                                c = hucre("A0", e["seq"], e["seed"], e["dr"], e["rejim"], aday, tag="kapid_fe")
                                gp.apply_async(_frontend_worker, ((c, SPOOL_FE),), callback=fe_bitti,
                                               error_callback=fe_hata)
                    # arka-uç
                    while kuyruk and S["be_calisan"] < a.jobs and mem >= a.bellek_kabul:
                        e, v = kuyruk.popleft()
                        S["be_calisan"] += 1
                        job = dict(parca=e["parca"], v=v, seq=e["seq"], seed=e["seed"], dr=e["dr"], rejim=e["rejim"],
                                   aday=aday, dump=e["dump"], base=e["base"], etiket=etiket,
                                   kaynak="sakli" if e["fe_dir"] is None else "fe")
                        cp.apply_async(_isci, (job,), callback=be_bitti, error_callback=hata)
                        mem -= 1500.0                    # yeni işçinin kabaca yer kaplayacağı varsayılır
                    bos = (S["be_calisan"] == 0 and S["fe_calisan"] == 0 and not kuyruk
                           and all(e["durum"] == "yok" and e["kaynak"] == "sakli" for e in kalan))
                else:
                    bos = False
                # ertelenen A0 denetimleri (hızlı ızgaranın kayıtları geldikçe)
                if bekleyen_a0 and time.time() - S["t_ref"] > 120.0:
                    S["t_ref"] = time.time()
                    ref = basvuru_a0(ontest)
                    for k in list(bekleyen_a0):
                        e = harita[k]
                        sonuc = a0_denetle(bekleyen_a0[k], ref.get((e["rejim"], e["seq"], int(e["seed"]), float(e["dr"]))))
                        if sonuc is None:
                            continue
                        del bekleyen_a0[k]
                        print(f"[dogrula] ertelenen A0 denetimi {k}: {'bit-aynı' if sonuc else 'FARKLI'}", flush=True)
                        if not sonuc:
                            S["sorun"] += 1
                            print(f"[dogrula] DUR: {k} A0: başvuru kaydıyla bit-aynı değil", flush=True)
                            dur.set()
            if bos:
                if not a.bekle:
                    print(f"[dogrula] saklı dökümü henüz olmayan {len(kalan)} hücre var; --bekle verilmedi, çıkılıyor",
                          flush=True)
                    cikis = 4
                    break
                time.sleep(30.0)
            else:
                time.sleep(3.0)
    finally:
        if cikis == 3:
            cp.terminate()
            if gp is not None:
                gp.terminate()
        else:
            cp.close()
            cp.join()
            if gp is not None:
                gp.close()
                gp.join()
        fout.flush()
        fout.close()
    # son ertelenen denetimler
    if bekleyen_a0:
        ref = basvuru_a0(ontest)
        for k in list(bekleyen_a0):
            e = harita[k]
            sonuc = a0_denetle(bekleyen_a0[k], ref.get((e["rejim"], e["seq"], int(e["seed"]), float(e["dr"]))))
            print(f"[dogrula] ertelenen A0 denetimi {k}: "
                  f"{'başvuru hâlâ yok (analiz betiği denetler)' if sonuc is None else 'bit-aynı' if sonuc else 'FARKLI'}",
                  flush=True)
            if sonuc is False:
                S["sorun"] += 1
    print(f"[dogrula] {etiket} bitti {time.strftime('%F %T')}: {S['be_biten']}/{n_be} arka-uç koşumu, sorun {S['sorun']}, "
          f"ön-uç hatası {S['fe_hata']}, {(time.time() - t_bas) / 60:.1f} dk", flush=True)
    if S["sorun"] or S["fe_hata"] or dur.is_set():
        return 2
    return cikis


# ----------------------------------------------------------------------------- ön test kapıları
def ontest_kapilari(aday) -> bool:
    """Sapma 6 ön testi (tohum 0, V1_01_easy). K1-K4 sert; K5 yumuşak (yalnız rapor)."""
    C = oku([cikti_yolu(True)])
    ref = basvuru_a0(True)
    G = kk.oku()                                     # geliştirme / ilk ön test kayıtları (aynı döküm, aynı bayraklar)
    ok = True
    s = ONTEST_SEQ

    def pen(k):
        return {int(r["window"]): r for r in C.get(k, []) if "arm" in r}

    for e in ontest_plani():
        for v in e["varyantlar"]:
            k = anahtar(e, v)
            p = pen(k)
            if len(p) < 3:
                print(f"[kapı] EKSİK {k}: {len(p)} pencere; hata {[r.get('error') for r in C.get(k, []) if 'arm' not in r]}",
                      flush=True)
                ok = False
                continue
            # K4: alanlar dolu
            for w, r in sorted(p.items()):
                dolu = (r.get("status") == "ok" and r.get("m_valid") and r.get("n_valid")
                        and abs(int(r.get("m_n_compare") or 0) - 400) <= 20
                        and all(r.get(x) is not None and math.isfinite(r.get(x)) for x in BIT))
                if not dolu:
                    print(f"[kapı] K4 EKSİK ALAN {k} w{w}: status={r.get('status')} m_valid={r.get('m_valid')} "
                          f"n_valid={r.get('n_valid')} n_compare={r.get('m_n_compare')} err={r.get('error')}", flush=True)
                    ok = False
            ozet, sorun, bekliyor = denetle(list(p.values()), e, v, aday, ref)
            if sorun or bekliyor:
                ok = False
            ek = ""
            r0 = p[min(p)]
            # K2: yeniden üretilen ön-uç istatistikleri
            if v == "A0" and (r0.get("kapid") or {}).get("kaynak") == "fe":
                m = ref.get((e["rejim"], e["seq"], 0, float(e["dr"]))) or {}
                fe_ok = hk.fe_ayni(r0.get("frontend_stats"), m.get("fe"))
                ok &= bool(fe_ok)
                ek += f"; K2 ön-uç istatistikleri başvuruyla {'aynı' if fe_ok else 'FARKLI'}"
            # K3: geliştirme kayıtlarıyla eşdeğerlik (R0, D/r0 8)
            if e["rejim"] == "R0" and e["dr"] == 8.0 and v in KAPILI:
                gv = aday if v == "KAPI" else "A1g"
                g = {int(r["window"]): r for r in G.get((gv, s), []) if "arm" in r}
                n = sum(all(p[w].get(x) == (g.get(w) or {}).get(x) for x in BIT) for w in p)
                kg = (g[min(g)].get("kapi") or {}) if g else {}
                yor = kg.get("sha_traj") is not None and kg.get("sha_traj") == r0["kapid"].get("sha_traj")
                rt = kg.get("sha_rt") is not None and kg.get("sha_rt") == r0["kapid"].get("sha_rt")
                iyi = n == 3 and yor and rt
                ok &= iyi
                ek += (f"; K3 geliştirme kaydıyla ({gv}, saklı döküm) bit-aynı {n}/3, yörünge sha256 "
                       f"{'aynı' if yor else 'FARKLI'}, gerçek-zaman {'aynı' if rt else 'FARKLI'}")
            print(f"[kapı] {k}: {ozet}{ek}{' SORUN: ' + '; '.join(sorun) if sorun else ''}", flush=True)

    # K5 (yumuşak): ATE oranları (3 pencere ortalaması) ve süreler
    def ort(k):
        x = [kk._ate(r) for r in pen(k).values()]
        return float(sum(x) / len(x)) if x and all(t is not None for t in x) else float("nan")

    for e in ontest_plani():
        a0 = ort(anahtar(e, "A0"))
        for v in e["varyantlar"]:
            if v == "A0":
                continue
            k = anahtar(e, v)
            sure = (pen(k)[min(pen(k))]["kapid"].get("sure_s") if pen(k) else None)
            print(f"[kapı] K5 (yumuşak) {e['rejim']} D/r0 {e['dr']:g} {v:4s} / A0 ATE (3 pencere ortalaması): "
                  f"{ort(k) / a0:.3f}; arka-uç {sure} s", flush=True)
    print(f"[kapı] Sapma 6 ön testi (K1-K4): {'GEÇTİ' if ok else 'KALDI'}", flush=True)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kuru", action="store_true")
    ap.add_argument("--ontest", action="store_true")
    ap.add_argument("--parca", choices=tuple(PARCA), default=None)
    ap.add_argument("--aday", choices=("D1", "D2"), default=None, help="yalnız --kuru için geçici aday (SECIM.json yoksa)")
    ap.add_argument("--dokumler", default="", help="saklı döküm dizinleri (virgüllü)")
    ap.add_argument("--tohumlar", default=",".join(str(t) for t in TOHUMLAR))
    ap.add_argument("--yalniz", default="", help="yalnız bu sekanslar (virgüllü)")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--gpu_jobs", type=int, default=1)
    ap.add_argument("--bekle", action="store_true", help="saklı dökümler gelene kadar bekle (parça I)")
    ap.add_argument("--bellek_alt", type=float, default=3600.0, help="[MB] MemAvailable bunun altına inerse bırak")
    ap.add_argument("--bellek_kabul", type=float, default=5500.0, help="[MB] bunun altındayken yeni arka-uç gönderme")
    ap.add_argument("--bellek_kabul_fe", type=float, default=7000.0, help="[MB] bunun altındayken yeni ön-uç gönderme")
    ap.add_argument("--yalniz_kapi", action="store_true", help="koşmadan ön test kapılarını değerlendir")
    a = ap.parse_args()

    os.makedirs(EV, exist_ok=True)
    fp = code_fingerprint()
    aday, gecici = secim_oku(a.aday if a.kuru else None)
    print(f"[dogrula] code_fp {fp} (ana ağaç {ANA_FP}); kök {ROOT}; seçilen aday {aday}"
          f"{' (GEÇİCİ: SECIM.json yok)' if gecici else ' (SECIM.json sha256 ' + _sha(SECIM)[:16] + ')'}", flush=True)
    if aday is None:
        sys.exit("[dogrula] DUR: evidence/KAPI/SECIM.json yok (yapılandırma dondurulmadan koşulmaz)")
    if fp == ANA_FP:
        sys.exit("[dogrula] DUR: bu ağaçta kapı kodu yok (parmak izi ana ağaçla aynı)")
    if not (a.ontest or a.parca or a.kuru):
        sys.exit("[dogrula] --ontest ya da --parca gerekli")
    toh = tuple(int(x) for x in a.tohumlar.split(",") if x != "")
    sec = [s for s in a.yalniz.split(",") if s]
    assert all(s in SEQS for s in sec), "bilinmeyen sekans"
    assert a.ontest or all(t in TOHUMLAR for t in toh), "doğrulayıcı ızgara yalnız tohum 1-4"
    if a.ontest or (a.kuru and not a.parca):
        E, etiket, cikti = ontest_plani(), "ontest", cikti_yolu(True)
    else:
        E, etiket, cikti = plan([a.parca], toh, [s for s in SEQS if not sec or s in sec]), a.parca, cikti_yolu(False, a.parca)
    if a.kuru and not a.ontest and not a.parca:
        E = ontest_plani() + plan(list(PARCA))
    print(f"[dogrula] kollar: " + "; ".join(f"{v} = {varyant(v, aday)[0]} {' '.join(varyant(v, aday)[1])}"
                                            for v in ("A0", "A3S", "KAPI", "A1g")), flush=True)
    for dr in sorted({e["dr"] for e in E}):
        print(f"    D/r0 {dr:4g}: oracle bireysel sigma {s_oracle(dr) if dr > 0 else 1.0:.6f} px; tau0/dt R5 "
              f"{tau0('R5', dr) if dr > 0 else float('nan'):.5f}")
    n_fe = sum(1 for e in E if e["kaynak"] in ("fe", "ikisi"))
    print(f"[dogrula] {etiket}: {len(E)} hücre, {sum(len(e['varyantlar']) for e in E)} arka-uç koşumu "
          f"({3 * sum(len(e['varyantlar']) for e in E)} kayıt); ön-uç gerekebilecek hücre {n_fe}; çıktı {cikti}", flush=True)

    if a.kuru:
        import config
        c = hucre("A0", ONTEST_SEQ, 0, 8.0, "R0", aday)
        b = config.backend_args(_build_cfg(c, "/tmp/_kapid_kuru"), argv=[])
        assert (b.pnp_reproj, b.px_noise, b.stereo_sigma) == (kk.PNP0, kk.PX0, kk.ST0)
        assert int(getattr(b, "turb_gate", -1)) == 0
        # hücre tanımı ana ızgara / hızlı ızgara ile aynı mı (etiket, kollar ve bayraklar dışında)
        from turb.experiments import cells_FAZ2
        ana = {(x.sequence, x.seed, x.D_over_r0): x for x in cells_FAZ2()}
        for e in plan(["II", "III"]):
            x, y = asdict(ana[(e["seq"], e["seed"], e["dr"])]), asdict(hucre("A0", e["seq"], e["seed"], e["dr"], "R0", aday))
            fark = {k for k in x if x[k] != y[k]} - {"arms", "backend_argv", "tag"}
            assert not fark, (e, fark)
        for e in plan(["I"]):
            x = asdict(hk._hucre(e["seq"], e["seed"], "R5", e["dr"], hk.KOLLAR))
            y = asdict(hucre("A0", e["seq"], e["seed"], e["dr"], "R5", aday))
            fark = {k for k in x if x[k] != y[k]} - {"arms", "backend_argv", "tag"}
            assert not fark, (e, fark)
        ref = basvuru_a0(True)
        say = collections.Counter(k[0] + ("@%g" % k[3]) + ("/s0" if k[2] == 0 else "/s1-4") for k in ref)
        print(f"[dogrula] hücre tanımları ana ızgara (II, III) ve hızlı ızgara (I) ile aynı (kol, bayrak, etiket dışında); "
              f"A0 başvuru hücreleri {dict(sorted(say.items()))}")
        dz = [d for d in a.dokumler.split(",") if d]
        for e in E[:4] + E[-2:]:
            print(f"    {anahtar(e)} {e['varyantlar']} kaynak {e['kaynak']}; saklı ad(lar) "
                  f"{sakli_adlar(e['rejim'], e['seq'], e['seed'], e['dr'])}; bulunan "
                  f"{sakli_dokum(e['rejim'], e['seq'], e['seed'], e['dr'], dz)}")
        if dz:
            var = sum(1 for e in E if sakli_dokum(e["rejim"], e["seq"], e["seed"], e["dr"], dz))
            print(f"[dogrula] saklı dökümü şu an bulunan hücre: {var}/{len(E)}")
        print("[dogrula] koşum BAŞLATILMADI")
        return 0
    if a.yalniz_kapi:
        return 0 if ontest_kapilari(aday) else 2
    if gecici:
        sys.exit("[dogrula] DUR: SECIM.json yok")
    env = {k: os.environ.get(k) for k in ("CUDA_DEVICE_ORDER", "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS")}
    if n_fe and (env["CUDA_DEVICE_ORDER"] != "PCI_BUS_ID" or env["CUDA_VISIBLE_DEVICES"] != "1"):
        sys.exit("[dogrula] DUR: yalnız RTX 3090 (CUDA_DEVICE_ORDER=PCI_BUS_ID, CUDA_VISIBLE_DEVICES=1)")
    os.makedirs(SPOOL, exist_ok=True)
    os.makedirs(SPOOL_FE, exist_ok=True)

    C = oku([cikti])
    for e in E:
        e["varyantlar"] = [v for v in e["varyantlar"] if not bitti(C.get(anahtar(e, v), []))]
    E = [e for e in E if e["varyantlar"]]
    print(f"[dogrula] kalan: {len(E)} hücre, {sum(len(e['varyantlar']) for e in E)} arka-uç koşumu; ortam {env}; "
          f"jobs {a.jobs} gpu_jobs {a.gpu_jobs}", flush=True)
    kod = kos(E, a, aday, cikti, etiket, a.ontest) if E else 0
    if a.ontest and kod == 0 and not ontest_kapilari(aday):
        kod = 2
    return kod


if __name__ == "__main__":
    sys.exit(main())
