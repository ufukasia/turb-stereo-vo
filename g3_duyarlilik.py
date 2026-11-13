#!/usr/bin/env python
"""G3 incelemesi: dosya yolu / kanca yolu farkının kaynağı ve arka-ucun küçük kalibrasyon farkına duyarlılığı.

Tanı aracı (ön kayıt dışı, karar üretmez); beklenenler evidence/ONARIM/ON_TEST.md §2b'de sonuçtan önce yazılıdır.
Parmak izli kodda değişiklik yok; bu dosya parmak izi kümesinin dışındadır.

    python g3_duyarlilik.py --dizin G fe                 # T-F için dosya yolu ön-ucu (ekran kartı; bellek boşalınca)
    python g3_duyarlilik.py --dizin G fark               # T-F döküm karşılaştırması (yalnız işlemci)
    python g3_duyarlilik.py --dizin G be T-K|T-D1|T-D1k|T-D2|T-E      # yalnız arka-uç, A0

G: çalışma dizini; G/dumps/<SEKANS>__landmarks_<SEKANS>.pkl kanca dökümleri (A3S ızgarasından kopya),
G/dumps/DOSYA__V1_01_easy.pkl dosya yolu dökümü (fe adımı yazar). Çıktı: evidence/ONARIM/g3_duyarlilik.jsonl
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import pickle
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

import a3s_kos  # noqa: E402
import onarim_kos  # noqa: E402
from turb.runner import _build_cfg, _backend_worker, _frontend_worker, _js, code_fingerprint  # noqa: E402

CIKTI = os.path.join(ROOT, "evidence", "ONARIM", "g3_duyarlilik.jsonl")
ONTEST = os.path.join(ROOT, "evidence", "ONARIM", "on_test.jsonl")
ALAN = a3s_kos.BIT_ALANLARI
KALIB = ("rect_K", "baseline", "T_BS", "R1")
SEQ_D = "V1_01_easy"                      # dosya yolu ön testi bu sekansta yapıldı


def yaz(d):
    d = dict(d, zaman=time.strftime("%F %T"), code_fp=code_fingerprint())
    with open(CIKTI, "a") as f:
        f.write(json.dumps(d, ensure_ascii=False, default=_js) + "\n")


def kanca_dokumu(G, seq):
    with open(os.path.join(G, "dumps", f"{seq}__landmarks_{seq}.pkl"), "rb") as f:
        return pickle.load(f)


def dosya_dokumu(G):
    with open(os.path.join(G, "dumps", f"DOSYA__{SEQ_D}.pkl"), "rb") as f:
        return pickle.load(f)


def kalibrasyon(base_path, sequence):
    """Ön-ucun arka-uca verdiği kalibrasyon (işlemcide; StereoMatcher'ın yalnız kalibrasyon adımları)."""
    import config
    import stereo_front_end as sfe
    cfg = config.load(str(config.DEFAULT_PATH))
    m = object.__new__(sfe.StereoMatcher)
    m.base_path, m.sequence = Path(base_path), sequence
    m.STEREO_EXTRINSICS = cfg["frontend"]["stereo_extrinsics"]
    m.load_calibration()
    m.compute_rectification()
    return m.get_calibration_export()


def kalib_fark(a, b):
    from gtsam_backend import T_body_rectcam0
    ka, kb = np.asarray(a["rect_K"], float), np.asarray(b["rect_K"], float)
    tbc = lambda x: T_body_rectcam0(np.asarray(x["T_BS"], np.float64), np.asarray(x["R1"], np.float64))
    return dict(rect_K_a=ka.tolist(), rect_K_b=kb.tolist(), rect_K_goreli=((kb - ka) / ka).tolist(),
                taban_goreli=float(b["baseline"] / a["baseline"] - 1.0),
                T_B_C_enbuyuk=float(np.abs(tbc(a) - tbc(b)).max()),      # arka-ucun kullandığı bileşke
                T_BS_enbuyuk=float(np.abs(np.asarray(a["T_BS"]) - np.asarray(b["T_BS"])).max()),
                R1_enbuyuk=float(np.abs(np.asarray(a["R1"]) - np.asarray(b["R1"])).max()))


# ----------------------------------------------------------------------------- arka-uç sınamaları
def ana_ref(seq):
    ana = a3s_kos.ana_izgara()[(seq, 0, 8.0)]
    return {x[1]: v for x, v in ana.items() if x != "fe" and x[0] == "A0"}


def ontest_ref():
    ref = {}
    with open(ONTEST) as f:
        for line in f:
            r = json.loads(line)
            c = r.get("cell") or {}
            if r.get("arm") == "A0" and c.get("sequence") == onarim_kos.ad_mod(SEQ_D, "R0", 8, "N"):
                ref[int(r["window"])] = {k: r.get(k) for k in ALAN + ("status",)}
    return ref


def be_kos(G, ad, cell, dokum, ref, not_=""):
    import config
    d = os.path.join(G, "spool", ad)
    os.makedirs(d, exist_ok=True)
    cfg = _build_cfg(cell, d)
    dp = config.dump_path(cfg)
    with open(dp, "wb") as f:
        pickle.dump(dokum, f)
    meta = dict(cell=cell, ok=True, dir=d, dump=dp, gt=config.gt_path(cfg), t_frontend=None, stats={}, gpu=None)
    t0 = time.time()
    recs = _backend_worker(meta)
    satir, hepsi_bit = [], True
    for r in recs:
        if "arm" not in r:
            yaz(dict(sinama=ad, hata=r.get("error")))
            print(f"[{ad}] HÜCRE HATASI: {r.get('error')}")
            return 2
        m = ref.get(int(r["window"])) or {}
        bit = all(r.get(k) == m.get(k) for k in ALAN)
        hepsi_bit &= bit
        yuzde = lambda k: (None if not (r.get(k) and m.get(k)) else 100.0 * (r[k] / m[k] - 1.0))
        satir.append(dict(pencere=int(r["window"]), status=r.get("status"), bit_ayni=bit,
                          olculen={k: r.get(k) for k in ALAN}, referans={k: m.get(k) for k in ALAN},
                          ate_fark_yuzde=yuzde("m_ate_se3"), anees_poz_fark_yuzde=yuzde("n_nees_pos_mean")))
    en = max((abs(s["ate_fark_yuzde"]) for s in satir if s["ate_fark_yuzde"] is not None), default=None)
    yaz(dict(sinama=ad, hucre=cell.key(), not_=not_, sure_s=round(time.time() - t0, 1), hepsi_bit_ayni=hepsi_bit,
             en_buyuk_ate_fark_yuzde=en, pencereler=satir))
    print(f"[{ad}] {cell.key()}  ({time.time() - t0:.0f} s)  hepsi bit-aynı: {hepsi_bit}")
    yz = lambda v: "—" if v is None else f"{v:+.2f} %"
    for s in satir:
        o, m = s["olculen"], s["referans"]
        print(f"   pencere {s['pencere']:>4}  ATE {o['m_ate_se3']:.6f} / {m.get('m_ate_se3')}  "
              f"ANEES_poz {o['n_nees_pos_mean']:.3f} / {m.get('n_nees_pos_mean')}  bit-aynı {s['bit_ayni']}  "
              f"ATE farkı {yz(s['ate_fark_yuzde'])}  ANEES farkı {yz(s['anees_poz_fark_yuzde'])}")
    return 0


def be(G, ad):
    ana_hucre = lambda s: a3s_kos._hucre(s, 0, 8, ("A0",), a3s_kos.K_A3S, "g3")
    if ad == "T-K":
        return be_kos(G, ad, ana_hucre("MH01"), kanca_dokumu(G, "MH01"), ana_ref("MH01"), "kanca dökümü değişmeden")
    if ad == "T-D2":
        dk = kanca_dokumu(G, "MH01")
        kd = kalibrasyon(onarim_kos.KOK, onarim_kos.ad_mod("MH01", "R0", 8, "N"))
        fark = kalib_fark(dk, kd)
        for k in KALIB:
            dk[k] = kd[k]
        return be_kos(G, ad, ana_hucre("MH01"), dk, ana_ref("MH01"), "dosya yolu kalibrasyonu: " + json.dumps(fark))
    if ad == "T-E":
        dk = kanca_dokumu(G, "MH01")
        fx, fy, cx, cy = dk["rect_K"]
        dk["rect_K"] = [fx * (1.0 + 1e-12), fy * (1.0 + 1e-12), cx, cy]
        assert dk["rect_K"][0] != fx
        return be_kos(G, ad, ana_hucre("MH01"), dk, ana_ref("MH01"),
                      f"fx, fy × (1 + 1e-12): {fx!r} -> {dk['rect_K'][0]!r}")
    if ad == "T-D1k":
        dk = kanca_dokumu(G, SEQ_D)
        dd = kalibrasyon(onarim_kos.KOK, onarim_kos.ad_mod(SEQ_D, "R0", 8, "N"))
        fark = kalib_fark(dk, dd)
        for k in KALIB:
            dk[k] = dd[k]
        return be_kos(G, ad, ana_hucre(SEQ_D), dk, ana_ref(SEQ_D), "dosya yolu kalibrasyonu: " + json.dumps(fark))
    if ad == "T-D1":
        cell = dataclasses.replace(onarim_kos.hucre(SEQ_D, "R0", 8, "N"), arms=("A0",))
        return be_kos(G, ad, cell, dosya_dokumu(G), ontest_ref(), "dosya yolu dökümü yeniden koşuldu")
    raise SystemExit(f"bilinmeyen sınama: {ad}")


# ----------------------------------------------------------------------------- T-F: ön-uç + döküm karşılaştırması
def bos_bellek_mib():
    o = subprocess.run(["nvidia-smi", "--query-gpu=memory.total,memory.used", "--format=csv,noheader,nounits"],
                       capture_output=True, text=True).stdout.strip().splitlines()
    for l in o:                                           # nvidia-smi iki kartı da listeler; 3090 = 24 GB satırı
        tot, kul = (int(x) for x in l.split(","))
        if tot > 20000:
            return tot - kul
    raise RuntimeError("RTX 3090 satırı bulunamadı")


def fe(G, esik_mib, bekle_s):
    cell = onarim_kos.hucre(SEQ_D, "R0", 8, "N")
    t0 = time.time()
    while bos_bellek_mib() < esik_mib:
        if time.time() - t0 > bekle_s:
            print(f"[fe] {bekle_s} s içinde {esik_mib} MiB boş bellek oluşmadı; çıkılıyor")
            return 3
        time.sleep(15)
    print(f"[fe] boş bellek {bos_bellek_mib()} MiB (bekleme {time.time() - t0:.0f} s); ön-uç başlıyor", flush=True)
    meta = _frontend_worker((cell, os.path.join(G, "spool_fe")))
    if not meta.get("ok"):
        yaz(dict(sinama="T-F/fe", hata=meta.get("error")))
        print("[fe] HATA:", meta.get("error"))
        return 2
    hedef = os.path.join(G, "dumps", f"DOSYA__{SEQ_D}.pkl")
    os.replace(meta["dump"], hedef)
    shutil.rmtree(meta["dir"], ignore_errors=True)
    ref = None
    with open(ONTEST) as f:
        for line in f:
            r = json.loads(line)
            if (r.get("cell") or {}).get("sequence") == cell.sequence and r.get("frontend_stats"):
                ref = r["frontend_stats"]
                break
    ayni = a3s_kos.fe_ayni(meta["stats"], ref)
    yaz(dict(sinama="T-F/fe", hucre=cell.key(), sure_s=meta.get("t_frontend"), gpu=meta.get("gpu"),
             istatistik_ontestle_ayni=ayni, istatistik=meta["stats"]))
    print(f"[fe] tamam ({meta.get('t_frontend'):.0f} s, {meta.get('gpu')}); istatistikler ön testle aynı: {ayni}; -> {hedef}")
    return 0 if ayni else 2


def fark(G):
    a, b = kanca_dokumu(G, SEQ_D), dosya_dokumu(G)
    fa, fb = a["frames"], b["frames"]
    S = dict(sinama="T-F/fark", anahtarlar_ayni=sorted(a) == sorted(b), n_kare=[len(fa), len(fb)],
             has_keypoints=[a.get("has_keypoints"), b.get("has_keypoints")], kalibrasyon=kalib_fark(a, b))
    # işlemcide hesaplanan kalibrasyon dökümdekiyle aynı mı (T-D2 / T-D1k bu hesabı kullanır)
    ik = kalibrasyon(a3s_kos.EUROC_FAZ2, SEQ_D)
    id_ = kalibrasyon(onarim_kos.KOK, onarim_kos.ad_mod(SEQ_D, "R0", 8, "N"))
    esit = lambda x, y: all(np.array_equal(np.asarray(x[k], float), np.asarray(y[k], float)) for k in KALIB)
    S["islemci_kalibrasyonu_dokumle_ayni"] = dict(kanca=esit(ik, a), dosya=esit(id_, b))
    n = min(len(fa), len(fb))
    S["zaman_damgasi_ayni"] = all(fa[i]["timestamp"] == fb[i]["timestamp"] for i in range(n)) and len(fa) == len(fb)
    S["kare_anahtarlari_ayni"] = all(sorted(fa[i]) == sorted(fb[i]) for i in range(n))
    sekil = ayrik = nan = 0
    ek_fark = {}
    goreli, n_goz, n_xyz, n_xyz_bit = [], 0, 0, 0
    for i in range(n):
        oa, ob = np.asarray(fa[i]["obs"]), np.asarray(fb[i]["obs"])
        if oa.shape != ob.shape:
            sekil += 1
            continue
        n_goz += oa.shape[0]
        if not np.array_equal(oa[:, :5], ob[:, :5]):
            ayrik += 1
        na, nb = np.isnan(oa[:, 5:8]), np.isnan(ob[:, 5:8])
        if not np.array_equal(na, nb):
            nan += 1
            continue
        v = ~na.any(axis=1)
        pa, pb = oa[v, 5:8], ob[v, 5:8]
        n_xyz += int(v.sum())
        n_xyz_bit += int((pa == pb).all(axis=1).sum())
        if len(pa):
            goreli.append(np.linalg.norm(pb - pa, axis=1) / np.linalg.norm(pa, axis=1))
        for k in fa[i]:
            if k in ("timestamp", "obs") or k not in fb[i]:
                continue
            xa, xb = fa[i][k], fb[i][k]
            esit = (np.array_equal(np.asarray(xa), np.asarray(xb), equal_nan=True)
                    if isinstance(xa, (np.ndarray, list, tuple, float, int)) else xa == xb)
            if not esit:
                ek_fark[k] = ek_fark.get(k, 0) + 1
    g = np.concatenate(goreli) if goreli else np.zeros(0)
    S.update(sekli_farkli_kare=sekil, ayrik_sutunu_farkli_kare=ayrik, nan_deseni_farkli_kare=nan,
             n_gozlem=n_goz, n_xyz=n_xyz, n_xyz_bit_ayni=n_xyz_bit, ek_alan_farkli_kare=ek_fark,
             xyz_goreli_fark=dict(en_buyuk=float(g.max()) if len(g) else None,
                                  medyan=float(np.median(g)) if len(g) else None,
                                  y99=float(np.percentile(g, 99)) if len(g) else None))
    S["gecti"] = bool(S["anahtarlar_ayni"] and S["zaman_damgasi_ayni"] and S["kare_anahtarlari_ayni"] and sekil == 0
                      and ayrik == 0 and nan == 0 and len(g) and g.max() < 1e-6)
    yaz(S)
    print(json.dumps(S, indent=1, ensure_ascii=False, default=_js))
    return 0 if S["gecti"] else 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dizin", required=True)
    ap.add_argument("adim", choices=("fe", "fark", "be"))
    ap.add_argument("sinama", nargs="?")
    ap.add_argument("--esik_mib", type=int, default=6500)
    ap.add_argument("--bekle_s", type=int, default=7200)
    a = ap.parse_args()
    assert code_fingerprint() == a3s_kos.BEKLENEN_FP, code_fingerprint()
    if a.adim == "fe":
        return fe(a.dizin, a.esik_mib, a.bekle_s)
    if a.adim == "fark":
        return fark(a.dizin)
    return be(a.dizin, a.sinama)


if __name__ == "__main__":
    sys.exit(main())
