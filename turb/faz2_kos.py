"""
Faz 2 EuRoC uzun koşum başlatıcısı (ON_TEST_FAZ2_EUROC.md §6). TEK KART: RTX 3090.

    # ortam (zorunlu; betik denetler):
    export CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=4
    export SLAM_DATA_ROOT="/data/sdata"
    python -m turb.faz2_kos --kuru                 # hiçbir şey koşmadan denetim + plan
    python -m turb.faz2_kos --jobs 10              # koşum (devam ettirilebilir: tamamlanan hücreler atlanır)
    python -m turb.faz2_kos --denetle              # koşum sırasında/sonunda kayıt denetimi
    python -m turb.faz2_kos --deney FAKT --kuru    # H2 faktöriyel (180 hücre, A0)
    python -m turb.faz2_kos --deney H6 --kirilma 12 --kuru   # H6 nedensel (405 hücre, A0)

Güvenlik kapıları (koşumdan önce, her başlatmada):
  G1  torch'un gördüğü tek kart RTX 3090 değilse DUR (PyEco ön-ucu GPU modelleri arasında bit-aynı değil;
      bir karşılaştırma ailesi tek kartta koşulmalı).
  G2  mevcut kayıtlarda fe_gpu RTX 3090 dışında bir değer taşıyorsa DUR.
  G3  mevcut kayıtlarda güncel code_fp dışında bir parmak izi varsa DUR (--farkli_koda_izin yoksa).
  G4  FAZ2_K_A0C tanımlı, hücre sayısı beklenen (9 x 7 x tohum).
Her başlatma `manifest_*.json` yazar (code_fp, config ve pencere dosyası sha256, k, hücre sayısı, ortam).
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BEKLENEN_GPU = "NVIDIA GeForce RTX 3090"
DENEYLER = {"FAZ2": ("FAZ2", "/tmp/turb_spool_faz2"),
            "FAKT": ("FAZ2_FAKT", "/tmp/turb_spool_faz2_fakt"),
            "H6": ("FAZ2_H6", "/tmp/turb_spool_faz2_h6")}


def _sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _kayitlar(path):
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as f:
        for line in f:
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return out


def denetle(path, sert=False) -> bool:
    """Kayıt denetimi: GPU tekliği, kod tekliği, durum, D/r0=0'da A1==A0, NEES örnek sayısı, temiz ATE."""
    import numpy as np
    from .runner import code_fingerprint
    R = _kayitlar(path)
    print(f"[denetim] {path}: {len(R)} kayıt")
    if not R:
        return True
    ok = True
    gpus = Counter(r.get("fe_gpu") for r in R if r.get("status") != "cell_fail")
    fps = Counter(r.get("code_fp") for r in R)
    st = Counter(r.get("status") for r in R)
    print(f"  fe_gpu: {dict(gpus)}\n  code_fp: {dict(fps)} (güncel {code_fingerprint()})\n  durum: {dict(st)}")
    if set(gpus) - {BEKLENEN_GPU}:
        print("  !! fe_gpu tek kart değil"); ok = False
    if len(fps) > 1:
        print("  !! birden çok code_fp"); ok = False
    bad = st.get("fail", 0) + st.get("cell_fail", 0)
    if bad / max(len(R), 1) > 0.05:
        print(f"  !! başarısız oranı {bad/len(R):.1%} > %5"); ok = False
    # D/r0=0: A1 == A0 bit-aynı (sigma_oracle = 0; gürültü iki kolda aynı döküm)
    by = defaultdict(dict)
    for r in R:
        if r.get("status") == "ok":
            c = r["cell"]
            by[(c["sequence"], c["D_over_r0"], c["seed"], r["window"])][r["arm"]] = r
    pairs = [(v["A0"], v["A1"]) for k, v in by.items() if k[1] == 0 and "A0" in v and "A1" in v]
    neq = [a["key"] for a, b in pairs if any(a.get(f) != b.get(f) for f in
                                            ("m_ate_se3", "ms_ate_se3", "n_nees_pos_median"))]
    print(f"  D/r0=0 A1==A0: {len(pairs) - len(neq)}/{len(pairs)} bit-aynı")
    if neq:
        print("  !! A1 != A0 (D/r0=0):", neq[:3]); ok = False
    nn = [r["n_n"] for r in R if r.get("n_valid") and r.get("n_n")]
    if nn:
        print(f"  NEES örnek/pencere: medyan {np.median(nn):.0f}, min {min(nn)}, maks {max(nn)}")
        if np.median(nn) < 30:
            print("  !! NEES örnek sayısı düşük"); ok = False
    clean = [r["m_ate_se3"] for r in R if r.get("status") == "ok" and r["cell"]["D_over_r0"] == 0
             and r["arm"] == "A0"]
    if clean:
        print(f"  temiz A0 ATE: medyan {np.median(clean):.4f} m, maks {max(clean):.4f} m")
        if max(clean) > 0.15:
            print("  !! temiz ATE > 0.15 m"); ok = False
    q = [(r["cell"]["D_over_r0"], (r.get("quality") or {}).get("psnr")) for r in R if r.get("status") == "ok"]
    turb_noq = [x for x in q if x[0] > 0 and x[1] is None]
    if turb_noq:
        print(f"  !! türbülanslı {len(turb_noq)} kayıtta PSNR yok (türbülans uygulanmamış olabilir)"); ok = False
    div = sum(1 for r in R if r.get("m_diverged"))
    print(f"  ıraksama (max_err > 5 m): {div}")
    print(f"[denetim] {'TAMAM' if ok else 'SORUN VAR — koşumu durdurup incele'}")
    if sert and not ok:
        sys.exit(2)
    return ok


def main():
    ap = argparse.ArgumentParser(description="Faz 2 EuRoC uzun koşum (tek kart RTX 3090)")
    ap.add_argument("--jobs", type=int, default=10, help="arka-uç (CPU) işçi sayısı")
    ap.add_argument("--gpu_jobs", type=int, default=1, help="eşzamanlı ön-uç (1: bit-tekrarlanabilir)")
    ap.add_argument("--chunk", type=int, default=24)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--deney", default="FAZ2", choices=sorted(DENEYLER))
    ap.add_argument("--kirilma", type=float, default=None, help="H6: ana ızgaradan kırılma noktası (D/r0)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--spool", default=None)
    ap.add_argument("--kuru", action="store_true", help="koşmadan denetle ve planı yaz")
    ap.add_argument("--denetle", action="store_true", help="yalnız kayıt denetimi")
    ap.add_argument("--farkli_koda_izin", action="store_true")
    ap.add_argument("--parca", default=None, help="I/N: yalnız sha1(key) %% N == I hücreleri (K'dan bağımsız, kararlı bölme); "
                    "çıktı records_pI.jsonl, spool _pI")
    a = ap.parse_args()
    _d, _sp = DENEYLER[a.deney]
    _ek = ""
    if a.parca:
        _pi, _pn = (int(x) for x in a.parca.split("/"))
        assert 0 <= _pi < _pn, "--parca I/N: 0 <= I < N"
        _ek = f"_p{_pi}"
    a.out = a.out or os.path.join(ROOT, "evidence", _d, f"records{_ek}.jsonl")
    a.spool = a.spool or (_sp + _ek)
    sys.path.insert(0, ROOT)
    os.chdir(ROOT)

    if a.denetle:
        denetle(a.out)
        return

    # G1 — ortam ve kart
    env = {k: os.environ.get(k) for k in ("CUDA_DEVICE_ORDER", "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS",
                                          "SLAM_DATA_ROOT")}
    import torch
    n = torch.cuda.device_count()
    name = torch.cuda.get_device_name(0) if n else None
    print(f"[G1] ortam {env}; görünen kart sayısı {n}: {name}")
    if n != 1 or name != BEKLENEN_GPU:
        sys.exit(f"[G1] DUR: tek kart {BEKLENEN_GPU} bekleniyordu "
                 f"(CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1)")
    # G4 — ızgara
    from .experiments import (cells_FAZ2, cells_FAZ2_faktoriyel, cells_FAZ2_h6, FAZ2_K_A0C,
                              FAZ2_WINDOWS_FILE, FAZ2_SEQS, FAZ2_DR)
    from .runner import code_fingerprint, run_experiment_pipelined
    import config as _cfg
    if a.deney == "FAZ2":
        cells = cells_FAZ2(a.seeds)
        beklenen = len(FAZ2_SEQS) * len(FAZ2_DR) * a.seeds
    elif a.deney == "FAKT":
        cells = cells_FAZ2_faktoriyel(a.seeds)
        beklenen = len(FAZ2_SEQS) * 2 * 2 * a.seeds
    else:
        if a.kirilma is None:
            sys.exit("[G4] DUR: H6 için --kirilma gerekli (ana ızgaradan, ON_KAYIT §2 tanımı)")
        cells = cells_FAZ2_h6(a.kirilma, a.seeds)
        beklenen = len(FAZ2_SEQS) * 3 * 3 * a.seeds
    if len(cells) != beklenen:
        sys.exit(f"[G4] DUR: {len(cells)} hücre, beklenen {beklenen}")
    if a.parca:
        import hashlib as _hl
        cells = [c for c in cells if int(_hl.sha1(c.key().encode()).hexdigest(), 16) % _pn == _pi]
        print(f"[parça] {a.parca}: {len(cells)} hücre")
    # G2/G3 — mevcut kayıtlar
    R = _kayitlar(a.out)
    gpus = {r.get("fe_gpu") for r in R if r.get("status") != "cell_fail"}
    if gpus - {BEKLENEN_GPU}:
        sys.exit(f"[G2] DUR: mevcut kayıtlarda başka kart: {gpus}")
    cur = code_fingerprint()
    fps = {r.get("code_fp") for r in R}
    if fps - {cur} and not a.farkli_koda_izin:
        sys.exit(f"[G3] DUR: mevcut kayıtlarda farklı code_fp {fps - {cur}} (güncel {cur})")
    done = {(r.get("key"), r.get("window"), r.get("arm")) for r in R}
    todo = [c for c in cells if not all((c.key(), int(w), x) in done for w in c.windows for x in c.arms)]
    man = {"deney": a.deney, "kirilma": a.kirilma, "zaman": f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}", "code_fp": cur, "gpu": name, "env": env,
           "config": str(_cfg.DEFAULT_PATH), "config_sha256": _sha256(str(_cfg.DEFAULT_PATH)),
           "pencere_dosyasi": FAZ2_WINDOWS_FILE, "pencere_sha256": _sha256(FAZ2_WINDOWS_FILE),
           "k_a0c": FAZ2_K_A0C, "hucre": len(cells), "kol_kosumu": sum(len(c.arms) for c in cells),
           "kalan_hucre": len(todo), "jobs": a.jobs, "gpu_jobs": a.gpu_jobs, "chunk": a.chunk,
           "out": a.out, "spool": a.spool, "kuru": a.kuru, "parca": a.parca}
    print(json.dumps(man, indent=1, ensure_ascii=False))
    if a.kuru:
        print("[kuru] koşum BAŞLATILMADI")
        return
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(os.path.join(os.path.dirname(a.out), f"manifest_{datetime.datetime.now():%Y%m%d_%H%M%S}.json"),
              "w") as f:
        json.dump(man, f, indent=1, ensure_ascii=False)
    run_experiment_pipelined(cells, a.out, cpu_jobs=a.jobs, gpu_jobs=a.gpu_jobs, chunk=a.chunk,
                             resume=True, spool=a.spool)
    denetle(a.out)


if __name__ == "__main__":
    main()
