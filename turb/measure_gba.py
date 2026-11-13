"""
GBA / loop-closure MALİYET ÖLÇÜMÜ — "kaç saat ekliyor, karşılığında ne veriyor".

Aynı ön-uç dump'ı üzerinde arka-uç üç yapılandırmayla koşulur ve HEM süre HEM
doğruluk ölçülür. Maliyet tek başına anlamsızdır: soru "GBA pahalı mı" değil,
"ödediği süre karşılığında ATE'yi düzeltiyor mu"dur.

    kapali   : loop.enable = false                  (loop closure hiç yok)
    pgo      : loop.enable = true,  gba = false      (yalnız poz-graf optimizasyonu)
    pgo_gba  : loop.enable = true,  gba = true       (şu anki koşum ayarı)

Ölçüm KOŞULLARI (adil karşılaştırma için zorunlu):
  * aynı dump, aynı sekans, aynı kare aralığı, aynı kol;
  * makine BOŞ olmalı (kuyruk SIGSTOP ile dondurulmuş) — aksi halde ölçülen süre
    yarışan işçilerin yükünü içerir ve oran bozulur;
  * her yapılandırma `--repeat` kez koşulur, MEDYAN alınır.

    ./turbvio/bin/python -m turb.measure_gba --sequences MH01 MH03 --arms A0 A3
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import statistics
import tempfile
import time

CFGS = {
    "kapali":  {"enable": False, "pgo": False, "gba": False},
    "pgo":     {"enable": True,  "pgo": True,  "gba": False},
    "pgo_gba": {"enable": True,  "pgo": True,  "gba": True},
}


def _find_dump(cache, seq, dr):
    pat = os.path.join(cache, f"euroc_{seq}_dr{dr:g}_s0_*.pkl")
    hits = [p for p in glob.glob(pat) if not p.endswith("_quality.json")]
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser(description="GBA/loop-closure maliyet ölçümü")
    ap.add_argument("--sequences", nargs="+", default=["MH01"])
    ap.add_argument("--arms", nargs="+", default=["A0"])
    ap.add_argument("--D_over_r0", type=float, default=0.0)
    ap.add_argument("--max_frames", type=int, default=2000)
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--cache", default="/data/turb_dump_cache")
    ap.add_argument("--out", default="evidence/gba_maliyeti.json")
    a = ap.parse_args()

    import config
    from gtsam_backend import main as backend_main
    from turb.metrics import trajectory_metrics

    rows = []
    tmp = tempfile.mkdtemp(prefix="gbacost_")
    try:
        for seq in a.sequences:
            dump = _find_dump(a.cache, seq, a.D_over_r0)
            if dump is None:
                print(f"[gba] {seq}: önbellekte dump yok, atlanıyor")
                continue
            print(f"[gba] {seq}: {os.path.basename(dump)}")
            for arm in a.arms:
                for name, lp in CFGS.items():
                    for rep in range(a.repeat):
                        cfg = config.load()
                        cfg["dataset"].update(sequence=seq, base_path=os.environ.get(
                            "TURB_EUROC", "/data/sdata/"),
                            max_frames=a.max_frames)
                        cfg["paths"]["outputs"] = tmp
                        cfg["loop"].update(lp)
                        out = os.path.join(tmp, f"{seq}_{arm}_{name}_{rep}.csv")
                        argv = ["--dump", dump, "--out", out,
                                "--max_frames", str(a.max_frames)]
                        if arm != "A0":
                            argv += ["--turb_arm", arm,
                                     "--turb_D_over_r0", str(max(a.D_over_r0, 5.0))]
                        t0 = time.time()
                        err = None
                        try:
                            backend_main(config.backend_args(cfg, argv=argv))
                        except BaseException as e:
                            err = f"{type(e).__name__}: {e}"
                        dt = time.time() - t0
                        m = {}
                        if err is None and os.path.exists(out):
                            m = trajectory_metrics(out, config.gt_path(cfg))
                        rows.append({"seq": seq, "arm": arm, "cfg": name, "rep": rep,
                                     "saniye": dt, "ate_se3": m.get("ate_se3"),
                                     "ate_noalign": m.get("ate_noalign"),
                                     "rpe_trans": m.get("rpe_trans"),
                                     "diverged": m.get("diverged"), "hata": err})
                        print(f"   {seq} {arm:4s} {name:8s} rep{rep}: {dt:7.1f} s  "
                              f"ATE={m.get('ate_se3', float('nan')):.4f} m"
                              + (f"  HATA {err}" if err else ""), flush=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(rows, f, indent=1)

    # ------------------------------- özet -------------------------------
    def med(sel, key):
        v = [r[key] for r in rows if sel(r) and r.get(key) is not None
             and r.get("hata") is None]
        return statistics.median(v) if v else float("nan")

    print("\n================ ÖZET ================")
    print(f"{'kol':5s} {'yapılandırma':12s} {'süre(s)':>9s} {'süre farkı':>11s} "
          f"{'ATE(m)':>9s} {'ATE farkı':>10s}")
    base_by_arm = {}
    for arm in a.arms:
        base_t = med(lambda r: r["arm"] == arm and r["cfg"] == "kapali", "saniye")
        base_a = med(lambda r: r["arm"] == arm and r["cfg"] == "kapali", "ate_se3")
        base_by_arm[arm] = (base_t, base_a)
        for name in CFGS:
            t = med(lambda r: r["arm"] == arm and r["cfg"] == name, "saniye")
            e = med(lambda r: r["arm"] == arm and r["cfg"] == name, "ate_se3")
            print(f"{arm:5s} {name:12s} {t:9.1f} {t - base_t:+11.1f} "
                  f"{e:9.4f} {e - base_a:+10.4f}")

    # GBA'nın SAF maliyeti = pgo_gba - pgo
    print("\n--- GBA'nın saf ek maliyeti (pgo_gba - pgo) ---")
    for arm in a.arms:
        t_g = med(lambda r: r["arm"] == arm and r["cfg"] == "pgo_gba", "saniye")
        t_p = med(lambda r: r["arm"] == arm and r["cfg"] == "pgo", "saniye")
        e_g = med(lambda r: r["arm"] == arm and r["cfg"] == "pgo_gba", "ate_se3")
        e_p = med(lambda r: r["arm"] == arm and r["cfg"] == "pgo", "ate_se3")
        pct = 100 * (t_g - t_p) / t_p if t_p else float("nan")
        print(f"  {arm}: +{t_g - t_p:.1f} s/koşum (%{pct:+.1f})  ·  "
              f"ATE {e_p:.4f} -> {e_g:.4f} ({e_g - e_p:+.4f} m)")
    print(f"\n({a.out})")


if __name__ == "__main__":
    main()
