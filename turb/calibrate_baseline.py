"""
A0C kalibrasyonu — temiz veride NEES'i tutarlılık bandına oturtan katsayı.

Hattın gürültü modeli ATE için optimize edildiğinden, TÜRBÜLANS OLMASA BİLE
kestirici aşırı-güvenlidir (ölçüldü: pencere-göreli NEES ~ 1022, bant [0.35, 7.81]).
Bu betik, izotropik sigma'ları k katsayısıyla ölçekleyip TEMİZ sekanslarda NEES'i
ölçer ve bandın ortasına (~3) en yakın k'yı seçer.

Sonuç iki şekilde kullanılır:
  (1) A0C kolunun `cal_factor`'ü olarak — adil kalibre edilmiş kontrol temel çizgisi;
  (2) BULGU olarak: k çok büyük çıkarsa, "skaler şişirme tutarlılığı ancak doğruluğu
      feda ederek satın alabilir" demektir ve bu, korelasyon YAPISININ gerekliliğinin
      doğrudan kanıtıdır.

    ./turbvio/bin/python -m turb.calibrate_baseline --factors 1 2 4 8 16 32
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile

from .experiments import EUROC, EUROC_SEQS
from .runner import Cell, _build_cfg, _fe_key, _cache_dir


def main():
    ap = argparse.ArgumentParser(description="A0C temel çizgi kalibrasyonu")
    ap.add_argument("--factors", nargs="+", type=float,
                    default=[1, 2, 4, 8, 16, 32, 64])
    ap.add_argument("--sequences", nargs="+", default=["MH01", "MH03", "MH05"])
    ap.add_argument("--out", default="evidence/A0C_calibration.json")
    a = ap.parse_args()

    import config
    from gtsam_backend import main as backend_main
    from turb.metrics import trajectory_metrics, nees_from_dump
    import pickle

    cache = _cache_dir()
    rows = []
    for seq in a.sequences:
        cell = Cell(dataset="euroc", base_path=EUROC, sequence=seq, seed=0,
                    D_over_r0=0.0, arms=("A0C",), nees=True, cov_every=10)
        cpath = os.path.join(cache, _fe_key(cell) + ".pkl")
        if not os.path.exists(cpath):
            print(f"[kal] {seq}: önbellekte dump yok, atlanıyor")
            continue
        tmp = tempfile.mkdtemp(prefix="a0ccal_")
        try:
            cfg = _build_cfg(cell, tmp)
            shutil.copy(cpath, config.dump_path(cfg))
            gt = config.gt_path(cfg)
            with open(config.dump_path(cfg), "rb") as f:
                f_ts = [int(fr["timestamp"]) for fr in pickle.load(f)["frames"]]
            for k in a.factors:
                out = os.path.join(tmp, f"t{k:g}.csv")
                cov = os.path.join(tmp, f"c{k:g}.npz")
                argv = ["--out", out, "--start_frame", "0",
                        "--max_frames", str(int(max(cell.windows) + cell.window_len)),
                        "--turb_arm", "A0C", "--turb_cal_factor", str(k),
                        "--cov_out", cov, "--cov_every", "10"]
                try:
                    backend_main(config.backend_args(cfg, argv=argv))
                except BaseException as e:
                    print(f"[kal] {seq} k={k:g}: HATA {type(e).__name__}")
                    rows.append({"seq": seq, "k": k, "status": "fail"})
                    continue
                for w in cell.windows:
                    i0, i1 = int(w), int(min(w + cell.window_len, len(f_ts)) - 1)
                    t0, t1 = float(f_ts[i0]), float(f_ts[i1])
                    m = trajectory_metrics(out, gt, t_min=t0, t_max=t1)
                    n = nees_from_dump(cov, gt, t_min=t0, t_max=t1) \
                        if os.path.exists(cov) else {}
                    rows.append({"seq": seq, "k": k, "window": w, "status": "ok",
                                 "ate": m.get("ate_se3"),
                                 "nees": n.get("nees_pos_median"),
                                 "cov": n.get("coverage_pos")})
                r = [x for x in rows if x.get("seq") == seq and x.get("k") == k
                     and x.get("status") == "ok"]
                if r:
                    import numpy as np
                    print(f"[kal] {seq} k={k:5g}: ATE={np.nanmean([x['ate'] for x in r]):.4f} m  "
                          f"NEES_med={np.nanmedian([x['nees'] for x in r]):9.2f}  "
                          f"kapsama={np.nanmean([x['cov'] for x in r]):.2f}", flush=True)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(rows, f, indent=1)

    # bandın ortasına (3.0) en yakın k
    import numpy as np
    # SEÇİM ÖLÇÜTÜ: "banda giren EN KÜÇÜK k".
    # Bandın ortasına (3.0) en yakın k'yı seçmek yanlıştır: kabul ölçütü BANTTIR,
    # banda girdikten sonra daha derine inmek istatistiksel kazanç sağlamaz ama
    # doğruluğu harcar. Ölçüldü (MH01): k=16 -> NEES 6.31 / ATE 0.0496;
    # k=32 -> NEES 4.04 / ATE 0.1420. İkisi de bandın içinde, ama k=32 üç kat
    # daha pahalı. Kalibre temel çizgi, kontrolü HAKSIZ yere zayıflatmamalıdır.
    LO, HI = 0.35, 7.81
    best = None
    for k in sorted(a.factors):
        v = [x["nees"] for x in rows if x.get("k") == k and x.get("status") == "ok"
             and x.get("nees") is not None and np.isfinite(x.get("nees"))]
        if not v:
            continue
        med = float(np.median(v))
        at = [x["ate"] for x in rows if x.get("k") == k and x.get("status") == "ok"
              and x.get("ate") is not None]
        print(f"  k={k:5g} -> NEES medyan {med:9.2f}  ATE {np.nanmean(at):.4f} m"
              + ("   <-- BANT İÇİ" if LO <= med <= HI else ""))
        if best is None and LO <= med <= HI:
            best = k
    print(f"\n✓ Önerilen cal_factor = {best}  ({a.out})")


if __name__ == "__main__":
    main()
