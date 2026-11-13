"""
R4 "etki uygulandı" denetimi: üretilen sekansın cam0 karelerinin rektifiye TEMİZ kareye göre PSNR'ı
(her 50. kare). Kullanım:  python -m turb.rakip_etki <sekans dizini> ... --out <json>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import cv2
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from turb.uret_sekans import asl_oku, Rektifikator  # noqa: E402


def psnr(a, b):
    m = np.mean((a.astype(np.float64) - b.astype(np.float64)) ** 2)
    return float("inf") if m == 0 else float(10 * np.log10(255.0 ** 2 / m))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dizinler", nargs="+")
    ap.add_argument("--adim", type=int, default=50)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = {}
    for d in map(Path, a.dizinler):
        meta = json.loads((d / "turb_meta.json").read_text())
        src = Path(meta["kaynak"]["yol"])
        cal = asl_oku(src)
        sd = "sensor_yaml" if "T_BS" in meta["rektifikasyon"]["taban_matrisi_kaynagi"] else "gomulu"
        rk = Rektifikator(cal, "rektifiye", sd)       # üretimdeki rektifikasyonun aynısı
        ts_out = sorted(int(l.split(",")[0]) for l in open(d / "mav0/cam0/data.csv") if l[0].isdigit())
        idx = {t: i for i, t in enumerate(cal["ts"])}
        ps = []
        for t in ts_out[::a.adim]:
            i0 = cv2.imread(str(src / "mav0/cam0/data" / cal["f0"][t]), cv2.IMREAD_GRAYSCALE)
            i1 = cv2.imread(str(src / "mav0/cam1/data" / cal["f1"][t]), cv2.IMREAD_GRAYSCALE)
            r0, _ = rk.uygula(i0, i1)
            g = cv2.imread(str(d / "mav0/cam0/data" / f"{t}.png"), cv2.IMREAD_UNCHANGED)
            ps.append(psnr(g, r0))
        res[d.name] = {"n": len(ps), "psnr_medyan_db": float(np.median(ps)),
                       "psnr_min_db": float(np.min(ps)), "psnr_maks_db": float(np.max(ps)),
                       "D_over_r0": meta["parametreler"]["D_over_r0"],
                       "sigma_dn": (meta.get("sensor_gurultusu") or {}).get("sigma_dn", 0.0)}
        print(d.name, res[d.name], flush=True)
    Path(a.out).write_text(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
