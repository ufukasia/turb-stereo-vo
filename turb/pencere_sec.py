"""
Faz 2 pencere seçimi — YALNIZ GT hareketine bakar (ON_KAYIT_FAZ2 §1).

    python -m turb.pencere_sec                     # evidence/pencereler_faz2.json yazar
    python -m turb.pencere_sec --kontrol           # yazmadan yazdır

NEDEN: eski sabit başlangıçlar (500/900/1300/1600) MH01/500'de neredeyse hareketsiz bir
pencere veriyordu (0.075 m; denetim B6) ve 1300/1600 pencereleri 100 kare örtüşüyordu.

KURAL (deterministik, koşum sonucu görmez):
  1. Kare indeksi = sıralı cam0 zaman damgaları (`stereo_front_end.process_sequence` ile aynı
     sıralama; ön-uç dökümündeki kare k bu listenin k. elemanıdır).
  2. Aday pencere [s, s+L) (L = 400) geçerlidir ancak ve ancak:
       a) s >= ISINMA (100 kare): arka-uç kare 0'dan koşar, ilk 100 kare ortak ısınmadır;
       b) pencerenin İLK ve SON karesinin zaman damgası GT kapsamı içindedir;
       c) GT yol uzunluğu >= 5.0 m. Yol uzunluğu, GT konumunun pencerenin KAMERA zaman
          damgalarına doğrusal enterpolasyonuyla (20 Hz) elde edilen polilinin uzunluğudur
          (`turb.metrics.trajectory_metrics` içindeki `path_length_gt` ile aynı tanım;
          200 Hz GT titreşimini yol uzunluğuna katmamak için kamera hızında örneklenir).
       d) (SEÇENEKLİ, varsayılan kapalı = ön-kayıt metni) GT konumlarının eksen-hizalı sınır
          kutusu köşegeni ("kapsam") >= --min_kapsam [m]. Gerekçe: MH01/MH02'nin başındaki IMU
          uyarım sallaması 5 m yolu yerinde salınarak (kapsam < 1 m) doldurur.
  3. Açgözlü, en erken başlangıç: s = ISINMA'dan başlayarak birer kare ilerlenir; ilk geçerli
     s kabul edilir, sonraki arama s + L'den (örtüşme yok) devam eder. En çok 3 pencere.
  4. Hiçbir koşum çıktısına, ATE'ye veya NEES'e bakılmaz; girdi yalnız cam0/data.csv ve
     state_groundtruth_estimate0/data.csv.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SEKANSLAR = ["MH01", "MH02", "MH03", "MH04", "MH05",
             "V1_01_easy", "V1_02_medium", "V2_01_easy", "V2_02_medium"]
PENCERE_UZ = 400
ISINMA = 100
MIN_YOL = 5.0
MIN_KAPSAM = 0.0        # ön-kayıt kuralı: kapsam ölçütü YOK; seçenek olarak --min_kapsam 1.0
EN_COK = 3
VARSAYILAN_KOK = os.environ.get("SLAM_DATA_ROOT", "/data/sdata")
CIKTI = os.path.join(ROOT, "evidence", "pencereler_faz2.json")


def _sha1(path: str) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def kare_zamanlari(kok: str, seq: str) -> np.ndarray:
    """Ön-uçla aynı sıra: cam0 zaman damgaları (dize olarak sıralı, sonra int)."""
    df = pd.read_csv(os.path.join(kok, seq, "mav0", "cam0", "data.csv"))
    df.columns = [c.strip() for c in df.columns]
    ts = sorted(set(df["#timestamp [ns]"].astype(str)))
    return np.array([int(t) for t in ts], dtype=np.int64)


def gt_konum(kok: str, seq: str):
    from .metrics import load_gt
    t, xyz, _ = load_gt(os.path.join(kok, seq, "mav0", "state_groundtruth_estimate0", "data.csv"))
    return t, xyz


def yol_uzunlugu(f_ts, gt_t, gt_xyz, i0: int, i1: int) -> float:
    """[i0, i1) karelerinin zaman damgalarında GT polilinin uzunluğu [m]."""
    tt = f_ts[i0:i1].astype(np.float64)
    p = np.stack([np.interp(tt, gt_t, gt_xyz[:, k]) for k in range(3)], axis=1)
    return float(np.sum(np.linalg.norm(np.diff(p, axis=0), axis=1)))


def kapsam(f_ts, gt_t, gt_xyz, i0: int, i1: int) -> float:
    """Pencere GT konumlarının sınır kutusu köşegeni [m]."""
    tt = f_ts[i0:i1].astype(np.float64)
    p = np.stack([np.interp(tt, gt_t, gt_xyz[:, k]) for k in range(3)], axis=1)
    return float(np.linalg.norm(p.max(0) - p.min(0)))


def sec(kok: str, seq: str, L: int = PENCERE_UZ, isinma: int = ISINMA,
        min_yol: float = MIN_YOL, en_cok: int = EN_COK, min_kapsam: float = MIN_KAPSAM):
    f_ts = kare_zamanlari(kok, seq)
    gt_t, gt_xyz = gt_konum(kok, seq)
    n = len(f_ts)
    secilen, s = [], int(isinma)
    while s + L <= n and len(secilen) < en_cok:
        t0, t1 = float(f_ts[s]), float(f_ts[s + L - 1])
        if t0 >= gt_t[0] and t1 <= gt_t[-1]:
            y = yol_uzunlugu(f_ts, gt_t, gt_xyz, s, s + L)
            kp = kapsam(f_ts, gt_t, gt_xyz, s, s + L)
            if y >= min_yol and kp >= min_kapsam:
                secilen.append({"start": int(s), "end_excl": int(s + L), "path_m": round(y, 4),
                                "extent_m": round(kp, 4),
                                "t_start_ns": int(f_ts[s]), "t_end_ns": int(f_ts[s + L - 1])})
                s += L
                continue
        s += 1
    # Bilgi amaçlı: sekansın bütün geçerli 400 karelik başlangıçlarında yol dağılımı
    yollar = [yol_uzunlugu(f_ts, gt_t, gt_xyz, k, k + L)
              for k in range(int(isinma), n - L + 1, 20)
              if f_ts[k] >= gt_t[0] and f_ts[k + L - 1] <= gt_t[-1]]
    return {
        "n_frames": int(n),
        "windows": secilen,
        "starts": [w["start"] for w in secilen],
        "path_max_m": round(float(max(yollar)), 3) if yollar else None,
        "path_median_m": round(float(np.median(yollar)), 3) if yollar else None,
        "cam0_sha1": _sha1(os.path.join(kok, seq, "mav0", "cam0", "data.csv")),
        "gt_sha1": _sha1(os.path.join(kok, seq, "mav0", "state_groundtruth_estimate0", "data.csv")),
    }


def yukle(path: str = CIKTI) -> dict:
    """{sekans: (başlangıç, ...)} — runner/deney tanımları bunu kullanır."""
    with open(path) as f:
        d = json.load(f)
    return {s: tuple(int(x) for x in v["starts"]) for s, v in d["sequences"].items()}


def main():
    ap = argparse.ArgumentParser(description="Faz 2 GT-hareket pencere seçimi")
    ap.add_argument("--kok", default=VARSAYILAN_KOK)
    ap.add_argument("--out", default=CIKTI)
    ap.add_argument("--kontrol", action="store_true", help="dosya yazma, yalnız yazdır")
    ap.add_argument("--min_kapsam", type=float, default=MIN_KAPSAM,
                    help="ek kapsam ölçütü [m] (varsayılan 0 = ön-kayıt metni)")
    a = ap.parse_args()
    out = {"rule": {"window_len": PENCERE_UZ, "warmup_min_start": ISINMA, "min_path_m": MIN_YOL,
                    "min_extent_m": a.min_kapsam,
                    "max_windows": EN_COK, "selection": "greedy earliest start, step 1 frame, "
                    "non-overlapping (next search from s+L); window first/last frame inside GT; "
                    "path = GT position interpolated at cam0 timestamps"},
           "generated": f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}",
           "script_sha1": _sha1(os.path.abspath(__file__)),
           "data_root": a.kok, "sequences": {}}
    for s in SEKANSLAR:
        r = sec(a.kok, s, min_kapsam=a.min_kapsam)
        out["sequences"][s] = r
        print(f"{s:14s} n={r['n_frames']:5d}  pencereler="
              + ", ".join(f"{w['start']}({w['path_m']:.2f} m, kapsam {w['extent_m']:.2f})" for w in r["windows"])
              + f"   [400-kare yol: medyan {r['path_median_m']} m, maks {r['path_max_m']} m]")
    if not a.kontrol:
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        with open(a.out, "w") as f:
            json.dump(out, f, indent=1)
        print(f"-> {a.out}")


if __name__ == "__main__":
    main()
