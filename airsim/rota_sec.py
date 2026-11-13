"""
Yükseklik haritasında (harita.py) engelsiz biçerdöver rotası arar ve rota.json yazar.

Rota: L uzunluğunda 2 şerit (gidiş-dönüş), şerit aralığı G, köşe yumuşatmalı (kayit.py ile aynı).
Koşul: üretilen GERÇEK yörüngede (kayit.yorunge, arazi takibi dahil) her pozun çevresindeki
yarıçap r içinde yüzey (DSM, ağaç tepeleri dahil) uçuş yüksekliğinin en az `bosluk` m altında.

    python rota_sec.py --harita africa/harita_africa.npz --irtifa 5 --L 150 --G 12 --cikti africa/rota_F1.json
"""
import argparse
import json
import math
import os
import sys
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kayit import yorunge, _yol_noktalari  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--harita", required=True)
    ap.add_argument("--irtifa", type=float, default=5.0)
    ap.add_argument("--L", type=float, default=150.0)
    ap.add_argument("--G", type=float, default=12.0)
    ap.add_argument("--kose", type=float, default=6.0)
    ap.add_argument("--hiz", type=float, default=3.0)
    ap.add_argument("--yaricap", type=float, default=4.0, help="poz çevresinde engel arama yarıçapı [m]")
    ap.add_argument("--bosluk", type=float, default=2.0, help="engel ile uçuş yüksekliği arası en az [m]")
    ap.add_argument("--cikti", required=True)
    a = ap.parse_args()
    d = np.load(a.harita)
    dsm_h, dtm, x0, y0, res = -d["dsm"], d["dtm"], float(d["x0"]), float(d["y0"]), float(d["res"])
    nx, ny = dsm_h.shape
    # engel haritası: yarıçap içindeki en yüksek yüzey (max filtre)
    import cv2
    k = int(2 * a.yaricap / res) | 1
    dsm_max = cv2.dilate(dsm_h.astype(np.float32), cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k)))
    # DTM'yi rota.json için 1 m'ye seyrelt
    zemin = {"x0": x0, "y0": y0, "res": res * 2, "z": dtm[::2, ::2].round(3).tolist()}

    def dene(xs, ys, th):
        base = _yol_noktalari("tarama", a.L, a.G, 2)          # (0,0)'dan başlar, +x boyunca
        c, s = math.cos(th), math.sin(th)
        pts = np.c_[xs + c * base[:, 0] - s * base[:, 1], ys + s * base[:, 0] + c * base[:, 1]]
        ka = SimpleNamespace(senaryo="dosya", rota=None, hz=20.0, hiz=a.hiz, kose=a.kose, sure=0.0,
                             yaw_hareket=1, tohum=0, irtifa=a.irtifa, irtifa_salinim=0.3, z_zemin=0.0,
                             uzunluk=a.L, genislik=a.G, seritler=2)
        import kayit
        kayit._rota_oku = lambda _: (pts, zemin)
        traj = yorunge(ka)
        P = np.array([p for p, _ in traj])
        gi = np.round((P[:, 0] - x0) / res).astype(int); gj = np.round((P[:, 1] - y0) / res).astype(int)
        if gi.min() < 0 or gj.min() < 0 or gi.max() >= nx or gj.max() >= ny:
            return None
        ucus_h = -P[:, 2]
        marj = ucus_h - dsm_max[gi, gj]                        # >0: engelin üstünde
        agl = ucus_h - (-dtm[gi, gj])
        return {"min_marj": float(marj.min()), "agl_min": float(agl.min()), "agl_max": float(agl.max()),
                "pts": pts.tolist(), "n": len(traj)}

    aday = []
    for th_deg in range(0, 180, 10):
        th = math.radians(th_deg)
        for xs in np.arange(x0 + 20, x0 + nx * res - 20, 10.0):
            for ys in np.arange(y0 + 20, y0 + ny * res - 20, 10.0):
                # hızlı ön eleme: düz çizgi üzerinde engel
                t = np.linspace(0, a.L, 40)
                px = xs + math.cos(th) * t; py = ys + math.sin(th) * t
                gi = np.round((px - x0) / res).astype(int); gj = np.round((py - y0) / res).astype(int)
                if gi.min() < 0 or gj.min() < 0 or gi.max() >= nx or gj.max() >= ny:
                    continue
                if (dsm_max[gi, gj] - (-dtm[gi, gj])).max() > a.irtifa - a.bosluk:
                    continue
                aday.append((xs, ys, th_deg))
    print(f"ön elemeden geçen: {len(aday)}")
    iyi = []
    for xs, ys, th_deg in aday:
        r = dene(xs, ys, math.radians(th_deg))
        if r and r["min_marj"] >= a.bosluk:
            r.update(x=float(xs), y=float(ys), th=th_deg); iyi.append(r)
    print(f"tam yörünge kontrolünü geçen: {len(iyi)}")
    if not iyi:
        sys.exit("uygun rota yok")
    iyi.sort(key=lambda r: -r["min_marj"])
    for r in iyi[:8]:
        print(f"  x={r['x']:.0f} y={r['y']:.0f} yön={r['th']}°  engel payı={r['min_marj']:.1f} m  "
              f"AGL {r['agl_min']:.1f}–{r['agl_max']:.1f} m")
    b = iyi[0]
    json.dump({"yol": b["pts"], "zemin": zemin,
               "secim": {k: b[k] for k in ("x", "y", "th", "min_marj", "agl_min", "agl_max", "n")},
               "parametre": vars(a), "adaylar": [{k: r[k] for k in ("x", "y", "th", "min_marj")} for r in iyi[:50]]},
              open(a.cikti, "w"))
    print("yazıldı:", a.cikti)


if __name__ == "__main__":
    main()
