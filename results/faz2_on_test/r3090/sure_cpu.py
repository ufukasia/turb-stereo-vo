"""Tam ızgara CPU-saati tahmini (3090 ön-testi; CPU-s = duvar süresi x örneklenen çekirdek oranı)."""
import json, numpy as np
W = json.load(open("/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb/evidence/pencereler_faz2_kapsam1m.json"))["sequences"]
# MH02 (1800 kare) ölçümleri: t_be [s] x çekirdek oranı (cpu.log: 1.18-1.35; kol başına)
MH = {0: {"A0": 1542*1.32, "A0C": 1503*1.27, "A1": 1450*1.29, "AG": 1037*1.35, "A0C+A3": 1502*0.97,
          "A3": 5798*1.19, "A5": 1345*1.26},
      8: {"A0": 403*1.25, "A0C": 345*1.25, "A1": 673*1.25, "AG": 626*1.25, "A0C+A3": 544*1.25,
          "A3": 3698*1.23, "A5": 4030*1.20}}
PER = {d: {a: v / 1800 for a, v in x.items()} for d, x in MH.items()}   # CPU-s / kare
V_OLCEK = np.median([750, 850, 466, 1139]) / 1300 / PER[0]["A0C"]        # V temiz A0C CPU/kare / MH02 temiz
YUKSEK = 1.5        # D/r0 >= 12: 4060 Ti'da A0 16/8 süre oranı ~1.9, A3 ~1.25-1.5 (varsayım, kaba)
FE = {0: 0.05, 8: 0.11}   # ön-uç duvar s/kare, en hafif yükte ölçülenler (3090: V2_01 D/r0=0 47 s/1300, D/r0=4 101 s/1300; MH02 D/r0=8 203 s/1800)
ARMS = ("A0", "A0C", "A1", "A3", "AG", "A0C+A3")
tot = a35 = fe = 0.0
for s, v in W.items():
    nf = v["windows"][-1]["end_excl"]; sc = 1.0 if s.startswith("MH") else V_OLCEK
    for dr in (0, 4, 8, 12, 16, 20, 24):
        base = PER[0] if dr == 0 else PER[8]
        f = 1.0 if dr <= 8 else YUKSEK
        fe += 5 * nf * (FE[0] + (FE[8] - FE[0]) * dr / 8)
        for sd in range(5):
            for a in ARMS + (("A5",) if s in ("MH02", "MH04") else ()):
                t = base[a] * nf * sc * f
                tot += t
                if a in ("A3", "A5"):
                    a35 += t
print(f"V/MH CPU ölçeği {V_OLCEK:.2f}")
print(f"arka-uç toplam {tot/3600:.0f} CPU-sa (A3+A5 payı {a35/tot:.0%}); ön-uç toplam duvar (tek akış) {fe/3600:.1f} sa")
print("| kullanılabilir CPU (donanım iş parçacığı) | arka-uç işçisi (≈1.3 iş parçacığı/işçi) | süre |")
for thr, w in ((24, 16), (16, 12), (12, 8)):
    print(f"| {thr} | {w} | ~{tot/3600/thr:.0f} sa |")
