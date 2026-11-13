"""ON_KAYIT_FAZ2 §2 "kırılma noktası": A0'ın medyan ATE'sinin D/r0=0 medyanının 2 katını ilk aştığı D/r0,
ya da ıraksama oranının %20'yi ilk aştığı D/r0; hangisi önce gelirse.

    python -m turb.kirilma_hesap                      # evidence/FAZ2/records.jsonl
    python -m turb.kirilma_hesap --kismi              # eksik ızgarada da tablo bas (karar verilmez)

ATE: gerçek zamanlı yörünge, SE(3) (m_ate_se3), §2 birincil metrik. Birim: (sekans, tohum, pencere) A0 kaydı.
Iraksama: m_diverged ya da status != ok ya da m_valid değil. Geçersiz ATE medyana +inf olarak girer (tutucu).
Karar yalnız ızgara tamamsa (her D/r0'da 9 sekans × 5 tohum × bütün pencereler) verilir.
"""
import argparse, json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def hesapla(records_path, beklenen_sekans=9, beklenen_tohum=5):
    ate, irak, sek = defaultdict(list), defaultdict(list), defaultdict(set)
    for l in open(records_path):
        r = json.loads(l)
        if r.get("arm") != "A0":
            continue
        c = r.get("cell") or {}
        dr = float(c.get("D_over_r0"))
        ok = r.get("status") == "ok" and r.get("m_valid") and not r.get("m_diverged")
        irak[dr].append(0 if ok else 1)
        v = r.get("m_ate_se3")
        ate[dr].append(float(v) if (r.get("m_valid") and v is not None and math.isfinite(v)) else math.inf)
        sek[dr].add((c.get("sequence"), int(c.get("seed"))))
    tablo = []
    for dr in sorted(ate):
        tablo.append(dict(D_over_r0=dr, n=len(ate[dr]), medyan_ate=float(np.median(ate[dr])),
                          irak_oran=float(np.mean(irak[dr])), sekans_tohum=len(sek[dr])))
    tam = bool(tablo) and all(t["sekans_tohum"] == beklenen_sekans * beklenen_tohum for t in tablo) \
        and len(tablo) == 7
    taban = next((t["medyan_ate"] for t in tablo if t["D_over_r0"] == 0.0), None)
    k = None
    if taban is not None:
        for t in tablo:
            t["oran"] = t["medyan_ate"] / taban
            if t["D_over_r0"] > 0 and k is None and (t["medyan_ate"] > 2 * taban or t["irak_oran"] > 0.20):
                k = t["D_over_r0"]
                t["kirilma_nedeni"] = "ATE>2x" if t["medyan_ate"] > 2 * taban else "irak>%20"
    return dict(tablo=tablo, tam=tam, kirilma=k, taban_medyan=taban)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--girdi", default=str(ROOT / "evidence/FAZ2/records.jsonl"))
    ap.add_argument("--kismi", action="store_true")
    ap.add_argument("--cikti", default=str(ROOT / "evidence/FAZ2/kirilma.json"))
    a = ap.parse_args()
    s = hesapla(a.girdi)
    for t in s["tablo"]:
        print(f"D/r0={t['D_over_r0']:4.0f} n={t['n']:4d} sek×toh={t['sekans_tohum']:2d} medyan ATE={t['medyan_ate']:.4f} "
              f"oran={t.get('oran', float('nan')):.2f} ıraksama={t['irak_oran']:.2%} {t.get('kirilma_nedeni', '')}")
    print(f"ızgara tam: {s['tam']}; kırılma: {s['kirilma']}")
    if not s["tam"]:
        if a.kismi:
            return 0
        print("KARAR YOK: ızgara eksik"); return 2
    if s["kirilma"] is None:
        print("24'e kadar kırılma yok (§3 H1: 'test edilen aralıkta eşik yok')")
    json.dump(s, open(a.cikti, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
