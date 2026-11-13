"""
Rakip ön-testleri R2 / R3b / R4 / R5 sürücüsü (ON_TEST_RAKIPLER.md §2). Aynı anda tek rakip süreci.

    python -m turb.rakip_on_test --adim r4_temiz,r2a,r3b,r4,r5

Kayıtlar: evidence/FAZ2_rakip/on_test/on_test.jsonl ; koşum dosyaları evidence/FAZ2_rakip/on_test/kosum/.
Sekanslar `evidence/FAZ2_rakip/on_test/uret_duman.sh` ile üretildi (/data/sdata_turb_rakip/on_test).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from turb import rakip_kos as RK  # noqa: E402

VERI = Path("/data/sdata")
OT = Path("/data/sdata_turb_rakip/on_test")
EV = ROOT / "evidence/FAZ2_rakip/on_test"
KAYIT = EV / "on_test.jsonl"
KOSUM = EV / "kosum"
SON = {"V1_01_easy": 1300, "MH02": 1800}


def ad(s, dr, sigma):
    n = f"_n{sigma:g}_ns0" if sigma > 0 else ""
    return OT / f"{s}__dr{dr:g}_tc0.15_t01_s0_speckle_global_sc0_sy{n}__k0-{SON[s]}"


def bekle(seq: Path):
    """Sekans üretimi bitene kadar bekle (üretici sonunda .yaziliyor -> ad yeniden adlandırır)."""
    import time
    while not (Path(seq) / "turb_meta.json").exists():
        time.sleep(15)


def kos(rk, seq, s, dr, sigma, protokol, tekrar=0, wins=None, etiket=None, kfp=None):
    if protokol != "pencere":
        bekle(seq)
    return RK.tek_kosum(rk, seq, VERI / s, s, dr, 0, sigma, protokol, KOSUM, KAYIT, tekrar=tekrar, wins=wins,
                        etiket=etiket or f"{s}_dr{dr:g}_n{sigma:g}_{protokol}", kimlik_fp=kfp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--adim", required=True)
    ap.add_argument("--rakipler", default=",".join(RK.RAKIP_LISTE))
    a = ap.parse_args()
    rks = a.rakipler.split(",")
    kfp = RK.kimlik_ozeti()
    for adim in a.adim.split(","):
        print(f"===== {adim}", flush=True)
        if adim == "r4_temiz":            # R2 (b) = R4'ün D/r0=0 hücreleri
            for s in ("V1_01_easy", "MH02"):
                for rk in rks:
                    kos(rk, ad(s, 0, 1), s, 0, 1, "onek", kfp=kfp)
        elif adim == "r2a":               # R2 (a): pencere başına kısa sekans (görünüm; aynı baytlar)
            for s in ("V1_01_easy", "MH02"):
                bekle(ad(s, 0, 1))
                for w in RK.pencereler(s):
                    g = RK.pencere_gorunumu(ad(s, 0, 1), EV / "gorunum" / f"{s}_w{w}", w)
                    for rk in rks:
                        kos(rk, g, s, 0, 1, "pencere", wins=[w], etiket=f"{s}_dr0_n1_pencere_w{w}", kfp=kfp)
        elif adim == "r3b":               # rektifiye temiz, gürültüsüz (V1_01 0:1300)
            for rk in rks:
                kos(rk, ad("V1_01_easy", 0, 0), "V1_01_easy", 0, 0, "onek_sigma0", kfp=kfp)
        elif adim == "r4":                # duman: D/r0 8 ve 24
            for dr in (8, 24):
                for s in ("V1_01_easy", "MH02"):
                    for rk in rks:
                        kos(rk, ad(s, dr, 1), s, dr, 1, "onek", kfp=kfp)
        elif adim == "r5":                # belirlenim: V1_01 D/r0=8 iki ek tekrar
            for t in (1, 2):
                for rk in rks:
                    kos(rk, ad("V1_01_easy", 8, 1), "V1_01_easy", 8, 1, "onek", tekrar=t, kfp=kfp)
        else:
            raise SystemExit(f"bilinmeyen adım {adim}")
    print("BITTI", flush=True)


if __name__ == "__main__":
    main()
