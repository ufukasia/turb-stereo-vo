"""AirSim ana ızgarası (ON_KAYIT_FAZ2 v0.3 §4): 120 hücre, RTX 3090, K paralel işçi.

    python ana_kuyruk.py --isciler 4 --uret_isci 6

- Hücre sırası tohum-öncelikli: önce her (sekans, pencere) için bir tohumun bütün C0 düzeyleri biter -> erken tam çiftler.
- Devam ettirilebilir: cikti_izgara/ana_v03/w*.jsonl içinde durum=tamam olan hücreler atlanır.
- Durdurma: cikti_izgara/ana_v03/DUR dosyası oluşturulursa işçiler yeni hücre almaz.
"""
import argparse, json, subprocess, sys, threading, time
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "cikti_izgara" / "ana_v03"
PY = str(HERE.parent / "kod/turbvio/bin/python")
C0S = [0.0, 1e-13, 1e-12, 3e-12, 1e-11]
PENS = ["20:857", "1157:1994"]
TOHUM = {"AF_F5_h5": list(range(10, 15)), "AF_F15_h5": list(range(10, 17))}


def hucreler():
    out = []
    for i in range(max(len(v) for v in TOHUM.values())):
        for seq, ths in TOHUM.items():
            if i < len(ths):
                for pen in PENS:
                    for c in C0S:
                        out.append((seq, c, ths[i], pen))
    return out


def yapilanlar():
    d = set()
    for f in OUT.glob("w*.jsonl"):
        for l in open(f):
            r = json.loads(l)
            if r.get("durum") == "tamam":
                d.add((r["sekans"], float(r["C0"]), int(r["tohum"]), r["pencere"]))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--isciler", type=int, default=4)
    ap.add_argument("--uret_isci", type=int, default=6)
    ap.add_argument("--gpu", default="1")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tum = hucreler()
    assert len(tum) == 120 and len(set(tum)) == 120
    d = yapilanlar()
    kuyruk = [h for h in tum if h not in d]
    kilit = threading.Lock()
    log = open(OUT / "ilerleme.log", "a")

    def yaz(m):
        with kilit:
            log.write(time.strftime("%Y-%m-%d %H:%M:%S ") + m + "\n"); log.flush()

    yaz(f"başla: {len(tum)} hücre, {len(d)} yapılmış, {len(kuyruk)} kuyrukta; isciler={a.isciler} uret_isci={a.uret_isci} gpu={a.gpu}")

    def isci(k):
        while True:
            if (OUT / "DUR").exists():
                yaz(f"w{k}: DUR dosyası, çıkış"); return
            with kilit:
                if not kuyruk:
                    return
                seq, c, th, pen = kuyruk.pop(0)
            t0 = time.time()
            r = subprocess.run([PY, "turb_izgara.py", "--sekanslar", seq, "--C0", repr(c), "--tohumlar", str(th),
                                "--pencereler", pen, "--gurultu", "1.0", "--gpu", a.gpu, "--isci", str(a.uret_isci),
                                "--cikti", str(OUT / f"w{k}.jsonl")], cwd=HERE, capture_output=True, text=True)
            son = (r.stdout.strip().splitlines() or [""])[-1]
            yaz(f"w{k}: {seq} C0={c:g} th={th} pen={pen} rc={r.returncode} {time.time() - t0:.0f}s {son[:160]}")
            if r.returncode != 0:
                yaz(f"w{k}: HATA {r.stderr[-800:]}")

    th = [threading.Thread(target=isci, args=(k,)) for k in range(a.isciler)]
    for t in th:
        t.start()
    for t in th:
        t.join()
    yaz(f"ANA_V03_BITTI yapılmış={len(yapilanlar())}/120")


if __name__ == "__main__":
    main()
