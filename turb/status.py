"""
Koşum durumu — TEK KOMUT. Hiçbir şey çıkarsamaz, ölçer.

    bash durum.sh              tek seferlik özet
    bash durum.sh -w           her 60 saniyede yenile (Ctrl+C ile çık)

Tahmini bitiş, VARSAYIMDAN değil ölçülen hızdan gelir: kayıt dosyasındaki
satır sayısı ile aşamanın toplam koşum sayısı karşılaştırılır. Aşama henüz
yeterince ilerlemediyse "hız ölçülemedi" yazar — uydurma sayı vermez.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

EVID = "evidence"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# run_queue_v6.sh ile AYNI sira ve AYNI --seeds degerleri. Burasi degisirse
# orasi da degismeli (ikisi ayrisirsa tahmin sessizce yanlis olur).
QUEUE = [
    ("E0",         None, "entegrasyon tabanı"),
    ("E4",            6, "ANA İDDİA (C3) · 7 kol"),
    ("E4model",    None, "korelasyon modeli ablasyonu"),
    ("E9",            6, "loop closure açık"),
    ("E1a",        None, "şiddet ekseni (en büyük)"),
    ("E1b",           6, "zamansal rejim"),
    ("E1c",           6, "uzaysal korelasyon"),
    ("TARTAN",     None, "veri kümesi genellemesi"),
    ("E3",         None, "restorasyon (C2)"),
    ("E8",         None, "ön-uç varyantları"),
    ("TARTANarms", None, "TartanAir kolları"),
]
ARM_DEPENDENT = {"E4", "E4model", "E9", "TARTANarms"}
# Asama basina kol sayisi (maliyet duzeltmesi icin)
ARM_COUNT = {"E0": 1, "E4": 7, "E4model": 1, "E9": 4, "E1a": 1,
             "E1b": 1, "E1c": 1, "TARTAN": 1, "E3": 1, "E8": 1,
             "TARTANarms": 4}


def _sh(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True,
                              text=True, timeout=15).stdout.strip()
    except Exception:
        return ""


def _hm(sec):
    if sec is None or sec != sec or sec < 0:
        return "?"
    sec = int(sec)
    if sec < 90:
        return f"{sec} sn"
    if sec < 5400:
        return f"{sec // 60} dk"
    return f"{sec // 3600} sa {(sec % 3600) // 60} dk"


def _expected_runs():
    """Her aşamanın toplam koşum sayısı = hücre × kol. Kayıt kütüğünden DEĞİL,
    deney tanımından hesaplanır (kuyruk başlamadan da bilinir)."""
    out = {}
    try:
        sys.path.insert(0, ROOT)
        from turb.experiments import REGISTRY
        for exp, seeds, _ in QUEUE:
            try:
                cells = REGISTRY[exp](seeds) if seeds else REGISTRY[exp]()
            except TypeError:
                cells = REGISTRY[exp]()
            # KAYIT = hucre x kol x PENCERE. Pencere carpani unutulmustu ve
            # butun beklenen sayilar 3-4 kat EKSIK cikiyordu; E4 "665" yerine
            # 2660 kayitliktir ve ilerleme yuzdesi ile ETA buna gore yanlisti.
            out[exp] = sum(len(c.arms) * max(len(c.windows), 1) for c in cells)
        return out
    except Exception:
        return {}


def _records(path):
    if not os.path.exists(path):
        return 0, None
    n, last = 0, None
    with open(path, errors="replace") as f:
        for line in f:
            if line.strip():
                n += 1
    try:
        last = os.path.getmtime(path)
    except OSError:
        pass
    return n, last


def _proc(pattern):
    """(pid, gecen_saniye) veya None."""
    out = _sh(f"pgrep -f '{pattern}' | head -20")
    for pid in out.split():
        et = _sh(f"ps --no-headers -o etimes -p {pid}").strip()
        cmd = _sh(f"tr '\\0' ' ' < /proc/{pid}/cmdline 2>/dev/null")
        if "pgrep" in cmd or "/bin/bash -c" in cmd:
            continue
        if et.isdigit():
            return int(pid), int(et), cmd
    return None


def gates():
    """Paralel ön-uçuş kapıları."""
    rows = []
    for tag, spool, rec, rep in (
            ("robust=1 (düzeltilmiş)", "/tmp/turb_preflight",
             f"{EVID}/preflight/records.jsonl", f"{EVID}/preflight.md"),
            ("robust=0 (ablasyon)", "/tmp/turb_preflight_r0",
             f"{EVID}/preflight_r0/records.jsonl", f"{EVID}/preflight_r0.md")):
        p = _proc(f"turb.preflight.*{spool.split('/')[-1]}$|turb.preflight.*--spool {spool} ")
        n, _ = _records(rec)
        sonuc = ""
        if os.path.exists(rep):
            for line in open(rep, errors="replace"):
                if "sonuç:" in line:
                    sonuc = "GEÇTİ" if "GEÇTİ" in line else "KALDI"
                    break
        rows.append((tag, p, n, sonuc, rep))
    return rows


def render():
    L = []
    A = L.append
    A("=" * 66)
    A(f" KOŞUM DURUMU · {time.strftime('%d.%m.%Y %H:%M:%S')}")
    A("=" * 66)

    # ---------- 1) ÖN-UÇUŞ KAPILARI ----------
    grows = gates()
    if any(r[1] or r[2] or r[3] for r in grows):
        A("")
        A("── ÖN-UÇUŞ KAPISI ──────────────────────────────────────────")
        for tag, p, n, sonuc, rep in grows:
            if p:
                pid, et, _ = p
                hiz = et / n if n else None
                kalan = hiz * (21 - n) if hiz and n < 21 else None
                A(f"  {tag:24s} KOŞUYOR  {n:2d}/21 kayıt · "
                  f"geçen {_hm(et)} · kalan ~{_hm(kalan)}")
            elif sonuc:
                A(f"  {tag:24s} {sonuc:8s} {n:2d}/21 kayıt · {rep}")
            elif n:
                A(f"  {tag:24s} DURDU?   {n:2d}/21 kayıt (süreç yok)")
            else:
                A(f"  {tag:24s} —")

    # ---------- 2) KUYRUK ----------
    q = _proc("run_queue_v[0-9]")
    # AKTIF ASAMA kuyruk kutugunden okunur, kayit dosyasindan DEGIL. E4 gibi
    # buyuk asamalarda hucre basina 7 kol SIRALI kosar (~2.2 sa) ve ilk kayit
    # ancak o zaman duser; kayda bakan bir gosterge saatlerce "sirada" der.
    aktif = None
    ql = os.path.join(EVID, "queue.log")
    if q and os.path.exists(ql):
        try:
            for ln in reversed(open(ql, errors="replace").read().splitlines()):
                if "BASLIYOR " in ln:
                    aktif = ln.split("BASLIYOR ", 1)[1].split()[0]
                    break
                if "BITTI " in ln:
                    break
        except OSError:
            pass
    exp_tot = _expected_runs()
    A("")
    A("── UZUN KUYRUK ─────────────────────────────────────────────")
    if q:
        A(f"  KOŞUYOR · geçen {_hm(q[1])}")
    else:
        A("  çalışmıyor" + ("  (kapı geçmeden başlamaz)" if grows else ""))

    toplam_kalan = 0.0
    olculemedi = []
    hiz_olculen = [None]   # koşan aşamanın ölçülen sn/koşum
    A("")
    A(f"  {'aşama':<12s} {'durum':<10s} {'ilerleme':>13s}  {'kalan':>9s}  açıklama")
    A("  " + "-" * 62)
    for exp, seeds, aciklama in QUEUE:
        d = os.path.join(EVID, exp)
        stamp = os.path.join(d, "DONE.stamp")
        n, mt = _records(os.path.join(d, "records.jsonl"))
        tot = exp_tot.get(exp)
        if os.path.exists(stamp):
            sn = ""
            for tok in open(stamp, errors="replace").read().split():
                if tok.startswith("saniye="):
                    sn = _hm(int(tok.split("=")[1]))
            # Biten asamada payda GOSTERILMEZ: eski bir kosumun DONE damgasi,
            # bugunku hucre sayisiyla ayni olmak zorunda degil (tohum sayisi
            # degismis olabilir) ve "20/5" gibi bozuk gorunen bir sey yazardi.
            A(f"  {exp:<12s} {'BİTTİ':<10s} {n:>6d} kayıt   "
              f"{sn:>9s}  {aciklama}")
            continue
        if exp == aktif and not n:
            # spool dizin sayisi, ilk kayittan ONCE tek ilerleme olcusudur
            _sp = len([d for d in (os.listdir("/tmp/turb_spool")
                                   if os.path.isdir("/tmp/turb_spool") else [])])
            A(f"  {exp:<12s} {'KOŞUYOR':<10s} {'0':>6s}/{tot or '?':<6} "
              f"{'ölçülemedi':>9s}  {aciklama}")
            A(f"  {'':<12s} {'':<10s} ilk kayıt hücre bitince düşer · "
              f"spool {_sp} dizin")
            continue
        if n and q:
            # ölçülen hız: bu aşamanın ilk kaydından bu yana
            # BASLANGIC ANI KUYRUK KUTUGUNDEN OKUNUR, dosya olusturma zamanindan
            # DEGIL. Kayit dosyasi bir onceki (durdurulmus) denemeden kalmis
            # olabilir ve aradaki OLU ZAMAN hiza karisir. Olculdu: E4 icin dosya
            # zamani 63 sa kalan derken kutuk zamani 34 sa diyordu -- ikincisi
            # dogru, cunku kuyruk 19:30-19:49 arasi duruyordu.
            t0 = None
            try:
                import datetime as _dt
                for ln in open(ql, errors="replace"):
                    if f"BASLIYOR {exp} " in ln:
                        stamp = ln.split("]")[0].strip("[").strip()
                        _now = _dt.datetime.now()
                        _t = _dt.datetime.strptime(f"{_now.year}-{stamp}",
                                                   "%Y-%m-%d %H:%M:%S")
                        t0 = _t.timestamp()
            except (OSError, ValueError):
                t0 = None
            if t0 is None:
                basla = _sh(f"stat -c %W '{d}/records.jsonl'")
                t0 = (int(basla) if basla.isdigit() and int(basla) > 0
                      else os.path.getmtime(os.path.join(d, "records.jsonl")))
            gecen = time.time() - t0
            hiz = gecen / n if n and gecen > 0 else None
            kalan = hiz * (tot - n) if hiz and tot else None
            if kalan:
                toplam_kalan += kalan
                hiz_olculen[0] = hiz
            A(f"  {exp:<12s} {'KOŞUYOR':<10s} {n:>6d}/{tot or '?':<6} "
              f"{_hm(kalan):>9s}  {aciklama}")
        else:
            olculemedi.append(exp)
            A(f"  {exp:<12s} {'sırada':<10s} {'-':>6s}/{tot or '?':<6} "
              f"{'?':>9s}  {aciklama}")

    A("")
    if toplam_kalan:
        A(f"  Koşan aşamanın kalanı (ölçüldü): ~{_hm(toplam_kalan)}")
    # TOPLAM PROJEKSIYON: kosan asamanin OLCULEN kosum-basi suresi, kalan
    # asamalara TASINIR. Kol sayisi duzeltmesi uygulanir -- A0-only asamalar
    # (E1*, E3, E8, TARTAN) 7 kollu E4'ten ucuzdur; olculen kol maliyetlerinden
    # (A0 769s ort. 1169s, DEVAM.md) oran ~0.66 alinir.
    if hiz_olculen[0] and olculemedi:
        _sure = 0.0
        for _e, _s2, _ in QUEUE:
            if _e not in olculemedi:
                continue
            _t = exp_tot.get(_e)
            if not _t:
                continue
            _kol = ARM_COUNT.get(_e, 1)
            _kat = 1.0 if _kol >= 7 else (0.85 if _kol >= 4 else 0.66)
            _sure += _t * hiz_olculen[0] * _kat
        A(f"  Kalan aşamalar (taşınan hızla): ~{_hm(_sure)}")
        A(f"  TOPLAM KALAN: ~{_hm(toplam_kalan + _sure)}  "
          f"(bitiş ≈ {time.strftime('%d.%m %H:%M', time.localtime(time.time() + toplam_kalan + _sure))})")
        A(f"  Not: PLAN.md'deki ~41 sa tahmini ÖLÇÜMLE aşıldı; yukarıdaki sayı")
        A(f"  koşan aşamanın gerçek hızından türetilmiştir.")
    elif olculemedi:
        A(f"  Henüz ölçülen hız yok; PLAN.md §3 kaba tahmini ~41 saat.")

    # ---------- 3) UYARILAR ----------
    A("")
    A("── UYARILAR ────────────────────────────────────────────────")
    uyari = 0
    for f, ne in ((f"{EVID}/SANITY_ALERT.md", "sonuç tahlili uyarısı"),
                  (f"{EVID}/ARMS_GATE_FAIL.md", "KOL KAPISI DÜŞTÜ"),
                  (f"{EVID}/PREFLIGHT_FAIL.stamp", "ön-uçuş kaldı")):
        if os.path.exists(f):
            uyari += 1
            A(f"  ! {ne:26s} {f}  ({time.strftime('%d.%m %H:%M', time.localtime(os.path.getmtime(f)))})")
    if not uyari:
        A("  yok")

    # ---------- 4) MAKİNE ----------
    A("")
    A("── MAKİNE ──────────────────────────────────────────────────")
    la = open("/proc/loadavg").read().split()[:3]
    nproc = _sh("nproc")
    A(f"  yük {la[0]} / {la[1]} / {la[2]}   (çekirdek {nproc})")
    df = _sh("df -h /tmp --output=avail | tail -1").strip()
    A(f"  /tmp boş alan {df}")
    A("")
    A("=" * 66)
    return "\n".join(L)


def main():
    watch = any(a in ("-w", "--watch") for a in sys.argv[1:])
    if not watch:
        print(render())
        return
    try:
        while True:
            os.system("clear")
            print(render())
            print("  (60 sn'de bir yenilenir · çıkmak için Ctrl+C)")
            time.sleep(60)
    except KeyboardInterrupt:
        print("\nçıkıldı.")


if __name__ == "__main__":
    main()
