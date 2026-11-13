#!/usr/bin/env python
"""Sayısal ikiz denetimi N′ (ON_KAYIT_FAZ2.md Sapma 3 eki, madde b).

Her N hücresinin (onarım yok) ön-uç dökümü saklanır; dökümde YALNIZ `rect_K`'nın fx ve fy değerleri (1 + EPS) ile
çarpılır, A0 ve A1 yeniden koşulur. Amaç: arka-ucun 1e-12 düzeyindeki girdi farkına verdiği yanıtı (sayısal gürültü
tabanı) O / N ve C / N oranlarının yanına koymak. Parmak izli kodda değişiklik yok; bu dosya parmak izi kümesinin
dışındadır.

    python ikiz_kos.py sakla --dokumler D --pid PID          # ızgara (PID) sürerken N dökümlerini sabit bağlantıyla sakla
    python ikiz_kos.py kos   --dokumler D [--jobs 3] [--pid PID] [--yalniz V1_01_easy] [--eps 1e-12]
    python ikiz_kos.py plan  --dokumler D

`kos` bir N hücresini ancak ızgaradaki kayıtları tamamlanınca işler (döküm kesinleşmiştir; ızgaranın arka-uç işçisi
spool dizinini siler, sabit bağlantı kalır). --pid verilirse ızgara süreci yaşadıkça yeni dökümleri bekler.
--eps 0 denetim koşumudur (ızgaranın N kayıtlarıyla bit-aynı olmalı); kayıtları ikiz_denetim.jsonl'a yazılır.
Çıktı: evidence/ONARIM/ikiz.jsonl (kayıt biçimi ızgarayla aynı; ek alanlar `onarim.mod = "Ni"` ve `ikiz`).
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import multiprocessing as mp
import os
import pickle
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

import onarim_kos  # noqa: E402
from turb.runner import _build_cfg, _backend_worker, _js, code_fingerprint  # noqa: E402

EPS = 1e-12
KOLLAR = ("A0", "A1")
SPOOL_IZGARA = "/tmp/turb_spool_onarim"
SPOOL_IKIZ = "/tmp/turb_spool_onarim_ikiz"
IZGARA = os.path.join(onarim_kos.EV, "records.jsonl")
BIT_ALANLARI = ("m_ate_se3", "n_nees_pos_mean", "n_nees_rot_mean", "m_rpe_trans")


def n_plani(seqs):
    return [p for p in onarim_kos.plan(onarim_kos.KOSULLAR, seqs) if p[3] == "N"]


def dokum_yolu(D, p):
    return os.path.join(D, onarim_kos.ad_mod(*p) + ".pkl")


def ikiz_hucre(p):
    c = onarim_kos.hucre(*p)
    return dataclasses.replace(c, arms=KOLLAR, tag=c.tag + "_ikiz")


def yasiyor(pid):
    if not pid:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def tamam(c, done):
    return all((c.key(), int(w), x) in done for w in c.windows for x in c.arms)


# ----------------------------------------------------------------------------- döküm saklama
def sakla_bir_tur(D, P):
    n = 0
    for p in P:
        c = onarim_kos.hucre(*p)
        src = os.path.join(SPOOL_IZGARA, c.key().replace("/", "_"), f"landmarks_{c.sequence}.pkl")
        dst = dokum_yolu(D, p)
        try:
            st = os.stat(src)
        except FileNotFoundError:
            continue
        if os.path.exists(dst) and os.stat(dst).st_ino == st.st_ino:
            continue
        tmp = dst + ".yeni"
        if os.path.exists(tmp):
            os.remove(tmp)
        try:
            os.link(src, tmp)               # aynı inode: yazım sürse de bağlantı son içeriği gösterir
        except FileNotFoundError:
            continue
        os.replace(tmp, dst)                # ön-uç yeniden denerse (yeni inode) bağlantı yenilenir
        n += 1
        print(f"[sakla] {time.strftime('%T')} {os.path.basename(dst)} (inode {st.st_ino})", flush=True)
    return n


def sakla(D, P, pid, adim=3.0):
    os.makedirs(D, exist_ok=True)
    while True:
        son = not yasiyor(pid)
        sakla_bir_tur(D, P)
        if son:
            break
        time.sleep(adim)
    var = sum(os.path.exists(dokum_yolu(D, p)) for p in P)
    print(f"[sakla] bitti: {var}/{len(P)} döküm saklı")


# ----------------------------------------------------------------------------- ikiz koşumu
def _ikiz_isci(args):
    p, yol, eps = args
    import config
    cell = ikiz_hucre(p)
    d = os.path.join(SPOOL_IKIZ, cell.key().replace("/", "_") + (f"_eps{eps:g}"))
    os.makedirs(d, exist_ok=True)
    with open(yol, "rb") as f:
        dk = pickle.load(f)
    fx, fy, cx, cy = (float(x) for x in dk["rect_K"])
    yeni = [fx * (1.0 + eps), fy * (1.0 + eps), cx, cy]
    if eps != 0.0:
        assert yeni[0] != fx and yeni[1] != fy, "pertürbasyon çift duyarlıkta kayboldu"
        dk["rect_K"] = yeni                 # T-E ile aynı biçim (düz liste); eps = 0'da döküm değişmez
    cfg = _build_cfg(cell, d)
    dp = config.dump_path(cfg)
    with open(dp, "wb") as f:
        pickle.dump(dk, f)
    del dk
    meta = dict(cell=cell, ok=True, dir=d, dump=dp, gt=config.gt_path(cfg), t_frontend=None, stats={}, gpu=None)
    t0 = time.time()
    recs = _backend_worker(meta)
    for r in recs:
        r["onarim"] = {"sekans": p[0], "rejim": p[1], "dr": p[2], "mod": "Ni"}
        r["ikiz"] = {"eps": eps, "fx": [fx, yeni[0]], "fy": [fy, yeni[1]], "sure_s": round(time.time() - t0, 1),
                     "dokum": os.path.basename(yol), "dokum_bayt": os.path.getsize(yol)}
    return recs


def izgara_kayitlari():
    """-> {(sekans adı, kol, pencere): kayıt} (ızgaranın N hücreleri)"""
    out = {}
    if os.path.exists(IZGARA):
        with open(IZGARA) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if "arm" in r and (r.get("onarim") or {}).get("mod") == "N":
                    out[(r["cell"]["sequence"], r["arm"], int(r["window"]))] = r
    return out


def kos(D, P, out, jobs, pid, eps):
    fp = code_fingerprint()
    print(f"[ikiz] code_fp {fp} (beklenen {onarim_kos.BEKLENEN_FP}); eps {eps:g}; kollar {KOLLAR}; çıktı {out}", flush=True)
    if fp != onarim_kos.BEKLENEN_FP:
        sys.exit("[ikiz] DUR: kod parmak izi ana ızgaradan farklı")
    os.makedirs(SPOOL_IKIZ, exist_ok=True)
    t_bas, gonderilen, bekleyen, sayac = time.time(), set(), [], {"n": 0}
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=jobs, maxtasksperchild=1) as cp, open(out, "a") as fout:
        def yaz(recs):
            ref = izgara_kayitlari()
            for r in recs:
                if "arm" in r:
                    m = ref.get((r["cell"]["sequence"], r["arm"], int(r["window"]))) or {}
                    r["ikiz"]["n_ile_bit_ayni"] = all(r.get(k) == m.get(k) for k in BIT_ALANLARI) if m else None
                fout.write(json.dumps(r, default=_js) + "\n")
            fout.flush()
            sayac["n"] += 1
            ozet = [(r.get("arm"), r.get("window"), r.get("status"), (r.get("ikiz") or {}).get("n_ile_bit_ayni"))
                    for r in recs]
            print(f"[ikiz] {sayac['n']} hücre yazıldı, {(time.time() - t_bas) / 60:.1f} dk: "
                  f"{recs[0].get('key') if recs else ''} {ozet if eps == 0.0 else [o[:3] for o in ozet]}", flush=True)

        def hata(e):
            print(f"[ikiz] İŞ HATASI: {type(e).__name__}: {e}", flush=True)

        while True:
            canli = yasiyor(pid)            # kayıtları okumadan ÖNCE: ızgara arada biterse son tur yine de okur
            izgara_done = onarim_kos.okunan(IZGARA)
            ikiz_done = onarim_kos.okunan(out)
            kalan = [p for p in P if ikiz_hucre(p).key() not in gonderilen and not tamam(ikiz_hucre(p), ikiz_done)]
            simdi = [p for p in kalan if tamam(onarim_kos.hucre(*p), izgara_done) and os.path.exists(dokum_yolu(D, p))]
            for p in simdi:
                gonderilen.add(ikiz_hucre(p).key())
                print(f"[ikiz] gönderildi: {onarim_kos.ad_mod(*p)}", flush=True)
                bekleyen.append(cp.apply_async(_ikiz_isci, ((p, dokum_yolu(D, p), eps),), callback=yaz,
                                               error_callback=hata))
            kalan = [p for p in kalan if p not in simdi]
            if not kalan or not canli:
                if kalan:
                    print(f"[ikiz] ızgara süreci yok; işlenemeyen {len(kalan)} hücre: "
                          f"{[onarim_kos.ad_mod(*p) for p in kalan]}", flush=True)
                break
            time.sleep(60)
        for b in bekleyen:
            b.wait()
    print(f"[ikiz] bitti, {(time.time() - t_bas) / 60:.1f} dk")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("komut", choices=["plan", "sakla", "kos"])
    ap.add_argument("--dokumler", required=True)
    ap.add_argument("--pid", type=int, default=0, help="onarım ızgarası süreci")
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--yalniz", default="")
    ap.add_argument("--eps", type=float, default=EPS)
    a = ap.parse_args()
    seqs = [s for s in a.yalniz.split(",") if s] or onarim_kos.SEQS
    P = n_plani(seqs)
    out = os.path.join(onarim_kos.EV, "ikiz.jsonl" if a.eps != 0.0 else "ikiz_denetim.jsonl")
    if a.komut == "plan":
        izgara_done, ikiz_done = onarim_kos.okunan(IZGARA), onarim_kos.okunan(out)
        for p in P:
            c, ci = onarim_kos.hucre(*p), ikiz_hucre(p)
            print(f"  {onarim_kos.ad_mod(*p):28s} döküm={int(os.path.exists(dokum_yolu(a.dokumler, p)))}  "
                  f"ızgara N tamam={int(tamam(c, izgara_done))}  ikiz tamam={int(tamam(ci, ikiz_done))}  {ci.key()}")
        print(f"[plan] {len(P)} N hücresi, kollar {KOLLAR}, eps {a.eps:g}; çıktı {out}")
    elif a.komut == "sakla":
        sakla(a.dokumler, P, a.pid)
    else:
        kos(a.dokumler, P, out, a.jobs, a.pid, a.eps)


if __name__ == "__main__":
    main()
