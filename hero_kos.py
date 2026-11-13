#!/usr/bin/env python
"""Makale görseli için tek hücrenin (MH05, tohum 3, R5, D/r0 8) A0 ve KAPI koşumunu yörünge dosyalarını saklayarak
yineler. Saklı ön-uç dökümünü kullanır; kayıtlar doğrulayıcı ızgara kayıtlarıyla bit-aynı olmalıdır (denetlenir).
Kod parmak izi kümesinin dışındadır.   python hero_kos.py CIKTI_DIZINI DOKUM_DIZINI"""
import json, multiprocessing as mp, os, pickle, shutil, sys
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
SEQ, SEED, DR, REJ = "MH05", 3, 8.0, "R5"


def isci(a):
    v, cikti, dump = a
    os.chdir(ROOT)
    import config
    import kapi_dogrula as kd
    from turb.runner import _build_cfg, _run_arms
    cell = kd.hucre(v, SEQ, SEED, DR, REJ, "D1")
    d = os.path.join(cikti, f"spool_{v}")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    cfg = _build_cfg(cell, d)
    dp = config.dump_path(cfg)
    os.makedirs(os.path.dirname(dp), exist_ok=True)
    os.symlink(os.path.abspath(dump), dp)
    with open(dp, "rb") as f:
        f_ts = [int(fr["timestamp"]) for fr in pickle.load(f)["frames"]]
    recs = _run_arms(cell, cfg, d, config.gt_path(cfg), f_ts, {})
    kol = cell.arms[0]
    for ad in (f"traj_{kol}.csv", f"traj_{kol}_realtime.csv"):
        shutil.copy(os.path.join(d, ad), os.path.join(cikti, ad.replace(kol, v)))
    shutil.copy(config.gt_path(cfg), os.path.join(cikti, "gt.csv"))
    json.dump(dict(f_ts=f_ts), open(os.path.join(cikti, "kare_zamanlari.json"), "w"))
    return v, recs


if __name__ == "__main__":
    cikti, dizin = os.path.abspath(sys.argv[1]), sys.argv[2]
    os.makedirs(cikti, exist_ok=True)
    import kapi_dogrula as kd
    dump = kd.sakli_dokum(REJ, SEQ, SEED, DR, [dizin])
    assert dump, "döküm yok"
    with mp.get_context("spawn").Pool(2) as p:
        out = p.map(isci, [("A0", cikti, dump), ("KAPI", cikti, dump)])
    json.dump({v: r for v, r in out}, open(os.path.join(cikti, "kayitlar.json"), "w"), default=float)
    for v, recs in out:
        for r in recs:
            print(v, r["window"], r.get("status"), r.get("m_ate_se3"), r.get("ms_ate_se3"))
