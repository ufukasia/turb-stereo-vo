"""Faz 2 ön-test sürücüsü (scratch). Ön-uç bir kez (GPU, sıralı), kollar paralel (CPU).

    python surucu.py plan.json out.jsonl [--cpu 6]

plan.json: [{"cell": {...Cell alanları...}, "arms": [["A0", []], ["A0C", ["--turb_cal_factor","12"]], ...],
             "fe_dir": "<kalıcı ön-uç dizini (varsa yeniden kullanılır)>"}, ...]
Her (hücre, kol) için runner._run_arms çağrılır -> kayıtlar runner ile AYNI biçimde.
"""
import json, os, sys, time, pickle, multiprocessing as mp, shutil, traceback
ROOT = "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb"


def _setup():
    sys.path.insert(0, ROOT); os.chdir(ROOT)


def frontend(job):
    _setup()
    from turb.runner import Cell, _build_cfg
    import config, run_frontend
    cell = Cell(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in job["cell"].items()})
    d = job["fe_dir"]
    os.makedirs(d, exist_ok=True)
    cfg = _build_cfg(cell, d)
    meta_p = os.path.join(d, "fe_meta.json")
    if os.path.exists(meta_p) and os.path.exists(config.dump_path(cfg)):
        return json.load(open(meta_p))
    import torch
    t0 = time.time()
    stats = run_frontend.main(cfg)
    meta = {"t_frontend": time.time() - t0, "gpu": torch.cuda.get_device_name(0),
            "stats": {k: float(v) for k, v in (stats or {}).items() if isinstance(v, (int, float))}}
    json.dump(meta, open(meta_p, "w"))
    return meta


def backend(args):
    job, arm, extra, tag = args
    _setup()
    try:
        import cv2; cv2.setNumThreads(1)
    except Exception:
        pass
    from turb.runner import Cell, _build_cfg, _run_arms
    import config
    c = dict(job["cell"]); c["arms"] = (arm,)
    c["backend_argv"] = tuple(c.get("backend_argv") or ()) + tuple(extra)
    cell = Cell(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in c.items()})
    d = job["fe_dir"]
    cfg = _build_cfg(cell, d)
    ad = os.path.join(d, "kol_" + tag)
    os.makedirs(ad, exist_ok=True)
    meta = json.load(open(os.path.join(d, "fe_meta.json")))
    qual = None
    qp = os.path.join(d, "turb_quality.json")
    if os.path.exists(qp):
        qual = json.load(open(qp))
    with open(config.dump_path(cfg), "rb") as f:
        f_ts = [int(fr["timestamp"]) for fr in pickle.load(f)["frames"]]
    base = {"t_frontend": meta["t_frontend"], "frontend_stats": meta["stats"], "quality": qual,
            "fe_gpu": meta.get("gpu"), "kol_etiket": tag, "n_dump_frames": len(f_ts)}
    log = open(os.path.join(ad, "be.log"), "w")
    so, se = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = log
    try:
        recs = _run_arms(cell, cfg, ad, config.gt_path(cfg), f_ts, base, verbose=True)
    finally:
        sys.stdout, sys.stderr = so, se
        log.close()
    return recs


def main():
    plan = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    cpu = int(sys.argv[sys.argv.index("--cpu") + 1]) if "--cpu" in sys.argv else 6
    ctx = mp.get_context("spawn")
    t0 = time.time()
    # Ön-uç: tek GPU işçisi, sıralı (bit-tekrarlanabilirlik: GPU paylaşılmaz)
    with ctx.Pool(1, maxtasksperchild=1) as gp, ctx.Pool(cpu, maxtasksperchild=4) as cp, \
            open(out, "a") as fo:
        pend = []
        for job in plan:
            m = gp.apply(frontend, (job,))
            print(f"[fe] {job['cell']['sequence']} dr{job['cell']['D_over_r0']} "
                  f"t={m['t_frontend']:.0f}s gpu={m.get('gpu')}  ({time.time()-t0:.0f}s)", flush=True)
            for arm, extra, *tg in job["arms"]:
                tag = tg[0] if tg else arm
                pend.append((job, cp.apply_async(backend, ((job, arm, extra, tag),))))
        for job, r in pend:
            try:
                recs = r.get()
            except Exception as e:
                recs = [{"status": "driver_fail", "error": repr(e), "cell": job["cell"]}]
            for x in recs:
                fo.write(json.dumps(x, default=lambda o: o.item() if hasattr(o, "item") else
                                    (list(o) if isinstance(o, (tuple, set)) else str(o))) + "\n")
            fo.flush()
            x = recs[0]
            print(f"[be] {x.get('key')} {x.get('kol_etiket', x.get('arm'))} "
                  f"t_be={x.get('t_backend', float('nan')):.0f}s "
                  f"ok={sum(1 for y in recs if y.get('status') == 'ok')}/{len(recs)}  "
                  f"({time.time()-t0:.0f}s)", flush=True)
    print("BITTI", flush=True)


if __name__ == "__main__":
    main()
