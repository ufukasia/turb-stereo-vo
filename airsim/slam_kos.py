"""
AirSim sekansında stereo VO/SLAM koşucusu — config.yaml'A DOKUNMAZ.

config.yaml (şu an KITTI/Tabu ayarlarını taşıyor) bellekte kopyalanır, yalnız
AirSim'e özgü alanlar değiştirilir ve kod/ hattının kendi fonksiyonları çağrılır:
  run_frontend.main(cfg) -> gtsam_backend.main(args) -> eval_trajectory.evaluate

AirSim'e özgü değişiklikler (TASARIM.md §4.3):
  - türbülans KAPALI (temiz taban çizgisi)
  - max_init_depth >= 2 x irtifa (EuRoC'un 19 m'si 15-20 m irtifada bütün noktaları keserdi)
  - max_landmark_depth >= 4 x irtifa
  - döngü kapanışı varsayılan KAPALI (saf VO); --loop 1 ile açılır

    turbvio/bin/python slam_kos.py --veri <sdata_airsim> --seq NH_tarama_h15 --irtifa 15
"""
import argparse
import copy
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
KOD = os.path.join(os.path.dirname(HERE), "kod")
sys.path.insert(0, KOD)
os.chdir(KOD)

import config                                   # noqa: E402
import run_frontend                             # noqa: E402
from gtsam_backend import main as backend_main  # noqa: E402
from eval_trajectory import evaluate            # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--veri", required=True)
    ap.add_argument("--seq", required=True)
    ap.add_argument("--irtifa", type=float, default=15.0)
    ap.add_argument("--loop", type=int, default=0)
    ap.add_argument("--max_frames", type=int, default=0)
    ap.add_argument("--cikti", default=os.path.join(HERE, "cikti"))
    ap.add_argument("--extractor", default=None)
    ap.add_argument("--max_keypoints", type=int, default=None)
    ap.add_argument("--atla_frontend", action="store_true")
    ap.add_argument("--set", nargs="*", default=[], help="bolum.anahtar=deger (yaml olarak çözülür)")
    a = ap.parse_args()

    os.environ.pop("SLAM_DATA_ROOT", None)
    import yaml
    # config.load() dosyadaki sekansı (KITTI) diskte arayıp düşüyor; ham YAML okunur,
    # alanlar burada değiştirilir, doğrulama aşağıda AirSim sekansı için yapılır.
    cfg = yaml.safe_load(open(os.environ.get("SLAM_CONFIG", os.path.join(KOD, "config.yaml"))))
    cfg["dataset"]["base_path"] = a.veri
    cfg["dataset"]["sequence"] = a.seq
    cfg["dataset"]["max_frames"] = a.max_frames if a.max_frames > 0 else 10**7
    config._check_sequence(cfg)
    cfg["frontend"]["visualize"] = False
    cfg["frontend"]["save_visualizations"] = False
    if a.extractor:
        cfg["frontend"]["extractor"] = a.extractor
    if a.max_keypoints:
        cfg["frontend"]["max_keypoints"] = a.max_keypoints
    cfg["frontend"]["max_landmark_depth"] = max(cfg["frontend"]["max_landmark_depth"], 4 * a.irtifa)
    cfg["backend"]["max_init_depth"] = max(cfg["backend"]["max_init_depth"], 2 * a.irtifa)
    cfg["loop"]["enable"] = bool(a.loop)
    cfg["loop"]["pgo"] = bool(a.loop)
    cfg.setdefault("turbulence", {})["enabled"] = False
    for kv in a.set:
        k, v = kv.split("=", 1); b, key = k.split(".", 1)
        cfg[b][key] = yaml.safe_load(v)
        print(f"[geçersiz kılma] {b}.{key} = {cfg[b][key]}")
    os.makedirs(a.cikti, exist_ok=True)
    cfg["paths"]["outputs"] = a.cikti

    t0 = time.time()
    if not a.atla_frontend:
        st = run_frontend.main(cfg)
    else:
        st = {}
    t1 = time.time()
    out_csv = os.path.join(a.cikti, f"traj_{a.seq}{'_loop' if a.loop else ''}.csv")
    out_png = os.path.join(a.cikti, f"traj_{a.seq}{'_loop' if a.loop else ''}.png")
    bargs = config.backend_args(cfg, argv=["--out", out_csv])
    backend_main(bargs)
    t2 = time.time()
    ate, axis, n = evaluate(traj_csv=out_csv, gt_csv=config.gt_path(cfg), out_png=out_png,
                            dataset_name=f"AirSim {a.seq}")
    ozet = {"seq": a.seq, "ate_m": float(ate), "eksen_rmse": [float(v) for v in axis], "n": int(n),
            "frontend_s": t1 - t0, "backend_s": t2 - t1,
            "frontend_istatistik": {k: (float(v) if isinstance(v, (int, float)) else str(v))
                                    for k, v in (st or {}).items()},
            "max_init_depth": cfg["backend"]["max_init_depth"], "loop": a.loop}
    json.dump(ozet, open(os.path.join(a.cikti, f"ozet_{a.seq}{'_loop' if a.loop else ''}.json"), "w"),
              indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in ozet.items() if k != "frontend_istatistik"}, indent=2))


if __name__ == "__main__":
    main()
