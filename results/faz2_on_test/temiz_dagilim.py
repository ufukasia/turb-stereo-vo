"""Temiz dağılım ön-testi: MH02 + V1_01, D/r0 {0, 8}, tohum 0-4, sigma=1 DN, A0, uzun-koşum yolu."""
import sys, os
ROOT = "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from turb.runner import Cell, run_experiment_pipelined, windows_for
from turb.experiments import EUROC_FAZ2
if __name__ == "__main__":
    F = sys.argv[1]
    cells = []
    for s in ("V1_01_easy", "MH02"):
        w = windows_for(s, os.path.join(ROOT, "evidence/pencereler_faz2.json"))
        for dr in (0.0, 8.0):
            for sd in range(5):
                cells.append(Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence=s, seed=sd, D_over_r0=dr,
                                  arms=("A0",), windows=w, window_len=400, nees=True, cov_every=10,
                                  psf_legacy_scale=False, backend_argv=("--turb_cal_factor", "12.0"),
                                  tag="faz2_temizdag", cache_dump=False, sensor_noise_dn=1.0))
    print(len(cells), cells[0].key())
    run_experiment_pipelined(cells, f"{F}/temiz_dagilim.jsonl", cpu_jobs=6, gpu_jobs=1, chunk=10,
                             resume=True, spool=f"{F}/spool_td")
    print("BITTI")
