"""3090 belirlenimlilik: cells_FAZ2 (cov_every=2, sigma=1) V2_01 D/r0 {0,4} tohum 0, A0/A0C, uzun-koşum yolu."""
import sys, os, dataclasses
ROOT = "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from turb.experiments import cells_FAZ2
from turb.runner import run_experiment_pipelined
if __name__ == "__main__":
    out, spool = sys.argv[1], sys.argv[2]
    cs = [dataclasses.replace(c, arms=("A0", "A0C")) for c in cells_FAZ2()
          if c.sequence == "V2_01_easy" and c.seed == 0 and c.D_over_r0 in (0.0, 4.0)]
    print([c.key() for c in cs])
    run_experiment_pipelined(cs, out, cpu_jobs=4, gpu_jobs=1, chunk=2, resume=False, spool=spool)
    print("BITTI")
