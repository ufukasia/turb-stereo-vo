"""Uçtan uca: experiments.cells_FAZ2 -> runner.run_experiment_pipelined (gerçek uzun koşum yolu), 2 hücre."""
import sys, os, dataclasses
ROOT = "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb"
sys.path.insert(0, ROOT); os.chdir(ROOT)
from turb.experiments import cells_FAZ2
from turb.runner import run_experiment_pipelined
if __name__ == "__main__":
    F = sys.argv[1]
    cs = [c for c in cells_FAZ2() if c.sequence == "V2_01_easy" and c.seed == 0 and c.D_over_r0 in (0.0, 4.0)]
    cs = [dataclasses.replace(c, arms=("A0", "A0C")) for c in cs]
    print([c.key() for c in cs])
    run_experiment_pipelined(cs, f"{F}/hat_testi.jsonl", cpu_jobs=4, gpu_jobs=1, chunk=2,
                             resume=False, spool=f"{F}/spool")
    print("BITTI")
