"""DOGRULAMA SINAVI 2 — robust=0'da A1 ile A5 OZDES olmali mi?

Sinav 1 sunu gosterdi: Huber granulerligi A5'in temiz-veri ATE cezasinin
buyuk kismini acikliyor (0.1388 -> 0.0792), ama %29 artik fark kaldi.
Aciklama adayi: temel cizgi (gtsam_backend.py:451) HER ZAMAN Huber kullanir,
--turb_robust'tan bagimsiz. Yani robust=0'da A5'in aykiri-deger reddi YOK,
A0'in VAR. Artik fark bundansa:

  ON-KAYITLI TAHMIN: robust=0 ve D/r0=0'da  A1 ~ A5  (ikisi de Huber'siz,
  ikisinin de Sigma'si temel cizgi kosegeni) ve ikisi de A0'dan biraz kotu.

  A1 ile A5 birbirinden AYRILIRSA baska bir kusur var demektir ve E4 baslamadan
  once bulunmalidir.


Sifir turbulansta Sigma temel cizgiye indirgeniyor (olculdu: kosegen farki
4e-10, kosegen-disi TAM SIFIR). O halde A0 ile A5 arasindaki fark kovaryanstan
DEGIL, faktor yapisindan gelir. Aday: A0 gozlem basina Huber alir, A5 BLOK
basina (GTSAM m-kestiriciyi beyazlatilmis artigin NORMUNA uygular) -> yamadaki
5 gozlemden biri kotuyse digerleri de bastirilir.

Sinav: ayni TEMIZ hucre, A0 ve A5, robust=1 ve robust=0.
  granulerlik ise -> A5(robust=0) ~ A0  ve  A5(robust=1) belirgin kotu.

    ./turbvio/bin/python -m turb._huber_probe

NOT: multiprocessing 'spawn' __main__'i yeniden import eder; asil is
__main__ korumasinin ICINDE olmak ZORUNDA (ilk surumde degildi, cocuk surec
sonsuz kuruluma girdi).
"""
import json
import os


def main():
    from .runner import Cell, run_experiment_pipelined
    from .experiments import EUROC
    cells = []
    for rb in (0,):
        cells.append(Cell(dataset="euroc", base_path=EUROC, sequence="MH05",
                          seed=0, D_over_r0=0.0, arms=("A0", "A1", "A5"), nees=True,
                          cov_every=10, windows=(500,), window_len=400,
                          warmup=50, cache_dump=True, tag=f"huber{rb}",
                          backend_argv=("--turb_robust", str(rb))))
    out = "evidence/huber_probe2/records.jsonl"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    run_experiment_pipelined(cells, out, cpu_jobs=4, gpu_jobs=1, chunk=2,
                             resume=False, spool="/tmp/turb_huber_spool2")
    print("\n===== SONUC =====")
    for line in open(out):
        r = json.loads(line)
        tag = r["cell"]["tag"]
        if r.get("status") != "ok":
            print(f"  {tag:8s} {r['arm']:4s}  BASARISIZ"); continue
        print(f"  {tag:8s} {r['arm']:4s}  ATE={r['m_ate_se3']:.4f} m  "
              f"NEES={r.get('n_nees_pos_median', float('nan')):.1f}")


if __name__ == "__main__":
    main()
