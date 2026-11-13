"""
E4/E5'in ÖNBELLEKTEN karşılanan hücrelerini hemen koş (C3'ün erken cevabı).

GEREKÇE: E4'ün ön-uç dump'ları E1a tarafından üretilip kalıcı önbelleğe yazılıyor
(turb/runner.py `_fe_key`). Önbellekte hazır olan hücreler için E4'ün GPU fazı
YOKTUR — yalnızca CPU arka-uç kollarıdır. Bu yüzden E1a hâlâ koşarken, GPU'yu hiç
rahatsız etmeden E4'ün hazır kısmı başlatılabilir ve yöntemin (C3) istatistiksel
cevabı ~10 saat erken alınır. Kalan hücreler kuyruk sırası gelince `--resume` ile
tamamlanır; hiçbir iş tekrar edilmez.

    ./turbvio/bin/python -m turb.run_e4_cached --jobs 6
"""
from __future__ import annotations

import argparse
import os

from .experiments import cells_E4
from .runner import run_experiment_pipelined, _fe_key, _cache_dir


def main():
    ap = argparse.ArgumentParser(description="E4'ün önbellekli hücrelerini koş")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--seeds", type=int, default=6)
    ap.add_argument("--out", default="evidence/E4/records.jsonl")
    ap.add_argument("--chunk", type=int, default=12)
    a = ap.parse_args()

    cells = cells_E4(a.seeds)
    cache = _cache_dir()
    ready = [c for c in cells
             if os.path.exists(os.path.join(cache, _fe_key(c) + ".pkl"))]
    print(f"[E4-cached] {len(ready)}/{len(cells)} hücre önbellekte hazır "
          f"(GPU fazı yok) -> {a.out}", flush=True)
    if not ready:
        print("[E4-cached] önbellek henüz boş; E1a ilerleyince tekrar deneyin.")
        return
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    # gpu_jobs=1: önbellek isabet etmezse tek bir ön-uç koşar (E1a'yı boğmaz).
    run_experiment_pipelined(ready, a.out, cpu_jobs=a.jobs, gpu_jobs=1,
                             chunk=a.chunk, resume=True)


if __name__ == "__main__":
    main()
