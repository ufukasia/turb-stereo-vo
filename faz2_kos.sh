#!/bin/bash
# Faz 2 EuRoC uzun koşumu — TEK KART RTX 3090 (ON_TEST_FAZ2_EUROC.md §6).
# Kullanım:  ./faz2_kos.sh kuru    [--deney FAZ2|FAKT|H6 --kirilma K]   (denetim + plan, koşmaz)
#            ./faz2_kos.sh baslat İŞÇİ [--deney ... --kirilma K]       (arka planda; devam ettirilebilir)
#            ./faz2_kos.sh denetle [--deney ...]                       (kayıt denetimi)
set -euo pipefail
cd "$(dirname "$0")"
export CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=4
export SLAM_DATA_ROOT="/data/sdata"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
PY="python"
mkdir -p evidence/FAZ2
case "${1:-kuru}" in
  kuru)    "$PY" -m turb.faz2_kos --kuru "${@:2}" ;;
  denetle) "$PY" -m turb.faz2_kos --denetle "${@:2}" ;;
  baslat)
    J="${2:-10}"
    "$PY" -m turb.faz2_kos --kuru "${@:3}" >/dev/null || { echo "kuru denetim KALDI"; exit 1; }
    nohup "$PY" -m turb.faz2_kos --jobs "$J" "${@:3}" > "evidence/FAZ2/run_$(date +%Y%m%d_%H%M%S)_${3:-}${4:-}.log" 2>&1 &
    echo "başladı PID $! (işçi $J); ilk hücreler bitince: ./faz2_kos.sh denetle" ;;
  *) echo "kullanım: $0 kuru|baslat [işçi]|denetle"; exit 1 ;;
esac
