#!/bin/bash
# Faz 2 rakip ızgarası (ORB-SLAM3 stereo, AirSLAM, DROID-SLAM) — ON_TEST_RAKIPLER.md §4.
# Tek rakip süreci, GPU yalnız RTX 3090; üretilen sekanslar hücre bitince silinir; devam ettirilebilir.
# Kullanım:  ./rakip_kuyruk.sh kuru          (denetim + plan, koşmaz)
#            ./rakip_kuyruk.sh baslat        (arka planda başlatır)
#            ./rakip_kuyruk.sh denetle       (kayıt denetimi)
#            ./rakip_kuyruk.sh analiz        (eşli oran özeti -> evidence/FAZ2_rakip/rakip_kayitlari_ozet.json)
set -euo pipefail
cd "$(dirname "$0")"
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export SLAM_DATA_ROOT="/data/sdata"
PY="python"
mkdir -p evidence/FAZ2_rakip
case "${1:-kuru}" in
  kuru)    "$PY" -m turb.rakip_kuyruk --kuru ;;
  denetle) "$PY" -m turb.rakip_kuyruk --denetle ;;
  analiz)  "$PY" -m turb.rakip_analiz ;;
  baslat)
    "$PY" -m turb.rakip_kuyruk --kuru >/dev/null || { echo "kuru denetim KALDI"; exit 1; }
    nohup "$PY" -m turb.rakip_kuyruk --sil > "evidence/FAZ2_rakip/kuyruk_$(date +%Y%m%d_%H%M%S).log" 2>&1 &
    echo "başladı PID $!; ilk hücreler bitince: ./rakip_kuyruk.sh denetle" ;;
  *) echo "kullanım: $0 kuru|baslat|denetle|analiz"; exit 1 ;;
esac
