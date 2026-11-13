#!/bin/bash
# (2026-11-08 13:10: --stereo_dis sensor_yaml ile yeniden; gömülü matrisli ilk üretim GEÇERSİZ, bkz. ON_TEST_RAKIPLER.md §4)
# Duman / R2 / R3b sekansları (önek: kare 0 -> son pencere sonu). Süre ve yük ortalaması kaydedilir.
cd "<project>/pyeco_turb"
PY="<project>/kod/turbvio/bin/python"
K="<data-disk>/sdata_turb_rakip/on_test"
export OMP_NUM_THREADS=2
u() {  # sekans son dr sigma
  echo "$(date +%T) yuk=$(cut -d' ' -f1-3 /proc/loadavg) BASLA $1 0:$2 dr$3 sigma$4"
  "$PY" -m turb.uret_sekans --girdi "<data-disk>/sdata/$1" --cikti_kok "$K" --D_over_r0 $3 --tohum 0 \
     --kareler 0:$2 --gurultu_dn $4 --isci 4 --stereo_dis sensor_yaml 2>&1 | grep uret_sekans
  echo "$(date +%T) yuk=$(cut -d' ' -f1-3 /proc/loadavg) BITTI"
}
u V1_01_easy 1300 0 1
u V1_01_easy 1300 0 0
u MH02 1800 0 1
u V1_01_easy 1300 8 1
u MH02 1800 8 1
u V1_01_easy 1300 24 1
u MH02 1800 24 1
echo HEPSI_BITTI
