#!/bin/bash
cd "../airsim_sim"
PY="python"
until grep -q "ZINCIR_BITTI" <tmp>/kayit_zinciri_F.log; do
  if ! pgrep -f "[k]ayit_zinciri_F.sh" >/dev/null; then echo "F ZINCIRI HATAYLA BITTI"; exit 1; fi
  sleep 20
done
setsid nohup ./kayit_zinciri_N.sh > <tmp>/kayit_zinciri_N.log 2>&1 < /dev/null &
"$PY" turb_izgara.py --sekanslar AF_F15_h5 AF_F5_h5 --C0 0 1e-12 1e-11 --tohumlar 0 --cikti cikti_izgara/pilot.jsonl --isci 14 --gpu 1
until grep -q "ZINCIR_N_BITTI" <tmp>/kayit_zinciri_N.log; do
  if ! pgrep -f "[k]ayit_zinciri_N.sh" >/dev/null; then echo "N ZINCIRI HATAYLA BITTI"; exit 1; fi
  sleep 20
done
"$PY" turb_izgara.py --sekanslar AF_N15_h30 --C0 0 1e-11 --tohumlar 0 --cikti cikti_izgara/pilot.jsonl --isci 14 --gpu 1
echo PILOT_BITTI
