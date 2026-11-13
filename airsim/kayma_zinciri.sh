#!/bin/bash
# Ara-değerleme kontrolü (GPU1, temiz referans GPU1'de 0.108 m): sabit alt-piksel kayma 0.5 ve 0.1 px.
cd "$(dirname "$0")"
PY=../kod/turbvio/bin/python
until grep -q "PILOT_BITTI\|HATAYLA" <tmp>/pilot_duzenle.log || ! pgrep -f "[p]ilot_duzenle.sh" >/dev/null; do sleep 30; done
for k in 0.5 0.1; do
  $PY turb_izgara.py --sekanslar AF_F15_h5 --C0 0 --sabit_kayma $k --gpu 1 --isci 10 --cikti cikti_izgara/kayma.jsonl
done
echo KAYMA_BITTI
