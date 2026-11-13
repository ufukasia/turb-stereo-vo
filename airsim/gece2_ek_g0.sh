#!/bin/bash
# F5 W2 ailesini GPU0'da tamamla (temiz ve 1e-11 tohumları GPU1'de koşmuştu; GPU'lar arası karşılaştırma yasak).
cd "$(dirname "$0")"
until grep -q GECE2_BITTI_0 <tmp>/gece2_g0.log || ! pgrep -f "[g]ece2_zinciri.sh 1157" >/dev/null; do sleep 30; done
../kod/turbvio/bin/python turb_izgara.py --sekanslar AF_F5_h5 --C0 0 1e-11 --tohumlar 0 1 2 --gurultu 1.0 --pencereler 1157:1994 --gpu 0 --isci 10 --cikti cikti_izgara/gece2_gpu0.jsonl
echo GECE2EK_BITTI
