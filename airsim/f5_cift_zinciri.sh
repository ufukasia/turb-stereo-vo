#!/bin/bash
cd "$(dirname "$0")"
until grep -q N12_BITTI <tmp>/n12_zinciri.log || ! pgrep -f "[n]12_zinciri.sh" >/dev/null; do sleep 30; done
../kod/turbvio/bin/python turb_izgara.py --sekanslar AF_F5_h5 --C0 0 1e-11 --tohumlar 0 1 2 --gurultu 1.0 \
   --pencereler 20:857 1157:1994 --gpu 1 --isci 10 --cikti cikti_izgara/f5_cift.jsonl
echo F5CIFT_BITTI
