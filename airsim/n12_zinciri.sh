#!/bin/bash
cd "$(dirname "$0")"
until grep -q KAYMA_BITTI <tmp>/kayma_zinciri.log || ! pgrep -f "[k]ayma_zinciri.sh" >/dev/null; do sleep 30; done
../kod/turbvio/bin/python turb_izgara.py --sekanslar AF_N15_h12 --C0 0 --gpu 1 --isci 10 --cikti cikti_izgara/n12.jsonl
echo N12_BITTI
