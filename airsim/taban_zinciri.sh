#!/bin/bash
# Gürültü tabanı ön testi (africa/ON_TEST_GURULTU_TABANI.md). Hepsi GPU0'da.
cd "$(dirname "$0")"
PY=../kod/turbvio/bin/python
until grep -q KONTROL_BITTI <tmp>/kontrol_tekrar.log; do sleep 30; done
echo "[F5 pencereler]"
$PY turb_izgara.py --sekanslar AF_F5_h5 --C0 0 --pencereler 20:857 1157:1994 --gpu 0 --isci 10 --cikti cikti_izgara/taban.jsonl
echo "[F15 gurultu tabani]"
$PY turb_izgara.py --sekanslar AF_F15_h5 --C0 0 1e-11 --tohumlar 0 1 2 --gurultu 1.0 --gpu 0 --isci 10 --cikti cikti_izgara/taban.jsonl
echo TABAN_BITTI
