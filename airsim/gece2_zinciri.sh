#!/bin/bash
# Ek ön test 2: pencere başına bir GPU (aile içi karşılaştırmalar aynı GPU'da).
cd "$(dirname "$0")"
PY=../kod/turbvio/bin/python
W=$1; G=$2
$PY turb_izgara.py --sekanslar AF_F15_h5 --C0 0 1e-12 1e-11 --tohumlar 0 1 2 --gurultu 1.0 --pencereler $W --gpu $G --isci 10 --cikti cikti_izgara/gece2_gpu$G.jsonl
$PY turb_izgara.py --sekanslar AF_F5_h5 --C0 1e-12 3e-12 --tohumlar 0 1 2 --gurultu 1.0 --pencereler $W --gpu $G --isci 10 --cikti cikti_izgara/gece2_gpu$G.jsonl
echo GECE2_BITTI_$G
