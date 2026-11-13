#!/bin/bash
# Gürültü tabanı kontrolleri (GPU0): PyEco tekrar-koşum belirlenimi ve "boş pertürbasyon" etkisi.
# 1) temiz F15 tam tekrar  2) temiz F15 yalnız arka-uç tekrarı (aynı ön-uç dökümü)
# 3) C0=1e-16 (hat çalışır, kayma ~1e-4 px)  4) C0=1e-12 tohum 1
cd "$(dirname "$0")"
PY=../kod/turbvio/bin/python
V=/data/sdata_airsim
K=cikti_kontrol; mkdir -p $K
export GPU=0
echo "[1] temiz tam tekrar"; ./pyeco_kos.sh $V AF_F15_h5 $K/F15_temiz_rep2
echo "[2] yalniz arka-uc tekrar"; DUMP_KAYNAK=$PWD/cikti_izgara/AF_F15_h5/landmarks_AF_F15_h5.pkl ./pyeco_kos.sh $V AF_F15_h5 $K/F15_temiz_be2
echo "[3-4] izgara kontrolleri"; $PY turb_izgara.py --sekanslar AF_F15_h5 --C0 1e-16 --tohumlar 0 --gpu 0 --isci 8 --cikti $K/kontrol.jsonl
$PY turb_izgara.py --sekanslar AF_F15_h5 --C0 1e-12 --tohumlar 1 --gpu 0 --isci 8 --cikti $K/kontrol.jsonl
echo KONTROL_BITTI
