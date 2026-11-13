#!/bin/bash
# Negatif kontrol: aynı rota, 15° kamera dik aşağı, 30 m AGL; taze sim.
cd "../airsim_sim"
PY=$HOME/anaconda3/envs/colosseum/bin/python
ENV=/data/Africa_Savannah/Africa_Savannah/LinuxNoEditor/Africa_001.sh
$PY kayit.py --out x --seq x --fov 15 --pitch -90 --taban 0.40 --settings_yaz
./sim_baslat.sh $ENV 0 || exit 1
sleep 10
$PY kayit.py --out /data/sdata_airsim --seq AF_N15_h30 --senaryo dosya --rota africa/rota_F1.json \
   --irtifa 30 --irtifa_salinim 0.3 --hiz 3 --kose 6 --sure 0 --fov 15 --pitch -90 --taban 0.40 --dondur 1 \
   --temizle 'Water_Plane|Poacher|Crocodile|Hippopotamus|Rhinoceros' || exit 1
$PY dogrula.py /data/sdata_airsim/AF_N15_h30 5
echo ZINCIR_N_BITTI
