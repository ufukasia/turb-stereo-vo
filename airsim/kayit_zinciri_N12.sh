#!/bin/bash
# Negatif kontrol: aynı rota, 15° kamera dik aşağı, 12 m AGL (30 m: tüm noktalar 40·B=16 m eşiğinin ötesinde, ölçek %3.2 kaydı); taze sim.
cd "../airsim_sim"
PY=$HOME/anaconda3/envs/colosseum/bin/python
ENV=/data/Africa_Savannah/Africa_Savannah/LinuxNoEditor/Africa_001.sh
$PY kayit.py --out x --seq x --fov 15 --pitch -90 --taban 0.40 --settings_yaz
./sim_baslat.sh $ENV 0 || exit 1
sleep 10
$PY kayit.py --out /data/sdata_airsim --seq AF_N15_h12 --senaryo dosya --rota africa/rota_F1.json \
   --irtifa 12 --irtifa_salinim 0.3 --hiz 3 --kose 6 --sure 0 --fov 15 --pitch -90 --taban 0.40 --dondur 1 \
   --temizle 'Water_Plane|Poacher|Crocodile|Hippopotamus|Rhinoceros' || exit 1
$PY dogrula.py /data/sdata_airsim/AF_N15_h12 5; ../kod/turbvio/bin/python dogrula_d5.py /data/sdata_airsim/AF_N15_h12 5
echo ZINCIR_N12_BITTI
