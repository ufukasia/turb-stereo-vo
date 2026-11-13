#!/bin/bash
# F rejimi temiz kayıtları: aynı rota, FOV 15 ve 5; her biri taze açılmış simde.
cd "../airsim_sim"
PY=$HOME/anaconda3/envs/colosseum/bin/python
ENV=/data/Africa_Savannah/Africa_Savannah/LinuxNoEditor/Africa_001.sh
for FOV in 15 5; do
  $PY kayit.py --out x --seq x --fov $FOV --pitch -3 --taban 0.40 --settings_yaz
  ./sim_baslat.sh $ENV 0 || exit 1
  sleep 10
  $PY kayit.py --out /data/sdata_airsim --seq AF_F${FOV}_h5 --senaryo dosya --rota africa/rota_F1.json \
     --irtifa 5 --irtifa_salinim 0.3 --hiz 3 --kose 6 --sure 0 --fov $FOV --pitch -3 --taban 0.40 --dondur 1 \
     --temizle 'Water_Plane|Poacher|Crocodile|Hippopotamus|Rhinoceros' || exit 1
  $PY dogrula.py /data/sdata_airsim/AF_F${FOV}_h5 20
done
echo ZINCIR_BITTI
