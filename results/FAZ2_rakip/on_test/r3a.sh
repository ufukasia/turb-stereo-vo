#!/bin/bash
# R3a: ham V1_01_easy bütün sekans (özgün görüntü + özgün kalibrasyon), üç rakip sırayla (aynı anda tek rakip)
cd "<project>/pyeco_turb"
PY="<project>/kod/turbvio/bin/python"
S="<data-disk>/sdata/V1_01_easy"
for r in orbslam3 airslam droid; do
  "$PY" -m turb.rakip_kos kos --rakip $r --seq "$S" --kaynak "$S" --sekans V1_01_easy --dr 0 --tohum 0 --sigma 0 \
     --protokol butun_ham --out_dir evidence/FAZ2_rakip/on_test/kosum --kayit evidence/FAZ2_rakip/on_test/r3a.jsonl --etiket V1_01_ham_butun
done
echo BITTI
