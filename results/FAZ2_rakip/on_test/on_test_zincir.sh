#!/bin/bash
# R3a ve duman üretimi bittikten sonra R2/R3b/R4/R5 (aynı anda tek rakip süreci)
cd "<project>/pyeco_turb"
D=evidence/FAZ2_rakip/on_test
# sekans beklemesi rakip_on_test.bekle() içinde (üretimle eşzamanlı)
echo "$(date +%T) basla yuk=$(cut -d' ' -f1-3 /proc/loadavg)"
"<project>/kod/turbvio/bin/python" -m turb.rakip_on_test --adim r4_temiz,r2a,r3b,r4,r5
echo "$(date +%T) zincir bitti yuk=$(cut -d' ' -f1-3 /proc/loadavg)"
