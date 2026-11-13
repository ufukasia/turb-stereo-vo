#!/bin/bash
# PyEco (pyeco_turb, türbülans kapalı) koşucusu — pyeco_turb/config.yaml'a DOKUNMAZ.
# Kullanım: pyeco_kos.sh <veri_kök> <sekans> <çıktı_dizini> [bolum.anahtar=deger ...]
# Çıktı dizinine türetilmiş config yazılır (SLAM_CONFIG), ön-uç + arka-uç koşulur.
set -e
VERI="$(realpath -m "$1")"; SEQ="$2"; OUT="$(realpath -m "$3")"; shift 3
PE="../pyeco_turb"
P="python"
mkdir -p "$OUT"
"$P" - "$PE/config.yaml" "$OUT/config.yaml" "$VERI" "$SEQ" "$OUT" "$@" <<'PY'
import sys, yaml
src, dst, veri, seq, out, *kv = sys.argv[1:]
c = yaml.safe_load(open(src))
c["dataset"]["base_path"] = veri.rstrip("/") + "/"; c["dataset"]["sequence"] = seq
c["paths"]["outputs"] = out
c.setdefault("turbulence", {})["enabled"] = False
for s in kv:
    k, v = s.split("=", 1); b, key = k.split(".", 1)
    assert key in c[b], f"bilinmeyen anahtar {k}"
    c[b][key] = yaml.safe_load(v); print(f"[geçersiz kılma] {k} = {c[b][key]}")
yaml.safe_dump(c, open(dst, "w"), sort_keys=False, allow_unicode=True)
PY
cd "$PE"
export SLAM_CONFIG="$OUT/config.yaml" CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=${GPU:-1} OMP_NUM_THREADS=4
# DUMP_KAYNAK: başka bir koşumun ön-uç dökümü (kollar ön-ucu paylaşır; D/r0=0'da kanca kimliktir)
if [ -n "$DUMP_KAYNAK" ]; then
  ln -sf "$DUMP_KAYNAK" "$OUT/landmarks_${SEQ}.pkl"; echo "on-uc atlandi: $DUMP_KAYNAK" > "$OUT/fe.log"
else
  "$P" -u run_frontend.py > "$OUT/fe.log" 2>&1
fi
# COV_OUT: pencere-göreli NEES için kovaryans dökümü (COV_EVERY ile seyreltme)
if [ -n "$COV_OUT" ]; then
  "$P" -u run_backend.py --cov_out "$COV_OUT" --cov_every "${COV_EVERY:-1}" > "$OUT/be.log" 2>&1
else
  "$P" -u run_backend.py > "$OUT/be.log" 2>&1
fi
grep -E "ATE|ölçek|scale" "$OUT/be.log" | tail -4
