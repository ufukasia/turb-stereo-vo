#!/bin/bash
# AirSim kolları (H3 + HA4) — başlatıcı. ON_TEST: africa/ON_TEST_KOLLAR.md; sapma taslağı: africa/SAPMA_TASLAK_KOLLAR.md
# Yalnız RTX 3090 (GPU1). Izgara arka-uç (CPU) koşar; GPU yalnız `kalibre`deki 3 yeniden üretimin ön-ucunda kullanılır.
#
#   ./kollar_kos.sh kuru                 # denetim + plan (koşmaz)
#   ./kollar_kos.sh kalibre [İŞÇİ]       # 1. adım: AirSim k seçimi (temiz pilot tohum 0-2), arka planda
#   ./kollar_kos.sh baslat  [İŞÇİ]       # 2. adım: ızgara (k_airsim.json gerekli), arka planda; ana ızgarayı bekler
#   ./kollar_kos.sh denetle              # kayıt denetimi (code_fp, A0 bit-özdeşliği, durum, NEES n)
#   ./kollar_kos.sh analiz               # H3 + HA4 + kollar -> cikti_kollar/izgara/analiz.json
# Durdurma: ilgili PID'i öldür (kayıtlar ana süreçte satır satır yazılır; yeniden başlatınca devam eder).
set -euo pipefail
cd "$(dirname "$0")"
export CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1 OMP_NUM_THREADS=4
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
unset SLAM_CONFIG
PY="python"
mkdir -p cikti_kollar/izgara cikti_kollar/kalibre
touch cikti_kollar/.trackerignore cikti_kollar/izgara/.trackerignore cikti_kollar/kalibre/.trackerignore
case "${1:-kuru}" in
  kuru)
    "$PY" kol_kos.py kalibre --kuru
    if [ -f cikti_kollar/k_airsim.json ]; then "$PY" kol_kos.py izgara --kuru
    else echo "[kuru] k_airsim.json yok -> önce: ./kollar_kos.sh kalibre"; fi ;;
  kalibre)
    J="${2:-3}"
    "$PY" kol_kos.py kalibre --kuru >/dev/null || { echo "kuru denetim KALDI"; exit 1; }
    nohup "$PY" kol_kos.py kalibre --isci "$J" > "cikti_kollar/kalibre/run_$(date +%Y%m%d_%H%M%S).log" 2>&1 &
    echo "kalibrasyon başladı PID $! (işçi $J); bitince: cat cikti_kollar/k_airsim.json" ;;
  baslat)
    J="${2:-3}"
    [ -f cikti_kollar/k_airsim.json ] || { echo "DUR: cikti_kollar/k_airsim.json yok (önce kalibre)"; exit 1; }
    "$PY" -c "import json,sys; d=json.load(open('cikti_kollar/k_airsim.json')); sys.exit(1 if d.get('sinirda') else 0)" \
      || { echo "DUR: k ızgara sınırında (k_airsim.json)"; exit 1; }
    "$PY" kol_kos.py izgara --kuru >/dev/null || { echo "kuru denetim KALDI"; exit 1; }
    nohup "$PY" kol_kos.py izgara --isci "$J" --bekle > "cikti_kollar/izgara/run_$(date +%Y%m%d_%H%M%S).log" 2>&1 &
    echo "ızgara başladı PID $! (işçi $J); ilk kayıtlar gelince: ./kollar_kos.sh denetle" ;;
  denetle)
    "$PY" - <<'PY'
import json, collections
R = [json.loads(l) for l in open("cikti_kollar/izgara/kayit.jsonl")] if __import__("os").path.exists("cikti_kollar/izgara/kayit.jsonl") else []
print(f"{len(R)} kayıt; durum {dict(collections.Counter((r.get('tur'), r.get('status')) for r in R))}")
K = [r for r in R if r.get("tur") == "kol" and r.get("status") == "ok"]
print("code_fp:", dict(collections.Counter(r.get("code_fp") for r in K)), "(beklenen 08e28b83407c)")
b = [r.get("a0_bit_ayni") for r in K if r["airsim"]["kol"] == "A0"]
print(f"A0 ana ızgarayla bayt aynı: {sum(1 for x in b if x)}/{len(b)}", "" if all(b) else "  !! DUR: A0 farklı")
nn = [r.get("n_n") for r in K if r.get("n_n")]
print(f"NEES örnek/koşum: min {min(nn) if nn else '-'}, medyan {sorted(nn)[len(nn)//2] if nn else '-'}")
for r in R:
    if r.get("status") != "ok":
        print("  !! başarısız:", r.get("airsim", {}).get("turb_sekans"), r.get("airsim", {}).get("kol"), r.get("error"))
PY
    ;;
  analiz) "$PY" analiz_kollar.py hepsi --etiket izgara --cikti cikti_kollar/izgara/analiz.json | tail -60 ;;
  *) echo "kullanım: $0 kuru|kalibre [işçi]|baslat [işçi]|denetle|analiz"; exit 1 ;;
esac
