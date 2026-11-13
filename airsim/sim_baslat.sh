#!/bin/bash
# Kullanim: sim_baslat.sh <ortam .sh yolu> [gpu_index]
# ~/sim_envs altindaki calisan butun UE ortamlarini durdurur, yenisini goruntusuz baslatir.
SCRIPT="$1"; GPU="${2:-1}"
for p in $(pgrep -f "$HOME/sim_envs/.*/Binaries/Linux/" ); do kill $p 2>/dev/null; done
sleep 4
cd "$(dirname "$SCRIPT")"
setsid nohup "$SCRIPT" -RenderOffscreen -nosound -graphicsadapter=$GPU > <tmp>/sim_env.log 2>&1 < /dev/null &
for i in $(seq 1 36); do ss -ltn | grep -q 41451 && { echo "API hazir"; exit 0; }; sleep 5; done
echo "API acilmadi"; tail -20 <tmp>/sim_env.log; exit 1
