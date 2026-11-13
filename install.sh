#!/usr/bin/env bash
# ==========================================================================================
#  PyEco-SLAM — one-command install (Ubuntu 20.04 / 22.04 / 24.04, x86_64)
#
#      bash install.sh            # CUDA 12.4 wheels (NVIDIA driver >= 550 recommended)
#      bash install.sh --cpu      # machine without a GPU (for trying things out; ALIKED/LightGlue are slow on CPU)
#
#  What it does:
#    1. System libraries (libGL/libglib for the OpenCV window) — via apt if sudo is available.
#    2. Installs Miniforge (conda-forge), so Python 3.12 comes independently of the distribution;
#       Ubuntu 20.04's 3.8, 22.04's 3.10 or 24.04's 3.12 make no difference.
#    3. Creates the `pyeco` environment: requirements.txt + torch/torchvision (CUDA or CPU index).
#    4. Pins numpy to 1.26.4 LAST (the gtsam wheel segfaults with numpy 2.x).
#    5. Verifies: numpy / cv2 / torch / gtsam / kornia import and CUDA visibility.
#
#  Environment variables: MINIFORGE_DIR (default ~/miniforge3), ENV_NAME (default pyeco).
#  Safe to re-run (idempotent). Reference machine: Ubuntu 24.04, RTX 3090,
#  Python 3.12.2 (conda-forge), torch 2.5.1+cu124, gtsam 4.2.1.
# ==========================================================================================
set -euo pipefail

CPU_ONLY=0
for a in "$@"; do
  case "$a" in
    --cpu) CPU_ONLY=1 ;;
    -h|--help) sed -n 2,20p "$0"; exit 0 ;;
    *) echo "unknown option: $a"; exit 2 ;;
  esac
done

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MINIFORGE_DIR="${MINIFORGE_DIR:-$HOME/miniforge3}"
ENV_NAME="${ENV_NAME:-pyeco}"
PY_VER="3.12"
TORCH_INDEX="https://download.pytorch.org/whl/cu124"
[ "$CPU_ONLY" = 1 ] && TORCH_INDEX="https://download.pytorch.org/whl/cpu"

say()  { printf '\n\033[1;34m[PyEco] %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33m[PyEco] WARNING: %s\033[0m\n' "$*"; }

# ---------------------------------------------------------------- 0) prerequisites
if [ "$(uname -s)" != "Linux" ] || [ "$(uname -m)" != "x86_64" ]; then
  echo "This script targets Linux x86_64 (found: $(uname -s) $(uname -m))."; exit 1
fi
if [ -r /etc/os-release ]; then
  . /etc/os-release
  say "Distribution: ${PRETTY_NAME:-unknown}"
  case "${VERSION_ID:-}" in
    20.04|22.04|24.04) ;;
    *) warn "Targeted versions are Ubuntu 20.04/22.04/24.04; continuing anyway." ;;
  esac
fi

# ---------------------------------------------------------------- 1) system libraries
NEED_PKGS=()
ldconfig -p | grep -q 'libGL.so.1'        || NEED_PKGS+=(libgl1)
ldconfig -p | grep -q 'libglib-2.0.so.0'  || NEED_PKGS+=(libglib2.0-0)
ldconfig -p | grep -q 'libgomp.so.1'      || NEED_PKGS+=(libgomp1)
command -v curl >/dev/null || NEED_PKGS+=(curl)
command -v git  >/dev/null || NEED_PKGS+=(git)
if [ "${#NEED_PKGS[@]}" -gt 0 ]; then
  say "Missing system packages: ${NEED_PKGS[*]}"
  if command -v sudo >/dev/null && command -v apt-get >/dev/null; then
    sudo apt-get update -qq && sudo apt-get install -y -qq "${NEED_PKGS[@]}"
  else
    warn "no sudo/apt; install these packages manually: ${NEED_PKGS[*]}"
  fi
fi

# ---------------------------------------------------------------- 2) Miniforge
if [ ! -x "$MINIFORGE_DIR/bin/conda" ]; then
  say "Installing Miniforge → $MINIFORGE_DIR"
  TMP_SH="$(mktemp /tmp/miniforge.XXXXXX.sh)"
  curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -o "$TMP_SH"
  bash "$TMP_SH" -b -p "$MINIFORGE_DIR"
  rm -f "$TMP_SH"
else
  say "Miniforge present: $MINIFORGE_DIR"
fi
CONDA="$MINIFORGE_DIR/bin/conda"

# ---------------------------------------------------------------- 3) environment + packages
if [ ! -x "$MINIFORGE_DIR/envs/$ENV_NAME/bin/python" ]; then
  say "Creating environment: $ENV_NAME (Python $PY_VER)"
  "$CONDA" create -y -q -n "$ENV_NAME" "python=$PY_VER" pip
else
  say "Environment present: $ENV_NAME"
fi
PY="$MINIFORGE_DIR/envs/$ENV_NAME/bin/python"

say "torch / torchvision ($([ "$CPU_ONLY" = 1 ] && echo CPU || echo CUDA 12.4) wheels)"
"$PY" -m pip install -q --upgrade pip
"$PY" -m pip install -q "torch==2.5.1" "torchvision==0.20.1" --index-url "$TORCH_INDEX"

say "requirements.txt"
# torch lines are already satisfied; pip only verifies the pins in requirements.txt.
"$PY" -m pip install -q -r "$HERE/requirements.txt"

say "Pinning numpy 1.26.4 (gtsam ABI)"
"$PY" -m pip install -q "numpy==1.26.4"

# ---------------------------------------------------------------- 4) verification
say "Verification"
"$PY" - <<'PYEOF'
import importlib, sys
mods = ["numpy", "scipy", "pandas", "yaml", "cv2", "torch", "torchvision", "kornia", "gtsam", "matplotlib"]
bad = []
for m in mods:
    try:
        mod = importlib.import_module(m)
        print(f"  {m:12s} {getattr(mod, '__version__', 'ok')}")
    except Exception as e:  # noqa: BLE001
        bad.append(m); print(f"  {m:12s} ERROR: {type(e).__name__}: {e}")
import numpy, torch
assert numpy.__version__.startswith("1."), "numpy 2.x was installed — incompatible with gtsam"
print(f"  CUDA visible: {torch.cuda.is_available()}" + (f" ({torch.cuda.get_device_name(0)})" if torch.cuda.is_available() else " — CPU mode"))
sys.exit(1 if bad else 0)
PYEOF

cat <<EOF

Done. Usage:

    source $MINIFORGE_DIR/bin/activate $ENV_NAME
    export SLAM_DATA_ROOT=/data/euroc              # EuRoC ASL layout (see README)
    python run_frontend.py && python run_backend.py

Pangolin is optional; see the note in requirements.txt if you want the 3-D trajectory window.
EOF
