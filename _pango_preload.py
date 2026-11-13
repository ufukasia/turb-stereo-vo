"""pypangolin'i güvenli import edilebilir hale getirir.
Anaconda'nın eski libstdc++/libgcc/glib'i sistem Pangolin'ini gölgeleyip import'u düşürebilir;
sistem kütüphaneleri varsa LD_PRELOAD ile öne alınır ve süreç bir kez yeniden başlatılır.
"""
import os
import sys


def ensure():
    if os.environ.get("_PANGO_PRELOAD_DONE"):
        return

    candidates = [
        "/usr/lib/x86_64-linux-gnu/libGL.so.1",
        "/usr/lib/x86_64-linux-gnu/libstdc++.so.6",
        "/lib/x86_64-linux-gnu/libgcc_s.so.1",
        "/lib/x86_64-linux-gnu/libglib-2.0.so.0",
        "/lib/x86_64-linux-gnu/libgobject-2.0.so.0",
    ]
    preload = " ".join(p for p in candidates if os.path.exists(p))
    if preload:
        os.environ["_PANGO_PRELOAD_DONE"] = "1"
        os.environ["LD_PRELOAD"] = (preload + " " + os.environ.get("LD_PRELOAD", "")).strip()
        os.execv(sys.executable, [sys.executable] + sys.argv)

