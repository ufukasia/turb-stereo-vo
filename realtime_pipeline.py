"""Ön uç ∥ arka uç boru hattı: arka uç ayrı süreçte koşar, ön uçla örtüşür.
GTSAM sarmalayıcısı GIL'i bırakmadığı için thread değil süreç kullanılır; `spawn` ile çocuk
ebeveynin CUDA bağlamını miras almaz. Aynı backend aynı gözlemleri aynı sırada aldığı için
yörünge ve ATE bit düzeyinde aynıdır. Tek görünür fark: canlı çizim kuyruk derinliği kadar gecikir.
"""

import os
import sys
import time
import importlib.util
from pathlib import Path

# Çocuk saf CPU'dur (GTSAM+numpy); BLAS tek thread'e sabitlenir, çok thread burada yavaşlatıyor.
_CHILD_THREAD_ENV = {
    'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
    'NUMEXPR_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1',
}

QUEUE_MAX = 24   # geri-basınç: ön uç arka uçtan bu kadar kare ileri gidebilir (~250 KB)


def _load_rtv():
    """`run_realtime_vision.py`'yi ana bloğunu çalıştırmadan yükle; çocukta yalnız RealtimeBackend gerekir."""
    here = Path(__file__).resolve().parent
    src = (here / 'run_realtime_vision.py').read_text()
    src = src.split("if __name__ ==")[0]
    spec = importlib.util.spec_from_loader('rtv_pipe', loader=None)
    mod = importlib.util.module_from_spec(spec)
    mod.__dict__['__file__'] = str(here / 'run_realtime_vision.py')
    sys.path.insert(0, str(here))
    os.chdir(here)
    exec(compile(src, 'run_realtime_vision.py', 'exec'), mod.__dict__)
    return mod


def backend_worker(calib, bargs, gt, gt_path, out_csv, q_in, q_out):
    """Arka uç süreci: q_in'den ('f', ts, obs, kps, desc, want_traj) al ve step() çalıştır; None gelince finalize + ATE.
    q_out'a ('traj', xyz_list), ('done', ate, be_seconds, n_frames) ya da ('err', mesaj) yazar."""
    for k, v in _CHILD_THREAD_ENV.items():
        os.environ[k] = v
    os.environ.setdefault('_PANGO_PRELOAD_DONE', '1')   # çocukta execv/pango yeniden başlatma yok
    try:
        rtv = _load_rtv()
        be = rtv.RealtimeBackend(calib, bargs)
        gt_t, gt_pos, gt_quat = gt
        be_time, n = 0.0, 0
        while True:
            item = q_in.get()
            if item is None:
                break
            _, ts, obs, kps, desc, want_traj = item
            t0 = time.time()
            be.step(ts, obs)
            if be.loop_closer is not None and kps is not None:
                be.loop_add(kps, desc)
            be_time += time.time() - t0
            n += 1
            if want_traj:
                try:
                    q_out.put(('traj', be.display_trajectory(gt_t, gt_pos, gt_quat)))
                except Exception:      # noqa: BLE001 — çizim verisi koşumu ASLA düşürmez
                    pass
        ate = be.finalize(gt_path, out_csv)
        q_out.put(('done', ate, be_time, n))
    except BaseException as e:         # noqa: BLE001
        msg = f'{type(e).__name__}: {str(e).strip().splitlines()[0] if str(e).strip() else ""}'
        try:
            q_out.put(('err', msg))
        except Exception:              # noqa: BLE001
            pass
