"""On uc || arka uc boru hatti: on uc (GPU) ana surecte kareleri uretir, asil arka uc (gtsam_backend.main) ayri
surecte kuyruktan okur; iki asama ortusur, duvar suresi kare basina max(on uc, arka uc) olur.

Arka uca giden kare sozlugu dokume yazilanin aynisidir (timestamp, obs, s_d, kps, desc); arka uc kareleri ayni sirada
aldigi icin yorunge, dokum + arka uc sirali calismasiyla ayni olur. Cocuk surec `spawn` ile baslar (CUDA baglami
miras alinmaz), yalniz CPU kullanir.

    python boru_hatti.py --out traj.csv [config --set ...]   (ayarlar config.yaml'dan)
"""
import os
import sys
import time
import multiprocessing as mp

KUYRUK = 32   # geri basinc: on uc arka uctan en cok bu kadar kare ileri gidebilir


class AkisListesi:
    """Kuyruktan gelen kareleri liste gibi sunar: len() bilinen kare sayisi, yineleme sirayla kuyruktan okur."""

    def __init__(self, q, n):
        self.q, self.n = q, n

    def __len__(self):
        return self.n

    def __getitem__(self, dilim):
        if isinstance(dilim, slice) and (dilim.start or 0) == 0 and dilim.step is None \
                and (dilim.stop is None or dilim.stop >= self.n):
            return self
        raise TypeError("AkisListesi yalniz bastan sona yinelenir")

    def __iter__(self):
        while True:
            fr = self.q.get()
            if fr is None:
                return
            yield fr


def _arka_uc(cfg, out_csv, q, sonuc_q):
    try:
        import config
        import gtsam_backend
        baslik, n = q.get()
        args = config.backend_args(cfg, argv=["--out", out_csv])
        args._veri = dict(baslik, frames=AkisListesi(q, n))
        t0 = time.time()
        gtsam_backend.main(args)
        sonuc_q.put(('bitti', time.time() - t0))
    except Exception as e:                                      # noqa: BLE001
        import traceback
        sonuc_q.put(('hata', f"{type(e).__name__}: {e}\n{traceback.format_exc()}"))


def calistir(cfg, out_csv):
    import config
    from stereo_front_end import StereoMatcher
    ds, fe = cfg["dataset"], cfg["frontend"]
    m = StereoMatcher(base_path=ds["base_path"], sequence=ds["sequence"], matcher_type=fe["matcher"],
                      filter_method=fe["filter"], y_axis_threshold=fe["y_axis_threshold"],
                      extractor=fe["extractor"], max_track_length=fe["max_track_length"],
                      max_keypoints=fe["max_keypoints"], detection_threshold=fe["detection_threshold"],
                      stereo_extrinsics=fe.get("stereo_extrinsics", "v1_gomulu"))
    config.apply_frontend(m, cfg)
    n = len(sorted(set(m.cam0_df['timestamp'].astype(str))))
    if ds["max_frames"]:
        n = min(n, int(ds["max_frames"]))

    ctx = mp.get_context('spawn')
    q, sonuc_q = ctx.Queue(maxsize=KUYRUK), ctx.Queue()
    p = ctx.Process(target=_arka_uc, args=(cfg, out_csv, q, sonuc_q), daemon=True)
    p.start()
    import queue as _queue

    def _koy(x):
        # arka uc dusmusse on uc dolu kuyrukta sonsuza dek beklemesin: hatasini yukselt
        while True:
            try:
                q.put(x, timeout=1.0)
                return
            except _queue.Full:
                if not p.is_alive():
                    durum_, deger_ = sonuc_q.get(timeout=5) if not sonuc_q.empty() else ('hata', 'arka uc sureci dustu')
                    raise RuntimeError(f"arka uc: {deger_}")

    baslik = dict(m.get_calibration_export(), has_keypoints=bool(getattr(m, 'EXPORT_KEYPOINTS', False)))
    _koy((baslik, n))
    sayac = {'n': 0, 'bos': 0}

    def on_frame(ts, obs, _idx):
        # dokum sozlugunun aynisi (stereo_front_end._save_landmark_export)
        fr = {'timestamp': str(ts), 'obs': obs}
        if m.EXPORT_S_D and ts in m._sd_by_ts:
            fr['s_d'] = m._sd_by_ts[ts]
        if getattr(m, 'EXPORT_KEYPOINTS', False) and m.keypoint_export and str(m.keypoint_export[-1][0]) == str(ts):
            _, kps, d = m.keypoint_export[-1]
            fr['kps'] = kps
            fr['desc'] = d
        _koy(fr)
        sayac['n'] += 1

    t0 = time.time()
    dokum = os.environ.get("BORU_DOKUM") or os.path.join(os.path.dirname(os.path.abspath(out_csv)), "boru_dokum.pkl")
    m.process_sequence(max_frames=ds["max_frames"], visualize=False, export_landmarks_path=dokum, on_frame=on_frame)
    t_on = time.time() - t0
    _koy(None)
    durum, deger = sonuc_q.get()
    p.join()
    t_top = time.time() - t0
    if durum == 'hata':
        raise RuntimeError(deger)
    k = max(sayac['n'], 1)
    print(f"BORU HATTI: {sayac['n']} kare, on uc {t_on:.1f} s ({t_on / k * 1e3:.1f} ms/kare), "
          f"arka uc {deger:.1f} s, toplam duvar {t_top:.1f} s ({t_top / k * 1e3:.1f} ms/kare)")
    return {'n': sayac['n'], 't_on': t_on, 't_arka': deger, 't_toplam': t_top}


if __name__ == "__main__":
    import config
    out = sys.argv[sys.argv.index("--out") + 1]
    calistir(config.load(), out)
