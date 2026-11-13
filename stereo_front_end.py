"""Stereo ön uç: EuRoC/ASL biçimli stereo çifti okur, rektifiye eder, öznitelik çıkarır
(SuperPoint/ALIKED/SIFT + LightGlue), sol-sağ eşleştirip üçgenler ve kareler arası track'leri
yönetir. Stereo yalnız derinlikli track sayısı eşiğin altına düşünce koşar.
Girdi: <base_path>/<sequence>/mav0/{cam0,cam1}/{sensor.yaml,data.csv,data/}.
Çıktı: kare başına istatistik CSV'si, isteğe bağlı landmark pickle'ı ve görselleştirme PNG'leri."""
import cv2
import numpy as np
import pandas as pd
import yaml
import os
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
import torch
import argparse
from scipy.spatial import cKDTree
from types import SimpleNamespace
from lightglue import ALIKED, SuperPoint, SIFT, LightGlue
import traceback
import threading, queue as _queue
from collections import Counter


class _AsyncImageWriter:
    """Hata-ayıklama PNG'lerini ayrı thread'lerde yazar; cv2.imwrite GIL'i bıraktığı için
    ana döngü beklemez. Kuyruk sınırlı: disk yavaşsa ana döngü bloklanır, RAM şişmez."""

    def __init__(self, n_workers=2, maxsize=8):
        self.q = _queue.Queue(maxsize=maxsize)
        self.n_dropped = 0
        self._threads = [threading.Thread(target=self._run, daemon=True)
                         for _ in range(max(1, int(n_workers)))]
        for t in self._threads:
            t.start()

    def _run(self):
        while True:
            item = self.q.get()
            if item is None:
                self.q.task_done()
                return
            path, img = item
            try:
                cv2.imwrite(path, img)
            except Exception:      # noqa: BLE001 — görsel yazımı koşumu ASLA düşürmez
                pass
            finally:
                self.q.task_done()

    def submit(self, path, img):
        # Çağıran kopya vermeli; vis_img sonraki karede yeniden kullanılıyor.
        self.q.put((path, img))

    def close(self, timeout=60.0):
        """Kuyruğu boşalt, thread'leri durdur."""
        try:
            self.q.join()
        except Exception:          # noqa: BLE001
            pass
        for _ in self._threads:
            try:
                self.q.put_nowait(None)
            except Exception:      # noqa: BLE001
                pass
        for t in self._threads:
            t.join(timeout=timeout / max(1, len(self._threads)))


# SIFT sarmalayıcısı: hat image_size'ı float verir (SIFT int ister) ve cv2.SIFT eşiği
# kurulumda donar; iki uyarlama da çağrı yerlerine değil çıkarıcı nesnesine konuldu.
class _SiftConf(SimpleNamespace):
    """conf.detection_threshold yazıldığında cv2 dedektörünü de günceller."""

    def __setattr__(self, ad, deger):
        super().__setattr__(ad, deger)
        if ad == "detection_threshold":
            d = self.__dict__.get("_cv2_dedektor")
            if d is not None and hasattr(d, "setContrastThreshold"):
                d.setContrastThreshold(float(deger))


class _SIFTHat(SIFT):
    """Hattın çağrı sözleşmesine uyan SIFT."""

    def __init__(self, **conf):
        super().__init__(**conf)
        self.conf = _SiftConf(**vars(self.conf))
        self.conf.__dict__["_cv2_dedektor"] = self.sift

    #: SIFT çıktısındaki alanlar ve son boyutları (boş çıktı bunlarla kurulur).
    _ALANLAR = {"keypoints": 2, "scales": None, "oris": None,
                "descriptors": 128, "keypoint_scores": None}

    def _bos_cikti(self, ornek):
        """Sıfır anahtar noktalı ama şekli geçerli çıktı."""
        B = int(ornek.shape[0])
        d, dt = ornek.device, ornek.dtype
        out = {}
        for ad, son in self._ALANLAR.items():
            bicim = (B, 0) if son is None else (B, 0, son)
            out[ad] = torch.zeros(bicim, device=d, dtype=dt)
        return out

    def forward(self, data):
        sz = data.get("image_size")
        if sz is not None and torch.is_floating_point(sz):
            data = {**data, "image_size": sz.long()}
        try:
            return super().forward(data)
        except ValueError as e:
            # lightglue/sift.py:filter_dog_point --- cv2 SIFT hiç aday bulamadığında
            # np.maximum.at boş diziye yayınlanamıyor ve ÇIKARIM ÇÖKÜYOR. Karanlık
            # basamaklarda (θ=0.04, ortalama parlaklık <13 DN) aday sayısı fiilen
            # sıfırdır; çökme yüzünden kare hiç işlenmiyor, yörünge ilk pozda
            # donuyor ve arz açığı kolu n=0'ı GÖREMEDEN devre dışı kalıyordu.
            # Boş ama geçerli çıktı döndür: kol n=0 okuyup eşiği θ_min'e indirir,
            # taban çizgi kareyi anahtar noktasız geçer.
            if "broadcastable" not in str(e):
                raise
            return self._bos_cikti(data["image"])


# Arz açığı kapılı dedektör kolu: aday sayısı hedefin altına düşünce eşik düşer,
# arz yeterliyken çıktı taban çizgiyle aynı kalır. İkinci çıkarım yalnız açık karelerde.
class ArzKapiliCikarici(torch.nn.Module):
    """Aday arzı referansın altına düştüğünde eşiği düşüren çıkarıcı sarmalayıcısı."""

    def __init__(self, ic, referans_n, thr_max, thr_min, pay=1.0, us=0.0):
        super().__init__()
        self.ic = ic
        # us<=0: doğrusal th = th_max - f*(th_max-th_min); us>0: çarpımsal th = th_max*(n/T)^us
        self.us = float(us)
        # Kol, teslim edilen sayı pay*referans'ın altına düşünce açılır (pay<1 ölü bant verir).
        self.referans_n = max(1, int(round(float(pay) * int(referans_n))))
        self.thr_max = float(thr_max)
        self.thr_min = float(thr_min)
        self.n_tetik = 0          # kolun açıldığı kare sayısı
        self.n_kare = 0
        self.n_geridusme = 0      # kütüphane geri düşmesi yakalanan kare sayısı

    @property
    def conf(self):
        """Eski kod self.aliked.conf okuyor; iç çıkarıcıya devret."""
        return self.ic.conf

    def _esik_yaz(self, v):
        self.ic.conf.detection_threshold = float(v)

    def forward(self, data):
        self.n_kare += 1
        self._esik_yaz(self.thr_max)
        cikti = self.ic(data)
        kp = cikti.get('keypoints')
        n = int(kp.shape[-2]) if kp is not None and kp.dim() >= 2 else 0
        # Geri düşme: kütüphane θ_max'ı aşan tepe bulamayınca eşiği sessizce ortalamaya indirir
        # ve havuz gürültüyle dolar. Tüm skorlar θ_max altındaysa kare aç sayılır (n=0).
        sc = cikti.get('keypoint_scores')
        if n > 0 and sc is not None and sc.numel() > 0 and float(sc.max()) < self.thr_max:
            self.n_geridusme += 1
            n = 0
        if n >= self.referans_n:
            return cikti                      # arz yeterli, taban çizgi
        if self.us > 0:
            # Çarpımsal kural; thr_min kelepçesi olmadan eşik gürültü tabanına iner.
            esik = self.thr_max * (max(n, 1) / float(self.referans_n)) ** self.us
            esik = max(self.thr_min, min(esik, self.thr_max))
        else:
            f = 1.0 - (n / float(self.referans_n))
            esik = self.thr_max - f * (self.thr_max - self.thr_min)
        self.n_tetik += 1
        self._esik_yaz(esik)
        try:
            return self.ic(data)
        finally:
            self._esik_yaz(self.thr_max)      # durum sızmasın


class TabanOranliCikarici(torch.nn.Module):
    """Eşiği karenin kendi skor dağılımından seçen sarmalayıcı: θ = clamp(s_(T), κ·s̄, θ_max),
    s_(T) = T'inci en büyük NMS tepesi, s̄ = harita ortalaması. Tek geçiş, aydınlıkta taban çizgi."""

    def __init__(self, ic, hedef_n, thr_max, kappa=2.0, gamma=1.0):
        super().__init__()
        self.ic = ic
        self.hedef_n = max(1, int(hedef_n))
        self.thr_max = float(thr_max)
        self.kappa = float(kappa)      # taban tepkinin kaç katı (alt kelepçe)
        self.gamma = float(gamma)      # θ ≤ gamma·s_max
        self.n_kare = 0
        self.n_tetik = 0               # θ'nın tavandan indiği kare sayısı
        self.n_aclik = 0               # taban kelepçesine dayanan kare sayısı
        self._son_esik = float(thr_max)

    @property
    def conf(self):
        return self.ic.conf

    def forward(self, data):
        from lightglue.aliked import simple_nms
        from kornia.color import grayscale_to_rgb
        self.n_kare += 1
        net = self.ic
        # Yoğun harita bir kez çıkarılır; eşik bunun üzerinden seçilir.
        img = data["image"]
        if img.shape[1] == 1:
            img = grayscale_to_rgb(img)
        _fm, sm = net.extract_dense_map(img)
        nms = simple_nms(sm.detach(), net.dkd.radius)
        duz = nms.reshape(nms.shape[0], -1)[0]
        s_ort = float(sm.detach().mean())
        s_max = float(duz.max())
        k = min(self.hedef_n, duz.numel())
        s_T = float(torch.topk(duz, k).values[-1])      # T'inci en büyük tepe

        alt = self.kappa * s_ort
        esik = min(max(s_T, alt), self.thr_max, self.gamma * s_max)
        if esik < self.thr_max:
            self.n_tetik += 1
        if s_T < alt:
            self.n_aclik += 1
        self._son_esik = esik

        eski = net.conf.detection_threshold
        net.conf.detection_threshold = float(esik)
        try:
            return net(data)
        finally:
            net.conf.detection_threshold = eski


def _cuda_oom(e):
    """CUDA bellek yetersizligi: altyapi hatasi; kare hatasi gibi yutulmaz, kosumu durdurur."""
    return isinstance(e, getattr(torch.cuda, 'OutOfMemoryError', ())) or 'out of memory' in str(e).lower()


def frame2tensor(frame, device):
    """Gri görüntüyü (HxW, 0-255) normalize tensöre çevir."""
    return torch.from_numpy(frame / 255.).float()[None, None].to(device)


def landmarks_to_obs(lms):
    """Aktif landmark listesini backend obs dizisine çevir:
    (N,8) float64 [id, u, v, has_depth, fresh_depth, X, Y, Z]; derinlik yoksa XYZ=NaN."""
    if not lms:
        return np.zeros((0, 8), dtype=np.float64)
    arr = np.full((len(lms), 8), np.nan, dtype=np.float64)
    for i, lm in enumerate(lms):
        arr[i, 0] = lm['id']
        arr[i, 1] = lm['kp'][0]
        arr[i, 2] = lm['kp'][1]
        arr[i, 3] = 1.0 if lm['has_depth'] else 0.0
        arr[i, 4] = 1.0 if lm['fresh_depth'] else 0.0
        if lm['p3d'] is not None:
            arr[i, 5:8] = lm['p3d']
    return arr

class StereoMatcher:
    # Varsayılan veri kökü $SLAM_DATA_ROOT (yoksa ./data); çağıranlar base_path'i açıkça verir.
    def __init__(self, base_path=None, sequence="MH01", matcher_type="lightglue", filter_method="y_axis_check", y_axis_threshold=8.0, extractor="superpoint", max_track_length=0, max_keypoints=1000, detection_threshold=None, stereo_extrinsics="gomulu"):
        if base_path is None:
            base_path = os.environ.get("SLAM_DATA_ROOT", "data")
        """
        StereoMatcher sınıfını başlat.

        Args:
            base_path (str): Veri setinin temel yolu. Varsayılan: "".
            sequence (str): İşlenecek EurocMAV sekans adı. Varsayılan: "MH_05_difficult".
            matcher_type (str): Eşleştirici türü.
                                Seçenekler: "lightglue", "classical_bf".
                                Varsayılan: "lightglue".
            filter_method (str): Filtreleme yöntemi.
                                 Seçenekler: "fundamental_matrix", "y_axis_check", "cross_check", "none".
                                 Varsayılan: "cross_check".
            y_axis_threshold (float): "y_axis_check" filtresi için eşik (piksel). Varsayılan: 8.0.
        """
        self.base_path = Path(base_path)
        self.sequence = sequence
        # cam0->cam1 donusumunun kaynagi: 'sensor_yaml' = inv(T_BS_cam1) @ T_BS_cam0 (sensorun kendi
        # kalibrasyonu); 'gomulu' = koda gomulu eski matris (EuRoC'ta 0.83 derece yunuslama + 0.10 derece
        # sapma hatali; teshis: dev/teshis_rekt.py). cam1/sensor.yaml'da T_cn_cnm1 varsa her iki kipte o kullanilir.
        self.STEREO_EXTRINSICS = str(stereo_extrinsics)
        self.results = []
        
        # Temporal takip ve görselleştirme durumu
        self.prev_vis_info = { # önceki karenin görselleştirme bilgileri
            'image': None,
            'kp_candidates': None,
            'filename0': None
        }
        self.tracked_features = [] # aktif track'ler: {id, kp, desc, p3d, age, last_seen_frame_idx, ...}
        self.current_frame_idx = -1
        self.next_track_id = 0
        self.MAX_TRACK_AGE_UNSEEN = 5  # kaybolan track kaç kare tutulur
        # Track ömür sınırı (kare); iSAM2'de dev klik oluşmasın diye turnover zorlar. 0 = kapalı.
        self.MAX_TRACK_LENGTH = int(max_track_length)
        self.REPLENISH_THRESHOLD = 50  # derinlikli track sayısı bunun altına düşerse stereo tetiklenir
        self.TARGET_FEATURE_COUNT = 200 # yeni nokta eklenirken hedef toplam

        # Ölçek-duyarlı stereo tetikleme; ablasyonda ATE'yi kötüleştirdi, varsayılan kapalı.
        self.SCALE_AWARE_TRIGGER = False
        # (A) yakın anchor kıtlığı
        self.SCALE_NEAR_TRIGGER = False
        self.SCALE_NEAR_DEPTH = 8.0       # m; bundan yakın anchor iyi koşullu sayılır
        self.SCALE_MIN_NEAR_ANCHORS = 12  # yakın anchor bundan azsa tetikle
        # (B) eski anchor kadansı (sabit kadansa dönüşüyor, adaptif değil)
        self.SCALE_STALE_TRIGGER = False
        self.SCALE_MAX_DEPTH_AGE = 8      # kare; anchor'ın son üçgenlemeden beri yaşı
        self.SCALE_MAX_STALE_FRAC = 0.75  # derinlikli track'lerin bu oranı eskiyse tetikle
        self.n_scale_triggers = 0
        # Izgara tabanlı nokta seçimi
        self.GRID_ROWS = 10
        self.GRID_COLS = 12

        # Backend'e gidecek landmark'lar için geometrik eşikler
        self.MIN_STEREO_DISPARITY = 0.0   # pozitif disparity zorunlu
        self.MAX_LANDMARK_DEPTH = 50.0    # m
        self.TEMPORAL_RANSAC_THRESHOLD = 2.0  # temporal RANSAC eşiği (px)
        # Parallaks açısı kapısı (apply_frontend ile açılır)
        self.PARALLAX_CHECK = False
        self.COS_MAX_PARALLAX = 0.9998    # cos(parallaks) bunun üstündeyse landmark elenir
        # Kare başına keypoint+descriptor dışa aktarımı
        self.EXPORT_KEYPOINTS = False
        self.keypoint_export = []         # (ts, kps Nx2, desc NxD) listesi

        # Track ayıklama; hepsi varsayılan kapalı, config.frontend.* ile açılır.
        # Sayaçlar cull_stats'ta birikir.
        # (1) RANSAC'ın N kez dışladığı track kalıcı silinir
        self.CULL_RANSAC_STRIKES = 0         # 0=kapalı
        # (2) güncel ve anchor descriptor arasındaki kosinüs mesafesi eşiği aşarsa sil
        self.CULL_DESC_DRIFT_MAX = 0.0       # 0=kapalı; 0..2
        # (3) ömür boyunca eşleşme oranı düşük (titreyen) track silinir
        self.CULL_FLICKER_RATIO = 0.0        # 0=kapalı
        self.CULL_FLICKER_MIN_LIFE = 10      # bu kadar kare yaşamadan yargılanmaz
        # (4) stereo üçgenleme reprojeksiyon kontrolü
        self.STEREO_REPROJ_CHECK = False
        self.STEREO_REPROJ_MAX_PX = 2.0
        # (5) göreli derinlik belirsizliği kapısı: σ_Z/Z = Z·σ_d/(f·b)
        self.MAX_DEPTH_SIGMA_REL = 0.0       # 0=kapalı
        self.STEREO_DISP_SIGMA_PX = 1.0      # disparity gürültüsü σ_d (px)
        self.cull_stats = Counter()

        # Yeni anahtarlar; varsayilanlar duzeltme oncesi davranisi verir (config.apply_frontend yazar).
        self.MIN_LANDMARK_DEPTH = 0.0            # m; ORB-SLAM3 minZ = taban (Frame.cc:841)
        self.TEMPORAL_MATCHER = 'bf_lowe'        # 'bf_lowe' | 'lightglue'
        self.TEMPORAL_LOWE_RATIO = 0.8           # Lowe (2004) oran testi
        self.TEMPORAL_RANSAC_CONF = 0.99
        self.TEMPORAL_RANSAC_ITERS = 2000
        self.STEREO_DEPTH_ASSOC_PX = 1.5         # stereo noktasi <-> track ayni anahtar nokta (NMS r=2 => >=3 px)
        self.STEREO_TRACK_RATE = 0.0             # >0: derinlikli track / son stereo karedeki sayi < oran -> stereo
        self.TRIGGER_CLOSE_DEPTH = 0.0           # >0: stereo tetigi yalniz bu derinlikten yakin track'leri sayar
        self.BOS_KARE_DUZELT = False             # basarisiz kareyi bos gozlemle aktar (hata duzeltmesi)
        self.n_bos_kare = 0
        self.n_stereo_geri_dusme = 0
        self._depth_at_last_stereo = 0
        self._prev_frame_feats = None            # LightGlue zamansal eslestirme icin onceki karenin ozellikleri
        # Aday 1: kapili zamansal LightGlue (karanlik kapisi acikken BF+Lowe yerine LightGlue). Kapaliyken ozdes.
        self.TEMPORAL_LG_KAPILI = False
        # Saglikla tetiklenen zamansal LightGlue: onceki karede gorulen izlerin devam orani bu esigin altindaysa bu
        # karede BF+Lowe yerine LightGlue (AirSLAM kareler arasi LightGlue). 0 = kapali (ozdes).
        self.TEMPORAL_LG_SAGLIK = 0.0
        self._son_devam_orani = 1.0
        self._kapi_acik = False
        # Aday 3: karenin S_d'si doku dosyasina yazilir (arka uc karanlik gurultu modeli icin). Kapaliyken ozdes.
        self.EXPORT_S_D = False
        self._sd_by_ts = {}
        self.n_lg_temporal = 0


        self.execution_time = None


        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.matcher_type = matcher_type
        self.filter_method = filter_method
        self.y_axis_threshold = y_axis_threshold

        # Çıkarıcı seçimi. self.aliked adı geriye uyumluluk için; içinde seçilen çıkarıcı durur.
        self.extractor = (extractor or "superpoint").lower()
        if self.extractor == "superpoint":
            if SuperPoint is None:
                raise ImportError("SuperPoint bu depoda dağıtılmaz (lisans; bkz. NOTICE). "
                                  "lightglue/superpoint.py dosyasını cvg/LightGlue'dan ekleyin ya da "
                                  "frontend.extractor: aliked kullanın.")
            thr = 0.0005 if detection_threshold is None else float(detection_threshold)
            self.aliked = SuperPoint(max_num_keypoints=int(max_keypoints), detection_threshold=thr).eval().to(self.device)
            lg_features = 'superpoint'
        elif self.extractor == "aliked":
            thr = 0.02 if detection_threshold is None else float(detection_threshold)
            self.aliked = ALIKED(max_num_keypoints=int(max_keypoints), detection_threshold=thr).eval().to(self.device)
            lg_features = 'aliked'
        elif self.extractor == "sift":
            # detection_threshold -> cv2.SIFT contrastThreshold; varsayılan COLMAP değeri.
            thr = 0.0066667 if detection_threshold is None else float(detection_threshold)
            self.aliked = _SIFTHat(max_num_keypoints=int(max_keypoints),
                                   detection_threshold=thr,
                                   backend="opencv").eval().to(self.device)
            lg_features = 'sift'
        else:
            raise ValueError(f"Geçersiz extractor: {self.extractor}. "
                             f"'superpoint', 'aliked' veya 'sift' olmalı.")
        print(f"Özellik çıkarıcı (extractor): {self.extractor}")

        if self.matcher_type == 'lightglue':
            print("LightGlue eşleştirici kullanılıyor.")
            self.lightglue_stereo = LightGlue(features=lg_features).eval().to(self.device)
            if hasattr(self.lightglue_stereo, 'conf') and hasattr(self.lightglue_stereo.conf, 'features'):
                self.extractor_name = self.lightglue_stereo.conf.features
            else:
                self.extractor_name = self.extractor
        elif self.matcher_type == 'classical_bf':
            print("Klasik BFMatcher kullanılıyor.")
            self.stereo_bf_matcher = cv2.BFMatcher(cv2.NORM_L2)  # float descriptor -> L2
            self.extractor_name = self.extractor
        else:
            raise ValueError(f"Geçersiz eşleştirici türü: {self.matcher_type}. 'lightglue' veya 'classical_bf' olmalı.")
        self.temporal_bf_matcher = cv2.BFMatcher(cv2.NORM_L2)

        # Adaptif ön-işleme: S_d = 1 - H/8 (normalize entropi) karanlık skoru.
        # CLAHE yalnız S_d > eşik iken, clipLimit S_d ile doğrusal; dedektör eşiği
        # thr = clamp(thr_max - S_d·(thr_max - thr_min)). Aydınlıkta taban çizgiyle aynı.
        # Varsayılan kapalı; config.frontend.adaptive_preprocessing ile açılır.
        self.ADAPTIVE_PREPROCESSING_ENABLED = False
        self.ADAPTIVE_CLAHE_THRESHOLD_S_D = 0.40   # S_d bunu aşarsa CLAHE
        self.ADAPTIVE_CLAHE_CLIP_LIMIT_MIN = 0.5
        self.ADAPTIVE_CLAHE_CLIP_LIMIT_MAX = 2.5
        self.ADAPTIVE_CLAHE_TILE = 8
        # Eşik tavanı = config'den gelen detection_threshold; karanlıkta tabana doğru gevşer.
        self.ADAPTIVE_DET_ENABLED = True
        self.ADAPTIVE_DET_MAX_THRESHOLD = float(thr)
        self.ADAPTIVE_DET_MIN_THRESHOLD = 1e-5
        # True: eşik kolu da fotometrik kolla aynı tau kapısına bağlı; False: her S_d>0 karede oynar.
        self.ADAPTIVE_DET_GATED = False
        self._clahe_cache = {}  # clipLimit (3 ondalık) -> cv2.CLAHE
        self._last_det_threshold = float(thr)  # son uygulanan eşik (rapor için)
        # [TÜRBÜLANS] bozulma enjektörü (turb.hooks.TurbulenceInjector); config.apply_frontend kurar.
        # None iken hiçbir ek hesap/kopya yapılmaz (E0: PyEco ile bit düzeyinde aynı).
        self.turb = None
        self._seq_frame_idx = 0   # sekans içi mutlak kare indeksi (türbülansın zaman ekseni)
        # Bütçe-eşlemeli CLAHE (BE-CLAHE; ARASTIRMA_KAYDI.md §2.1). Kapalıyken hiçbir yol değişmez.
        # Kapı aynı (sol S_d > tau). Açıkken c = C'de n0(c) >= kp olan EN KÜÇÜK clip (yoksa C[-1]);
        # n0 = ALIKED'in kendi eşiği (theta_max) üstündeki NMS tepeleri, top-k (n_limit) kesiminden önce.
        # Sol ve sağ aynı c'yi alır.
        self.BE_CLAHE = False
        self.BE_MERDIVEN = (0.5, 1.0, 1.5, 2.0, 2.5)
        self.BE_HEDEF = int(max_keypoints)
        self.BE_ARAMA = 'tam'          # 'tam' (kesin, artan sıra) | 'nedensel' (önceki c'den başla)
        self.BE_PROFIL = False         # True: her açık karede tüm C için n0 (yalnız kayıt; seçimi değiştirmez)
        self._be_onceki_i = 0
        # Seçilen basamağın CLAHE görüntüsü ve yoğun haritası saklanır; sol çıkarım haritayı yeniden hesaplamaz.
        # Kapalıyken davranış bayt-özdeş. Girdi tensörü birebir aynı değilse harita yeniden hesaplanır.
        self.BE_ONBELLEK = False
        self._be_secilen = None        # (clahe_img_u8, (t_rgb, feature_map, score_map)) arama sonucu
        self._be_onbellek = None       # sol çıkarımı bekleyen (t_rgb, feature_map, score_map)
        self._be_isabet = 0
        self._be_iska = 0

        print("Kalibrasyon parametreleri yükleniyor...")
        self.load_calibration()
        print("Görüntü zaman damgaları yükleniyor...")
        self.load_timestamps()

    def load_calibration(self):
        """Kamera kalibrasyonunu sensor.yaml dosyalarından yükle."""
        try:
            cam0_yaml = self.base_path / self.sequence / "mav0/cam0/sensor.yaml"
            cam1_yaml = self.base_path / self.sequence / "mav0/cam1/sensor.yaml"

            with open(cam0_yaml, 'r') as f:
                self.cam0_params = yaml.safe_load(f)
            with open(cam1_yaml, 'r') as f:
                self.cam1_params = yaml.safe_load(f)

            intrinsics1 = self.cam0_params['intrinsics']
            intrinsics2 = self.cam1_params['intrinsics']

            self.K1 = np.array([
                [intrinsics1[0], 0, intrinsics1[2]],
                [0, intrinsics1[1], intrinsics1[3]],
                [0, 0, 1]
            ])

            self.K2 = np.array([
                [intrinsics2[0], 0, intrinsics2[2]],
                [0, intrinsics2[1], intrinsics2[3]],
                [0, 0, 1]
            ])

            self.D1 = np.array(self.cam0_params['distortion_coefficients'])
            self.D2 = np.array(self.cam1_params['distortion_coefficients'])

            # cam0 -> cam1 dönüşümü. Varsayılan EuRoC; cam1/sensor.yaml'da T_cn_cnm1 varsa o kullanılır.
            baseline_matrix = np.array([
                [ 0.99999609,  0.0023145 , -0.00136141, -0.110073  ],
                [-0.00231508,  0.99999693, -0.00042232, -0.000399  ],
                [ 0.00136043,  0.00042547,  0.99999868, -0.000284  ],
                [ 0.        ,  0.        ,  0.        ,  1.        ]
            ])
            if 'T_cn_cnm1' in self.cam1_params:
                baseline_matrix = np.array(self.cam1_params['T_cn_cnm1'],
                                           dtype=np.float64).reshape(4, 4)
                print("T_cn_cnm1 sensor.yaml'dan okundu (EuRoC dışı veri kümesi).")
            elif self.STEREO_EXTRINSICS == 'sensor_yaml':
                # EuRoC: T_c1_c0 = inv(T_B_c1) @ T_B_c0, iki kameranin T_BS'sinden (ORB-SLAM3
                # Examples/Stereo/EuRoC.yaml Stereo.T_c1_c2 ile ayni donusumun tersi).
                _T0 = np.array(self.cam0_params['T_BS']['data'], dtype=np.float64).reshape(4, 4)
                _T1 = np.array(self.cam1_params['T_BS']['data'], dtype=np.float64).reshape(4, 4)
                baseline_matrix = np.linalg.inv(_T1) @ _T0
                print("cam0->cam1 donusumu iki T_BS'den hesaplandi.")
            elif self.STEREO_EXTRINSICS != 'gomulu':
                raise ValueError(f"stereo_extrinsics bilinmiyor: {self.STEREO_EXTRINSICS}")

            self.R = baseline_matrix[:3, :3]
            self.T = baseline_matrix[:3, 3]

            baseline_norm = np.linalg.norm(self.T)
            print("Baseline (cam0 to cam1) matris:\n", baseline_matrix)
            print(f"Baseline norm: {baseline_norm:.9f} [m]")

            print("Rektifikasyon haritaları hesaplanıyor...")
            self.compute_rectification()

        except Exception as e:
            print(f"Kalibrasyon yükleme hatası: {e}")
            raise

    def compute_rectification(self):
        """Stereo rektifikasyon haritalarını hesapla."""
        # Çözünürlük sensor.yaml'dan; yoksa EuRoC (752, 480).
        res = getattr(self, 'cam0_params', {}).get('resolution')
        img_size = (int(res[0]), int(res[1])) if res else (752, 480)

        self.R1, self.R2, self.P1, self.P2, self.Q, self.roi1, self.roi2 = cv2.stereoRectify(
            self.K1, self.D1, self.K2, self.D2,
            img_size, self.R, self.T,
            flags=cv2.CALIB_ZERO_DISPARITY,
            alpha=0
        )

        self.map1x, self.map1y = cv2.initUndistortRectifyMap(
            self.K1, self.D1, self.R1, self.P1,
            img_size, cv2.CV_32FC1)

        self.map2x, self.map2y = cv2.initUndistortRectifyMap(
            self.K2, self.D2, self.R2, self.P2,
            img_size, cv2.CV_32FC1)

    def load_timestamps(self):
        """Görüntü zaman damgalarını data.csv'den yükle."""
        try:
            cam0_csv = self.base_path / self.sequence / "mav0/cam0/data.csv"
            cam1_csv = self.base_path / self.sequence / "mav0/cam1/data.csv"

            self.cam0_df = pd.read_csv(cam0_csv)
            self.cam1_df = pd.read_csv(cam1_csv)

            # Damgalar metin olarak tutulur (sabit genişlik, sıralama için)
            self.cam0_df = self.cam0_df.rename(columns={'#timestamp [ns]': 'timestamp'})
            self.cam1_df = self.cam1_df.rename(columns={'#timestamp [ns]': 'timestamp'})
            self.cam0_df['timestamp'] = self.cam0_df['timestamp'].astype(str)
            self.cam1_df['timestamp'] = self.cam1_df['timestamp'].astype(str)

            print(f"Kamera 0'dan {len(self.cam0_df)} zaman damgası yüklendi")
            print(f"Kamera 1'den {len(self.cam1_df)} zaman damgası yüklendi")

        except Exception as e:
            print(f"CSV dosya yükleme hatası: {e}")
            raise

    def rectify_images(self, img1, img2):
        """Stereo çifti rektifiye et."""
        rect1 = cv2.remap(img1, self.map1x, self.map1y, cv2.INTER_LINEAR)
        rect2 = cv2.remap(img2, self.map2x, self.map2y, cv2.INTER_LINEAR)
        return rect1, rect2

    # Adaptif özellik çıkarımı: entropi -> CLAHE clipLimit + dedektör eşiği
    def calculate_shannon_entropy(self, image_cv2, roi=None):
        """8-bit gri görüntünün normalize Shannon entropisi [0,1]; H/8.
        roi=(x,y,w,h) verilirse yalnız o bölgede hesaplanır (teşhis amaçlı)."""
        if image_cv2 is None or image_cv2.size == 0:
            return 0.0

        if image_cv2.ndim > 2 and image_cv2.shape[2] > 1:
            gray_image = cv2.cvtColor(image_cv2, cv2.COLOR_BGR2GRAY)
        else:
            gray_image = image_cv2

        if roi is not None:
            x, y, w, h = (int(v) for v in roi)
            if w > 0 and h > 0:
                gray_image = gray_image[y:y + h, x:x + w]
        if gray_image.size == 0:
            return 0.0

        hist = cv2.calcHist([gray_image], [0], None, [256], [0, 256]).ravel().astype(np.float64)
        total = hist.sum()
        if total <= 0:
            return 0.0
        p = hist / total
        nz = p > 0
        entropy = float(-np.sum(p[nz] * np.log2(p[nz])))

        normalized_entropy = entropy / 8.0  # 8 = log2(256)
        return min(max(normalized_entropy, 0.0), 1.0)

    def _adaptive_preprocess_image(self, image, roi=None, image_name_for_log="image"):
        """Entropi güdümlü CLAHE. Dönüş: (görüntü, stats); clahe_clip=0 uygulanmadı demek."""
        norm_entropy = self.calculate_shannon_entropy(image, roi=None)
        s_d = 1.0 - norm_entropy
        # ROI varyantı yalnız teşhis için
        norm_entropy_roi = self.calculate_shannon_entropy(image, roi=roi) if roi is not None else float('nan')
        s_d_roi = 1.0 - norm_entropy_roi if roi is not None else float('nan')

        output_image = image
        applied_clip_limit = 0.0

        s_d_thr = self.ADAPTIVE_CLAHE_THRESHOLD_S_D
        if s_d > s_d_thr:
            # [S_d_thr, 1] -> [0, 1]
            s_d_capped = min(max(s_d, s_d_thr), 1.0)
            denominator = 1.0 - s_d_thr
            if denominator > 1e-6:
                clahe_strength_factor = (s_d_capped - s_d_thr) / denominator
            else:
                clahe_strength_factor = 1.0

            clip_min = self.ADAPTIVE_CLAHE_CLIP_LIMIT_MIN
            clip_max = self.ADAPTIVE_CLAHE_CLIP_LIMIT_MAX
            current_clip_limit = clip_min + clahe_strength_factor * (clip_max - clip_min)
            current_clip_limit = max(clip_min, min(current_clip_limit, clip_max))
            applied_clip_limit = current_clip_limit

            output_image = self._get_clahe(current_clip_limit).apply(image)

        stats = {
            's_d': s_d,
            'norm_entropy': norm_entropy,
            'clahe_clip': applied_clip_limit,
            's_d_roi': s_d_roi,
            'norm_entropy_roi': norm_entropy_roi,
        }
        return output_image, stats

    def _get_clahe(self, clip_limit):
        """clipLimit başına önbellekli cv2.CLAHE nesnesi."""
        key = round(float(clip_limit), 3)
        clahe = self._clahe_cache.get(key)
        if clahe is None:
            tile = int(self.ADAPTIVE_CLAHE_TILE)
            clahe = cv2.createCLAHE(clipLimit=key, tileGridSize=(tile, tile))
            self._clahe_cache[key] = clahe
        return clahe



    @staticmethod
    def _lg_scale_ori(n, sc, ori, device):
        """LightGlue add_scale_ori girdisi: dedektör ölçek/yönelim veriyorsa onu, yoksa (1, 0) ver."""
        if sc is not None and ori is not None and sc.shape[0] == n:
            return sc.reshape(1, n).float(), ori.reshape(1, n).float()
        if n > 0:
            return (torch.ones(1, n, device=device, dtype=torch.float32),
                    torch.zeros(1, n, device=device, dtype=torch.float32))
        return (torch.empty(1, 0, device=device, dtype=torch.float32),
                torch.empty(1, 0, device=device, dtype=torch.float32))

    def _adapt_detection_threshold(self, s_d):
        """Dedektör eşiğini S_d ile doğrusal gevşet: thr = clamp(thr_max - S_d·(thr_max - thr_min)).
        Çıkarıcılar conf.detection_threshold'u her çağrıda okur, o yüzden conf'a yazmak yeter."""
        thr_max = self.ADAPTIVE_DET_MAX_THRESHOLD
        thr_min = self.ADAPTIVE_DET_MIN_THRESHOLD
        if getattr(self, 'ADAPTIVE_DET_GATED', False):
            # Kapılı kip: S_d <= tau iken tavan; üstünde [tau,1] -> [0,1] (CLAHE ile aynı)
            tau = self.ADAPTIVE_CLAHE_THRESHOLD_S_D
            if s_d <= tau:
                f = 0.0
            else:
                payda = 1.0 - tau
                f = ((min(max(s_d, tau), 1.0) - tau) / payda) if payda > 1e-6 else 1.0
            new_threshold = thr_max - f * (thr_max - thr_min)
        else:
            new_threshold = thr_max - s_d * (thr_max - thr_min)
        new_threshold = max(thr_min, min(new_threshold, thr_max))

        if hasattr(self.aliked, 'conf') and hasattr(self.aliked.conf, 'detection_threshold'):
            self.aliked.conf.detection_threshold = new_threshold
            self._last_det_threshold = new_threshold
        else:
            print(f"  Uyarı: {self.extractor} conf.detection_threshold ayarlanamadı; "
                  f"hesaplanan eşik {new_threshold:.6f} UYGULANMADI.")
        return new_threshold

    def _kare_kapi_kaydet(self, timestamp, adapt):
        """Karenin karanlik kapisi (sol S_d > tau) ve S_d'si; yalniz aday bayraklari aciksa kullanilir."""
        s_d = adapt.get('s_d1', float('nan'))
        self._kapi_acik = bool(s_d == s_d and s_d > self.ADAPTIVE_CLAHE_THRESHOLD_S_D)
        if self.EXPORT_S_D and s_d == s_d:
            self._sd_by_ts[timestamp] = float(s_d)

    def _preprocess_stereo_pair(self, rect1, rect2):
        """Stereo çift için adaptif ön işleme + eşik uyarlaması. Dönüş: (out1, out2, adapt).
        Adaptasyon kapalıysa görüntüler değişmeden döner."""
        adapt = {
            's_d1': float('nan'), 's_d2': float('nan'),
            's_d1_roi': float('nan'), 's_d2_roi': float('nan'),
            'norm_entropy1': float('nan'), 'norm_entropy2': float('nan'),
            'clahe_clip1': 0.0, 'clahe_clip2': 0.0,
            'det_threshold': self._last_det_threshold,
            'adaptive': False,
        }
        if not self.ADAPTIVE_PREPROCESSING_ENABLED:
            return rect1, rect2, adapt
        if self.BE_CLAHE:
            return self._be_clahe_cift(rect1, rect2, adapt)

        out1, st1 = self._adaptive_preprocess_image(rect1, getattr(self, 'roi1', None), "rect1")
        adapt.update(s_d1=st1['s_d'], s_d1_roi=st1['s_d_roi'],
                     norm_entropy1=st1['norm_entropy'], clahe_clip1=st1['clahe_clip'],
                     adaptive=True)

        out2 = rect2
        if rect2 is not None:
            out2, st2 = self._adaptive_preprocess_image(rect2, getattr(self, 'roi2', None), "rect2")
            adapt.update(s_d2=st2['s_d'], s_d2_roi=st2['s_d_roi'],
                         norm_entropy2=st2['norm_entropy'], clahe_clip2=st2['clahe_clip'])

        # Eşik sol görüntünün S_d'siyle uyarlanır; track'ler sol kamera çerçevesinde.
        if self.ADAPTIVE_DET_ENABLED:
            adapt['det_threshold'] = self._adapt_detection_threshold(st1['s_d'])

        print(f"  ADAPTİF: S_d(sol)={adapt['s_d1']:.3f} (ROI {adapt['s_d1_roi']:.3f}), "
              f"clip={adapt['clahe_clip1']:.2f}, det_thr={adapt['det_threshold']:.6f}")
        return out1, out2, adapt

    # ------------------------------------------------------------------ BE-CLAHE
    def _ic_cikarici(self):
        """Sarmalayıcıların (ArzKapiliCikarici / TabanOranliCikarici) içindeki ALIKED ağı."""
        net = self.aliked
        while hasattr(net, 'ic'):
            net = net.ic
        return net

    def _theta_max(self):
        """Dedektörün kendi eşiği (arz kolundan önceki tavan)."""
        a = self.aliked
        if hasattr(a, 'thr_max'):
            return float(a.thr_max)
        return float(self.ADAPTIVE_DET_MAX_THRESHOLD)

    def n0_say(self, img_u8, sakla=False):
        """n0: ALIKED skor haritasında theta_max üstündeki NMS tepeleri (kenar bandı hariç, DKD ile aynı),
        n_limit (top-k) kesiminden ÖNCE. Kütüphanenin 'tepe yoksa ortalamaya düş' yolu sayılmaz (n0=0)."""
        from lightglue.aliked import simple_nms
        from kornia.color import grayscale_to_rgb
        net = self._ic_cikarici()
        with torch.no_grad():
            t = frame2tensor(img_u8, self.device)
            t_rgb = grayscale_to_rgb(t)
            _fm, sm = net.extract_dense_map(t_rgb)
            nms = simple_nms(sm, net.dkd.radius)
            r = net.dkd.radius
            nms[:, :, :r, :] = 0
            nms[:, :, :, :r] = 0
            nms[:, :, -r:, :] = 0
            nms[:, :, :, -r:] = 0
            n0 = int((nms > self._theta_max()).sum().item())
            return (n0, (t_rgb, _fm, sm)) if sakla else n0

    def _harita_kancasi(self):
        """ALIKED.extract_dense_map'i sarar: girdi, aramada seçilen basamağın tensörüyle birebir aynıysa
        saklanan haritayı döndürür (tek kullanımlık), değilse asıl hesaplamayı yapar."""
        net = self._ic_cikarici()
        if getattr(net, '_be_kanca', False):
            return
        asil = net.extract_dense_map

        def kanca(image):
            ob, self._be_onbellek = self._be_onbellek, None
            if ob is not None and ob[0].shape == image.shape and torch.equal(ob[0], image):
                self._be_isabet += 1
                return ob[1], ob[2]
            if ob is not None:
                self._be_iska += 1
            return asil(image)

        net.extract_dense_map = kanca
        net._be_kanca = True

    def _be_ara(self, img_u8):
        """Kural: C'de n0(c) >= kp olan en küçük c; yoksa C[-1]. Dönüş: (i, n0(c_i), {c: n0}, ileri_sayisi)."""
        C, hedef = self.BE_MERDIVEN, self.BE_HEDEF
        n0s = {}
        sakla = self.BE_ONBELLEK and not self.BE_PROFIL
        haritalar = {}

        def n0(i):
            if i not in n0s:
                img = self._get_clahe(C[i]).apply(img_u8)
                if sakla:
                    n0s[i], h = self.n0_say(img, sakla=True)
                    # Seçilemeyecek basamağın haritası hemen bırakılır (GPU belleği): yalnız n0 >= hedef
                    # olanlar ya da merdivenin son basamağı (geri düşüş) seçilebilir.
                    if n0s[i] >= hedef or i == len(C) - 1:
                        haritalar[i] = (img, h)
                else:
                    n0s[i] = self.n0_say(img)
            return n0s[i]

        if self.BE_ARAMA == 'nedensel':
            i = min(max(self._be_onceki_i, 0), len(C) - 1)
            if n0(i) >= hedef:
                while i > 0 and n0(i - 1) >= hedef:
                    i -= 1
            else:
                while i < len(C) - 1 and n0(i) < hedef:
                    i += 1
        else:  # 'tam': artan sırada kesin arama
            i = len(C) - 1
            for j in range(len(C)):
                if n0(j) >= hedef:
                    i = j
                    break
        ileri = len(n0s)
        if self.BE_PROFIL:
            for j in range(len(C)):
                n0(j)
        self._be_onceki_i = i
        n0i = n0s[i] if i in n0s else n0(i)
        self._be_secilen = haritalar.get(i) if sakla else None
        return i, n0i, {C[j]: n0s[j] for j in sorted(n0s)}, ileri

    def _be_clahe_cift(self, rect1, rect2, adapt):
        """BE-CLAHE: kapı sol S_d ile (çift başına tek karar); açıksa aynı c sol ve sağa."""
        import time as _t
        s_d1 = 1.0 - self.calculate_shannon_entropy(rect1, roi=None)
        adapt.update(s_d1=s_d1, norm_entropy1=1.0 - s_d1, adaptive=True,
                     be_c=0.0, be_n0=float('nan'), be_ileri=0, be_ms=0.0)
        if rect2 is not None:
            s_d2 = 1.0 - self.calculate_shannon_entropy(rect2, roi=None)
            adapt.update(s_d2=s_d2, norm_entropy2=1.0 - s_d2)
        out1, out2 = rect1, rect2
        if s_d1 > self.ADAPTIVE_CLAHE_THRESHOLD_S_D:
            if self.device == 'cuda':
                torch.cuda.synchronize()
            t0 = _t.perf_counter()
            i, n0c, prof, ileri = self._be_ara(rect1)
            c = self.BE_MERDIVEN[i]
            if self.device == 'cuda':
                torch.cuda.synchronize()
            adapt.update(be_c=c, be_n0=float(n0c), be_ileri=ileri,
                         be_ms=1000.0 * (_t.perf_counter() - t0), clahe_clip1=c)
            if self.BE_PROFIL:
                for cc, nn in prof.items():
                    adapt[f'be_n0_c{cc:g}'] = float(nn)
            if self._be_secilen is not None:
                self._harita_kancasi()
                out1, self._be_onbellek = self._be_secilen
                self._be_secilen = None
            else:
                out1 = self._get_clahe(c).apply(rect1)
            if rect2 is not None:
                out2 = self._get_clahe(c).apply(rect2)
                adapt['clahe_clip2'] = c
        else:
            self._be_onceki_i = 0
        if self.ADAPTIVE_DET_ENABLED:
            adapt['det_threshold'] = self._adapt_detection_threshold(s_d1)
        print(f"  BE-CLAHE: S_d(sol)={s_d1:.3f}, c={adapt['be_c']:.2f}, n0={adapt['be_n0']}, "
              f"ileri={adapt['be_ileri']}, {adapt['be_ms']:.1f} ms"
              + (f", onb={self._be_isabet}/{self._be_iska}" if self.BE_ONBELLEK else ""))
        return out1, out2, adapt

    def filter_matches_with_fundamental(self, kp1_coords, kp2_coords, matches, threshold=1.0):
        """Eşleşmeleri RANSAC temel matrisiyle filtrele. Dönüş: (DMatch listesi, (pts1, pts2))."""
        if len(matches) < 8: # findFundamentalMat için en az 8 nokta
            return [], (np.array([]), np.array([]))

        pts1_fm = np.float32([kp1_coords[m.queryIdx] for m in matches])
        pts2_fm = np.float32([kp2_coords[m.trainIdx] for m in matches])

        if len(pts1_fm) < 8:
             return [], (np.array([]), np.array([]))

        F, mask = cv2.findFundamentalMat(pts1_fm, pts2_fm, cv2.FM_RANSAC, threshold)

        if F is None or mask is None:
            return [], (np.array([]), np.array([]))

        filtered_dmatches = []
        good_pts1_list = []
        good_pts2_list = []
        
        mask_flat = mask.ravel()
        for i, m_fm in enumerate(matches):
            if mask_flat[i]:
                filtered_dmatches.append(m_fm)
                good_pts1_list.append(pts1_fm[i])
                good_pts2_list.append(pts2_fm[i])
        
        good_pts1_arr = np.array(good_pts1_list, dtype=np.float32)
        good_pts2_arr = np.array(good_pts2_list, dtype=np.float32)

        return filtered_dmatches, (good_pts1_arr, good_pts2_arr)

    def filter_matches_y_axis(self, kp1_coords, kp2_coords, matches, y_threshold=2.0, min_disparity=0.0):
        """Rektifiye çiftte |dy| < y_threshold ve disparity > min_disparity olan eşleşmeleri tut.
        Dönüş: (DMatch listesi, (pts1, pts2))."""
        if not matches:
            return [], (np.array([]), np.array([]))

        filtered_dmatches = []
        good_pts1_list = []
        good_pts2_list = []

        for m in matches:
            pt1 = kp1_coords[m.queryIdx]
            pt2 = kp2_coords[m.trainIdx]

            disparity = pt1[0] - pt2[0]
            if abs(pt1[1] - pt2[1]) < y_threshold and disparity > min_disparity:
                filtered_dmatches.append(m)
                good_pts1_list.append(pt1)
                good_pts2_list.append(pt2)

        good_pts1_arr = np.array(good_pts1_list, dtype=np.float32)
        good_pts2_arr = np.array(good_pts2_list, dtype=np.float32)
        return filtered_dmatches, (good_pts1_arr, good_pts2_arr)

    def triangulate_points(self, pts1, pts2):
        """Rektifiye piksel çiftlerini P1/P2 ile üçgenle. Noktalar zaten rektifiye;
        undistortPoints uygulanmaz (derinliği bozar)."""
        pts1 = np.asarray(pts1, dtype=np.float32).reshape(-1, 2).T  # (2, N)
        pts2 = np.asarray(pts2, dtype=np.float32).reshape(-1, 2).T  # (2, N)

        points_4d = cv2.triangulatePoints(self.P1, self.P2, pts1, pts2)
        points_3d = points_4d[:3] / points_4d[3]

        return points_3d.T

    def process_frame_pair(self, timestamp):
        """Bir stereo çifti işle: çıkarım, sol-sağ eşleştirme, filtre, üçgenleme."""
        try:
            filename0 = self.cam0_df.loc[self.cam0_df['timestamp'] == timestamp, 'filename'].values
            filename1 = self.cam1_df.loc[self.cam1_df['timestamp'] == timestamp, 'filename'].values

            if len(filename0) == 0 or len(filename1) == 0:
                print(f"Uyarı: {timestamp} için dosya adları bulunamadı")
                return None, None, None, None

            filename0 = filename0[0]
            filename1 = filename1[0]

            img1_path = self.base_path / self.sequence / "mav0/cam0/data" / filename0
            img2_path = self.base_path / self.sequence / "mav0/cam1/data" / filename1

            img1 = cv2.imread(str(img1_path), cv2.IMREAD_GRAYSCALE)
            img2 = cv2.imread(str(img2_path), cv2.IMREAD_GRAYSCALE)

            if img1 is None or img2 is None:
                print(f"Uyarı: {timestamp} için görüntüler yüklenemedi")
                return None, None, None, None

            rect1, rect2 = self.rectify_images(img1, img2)

            # [TÜRBÜLANS] rektifiye düzlemde, BE-CLAHE/ön işlemeden ÖNCE (atmosfer görüntü oluşumuna aittir)
            if self.turb is not None:
                rect1, rect2 = self.turb.apply_stereo(rect1, rect2, self._seq_frame_idx)

            # Adaptif ön işleme; rect1/rect2 bundan sonra işlenmiş görüntü (geometri değişmez)
            rect1, rect2, adapt_stats = self._preprocess_stereo_pair(rect1, rect2)
            self._kare_kapi_kaydet(timestamp, adapt_stats)

            # Özellik çıkarımı
            img1_tensor = frame2tensor(rect1, self.device)
            img2_tensor = frame2tensor(rect2, self.device)
            h0, w0 = img1_tensor.shape[-2:]
            h1, w1 = img2_tensor.shape[-2:]

            # Önceki karenin görselleştirme verisi
            actual_prev_frame_image_for_vis = self.prev_vis_info['image']
            actual_prev_frame_candidate_kps_for_vis = self.prev_vis_info['kp_candidates']

            image_size0_aliked = img1_tensor.new_tensor([[w0, h0]])
            image_size1_aliked = img2_tensor.new_tensor([[w1, h1]])
            data0_for_aliked = {'image': img1_tensor, 'image_size': image_size0_aliked}
            data1_for_aliked = {'image': img2_tensor, 'image_size': image_size1_aliked}

            with torch.no_grad():
                feats0_explicit = self.aliked(data0_for_aliked)
                feats1_explicit = self.aliked(data1_for_aliked)

            kp1_aliked_tensor = feats0_explicit['keypoints'][0]
            desc1_aliked_tensor = feats0_explicit['descriptors'][0]
            scores1_aliked_tensor = feats0_explicit.get('keypoint_scores')[0] if 'keypoint_scores' in feats0_explicit and feats0_explicit['keypoint_scores'] is not None else None
            
            kp2_aliked_tensor = feats1_explicit['keypoints'][0]
            desc2_aliked_tensor = feats1_explicit['descriptors'][0]
            scores2_aliked_tensor = feats1_explicit.get('keypoint_scores')[0] if 'keypoint_scores' in feats1_explicit and feats1_explicit['keypoint_scores'] is not None else None

            # LightGlue-SIFT 'scales'/'oris' ister; ALIKED/SuperPoint'te None kalır.
            def _so(feats, ad):
                v = feats.get(ad)
                return v[0] if v is not None else None
            sc1_tensor, or1_tensor = _so(feats0_explicit, 'scales'), _so(feats0_explicit, 'oris')
            sc2_tensor, or2_tensor = _so(feats1_explicit, 'scales'), _so(feats1_explicit, 'oris')

            if kp1_aliked_tensor.shape[0] == 0 or kp2_aliked_tensor.shape[0] == 0:
                print(f"Uyarı: {timestamp} için temel özellik tespiti (self.aliked ile) başarısız.")
                return None, None, None, None
            
            kp1_all_coords = kp1_aliked_tensor.cpu().numpy()
            desc1_all = desc1_aliked_tensor.cpu().numpy()
            kp2_all_coords = kp2_aliked_tensor.cpu().numpy()
            desc2_all = desc2_aliked_tensor.cpu().numpy()

            initial_matches_for_processing = [] # sol->sağ DMatch listesi

            # Sol->sağ eşleştirme
            if self.matcher_type == 'lightglue':
                # LightGlue image_size'i (w, h) bekler (lightglue.py normalize_keypoints; ALIKED de [w, h]
                # verir). (h, w) verilirse konum kodlamasi kayar; LG_IMAGE_SIZE_WH ile duzeltilir.
                if getattr(self, 'LG_IMAGE_SIZE_WH', False):
                    hw0_tensor_lg = img1_tensor.new_tensor([w0, h0])
                    hw1_tensor_lg = img2_tensor.new_tensor([w1, h1])
                else:
                    hw0_tensor_lg = img1_tensor.new_tensor([h0, w0])
                    hw1_tensor_lg = img2_tensor.new_tensor([h1, w1])
                lg_input_data = {
                    'image0': {'keypoints': kp1_aliked_tensor.unsqueeze(0), 'descriptors': desc1_aliked_tensor.unsqueeze(0), 'image_size': hw0_tensor_lg.unsqueeze(0)},
                    'image1': {'keypoints': kp2_aliked_tensor.unsqueeze(0), 'descriptors': desc2_aliked_tensor.unsqueeze(0), 'image_size': hw1_tensor_lg.unsqueeze(0)}
                }
                if scores1_aliked_tensor is not None: lg_input_data['image0']['scores'] = scores1_aliked_tensor.unsqueeze(0)
                if scores2_aliked_tensor is not None: lg_input_data['image1']['scores'] = scores2_aliked_tensor.unsqueeze(0)

                if self.lightglue_stereo.conf.add_scale_ori: # sift için True
                    (lg_input_data['image0']['scales'],
                     lg_input_data['image0']['oris']) = self._lg_scale_ori(
                        kp1_aliked_tensor.shape[0], sc1_tensor, or1_tensor, self.device)
                    (lg_input_data['image1']['scales'],
                     lg_input_data['image1']['oris']) = self._lg_scale_ori(
                        kp2_aliked_tensor.shape[0], sc2_tensor, or2_tensor, self.device)
                
                try:
                    with torch.no_grad():
                        res_lightglue = self.lightglue_stereo(lg_input_data)
                    matches0_lg_indices = res_lightglue['matches0'][0].cpu()
                    match_scores0 = res_lightglue.get('matching_scores0', [torch.zeros_like(matches0_lg_indices, dtype=torch.float)])[0].cpu()

                    for i0 in range(kp1_aliked_tensor.shape[0]):
                        i1 = matches0_lg_indices[i0].item()
                        if i1 != -1:
                            score_val = match_scores0[i0].item() if match_scores0.numel() > i0 else 0.0
                            initial_matches_for_processing.append(cv2.DMatch(_queryIdx=i0, _trainIdx=i1, _distance=1.0 - score_val))
                except KeyError as e:
                     print(f"LightGlue eşleştirme hatası (muhtemelen 'scales' eksik ve handle edilmedi): {e}")
                     traceback.print_exc()

            elif self.matcher_type == 'classical_bf':
                if desc1_all is not None and desc2_all is not None and len(desc1_all) > 0 and len(desc2_all) > 0:
                    desc1_bf = np.ascontiguousarray(desc1_all, dtype=np.float32)
                    desc2_bf = np.ascontiguousarray(desc2_all, dtype=np.float32)
                    
                    if hasattr(self, 'stereo_bf_matcher') and self.stereo_bf_matcher is not None:
                        raw_matches_bf = self.stereo_bf_matcher.knnMatch(desc1_bf, desc2_bf, k=2)
                        for match_pair in raw_matches_bf:
                            if len(match_pair) == 2:
                                m, n = match_pair
                                if m.distance < 0.8 * n.distance: # Lowe oran testi
                                    initial_matches_for_processing.append(m)
                    else:
                        print("Uyarı: Klasik BFMatcher (stereo) başlatılmamış.")
            
            if not initial_matches_for_processing:
                print(f"Uyarı: {timestamp} için {self.matcher_type} ile Sol->Sağ stereo eşleşme bulunamadı.")
                self.prev_vis_info.update({
                    'image': rect1.copy() if rect1 is not None else None,
                    'kp_candidates': None,
                    'filename0': filename0
                })
                return None, None, None, None

            # Filtreleme / çapraz kontrol
            good_stereo_matches_f = []
            filtered_points_stereo = (np.array([]), np.array([]))

            if self.filter_method == 'cross_check':
                print("Çapraz Kontrol (Cross-Check) filtresi uygulanıyor...")
                matches_RL_dmatch_temp = [] # sağ->sol eşleşmeleri

                if self.matcher_type == 'lightglue':
                    hw1_tensor_lg_rev = img2_tensor.new_tensor([h1, w1])
                    hw0_tensor_lg_rev = img1_tensor.new_tensor([h0, w0])
                    lg_input_data_rev = { # image0: sağ, image1: sol
                        'image0': {'keypoints': kp2_aliked_tensor.unsqueeze(0), 'descriptors': desc2_aliked_tensor.unsqueeze(0), 'image_size': hw1_tensor_lg_rev.unsqueeze(0)},
                        'image1': {'keypoints': kp1_aliked_tensor.unsqueeze(0), 'descriptors': desc1_aliked_tensor.unsqueeze(0), 'image_size': hw0_tensor_lg_rev.unsqueeze(0)}
                    }
                    if scores2_aliked_tensor is not None: lg_input_data_rev['image0']['scores'] = scores2_aliked_tensor.unsqueeze(0)
                    if scores1_aliked_tensor is not None: lg_input_data_rev['image1']['scores'] = scores1_aliked_tensor.unsqueeze(0)
                    if self.lightglue_stereo.conf.add_scale_ori:
                        (lg_input_data_rev['image0']['scales'],
                         lg_input_data_rev['image0']['oris']) = self._lg_scale_ori(
                            kp2_aliked_tensor.shape[0], sc2_tensor, or2_tensor, self.device)
                        (lg_input_data_rev['image1']['scales'],
                         lg_input_data_rev['image1']['oris']) = self._lg_scale_ori(
                            kp1_aliked_tensor.shape[0], sc1_tensor, or1_tensor, self.device)
                    try:
                        with torch.no_grad(): res_lightglue_rev = self.lightglue_stereo(lg_input_data_rev)
                        matches0_lg_indices_rev = res_lightglue_rev['matches0'][0].cpu() # kp2_idx -> kp1_idx
                        match_scores0_rev = res_lightglue_rev.get('matching_scores0', [torch.zeros_like(matches0_lg_indices_rev, dtype=torch.float)])[0].cpu()
                        for i2_rev in range(kp2_aliked_tensor.shape[0]):
                            i1_rev = matches0_lg_indices_rev[i2_rev].item()
                            if i1_rev != -1:
                                score_val_rev = match_scores0_rev[i2_rev].item() if match_scores0_rev.numel() > i2_rev else 0.0
                                matches_RL_dmatch_temp.append(cv2.DMatch(_queryIdx=i2_rev, _trainIdx=i1_rev, _distance=1.0 - score_val_rev))
                    except KeyError as e: print(f"LightGlue Sağ->Sol eşleştirme hatası: {e}"); traceback.print_exc()

                elif self.matcher_type == 'classical_bf':
                    if desc2_all is not None and desc1_all is not None and len(desc2_all) > 0 and len(desc1_all) > 0:
                        desc2_bf_rev = np.ascontiguousarray(desc2_all, dtype=np.float32)
                        desc1_bf_rev = np.ascontiguousarray(desc1_all, dtype=np.float32)
                        if hasattr(self, 'stereo_bf_matcher') and self.stereo_bf_matcher is not None:
                            raw_matches_bf_rev = self.stereo_bf_matcher.knnMatch(desc2_bf_rev, desc1_bf_rev, k=2)
                            for match_pair_rev in raw_matches_bf_rev:
                                if len(match_pair_rev) == 2:
                                    m_rev, n_rev = match_pair_rev
                                    if m_rev.distance < 0.8 * n_rev.distance:
                                        matches_RL_dmatch_temp.append(m_rev)

                # İki yönde tutarlı eşleşmeleri tut
                lookup_matches_RL = {m.queryIdx: m.trainIdx for m in matches_RL_dmatch_temp}
                for m_LR in initial_matches_for_processing:
                    kp1_idx_fwd = m_LR.queryIdx
                    kp2_idx_fwd = m_LR.trainIdx
                    if kp2_idx_fwd in lookup_matches_RL and lookup_matches_RL[kp2_idx_fwd] == kp1_idx_fwd:
                        good_stereo_matches_f.append(m_LR)

                if good_stereo_matches_f:
                    pts1_nofilter = np.float32([kp1_all_coords[m.queryIdx] for m in good_stereo_matches_f])
                    pts2_nofilter = np.float32([kp2_all_coords[m.trainIdx] for m in good_stereo_matches_f])
                    filtered_points_stereo = (pts1_nofilter, pts2_nofilter)

            else:
                if self.filter_method == 'y_axis_check' and initial_matches_for_processing:
                    print(f"--- {timestamp} için Y-Ekseni Filtresi Öncesi Y Farkları ---")
                    y_diffs = []
                    for m_debug in initial_matches_for_processing:
                        pt1_debug = kp1_all_coords[m_debug.queryIdx]
                        pt2_debug = kp2_all_coords[m_debug.trainIdx]
                        y_diffs.append(abs(pt1_debug[1] - pt2_debug[1]))
                    if y_diffs:
                        print(f"Min Y Farkı: {min(y_diffs):.2f}, Max Y Farkı: {max(y_diffs):.2f}, Ortalama Y Farkı: {np.mean(y_diffs):.2f}, Medyan Y Farkı: {np.median(y_diffs):.2f}")

                if self.filter_method == 'fundamental_matrix':
                    good_stereo_matches_f, filtered_points_stereo = self.filter_matches_with_fundamental(
                        kp1_all_coords, kp2_all_coords, initial_matches_for_processing
                    )
                elif self.filter_method == 'y_axis_check':
                    good_stereo_matches_f, filtered_points_stereo = self.filter_matches_y_axis(
                        kp1_all_coords, kp2_all_coords, initial_matches_for_processing, y_threshold=self.y_axis_threshold, min_disparity=self.MIN_STEREO_DISPARITY
                    )
                elif self.filter_method == 'none':
                    good_stereo_matches_f = initial_matches_for_processing
                    if initial_matches_for_processing:
                        pts1_nofilter = np.float32([kp1_all_coords[m.queryIdx] for m in initial_matches_for_processing])
                        pts2_nofilter = np.float32([kp2_all_coords[m.trainIdx] for m in initial_matches_for_processing])
                        filtered_points_stereo = (pts1_nofilter, pts2_nofilter)
                else:
                    raise ValueError(f"Geçersiz filtreleme yöntemi (ikincil filtreleme bloğu): {self.filter_method}")

            if not good_stereo_matches_f or (filtered_points_stereo[0] is not None and filtered_points_stereo[0].size == 0):
                filter_name_for_msg = "cross_check" if self.filter_method == "cross_check" else self.filter_method
                print(f"Uyarı: {timestamp} için '{filter_name_for_msg}' filtrelemesi sonrası geçerli stereo eşleşme bulunamadı")
                self.prev_vis_info.update({
                    'image': rect1.copy() if rect1 is not None else None,
                    'kp_candidates': None, 
                    'filename0': filename0
                })
                return None, None, None, None

            pts1_final, pts2_final = filtered_points_stereo

            desc1_final = np.array([desc1_all[m.queryIdx] for m in good_stereo_matches_f])

            points_3d = self.triangulate_points(pts1_final, pts2_final)

            # Geçerlilik: Z>0 ve makul derinlik; elenenler tüm hizalı dizilerden çıkarılır
            if points_3d.shape[0] > 0:
                depths = points_3d[:, 2]
                valid_depth_mask = (depths > max(0.0, self.MIN_LANDMARK_DEPTH)) & (depths < self.MAX_LANDMARK_DEPTH)
                # Parallaks kapısı: iki bakış ışını arası açı çok küçükse nokta uzak, derinlik güvenilmez
                if self.PARALLAX_CHECK:
                    b = float(np.linalg.norm(self.T))
                    C1 = np.array([b, 0.0, 0.0])
                    ray0 = points_3d
                    ray1 = points_3d - C1
                    n0 = np.linalg.norm(ray0, axis=1)
                    n1 = np.linalg.norm(ray1, axis=1)
                    denom = np.maximum(n0 * n1, 1e-12)
                    cos_par = np.sum(ray0 * ray1, axis=1) / denom
                    valid_depth_mask &= (cos_par <= self.COS_MAX_PARALLAX)
                # Stereo reprojeksiyon kontrolü: X'i P1/P2 ile geri yansıt, hata büyükse ele
                if self.STEREO_REPROJ_CHECK and np.any(valid_depth_mask):
                    Xh = np.hstack([points_3d, np.ones((points_3d.shape[0], 1))])
                    def _reproj_err(P, pts):
                        q = Xh @ np.asarray(P, dtype=np.float64).T
                        w = np.where(np.abs(q[:, 2]) > 1e-12, q[:, 2], 1e-12)
                        return np.linalg.norm(q[:, :2] / w[:, None] - pts, axis=1)
                    e1 = _reproj_err(self.P1, pts1_final)
                    e2 = _reproj_err(self.P2, pts2_final)
                    reproj_ok = (e1 <= self.STEREO_REPROJ_MAX_PX) & (e2 <= self.STEREO_REPROJ_MAX_PX)
                    self.cull_stats['stereo_reproj_rejected'] += int(np.sum(valid_depth_mask & ~reproj_ok))
                    # medyan artık (px), eşik seçimi için teşhis
                    self.cull_stats['stereo_reproj_med_px_sum'] += float(np.median(np.maximum(e1, e2)[valid_depth_mask]))
                    self.cull_stats['stereo_reproj_frames'] += 1
                    valid_depth_mask &= reproj_ok
                # Göreli derinlik belirsizliği: σ_Z/Z = Z·σ_d/(f·b)
                if self.MAX_DEPTH_SIGMA_REL > 0 and np.any(valid_depth_mask):
                    fb = float(self.P1[0, 0]) * float(np.linalg.norm(self.T))
                    sig_rel = np.abs(depths) * self.STEREO_DISP_SIGMA_PX / max(fb, 1e-9)
                    sig_ok = sig_rel <= self.MAX_DEPTH_SIGMA_REL
                    self.cull_stats['depth_sigma_rejected'] += int(np.sum(valid_depth_mask & ~sig_ok))
                    valid_depth_mask &= sig_ok
                if not np.all(valid_depth_mask):
                    pts1_final = pts1_final[valid_depth_mask]
                    pts2_final = pts2_final[valid_depth_mask]
                    desc1_final = desc1_final[valid_depth_mask]
                    points_3d = points_3d[valid_depth_mask]
                    good_stereo_matches_f = [m for m, keep in zip(good_stereo_matches_f, valid_depth_mask) if keep]

            # Temporal takip
            vis_prev_kps_continued_np, \
            vis_curr_kps_continued_np, \
            vis_newly_added_kps_curr_np, \
            motion_vectors_t, \
            num_continued_tracks = self.update_and_select_tracks(
                pts1_final,
                desc1_final,
                points_3d,
                kp1_all_coords,
                desc1_all
            )

            result = {
                'timestamp': timestamp,
                'kp1_count': len(kp1_all_coords),
                'kp2_count': len(kp2_all_coords),
                'stereo_matches': len(initial_matches_for_processing),  # filtre öncesi sol->sağ
                'current_stereo_triangulated_points': len(pts1_final),
                'temporal_matches': num_continued_tracks,
                'active_depth_tracks': len(self.get_active_landmarks(only_with_depth=True)),
                'ran_stereo': True
            }
            result.update(adapt_stats)

            # Görselleştirme verisi
            actual_prev_frame_image_for_vis = self.prev_vis_info['image']
            actual_prev_frame_candidate_kps_for_vis = self.prev_vis_info['kp_candidates']

            vis_data_temporal_tuple = None
            if actual_prev_frame_image_for_vis is not None and \
               vis_prev_kps_continued_np is not None and \
               vis_curr_kps_continued_np is not None and \
               vis_newly_added_kps_curr_np is not None:
                vis_data_temporal_tuple = (
                    actual_prev_frame_image_for_vis,
                    actual_prev_frame_candidate_kps_for_vis,
                    rect1,
                    vis_prev_kps_continued_np,
                    vis_curr_kps_continued_np,
                    vis_newly_added_kps_curr_np,
                    motion_vectors_t
                )
            vis_data = {
                'stereo': (rect1, rect2, kp1_all_coords, kp2_all_coords, good_stereo_matches_f, pts1_final),
                'temporal': vis_data_temporal_tuple
            }

            prev_filename_for_temporal_panel = self.prev_vis_info['filename0']

            self.prev_vis_info['image'] = rect1.copy() if rect1 is not None else None
            self.prev_vis_info['kp_candidates'] = pts1_final.copy() if pts1_final is not None and pts1_final.size > 0 else None
            self.prev_vis_info['filename0'] = filename0

            filenames_stereo = (filename0, filename1)
            filenames_temporal = (prev_filename_for_temporal_panel, filename0)
            
            matcher_name_print = "LightGlue" if self.matcher_type == "lightglue" else "BFMatcher"
            filter_name_print = self.filter_method if self.filter_method != "fundamental_matrix" else "F-Matrix"
            print(f"İşleme başarılı - {matcher_name_print} Eşleşme / {filter_name_print} Filtreli: {len(initial_matches_for_processing)}/{len(pts1_final)}")
            print(f"Temporal eşleşme sayısı (devam eden): {num_continued_tracks}")

            return (result, vis_data, filenames_stereo, filenames_temporal)

        except Exception as e:
            if _cuda_oom(e):
                raise
            print(f"Kare çifti işleme hatası: {e}")
            traceback.print_exc()
            if 'desc1_all' in locals() and desc1_all is not None:
                print(f"desc1_all (tanımlayıcılar) tanımlı, şekli: {desc1_all.shape}")
            elif 'desc1_aliked_tensor' in locals() and desc1_aliked_tensor is not None:
                print("desc1_all tanımsızdı, ancak desc1_aliked_tensor mevcut.")
            else:
                print("Ne desc1_all ne de desc1_aliked_tensor tanımlı (veya None).")

            return None, None, None, None

    def process_frame_temporal_only(self, timestamp):
        """Sadece sol görüntüyle temporal takip (stereo atlanır)."""
        try:
            filename0_series = self.cam0_df.loc[self.cam0_df['timestamp'] == timestamp, 'filename']
            if filename0_series.empty:
                print(f"Uyarı: {timestamp} için cam0 dosya adı bulunamadı")
                return None, None, None, None
            filename0 = filename0_series.iloc[0]

            img1_path = self.base_path / self.sequence / "mav0/cam0/data" / filename0
            img1 = cv2.imread(str(img1_path), cv2.IMREAD_GRAYSCALE)

            if img1 is None:
                print(f"Uyarı: {timestamp} için cam0 görüntüsü yüklenemedi")
                return None, None, None, None

            rect1 = cv2.remap(img1, self.map1x, self.map1y, cv2.INTER_LINEAR)

            # [TÜRBÜLANS] stereo yolundakiyle aynı kare indeksi ve aynı gerçekleme (yalnız cam0)
            if self.turb is not None:
                rect1, _ = self.turb.apply_stereo(rect1, None, self._seq_frame_idx)

            # Aynı ön işleme burada da uygulanır; yoksa ardışık karelerin descriptor uzayı ayrışır
            rect1, _, adapt_stats = self._preprocess_stereo_pair(rect1, None)
            self._kare_kapi_kaydet(timestamp, adapt_stats)

            img1_tensor = frame2tensor(rect1, self.device)
            h0, w0 = img1_tensor.shape[-2:]
            image_size0_aliked = img1_tensor.new_tensor([[w0, h0]])
            data0_for_aliked = {'image': img1_tensor, 'image_size': image_size0_aliked}

            with torch.no_grad():
                feats0_explicit = self.aliked(data0_for_aliked)

            kp1_all_coords_tensor = feats0_explicit['keypoints'][0]
            desc1_all_tensor = feats0_explicit['descriptors'][0]

            if kp1_all_coords_tensor.shape[0] == 0:
                print(f"Uyarı: {timestamp} için cam0'dan temel özellik tespiti (self.aliked ile) başarısız.")
                self.prev_vis_info['image'] = rect1.copy() if rect1 is not None else None
                self.prev_vis_info['kp_candidates'] = np.array([])
                self.prev_vis_info['filename0'] = filename0
                return None, None, None, None

            kp1_all_coords = kp1_all_coords_tensor.cpu().numpy()
            desc1_all = desc1_all_tensor.cpu().numpy()

            # Stereo verisi olmadan track güncelle
            vis_prev_kps_continued_np, \
            vis_curr_kps_continued_np, \
            vis_newly_added_kps_curr_np, \
            motion_vectors_t, \
            num_continued_tracks = self.update_and_select_tracks(
                current_stereo_good_kps=np.array([]),
                current_stereo_good_descs=np.array([]),
                current_stereo_good_p3ds=np.array([]),
                all_kps_in_current_rect1=kp1_all_coords,
                all_descs_in_current_rect1=desc1_all
            )

            result = {
                'timestamp': timestamp,
                'kp1_count': len(kp1_all_coords),
                'kp2_count': 0,
                'stereo_matches': 0,
                'current_stereo_triangulated_points': 0,
                'temporal_matches': num_continued_tracks,
                'active_depth_tracks': len(self.get_active_landmarks(only_with_depth=True)),
                'ran_stereo': False
            }
            result.update(adapt_stats)

            actual_prev_frame_image_for_vis = self.prev_vis_info['image']
            actual_prev_frame_candidate_kps_for_vis = self.prev_vis_info['kp_candidates']

            vis_data_stereo_tuple = (rect1, None, kp1_all_coords, None, [], np.array([]))
            vis_data_temporal_tuple = (
                actual_prev_frame_image_for_vis, actual_prev_frame_candidate_kps_for_vis,
                rect1, vis_prev_kps_continued_np, vis_curr_kps_continued_np,
                vis_newly_added_kps_curr_np,
                motion_vectors_t
            )
            vis_data = {'stereo': vis_data_stereo_tuple, 'temporal': vis_data_temporal_tuple}

            prev_filename_for_temporal_panel = self.prev_vis_info['filename0']
            self.prev_vis_info.update({'image': rect1.copy(), 'kp_candidates': np.array([]), 'filename0': filename0})
            
            filenames_stereo = (filename0, "N/A")
            filenames_temporal = (prev_filename_for_temporal_panel, filename0)
            print(f"Sadece temporal takip başarılı - Devam eden: {num_continued_tracks}")
            return result, vis_data, filenames_stereo, filenames_temporal

        except Exception as e:
            if _cuda_oom(e):
                raise
            print(f"Sadece temporal kare işleme hatası: {e}")
            traceback.print_exc()
            return None, None, None, None

    def update_and_select_tracks(self,
                                 current_stereo_good_kps,    # pts1_final
                                 current_stereo_good_descs,  # desc1_final
                                 current_stereo_good_p3ds,   # points_3d
                                 all_kps_in_current_rect1,   # rect1'deki tüm kp'ler
                                 all_descs_in_current_rect1  # rect1'deki tüm desc'ler
                                ):
        """Track'leri eşleştir, güncelle, eskileri sil, yenilerini ekle. Stereo'lu ve mono
        karelerde çağrılır; has_depth bayrağı stereo derinliği olan track'i işaretler."""
        self.current_frame_idx += 1

        # Keypoint dışa aktarımı (backend'de search-by-projection / loop closure için)
        if getattr(self, 'EXPORT_KEYPOINTS', False):
            if all_kps_in_current_rect1 is not None and len(all_kps_in_current_rect1) > 0:
                self._cur_all_kps = np.asarray(all_kps_in_current_rect1, dtype=np.float32).reshape(-1, 2)
                self._cur_all_descs = np.asarray(all_descs_in_current_rect1, dtype=np.float32)
            else:
                self._cur_all_kps = np.zeros((0, 2), dtype=np.float32)
                self._cur_all_descs = None

        continued_tracks_info_list = []
        vis_prev_kps_continued = []     # devam eden track'lerin önceki karedeki kp'leri
        vis_curr_kps_continued = []     # devam eden track'lerin bu karedeki kp'leri
        culled_ids = set()              # bu karede kalıcı silinecek track id'leri

        # Stereo noktalarına derinlik ataması için KDTree
        has_stereo = current_stereo_good_kps is not None and len(current_stereo_good_kps) > 0
        stereo_kdtree = cKDTree(np.asarray(current_stereo_good_kps, dtype=np.float32)) if has_stereo else None

        # 1. Aktif track'leri eşleştir (MAX_TRACK_AGE_UNSEEN içinde kaybolanlar da dahil)
        active_tracks_from_last_frame = [
            track for track in self.tracked_features
            if 1 <= (self.current_frame_idx - track['last_seen_frame_idx']) <= self.MAX_TRACK_AGE_UNSEEN
        ]

        if active_tracks_from_last_frame and all_descs_in_current_rect1 is not None and len(all_descs_in_current_rect1) > 0:
            prev_track_descs_np = np.array([t['desc'] for t in active_tracks_from_last_frame], dtype=np.float32)
            curr_all_descs_np = np.ascontiguousarray(all_descs_in_current_rect1, dtype=np.float32)

            if prev_track_descs_np.size > 0:
                candidate_matches = []
                bf_query = list(range(len(active_tracks_from_last_frame)))
                lg_bu_kare = (self.TEMPORAL_MATCHER == 'lightglue' or (self.TEMPORAL_LG_KAPILI and self._kapi_acik)
                              or (self.TEMPORAL_LG_SAGLIK > 0 and self._son_devam_orani < self.TEMPORAL_LG_SAGLIK))
                if lg_bu_kare and self._prev_frame_feats is not None:
                    # Onceki karenin TUM ozellikleri ile bu karenin ozellikleri LightGlue ile eslenir
                    # (AirSLAM map_builder.cc:122 kareler arasi takipte de LightGlue kullanir); bir onceki
                    # karede gorulen track'in anahtar noktasi o karenin tespit kumesinde birebir bulunur.
                    pk, pdsc = self._prev_frame_feats
                    pidx = {tuple(np.asarray(k, dtype=np.float32).tolist()): i for i, k in enumerate(pk)}
                    lg_m0 = self._lightglue_match(pk, pdsc, all_kps_in_current_rect1, curr_all_descs_np)
                    bf_query = []
                    for qi, t in enumerate(active_tracks_from_last_frame):
                        j = pidx.get(tuple(np.asarray(t['kp'], dtype=np.float32).tolist())) \
                            if t['last_seen_frame_idx'] == self.current_frame_idx - 1 else None
                        if j is None:
                            bf_query.append(qi)          # bir kare goremeyen track: BF+Lowe yolu
                            continue
                        k = int(lg_m0[j][0])
                        if k >= 0:
                            candidate_matches.append(cv2.DMatch(_queryIdx=qi, _trainIdx=k,
                                                                _distance=float(1.0 - lg_m0[j][1])))
                    self.n_lg_temporal += 1
                if bf_query:
                    raw_matches = self.temporal_bf_matcher.knnMatch(prev_track_descs_np[bf_query], curr_all_descs_np, k=2)
                    for match_pair in raw_matches:
                        if len(match_pair) == 2:
                            m, n = match_pair
                            if m.distance < self.TEMPORAL_LOWE_RATIO * n.distance: # Lowe oran testi
                                candidate_matches.append(cv2.DMatch(_queryIdx=bf_query[m.queryIdx], _trainIdx=m.trainIdx,
                                                                    _distance=m.distance + 2.0) if lg_bu_kare else m)

                # Her keypoint en fazla bir track'e; çakışmada en düşük mesafe kazanır
                best_match_for_trainidx = {}
                for m in candidate_matches:
                    existing = best_match_for_trainidx.get(m.trainIdx)
                    if existing is None or m.distance < existing.distance:
                        best_match_for_trainidx[m.trainIdx] = m

                for m in best_match_for_trainidx.values():
                    track_idx_in_active_list = m.queryIdx
                    matched_kp_idx_in_all_current = m.trainIdx

                    original_track = active_tracks_from_last_frame[track_idx_in_active_list]

                    # Ömür sınırını aşan track emekli edilir; keypoint yeni track'e kalır
                    if self.MAX_TRACK_LENGTH > 0 and original_track['age'] >= self.MAX_TRACK_LENGTH:
                        continue

                    new_kp = all_kps_in_current_rect1[matched_kp_idx_in_all_current]
                    new_desc = all_descs_in_current_rect1[matched_kp_idx_in_all_current]

                    # Descriptor sürüklenmesi: anchor descriptor'dan kosinüs mesafesi eşiği aşarsa kalıcı sil
                    if self.CULL_DESC_DRIFT_MAX > 0:
                        d0 = original_track.get('desc0')
                        if d0 is not None:
                            den = float(np.linalg.norm(d0)) * float(np.linalg.norm(new_desc))
                            cos_dist = 1.0 - float(np.dot(d0, new_desc)) / den if den > 1e-12 else 0.0
                            if cos_dist > self.CULL_DESC_DRIFT_MAX:
                                culled_ids.add(original_track['id'])
                                self.cull_stats['desc_drift'] += 1
                                continue

                    # Derinlik durumunu önceki track'ten devral
                    new_p3d = original_track['p3d']
                    new_has_depth = original_track.get('has_depth', original_track['p3d'] is not None)
                    init_kp = original_track['initial_stereo_kp']
                    init_desc = original_track['initial_stereo_desc']
                    new_depth_frame_idx = original_track.get('depth_frame_idx', -1)  # p3d hangi karenin çerçevesinde

                    # Bu karede stereo varsa yakın stereo noktasından derinliği tazele
                    if stereo_kdtree is not None:
                        dist, stereo_idx = stereo_kdtree.query(new_kp)
                        if dist < self.STEREO_DEPTH_ASSOC_PX: # px
                            new_p3d = current_stereo_good_p3ds[stereo_idx]
                            new_depth_frame_idx = self.current_frame_idx
                            if not new_has_depth:
                                # ilk derinlik: stereo anchor'ını kaydet
                                init_kp = new_kp.copy()
                                init_desc = current_stereo_good_descs[stereo_idx].copy()
                            new_has_depth = True

                    p3d_before = original_track['p3d'].copy() if original_track['p3d'] is not None else None
                    continued_tracks_info_list.append({
                        'id': original_track['id'],
                        'p3d_before_update': p3d_before,
                        'kp': new_kp, 'desc': new_desc, 'p3d': new_p3d, 'has_depth': new_has_depth,
                        'depth_frame_idx': new_depth_frame_idx,
                        'age': original_track['age'] + 1,
                        'last_seen_frame_idx': self.current_frame_idx,
                        'initial_stereo_kp': init_kp,
                        'initial_stereo_desc': init_desc,
                        # ömür istatistikleri (eski kayıtlarda yoksa türetilir)
                        'created_frame_idx': original_track.get('created_frame_idx', self.current_frame_idx - original_track['age']),
                        'n_outlier': original_track.get('n_outlier', 0),
                        'desc0': original_track.get('desc0', original_track['desc']),
                    })
                    vis_prev_kps_continued.append(original_track['kp'])
                    vis_curr_kps_continued.append(new_kp)

        # 1b. Temporal RANSAC: devam eden track'lerin hareketine F-matris uydur, aykırıları ele
        if len(continued_tracks_info_list) >= 8:
            prev_pts_ransac = np.ascontiguousarray(vis_prev_kps_continued, dtype=np.float32).reshape(-1, 1, 2)
            curr_pts_ransac = np.ascontiguousarray(vis_curr_kps_continued, dtype=np.float32).reshape(-1, 1, 2)
            ransac_mask = None
            # Dejenere geometride OpenCV iç hata verebilir: önce MAGSAC, sonra FM_RANSAC, olmazsa atla
            for _method in (getattr(cv2, 'USAC_MAGSAC', cv2.FM_RANSAC), cv2.FM_RANSAC):
                try:
                    _F, ransac_mask = cv2.findFundamentalMat(
                        prev_pts_ransac, curr_pts_ransac, _method,
                        self.TEMPORAL_RANSAC_THRESHOLD, self.TEMPORAL_RANSAC_CONF, self.TEMPORAL_RANSAC_ITERS
                    )
                    if ransac_mask is not None:
                        break
                except cv2.error:
                    ransac_mask = None
            if ransac_mask is None:
                print("  Temporal RANSAC atlandı (dejenere geometri).")
            if ransac_mask is not None:
                keep = ransac_mask.ravel().astype(bool)
                if keep.any(): # hepsi elenirse filtreleme yapma
                    num_rejected = int((~keep).sum())
                    # Dışlanan track'in strike sayacı artar; eşiğe ulaşan kalıcı silinir
                    if self.CULL_RANSAC_STRIKES > 0 and num_rejected > 0:
                        by_id = {t['id']: t for t in self.tracked_features}
                        for ct, k in zip(continued_tracks_info_list, keep):
                            if k:
                                continue
                            old = by_id.get(ct['id'])
                            if old is None:
                                continue
                            old['n_outlier'] = old.get('n_outlier', 0) + 1
                            self.cull_stats['ransac_strikes'] += 1
                            if old['n_outlier'] >= self.CULL_RANSAC_STRIKES:
                                culled_ids.add(ct['id'])
                                self.cull_stats['ransac_culled'] += 1
                    continued_tracks_info_list = [ct for ct, k in zip(continued_tracks_info_list, keep) if k]
                    vis_prev_kps_continued = [p for p, k in zip(vis_prev_kps_continued, keep) if k]
                    vis_curr_kps_continued = [p for p, k in zip(vis_curr_kps_continued, keep) if k]
                    if num_rejected > 0:
                        print(f"  Temporal RANSAC: {num_rejected} outlier elendi, {len(continued_tracks_info_list)} tutuldu.")

        if self.TEMPORAL_LG_SAGLIK > 0:
            _onceki = sum(1 for t in self.tracked_features if t['last_seen_frame_idx'] == self.current_frame_idx - 1)
            self._son_devam_orani = (len(continued_tracks_info_list) / _onceki) if _onceki > 0 else 1.0
        # 2. tracked_features = devam edenler + henüz eskimemiş kayıplar
        new_tracked_features_list = []
        for updated_track_dict in continued_tracks_info_list:
            new_tracked_features_list.append(updated_track_dict)

        updated_ids = {t['id'] for t in new_tracked_features_list}
        for old_track in self.tracked_features:
            if old_track['id'] not in updated_ids:
                if (self.current_frame_idx - old_track['last_seen_frame_idx']) < self.MAX_TRACK_AGE_UNSEEN:
                    new_tracked_features_list.append(old_track)

        # Titreme: ömrü yeterli track'in eşleşme oranı (age/ömür) düşükse kalıcı sil
        if self.CULL_FLICKER_RATIO > 0:
            for t in new_tracked_features_list:
                created = t.get('created_frame_idx')
                if created is None:
                    continue
                life = self.current_frame_idx - created + 1
                if life >= self.CULL_FLICKER_MIN_LIFE and (t['age'] / float(life)) < self.CULL_FLICKER_RATIO:
                    culled_ids.add(t['id'])
                    self.cull_stats['flicker_culled'] += 1

        # Kalıcı silme; vis listeleri de tutarlı tutulur
        if culled_ids:
            new_tracked_features_list = [t for t in new_tracked_features_list if t['id'] not in culled_ids]
            if continued_tracks_info_list:
                keep_ct = [ct['id'] not in culled_ids for ct in continued_tracks_info_list]
                continued_tracks_info_list = [ct for ct, k in zip(continued_tracks_info_list, keep_ct) if k]
                vis_prev_kps_continued = [p for p, k in zip(vis_prev_kps_continued, keep_ct) if k]
                vis_curr_kps_continued = [p for p, k in zip(vis_curr_kps_continued, keep_ct) if k]
            self.cull_stats['tracks_culled_total'] += len(culled_ids)
        self.tracked_features = new_tracked_features_list

        num_currently_active_tracks = len([t for t in self.tracked_features if t['last_seen_frame_idx'] == self.current_frame_idx])

        # 3. Hedefin altındaysa ızgara tabanlı doldurma: her karede, boş hücreye önce stereo sonra mono aday
        vis_newly_added_kps_curr = []

        if num_currently_active_tracks < self.TARGET_FEATURE_COUNT and \
           all_kps_in_current_rect1 is not None and len(all_kps_in_current_rect1) > 0:

            needed = self.TARGET_FEATURE_COUNT - num_currently_active_tracks
            print(f"  Aktif takip: {num_currently_active_tracks}. Hedef: {self.TARGET_FEATURE_COUNT}. İhtiyaç: {needed} nokta. (stereo={'var' if has_stereo else 'yok'})")

            # Hücre boyutu gerçek görüntü boyutundan (rektifikasyon haritası)
            img_h, img_w = self.map1x.shape[:2]
            cell_h = img_h / self.GRID_ROWS
            cell_w = img_w / self.GRID_COLS

            # Devam eden track'lerin doldurduğu hücreler ve kullanılan kp'ler
            occupied_cells = set()
            used_kp_tuples = set()
            for kp in vis_curr_kps_continued:
                occupied_cells.add((int(kp[1] // cell_h), int(kp[0] // cell_w)))
                used_kp_tuples.add(tuple(map(int, kp)))

            # Stereo (derinlikli) ve mono adaylar için ayrı bucket'lar
            stereo_buckets = [[[] for _ in range(self.GRID_COLS)] for _ in range(self.GRID_ROWS)]
            mono_buckets   = [[[] for _ in range(self.GRID_COLS)] for _ in range(self.GRID_ROWS)]

            if has_stereo:
                for i in range(len(current_stereo_good_kps)):
                    kp = current_stereo_good_kps[i]
                    if tuple(map(int, kp)) in used_kp_tuples:
                        continue
                    r, c = int(kp[1] // cell_h), int(kp[0] // cell_w)
                    if (r, c) in occupied_cells:
                        continue
                    stereo_buckets[r][c].append(i)

            for j in range(len(all_kps_in_current_rect1)):
                kp = all_kps_in_current_rect1[j]
                if tuple(map(int, kp)) in used_kp_tuples:
                    continue
                r, c = int(kp[1] // cell_h), int(kp[0] // cell_w)
                if (r, c) in occupied_cells:
                    continue
                mono_buckets[r][c].append(j)

            # Hücrelerden dengeli seç: önce stereo, yoksa mono
            added_count = 0
            added_with_depth = 0
            points_added_in_last_pass = True
            while added_count < needed and points_added_in_last_pass:
                points_added_in_last_pass = False
                for r in range(self.GRID_ROWS):
                    for c in range(self.GRID_COLS):
                        new_track_obj = None
                        if stereo_buckets[r][c]:
                            idx = stereo_buckets[r][c].pop(0)
                            kp_new = current_stereo_good_kps[idx]
                            kp_tuple = tuple(map(int, kp_new))
                            if kp_tuple in used_kp_tuples:
                                continue
                            new_track_obj = {
                                'id': self.next_track_id, 'kp': kp_new,
                                'desc': current_stereo_good_descs[idx],
                                'p3d': current_stereo_good_p3ds[idx], 'has_depth': True,
                                'depth_frame_idx': self.current_frame_idx,
                                'age': 1,
                                'last_seen_frame_idx': self.current_frame_idx,
                                'initial_stereo_kp': kp_new.copy(),
                                'initial_stereo_desc': current_stereo_good_descs[idx].copy(),
                                'created_frame_idx': self.current_frame_idx, 'n_outlier': 0,
                                'desc0': current_stereo_good_descs[idx],
                            }
                            added_with_depth += 1
                        elif mono_buckets[r][c]:
                            idx = mono_buckets[r][c].pop(0)
                            kp_new = all_kps_in_current_rect1[idx]
                            kp_tuple = tuple(map(int, kp_new))
                            if kp_tuple in used_kp_tuples:
                                continue
                            new_track_obj = {
                                'id': self.next_track_id, 'kp': kp_new,
                                'desc': all_descs_in_current_rect1[idx],
                                'p3d': None, 'has_depth': False,
                                'depth_frame_idx': -1,
                                'age': 1,
                                'last_seen_frame_idx': self.current_frame_idx,
                                'initial_stereo_kp': None,
                                'initial_stereo_desc': None,
                                'created_frame_idx': self.current_frame_idx, 'n_outlier': 0,
                                'desc0': all_descs_in_current_rect1[idx],
                            }
                        else:
                            continue

                        self.tracked_features.append(new_track_obj)
                        vis_newly_added_kps_curr.append(kp_new)
                        used_kp_tuples.add(kp_tuple)
                        self.next_track_id += 1
                        added_count += 1
                        points_added_in_last_pass = True

                        if added_count >= needed:
                            break
                    if added_count >= needed:
                        break

            print(f"  Yeni eklenen (Izgara): {added_count} -> {added_with_depth} derinlikli + {added_count - added_with_depth} mono. Aktif: {num_currently_active_tracks + added_count}")

        num_continued_for_result = len(vis_prev_kps_continued)

        # Hareket vektörleri: iki karede de derinliği olan track'ler için
        motion_vectors = None
        if continued_tracks_info_list:
            depth_pairs = [(ct['p3d_before_update'], ct['p3d']) for ct in continued_tracks_info_list
                           if ct['p3d_before_update'] is not None and ct['p3d'] is not None]
            p3d_prev_list = [a for a, b in depth_pairs]
            p3d_curr_list = [b for a, b in depth_pairs]
            if p3d_prev_list and len(p3d_prev_list) == len(p3d_curr_list):
                p3d_prev_np = np.array(p3d_prev_list)
                p3d_curr_np = np.array(p3d_curr_list)
                if p3d_prev_np.ndim == 2 and p3d_curr_np.ndim == 2 and \
                   p3d_prev_np.shape[0] > 0 and p3d_prev_np.shape == p3d_curr_np.shape:
                    motion_vectors = p3d_curr_np - p3d_prev_np

        if all_kps_in_current_rect1 is not None and len(all_kps_in_current_rect1) > 0 and \
                all_descs_in_current_rect1 is not None and len(all_descs_in_current_rect1) > 0:
            self._prev_frame_feats = (np.asarray(all_kps_in_current_rect1, dtype=np.float32).reshape(-1, 2),
                                      np.ascontiguousarray(all_descs_in_current_rect1, dtype=np.float32))
        else:
            self._prev_frame_feats = None

        return (np.array(vis_prev_kps_continued, dtype=np.float32) if vis_prev_kps_continued else np.array([]),
                np.array(vis_curr_kps_continued, dtype=np.float32) if vis_curr_kps_continued else np.array([]),
                np.array(vis_newly_added_kps_curr, dtype=np.float32) if vis_newly_added_kps_curr else np.array([]),
                motion_vectors,
                num_continued_for_result)

    def _lightglue_match(self, kps0, desc0, kps1, desc1):
        """LightGlue (cvg varsayilan esikleriyle, lightglue.py:339-341) iki tespit kumesi arasinda.
        Donus: (N0, 2) [eslesen indeks (-1 yok), eslesme skoru]."""
        h, w = self.map1x.shape[:2]
        dev = self.device
        t = lambda a: torch.from_numpy(np.ascontiguousarray(a, dtype=np.float32))[None].to(dev)
        size = torch.tensor([[w, h]], dtype=torch.float32, device=dev)
        with torch.no_grad():
            res = self.lightglue_stereo({'image0': {'keypoints': t(kps0), 'descriptors': t(desc0), 'image_size': size},
                                         'image1': {'keypoints': t(kps1), 'descriptors': t(desc1), 'image_size': size}})
        m0 = res['matches0'][0].cpu().numpy()
        sc = res.get('matching_scores0', torch.zeros_like(res['matches0'], dtype=torch.float))[0].cpu().numpy()
        return np.stack([m0.astype(np.float64), sc.astype(np.float64)], 1)

    def get_active_landmarks(self, only_with_depth=False):
        """Bu karede görülen track'leri backend için döndür: {id, kp, p3d, has_depth, fresh_depth, age}.
        only_with_depth=True ise yalnız derinlikli olanlar."""
        landmarks = []
        for t in self.tracked_features:
            if t['last_seen_frame_idx'] != self.current_frame_idx:
                continue
            has_depth = t.get('has_depth', t['p3d'] is not None)
            if only_with_depth and not has_depth:
                continue
            # fresh_depth: p3d bu karede üçgenlendi (backend yalnız o karede init eder)
            fresh_depth = has_depth and (t.get('depth_frame_idx', -1) == self.current_frame_idx)
            landmarks.append({
                'id': t['id'],
                'kp': np.asarray(t['kp'], dtype=np.float32),
                'p3d': (np.asarray(t['p3d'], dtype=np.float32) if t['p3d'] is not None else None),
                'has_depth': has_depth,
                'fresh_depth': fresh_depth,
                'age': t['age'],
            })
        return landmarks

    def _draw_text_with_background(self, image, text, baseline_origin_for_text, font, font_scale, text_color, text_thickness, bg_color, padding=5, line_type=cv2.LINE_AA):
        """Metni dolu arka plan kutusuyla çiz; origin cv2.putText baseline noktası."""
        if image is None or text is None or len(text.strip()) == 0:
            return

        (text_w, text_h_above_baseline), baseline_from_bottom = cv2.getTextSize(text, font, font_scale, text_thickness)
        
        x_bl = baseline_origin_for_text[0]
        y_bl = baseline_origin_for_text[1]

        rect_x1 = max(0, x_bl - padding)
        rect_y1 = max(0, y_bl - text_h_above_baseline - padding)
        rect_x2 = min(image.shape[1], x_bl + text_w + padding)
        rect_y2 = min(image.shape[0], y_bl + baseline_from_bottom + padding)

        if rect_x1 < rect_x2 and rect_y1 < rect_y2:
            cv2.rectangle(image, (rect_x1, rect_y1), (rect_x2, rect_y2), bg_color, cv2.FILLED)
        
        cv2.putText(image, text, baseline_origin_for_text, font, font_scale, text_color, text_thickness, line_type)

    def visualize_combined_tracking(self, stereo_data, temporal_data, filenames_stereo, filenames_temporal, frame_stats, save_path=None, display_live=False):
        """Stereo ve temporal eşleşmeleri 1x4 panelde çiz; isteğe göre kaydet / göster."""
        if stereo_data is None:
            return
        try:
            rect1, rect2, kp1_all_coords, kp2_all_coords, good_filtered_stereo_matches, current_good_stereo_pts1_coords_for_temporal_np = stereo_data

            # Renkler (BGR)
            color_all_stereo_candidates = (0, 0, 255)    # kırmızı: tüm adaylar
            color_good_stereo_matches = (0, 255, 0)      # yeşil: filtre sonrası stereo eşleşmeler
            color_match_line = (0, 255, 0)
            color_all_temporal_candidates = (0, 255, 0)
            color_temporally_tracked = (255, 162, 0)    # sarı: devam eden track'ler
            color_newly_added_temporal = (255, 255, 0)   # camgöbeği: yeni eklenenler

            h, w = rect1.shape[:2]
            combined_h = h
            combined_w = 4 * w
            vis_img = np.zeros((combined_h, combined_w, 3), dtype=np.uint8)

            overlay_alpha = 0.7

            font = cv2.FONT_HERSHEY_SIMPLEX
            info_text_font_scale = 0.7
            info_text_thickness = 1

            # 1. Temel görüntüler
            if rect1 is not None:
                vis_img[:, :w] = cv2.cvtColor(rect1, cv2.COLOR_GRAY2BGR)
            else:
                vis_img[:, :w].fill(0)
                cv2.putText(vis_img, "No Cam0 Data", (50, h // 2), font, info_text_font_scale, (255,255,255), info_text_thickness)

            if rect2 is None:
                vis_img[:, w:2*w].fill(0)
                cv2.putText(vis_img, "No Cam1 Data", (w + 50, h // 2), font, info_text_font_scale, (255,255,255), info_text_thickness)
            else:
                vis_img[:, w:2*w] = cv2.cvtColor(rect2, cv2.COLOR_GRAY2BGR)

            panel3_base_exists = False
            panel4_base_exists = False
            if temporal_data is not None and len(temporal_data) == 7 and \
               temporal_data[0] is not None and temporal_data[2] is not None:
                # temporal_data: (prev_img, prev_candidates, curr_img, prev_tracked, curr_continued, curr_new, motion_vectors)
                prev_img_t_panel_base = temporal_data[0]
                curr_img_rect1_t_panel_base = temporal_data[2]
                if prev_img_t_panel_base is not None and curr_img_rect1_t_panel_base is not None:
                    vis_img[:, 2*w:3*w] = cv2.cvtColor(prev_img_t_panel_base.astype(np.uint8), cv2.COLOR_GRAY2BGR)
                    vis_img[:, 3*w:4*w] = cv2.cvtColor(curr_img_rect1_t_panel_base.astype(np.uint8), cv2.COLOR_GRAY2BGR)
                    panel3_base_exists = True
                    panel4_base_exists = True

            # Panel başlık metni parametreleri
            title_font_scale = 1
            text_color_fg = (0, 255, 255)
            text_thickness = 4
            bg_color_text = (0, 0, 0)
            title_text_padding = 5
            top_margin_for_box = 10
            inter_line_spacing_between_boxes = 5

            # 2. Noktalar ve çizgiler tek overlay'e çizilir, sonda bir kez harmanlanır.
            # Tampon self üzerinde yeniden kullanılır.
            _ov = getattr(self, '_viz_overlay', None)
            if _ov is None or _ov.shape != vis_img.shape or _ov.dtype != vis_img.dtype:
                _ov = np.zeros_like(vis_img)
                self._viz_overlay = _ov
            else:
                _ov[:] = 0

            # Panel 1-2: tüm adaylar
            if kp1_all_coords is not None:
                for pt_coord in kp1_all_coords:
                    cv2.circle(_ov, tuple(map(int, pt_coord)), 10, color_all_stereo_candidates, -1)

            if kp2_all_coords is not None:
                for pt_coord in kp2_all_coords:
                    cv2.circle(_ov, tuple(map(int, (pt_coord[0] + w, pt_coord[1]))), 10, color_all_stereo_candidates, -1)

            # Panel 1-2: iyi stereo eşleşmeler (nokta + çizgi)
            overlay_good_stereo = _ov
            if good_filtered_stereo_matches and kp1_all_coords is not None and kp2_all_coords is not None:
                for m_vis in good_filtered_stereo_matches:
                    if m_vis.queryIdx < len(kp1_all_coords) and m_vis.trainIdx < len(kp2_all_coords):
                        pt1_coord_vis = kp1_all_coords[m_vis.queryIdx]
                        pt2_coord_vis = kp2_all_coords[m_vis.trainIdx]

                        pt1_on_vis_img = tuple(map(int, pt1_coord_vis))
                        pt2_on_vis_img = tuple(map(int, (pt2_coord_vis[0] + w, pt2_coord_vis[1])))

                        cv2.circle(overlay_good_stereo, pt1_on_vis_img, 5, color_good_stereo_matches, -1)
                        cv2.circle(overlay_good_stereo, pt2_on_vis_img, 5, color_good_stereo_matches, -1)
                        cv2.line(overlay_good_stereo, pt1_on_vis_img, pt2_on_vis_img, color_match_line, 2)

            # Panel 3-4: temporal
            if temporal_data is not None and panel3_base_exists and panel4_base_exists:
                overlay_temporal_panel3_candidates = _ov
                overlay_temporal_panel3_tracked = _ov
                overlay_temporal_panel4_candidates = _ov
                overlay_temporal_panel4_continued = _ov
                overlay_temporal_panel4_newly_added = _ov
                overlay_temporal_lines = _ov

                prev_kps_candidates_for_panel3_green = temporal_data[1]
                prev_tm_pts_t_for_panel3 = temporal_data[3]
                curr_tm_pts_t_continued_for_panel4_yellow = temporal_data[4]
                curr_tm_pts_t_newly_added_for_panel4_cyan = temporal_data[5]

                # Panel 3 (önceki kare): adaylar yeşil, takip edilenler sarı
                if prev_kps_candidates_for_panel3_green is not None:
                    for pt_coord_prev_good_stereo in prev_kps_candidates_for_panel3_green:
                        cv2.circle(overlay_temporal_panel3_candidates, tuple(map(int, (pt_coord_prev_good_stereo[0] + 2*w, pt_coord_prev_good_stereo[1]))), 10, color_all_temporal_candidates, -1)
                if prev_tm_pts_t_for_panel3 is not None and prev_tm_pts_t_for_panel3.size > 0:
                    for pt_coord_tracked_prev in prev_tm_pts_t_for_panel3:
                        cv2.circle(overlay_temporal_panel3_tracked, tuple(map(int, (pt_coord_tracked_prev[0] + 2*w, pt_coord_tracked_prev[1]))), 5, color_temporally_tracked, -1)

                # Panel 4 (bu kare): adaylar yeşil, devam edenler sarı, yeniler camgöbeği
                if current_good_stereo_pts1_coords_for_temporal_np is not None and current_good_stereo_pts1_coords_for_temporal_np.size > 0:
                    for pt_coord_curr_good_stereo in current_good_stereo_pts1_coords_for_temporal_np:
                        cv2.circle(overlay_temporal_panel4_candidates, tuple(map(int, (pt_coord_curr_good_stereo[0] + 3*w, pt_coord_curr_good_stereo[1]))), 10, color_all_temporal_candidates, -1)
                if curr_tm_pts_t_continued_for_panel4_yellow is not None and curr_tm_pts_t_continued_for_panel4_yellow.size > 0:
                    for pt_coord_tracked_curr in curr_tm_pts_t_continued_for_panel4_yellow:
                        cv2.circle(overlay_temporal_panel4_continued, tuple(map(int, (pt_coord_tracked_curr[0] + 3*w, pt_coord_tracked_curr[1]))), 5, color_temporally_tracked, -1)
                if curr_tm_pts_t_newly_added_for_panel4_cyan is not None and curr_tm_pts_t_newly_added_for_panel4_cyan.size > 0:
                    for pt_coord_newly_added in curr_tm_pts_t_newly_added_for_panel4_cyan:
                        cv2.circle(overlay_temporal_panel4_newly_added, tuple(map(int, (pt_coord_newly_added[0] + 3*w, pt_coord_newly_added[1]))), 5, color_newly_added_temporal, -1)

                # Temporal eşleşme çizgileri
                if prev_tm_pts_t_for_panel3 is not None and curr_tm_pts_t_continued_for_panel4_yellow is not None and \
                   len(prev_tm_pts_t_for_panel3) == len(curr_tm_pts_t_continued_for_panel4_yellow):
                    for idx, pt_coord_tracked_prev in enumerate(prev_tm_pts_t_for_panel3):
                            pt_prev_on_vis_img = tuple(map(int, (pt_coord_tracked_prev[0] + 2*w, pt_coord_tracked_prev[1])))
                            pt_curr_on_vis_img = tuple(map(int, (curr_tm_pts_t_continued_for_panel4_yellow[idx][0] + 3*w, curr_tm_pts_t_continued_for_panel4_yellow[idx][1])))
                            cv2.line(overlay_temporal_lines, pt_prev_on_vis_img, pt_curr_on_vis_img, color_temporally_tracked, 2)

            # Tek harmanlama, yerinde
            cv2.addWeighted(vis_img, 1.0, _ov, overlay_alpha, 0.0, dst=vis_img)

            # 3. Metinler en son
            num_good_stereo_matches_display = len(good_filtered_stereo_matches) if good_filtered_stereo_matches is not None else 0
            num_temporal_tracks_display = frame_stats.get('temporal_matches', 0) if frame_stats is not None else 0

            # Panel 1 metinleri
            text1_content = f'Cam0 Curr {filenames_stereo[0]}'
            (tw1, th1_ab), bl1_fb = cv2.getTextSize(text1_content, font, title_font_scale, text_thickness)
            baseline_y1 = top_margin_for_box + th1_ab + title_text_padding
            self._draw_text_with_background(vis_img, text1_content, (10, baseline_y1), font, title_font_scale, text_color_fg, text_thickness, bg_color_text, title_text_padding)

            text_stereo_count_content = f'Stereo Matches: {num_good_stereo_matches_display}'
            (tw_sc, th_sc_ab), bl_sc_fb = cv2.getTextSize(text_stereo_count_content, font, title_font_scale, text_thickness)
            top_of_box_stereo_count_y = baseline_y1 + bl1_fb + title_text_padding + inter_line_spacing_between_boxes
            baseline_y_stereo_count = top_of_box_stereo_count_y + th_sc_ab + title_text_padding
            self._draw_text_with_background(vis_img, text_stereo_count_content, (10, baseline_y_stereo_count), font, title_font_scale, text_color_fg, text_thickness, bg_color_text, title_text_padding)

            matcher_label_short = "LG" if self.matcher_type == "lightglue" else "BF"
            filter_label = {"fundamental_matrix": "F", "y_axis_check": "Y", "cross_check": "CC", "none": "N"}.get(self.filter_method, "")
            extractor_label_display = self.extractor_name.upper()
            text2_content = f'{extractor_label_display}+{matcher_label_short} (Filtre:{filter_label})'
            (tw2, th2_ab), bl2_fb = cv2.getTextSize(text2_content, font, title_font_scale, text_thickness)
            top_of_box2_y = baseline_y_stereo_count + bl_sc_fb + title_text_padding + inter_line_spacing_between_boxes
            baseline_y2 = top_of_box2_y + th2_ab + title_text_padding
            self._draw_text_with_background(vis_img, text2_content, (10, baseline_y2), font, title_font_scale, text_color_fg, text_thickness, bg_color_text, title_text_padding)

            # Panel 2 metni
            text3_content = f'Cam1 Curr {filenames_stereo[1] if filenames_stereo[1] else "N/A"}'
            (tw3, th3_ab), bl3_fb = cv2.getTextSize(text3_content, font, title_font_scale, text_thickness)
            baseline_y3 = top_margin_for_box + th3_ab + title_text_padding
            self._draw_text_with_background(vis_img, text3_content, (10 + w, baseline_y3), font, title_font_scale, text_color_fg, text_thickness, bg_color_text, title_text_padding)

            # Panel 3-4 metinleri
            if temporal_data is not None and panel3_base_exists and panel4_base_exists:
                filename_prev_t, filename_curr_t = filenames_temporal
                text4_content = f'Prev: {filename_prev_t if filename_prev_t else "N/A"}'
                (tw4, th4_ab), bl4_fb = cv2.getTextSize(text4_content, font, title_font_scale, text_thickness)
                baseline_y4 = top_margin_for_box + th4_ab + title_text_padding
                self._draw_text_with_background(vis_img, text4_content, (10 + 2*w, baseline_y4), font, title_font_scale, text_color_fg, text_thickness, bg_color_text, title_text_padding)

                text5_content = f'Curr: {filename_curr_t if filename_curr_t else "N/A"}'
                if text5_content:
                    (tw5, th5_ab), bl5_fb = cv2.getTextSize(text5_content, font, title_font_scale, text_thickness)
                    baseline_y5 = top_margin_for_box + th5_ab + title_text_padding
                    self._draw_text_with_background(vis_img, text5_content, (10 + 3*w, baseline_y5), font, title_font_scale, text_color_fg, text_thickness, bg_color_text, title_text_padding)

                    text_temporal_count_content = f'Temporal Tracks: {num_temporal_tracks_display}'
                    (tw_tc, th_tc_ab), bl_tc_fb = cv2.getTextSize(text_temporal_count_content, font, title_font_scale, text_thickness)
                    top_of_box_temporal_count_y = baseline_y5 + bl5_fb + title_text_padding + inter_line_spacing_between_boxes
                    baseline_y_temporal_count = top_of_box_temporal_count_y + th_tc_ab + title_text_padding
                    self._draw_text_with_background(vis_img, text_temporal_count_content, (10 + 3*w, baseline_y_temporal_count), font, title_font_scale, text_color_fg, text_thickness, bg_color_text, title_text_padding)

            if save_path:
                if vis_img is None or vis_img.size == 0:
                    print(f"Uyarı: Kaydedilecek vis_img boş veya hatalı: {save_path}")
                elif getattr(self, '_img_writer', None) is not None:
                    # Asenkron yazım; kopya gerekli, imshow aynı tampona bakıyor
                    self._img_writer.submit(save_path, vis_img.copy())
                else:
                    success = cv2.imwrite(save_path, vis_img)
                    if not success:
                        print(f"HATA: Görüntü kaydedilemedi: {save_path}")

            if display_live:
                cv2.imshow("Live Tracking", vis_img)
        except Exception as e:
            print(f"Görselleştirme hatası: {e}")
            print(f"Satır: {e.__traceback__.tb_lineno}")

    def process_sequence(self, max_frames=None, visualize=False, viz_interval=1, export_landmarks_path=None, save_visualizations=True, on_frame=None):
        """Tüm sekansı işle. export_landmarks_path verilirse kare başına aktif landmark'lar
        pickle'a yazılır; on_frame verilirse her karede (ts, obs) çağrılır."""
        self.landmark_export = [] # (timestamp, aktif landmark listesi)
        try:
            timestamps = sorted(set(self.cam0_df['timestamp'].astype(str)))

            if max_frames:
                timestamps = timestamps[:max_frames]

            viz_window_name = "Live Tracking"
            if visualize:
                if save_visualizations:
                    self.viz_dir = f"visualizations_{self.extractor_name.lower()}_{self.matcher_type}_{self.filter_method}"
                    os.makedirs(self.viz_dir, exist_ok=True)
                    # PNG yazımı arka plan thread'inde; VIZ_ASYNC_WRITE=0 senkron yapar
                    if os.environ.get('VIZ_ASYNC_WRITE', '1') != '0':
                        self._img_writer = _AsyncImageWriter(
                            n_workers=int(os.environ.get('VIZ_WRITER_THREADS', '2')))
                cv2.namedWindow(viz_window_name, cv2.WINDOW_NORMAL)

            total_frames = len(timestamps)
            print(f"{total_frames} kare işleniyor...")

            start_time = datetime.now()

            for i, timestamp in enumerate(timestamps):
                # [TÜRBÜLANS] sekans içi mutlak kare indeksi: işlem başarısız olsa da ilerler
                self._seq_frame_idx = i
                print(f"\nKare işleniyor {i+1}/{total_frames} (Timestamp: {timestamp}, Dahili Frame Idx: {self.current_frame_idx})")

                # Stereo tetiği derinlikli track sayısına bakar (toplam değil)
                num_active_depth_tracks_at_prev_step = 0
                num_near_anchors = 0           # Z <= SCALE_NEAR_DEPTH olan anchor sayısı
                num_stale_anchors = 0          # yaşı SCALE_MAX_DEPTH_AGE'i aşan anchor sayısı
                if self.current_frame_idx >= 0:
                    cfi = self.current_frame_idx
                    for t in self.tracked_features:
                        if t['last_seen_frame_idx'] != cfi or not t.get('has_depth', False):
                            continue
                        # tetik yalniz YAKIN (Z < ThDepth*b) derinlikli track'leri sayar: ORB-SLAM3
                        # NeedNewKeyFrame (Tracking.cc:3116) anahtar kare kararinda yalniz yakin noktalari sayar.
                        if self.TRIGGER_CLOSE_DEPTH > 0:
                            _p = t.get('p3d', None)
                            if _p is None or not (0.0 < float(_p[2]) < self.TRIGGER_CLOSE_DEPTH):
                                continue
                        num_active_depth_tracks_at_prev_step += 1
                        p3d = t.get('p3d', None)
                        if p3d is not None and 0.0 < float(p3d[2]) <= self.SCALE_NEAR_DEPTH:
                            num_near_anchors += 1
                        if (cfi - t.get('depth_frame_idx', cfi)) > self.SCALE_MAX_DEPTH_AGE:
                            num_stale_anchors += 1

                print(f"  Önceki adımdaki aktif derinlikli track: {num_active_depth_tracks_at_prev_step}, Stereo Tetik Eşiği: {self.REPLENISH_THRESHOLD}")

                # Çekirdek tetik: ilk kare ya da derinlikli track sayısı eşiğin altında
                run_full_stereo = (self.current_frame_idx == -1) or \
                                  (num_active_depth_tracks_at_prev_step < self.REPLENISH_THRESHOLD)
                # AirSLAM anahtar kare kurali (map_builder.cc:451, tracking_point_rate=0.65): son stereo
                # karedeki derinlikli track'lerin orani bu degerin altina dusunce yeni stereo (anahtar) kare.
                if (not run_full_stereo and self.STEREO_TRACK_RATE > 0 and self._depth_at_last_stereo > 0 and
                        num_active_depth_tracks_at_prev_step < self.STEREO_TRACK_RATE * self._depth_at_last_stereo):
                    run_full_stereo = True

                # Ölçek-duyarlı ek tetik (varsayılan kapalı)
                if (self.SCALE_AWARE_TRIGGER and not run_full_stereo
                        and num_active_depth_tracks_at_prev_step > 0):
                    stale_frac = num_stale_anchors / num_active_depth_tracks_at_prev_step
                    near_scarce = (self.SCALE_NEAR_TRIGGER and
                                   num_near_anchors < self.SCALE_MIN_NEAR_ANCHORS)
                    stale_dominant = (self.SCALE_STALE_TRIGGER and
                                      stale_frac > self.SCALE_MAX_STALE_FRAC)
                    if near_scarce or stale_dominant:
                        run_full_stereo = True
                        self.n_scale_triggers += 1
                        reason = ("yakın-anchor kıtlığı" if near_scarce else "") + \
                                 (" + " if (near_scarce and stale_dominant) else "") + \
                                 ("eski-anchor baskın" if stale_dominant else "")
                        print(f"  #6 ÖLÇEK-DUYARLI TETİK ({reason}): "
                              f"yakın={num_near_anchors}/{self.SCALE_MIN_NEAR_ANCHORS}, "
                              f"eski-oran={stale_frac:.2f}/{self.SCALE_MAX_STALE_FRAC}")

                processed_data_tuple = None
                if run_full_stereo:
                    print(f"  Tam stereo işleme çalıştırılıyor...")
                    processed_data_tuple = self.process_frame_pair(timestamp)
                    if (self.BOS_KARE_DUZELT and processed_data_tuple is not None and processed_data_tuple[0] is None):
                        # stereo eslesmesi bu karede basarisiz (or. y-filtresi sonrasi eslesme yok): kare
                        # dusurulmez, yalniz zamansal takiple islenir.
                        self.n_stereo_geri_dusme += 1
                        processed_data_tuple = self.process_frame_temporal_only(timestamp)
                else:
                    print(f"  Sadece temporal takip çalıştırılıyor...")
                    processed_data_tuple = self.process_frame_temporal_only(timestamp)

                # HATA DUZELTMESI: islem fonksiyonlari basarisizlikta None degil (None, None, None, None) dondurur;
                # duzeltme oncesi bu dal hic calismaz, kare indeksi ilerlemez ve onceki karenin gozlemleri bu karenin damgasiyla
                # disa aktarilirdi (bayat gozlem). Simdi: basarisiz kare track'leri yaslandirir ve BOS gozlemle aktarilir
                # (arka uc o kareyi sabit-hiz ongorusuyle pozlar).
                _bos = (self.BOS_KARE_DUZELT and processed_data_tuple is not None and processed_data_tuple[0] is None)
                if _bos:
                    self.n_bos_kare += 1
                    _ = self.update_and_select_tracks(np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
                    if export_landmarks_path is not None:
                        self.landmark_export.append((timestamp, []))
                        if getattr(self, 'EXPORT_KEYPOINTS', False):
                            self.keypoint_export.append((timestamp, np.zeros((0, 2), dtype=np.float32), None))
                    if on_frame is not None:
                        on_frame(str(timestamp), landmarks_to_obs([]), self.current_frame_idx)
                    continue
                if processed_data_tuple is None:
                    print(f"  Uyarı: {timestamp} için işlem başarısız oldu (stereo veya temporal).")
                    print(f"  self.current_frame_idx: {self.current_frame_idx} (işlem öncesi)")
                    # İndeksi ilerletip track'leri yaşlandırmak için boş çağrı
                    _ = self.update_and_select_tracks(
                        np.array([]), np.array([]), np.array([]),
                        np.array([]), np.array([])
                    )
                    print(f"  self.current_frame_idx: {self.current_frame_idx} (update_and_select_tracks sonrası)")

                    # prev_vis_info'yu güvenli duruma getir (içeride güncellenmediyse)
                    current_filename0_for_failed_frame = "UNKNOWN"
                    filename0_series_fail = self.cam0_df.loc[self.cam0_df['timestamp'] == timestamp, 'filename']
                    if not filename0_series_fail.empty:
                        current_filename0_for_failed_frame = filename0_series_fail.iloc[0]

                    if self.prev_vis_info['filename0'] != current_filename0_for_failed_frame :
                        self.prev_vis_info['image'] = None
                        self.prev_vis_info['kp_candidates'] = None
                        self.prev_vis_info['filename0'] = current_filename0_for_failed_frame
                    continue

                result, vis_data, filenames_stereo, filenames_temporal = processed_data_tuple
                self.results.append(result)
                if result is not None and result.get('ran_stereo'):
                    self._depth_at_last_stereo = int(result.get('active_depth_tracks', 0))

                # Backend için bu karenin aktif landmark gözlemleri
                _active = None
                if export_landmarks_path is not None or on_frame is not None:
                    _active = self.get_active_landmarks()
                if export_landmarks_path is not None:
                    self.landmark_export.append((timestamp, _active))
                    if getattr(self, 'EXPORT_KEYPOINTS', False):
                        self.keypoint_export.append((
                            timestamp,
                            getattr(self, '_cur_all_kps', np.zeros((0, 2), dtype=np.float32)),
                            getattr(self, '_cur_all_descs', None)))
                # Gerçek zamanlı kanca: bu karenin (N,8) gözlemlerini backend'e ver
                if on_frame is not None:
                    on_frame(str(timestamp), landmarks_to_obs(_active), self.current_frame_idx)

                if visualize and i % viz_interval == 0 and vis_data is not None:
                    save_path = f"{self.viz_dir}/combined_tracking_{timestamp}.png" if save_visualizations else None
                    self.visualize_combined_tracking(
                        vis_data['stereo'],
                        vis_data['temporal'],
                        filenames_stereo,
                        filenames_temporal,
                        result,
                        save_path=save_path,
                        display_live=True
                    )

                if visualize:
                    # 'q' ile çıkış
                    if cv2.waitKey(1) & 0xFF == ord('q'):
                        print("Kullanıcı tarafından çıkış yapılıyor...")
                        break 

            end_time = datetime.now()
            self.execution_time = (end_time - start_time).total_seconds()
            # Süre ölçümü PNG kuyruğu boşalmadan önce alınır
            self._flush_img_writer()

            print("\nİşlem tamamlandı!")
            self.save_results()
            if export_landmarks_path is not None:
                self._save_landmark_export(export_landmarks_path)

        except KeyboardInterrupt:
            print("\nKullanıcı tarafından işlem kesildi (Ctrl+C).")
            self._flush_img_writer()
            self.save_results()
            if export_landmarks_path is not None:
                self._save_landmark_export(export_landmarks_path)
        except Exception as e:
            self._flush_img_writer()
            print(f"Sekans işleme hatası: {e}")
            raise

    def _flush_img_writer(self):
        """PNG kuyruğunu boşalt ve yazıcı thread'lerini kapat; her çıkış yolundan çağrılır."""
        w = getattr(self, '_img_writer', None)
        if w is not None:
            n = w.q.qsize()
            if n:
                print(f"  Bekleyen {n} görsel diske yazılıyor...")
            w.close()
            self._img_writer = None

    def get_calibration_export(self):
        """Backend için kalibrasyon: rektifiye intrinsics (P1), T_BS, R1, baseline."""
        return {
            'rect_K': [float(self.P1[0, 0]), float(self.P1[1, 1]),
                       float(self.P1[0, 2]), float(self.P1[1, 2])],
            'T_BS': np.array(self.cam0_params['T_BS']['data'], dtype=np.float64).reshape(4, 4),
            'R1': np.asarray(self.R1, dtype=np.float64),
            'baseline': float(np.linalg.norm(self.T)),
        }

    def _save_landmark_export(self, path):
        """Kare başına landmark gözlemlerini pickle'a yaz:
        {rect_K, T_BS, R1, baseline, frames: [{timestamp, obs (N,8)}], has_keypoints}."""
        import pickle
        frames = [{'timestamp': str(ts), 'obs': landmarks_to_obs(lms)}
                  for ts, lms in self.landmark_export]
        if self.EXPORT_S_D:
            for fr, (ts, _) in zip(frames, self.landmark_export):
                if ts in self._sd_by_ts:
                    fr['s_d'] = self._sd_by_ts[ts]
        # export_keypoints açıkken kare başına tüm kp + descriptor eklenir
        if getattr(self, 'EXPORT_KEYPOINTS', False) and self.keypoint_export:
            kp_by_ts = {str(ts): {'kps': np.asarray(kps, dtype=np.float32),
                                  'desc': (None if d is None else np.asarray(d, dtype=np.float32))}
                        for ts, kps, d in self.keypoint_export}
            for fr in frames:
                kd = kp_by_ts.get(fr['timestamp'])
                if kd is not None:
                    fr['kps'] = kd['kps']
                    fr['desc'] = kd['desc']
        export = dict(self.get_calibration_export(), frames=frames)
        export['has_keypoints'] = bool(getattr(self, 'EXPORT_KEYPOINTS', False))
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        with open(path, 'wb') as f:
            pickle.dump(export, f)
        n_obs = sum(fr['obs'].shape[0] for fr in frames)
        extra = ""
        if export['has_keypoints']:
            n_kp = sum(fr.get('kps', np.zeros((0, 2))).shape[0] for fr in frames)
            extra = f", {n_kp} keypoint (export_keypoints AÇIK)"
        print(f"Landmark export kaydedildi: {path} ({len(frames)} kare, {n_obs} gözlem{extra})")

    def save_results(self):
        """Kare başına sonuçları CSV'ye yaz."""
        if not self.results:
            print("Kaydedilecek sonuç yok!")
            return

        # SLAM_FRAME_STATS_DIR verilmişse CSV oraya
        results_dir = (os.environ.get("SLAM_FRAME_STATS_DIR")
                       or f"results_{self.extractor_name.lower()}_{self.matcher_type}_{self.filter_method}")
        os.makedirs(results_dir, exist_ok=True)

        # dict olmayan bozuk kayıtları süz
        self.results = [r for r in self.results if isinstance(r, dict)]
        df = pd.DataFrame(self.results)
        output_file = f"{results_dir}/stereo_matching_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        df.to_csv(output_file, index=False)
        print(f"Sonuçlar kaydedildi: {output_file}")

    def compute_statistics(self):
        """Sekans istatistiklerini hesapla, yazdır ve sözlük olarak döndür."""
        if not self.results:
            print("İstatistik hesaplanacak sonuç yok!")
            return None

        self.results = [r for r in self.results if isinstance(r, dict)]
        df = pd.DataFrame(self.results)

        # Eksik kolona dayanıklı erişim (NaN)
        def _m(col):
            return df[col].mean() if col in df.columns else float('nan')

        def _s(col):
            return df[col].std() if col in df.columns else float('nan')

        stats = {
            'avg_kp1_count': _m('kp1_count'),
            'std_kp1_count': _s('kp1_count'),
            'avg_kp2_count': _m('kp2_count'),
            'std_kp2_count': _s('kp2_count'),
            'avg_stereo_matches': _m('stereo_matches'),
            'std_stereo_matches': _s('stereo_matches'),
            'avg_current_stereo_triangulated_points': _m('current_stereo_triangulated_points'),
            'std_current_stereo_triangulated_points': _s('current_stereo_triangulated_points'),
            'avg_temporal_matches': _m('temporal_matches'),
            'std_temporal_matches': _s('temporal_matches'),
            'total_frames_processed': len(df),
            'total_execution_time_seconds': self.execution_time if self.execution_time is not None else 0
        }

        # Backend'e giden landmark sayısı ve stereo'nun ne kadar seyrek koştuğu
        if 'active_depth_tracks' in df.columns:
            stats['avg_active_depth_tracks'] = df['active_depth_tracks'].mean()
            stats['std_active_depth_tracks'] = df['active_depth_tracks'].std()
        if 'ran_stereo' in df.columns:
            stereo_run_ratio = float(df['ran_stereo'].mean())
            stats['stereo_run_ratio'] = stereo_run_ratio
            stats['stereo_skip_ratio'] = 1.0 - stereo_run_ratio
        stats['scale_aware_triggers'] = int(self.n_scale_triggers)
        stats['bos_kare'] = int(getattr(self, 'n_bos_kare', 0))
        stats['stereo_geri_dusme'] = int(getattr(self, 'n_stereo_geri_dusme', 0))
        stats['lg_zamansal'] = int(getattr(self, 'n_lg_temporal', 0))
        if len(df) > 0:
            stats['scale_aware_trigger_ratio'] = float(self.n_scale_triggers) / len(df)

        # Adaptif ön işleme metrikleri; kapalıyken sütunlar NaN, blok atlanır
        stats['adaptive_preprocessing'] = 1.0 if self.ADAPTIVE_PREPROCESSING_ENABLED else 0.0
        if 's_d1' in df.columns and df['s_d1'].notna().any():
            stats['avg_s_d1'] = float(df['s_d1'].mean())
            stats['std_s_d1'] = float(df['s_d1'].std())
            stats['max_s_d1'] = float(df['s_d1'].max())
        if 's_d1_roi' in df.columns and df['s_d1_roi'].notna().any():
            # ROI içi S_d yalnız teşhis; karar tüm görüntüyle verilir
            stats['avg_s_d1_roi'] = float(df['s_d1_roi'].mean())
            stats['avg_s_d1_roi_bias'] = float((df['s_d1'] - df['s_d1_roi']).mean())
        if 'clahe_clip1' in df.columns and df['clahe_clip1'].notna().any():
            applied = df['clahe_clip1'] > 0.0
            stats['clahe_applied_ratio'] = float(applied.mean())
            stats['avg_clahe_clip1'] = float(df.loc[applied, 'clahe_clip1'].mean()) if applied.any() else 0.0
        if 'det_threshold' in df.columns and df['det_threshold'].notna().any():
            stats['avg_det_threshold'] = float(df['det_threshold'].mean())
            stats['min_det_threshold'] = float(df['det_threshold'].min())

        # Ayıklama sayaçları (kapalıyken sözlük boş)
        for k, v in sorted(self.cull_stats.items()):
            stats[f'cull_{k}'] = float(v)

        print("\nEşleştirme İstatistikleri:")
        print("-" * 50)
        for key, value in stats.items():
            print(f"{key}: {value:.2f}")
        print("-" * 50)
        if self.execution_time is not None:
            print(f"Toplam Çalışma Süresi: {self.execution_time:.2f} saniye")

        return stats

def main():
    try:
        parser = argparse.ArgumentParser(description="Özellik tabanlı Stereo Eşleştirme Scripti.")
        # default=None: verilmeyen argüman kod içi varsayılana bırakılır
        parser.add_argument('--matcher_type', type=str, default=None,
                            choices=['lightglue', 'classical_bf'],
                            help="Kullanılacak stereo eşleştirici türü: 'lightglue' veya 'classical_bf'. "
                                 "Eğer belirtilmezse, kod içindeki varsayılan kullanılır.")
        parser.add_argument('--sequence', type=str, default='MH_05_difficult', help="İşlenecek EurocMAV sekansı.")
        parser.add_argument('--base_path', type=str, default="", help="EurocMAV veri setinin temel yolu.")
        parser.add_argument('--filter_method', type=str, default=None,
                            choices=['fundamental_matrix', 'y_axis_check', 'cross_check', 'none'],
                            help="Stereo eşleştirmeler için kullanılacak filtreleme yöntemi. "
                                 "Eğer belirtilmezse, kod içindeki varsayılan kullanılır. "
                                 "Seçenekler: 'fundamental_matrix', 'y_axis_check', 'cross_check', 'none'.")
        parser.add_argument('--y_axis_threshold', type=float, default=None,
                            help="Y ekseni kontrol filtresi için piksel cinsinden eşik değeri. "
                                 "Eğer belirtilmezse, kod içindeki varsayılan kullanılır. "
                                 "Örn: 8.0.")
        parser.add_argument('--extractor', type=str, default=None,
                            choices=['superpoint', 'aliked'],
                            help="Özellik çıkarıcı. Belirtilmezse kod içi varsayılan (superpoint) kullanılır.")
        parser.add_argument('--max_track_length', type=int, default=None,
                            help="Track yaş sınırı (kare). Belirtilmezse kod içi varsayılan (60) kullanılır. 0=kapalı.")

        args = parser.parse_args()

        device_info = "GPU" if torch.cuda.is_available() else "CPU"
        print(f"Bu kod şu anda {device_info} üzerinde çalışıyor.")

        # Yalnız kullanıcının verdiği argümanlar kurucuya geçer
        init_args_from_cli = {k: v for k, v in vars(args).items() if v is not None}

        base_path_arg = init_args_from_cli.pop('base_path', args.base_path)
        sequence_arg = init_args_from_cli.pop('sequence', args.sequence)

        matcher = StereoMatcher(base_path=base_path_arg, sequence=sequence_arg, **init_args_from_cli)

        matcher.process_sequence(
            max_frames=None,
            visualize=True,
            viz_interval=1
        )

        stats = matcher.compute_statistics()

    except Exception as e:
        print(f"Ana program hatası: {e}")
        raise

if __name__ == "__main__":
    main()
