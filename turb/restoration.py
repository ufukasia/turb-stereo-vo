"""
Türbülans giderme (restorasyon) kolları — C2 / E3.

TEZ (rapor §3-C2): restorasyon yöntemleri ALGISAL kaliteye (PSNR/SSIM) göre
optimize edilir; kareler arası GEOMETRİK tutarlılığı garanti etmezler. Dolayısıyla
"önce restore et, sonra VIO çalıştır" hattında PSNR iyileşmesi ATE iyileşmesine
DÖNÜŞMEYEBİLİR — hatta ters yönde ilişkilenebilir. Bu bir hipotezdir; ölçüm sonucu
ne olursa olsun raporlanır.

Kollar:
    (a) none            bozuk -> VIO                                  (referans)
    (b) lucky           şanslı-kare seçimi (keskinlik ölçütü + zamansal medyan)
    (c) temporal_median kayan pencere zamansal medyan (klasik dalgalanma bastırma)
    (d) elastic         zamansal medyan referansına elastik kayıt (optik akış)
                        + keskinleştirme — klasik "türbülans giderme" hattının özü
    (e) clean           temiz -> VIO                                  (üst sınır)

(d) kolu, klasik literatürdeki "registration + fusion" ailesinin (Zhu & Milanfar,
Mao et al.) hafif bir uygulamasıdır: kareler ortak bir referansa hizalanır, böylece
DALGALANMA bastırılır. Kritik nokta: bu hizalama her kareyi KENDİ referansına
çektiği için sahnenin gerçek hareketiyle karışır — C2'nin ölçtüğü tam da budur.

Derin modeller (TMT / MambaTM) isteğe bağlıdır: ağırlıkları yerelde mevcutsa
`--deep <ad>` ile eklenir; yoksa klasik kollar tek başına C2 için yeterlidir
(hipotez, "restorasyon geometrik tutarlılığı garanti etmez" biçiminde geneldir).
"""
from __future__ import annotations

from collections import deque

import cv2
import numpy as np


class Restorer:
    """Kare-akışı üzerinde çalışan restorasyon kolu (nedensel, kayan pencere)."""

    def __init__(self, method: str = "none", window: int = 7,
                 sharpen: float = 0.6, flow_levels: int = 3):
        self.method = method
        self.window = int(window) | 1
        self.sharpen = float(sharpen)
        self.flow_levels = int(flow_levels)
        self._buf = deque(maxlen=self.window)

    def reset(self):
        self._buf.clear()

    # ------------------------------------------------------------------ ana
    def __call__(self, img: np.ndarray) -> np.ndarray:
        if self.method == "none":
            return img
        self._buf.append(img.astype(np.float32))
        if len(self._buf) < 3:
            return img
        if self.method == "temporal_median":
            return self._median()
        if self.method == "lucky":
            return self._lucky()
        if self.method == "elastic":
            return self._elastic()
        raise ValueError(f"bilinmeyen restorasyon kolu: {self.method}")

    # -------------------------------------------------------------- yöntemler
    def _median(self) -> np.ndarray:
        m = np.median(np.stack(self._buf, axis=0), axis=0)
        return np.clip(m, 0, 255).astype(np.uint8)

    @staticmethod
    def _sharpness(a: np.ndarray) -> float:
        return float(cv2.Laplacian(a, cv2.CV_32F).var())

    def _lucky(self) -> np.ndarray:
        """Şanslı-kare: pencere içindeki EN KESKİN kare seçilir.

        Klasik "lucky imaging" ilkesi. Geometrik açıdan kritik kusuru: seçilen kare
        MEVCUT kare olmayabilir -> zaman ekseninde kayma (VIO için zehirli).
        """
        idx = int(np.argmax([self._sharpness(a) for a in self._buf]))
        return np.clip(self._buf[idx], 0, 255).astype(np.uint8)

    def _elastic(self) -> np.ndarray:
        """Zamansal medyan referansına elastik kayıt + keskinleştirme."""
        cur = self._buf[-1]
        ref = np.median(np.stack(self._buf, axis=0), axis=0).astype(np.float32)
        flow = cv2.calcOpticalFlowFarneback(
            cv2.convertScaleAbs(cur), cv2.convertScaleAbs(ref), None,
            0.5, self.flow_levels, 21, 3, 5, 1.2, 0)
        h, w = cur.shape
        gx, gy = np.meshgrid(np.arange(w, dtype=np.float32),
                             np.arange(h, dtype=np.float32))
        warped = cv2.remap(cur, gx + flow[..., 0], gy + flow[..., 1],
                           cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT_101)
        if self.sharpen > 0:
            blur = cv2.GaussianBlur(warped, (0, 0), 1.5)
            warped = cv2.addWeighted(warped, 1 + self.sharpen, blur, -self.sharpen, 0)
        return np.clip(warped, 0, 255).astype(np.uint8)


# =============================================================================
class RestorationHook:
    """StereoMatcher'a takılan kanca: türbülanstan SONRA restorasyon uygular.

    Her kamera için ayrı tampon tutulur (stereo tutarlılığı korunur).
    PSNR/SSIM ölçümü için temiz kare de verilebilir.
    """

    def __init__(self, method: str, window: int = 7, sharpen: float = 0.6):
        self.method = method
        self.r0 = Restorer(method, window, sharpen)
        self.r1 = Restorer(method, window, sharpen)
        self.psnr_sum = 0.0
        self.ssim_sum = 0.0
        self.n = 0

    def apply(self, img, cam: int):
        return (self.r0 if cam == 0 else self.r1)(img)

    def score(self, clean, restored):
        from .metrics import psnr, ssim
        self.psnr_sum += psnr(clean, restored)
        self.ssim_sum += ssim(clean, restored)
        self.n += 1

    def summary(self):
        return {"method": self.method, "n": self.n,
                "psnr": self.psnr_sum / self.n if self.n else float("nan"),
                "ssim": self.ssim_sum / self.n if self.n else float("nan")}
