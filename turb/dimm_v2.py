"""DIMM kestiricisinin AMPİRİK ÇEKİRDEKLİ sürümü — hazırlık, henüz devrede DEĞİL.

KUSUR (2026-08-22'de E4'ün 672 kaydından ölçülerek bulundu).
`turb/dimm.py::estimate` yapı fonksiyonuna FİZİKSEL Kolmogorov tensörünü
(`correlation.rho_long/rho_tran`) uyduruyor. Oysa kovaryans yolu 2026-08-19'da
E2 ile kalibre edilmiş AMPİRİK çekirdeği kullanıyor:

    rho(dtheta) = A * exp(-(dtheta/theta_c)^p),   A = 0.75,  p = 0.40

Fiziksel çekirdek küçük ayrımlarda çok daha hızlı düşer. Ampirik veriye fiziksel
model uydurunca theta_c aşağı itiliyor ve ALT SINIRA ÇAKILIYOR. Ölçülen:

    theta_c_hat = 33.333 px  (D/r0 = 3, 5, 8'in ÜÇÜNDE DE, hücreden hücreye sabit)
    tc_lo       = max(2 * ilk_bin_merkezi, 4) = max(2*16.667, 4) = 33.333 px
    oracle      = 112.8 px

Yani kestirim hiçbir şey ÖLÇMÜYOR, sınırı raporluyor. Sonuçları:
  * A3'ün yamaları 3.3 kat küçük -> yapı bilgisi kaybı
  * sigma 1.25-2.01 kat fazla kestiriliyor -> A3 tutarlılığı YAPIYLA değil
    ŞİŞİRMEYLE satın alıyor (A0C'nin yaptığının aynısı) ve NEES'i yanıltıcı
    biçimde iyi görünüyor (38.9; oracle kolu A5'inki 100.3)

A3, "oracle bilgisi olmadan işe yarıyor mu?" sorusunun cevabıdır; hakemin ilk
soracağı şeydir. Sakat bir kestiriciyle raporlamak yöntemi olduğundan KÖTÜ
gösterir.

DÜZELTME. Yapı fonksiyonu modeli, kovaryans yolunun kullandığı çekirdekle AYNI
olmalıdır:

    D(sep) = 2 * s^2 * (1 - A * exp(-(sep/theta_c)^p))

A ve p, E2'de kalibre edildikleri değerlerde SABİT tutulur (serbest bırakmak
12 kutuda üç parametreyi dejenere yapar). Ayrıca alt sınır gevşetilir: eski
sınırın gerekçesi ölçülebilirlikti ve doğruydu, ama 2x ilk-bin-merkezi fazla
agresif — asıl sorun sınır değil MODEL UYUŞMAZLIĞIYDI.

DEVREYE ALMA: koşan kuyruk bitince `turb/dimm.py` bununla değiştirilecek,
`dimm.py` kod parmak izi listesine EKLENECEK (şu an listede yok — A3/A4
sonuçlarını etkilediği hâlde; bu ayrı bir eksikti) ve kola bağlı aşamalar
(E4, E4model, E9, TARTANarms) yeniden koşulacak.
"""
from __future__ import annotations

import numpy as np

from .covariance import EMP_A, EMP_P
from .dimm import DifferentialMotionMonitor as _Base
from .optics import C_GTILT
from .dimm import r0_from_sigma_px


class DifferentialMotionMonitorV2(_Base):
    """Ampirik çekirdekli DIMM. add_frame/structure_function tabandan gelir.

    `last_pinned`: son kestirim SINIRA yapıştı mı. True ise dönen theta_c bir
    ÖLÇÜM DEĞİLDİR ve kullanan taraf bunu bilmelidir (v1'de böyle bir bayrak
    yoktu; sınır sessizce ölçüm diye raporlanıyordu).
    """

    last_pinned = False
    last_bounds = (float("nan"), float("nan"))

    def estimate(self, theta_c_px: float | None = None):
        """(sigma_px, theta_c_px, r0_m). Yetersiz veri -> (nan, nan, nan)."""
        sep, Dv, C = self.structure_function()
        ok = np.isfinite(Dv) & (C > 10)
        if ok.sum() < 3 or C.sum() < self.min_pairs:
            return np.nan, np.nan, np.nan

        # KOVARYANS YOLUYLA AYNI ÇEKİRDEK (turb/covariance.py EMP_A, EMP_P).
        # rho(0) = A < 1: nugget (dedektör gürültüsü + korelasyonsuz kalıntı)
        # genliğin İÇİNDEDİR, bu yüzden AYRI bir gürültü terimi EKLENMEZ —
        # eklemek nugget'i iki kez sayardı (2026-08-20'de kovaryans yolunda
        # bulunan ve düzeltilen hatanın aynısı).
        # BILESEN SAYIMI (2026-08-23'te bulunan IKINCI kusur, ilkinden bagimsiz).
        # add_frame `dd = sum((r_a - r_b)**2, axis=1)` ile IKI EKSENI toplar:
        #     D(sep) = E|dr|^2 = 2 * [2 sigma^2 (1-rho)] = 4 sigma^2 (1-rho)
        # Eski model katsayisi 2.0 idi -> s2 = 2 sigma^2 uyduruyor ->
        # sigma_hat = sqrt(2) * sigma_gercek, yani SISTEMATIK %41 fazla kestirim.
        # E4'te olculen sigma_hat/sigma_oracle = 1.25-2.01 araliginin bir kismi
        # tam olarak budur (sqrt(2) = 1.414); kalani cekirdek uyusmazligindandi.
        # Bu, r0 kestirimini de (sigma^-6/5 ile) yanli kilar -> E6 eksenini etkiler.
        def model(sep_px, s2, tc_px):
            x = np.maximum(sep_px, 0.0) / max(tc_px, 1e-6)
            return 4.0 * s2 * (1.0 - EMP_A * np.exp(-(x ** EMP_P)))

        from scipy.optimize import curve_fit
        x, y, w = sep[ok], Dv[ok], np.sqrt(C[ok])
        try:
            if theta_c_px is None:
                # ÖLÇÜLEBİLİRLİK SINIRI korunur ama gevşetilir: alt sınır artık
                # 2x değil 1x ilk ayrım (0.5x güvenlik payıyla), üst sınır
                # 4x en büyük ayrım. Gerçek theta_c (~113 px) bu aralığın
                # ORTASINDA kalmalı; sınıra dayanırsa P14 kırmızı verir.
                # OLCULEBILIRLIK: en uzun ayrimdan UZUN bir korelasyon boyu
                # olculemez -- yapi fonksiyonu o araligin otesinde duz olur ve
                # uydurma gurultuye oturur (olculdu: tc_hi = 4*max_sep iken
                # gercek 200 px icin 1533 px donuyordu, %667 hata).
                # Ust sinir 1.5 * en uzun ayrim.
                tc_lo = max(0.5 * float(np.min(x)), 2.0)
                tc_hi = max(1.5 * float(np.max(x)), tc_lo * 4.0)
                # COK BASLANGICLI UYDURMA. Tek p0 ile yuksek sigma'da yerel
                # minimuma dusuyordu (olculdu: theta_c=112.8 gercekte iken 63.9,
                # theta_c=200 iken 282 veriyordu). Birkac baslangic denenip
                # agirlikli kalinti karesi EN KUCUK olan secilir.
                best = None
                _s0 = max(np.nanmax(y) / 4.0, 0.1)
                for _f in (0.25, 0.5, 1.0, 2.0, 4.0):
                    _tc0 = float(np.clip(np.nanmedian(x) * _f, tc_lo, tc_hi))
                    try:
                        _po, _ = curve_fit(model, x, y, p0=[_s0, _tc0],
                                           sigma=1.0 / w,
                                           bounds=([1e-6, tc_lo], [1e4, tc_hi]),
                                           maxfev=20000)
                    except Exception:
                        continue
                    _r = float(np.sum(w * (model(x, *_po) - y) ** 2))
                    if best is None or _r < best[0]:
                        best = (_r, _po)
                if best is None:
                    return np.nan, np.nan, np.nan
                popt = best[1]
                s2, tc = float(popt[0]), float(popt[1])
                self.last_bounds = (tc_lo, tc_hi)
                # SINIRA YAPISMA TESPITI — v1'in ASIL KUSURU BUYDU.
                # v1 her hucrede alt siniri (33.3 px) donduruyordu ve bunu
                # OLCUM gibi raporluyordu; hicbir sey uyarmiyordu. Sinira
                # yapisan bir kestirim bir OLCUM DEGIL, "olculemedi"dir.
                self.last_pinned = bool(tc <= tc_lo * 1.01 or tc >= tc_hi * 0.99)
            else:
                f = lambda xx, s2_: model(xx, s2_, theta_c_px)     # noqa: E731
                popt, _ = curve_fit(f, x, y, p0=[max(np.nanmax(y) / 2.0, 0.1)],
                                    sigma=1.0 / w, bounds=([1e-6], [1e4]),
                                    maxfev=20000)
                s2, tc = float(popt[0]), float(theta_c_px)
                self.last_bounds = (np.nan, np.nan)
                self.last_pinned = False
        except Exception:
            return np.nan, np.nan, np.nan

        sigma_px = float(np.sqrt(max(s2, 0.0)))
        r0 = r0_from_sigma_px(sigma_px, self.D, self.lam, self.px_per_rad)
        return sigma_px, tc, r0
