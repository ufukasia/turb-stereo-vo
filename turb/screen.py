"""
Kolmogorov faz ekranı üretimi (FFT + subharmonik alçak-frekans telafisi).

Neden subharmonik ŞART: Kolmogorov PSD'si Phi_phi(f) ~ f^{-11/3} olduğundan faz
varyansı en alçak uzaysal frekanslarda yoğunlaşır. Düz FFT ekranı, ekran boyu
L_scr'nin altındaki frekansları hiç içermez ve yapı fonksiyonunu sistematik olarak
eksik üretir (ölçüldü: r = 1 m'de %58 eksik). Lane/Johansson-Gavel subharmonik
şeması DC hücresini özyinelemeli 3x3 alt hücrelere bölerek bu enerjiyi geri kazandırır.
Bu proje TILT istatistiğini ölçtüğü için düzeltme isteğe bağlı değildir.

TASARIM: `KolmogorovRealization` gerçeklemeyi SPEKTRAL katsayı düzeyinde tutar.
Bu sayede fazın herhangi bir DOĞRUSAL SÜZGEÇLENMİŞ hâli (ör. açıklık-ortalamalı
gradyan = G-tilt) doğrudan frekans uzayında, artefaktsız üretilebilir:

    * Sonlu farkla türev alıp sonra ayrık diskle ortalamak YANLIŞTIR: saf Kolmogorov
      fazının gradyan varyansı iç ölçek olmadan ıraksar; ayrık türev şablonu Nyquist
      civarını yükseltir ve tilt varyansı ızgara adımına bağımlı hale gelir.
    * Doğru yol, tam analitik transfer fonksiyonunu uygulamaktır:
          H_x(f) = 2*pi*i*f_x * A(f) / k ,   A(f) = 2 J1(pi D f)/(pi D f)
      Bu, açıklık üzerinde ortalanmış dalga cephesi eğimini VERİR ve integrali
      yakınsaktır. Sayısal doğrulama (quad): katsayı 0.1706 — literatürdeki G-tilt
      değeri 0.170 ile %0.3 uyum (rapor Ek B/5'teki Z-tilt/G-tilt belirsizliği
      böylece kapanır: 0.170 G-tilt'e aittir).
    * Subharmonik bileşenler ekran üzerinde periyodik OLMADIĞINDAN spektral türev
      onlara uygulanamaz; bunlar düzlem dalga oldukları için transfer fonksiyonu
      ANALİTİK olarak uygulanır (aşağıda `filtered`).

Referans: J. D. Schmidt, "Numerical Simulation of Optical Wave Propagation", §9.
"""
from __future__ import annotations

import numpy as np


def kolmogorov_psd(f: np.ndarray, r0: float, L0: float | None = None,
                   l0: float | None = None) -> np.ndarray:
    """Faz PSD'si [rad^2 m^2]; f çevrim/m.

    Kolmogorov:       Phi = 0.023 r0^{-5/3} f^{-11/3}
    von Karman (L0):  Phi = 0.023 r0^{-5/3} exp(-(f/fm)^2) / (f^2 + 1/L0^2)^{11/6}
    """
    f = np.asarray(f, dtype=np.float64)
    if L0 is None and l0 is None:
        with np.errstate(divide="ignore", invalid="ignore"):
            psd = 0.023 * r0 ** (-5.0 / 3.0) * f ** (-11.0 / 3.0)
        psd[~np.isfinite(psd)] = 0.0
        return psd
    f0 = 0.0 if L0 is None else 1.0 / L0
    fm = np.inf if l0 is None else 5.92 / l0 / (2 * np.pi)
    with np.errstate(divide="ignore", invalid="ignore"):
        psd = (0.023 * r0 ** (-5.0 / 3.0) * np.exp(-((f / fm) ** 2))
               / (f ** 2 + f0 ** 2) ** (11.0 / 6.0))
    psd[~np.isfinite(psd)] = 0.0
    return psd


class KolmogorovRealization:
    """Tek bir atmosfer gerçeklemesinin spektral gösterimi.

    Alanlar:
        cn      (N,N) kompleks ızgara katsayıları (kaydırılmış frekans ızgarasında)
        FX, FY  (N,N) frekans ızgarası [çevrim/m]
        subs    [(fx, fy, c)] subharmonik düzlem dalga modları
    """

    def __init__(self, N: int, delta: float, r0: float, rng: np.random.Generator,
                 n_subharmonics: int = 8, L0: float | None = None,
                 l0: float | None = None):
        self.N, self.delta, self.r0 = int(N), float(delta), float(r0)
        delta_f = 1.0 / (N * delta)
        f1 = (np.arange(-N // 2, N // 2)) * delta_f
        self.FX, self.FY = np.meshgrid(f1, f1)
        F = np.sqrt(self.FX ** 2 + self.FY ** 2)

        psd = kolmogorov_psd(F, r0, L0, l0)
        psd[N // 2, N // 2] = 0.0
        # (randn + i randn) formunda kompleks katsayının gücü 2*PSD*df^2'dir; reel
        # kısım alındığında yarısı kalır -> reel alanın spektral yoğunluğu tam PSD.
        self.cn = ((rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N)))
                   * np.sqrt(psd) * delta_f)

        self.subs = []
        D_scr = N * delta
        for p in range(1, int(n_subharmonics) + 1):
            df = 1.0 / (3.0 ** p * D_scr)
            fx = np.array([-1.0, 0.0, 1.0]) * df
            SFX, SFY = np.meshgrid(fx, fx)
            SF = np.sqrt(SFX ** 2 + SFY ** 2)
            spsd = kolmogorov_psd(SF, r0, L0, l0)
            spsd[1, 1] = 0.0
            c = ((rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))
                 * np.sqrt(spsd) * df)
            for i in range(3):
                for j in range(3):
                    if c[i, j] != 0:
                        self.subs.append((SFX[i, j], SFY[i, j], c[i, j]))

        x = (np.arange(-N // 2, N // 2)) * delta
        self._x = x

    # ------------------------------------------------------------------ alanlar
    def _spectral_sum(self, C: np.ndarray) -> np.ndarray:
        N = self.N
        return np.fft.ifftshift(np.fft.ifft2(np.fft.ifftshift(C))) * (N ** 2)

    def filtered(self, H_grid: np.ndarray, H_sub) -> np.ndarray:
        """Doğrusal süzgeç uygulanmış reel alan.

        H_grid: (N,N) kompleks transfer fonksiyonu (ızgara modları için)
        H_sub:  çağrılabilir H(fx, fy) -> kompleks (subharmonik modlar için)
        """
        out = self._spectral_sum(self.cn * H_grid).real
        if self.subs:
            lo = np.zeros((self.N, self.N), dtype=np.complex128)
            for fx, fy, c in self.subs:
                ex = np.exp(2j * np.pi * fx * self._x)
                ey = np.exp(2j * np.pi * fy * self._x)
                lo += (c * H_sub(fx, fy)) * np.outer(ey, ex)
            out = out + lo.real
        return out

    def phase(self) -> np.ndarray:
        """Faz ekranı [rad] (ortalaması çıkarılmış)."""
        one = np.ones((self.N, self.N), dtype=np.complex128)
        ph = self.filtered(one, lambda fx, fy: 1.0 + 0j)
        return ph - ph.mean()


def phase_screen(N: int, delta: float, r0: float, rng: np.random.Generator,
                 n_subharmonics: int = 8, L0: float | None = None,
                 l0: float | None = None) -> np.ndarray:
    """Kısayol: N x N Kolmogorov faz ekranı [rad]."""
    return KolmogorovRealization(N, delta, r0, rng, n_subharmonics, L0, l0).phase()


def structure_function_shift(phz: np.ndarray, delta: float, shifts) -> tuple:
    """Ampirik yapı fonksiyonu D_phi(r) = <[phi(x+r)-phi(x)]^2> (eksen boyu kaydırma).

    Kolmogorov teorisi: D_phi(r) = 6.88 (r/r0)^{5/3}.
    """
    r, D = [], []
    for s in shifts:
        r.append(s * delta)
        D.append(float(np.mean((phz[:, s:] - phz[:, :-s]) ** 2)))
    return np.array(r), np.array(D)
