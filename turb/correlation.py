"""
Açıklık-ortalamalı tilt alanının UZAYSAL KORELASYON YAPISI.

Bu modül, C3 katkısının (türbülans-farkında ölçüm kovaryansı) matematiksel çekirdeğidir.

İki model tutulur ve makalede KARŞILAŞTIRILIR:

  (M1) "iso53" — proje raporu §4.2'deki mühendislik modeli
           rho(dtheta) = exp[-(dtheta/theta_c)^{5/3}] ,  Sigma = sigma^2 rho (x) I_2
       Yaygın kabul; iki varsayım içerir: (i) 5/3 gerilmiş üstel biçim,
       (ii) x/y eksenlerinin ayrışması (Kronecker (x) I_2).

  (M2) "physical" — Kolmogorov spektrumundan TÜRETİLEN tam tensör
           C_ab(dx) = (1/k^2) integral (2 pi f_a)(2 pi f_b) A(f)^2 Phi_phi(f)
                                        exp(2 pi i f . dx) d^2f
       A(f) = 2 J1(pi D f)/(pi D f) açıklık ortalaması. Bu, boyuna (longitudinal) ve
       enine (transverse) korelasyonların FARKLI olduğunu ve ayrımın eksene göre
       dönmesiyle C_xy != 0 doğduğunu verir.

ÖLÇÜLEN AYRIM (makalenin bulgusu): M2'nin kuyruğu M1'inkinden çok daha ağırdır —
4 açıklık ayrımında bile rho ~ 0.34 iken iso53 modeli ~0.02 verir. Yani yaygın
mühendislik modeli, uzak öznitelikler arasındaki korelasyonu ciddi biçimde
KÜÇÜMSER. E2 bu iki modeli ampirik rho_hat'e karşı sınar.

Korelasyon YALNIZCA dx/D'ye bağlıdır (r0'dan bağımsız): r0 yalnız genlik ölçeğidir.
Bu yüzden tek bir evrensel tablo yeterlidir ve önbelleğe alınır.
"""
from __future__ import annotations

import os
from functools import lru_cache

import numpy as np
from scipy.special import j1

_CACHE = os.path.join(os.path.dirname(__file__), "_tilt_corr_table.npz")


# =============================================================================
def _aperture_ft(f: np.ndarray, D: float) -> np.ndarray:
    x = np.pi * D * f
    out = np.ones_like(x)
    nz = np.abs(x) > 1e-12
    out[nz] = 2.0 * j1(x[nz]) / x[nz]
    return out


def _compute_table(n_sep: int = 241, sep_max: float = 24.0,
                   n_f: int = 4000, n_th: int = 721):
    """Boyuna/enine tilt korelasyonlarını dx/D ızgarasında sayısal integralle üret.

    Normalize edilmiş (rho(0) = 1) döndürülür; D = 1 alınır (ölçekten bağımsız).
    """
    D = 1.0
    f = np.logspace(-4, 3.5, n_f)                  # çevrim/D
    th = np.linspace(0.0, 2.0 * np.pi, n_th)
    A2 = _aperture_ft(f, D) ** 2
    # Phi ~ f^{-11/3}; sabit katsayılar normalizasyonda düşer
    W = (2 * np.pi * f) ** 2 * A2 * f ** (-11.0 / 3.0) * f      # son f: d^2f = f df dtheta
    cos2 = np.cos(th) ** 2
    sin2 = np.sin(th) ** 2
    seps = np.linspace(0.0, sep_max, n_sep)
    C_long = np.empty(n_sep)
    C_tran = np.empty(n_sep)
    for i, s in enumerate(seps):
        a = 2 * np.pi * np.outer(f, np.cos(th)) * s              # (n_f, n_th)
        ca = np.cos(a)
        gl = np.trapz(ca * cos2[None, :], th, axis=1)
        gt = np.trapz(ca * sin2[None, :], th, axis=1)
        C_long[i] = np.trapz(W * gl, f)
        C_tran[i] = np.trapz(W * gt, f)
    return seps, C_long / C_long[0], C_tran / C_tran[0]


@lru_cache(maxsize=1)
def _table():
    if os.path.exists(_CACHE):
        d = np.load(_CACHE)
        return d["seps"], d["c_long"], d["c_tran"]
    seps, cl, ct = _compute_table()
    try:
        np.savez_compressed(_CACHE, seps=seps, c_long=cl, c_tran=ct)
    except OSError:
        pass
    return seps, cl, ct


# =============================================================================
def rho_long(sep_over_D: np.ndarray) -> np.ndarray:
    """Boyuna (ayrım yönüne PARALEL bileşen) korelasyonu."""
    seps, cl, _ = _table()
    s = np.abs(np.asarray(sep_over_D, dtype=np.float64))
    return np.interp(s, seps, cl, right=cl[-1])


def rho_tran(sep_over_D: np.ndarray) -> np.ndarray:
    """Enine (ayrım yönüne DİK bileşen) korelasyonu."""
    seps, _, ct = _table()
    s = np.abs(np.asarray(sep_over_D, dtype=np.float64))
    return np.interp(s, seps, ct, right=ct[-1])


def tensor_blocks(dx: np.ndarray, dy: np.ndarray, D: float):
    """(N,N) ayrım dizileri için 2x2 korelasyon tensörünün bileşenleri.

    Ayrım birim vektörü u = (dx,dy)/|d| olmak üzere
        C = rho_L(u u^T) + rho_T(I - u u^T)
    Dönüş: (Cxx, Cxy, Cyy) — hepsi giriş dizileriyle aynı biçimde.
    """
    dx = np.asarray(dx, dtype=np.float64)
    dy = np.asarray(dy, dtype=np.float64)
    r = np.hypot(dx, dy)
    s = r / D
    rl = rho_long(s)
    rt = rho_tran(s)
    with np.errstate(invalid="ignore", divide="ignore"):
        ux = np.where(r > 1e-12, dx / np.maximum(r, 1e-30), 1.0)
        uy = np.where(r > 1e-12, dy / np.maximum(r, 1e-30), 0.0)
    Cxx = rl * ux * ux + rt * (1.0 - ux * ux)
    Cyy = rl * uy * uy + rt * (1.0 - uy * uy)
    Cxy = (rl - rt) * ux * uy
    return Cxx, Cxy, Cyy


def rho_iso53(dtheta: np.ndarray, theta_c: float, exponent: float = 5.0 / 3.0) -> np.ndarray:
    """Mühendislik modeli: exp[-(dtheta/theta_c)^{5/3}] (rapor §4.2)."""
    d = np.abs(np.asarray(dtheta, dtype=np.float64))
    with np.errstate(over="ignore"):
        return np.exp(-np.power(d / max(theta_c, 1e-30), exponent))


# =============================================================================
def covariance_matrix(theta_x: np.ndarray, theta_y: np.ndarray, sigma_px: float,
                      z_layer: float, D: float, model: str = "physical",
                      theta_c: float | None = None) -> np.ndarray:
    """Bir karedeki N öznitelik için 2N x 2N ölçüm kovaryansı Sigma_yama [px^2].

    Sıralama: [u_1, v_1, u_2, v_2, ...] (öznitelik-içi bloklar bitişik).

    Args:
        theta_x, theta_y: özniteliklerin açısal konumları [rad] (N,)
        sigma_px: tek eksen tilt jitteri [px]
        z_layer: eşdeğer katman mesafesi [m]  (ayrım dx = z * dtheta)
        D: açıklık [m]
        model: "physical" (tam tensör) | "iso53" (rapor modeli) | "diag"
        theta_c: iso53 için korelasyon açısı [rad]
    """
    tx = np.asarray(theta_x, dtype=np.float64).ravel()
    ty = np.asarray(theta_y, dtype=np.float64).ravel()
    N = tx.size
    S = np.zeros((2 * N, 2 * N), dtype=np.float64)
    if N == 0:
        return S
    dth_x = tx[:, None] - tx[None, :]
    dth_y = ty[:, None] - ty[None, :]

    if model == "diag":
        np.fill_diagonal(S[0::2, 0::2], 1.0)
        np.fill_diagonal(S[1::2, 1::2], 1.0)
        return S * sigma_px ** 2

    if model == "iso53":
        if theta_c is None:
            raise ValueError("iso53 modeli theta_c ister")
        R = rho_iso53(np.hypot(dth_x, dth_y), theta_c)
        S[0::2, 0::2] = R
        S[1::2, 1::2] = R
        return S * sigma_px ** 2

    # physical: tam 2x2 tensör
    Cxx, Cxy, Cyy = tensor_blocks(z_layer * dth_x, z_layer * dth_y, D)
    S[0::2, 0::2] = Cxx
    S[0::2, 1::2] = Cxy
    S[1::2, 0::2] = Cxy
    S[1::2, 1::2] = Cyy
    return S * sigma_px ** 2


def fit_rho_model(sep_theta: np.ndarray, rho_emp: np.ndarray,
                  exponent: float | None = 5.0 / 3.0):
    """Ampirik rho_hat'e exp[-(dtheta/theta_c)^p] uydur.

    exponent None ise p de serbest bırakılır. Dönüş: (theta_c, p, R2).
    """
    from scipy.optimize import curve_fit
    x = np.asarray(sep_theta, dtype=np.float64)
    y = np.asarray(rho_emp, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[m], y[m]
    if x.size < 3:
        return np.nan, np.nan, np.nan
    if exponent is None:
        f = lambda t, tc, p: np.exp(-np.power(t / tc, p))       # noqa: E731
        p0 = [np.median(x), 1.0]
    else:
        f = lambda t, tc: np.exp(-np.power(t / tc, exponent))   # noqa: E731
        p0 = [np.median(x)]
    try:
        popt, _ = curve_fit(f, x, y, p0=p0, maxfev=20000)
    except Exception:
        return np.nan, np.nan, np.nan
    res = y - f(x, *popt)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1.0 - np.sum(res ** 2) / ss_tot if ss_tot > 0 else np.nan
    tc = float(popt[0])
    p = float(popt[1]) if exponent is None else float(exponent)
    return tc, p, float(r2)


def r2_of_model(sep_theta: np.ndarray, rho_emp: np.ndarray, rho_model: np.ndarray) -> float:
    """Serbest parametresiz bir modelin açıkladığı varyans oranı."""
    y = np.asarray(rho_emp, dtype=np.float64)
    m = np.asarray(rho_model, dtype=np.float64)
    ok = np.isfinite(y) & np.isfinite(m)
    y, m = y[ok], m[ok]
    ss_tot = np.sum((y - y.mean()) ** 2)
    if ss_tot <= 0:
        return np.nan
    return float(1.0 - np.sum((y - m) ** 2) / ss_tot)
