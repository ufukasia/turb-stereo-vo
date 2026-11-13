"""
YAMA KOVARYANSI MODEL DOĞRULAMASI — SLAM'den bağımsız, dakikalar süren kesin test.

SORU: `turb.covariance.patch_covariance`'ın ürettigi Sigma, simulatorun o
oznitelik konumlarinda URETTIGI GERCEK yer degistirmelerin kovaryansi midir?

Bu, tum zincirin en kritik halkasidir. Sigma yanlissa faktor dogru whitening
yapmaz ve NEES patlar — ki olculdu: ayni sigma olcegiyle skaler sisirme (A1)
NEES 32 verirken yapisal kovaryans (A5) 254 veriyor.

YONTEM (dairesel degil):
  1. Simulatorden BAGIMSIZ gerceklemeler uret.
  2. Her gerceklemede, sabit bir oznitelik kumesinin GERCEK yer degistirmesini
     olc (simulatorun displacement alani).
  3. Bu orneklerin AMPIRIK kovaryansini hesapla.
  4. patch_covariance'in ONGORDUGU Sigma ile karsilastir (taban gurultu haric).
  5. Ayrica NEES benzeri bir olcut: d^T Sigma^-1 d ortalamasi = boyut olmali.

    ./turbvio/bin/python -m turb.test_cov_model --D_over_r0 5 --n_real 300
"""
from __future__ import annotations

import argparse

import numpy as np

from . import covariance as _cov
from . import optics as _optics
from .simulator import TurbParams, TurbulenceSimulator


def main():
    ap = argparse.ArgumentParser(description="Yama kovaryansı model doğrulaması")
    ap.add_argument("--D_over_r0", type=float, default=5.0)
    ap.add_argument("--theta_c_over_fov", type=float, default=0.15)
    ap.add_argument("--n_real", type=int, default=200)
    ap.add_argument("--n_feat", type=int, default=5)
    ap.add_argument("--image_width", type=int, default=752)
    ap.add_argument("--shrinkage", type=float, default=0.05)
    ap.add_argument("--blur", default="none",
                    help="none | speckle_global — bagimsiz bilesenin etkisi")
    ap.add_argument("--sigma_ind", type=float, default=0.0,
                    help="modele eklenecek BAGIMSIZ gurultu tabani [px]")
    a = ap.parse_args()

    W = a.image_width
    px_per_rad = 170648.0
    D = 0.10
    theta_c = a.theta_c_over_fov * (W / px_per_rad)
    z = _optics.layer_distance_for_theta_c(D, theta_c)
    sigma_px = _optics.sigma_tilt_px(a.D_over_r0, 550e-9, D, px_per_rad)
    tc_px = theta_c * px_per_rad
    print(f"D/r0={a.D_over_r0:g}  sigma_teori={sigma_px:.3f} px  "
          f"theta_c={tc_px:.1f} px  z_layer={z:.1f} m  n_feat={a.n_feat}")

    # Bir YAMA icindeki tipik oznitelik dagilimi (hucre boyu = theta_c)
    rng = np.random.default_rng(0)
    uv = rng.uniform(0, tc_px, size=(a.n_feat, 2)) + np.array([W / 3, 200.0])
    theta = uv / px_per_rad

    # --- 1) GERCEK yer degistirmeleri bagimsiz gerceklemelerden ornekle -------
    samples = []
    for k in range(a.n_real):
        p = TurbParams(enabled=True, seed=k, D_over_r0=a.D_over_r0,
                       theta_c_over_fov=a.theta_c_over_fov, blur=a.blur,
                       tilt=True, px_per_rad=px_per_rad)
        sim = TurbulenceSimulator(p, image_size=(W, 480))
        dx, dy = sim.displacement_at(0, 0, uv[:, 0], uv[:, 1])   # iki AYRI dizi
        # SIRA KRİTİK: patch_covariance blok yerleşimi [u1,v1,u2,v2,...] bekler.
        # `displacement_at` (dx, dy) DEMETİ döndürür; düz reshape [dx..,dy..] verir
        # ve bileşenleri karıştırır. Bu hata testin ilk sürümünde modeli haksız yere
        # "uyuşmuyor" gösterdi — testin kendisi de doğrulanmalıdır.
        samples.append(np.stack([dx, dy], axis=1).reshape(-1))
    Xs = np.array(samples)                                          # (n_real, 2N)
    Xs = Xs - Xs.mean(axis=0, keepdims=True)
    S_emp = Xs.T @ Xs / (len(Xs) - 1)

    # --- 2) MODELIN ongordugu Sigma (taban gurultu HARIC, shrinkage HARIC) ----
    S_mod_full = _cov.patch_covariance(
        theta, np.zeros(a.n_feat, dtype=bool), sigma_px,
        base_sigma_mono=a.sigma_ind, base_sigma_stereo=a.sigma_ind, z_layer=z, D=D,
        stereo_baseline=0.11, model="physical", theta_c=theta_c, shrinkage=0.0)

    # --- 3) Karsilastir --------------------------------------------------------
    de = np.sqrt(np.diag(S_emp)); dm = np.sqrt(np.diag(S_mod_full))
    print(f"\nkosegen sigma [px]  ampirik medyan {np.median(de):.3f}  "
          f"model {np.median(dm):.3f}  oran {np.median(de)/np.median(dm):.3f}")
    Re = S_emp / np.outer(de, de)
    Rm = S_mod_full / np.outer(dm, dm)
    iu = np.triu_indices_from(Re, 1)
    print(f"kosegen-disi korelasyon: ampirik medyan {np.median(Re[iu]):+.3f}  "
          f"model {np.median(Rm[iu]):+.3f}")
    print(f"korelasyon farki: ortalama |fark| {np.mean(np.abs(Re[iu]-Rm[iu])):.3f}  "
          f"max {np.max(np.abs(Re[iu]-Rm[iu])):.3f}")

    # --- 4) MAHALANOBIS testi: d^T Sigma^-1 d ortalamasi = 2N olmali ----------
    print("\nMahalanobis (E[d^T Sigma^-1 d] = boyut = %d olmali):" % (2 * a.n_feat))
    for name, S in (("model (tam yapi)", S_mod_full),
                    ("model + shrinkage %.2f" % a.shrinkage,
                     _cov.patch_covariance(theta, np.zeros(a.n_feat, bool), sigma_px,
                                           0.0, 0.0, z, D, 0.11, model="physical",
                                           theta_c=theta_c, shrinkage=a.shrinkage)),
                    ("SKALER (kosegen, ayni sigma)",
                     np.eye(2 * a.n_feat) * (sigma_px ** 2 + a.sigma_ind ** 2)),
                    ("ampirik (ust sinir)", S_emp)):
        try:
            Si = np.linalg.inv(S + 1e-12 * np.eye(len(S)))
            m = float(np.mean(np.einsum("ij,jk,ik->i", Xs, Si, Xs)))
            print(f"  {name:34s} {m:10.2f}   (oran {m/(2*a.n_feat):6.2f})")
        except np.linalg.LinAlgError:
            print(f"  {name:34s} tekil")

    ok = (abs(np.median(de) / np.median(dm) - 1) < 0.25
          and np.mean(np.abs(Re[iu] - Rm[iu])) < 0.15)
    print(f"\n{'✓ MODEL TUTARLI' if ok else '✗ MODEL UYUŞMUYOR'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
