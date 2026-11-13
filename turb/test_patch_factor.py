"""
YAMA FAKTORU DOGRULAMASI — Sigma dogru oldugu KANITLANDIKTAN sonraki halka.

test_cov_model gosterdi ki patch_covariance'in urettigi Sigma, simulatorun gercek
yer degistirme istatistigiyle uyusuyor (Mahalanobis orani 1.04; skaler model 1.27).
Yani kusur modelde DEGIL. Sirada faktor var: Sigma dogru olsa bile faktor onu yanlis
kullaniyor olabilir.

SORU 1 (whitening): faktorun hatasi gercekten Sigma^{-1/2} r midir?
  GTSAM'de noiseModel.Gaussian.Covariance(S) whitened hatayi verir; unwhitened
  artik r icin  e_w = R r  olmali (R^T R = S^-1).
SORU 2 (Jakobiyen): analitik Jakobiyen sayisal turevle ayni mi?
SORU 3 (ölçek): tek bir gozlem icin yama faktoru, standart GTSAM faktoruyle AYNI
  ki-kare degerini vermeli (yapı yokken ikisi ozdes olmali).

    ./turbvio/bin/python -m turb.test_patch_factor
"""
from __future__ import annotations

import numpy as np
import gtsam
from gtsam import Pose3, Rot3, Point3
from gtsam.symbol_shorthand import L, X

from .factors import make_patch_factor


def _setup(n=3, stereo=False, seed=0):
    rng = np.random.default_rng(seed)
    fx = fy = 400.0
    cx, cy = 376.0, 240.0
    base = 0.11
    pose = Pose3(Rot3.Rodrigues(0.02, -0.01, 0.03), Point3(0.1, -0.2, 0.05))
    lms = [np.array(rng.normal(0, 0.5, 3) + np.array([0.0, 0.0, 5.0]), dtype=float)
           for _ in range(n)]
    keys_l = [L(i) for i in range(n)]
    is_st = np.array([stereo] * n, dtype=bool)
    dim = sum(3 if s else 2 for s in is_st)
    A = rng.normal(size=(dim, dim))
    S = A @ A.T + dim * np.eye(dim) * 0.3          # rastgele ama iyi kosullu SPD
    meas = []
    Tcw = np.linalg.inv(pose.matrix())
    for i, p in enumerate(lms):
        pc = Tcw @ np.array([p[0], p[1], p[2], 1.0])
        u = fx * pc[0] / pc[2] + cx
        v = fy * pc[1] / pc[2] + cy
        if is_st[i]:
            meas += [u, u - fx * base / pc[2], v]
        else:
            meas += [u, v]
    meas = np.array(meas) + rng.normal(0, 0.3, dim)   # kucuk artik
    return pose, lms, keys_l, meas, is_st, S, (fx, fy, cx, cy, base)


def run(n=3, stereo=False, seed=0):
    pose, lms, keys_l, meas, is_st, S, K = _setup(n, stereo, seed)
    fx, fy, cx, cy, base = K
    fac, errfn = make_patch_factor(X(0), keys_l, meas, is_st, S, fx, fy, cx, cy, base)

    vals = gtsam.Values()
    vals.insert(X(0), pose)
    for k, p in zip(keys_l, lms):
        vals.insertPoint3(k, p)      # insert() dinamik vektor saklar -> atPoint3 patlar

    r = np.asarray(errfn(None, vals, None))   # unwhitened artik
    # GTSAM faktorunun ki-kare degeri = 0.5 * r^T S^-1 r olmali
    chi2_gtsam = 2.0 * fac.error(vals)
    chi2_ref = float(r @ np.linalg.solve(S, r))
    d1 = abs(chi2_gtsam - chi2_ref) / max(abs(chi2_ref), 1e-12)

    # --- Jakobiyen: analitik vs sayisal ------------------------------------
    lin = fac.linearize(vals)
    # JacobianFactor'un tam A matrisi; sutunlar anahtar sirasina gore dilimlenir.
    A_all = np.asarray(lin.jacobian()[0])
    keys_order = list(lin.keys())
    dims = [6 if k == X(0) else 3 for k in keys_order]
    offs = np.cumsum([0] + dims)
    Js = {k: A_all[:, offs[i]:offs[i + 1]] for i, k in enumerate(keys_order)}
    # sayisal
    eps = 1e-6
    num = {}
    d = len(r)
    for k in [X(0)] + list(keys_l):
        dim_k = 6 if k == X(0) else 3
        Jn = np.zeros((d, dim_k))
        for j in range(dim_k):
            dv = np.zeros(dim_k); dv[j] = eps
            v2 = gtsam.Values(vals)
            if k == X(0):
                v2.update(k, vals.atPose3(k).retract(dv))
            else:
                v2.update(k, gtsam.Point3(*(vals.atPoint3(k) + dv)))
            Jn[:, j] = (np.asarray(errfn(None, v2, None)) - r) / eps
        num[k] = Jn
    # GTSAM'in A'si WHITENED oldugundan karsilastirma icin beyazlatilir
    Rw = np.linalg.cholesky(np.linalg.inv(S)).T
    errs = {}
    for k in Js:
        errs[k] = float(np.max(np.abs(Js[k] - Rw @ num[k])))
    return d1, errs, chi2_ref


def main():
    print("SORU 1+2: whitening ve Jakobiyen (rastgele SPD Sigma ile)")
    worst_chi, worst_j = 0.0, 0.0
    for stereo in (False, True):
        for n in (1, 3, 5):
            for seed in (0, 1):
                d1, errs, chi2 = run(n, stereo, seed)
                worst_chi = max(worst_chi, d1)
                worst_j = max(worst_j, max(errs.values()))
                print(f"  n={n} stereo={int(stereo)} seed={seed}: "
                      f"chi2 bagil fark={d1:.2e}  max|J_analitik-J_sayisal|="
                      f"{max(errs.values()):.2e}")
    print(f"\n  EN KOTU: chi2 {worst_chi:.2e}   Jakobiyen {worst_j:.2e}")

    print("\nSORU 3: tek MONO gozlem, kosegen Sigma -> standart GTSAM faktoruyle ayni mi")
    pose, lms, keys_l, meas, is_st, _S, K = _setup(1, False, 3)
    fx, fy, cx, cy, base = K
    sig = 0.7
    S = np.eye(2) * sig ** 2
    fac, _ = make_patch_factor(X(0), keys_l, meas, is_st, S, fx, fy, cx, cy, base)
    ref = gtsam.GenericProjectionFactorCal3_S2(
        np.array(meas[:2]), gtsam.noiseModel.Isotropic.Sigma(2, sig),
        X(0), keys_l[0], gtsam.Cal3_S2(fx, fy, 0.0, cx, cy), False, False)
    vals = gtsam.Values(); vals.insert(X(0), pose); vals.insertPoint3(keys_l[0], lms[0])
    e1, e2 = fac.error(vals), ref.error(vals)
    print(f"  yama faktoru {e1:.10f}  vs  GTSAM {e2:.10f}  "
          f"bagil fark {abs(e1-e2)/max(e2,1e-12):.2e}")

    ok = worst_chi < 1e-8 and worst_j < 1e-4 and abs(e1 - e2) / max(e2, 1e-12) < 1e-8
    print(f"\n{'✓ FAKTÖR DOĞRU' if ok else '✗ FAKTÖRDE KUSUR'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
