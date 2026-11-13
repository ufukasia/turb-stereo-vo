"""DIMM v1 (fiziksel cekirdek) vs v2 (ampirik cekirdek) — SENTETIK KESIN SINAV.

GERCEK theta_c ve sigma bilinen veri uretilir; iki kestirici de ayni veriyi
gorur. Veri, E2'de OLCULEN ampirik cekirdekle uretilir (gercek simulatorun
davranisi budur): rho(d) = A exp(-(d/theta_c)^p), A=0.75, p=0.40.

    ./turbvio/bin/python -m turb._dimm_probe
"""
import numpy as np


def main():
    from .dimm import DifferentialMotionMonitor as V1
    from .dimm_v2 import DifferentialMotionMonitorV2 as V2
    from .covariance import EMP_A, EMP_P

    PX_PER_RAD, D, LAM = 170648.0, 0.10, 550e-9
    W = 752
    rng = np.random.default_rng(7)

    print(f"{'theta_c':>8} {'sigma':>7} | {'v1 tc':>8} {'v1 s':>7} {'v1 hata':>9} "
          f"| {'v2 tc':>8} {'v2 s':>7} {'v2 hata':>9}")
    print("-" * 82)
    for tc_true in (60.0, 112.8, 200.0):
        for s_true in (1.0, 2.2):
            d1 = V1(PX_PER_RAD, D, LAM, base_sigma_px=0.5, window_frames=60)
            d2 = V2(PX_PER_RAD, D, LAM, base_sigma_px=0.5, window_frames=60)
            for _ in range(60):
                n = 120
                uv = rng.uniform(0, W, size=(n, 2))
                dd = np.hypot(uv[:, 0][:, None] - uv[:, 0][None, :],
                              uv[:, 1][:, None] - uv[:, 1][None, :])
                R = EMP_A * np.exp(-((dd / tc_true) ** EMP_P))
                np.fill_diagonal(R, 1.0)
                S = (s_true ** 2) * R + 1e-9 * np.eye(n)
                L = np.linalg.cholesky(S)
                res = np.column_stack([L @ rng.standard_normal(n),
                                       L @ rng.standard_normal(n)])
                d1.add_frame(uv, res)
                d2.add_frame(uv, res)
            s1, t1, _ = d1.estimate()
            s2, t2, _ = d2.estimate()
            e1 = abs(t1 - tc_true) / tc_true if np.isfinite(t1) else float("nan")
            e2 = abs(t2 - tc_true) / tc_true if np.isfinite(t2) else float("nan")
            print(f"{tc_true:>8.1f} {s_true:>7.2f} | {t1:>8.1f} {s1:>7.2f} "
                  f"{100*e1:>8.0f}% | {t2:>8.1f} {s2:>7.2f} {100*e2:>8.0f}%")
    print("\nv1'in alt siniri (tc_lo) = 33.3 px — bu degere yapisan her satir")
    print("kestirimin ARANMADIGINI, SINIRIN raporlandigini gosterir.")


if __name__ == "__main__":
    main()
