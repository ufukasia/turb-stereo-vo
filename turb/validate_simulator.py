"""
Simülatör doğrulaması — makalenin "simülasyon hattı fiziksel olarak doğrudur"
kanıtı (rapor §12: "Hepsi sentetik" ve "Kovaryans modeliniz uydurma" itirazlarına
karşı ilk savunma hattı).

Dört bağımsız sınama (hiçbiri simülatörün kendi varsayımını tekrar etmez):

  V1  Faz yapı fonksiyonu      D_phi(r)  =? 6.88 (r/r0)^{5/3}
  V2  Tilt varyansı            sigma_a^2 =? 0.170 (D/r0)^{5/3} (lam/D)^2   [G-tilt]
  V3  Uzaysal tilt korelasyonu rho_hat(dx) =? Kolmogorov tensör integrali
                               + mühendislik modeli exp[-(dtheta/theta_c)^{5/3}] uyumu
  V4  Zamansal korelasyon      donmuş akışta rho(dt) ve tau0 = 0.314 r0/v tutarlılığı

Çıktı: evidence/simulator_validation/ altında JSON + PNG.

    python -m turb.validate_simulator --n_seeds 16
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from . import optics, screen as _screen, correlation as _corr
from .fields import TiltField
from .simulator import TurbParams, TurbulenceSimulator


def _json_safe(o):
    if isinstance(o, dict):
        return {k: _json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_safe(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


# ---------------------------------------------------------------- V1
def v1_structure_function(N, ppd, nsh, D, r0, n_seeds, out):
    delta = D / ppd
    shifts = [1, 2, 4, 8, 16, 32, 64, 128]
    acc = np.zeros(len(shifts))
    for s in range(n_seeds):
        ph = _screen.phase_screen(N, delta, r0, np.random.default_rng(1000 + s),
                                  n_subharmonics=nsh)
        r, Dp = _screen.structure_function_shift(ph, delta, shifts)
        acc += Dp / n_seeds
    r = np.array(shifts) * delta
    th = 6.88 * (r / r0) ** (5.0 / 3.0)
    out["V1"] = {"r_m": r, "D_emp": acc, "D_theory": th, "ratio": acc / th,
                 "ratio_median": float(np.median(acc / th))}
    return r, acc, th


# ---------------------------------------------------------------- V2
def v2_tilt_variance(N, ppd, nsh, D, lam, n_seeds, out, D_over_r0_list=(1, 3, 5, 12)):
    rows = []
    for dr in D_over_r0_list:
        r0 = D / dr
        th = optics.sigma_tilt_rad(dr, lam, D, "gtilt")
        vals = []
        for s in range(n_seeds):
            tf = TiltField(N, D / ppd, r0, D, lam, np.random.default_rng(2000 + s),
                           n_subharmonics=nsh)
            vals.append(tf.tilt_std())
        vals = np.array(vals)
        rows.append({"D_over_r0": dr, "sigma_theory_rad": th,
                     "sigma_sim_rad_mean": float(vals.mean()),
                     "sigma_sim_rad_std": float(vals.std()),
                     "ratio": float(vals.mean() / th)})
    ratios = np.array([r["ratio"] for r in rows])
    out["V2"] = {"rows": rows, "ratio_mean": float(ratios.mean()),
                 "calibration_constant": float(1.0 / ratios.mean()),
                 "note": "calibration_constant, TurbParams.tilt_calibration olarak kullanılır"}
    return rows


# ---------------------------------------------------------------- V3
def v3_spatial_correlation(N, ppd, nsh, D, lam, n_seeds, out, r0=0.02, max_sep_D=6.0):
    delta = D / ppd
    kmax = int(max_sep_D * D / delta)
    seps = np.arange(0, kmax + 1)
    cl = np.zeros(len(seps))
    ct = np.zeros(len(seps))
    for s in range(n_seeds):
        tf = TiltField(N, delta, r0, D, lam, np.random.default_rng(3000 + s),
                       n_subharmonics=nsh)
        ax = tf.alpha_x - tf.alpha_x.mean()
        ay = tf.alpha_y - tf.alpha_y.mean()
        vx = np.mean(ax * ax)
        vy = np.mean(ay * ay)
        for i, k in enumerate(seps):
            if k == 0:
                cl[i] += 1.0 / n_seeds
                ct[i] += 1.0 / n_seeds
            else:
                # x yönünde ayrım: alpha_x boyuna, alpha_y enine
                cl[i] += np.mean(ax[:, k:] * ax[:, :-k]) / vx / n_seeds
                ct[i] += np.mean(ay[:, k:] * ay[:, :-k]) / vy / n_seeds
    sep_D = seps * delta / D
    th_l = _corr.rho_long(sep_D)
    th_t = _corr.rho_tran(sep_D)

    # mühendislik modeli uyumu (H-2 kabul ölçütü): dtheta = dx / z ; z serbest ölçek
    # olduğundan theta_c uydurulur -> biçim sınanır
    tc53, _, r2_53 = _corr.fit_rho_model(sep_D[1:], cl[1:], exponent=5.0 / 3.0)
    tcfree, pfree, r2_free = _corr.fit_rho_model(sep_D[1:], cl[1:], exponent=None)
    r2_phys = _corr.r2_of_model(sep_D[1:], cl[1:], th_l[1:])

    out["V3"] = {
        "sep_over_D": sep_D, "rho_long_sim": cl, "rho_tran_sim": ct,
        "rho_long_theory": th_l, "rho_tran_theory": th_t,
        "R2_physical_model": r2_phys,
        "R2_iso53_model": r2_53, "theta_c_over_D_iso53": tc53,
        "R2_free_exponent": r2_free, "free_exponent": pfree,
        "theta_c_over_D_free": tcfree,
    }
    return sep_D, cl, ct, th_l, th_t


# ---------------------------------------------------------------- V4
def v4_temporal(D, lam, out, D_over_r0=5.0, tau_ratios=(0.1, 0.5, 1.0, 2.0, 5.0),
                n_frames=40, seed=7):
    rows = []
    for tr in tau_ratios:
        p = TurbParams(enabled=True, seed=seed, D_over_r0=D_over_r0,
                       theta_c_over_fov=0.15, tau0_over_dt=tr, blur="none",
                       aperture_D=D, wavelength=lam)
        sim = TurbulenceSimulator(p)
        # merkez pikselin zamansal tilt izi
        u = np.array([sim.W / 2.0]); v = np.array([sim.H / 2.0])
        tr_series = np.array([sim.displacement_at(f, 0, u, v)[0][0] for f in range(n_frames)])
        tr_series -= tr_series.mean()
        denom = np.mean(tr_series ** 2)
        ac = [float(np.mean(tr_series[k:] * tr_series[:len(tr_series) - k]) / denom)
              for k in range(min(20, n_frames // 2))]
        # 1/e düşüş süresi [kare]
        te = next((k for k, a in enumerate(ac) if a < np.exp(-1.0)), np.nan)
        rows.append({"tau0_over_dt": tr, "tau0_s": sim.tau0, "v_wind": sim.v_wind,
                     "autocorr": ac, "decay_frames_1e": te})
    out["V4"] = {"rows": rows}
    return rows


# ---------------------------------------------------------------- çizim
def _plots(out, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    a = axes[0, 0]
    v1 = out["V1"]
    a.loglog(v1["r_m"], v1["D_emp"], "o-", label="simülasyon")
    a.loglog(v1["r_m"], v1["D_theory"], "k--", label=r"$6.88\,(r/r_0)^{5/3}$")
    a.set_xlabel("ayrım r [m]"); a.set_ylabel(r"$D_\varphi$ [rad$^2$]")
    a.set_title(f"V1 faz yapı fonksiyonu (medyan oran {v1['ratio_median']:.3f})")
    a.legend(); a.grid(True, which="both", alpha=.3)

    a = axes[0, 1]
    rows = out["V2"]["rows"]
    x = [r["D_over_r0"] for r in rows]
    a.plot(x, [r["sigma_sim_rad_mean"] * 1e6 for r in rows], "o-", label="simülasyon")
    a.plot(x, [r["sigma_theory_rad"] * 1e6 for r in rows], "k--",
           label=r"$\sqrt{0.170}\,(D/r_0)^{5/6}\lambda/D$")
    a.set_xlabel(r"$D/r_0$"); a.set_ylabel(r"$\sigma_\alpha$ [$\mu$rad]")
    a.set_title(f"V2 tilt varyansı (oran {out['V2']['ratio_mean']:.3f})")
    a.legend(); a.grid(alpha=.3)

    a = axes[1, 0]
    v3 = out["V3"]
    s = np.asarray(v3["sep_over_D"])
    a.plot(s, v3["rho_long_sim"], "o", ms=3, label="sim. boyuna")
    a.plot(s, v3["rho_tran_sim"], "s", ms=3, label="sim. enine")
    a.plot(s, v3["rho_long_theory"], "k-", label="Kolmogorov tensör (boyuna)")
    a.plot(s, v3["rho_tran_theory"], "k:", label="Kolmogorov tensör (enine)")
    a.plot(s, np.exp(-(s / v3["theta_c_over_D_iso53"]) ** (5 / 3)), "r--",
           label=r"iso53: $e^{-(\Delta\theta/\theta_c)^{5/3}}$")
    a.set_xlabel(r"ayrım $\Delta x / D$"); a.set_ylabel(r"$\rho$")
    a.set_title(f"V3 uzaysal korelasyon  (R²: fizik {v3['R2_physical_model']:.3f} / "
                f"iso53 {v3['R2_iso53_model']:.3f})")
    a.legend(fontsize=7); a.grid(alpha=.3)

    a = axes[1, 1]
    for r in out["V4"]["rows"]:
        a.plot(r["autocorr"], "o-", ms=3, label=f"$\\tau_0/\\Delta t$={r['tau0_over_dt']}")
    a.axhline(np.exp(-1), color="k", ls=":", lw=.8)
    a.set_xlabel("kare gecikmesi"); a.set_ylabel("tilt otokorelasyonu")
    a.set_title("V4 zamansal korelasyon (donmuş akış)")
    a.legend(fontsize=7); a.grid(alpha=.3)

    fig.tight_layout()
    p = os.path.join(outdir, "simulator_validation.png")
    fig.savefig(p, dpi=150)
    plt.close(fig)
    return p


# ---------------------------------------------------------------- ana
def main():
    ap = argparse.ArgumentParser(description="Türbülans simülatörü doğrulaması")
    ap.add_argument("--n_seeds", type=int, default=12)
    ap.add_argument("--screen_N", type=int, default=2048)
    ap.add_argument("--px_per_D", type=int, default=8)
    ap.add_argument("--subharmonics", type=int, default=3)
    ap.add_argument("--D", type=float, default=0.10)
    ap.add_argument("--lam", type=float, default=550e-9)
    ap.add_argument("--outdir", default="evidence/simulator_validation")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out = {"config": vars(args), "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}

    t0 = time.time()
    print("V1 faz yapı fonksiyonu ...")
    v1_structure_function(args.screen_N, args.px_per_D, args.subharmonics,
                          args.D, args.D / 5.0, args.n_seeds, out)
    print(f"   medyan oran = {out['V1']['ratio_median']:.4f}   ({time.time()-t0:.0f}s)")

    t0 = time.time()
    print("V2 tilt varyansı ...")
    v2_tilt_variance(args.screen_N, args.px_per_D, args.subharmonics,
                     args.D, args.lam, args.n_seeds, out)
    print(f"   oran = {out['V2']['ratio_mean']:.4f}  -> kalibrasyon "
          f"{out['V2']['calibration_constant']:.4f}   ({time.time()-t0:.0f}s)")

    t0 = time.time()
    print("V3 uzaysal korelasyon ...")
    v3_spatial_correlation(args.screen_N, args.px_per_D, args.subharmonics,
                           args.D, args.lam, args.n_seeds, out)
    print(f"   R2 fizik={out['V3']['R2_physical_model']:.4f}  "
          f"iso53={out['V3']['R2_iso53_model']:.4f}  "
          f"serbest üs={out['V3']['free_exponent']:.2f} (R2={out['V3']['R2_free_exponent']:.4f})"
          f"   ({time.time()-t0:.0f}s)")

    t0 = time.time()
    print("V4 zamansal korelasyon ...")
    v4_temporal(args.D, args.lam, out)
    for r in out["V4"]["rows"]:
        print(f"   tau0/dt={r['tau0_over_dt']:<4} -> 1/e düşüşü {r['decay_frames_1e']} kare")
    print(f"   ({time.time()-t0:.0f}s)")

    png = _plots(out, args.outdir)
    with open(os.path.join(args.outdir, "validation.json"), "w") as f:
        json.dump(_json_safe(out), f, indent=2)
    print(f"\n✓ {png}\n✓ {os.path.join(args.outdir,'validation.json')}")


if __name__ == "__main__":
    main()
