"""
E2 analizi — ön-uç karakterizasyonu ve korelasyon modeli doğrulaması.

    python -m turb.analyze_e2 --results evidence/E2/e2_results.json

Üretir: e2_summary.csv, e2_correlation.png, e2_frontend.png, e2_report.md
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd


def load(path):
    with open(path) as f:
        rows = json.load(f)
    return rows


def summarize(rows):
    recs = []
    for r in rows:
        recs.append({k: r[k] for k in (
            "sequence", "D_over_r0", "theta_c_over_fov", "seed", "blur",
            "sigma_px_theory", "theta_c_px", "meas_sigma_px", "sigma_tilt_implied_px",
            "A_physical", "A_iso53", "pos_err_rms_px", "pos_err_median_px",
            "gt_disp_rms_px", "repeatability", "matching_score", "desc_score",
            "n_kp_clean", "n_kp_dist", "R2_physical", "R2_iso53",
            "R2_physical_fitted", "R2_iso53_fitted", "R2_free", "free_exponent",
            "theta_c_px_physical_fit", "theta_c_px_iso53_fit") if k in r})
    return pd.DataFrame(recs)


def figures(rows, df, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # --- Şekil 1: ampirik korelasyon vs modeller ---
    sel = [r for r in rows if abs(r["theta_c_over_fov"] - 0.15) < 1e-9]
    if not sel:
        sel = rows
    drs = sorted({r["D_over_r0"] for r in sel})[:4]
    fig, axes = plt.subplots(1, len(drs), figsize=(4 * len(drs), 3.4), squeeze=False)
    for ax, dr in zip(axes[0], drs):
        g = [r for r in sel if r["D_over_r0"] == dr]
        sep = np.array(g[0]["sep_px"])
        emp = np.mean([np.array(r["rho_long_emp"]) for r in g], axis=0)
        cnt = np.mean([np.array(r["bin_counts"]) for r in g], axis=0)
        A = float(np.mean([r["A_physical"] for r in g]))
        m = cnt >= 50
        ax.plot(sep[m], emp[m], "o", ms=4, label=r"ampirik $\hat\rho$ (ölçülen)")
        ax.plot(sep[m], A * np.array(g[0]["rho_long_model_phys"])[m], "-",
                label=f"fizik tensörü (A={A:.2f})")
        Ai = float(np.mean([r["A_iso53"] for r in g]))
        ax.plot(sep[m], Ai * np.array(g[0]["rho_model_iso53"])[m], "--",
                label=r"iso53: $e^{-(\Delta\theta/\theta_c)^{5/3}}$")
        r2p = np.mean([r["R2_physical"] for r in g])
        r2i = np.mean([r["R2_iso53"] for r in g])
        ax.set_title(f"$D/r_0$={dr:g}   R²: fizik {r2p:.2f} / iso53 {r2i:.2f}",
                     fontsize=9)
        ax.set_xlabel("açısal ayrım [px]")
        ax.grid(alpha=.3)
    axes[0][0].set_ylabel(r"$\rho$")
    axes[0][0].legend(fontsize=7)
    fig.tight_layout()
    p1 = os.path.join(outdir, "e2_correlation.png")
    fig.savefig(p1, dpi=150); plt.close(fig)

    # --- Şekil 2: ön-uç metrikleri ---
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
    g = df.groupby("D_over_r0")
    axes[0].errorbar(g["meas_sigma_px"].mean().index, g["sigma_tilt_implied_px"].mean(),
                     yerr=g["sigma_tilt_implied_px"].std(), fmt="o-", capsize=3,
                     label="ölçülen (gürültü düzeltmeli)")
    axes[0].plot(g["sigma_px_theory"].mean().index, g["sigma_px_theory"].mean(), "k--",
                 label=r"teori $\sqrt{0.170}(D/r_0)^{5/6}\lambda/D$")
    axes[0].set_xlabel(r"$D/r_0$"); axes[0].set_ylabel(r"$\sigma_{tilt}$ [px]")
    axes[0].set_title("tilt genliği"); axes[0].legend(fontsize=7); axes[0].grid(alpha=.3)

    axes[1].errorbar(g["repeatability"].mean().index, g["repeatability"].mean(),
                     yerr=g["repeatability"].std(), fmt="o-", capsize=3)
    axes[1].set_xlabel(r"$D/r_0$"); axes[1].set_ylabel("tekrarlanabilirlik")
    axes[1].set_title("dedektör tekrarlanabilirliği"); axes[1].grid(alpha=.3)

    axes[2].errorbar(g["matching_score"].mean().index, g["matching_score"].mean(),
                     yerr=g["matching_score"].std(), fmt="o-", capsize=3,
                     label="eşleştirme skoru")
    axes[2].errorbar(g["pos_err_rms_px"].mean().index, g["pos_err_rms_px"].mean(),
                     yerr=g["pos_err_rms_px"].std(), fmt="s-", capsize=3,
                     label="konum hatası [px]")
    axes[2].set_xlabel(r"$D/r_0$"); axes[2].set_title("eşleştirme ve konum hatası")
    axes[2].legend(fontsize=7); axes[2].grid(alpha=.3)
    fig.tight_layout()
    p2 = os.path.join(outdir, "e2_frontend.png")
    fig.savefig(p2, dpi=150); plt.close(fig)
    return p1, p2


def main():
    ap = argparse.ArgumentParser(description="E2 analizi")
    ap.add_argument("--results", default="evidence/E2/e2_results.json")
    a = ap.parse_args()
    outdir = os.path.dirname(a.results)
    rows = load(a.results)
    df = summarize(rows)
    df.to_csv(os.path.join(outdir, "e2_summary.csv"), index=False)
    p1, p2 = figures(rows, df, outdir)

    g = df.groupby("D_over_r0").agg(
        sigma_theory=("sigma_px_theory", "mean"),
        sigma_measured=("sigma_tilt_implied_px", "mean"),
        ratio=("sigma_tilt_implied_px", lambda s: np.nan),
        A=("A_physical", "mean"),
        pos_err=("pos_err_rms_px", "mean"),
        repeat=("repeatability", "mean"),
        match=("matching_score", "mean"),
        R2_phys=("R2_physical", "mean"),
        R2_iso53=("R2_iso53", "mean"),
        R2_phys_fit=("R2_physical_fitted", "mean"),
        R2_iso53_fit=("R2_iso53_fitted", "mean"),
        exponent=("free_exponent", "mean"))
    g["ratio"] = g["sigma_measured"] / g["sigma_theory"]

    ok_phys = float(np.mean(df["R2_physical"] >= 0.85))
    ok_iso = float(np.mean(df["R2_iso53"] >= 0.85))
    rep = ["# E2 — ön-uç karakterizasyonu ve korelasyon modeli doğrulaması", "",
           f"Hücre sayısı: {len(df)}  (sekans × D/r₀ × θc × tohum)", "",
           "## H-2 kabul ölçütü (R² ≥ 0.85)", "",
           f"- **Fizik tensörü** (turb/correlation.py): hücrelerin **%{100*ok_phys:.0f}**'inde sağlandı",
           f"- **Rapor §4.2 modeli** exp[−(Δθ/θc)^{{5/3}}]: hücrelerin **%{100*ok_iso:.0f}**'inde sağlandı",
           f"- Serbest üs uydurması ortalama üs = **{df['free_exponent'].mean():.2f}** "
           f"(rapor 5/3 = 1.67 varsayıyordu)", "",
           "## D/r₀ ekseninde özet", "", g.round(3).to_markdown(), "",
           f"![korelasyon]({os.path.basename(p1)})", "",
           f"![ön-uç]({os.path.basename(p2)})", ""]
    with open(os.path.join(outdir, "e2_report.md"), "w") as f:
        f.write("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
