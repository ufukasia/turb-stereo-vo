"""
Kanıt analizi — JSONL kayıtlarından tablo ve şekil üretimi (rapor §6-§7).

    python -m turb.analyze --exp E1a
    python -m turb.analyze --exp E4

Üretilenler (evidence/<EXP>/ altına):
    summary.csv        hücre bazlı özet (ortalama + %95 bootstrap GA)
    tests.json         eşleştirilmiş üçlü ölçüt sonuçları (rapor §6.0)
    *.png              şekiller
    report.md          insan-okur özet

EŞLEŞTİRME BİRİMİ: (sekans, pencere, tohum). Temiz koşum belirlenimci olduğundan
her (sekans, pencere) için tek bir referans değeri vardır ve tüm tohumlarla
eşleştirilir; kollar arası karşılaştırmada ise eşleştirme (sekans, pencere, tohum)
üçlüsü üzerindendir — aynı türbülans gerçeklemesi, farklı ölçüm modeli.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd

from .runner import load_records
from .stats import paired_test, bootstrap_ci

METRICS = ["ate_se3", "ate_noalign", "ate_sim3", "rpe_trans", "rpe_rot_deg"]


# Bir pencerenin ATE'sinin anlamlı olması için asgari GT yol uzunluğu [m].
# GEREKÇE (turb.sanity I1 bulgusu): EuRoC MH01'in ilk ~25 saniyesinde platform
# neredeyse hareketsizdir; bu pencerede GT yolu 0.08 m'dir ve ATE ~ 4e-4 m çıkar.
# Böyle pencereler hiçbir kolu ayırt etmez, yalnız etki büyüklüklerini seyreltir ve
# NEES'i dejenere kılar. Zhang & Scaramuzza 2018 (IROS) yörünge değerlendirme
# protokolü de bölüt uzunluğuna alt sınır koyar.
MIN_PATH_M = 1.0


def _prep(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df["status"] == "ok"].copy()
    for c in ("m_" + m for m in METRICS):
        if c not in df:
            df[c] = np.nan
    if "m_path_length_gt" in df:
        n0 = len(df)
        df = df[df["m_path_length_gt"].fillna(0.0) >= MIN_PATH_M].copy()
        if len(df) < n0:
            print(f"[analiz] dejenere pencere filtresi: {n0 - len(df)}/{n0} kayıt "
                  f"düştü (GT yolu < {MIN_PATH_M} m)")
    df["seq"] = df["c_sequence"]
    df["dr"] = df["c_D_over_r0"].astype(float)
    df["tc"] = df["c_theta_c_over_fov"].astype(float)
    df["t0"] = df["c_tau0_over_dt"].astype(float)
    df["seed"] = df["c_seed"].astype(int)
    df["tag"] = df["c_tag"].fillna("")
    return df


# =============================================================================
def analyze_sweep(df: pd.DataFrame, axis: str, outdir: str, metric="m_ate_se3",
                  tag_filter: str = ""):
    """Bozulma süpürmesi: eksen boyunca metrik + temizle eşleştirilmiş test."""
    d = df[(df["arm"] == "A0") & (df["tag"] == tag_filter)]
    clean = d[d["dr"] == 0.0]
    if clean.empty:
        # Etiketli kollar (ör. "tiltonly", "fe_superpoint_lg") kendi temiz hücrelerini
        # koşmaz: D/r0 = 0'da hiçbir kanal enjekte edilmediği için temiz koşum
        # ETİKETTEN BAĞIMSIZDIR ve etiketsiz kümeden ödünç alınır. Bu ödünç alma
        # olmadan `clean_mean`/`ratio_to_clean` boş kalır ve eşleştirilmiş test
        # tamamen NaN üzerinde koşar (turb.sanity çalışırken fark edildi).
        clean = df[(df["arm"] == "A0") & (df["tag"] == "") & (df["dr"] == 0.0)]
        if len(clean):
            print(f"[analiz] '{tag_filter}' için temiz referans etiketsiz kümeden "
                  f"ödünç alındı ({len(clean)} hücre)")
    ref = {(r["seq"], r["window"]): r[metric] for _, r in clean.iterrows()}

    rows, tests = [], {}
    for val, g in d[d["dr"] > 0].groupby(axis):
        vals = g[metric].values.astype(float)
        mean, lo, hi = bootstrap_ci(vals)
        base = np.array([ref.get((r["seq"], r["window"]), np.nan)
                         for _, r in g.iterrows()], dtype=float)
        pr = paired_test(base, vals)
        rows.append({axis: val, "n": len(vals), "mean": mean, "ci_low": lo, "ci_high": hi,
                     "median": float(np.nanmedian(vals)),
                     "clean_mean": float(np.nanmean(base)),
                     "ratio_to_clean": float(np.nanmean(vals) / np.nanmean(base))
                     if np.nanmean(base) else np.nan,
                     "diverged_frac": float(np.mean(g["m_diverged"].astype(bool)))
                     if "m_diverged" in g else np.nan,
                     "accepted_degradation": pr.accepted or (pr.direction == "kötüleşme"
                                                             and pr.wilcoxon_ok)})
        tests[str(val)] = pr.to_dict()
    if not rows:
        return pd.DataFrame(columns=[axis, "n", "mean", "ci_low", "ci_high"]), tests
    return pd.DataFrame(rows).sort_values(axis), tests


ARM_ORDER = ("A0", "A0C", "A1", "A2", "A3", "A4", "A5")
ARM_PAIRS = (("A0", "A1"), ("A0", "A2"), ("A0", "A3"), ("A0", "A4"), ("A0", "A5"),
             ("A0", "A0C"), ("A0C", "A3"), ("A0C", "A5"),
             ("A1", "A3"), ("A1", "A5"), ("A3", "A5"))


def analyze_arms(df: pd.DataFrame, outdir: str, metric="m_ate_se3",
                 raw: pd.DataFrame | None = None,
                 inter_arms: tuple | None = None):
    """E4 — kollar arası eşleştirilmiş karşılaştırma (rapor §6 üçlü ölçütü).

    İKİ TABAN, AÇIKÇA AYRILMIŞ (turb.sanity I2/I3 bulgusundan sonra eklendi):

      * `kesisim`  — BÜTÜN kolların başarılı olduğu hücreler. Birincil sonuç budur.
        Gerekçe: bir kol zor hücrelerde ıraksayıp analizden düşerse, kalan
        hücreler onun için kolaydır; ikili dropna ile hesaplanan ortalama bu
        yüzden YANLIDIR (MNAR — Little & Rubin 2019). Ölçüldü: A4 %67, A5 %17
        ıraksama oranıyla düşüyordu.
      * `ikili`    — yalnız o ikilinin ortak hücreleri. Yardımcı/şeffaflık için.

    Ayrıca her kol için IRAKSAMA ORANI birinci sınıf sonuç olarak döner: bir kol
    "daha düşük ATE" gösterip daha sık ıraksıyorsa bu bir kazanç değildir.
    """
    d = df[df["dr"] > 0]
    out, tests = [], {}

    # --- kol başına başarısızlık/ıraksama oranı (ham kayıtlardan) -------------
    fail_rate = {}
    if raw is not None and len(raw):
        rr = raw[raw["c_D_over_r0"].astype(float) > 0] if "c_D_over_r0" in raw else raw
        for arm, g in rr.groupby("arm"):
            nfail = int((g["status"] != "ok").sum())
            ndiv = int(g.get("m_diverged", pd.Series(dtype=bool)).fillna(False).sum())
            fail_rate[arm] = {"n": int(len(g)), "n_fail": nfail,
                              "fail_frac": nfail / max(len(g), 1),
                              "n_diverged": ndiv,
                              "div_frac": ndiv / max(len(g), 1)}
    tests["_kol_basarisizlik"] = fail_rate

    for dr, g in d.groupby("dr"):
        piv = g.pivot_table(index=["seq", "window", "seed"], columns="arm",
                            values=metric, aggfunc="first")
        arms = [a for a in ARM_ORDER if a in piv.columns]
        # Kesişim, KARŞILAŞTIRILAN kollar üzerinden alınır. Analiz dışı bırakılmış
        # bir kol (ör. E4'te bozuk A4 sütunu) kesişime dahil edilirse kendi yüksek
        # ıraksama oranıyla BÜTÜN kolların örneklemini küçültür ve gücü düşürür.
        ia = [a for a in arms if inter_arms is None or a in inter_arms]
        inter = piv[ia].dropna() if ia else piv.iloc[:0]
        tests[f"dr{dr:g}/_kesisim_n"] = int(len(inter))
        tests[f"dr{dr:g}/_kesisim_kollari"] = ia
        tests[f"dr{dr:g}/_kol_n"] = {a: int(piv[a].notna().sum()) for a in arms}

        for a in arms:
            v_all = piv[a].values.astype(float)
            mean, lo, hi = bootstrap_ci(v_all)
            row = {"D_over_r0": dr, "arm": a, "n": int(piv[a].notna().sum()),
                   "mean": mean, "ci_low": lo, "ci_high": hi,
                   "median": float(np.nanmedian(v_all)),
                   "n_kesisim": int(len(inter))}
            if len(inter) >= 3 and a in inter.columns:
                vi = inter[a].values.astype(float)
                mi, li, hi2 = bootstrap_ci(vi)
                row.update({"mean_kesisim": mi, "ci_low_kesisim": li,
                            "ci_high_kesisim": hi2,
                            "median_kesisim": float(np.nanmedian(vi))})
            row.update({"fail_frac": fail_rate.get(a, {}).get("fail_frac", np.nan),
                        "div_frac": fail_rate.get(a, {}).get("div_frac", np.nan)})
            out.append(row)

        for a, b in ARM_PAIRS:
            if a not in piv.columns or b not in piv.columns:
                continue
            # BİRİNCİL: bütün kolların kesişimi
            if len(inter) >= 3 and a in inter.columns and b in inter.columns:
                pr = paired_test(inter[a].values.astype(float),
                                 inter[b].values.astype(float))
                tests[f"dr{dr:g}/{b}_vs_{a}"] = dict(pr.to_dict(), taban="kesisim",
                                                     n=int(len(inter)))
            # YARDIMCI: yalnız ikilinin ortak hücreleri
            m = piv[[a, b]].dropna()
            if len(m) >= 3:
                pr2 = paired_test(m[a].values.astype(float), m[b].values.astype(float))
                tests[f"dr{dr:g}/{b}_vs_{a}__ikili"] = dict(pr2.to_dict(),
                                                            taban="ikili", n=int(len(m)))
    return pd.DataFrame(out), tests


def analyze_restoration(df: pd.DataFrame, outdir: str, metric="m_ate_se3"):
    """E3 / C2 — ΔPSNR ile ΔATE arasındaki ilişki.

    Eşleştirme birimi (sekans, pencere, tohum, D/r0): aynı türbülans gerçeklemesinde
    restorasyonlu ve restorasyonsuz koşum karşılaştırılır. Delta'lar referans kol
    "none"a (bozuk -> VIO) GÖRE tanımlıdır:
        ΔPSNR = PSNR(restorasyon) - PSNR(bozuk)      [algısal eksen]
        ΔATE  = ATE(restorasyon)  - ATE(bozuk)       [navigasyon ekseni]
    C2'nin iddiası, bu iki eksenin BEKLENDİĞİ GİBİ ilişkilenmediğidir.
    """
    d = df[df["dr"] > 0].copy()
    d["method"] = d["c_restoration"].fillna("none")
    d["psnr"] = d["quality"].apply(
        lambda q: q.get("psnr") if isinstance(q, dict) else np.nan)
    d["ssim"] = d["quality"].apply(
        lambda q: q.get("ssim") if isinstance(q, dict) else np.nan)
    idx = ["seq", "window", "seed", "dr"]
    ref = d[d["method"] == "none"].set_index(idx)
    rows = []
    for meth, g in d[d["method"] != "none"].groupby("method"):
        g = g.set_index(idx)
        common = g.index.intersection(ref.index)
        if len(common) < 3:
            continue
        d_psnr = (g.loc[common, "psnr"] - ref.loc[common, "psnr"]).astype(float)
        d_ate = (g.loc[common, metric] - ref.loc[common, metric]).astype(float)
        m = np.isfinite(d_psnr) & np.isfinite(d_ate)
        if m.sum() < 3:
            continue
        pr = paired_test(ref.loc[common, metric].values.astype(float),
                         g.loc[common, metric].values.astype(float))
        corr = float(np.corrcoef(d_psnr[m], d_ate[m])[0, 1])
        from scipy import stats as _st
        rho, pval = _st.spearmanr(d_psnr[m], d_ate[m])
        rows.append({"method": meth, "n": int(m.sum()),
                     "dPSNR_mean": float(d_psnr[m].mean()),
                     "dATE_mean": float(d_ate[m].mean()),
                     "dATE_median": float(np.median(d_ate[m])),
                     "pearson_dPSNR_dATE": corr,
                     "spearman": float(rho), "spearman_p": float(pval),
                     "ATE_improved_frac": float(np.mean(d_ate[m] < 0)),
                     "PSNR_improved_frac": float(np.mean(d_psnr[m] > 0)),
                     "accepted_ATE_gain": pr.accepted})
    tab = pd.DataFrame(rows)
    if len(tab):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6, 4.5))
            for meth, g in d[d["method"] != "none"].groupby("method"):
                gg = g.set_index(idx)
                common = gg.index.intersection(ref.index)
                if len(common) < 3:
                    continue
                ax.scatter((gg.loc[common, "psnr"] - ref.loc[common, "psnr"]),
                           (gg.loc[common, metric] - ref.loc[common, metric]),
                           s=12, alpha=.6, label=meth)
            ax.axhline(0, color="k", lw=.8); ax.axvline(0, color="k", lw=.8)
            ax.set_xlabel(r"$\Delta$PSNR [dB]  (restorasyon − bozuk)")
            ax.set_ylabel(r"$\Delta$ATE [m]  (restorasyon − bozuk)")
            ax.set_title("C2: algısal kalite ile navigasyon hatası ilişkisi")
            ax.legend(fontsize=8); ax.grid(alpha=.3)
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, "e3_psnr_vs_ate.png"), dpi=150)
            plt.close(fig)
        except Exception:
            pass
    return tab


from .metrics import NEES_BAND as _NB
_NEES_LO, _NEES_HI = _NB


def analyze_nees(df: pd.DataFrame, intersection_only: bool = True):
    """E5 — tutarlılık tablosu.

    `intersection_only` (varsayılan): NEES yalnız BÜTÜN kolların geçerli pencere
    ürettiği (sekans, pencere, tohum, D/r0) hücrelerinde hesaplanır. Gerekçe
    (turb.sanity I4): kol başına pencere sayısı 7 ile 66 arasında değişiyordu;
    Bar-Shalom §5.4 NEES testi ancak AYNI koşullardaki koşumlar üzerinde geçerlidir
    ve farklı pencere kümelerinden gelen NEES medyanları karşılaştırılamaz.
    Tablo, filtrelenmemiş sayıyı `n_ham` sütununda saklar (şeffaflık).
    """
    if "n_nees_pos_mean" not in df:
        return pd.DataFrame()
    d = df[df["n_nees_pos_mean"].notna()]
    n_ham = {k: len(g) for k, g in d.groupby(["dr", "arm"])}
    if intersection_only and len(d):
        key = ["seq", "window", "seed", "dr"]
        arms = sorted(d["arm"].unique())
        if len(arms) > 1:
            sets = [set(map(tuple, d[d["arm"] == a][key].values)) for a in arms]
            common = set.intersection(*sets)
            if len(common) >= 3:
                mask = [tuple(v) in common for v in d[key].values]
                d = d[np.array(mask)]
            else:
                print(f"[analiz] NEES kesişimi yetersiz ({len(common)}); "
                      f"filtrelenmemiş küme kullanılıyor")
    rows = []
    for (dr, arm), g in d.groupby(["dr", "arm"]):
        rows.append({
            "D_over_r0": dr, "arm": arm, "n": len(g),
            "nees_pos_mean": float(np.nanmean(g["n_nees_pos_mean"])),
            "nees_pos_median": float(np.nanmedian(g["n_nees_pos_median"]))
            if "n_nees_pos_median" in g else np.nan,
            "nees_rot_mean": float(np.nanmean(g["n_nees_rot_mean"]))
            if "n_nees_rot_mean" in g else np.nan,
            "coverage_pos": float(np.nanmean(g["n_coverage_pos"])),
            "coverage_rot": float(np.nanmean(g["n_coverage_rot"]))
            if "n_coverage_rot" in g else np.nan,
            "in_band_cells": float(np.mean(
                (g["n_nees_pos_mean"] >= _NEES_LO) & (g["n_nees_pos_mean"] <= _NEES_HI))),
        })
    return pd.DataFrame(rows).sort_values(["D_over_r0", "arm"])


def analyze_r0(df: pd.DataFrame):
    """E6 — DIMM r0 kestirimi doğruluğu."""
    rows = []
    for _, r in df.iterrows():
        st = r.get("arm_stats")
        if not isinstance(st, dict):
            continue
        if not np.isfinite(st.get("sigma_hat_median", np.nan)):
            continue
        rows.append({"seq": r["seq"], "dr": r["dr"], "seed": r["seed"], "arm": r["arm"],
                     "sigma_hat": st["sigma_hat_median"],
                     "sigma_oracle": st.get("sigma_px_oracle", np.nan),
                     "r0_hat_cm": 100 * st.get("r0_hat_median", np.nan),
                     "theta_c_hat_px": st.get("theta_c_hat_median", np.nan),
                     "theta_c_oracle_px": st.get("theta_c_px_oracle", np.nan)})
    d = pd.DataFrame(rows)
    if d.empty:
        return d, {}
    d = d.drop_duplicates(subset=["seq", "dr", "seed", "arm"])
    stat = {}
    m = d[np.isfinite(d["sigma_hat"]) & np.isfinite(d["sigma_oracle"])]
    if len(m) >= 3:
        stat["corr_sigma"] = float(np.corrcoef(m["sigma_hat"], m["sigma_oracle"])[0, 1])
        stat["bias_sigma"] = float(np.mean(m["sigma_hat"] - m["sigma_oracle"]))
        stat["rel_err_sigma"] = float(np.mean(np.abs(m["sigma_hat"] - m["sigma_oracle"])
                                              / m["sigma_oracle"]))
        r0_true = 100 * 0.10 / m["dr"].replace(0, np.nan)
        stat["corr_r0"] = float(np.corrcoef(m["r0_hat_cm"], r0_true)[0, 1])
    return d, stat


# =============================================================================
def _fig_sweep(tab, axis, outdir, xlabel, title, fname):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    x = tab[axis].values
    ax.errorbar(x, tab["mean"], yerr=[tab["mean"] - tab["ci_low"],
                                      tab["ci_high"] - tab["mean"]],
                fmt="o-", capsize=3, label="ATE (SE3) ortalama ± %95 GA")
    if "clean_mean" in tab and np.isfinite(tab["clean_mean"]).any():
        ax.axhline(float(tab["clean_mean"].iloc[0]), color="k", ls="--", lw=1,
                   label="temiz (D/r₀ = 0)")
    ax.set_xlabel(xlabel); ax.set_ylabel("ATE [m]"); ax.set_title(title)
    ax.grid(alpha=.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, fname), dpi=150)
    plt.close(fig)


def _fig_arms(tab, outdir, fname="e4_arms.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    for arm, g in tab.groupby("arm"):
        ax.errorbar(g["D_over_r0"], g["mean"],
                    yerr=[g["mean"] - g["ci_low"], g["ci_high"] - g["mean"]],
                    fmt="o-", capsize=3, label=arm)
    ax.set_xlabel(r"$D/r_0$"); ax.set_ylabel("ATE [m]")
    ax.set_title("E4 — ölçüm modeli kolları")
    ax.grid(alpha=.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, fname), dpi=150)
    plt.close(fig)


# =============================================================================
def main():
    ap = argparse.ArgumentParser(description="Kanıt analizi")
    ap.add_argument("--exp", required=True)
    ap.add_argument("--records", default=None)
    ap.add_argument("--metric", default="m_ate_se3")
    a = ap.parse_args()

    rec = a.records or f"evidence/{a.exp}/records.jsonl"
    outdir = os.path.dirname(rec)
    raw = load_records(rec)
    df = _prep(raw)
    print(f"[{a.exp}] {len(df)} geçerli kayıt ({len(raw)} ham)")
    report = [f"# {a.exp} — kanıt özeti", "",
              f"Geçerli kayıt: {len(df)}", ""]

    if a.exp.startswith("E1a") or a.exp in ("TARTAN",):
        tab, tests = analyze_sweep(df, "dr", outdir, a.metric)
        tab.to_csv(os.path.join(outdir, "sweep_dr.csv"), index=False)
        _fig_sweep(tab, "dr", outdir, r"$D/r_0$", f"{a.exp}: bozulma şiddeti",
                   "sweep_dr.png")
        report += ["## Şiddet ekseni (D/r₀)", tab.to_markdown(index=False), ""]
        tt = df[df["tag"] == "tiltonly"]
        if len(tt):
            tab2, tests2 = analyze_sweep(df, "dr", outdir, a.metric,
                                         tag_filter="tiltonly")
            tab2.to_csv(os.path.join(outdir, "sweep_dr_tiltonly.csv"), index=False)
            report += ["## Yalnız-tilt kanalı (bulanıklık kapalı)",
                       tab2.to_markdown(index=False), ""]
            tests["tiltonly"] = tests2
        with open(os.path.join(outdir, "tests.json"), "w") as f:
            json.dump(tests, f, indent=1)

    elif a.exp == "E1b":
        tab, tests = analyze_sweep(df, "t0", outdir, a.metric)
        tab.to_csv(os.path.join(outdir, "sweep_tau0.csv"), index=False)
        _fig_sweep(tab, "t0", outdir, r"$\tau_0/\Delta t$", "E1-b: zamansal rejim",
                   "sweep_tau0.png")
        report += ["## Zamansal rejim (τ₀/Δt)", tab.to_markdown(index=False), ""]
        with open(os.path.join(outdir, "tests.json"), "w") as f:
            json.dump(tests, f, indent=1)

    elif a.exp == "E1c":
        tab, tests = analyze_sweep(df, "tc", outdir, a.metric)
        tab.to_csv(os.path.join(outdir, "sweep_thetac.csv"), index=False)
        _fig_sweep(tab, "tc", outdir, r"$\theta_c/\mathrm{FOV}$",
                   "E1-c: uzaysal korelasyon ölçeği", "sweep_thetac.png")
        report += ["## Uzaysal korelasyon (θ_c/FOV)", tab.to_markdown(index=False), ""]
        with open(os.path.join(outdir, "tests.json"), "w") as f:
            json.dump(tests, f, indent=1)

    elif a.exp == "E8":
        report += ["## E8-a — genellenebilirlik (ön-uç yapılandırmaları)", ""]
        all_tests = {}
        for tag in sorted(t for t in df["tag"].unique() if str(t).startswith("fe_")):
            tab, tests = analyze_sweep(df, "dr", outdir, a.metric, tag_filter=tag)
            tab.to_csv(os.path.join(outdir, f"sweep_{tag}.csv"), index=False)
            report += [f"### {tag}", tab.to_markdown(index=False) if len(tab)
                       else "(veri yok)", ""]
            all_tests[tag] = tests
        with open(os.path.join(outdir, "tests.json"), "w") as f:
            json.dump(all_tests, f, indent=1)

    elif a.exp == "E3":
        tab = analyze_restoration(df, outdir, a.metric)
        tab.to_csv(os.path.join(outdir, "restoration.csv"), index=False)
        report += ["## E3 — restorasyon kolları (C2)", "",
                   tab.to_markdown(index=False) if len(tab) else "(veri yok)", "",
                   "ΔPSNR > 0 ve ΔATE < 0 aynı anda gerçekleşiyorsa restorasyon hem "
                   "algısal hem geometrik olarak yardımcı olmuştur; işaretler ayrışıyorsa "
                   "C2 hipotezi (PSNR ≠ ATE) desteklenir.", ""]

    elif a.exp.startswith("E4") or a.exp.startswith("TARTANarms"):
        # E4'ün A4 sütunu bozuk/karışıktır (eski seçim mantığı) ve AYRI bir deneyde
        # (E4extra) yeniden koşulmaktadır; kesişime dahil edilmez.
        ia = tuple(x for x in ARM_ORDER if not (a.exp == "E4" and x == "A4"))
        tab, tests = analyze_arms(df, outdir, a.metric, raw=raw, inter_arms=ia)
        tab.to_csv(os.path.join(outdir, "arms.csv"), index=False)
        _fig_arms(tab, outdir)
        report += ["## Kollar (ATE)", tab.to_markdown(index=False), ""]
        nees = analyze_nees(df)
        if len(nees):
            nees.to_csv(os.path.join(outdir, "nees.csv"), index=False)
            report += ["## E5 tutarlılık (pencere-göreli NEES)",
                       nees.to_markdown(index=False), ""]

        # ---- C3'ün ANA İDDİASI: NEES üzerinde EŞLEŞTİRİLMİŞ test ----------
        # Yukarıdaki tablo kol başına toplulaştırır; hakem eşleştirilmiş bir test
        # ister ve haklıdır (hücreden hücreye NEES saçılımı, kollar arası farktan
        # büyüktür). İki ölçüt raporlanır:
        #   log10(NEES)  — fark, doğrudan KAT cinsinden okunur (10^fark)
        #   bant uzaklığı — bandın İÇİNDE 0, dışında en yakın kenara log uzaklık.
        #     Bu ikincisi gereklidir çünkü "daha düşük NEES" ancak bandın ÜSTÜNDE
        #     iken iyidir; aşırı şişirme bandın ALTINA düşürür ve o bir kazanç
        #     değildir. (A0C k=16 ile bandın içine giriyor; daha derini doğruluğu
        #     harcar — bkz. config.yaml cal_factor gerekçesi.)
        if "n_nees_pos_median" in df.columns:
            _lo_b, _hi_b = _NEES_LO, _NEES_HI
            dnl = df.copy()
            _v = pd.to_numeric(dnl["n_nees_pos_median"], errors="coerce")
            dnl["n_nees_log"] = np.log10(_v.where(_v > 0))
            dnl["n_nees_bantd"] = np.maximum.reduce([
                np.zeros(len(dnl)),
                np.log10(_v.where(_v > 0) / _hi_b).fillna(0).values,
                np.log10(_lo_b / _v.where(_v > 0)).fillna(0).values])
            dnl.loc[~np.isfinite(dnl["n_nees_log"]), "n_nees_bantd"] = np.nan
            for _met, _bas, _birim in (("n_nees_log", "log10(NEES)", "kat"),
                                       ("n_nees_bantd", "bant uzaklığı (log10)", "kat")):
                try:
                    _t, _ts = analyze_arms(dnl, outdir, _met, raw=raw, inter_arms=ia)
                except Exception as _e:
                    report += [f"*{_bas} eşleştirilmiş testi koşulamadı: "
                               f"{type(_e).__name__}: {_e}*", ""]
                    continue
                _rows = [(k, v) for k, v in _ts.items()
                         if isinstance(v, dict) and "mean_diff" in v
                         and v.get("taban") == "kesisim"]
                if not _rows:
                    continue
                report += [f"## C3 eşleştirilmiş tutarlılık testi — {_bas}", "",
                           "`fark` negatifse ikinci kol daha TUTARLIdır. "
                           "`kat` = 10^|fark| (etki büyüklüğü).", "",
                           "| karşılaştırma | n | fark | %95 GA | kat | Wilcoxon p |",
                           "|---|---|---|---|---|---|"]
                for k, v in sorted(_rows):
                    report.append(
                        f"| {k} | {v['n']} | {v['mean_diff']:+.3f} | "
                        f"[{v['ci_low']:+.3f}, {v['ci_high']:+.3f}] | "
                        f"{10 ** abs(v['mean_diff']):.2f}x | {v['wilcoxon_p']:.2e} |")
                report.append("")
                with open(os.path.join(outdir, f"tests_{_met}.json"), "w") as _f:
                    json.dump(_ts, _f, indent=1, default=str)
        r0d, r0s = analyze_r0(df)
        if len(r0d):
            r0d.to_csv(os.path.join(outdir, "r0_estimates.csv"), index=False)
            report += ["## E6 çevrimiçi r₀ kestirimi",
                       json.dumps(r0s, indent=1), ""]
        fr = tests.get("_kol_basarisizlik") or {}
        if fr:
            report += ["## Kol dayanıklılığı (BİRİNCİL SONUÇ — ATE ile birlikte okunmalı)",
                       "", "Bir kolun düşük ATE'si, yalnız kolay hücrelerde hayatta "
                       "kaldığı için düşük olabilir. Aşağıdaki oranlar ıraksayıp "
                       "analizden düşen hücreleri sayar.", "",
                       "| kol | n | başarısız | oran | ıraksak | oran |",
                       "|---|---|---|---|---|---|"]
            for k in sorted(fr):
                v = fr[k]
                report.append(f"| {k} | {v['n']} | {v['n_fail']} | "
                              f"%{100*v['fail_frac']:.0f} | {v['n_diverged']} | "
                              f"%{100*v['div_frac']:.0f} |")
            report.append("")
        report += ["## Eşleştirilmiş testler (üçlü ölçüt)", "",
                   "Taban `kesisim` = bütün kolların hayatta kaldığı hücreler "
                   "(BİRİNCİL); `ikili` = yalnız o ikilinin ortak hücreleri (yardımcı).",
                   "",
                   "| karşılaştırma | taban | n | ort. fark | %95 GA | Wilcoxon p | işaret | KABUL |",
                   "|---|---|---|---|---|---|---|---|"]
        for k, v in tests.items():
            if k.startswith("_") or "/_" in k or not isinstance(v, dict) \
                    or "mean_diff" not in v:
                continue
            report.append(f"| {k} | {v.get('taban','-')} | {v['n']} | {v['mean_diff']:+.4f} | "
                          f"[{v['ci_low']:+.4f}, {v['ci_high']:+.4f}] | "
                          f"{v['wilcoxon_p']:.2e} | {v['n_improved']}/{v['n_worsened']} | "
                          f"{'EVET' if v['accepted'] else 'hayır'} |")
        with open(os.path.join(outdir, "tests.json"), "w") as f:
            json.dump(tests, f, indent=1)

    with open(os.path.join(outdir, "report.md"), "w") as f:
        f.write("\n".join(report) + "\n")
    print("\n".join(report[:40]))
    print(f"\n✓ {outdir}/report.md")


if __name__ == "__main__":
    main()
