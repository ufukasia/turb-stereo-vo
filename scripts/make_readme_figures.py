#!/usr/bin/env python
"""Draw every README figure from the result files in results/ (nothing is typed in by hand).

    python scripts/make_readme_figures.py            # writes docs/figures/*.png (each 2200 px wide)

hero.png         frames and causal trajectories of the window with the largest gain (results/HERO)
degradation.png  ATE growth and failure rate of four stereo odometry systems
gate.png         paired ATE and ANEES ratios of the gated arms in the confirmatory grids
real.png         online tilt estimate against the scintillometer on real imagery
"""
import json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, ROOT)
import figstyle as fs   # noqa: E402

OUT = os.path.join(ROOT, "docs", "figures")
BLUE, SKY, VERM, ORANGE, GREEN, PINK, GREY = ("#0072B2", "#56B4E9", "#D55E00", "#E69F00", "#009E73", "#CC79A7",
                                              "#6E6E6E")
fs.SYS = {
    "gt":       (fs.GT, "-", 2.8, "Ground truth"),
    "gate":     (BLUE, "-", 2.0, "Adaptive gate (proposed)"),
    "oracle":   (SKY, (0, (4, 2)), 1.5, "Oracle gate"),
    "baseline": (VERM, (0, (1.5, 1.5)), 1.5, "Baseline (fixed noise and gates)"),
    "pyeco":    (VERM, (0, (1.5, 1.5)), 1.5, "PyEco-SLAM (baseline)"),
    "orb":      (GREEN, (0, (6, 2, 1.5, 2)), 1.5, "ORB-SLAM3"),
    "airslam":  (ORANGE, (0, (5, 2)), 1.5, "AirSLAM"),
    "droid":    (GREY, (0, (8, 3)), 1.5, "DROID-SLAM"),
    "generic":  (PINK, (0, (6, 2, 1, 2, 1, 2)), 1.5, "Gated generic"),
    "scalar":   (fs.NEUTRAL_BAR, "-", 1.5, "Scalar variant (fixed gates)"),
}
COND = ["R0-8", "R5-8", "R5-4", "R0-16", "R0-0"]
COND_LABEL = ["Quasi-static\nD/r0 = 8", "Fast\nD/r0 = 8", "Fast\nD/r0 = 4", "Quasi-static\nD/r0 = 16", "Clean\n"]


def J(*p):
    with open(os.path.join(ROOT, *p)) as f:
        return json.load(f)


# ------------------------------------------------------------------------------------------------ hero
def hero_window(arm):
    """Aligned causal trajectory of one arm in the selected window; asserts the stored window ATE."""
    from turb.metrics import load_gt, load_traj, associate, umeyama
    H = os.path.join(ROOT, "results", "HERO")
    meta = J("results", "HERO", "window.json")
    et, exyz, eq = load_traj(os.path.join(H, f"traj_{arm}_realtime.csv"))
    m = (et >= meta["t_start_ns"]) & (et <= meta["t_end_ns"])
    t, e, _, g, _ = associate(et[m], exyz[m], eq[m], *load_gt(os.path.join(H, "gt_window.csv")))
    R, tr, _ = umeyama(e, g, with_scale=False)
    al = (R @ e.T).T + tr
    err = np.linalg.norm(g - al, axis=1)
    ate = float(np.sqrt(np.mean(err ** 2)))
    assert abs(ate - meta["ate_m"][arm]) < 1e-9, (arm, ate, meta["ate_m"][arm])
    return (t - meta["t_start_ns"]) * 1e-9, al, g, err, ate


def hero(plt):
    meta = J("results", "HERO", "window.json")
    A0, KA = hero_window("A0"), hero_window("KAPI")
    fig = fs.new_fig(plt, fs.HEIGHTS["hero"])
    gs = fig.add_gridspec(2, 2, left=0.065, right=0.985, top=0.905, bottom=0.075, hspace=0.30, wspace=0.16,
                          height_ratios=[1.0, 0.92])
    for j, (name, title) in enumerate((("frame_clean.png", "(a) EuRoC MH05, clean frame"),
                                       ("frame_turbulent.png", "(b) Same frame, D/r0 = 8 (fast regime)"))):
        ax = fig.add_subplot(gs[0, j])
        ax.imshow(plt.imread(os.path.join(ROOT, "results", "HERO", name)), cmap="gray", vmin=0, vmax=1,
                  interpolation="lanczos")
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(title)
    ax = fig.add_subplot(gs[1, 0])
    ex = fs.topview(ax, A0[2])
    fs.plot(ax, "gt", A0[2][:, ex[0]], A0[2][:, ex[1]])
    fs.plot(ax, "baseline", A0[1][:, ex[0]], A0[1][:, ex[1]])
    fs.plot(ax, "gate", KA[1][:, ex[0]], KA[1][:, ex[1]])
    ax.set_title("(c) Causal trajectories in the window")
    ax = fig.add_subplot(gs[1, 1])
    fs.plot(ax, "baseline", A0[0], A0[3])
    fs.plot(ax, "gate", KA[0], KA[3])
    ax.axvline(meta["t_frame_s"], color="#888888", lw=0.8, ls=":")
    ax.set_xlim(0, A0[0][-1]); ax.set_ylim(0, None)
    ax.set_xlabel("time in the window [s]"); ax.set_ylabel("position error [m]")
    ax.set_title("(d) Position error after rigid alignment")
    fs.legend(fig, [fs.handle("gt"), fs.handle("gate"), fs.handle("baseline")], ncol=3)
    fs.save(fig, OUT, "hero.png")


# ------------------------------------------------------------------------------------------------ degradation
def degradation(plt):
    F = J("results", "FAZ2", "analiz.json")["H1"]["kirilma_tablo"]
    R = J("results", "FAZ2_rakip_ek", "rakip_ek_analiz.json")["HR1"]["rakipler"]
    fig = fs.new_fig(plt, fs.HEIGHTS["row"])
    gs = fig.add_gridspec(1, 2, left=0.075, right=0.985, top=0.80, bottom=0.17, wspace=0.22)
    a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    lv = [r["D_over_r0"] for r in F if r["D_over_r0"] > 0]
    fs.plot(a, "pyeco", lv, [r["oran"] for r in F if r["D_over_r0"] > 0], marker="o", ms=6, mec="white", mew=0.6)
    fs.plot(b, "pyeco", lv, [100 * r["irak_oran"] for r in F if r["D_over_r0"] > 0], marker="o", ms=6,
            mec="white", mew=0.6)
    for key, name, mk in (("orb", "orbslam3", "s"), ("airslam", "airslam", "^"), ("droid", "droid", "D")):
        d = R[name]
        ls = [l for l in ("4", "8", "16", "24") if d[l].get("medyan_oran")]
        fs.plot(a, key, [float(l) for l in ls], [d[l]["medyan_oran"] for l in ls], marker=mk, ms=6, mec="white",
                mew=0.6)
        ls = ["4", "8", "16", "24"]
        fs.plot(b, key, [float(l) for l in ls], [100 * d[l]["n_basarisiz"] / d[l]["n_pencere"] for l in ls],
                marker=mk, ms=6, mec="white", mew=0.6)
    a.axhline(2, color="#888888", lw=0.8, ls=":")
    a.set_yscale("log"); a.set_yticks([1, 2, 5, 10, 20, 50]); a.set_yticklabels(["1", "2", "5", "10", "20", "50"])
    a.set_ylim(1, 70)
    for ax in (a, b):
        ax.set_xticks([4, 8, 12, 16, 20, 24]); ax.set_xlim(2.5, 25.5); ax.set_xlabel("turbulence strength D/r0")
    a.set_ylabel("ATE / turbulence-free ATE"); b.set_ylabel("failed windows [%]"); b.set_ylim(-3, 103)
    a.set_title("(a) ATE growth"); b.set_title("(b) Failure rate")
    fs.legend(fig, [fs.handle("pyeco", marker="o"), fs.handle("orb", marker="s"), fs.handle("airslam", marker="^"),
                    fs.handle("droid", marker="D")], ncol=4)
    fs.save(fig, OUT, "degradation.png")


# ------------------------------------------------------------------------------------------------ gate
def _pts(ax, key, x, rows, marker):
    c = fs.SYS[key][0]
    med = np.array([r["medyan_oran"] for r in rows])
    lo = np.array([r["ga95"][0] for r in rows]); hi = np.array([r["ga95"][1] for r in rows])
    ax.errorbar(x, med, yerr=[med - lo, hi - med], fmt=marker, color=c, ms=7, mec="white", mew=0.6, lw=1.4,
                capsize=2.5, zorder=4 if key == "gate" else 3)


def gate(plt):
    K = J("results", "KAPI", "dogrula_analiz.json")
    s7 = os.path.join(ROOT, "results", "SAPMA7", "s7_analiz.json")
    G = json.load(open(s7)) if os.path.exists(s7) else None
    fam = {"R0-8": "B", "R5-8": "B", "R5-4": "B", "R0-16": "I1", "R0-0": "T"}
    fig = fs.new_fig(plt, fs.HEIGHTS["row"])
    gs = fig.add_gridspec(1, 2, left=0.075, right=0.985, top=0.80, bottom=0.21, wspace=0.22, width_ratios=[1.25, 1])
    a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    x5, x4 = np.arange(5.0), np.arange(4.0)
    _pts(a, "scalar", x4 - 0.27, [K["I4"][f"A3S/A0@{c}"] for c in COND[:4]], "s")
    _pts(a, "scalar", [4 - 0.27], [K["T2_ATE"]["A3S/A0@R0-0"]], "s")
    _pts(a, "gate", x5 - 0.09, [K[fam[c]][f"KAPI/A0@{c}"] for c in COND], "o")
    hs = [fs.handle("gate", marker="o")]
    if G is not None:
        _pts(a, "generic", x4 + 0.09, [G["HG1"][f"AGg/A0@{c}"] for c in COND[:4]], "^")
        _pts(a, "generic", [4 + 0.09], [G["HG3"]["AGg/A0@R0-0"]], "^")
        hs.append(fs.handle("generic", marker="^"))
    _pts(a, "oracle", x4 + 0.27, [K["I4"][f"A1g/A0@{c}"] for c in COND[:4]], "D")
    hs += [fs.handle("oracle", marker="D"), fs.handle("scalar", marker="s")]
    a.axhline(1, color="#888888", lw=0.8, ls=":")
    a.set_xticks(x5); a.set_xticklabels(COND_LABEL); a.set_ylabel("ATE / baseline ATE"); a.set_ylim(0.2, 1.2)
    a.set_title("(a) Accuracy, paired ratio with 95% CI")
    _pts(b, "gate", x4 - 0.18, [K["I5_poz"][f"KAPI/A0@{c}"] for c in COND[:4]], "o")
    if G is not None:   # S3: position ANEES ratios of Deviation 7
        _pts(b, "generic", x4, [G["S3"][f"AGg/A0@{c}"] for c in COND[:4]], "^")
    _pts(b, "oracle", x4 + 0.18, [K["I5_poz"][f"A1g/A0@{c}"] for c in COND[:4]], "D")
    b.axhline(1, color="#888888", lw=0.8, ls=":")
    b.set_yscale("log"); b.set_yticks([0.02, 0.05, 0.1, 0.2, 0.5, 1]); b.set_ylim(0.015, 1.3)
    b.set_yticklabels(["0.02", "0.05", "0.1", "0.2", "0.5", "1"])
    b.set_xticks(x4); b.set_xticklabels(COND_LABEL[:4]); b.set_ylabel("ANEES / baseline ANEES")
    b.set_title("(b) Consistency")
    for ax in (a, b):
        ax.tick_params(axis="x", labelsize=9.5)
    for h in hs:
        h.set_linestyle("none")
    fs.legend(fig, hs, ncol=len(hs))
    fs.save(fig, OUT, "gate.png")


# ------------------------------------------------------------------------------------------------ real imagery
def real(plt):
    S = J("results", "GERCEK_VERI", "saha_ozet.json")["ort15"]["2021-07"]
    with open(os.path.join(ROOT, "results", "GERCEK_VERI", "saha_temmuz.jsonl")) as f:
        rows = [json.loads(l) for l in f]
    rows = [r for r in rows if r["lk"] == 15 and r["ref"] == "ort" and r["durum"] == "ok"]
    cn2, s2 = np.array([r["cn2"] for r in rows]), np.array([r["sigma_px"] ** 2 for r in rows])
    assert len(cn2) == S["n_ok"], (len(cn2), S["n_ok"])
    fig = fs.new_fig(plt, fs.HEIGHTS["row"])
    gs = fig.add_gridspec(1, 2, left=0.075, right=0.985, top=0.80, bottom=0.17, wspace=0.22)
    a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    a.scatter(cn2 * 1e15, s2, s=22, color=BLUE, alpha=0.65, edgecolor="white", linewidth=0.4, zorder=3)
    k, c = np.polyfit(np.log(cn2), np.log(s2), 1)
    xs = np.array([cn2.min(), cn2.max()])
    a.plot(xs * 1e15, np.exp(c) * xs ** k, color=fs.GT, lw=1.5, ls=(0, (5, 2)), zorder=4)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xticks([2, 5, 10, 20, 40]); a.set_xticklabels(["2", "5", "10", "20", "40"])
    a.set_yticks([0.2, 0.5, 1, 2]); a.set_yticklabels(["0.2", "0.5", "1", "2"])
    from matplotlib.ticker import NullFormatter
    a.xaxis.set_minor_formatter(NullFormatter()); a.yaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel("scintillometer Cn² [10⁻¹⁵ m⁻²ᐟ³]"); a.set_ylabel("estimated tilt variance [px²]")
    a.set_title(f"(a) {S['n_ok']} one-minute records, 724 m path")
    names = ["Rank\ncorrelation", "Log-log\nslope", "Estimated /\nmeasured Cn²"]
    val = [S["G1_spearman"], S["G2_egim"], S["G3_oran"]]
    ci = [S["G1_ga"], S["G2_ga"], S["G3_ga"]]
    ref = [None, 1.0, 0.375]
    x = np.arange(3.0)
    b.errorbar(x, val, yerr=[[v - c_[0] for v, c_ in zip(val, ci)], [c_[1] - v for v, c_ in zip(val, ci)]],
               fmt="o", color=BLUE, ms=7, mec="white", mew=0.6, lw=1.4, capsize=2.5, zorder=4)
    for xi, r in zip(x, ref):
        if r is not None:
            b.plot([xi - 0.3, xi + 0.3], [r, r], color=fs.GT, lw=1.5, ls=(0, (5, 2)), zorder=3)
    b.set_xticks(x); b.set_xticklabels(names); b.set_xlim(-0.6, 2.6); b.set_ylim(0, 1.15)
    b.set_ylabel("value with 95% CI"); b.set_title("(b) Pre-declared criteria")
    from matplotlib.lines import Line2D
    fs.legend(fig, [Line2D([], [], color=BLUE, marker="o", ls="none", ms=7, mec="white", mew=0.6,
                           label="Online estimate (unchanged code)"),
                    Line2D([], [], color=fs.GT, lw=1.5, ls=(0, (5, 2)), label="Fit in (a); theory in (b)")], ncol=2)
    assert abs(k - S["G2_egim"]) < 5e-3, (k, S["G2_egim"])
    fs.save(fig, OUT, "real.png")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    plt = fs.style()
    only = sys.argv[1:]
    for name, fn in (("hero", hero), ("degradation", degradation), ("gate", gate), ("real", real)):
        if not only or name in only:
            fn(plt)
    fs.check_all(OUT)
