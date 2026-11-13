"""Makale sekilleri. Butun sayilar evidence/ kayitlarindan okunur.

TASARIM KURALLARI (dataviz):
  * CIFT EKSEN YOK -- farkli olcekli buyuklukler ayri panellere konur
  * kategorik renkler SABIT sirada: blue, orange, aqua, yellow
  * ince isaretler, geri planda izgara, secici dogrudan etiket
  * >=2 seri oldugunda lejant her zaman var
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID = "#1a1a19", "#5c5b55", "#d8d7d0"
plt.rcParams.update({
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.edgecolor": INK2, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK,
    "grid.color": GRID, "grid.linewidth": 0.4,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "legend.frameon": False,
})
ROOT = "."
OUT = sys.argv[1] if len(sys.argv) > 1 else "figures"
os.makedirs(OUT, exist_ok=True)


def load(e):
    p = os.path.join(ROOT, f"evidence/{e}/records.jsonl")
    return [json.loads(l) for l in open(p) if l.strip()] if os.path.exists(p) else []


def tidy(ax):
    ax.grid(True, axis="y", zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ---------------------------------------------------------------- Sekil 1 + 2
by = defaultdict(list)
for r in load("E1a") + load("E0"):
    if r.get("status") != "ok" or not r.get("m_valid"):
        continue
    fs = r.get("frontend_stats") or {}
    by[r["cell"].get("D_over_r0") or 0].append(
        (r["m_ate_se3"], r.get("m_diverged", False),
         fs.get("stereo_run_ratio"), fs.get("avg_active_depth_tracks")))
dr = np.array(sorted(by))
ate = np.array([np.median([z[0] for z in by[d]]) for d in dr])
div = np.array([100 * np.mean([z[1] for z in by[d]]) for d in dr])
ster = np.array([100 * np.median([z[2] for z in by[d] if z[2] is not None]) for d in dr])
trk = np.array([np.median([z[3] for z in by[d] if z[3] is not None]) for d in dr])

fig, ax = plt.subplots(2, 1, figsize=(3.4, 3.1), sharex=True)
ax[0].plot(dr, ate, "-o", color=C[0], lw=1.4, ms=4, zorder=3)
ax[0].set_yscale("log"); ax[0].set_ylabel("ATE medyan [m]")
ax[0].axvline(5, color=INK2, lw=0.7, ls="--", zorder=1)
ax[0].annotate("mimari eşik\n$D/r_0\\approx5$", xy=(5, ate[dr == 5][0]),
               xytext=(5.8, 0.09), color=INK2, fontsize=6.2,
               arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
tidy(ax[0])
ax[1].plot(dr, div, "-o", color=C[1], lw=1.4, ms=4, zorder=3)
ax[1].axvline(5, color=INK2, lw=0.7, ls="--", zorder=1)
ax[1].set_ylabel("ıraksama [%]"); ax[1].set_xlabel("$D/r_0$")
tidy(ax[1])
fig.align_ylabels(ax)
fig.savefig(f"{OUT}/fig1_bozulma.pdf"); plt.close(fig)

fig, ax = plt.subplots(2, 1, figsize=(3.4, 3.1), sharex=True)
ax[0].plot(dr, ster, "-o", color=C[2], lw=1.4, ms=4, zorder=3)
ax[0].axvline(5, color=INK2, lw=0.7, ls="--", zorder=1)
ax[0].set_ylabel("stereo tetiklenen kare [%]"); tidy(ax[0])
ax[1].plot(dr, trk, "-o", color=C[3], lw=1.4, ms=4, zorder=3, label="derinlikli track")
ax[1].axhline(50, color="#e34948", lw=1.0, ls=":", zorder=2,
              label="ikmal eşiği (50)")
ax[1].axvline(5, color=INK2, lw=0.7, ls="--", zorder=1)
ax[1].set_ylabel("aktif derinlikli track"); ax[1].set_xlabel("$D/r_0$")
ax[1].legend(loc="upper right"); tidy(ax[1])
fig.align_ylabels(ax)
fig.savefig(f"{OUT}/fig2_mekanizma.pdf"); plt.close(fig)

# ---------------------------------------------------------------- Sekil 3 (C2)
meth = ["none", "temporal_median", "lucky", "elastic"]
lbl = ["yok", "zamansal\nmedyan", "lucky", "elastik"]
clean_psnr = [np.nan, 11.91, 9.19, 12.84]          # turb/_c2_probe.py, TEMIZ veri
bw = defaultdict(list)
for r in load("E3win"):
    if r.get("status") != "ok" or not r.get("m_valid"):
        continue
    c = r["cell"]
    if abs((c.get("D_over_r0") or 0) - 3.0) > 1e-9:
        continue
    bw[c.get("tag")].append(r["m_ate_se3"])
wins = ["win_none", "win3_temporal_median", "win5_temporal_median", "win7_temporal_median"]
wl = ["yok", "3 kare", "5 kare", "7 kare"]
wv = [np.median(bw[t]) if bw.get(t) else np.nan for t in wins]

fig, ax = plt.subplots(1, 2, figsize=(3.4, 1.75))
xs = np.arange(3)
ax[0].bar(xs, clean_psnr[1:], color=C[0], width=0.62, zorder=3)
ax[0].axhline(40, color=INK2, lw=0.6, ls="--", zorder=2)
ax[0].set_xticks(xs); ax[0].set_xticklabels(lbl[1:], fontsize=6)
ax[0].set_ylabel("PSNR [dB]"); ax[0].set_title("(a) türbülans YOK", fontsize=7)
for i, v in enumerate(clean_psnr[1:]):
    ax[0].text(i, v + 0.6, f"{v:.1f}", ha="center", fontsize=6, color=INK)
tidy(ax[0])
ax[1].bar(np.arange(4), wv, color=[GRID, C[1], C[1], C[1]], width=0.62, zorder=3)
ax[1].set_xticks(np.arange(4)); ax[1].set_xticklabels(wl, fontsize=6)
ax[1].set_ylabel("ATE medyan [m]")
ax[1].set_title("(b) pencere boyu, $D/r_0{=}3$", fontsize=7)
tidy(ax[1])
fig.tight_layout(pad=0.4)
fig.savefig(f"{OUT}/fig3_c2.pdf"); plt.close(fig)

# ---------------------------------------------------------------- Sekil 4 (kollar)
def paired(exp):
    per, ate_ = defaultdict(dict), defaultdict(dict)
    for r in load(exp):
        c = r["cell"]
        if r.get("status") != "ok" or (c.get("D_over_r0") or 0) <= 0:
            continue
        k = (c["sequence"], c["seed"], c["D_over_r0"], r.get("window"))
        if r.get("m_valid"):
            ate_[k][r["arm"]] = r["m_ate_se3"]
        if r.get("n_valid") and r.get("n_nees_pos_median"):
            per[k][r["arm"]] = r["n_nees_pos_median"]
    return per, ate_


def rat(d, a, b):
    ks = [k for k in d if a in d[k] and b in d[k] and d[k][b] > 1e-12]
    return float(np.median([d[k][a] / d[k][b] for k in ks])) if len(ks) >= 4 else np.nan


arms = ["A1", "A5", "A3", "A0C"]
# ETIKETLER KISA: uzun cok satirli etiketler 3.4 inc sutunda CAKISIYORDU.
# Kolun ne oldugu sekil altyazisinda tanimlanir.
alab = ["A1", "A5", "A3", "A0C"]
p4, a4 = paired("E4"); p9, a9 = paired("E9")
g4 = [rat(p4, "A0", x) for x in arms]
g9 = [rat(p9, "A0", x) for x in arms]
c4 = [rat(a4, x, "A0") for x in arms]

fig, ax = plt.subplots(1, 2, figsize=(3.45, 1.75))
x = np.arange(4); w = 0.36
# DOGRUSAL eksen: log olcek oran farklarini okunaksiz kiliyordu (2x10^1 gibi
# etiketler). A0C zaten belirgin sekilde ayri; dogrudan etiket yeterli.
b1 = ax[0].bar(x - w / 2, g4, w, color=C[0], zorder=3, label="loop kapalı")
b2 = ax[0].bar(x + w / 2, g9, w, color=C[1], zorder=3, label="loop açık")
for bs in (b1, b2):
    for r_ in bs:
        h = r_.get_height()
        if np.isfinite(h):
            ax[0].text(r_.get_x() + r_.get_width() / 2, h + 2.5, f"{h:.0f}",
                       ha="center", fontsize=5.6, color=INK)
ax[0].set_ylabel("tutarlılık kazancı (×)")
ax[0].set_xticks(x); ax[0].set_xticklabels(alab)
ax[0].set_ylim(0, 118)
ax[0].legend(loc="upper left", fontsize=6)
ax[0].set_title("(a) A0'a göre NEES oranı", fontsize=7)
tidy(ax[0])
b3 = ax[1].bar(x, c4, 0.6, color=C[2], zorder=3)
for r_ in b3:
    h = r_.get_height()
    ax[1].text(r_.get_x() + r_.get_width() / 2, h + 0.012, f"{h:.2f}",
               ha="center", fontsize=5.6, color=INK)
ax[1].axhline(1.0, color=INK2, lw=0.7, ls="--", zorder=2)
ax[1].set_ylabel("ATE oranı (A0 = 1,0)")
ax[1].set_xticks(x); ax[1].set_xticklabels(alab)
ax[1].set_ylim(0.85, 1.13)
ax[1].set_title("(b) doğruluk bedeli", fontsize=7)
tidy(ax[1])
fig.tight_layout(pad=0.35)
fig.savefig(f"{OUT}/fig4_kollar.pdf"); plt.close(fig)

print("uretilen sekiller:")
for f in sorted(os.listdir(OUT)):
    print(f"  {OUT}/{f}  ({os.path.getsize(os.path.join(OUT,f))} bayt)")
