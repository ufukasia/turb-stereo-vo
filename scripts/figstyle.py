"""Shared matplotlib style for README figures (copy next to the repo's figure script and import it).

Every PNG: exactly 10.0 in x 220 dpi = 2200 px wide, white, DejaVu Sans, one legend top-left, no suptitle.
Heights only from HEIGHTS. save() asserts the width. Adapt SYS (series -> colour/dash/width/label) per project,
keeping the meaning of the colours: proposed = blue solid, rivals = orange/vermillion/green/pink/grey with dashes.
"""
from pathlib import Path

import numpy as np

W, DPI = 10.0, 220
HEIGHTS = {"row": 3.2, "two_rows": 5.6, "hero": 7.6}
GRID = "#EEEEEE"
GT = "#1A1A1A"
# key: (colour, linestyle, linewidth, default label)
SYS = {
    "gt":        (GT,        "-",                    2.8, "Ground truth"),
    "proposed":  ("#0072B2", "-",                    2.0, "Proposed"),
    "variant":   ("#56B4E9", (0, (4, 2)),            1.5, "Proposed (variant)"),
    "rival1":    ("#E69F00", (0, (5, 2)),            1.5, "Rival 1"),
    "rival2":    ("#D55E00", (0, (1.5, 1.5)),        1.5, "Rival 2 / scaled"),
    "rival3":    ("#009E73", (0, (6, 2, 1.5, 2)),    1.5, "Rival 3"),
    "rival4":    ("#CC79A7", (0, (6, 2, 1, 2, 1, 2)), 1.5, "Rival 4"),
    "rival5":    ("#6E6E6E", (0, (8, 3)),            1.5, "Rival 5"),
}
NEUTRAL_BAR = "#BBBBBB"


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 12, "axes.labelsize": 12, "xtick.labelsize": 11,
        "ytick.labelsize": 11, "legend.fontsize": 11.5, "axes.titlesize": 12.5, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "axes.titlepad": 8, "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#888888", "axes.linewidth": 0.8, "axes.labelcolor": "#222222", "axes.axisbelow": True,
        "xtick.color": "#444444", "ytick.color": "#444444", "xtick.direction": "out", "ytick.direction": "out",
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
        "legend.handlelength": 2.6, "legend.columnspacing": 1.6, "figure.facecolor": "white",
        "savefig.facecolor": "white", "savefig.transparent": False, "savefig.bbox": "standard",
        "lines.solid_capstyle": "round", "lines.dash_capstyle": "butt"})
    return plt


def new_fig(plt, h):
    assert h in HEIGHTS.values(), f"height {h} not in {HEIGHTS}"
    return plt.figure(figsize=(W, h), dpi=DPI)


def legend(fig, handles, ncol, y=1.0):
    """The ONE legend of the figure: top-left, horizontal, frameless. Leave room with gridspec top≈0.88."""
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, y), ncol=ncol, borderaxespad=0.4)


def handle(key, label=None, marker=None):
    from matplotlib.lines import Line2D
    c, ls, lw, lab = SYS[key]
    return Line2D([], [], color=c, ls=ls, lw=lw, marker=marker, ms=6, mec="white", mew=0.6, label=label or lab)


def plot(ax, key, x, y, **kw):
    c, ls, lw, _ = SYS[key]
    z = 1 if key == "gt" else (4 if key == "proposed" else 3)
    return ax.plot(x, y, color=c, ls=ls, lw=lw, zorder=z, **kw)


def save(fig, out_dir, name):
    import matplotlib.pyplot as plt
    from PIL import Image
    p = Path(out_dir) / name
    fig.savefig(p, dpi=DPI, metadata={"Software": "Matplotlib"})   # never bbox_inches="tight"
    plt.close(fig)
    w, h = Image.open(p).size
    assert w == int(W * DPI), (name, w)
    print(f"{name}: {w} x {h}")
    return p


def topview(ax, g, labels="xyz"):
    """Limits from ground truth with equal padding; equal aspect unless the path is elongated (ratio >= 6).
    Also fixes 'flat' panels: the short axis gets its own padding instead of a huge empty frame."""
    span = g.max(0) - g.min(0)
    ex = np.argsort(span)[::-1][:2]
    lo, hi = g.min(0), g.max(0)
    equal = span[ex[0]] / max(span[ex[1]], 1e-6) < 6
    px = 0.25 * max(span[ex[0]], span[ex[1]]) if equal else 0.06 * span[ex[0]]
    py = px if equal else 0.6 * span[ex[1]] + 0.05
    ax.set_xlim(lo[ex[0]] - px, hi[ex[0]] + px)
    ax.set_ylim(lo[ex[1]] - py, hi[ex[1]] + py)
    if equal:
        ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(f"{labels[ex[0]]} [m]" + ("" if equal else "  (axes not to scale)"))
    ax.set_ylabel(f"{labels[ex[1]]} [m]")
    return ex


def gaps(t, *ys, limit=1.0):
    """Insert NaN where timestamps jump more than `limit` s, so lost segments are not bridged by a line."""
    cut = np.where(np.diff(t) > limit)[0] + 1
    out = [np.insert(np.asarray(t, float), cut, np.nan)]
    for y in ys:
        out.append(np.insert(np.asarray(y, float), cut, np.nan, axis=0))
    return out


def check_all(fig_dir):
    """Assert every PNG in the figure folder is 2200 px wide."""
    from PIL import Image
    bad = {p.name: Image.open(p).size[0] for p in Path(fig_dir).glob("*.png") if Image.open(p).size[0] != 2200}
    assert not bad, bad
