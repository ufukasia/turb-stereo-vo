"""
Eşleştirilmiş istatistik protokolü — ÖNCEKİ MAKALEDEN DEĞİŞTİRİLMEDEN DEVRALINIR
(rapor §6.0).

Kabul ölçütü ÜÇLÜSÜ (üçü birden sağlanmalı):
    (i)   eşleştirilmiş farkların bootstrap %95 güven aralığı sıfırı KESMEYECEK
    (ii)  Wilcoxon işaretli sıra testi p < 0.05
    (iii) işaret testi aynı yönü gösterecek

Bu ölçütler ÖN-KAYITLIDIR: sonuçlara bakılarak gevşetilmez.

Monte-Carlo notu: hat belirlenimcidir (koşum-içi varyans = 0). Bu deneyde
belirsizliğin TEK kaynağı türbülans gerçeklemesidir; bu yüzden hücre başına M
tohum, doğrudan bozulmanın etkisinin dağılımını verir — çoğu SLAM çalışmasının
ayıramadığı iki gürültü kaynağı burada ayrıktır.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy import stats


@dataclass
class PairedResult:
    n: int
    mean_diff: float
    median_diff: float
    ci_low: float
    ci_high: float
    ci_excludes_zero: bool
    wilcoxon_p: float
    wilcoxon_ok: bool
    sign_p: float
    n_improved: int
    n_worsened: int
    sign_ok: bool
    accepted: bool
    direction: str
    rel_improvement: float

    def to_dict(self):
        return asdict(self)


def paired_test(baseline: np.ndarray, treatment: np.ndarray, n_boot: int = 10000,
                alpha: float = 0.05, seed: int = 12345,
                lower_is_better: bool = True) -> PairedResult:
    """Eşleştirilmiş karşılaştırma: treatment - baseline farkları üzerinde üçlü ölçüt.

    lower_is_better=True (ATE gibi) ise "iyileşme" = negatif fark.
    """
    b = np.asarray(baseline, dtype=np.float64)
    t = np.asarray(treatment, dtype=np.float64)
    m = np.isfinite(b) & np.isfinite(t)
    b, t = b[m], t[m]
    n = b.size
    if n < 3:
        return PairedResult(n, *(np.nan,) * 3, np.nan, False, np.nan, False,
                            np.nan, 0, 0, False, False, "yetersiz", np.nan)
    d = t - b

    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = d[idx].mean(axis=1)
    lo, hi = np.percentile(boot, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    ci_ok = (lo > 0) or (hi < 0)

    try:
        w_p = float(stats.wilcoxon(d, zero_method="wilcox").pvalue)
    except ValueError:
        w_p = 1.0
    n_imp = int(np.sum(d < 0)) if lower_is_better else int(np.sum(d > 0))
    n_wor = int(np.sum(d > 0)) if lower_is_better else int(np.sum(d < 0))
    n_eff = n_imp + n_wor
    s_p = float(stats.binomtest(n_imp, n_eff, 0.5).pvalue) if n_eff > 0 else 1.0
    sign_ok = (n_imp > n_wor) and (s_p < alpha)

    improved = (d.mean() < 0) if lower_is_better else (d.mean() > 0)
    accepted = bool(ci_ok and (w_p < alpha) and sign_ok and improved)
    rel = float(d.mean() / np.mean(b)) if np.mean(b) != 0 else np.nan

    return PairedResult(
        n=n, mean_diff=float(d.mean()), median_diff=float(np.median(d)),
        ci_low=float(lo), ci_high=float(hi), ci_excludes_zero=bool(ci_ok),
        wilcoxon_p=w_p, wilcoxon_ok=bool(w_p < alpha),
        sign_p=s_p, n_improved=n_imp, n_worsened=n_wor, sign_ok=bool(sign_ok),
        accepted=accepted, direction=("iyileşme" if improved else "kötüleşme"),
        rel_improvement=rel)


def bootstrap_ci(x: np.ndarray, n_boot: int = 10000, alpha: float = 0.05,
                 seed: int = 12345, stat=np.mean):
    """Tek örneklem bootstrap güven aralığı."""
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(n_boot, x.size))
    boot = stat(x[idx], axis=1)
    lo, hi = np.percentile(boot, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(stat(x)), float(lo), float(hi)


def summarize_cells(values_by_cell: dict, n_boot: int = 10000):
    """{hücre: [tohum değerleri]} -> hücre başına ortalama + %95 GA."""
    out = {}
    for k, v in values_by_cell.items():
        mean, lo, hi = bootstrap_ci(np.asarray(v, dtype=np.float64), n_boot=n_boot)
        out[k] = {"n": int(np.sum(np.isfinite(v))), "mean": mean,
                  "ci_low": lo, "ci_high": hi,
                  "median": float(np.nanmedian(v)) if len(v) else np.nan}
    return out


def holm_bonferroni(pvals: dict, alpha: float = 0.05) -> dict:
    """Çoklu karşılaştırma düzeltmesi (Holm-Bonferroni)."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    n = len(items)
    out, prev = {}, 0.0
    for i, (k, p) in enumerate(items):
        thr = alpha / (n - i)
        prev = max(prev, p)
        out[k] = {"p": p, "threshold": thr, "reject": bool(prev < thr)}
    return out
