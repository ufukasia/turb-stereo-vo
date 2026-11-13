"""Sapma 5 tasarım duyarlılığı: ana ızgaranın (bilinen veri) birim log-oran saçılımından, n = 36 için
Holm'un en sıkı basamağında (alfa / m) %80 güçle saptanabilir en küçük medyan oran. Sonuçtan önce hesaplandı."""
import sys, numpy as np
from scipy import stats
sys.path.insert(0, ".")
from turb.faz2_analiz import oku, eslesik
D = oku("evidence/FAZ2/records.jsonl")
N = 36; ARE = 0.955
def mde(sd, m, n=N, guc=0.80):
    za = stats.norm.ppf(1 - 0.05 / m / 2); zb = stats.norm.ppf(guc)
    d = (za + zb) * sd / np.sqrt(n * ARE)
    return float(np.exp(d))
def guc(sd, etki_log, m, n=N):
    za = stats.norm.ppf(1 - 0.05 / m / 2)
    return float(stats.norm.cdf(abs(etki_log) * np.sqrt(n * ARE) / sd - za))
def sd_log(kol, ref, dr, metrik, dr_ref=None):
    a, b, U, _ = eslesik(D, kol, ref, dr, metrik, dr_ref)
    lr = np.log(np.asarray(a) / np.asarray(b))
    # sağlam saçılım (MAD) ve klasik SD birlikte
    return len(lr), float(np.std(lr, ddof=1)), float(1.4826 * np.median(np.abs(lr - np.median(lr)))), float(np.median(lr))
print("karşılaştırma | n_ana | SD(log) | MAD-SD | ana medyan oran | n=36 MDE oran (SD; Holm m) | ana etki için güç")
for ad, kol, ref, metrik, m in (("F1 A1/A0 ANEES", "A1", "A0", "ANEES_poz", 10), ("F1 A3/A0 ANEES", "A3", "A0", "ANEES_poz", 10),
                                ("F1 AG/A0 ANEES", "AG", "A0", "ANEES_poz", 10), ("F1 A3/A0 ATE", "A3", "A0", "ATE", 10),
                                ("F1 A1/A0 ATE", "A1", "A0", "ATE", 10), ("F1 AG/A0 ATE", "AG", "A0", "ATE", 10),
                                ("F1 A0C/A0 ATE", "A0C", "A0", "ATE", 10), ("F2 A3/AG ANEES", "A3", "AG", "ANEES_poz", 2)):
    for dr in (4.0, 8.0):
        n, sd, mad, med = sd_log(kol, ref, dr, metrik)
        print(f"{ad} D/r0 {dr:g} | {n} | {sd:.3f} | {mad:.3f} | {np.exp(med):.3f} | {mde(sd, m):.3f} (1/{mde(sd, m):.3f} = {1/mde(sd, m):.3f}) | {guc(sd, med, m):.2f}")
for dr in (4.0, 8.0):
    n, sd, mad, med = sd_log("A0", "A0", dr, "ATE", dr_ref=0.0)
    print(f"bozulma A0 ATE D/r0 {dr:g}/temiz | {n} | {sd:.3f} | {mad:.3f} | {np.exp(med):.3f} | {mde(sd, 1):.3f} | {guc(sd, med, 1):.2f}")
