"""
Derinliğe bağlı türbülans modeli — fizik ve uygulama doğrulama kapısı (uzun koşumlardan ÖNCE).

    python -m turb.test_derinlik_model [--sekans /data/sdata_airsim/AF_F1_h5]
        [--dtm ../airsim_sim/africa/harita_africa.npz] [--json cikti.json]

Testler (beklenen değerler koşumdan önce TASARIM.md §5.3 / ON_TEST dosyasında yazılı):
  T0  C0 = 0: model alanları sıfır; üretici çıktısı bayt bayt aynı
  T1  ekran gerçeklemesi: bant-sınırlı warp varyansı (ampirik) = analitik σ²_düşük  (±%10)
  T2  homojen yol: σ²_tot ∝ L (log-log eğim 1 ± 0.05)
  T3  homojen yol: σ_tot = optics.sigma_tilt_rad(D/r0_küre) ile birebir (±%2, kareleme)
  T4  perde etkisi: Cn² yalnız kameraya yakın %10'luk dilimde / nesneye yakın %10'luk dilimde
      σ²_tot oranı = ∫(z/L)^(5/3) ağırlık oranı (±%5)
  T5  kısa-poz bulanıklık: Airy FWHM 1.029 λ/D; uzun-poz seeing 0.98 λ/r0; SE/LE(D/r0=10)
      0.55–0.80; ek varyans D/r0 -> 0 iken 0 ve monoton
  T6  gerçek kare (F rejimi): σ_tot derinlik dilimleriyle artar (Spearman ρ ≥ 0.9)
  T7  görüntü düzeyi: temiz↔türbülanslı optik akış RMS'i ≈ modelin warp RMS'i (oran 0.75–1.33),
      akış–model korelasyonu ≥ 0.7 (işaret ve ölçek doğrulaması)
  T8  belirlenimcilik: aynı tohum -> aynı sha256
  T9  stereo: cam0–cam1 warp korelasyonu (rapor; 0 ile 1 arası beklenir)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import tempfile
from pathlib import Path

import cv2
import numpy as np

from . import optics
from .derinlik_model import Arazi, DerinlikModeli, DerinlikParam, se_bulaniklik_var_rad2, uygula

R_ILERI = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])   # kamera -> NED, kuzeye bakış


def _izgara(w=256, h=192, fxl=972.0):
    u, v = np.meshgrid(np.arange(w, dtype=float), np.arange(h, dtype=float))
    return (u - (w - 1) / 2) / fxl, (v - (h - 1) / 2) / fxl


def _model(C0=1e-12, **kw):
    p = DerinlikParam(C0=C0, **kw)
    return DerinlikModeli(p, fx=3889.0, arazi=Arazi(None, 0.0))


def T0_sifir():
    m = _model(C0=0.0)
    ex, ey = _izgara()
    Z = np.full(ex.shape, 100.0)
    dx, dy, b, _ = m.kare(Z, ex, ey, R_ILERI, np.array([0, 0, -5.0]), np.zeros((m.p.M, 2)), np.zeros(m.p.M))
    return {"maks_alan": float(max(np.abs(dx).max(), np.abs(dy).max(), b.max())),
            "gecti": bool(max(np.abs(dx).max(), np.abs(dy).max(), b.max()) == 0.0)}


def T1_ekran(n_deneme=40):
    m = _model(C0=1e-12)
    m.cn_fn = lambda h, a, Z: np.ones_like(h)                 # homojen
    ex, ey = _izgara()
    Z = np.full(ex.shape, 200.0)
    rng = np.random.default_rng(1)
    emp, ana = [], []
    for _ in range(n_deneme):
        of = rng.uniform(0, 1, (m.p.M, 2)) * (m.p.N_ekran * m.delta)[:, None]
        dx, dy, _, t = m.kare(Z, ex, ey, R_ILERI, np.array([0, 0, -5.0]), of, rng.uniform(0, 6.28, m.p.M))
        emp.append(0.5 * (np.mean(dx ** 2) + np.mean(dy ** 2))); ana.append(np.mean(t["sigma_lo_px"] ** 2))
    oran = float(np.mean(emp) / np.mean(ana))
    return {"ampirik_var_px2": float(np.mean(emp)), "analitik_var_px2": float(np.mean(ana)),
            "oran": oran, "gecti": bool(abs(oran - 1) < 0.10)}


def T2_T3_homojen():
    m = _model(C0=1e-13)
    m.cn_fn = lambda h, a, Z: np.ones_like(h)
    ex, ey = np.zeros((1, 1)), np.zeros((1, 1))
    L = np.array([10, 30, 100, 300, 1000.0]); s2 = []; fark = []
    for Li in L:
        _, _, _, r0s = m.hat_integralleri(np.full((1, 1), Li), ex, ey, R_ILERI, np.array([0, 0, -5.0]))
        var = 0.170 * m.p.lam ** 2 * m.p.D ** (-1 / 3) * r0s[0, 0]
        s2.append(var)
        r0 = (0.423 * m.k ** 2 * m.p.C0 * 3 * Li / 8) ** (-3 / 5)        # ∫(z/L)^{5/3}dz = 3L/8
        ref = optics.sigma_tilt_rad(m.p.D / r0, m.p.lam, m.p.D)
        fark.append(math.sqrt(var) / ref - 1)
    egim = float(np.polyfit(np.log(L), np.log(s2), 1)[0])
    return ({"egim": egim, "gecti": bool(abs(egim - 1) < 0.05)},
            {"maks_bagil_fark": float(np.max(np.abs(fark))), "gecti": bool(np.max(np.abs(fark)) < 0.02)})


def T4_perde():
    m = _model(C0=1e-13, alt_ornek=256)          # süreksiz göstergeç profili: ince kareleme
    ex, ey = np.zeros((1, 1)), np.zeros((1, 1)); L = 200.0
    out = {}
    for ad, (a0, a1) in {"kamera": (0.0, 0.1), "nesne": (0.9, 1.0)}.items():
        m.cn_fn = lambda h, a, Z, a0=a0, a1=a1: ((a >= a0 * Z) & (a < a1 * Z)).astype(float)
        _, _, _, r0s = m.hat_integralleri(np.full((1, 1), L), ex, ey, R_ILERI, np.array([0, 0, -5.0]))
        out[ad] = float(r0s[0, 0])
    w = lambda s0, s1: ((1 - s0) ** (8 / 3) - (1 - s1) ** (8 / 3)) * 3 / 8     # ∫ (1-s)^{5/3} ds
    beklenen = w(0.0, 0.1) / w(0.9, 1.0)
    oran = out["kamera"] / out["nesne"]
    return {"oran": oran, "beklenen": beklenen, "gecti": bool(abs(oran / beklenen - 1) < 0.05)}


def T4b_kareleme(seq, dtm):
    """Üretim ayarı (alt_ornek=4) gerçek F karesinde, ince kareleme (64) referansına göre:
    σ_tot ve W_m toplamı medyan bağıl hata < %3."""
    from . import uret_derinlik as U
    S = U.sekans_oku(Path(seq)); i = 600; ts = S["ts"][i]
    Z = U._derinlik(S["m"] / "depth0/data" / f"{ts}.npy")
    h, w = Z.shape
    ex, ey, _ = U._egim(S["fx"], S["fy"], S["cx"], S["cy"], w, h, S["W"], S["H"])
    R, t = U._R(S["gt"][i, 4:8]), S["gt"][i, 1:4]
    res = {}
    for K in (4, 64):
        m = DerinlikModeli(DerinlikParam(C0=1e-12, alt_ornek=K), S["fx"], Arazi(dtm, 0.0))
        Wm, _, _, r0s = m.hat_integralleri(Z, ex, ey, R, t)
        res[K] = (r0s, Wm.sum(0))
    ok = np.isfinite(Z)
    e1 = float(np.median(np.abs(res[4][0][ok] / res[64][0][ok] - 1)))
    e2 = float(np.median(np.abs(res[4][1][ok] / res[64][1][ok] - 1)))
    return {"r0kure_medyan_bagil_hata": e1, "Wtoplam_medyan_bagil_hata": e2, "gecti": bool(e1 < 0.03 and e2 < 0.03)}


def T5_bulaniklik():
    """Sayısal PSF/FWHM zincirini kanonik değerlerle doğrular, sonra kısa-poz eki:
    (a) Airy FWHM = 1.029 λ/D (±%2); (b) uzun-poz seeing FWHM ≈ 0.98 λ/r0 (D/r0 = 30, ±%5);
    (c) kısa-poz / uzun-poz FWHM oranı D/r0 = 10'da 0.55–0.80 (moment tahmini ≈ 0.71);
    (d) ek varyans D/r0 -> 0 iken -> 0 ve monoton."""
    from .derinlik_model import _FWHM_KIRINIM, _fwhm
    lam, D = 5.5e-7, 4.8e-3
    q = np.array([1e-3, 0.1, 0.3, 1, 3, 10])
    s = np.sqrt(se_bulaniklik_var_rad2(q, lam, D)) / (lam / D)
    nu = np.linspace(0.0, 1.0, 4001); r = np.linspace(0, 40.0, 16001)
    mtf_d = (2 / np.pi) * (np.arccos(nu) - nu * np.sqrt(np.clip(1 - nu * nu, 0, 1)))
    le30 = _fwhm(mtf_d * np.exp(-3.44 * 30 ** (5 / 3) * nu ** (5 / 3)), nu, r) / 30      # λ/r0 birimi
    le10 = _fwhm(mtf_d * np.exp(-3.44 * 10 ** (5 / 3) * nu ** (5 / 3)), nu, r)
    se10 = math.sqrt(_FWHM_KIRINIM ** 2 + 8 * math.log(2) * se_bulaniklik_var_rad2(np.array([10.0]), lam, D)[0] / (lam / D) ** 2)
    oran = se10 / le10
    return {"sigma_ek_lamD": s.round(4).tolist(), "FWHM_kirinim_lamD": round(float(_FWHM_KIRINIM), 4),
            "LE_FWHM_lam_r0_D_r0_30": round(float(le30), 4), "SE_LE_orani_D_r0_10": round(float(oran), 3),
            "gecti": bool(abs(_FWHM_KIRINIM / 1.029 - 1) < 0.02 and abs(le30 / 0.98 - 1) < 0.05
                          and 0.55 < oran < 0.80 and s[0] < 1e-3 and np.all(np.diff(s) > 0))}


def _gercek_kare(seq, dtm, C0, i=600, tohum=0):
    from . import uret_derinlik as U
    S = U.sekans_oku(Path(seq))
    p = DerinlikParam(C0=C0, tohum=tohum)
    m = DerinlikModeli(p, S["fx"], Arazi(dtm, 0.0))
    ts = S["ts"][i]
    Z = U._derinlik(S["m"] / "depth0/data" / f"{ts}.npy")
    h, w = Z.shape
    ex, ey, fxl = U._egim(S["fx"], S["fy"], S["cx"], S["cy"], w, h, S["W"], S["H"])
    R, t = U._R(S["gt"][i, 4:8]), S["gt"][i, 1:4]
    of = np.zeros((p.M, 2)); km = np.zeros(p.M)
    r0 = m.kare(Z, ex, ey, R, t, of, km)
    Z1 = U._cam1_derinlik(Z, fxl, S["B"])
    r1 = m.kare(Z1, ex, ey, R, t, of, km, ofset_x=S["B"])
    img = cv2.imread(str(S["m"] / f"cam0/data/{ts}.png"), cv2.IMREAD_UNCHANGED)
    return S, Z, r0, r1, img


def T6_derinlik(seq, dtm):
    _, Z, (dx, dy, b, t), _, _ = _gercek_kare(seq, dtm, 1e-12)
    from scipy.stats import spearmanr
    kenar = [0, 25, 50, 100, 200, 400, 1e9]; med = []; mer = []
    for a, c in zip(kenar[:-1], kenar[1:]):
        mm = np.isfinite(Z) & (Z >= a) & (Z < c)
        if mm.sum() > 50:
            med.append(float(np.median(t["sigma_tot_px"][mm]))); mer.append(a)
    rho = float(spearmanr(mer, med)[0]) if len(med) > 2 else float("nan")
    return {"dilim_alt_m": mer, "sigma_tot_medyan_px": [round(x, 4) for x in med], "spearman": rho,
            "gecti": bool(rho >= 0.9)}


def T7_goruntu(seq, dtm, C0=1e-11):
    S, Z, (dx, dy, b, t), _, img = _gercek_kare(seq, dtm, C0)
    H, W = img.shape
    # yalnız warp (bulanıklık ayrı): akış ölçümü bulanıklıktan etkilenmesin diye blur = 0
    tur = uygula(img, dx, dy, np.zeros_like(b))
    flow = cv2.calcOpticalFlowFarneback(img, tur, None, 0.5, 4, 21, 5, 7, 1.5, 0)
    DX = cv2.resize(dx.astype(np.float32), (W, H), interpolation=cv2.INTER_CUBIC)
    DY = cv2.resize(dy.astype(np.float32), (W, H), interpolation=cv2.INTER_CUBIC)
    # dokulu bölgeler (akış güvenilir) ve kenar dışı
    gx = cv2.Sobel(img, cv2.CV_32F, 1, 0); gy = cv2.Sobel(img, cv2.CV_32F, 0, 1)
    dok = np.hypot(gx, gy) > 40
    dok[:30] = dok[-30:] = False; dok[:, :30] = dok[:, -30:] = False
    f = np.concatenate([flow[..., 0][dok], flow[..., 1][dok]])
    mdl = np.concatenate([DX[dok], DY[dok]])
    oran = float(np.sqrt(np.mean(f ** 2)) / np.sqrt(np.mean(mdl ** 2)))
    kor = float(np.corrcoef(f, mdl)[0, 1])
    return {"C0": C0, "model_warp_rms_px": float(np.sqrt(np.mean(mdl ** 2))), "akis_rms_px": float(np.sqrt(np.mean(f ** 2))),
            "oran": oran, "korelasyon": kor, "blur_medyan_px": float(np.median(b)),
            "gecti": bool(0.75 < oran < 1.33 and kor >= 0.7)}


def T8_T0_uretici(seq, dtm):
    from . import uret_derinlik as U
    sha = []
    with tempfile.TemporaryDirectory() as td:
        for k in range(2):
            out = U.main(["--girdi", seq, "--cikti_kok", td, "--C0", "1e-12", "--tohum", "3",
                          "--dtm", dtm, "--kareler", "600:606", "--isci", "3", "--ad", f"t{k}"])
            h = hashlib.sha256()
            for f in sorted((out / "mav0/cam0/data").glob("*.png")) + sorted((out / "mav0/cam1/data").glob("*.png")):
                h.update(f.read_bytes())
            sha.append(h.hexdigest())
        out0 = U.main(["--girdi", seq, "--cikti_kok", td, "--C0", "0", "--dtm", dtm, "--kareler", "600:606",
                       "--isci", "2", "--ad", "s0"])
        ayni = all((out0 / f"mav0/{c}/data" / f.name).read_bytes() == f.read_bytes()
                   for c in ("cam0", "cam1") for f in (Path(seq) / f"mav0/{c}/data").glob("*.png")
                   if (out0 / f"mav0/{c}/data" / f.name).exists())
        n0 = len(list((out0 / "mav0/cam0/data").glob("*.png")))
    return ({"sha_1": sha[0][:16], "sha_2": sha[1][:16], "gecti": sha[0] == sha[1]},
            {"C0_0_bayt_ayni": bool(ayni and n0 == 6), "gecti": bool(ayni and n0 == 6)})


def T9_stereo(seq, dtm):
    _, Z, (dx0, dy0, _, _), (dx1, dy1, _, _), _ = _gercek_kare(seq, dtm, 1e-12)
    ok = np.isfinite(Z)
    return {"cam0_cam1_dx_korelasyon": float(np.corrcoef(dx0[ok], dx1[ok])[0, 1]),
            "cam0_cam1_dy_korelasyon": float(np.corrcoef(dy0[ok], dy1[ok])[0, 1]), "gecti": True}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sekans", default="/data/sdata_airsim/AF_F1_h5")
    ap.add_argument("--dtm", default=str(Path(__file__).resolve().parents[2] / "airsim_sim/africa/harita_africa.npz"))
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    R = {}

    class _Yaz(dict):
        def __setitem__(self, k, v):
            super().__setitem__(k, v)
            print(f"  ... {k}: {'geçti' if v['gecti'] else 'KALDI'}", flush=True)
    R = _Yaz()
    R["T0_model_sifir"] = T0_sifir()
    R["T1_ekran_varyansi"] = T1_ekran()
    R["T2_L_olcekleme"], R["T3_kapali_form"] = T2_T3_homojen()
    R["T4_perde_etkisi"] = T4_perde()
    R["T4b_kareleme"] = T4b_kareleme(a.sekans, a.dtm)
    R["T5_kisa_poz"] = T5_bulaniklik()
    R["T6_derinlik_artisi"] = T6_derinlik(a.sekans, a.dtm)
    R["T7_goruntu_akis"] = T7_goruntu(a.sekans, a.dtm)
    R["T8_belirlenimcilik"], R["T0b_uretici_C0_0"] = T8_T0_uretici(a.sekans, a.dtm)
    R["T9_stereo_korelasyon"] = T9_stereo(a.sekans, a.dtm)
    for k, v in R.items():
        print(f"{'GEÇTİ ' if v['gecti'] else 'KALDI '} {k}: " + json.dumps({x: y for x, y in v.items() if x != 'gecti'}))
    if a.json:
        json.dump(R, open(a.json, "w"), indent=1)
    print("HEPSİ GEÇTİ" if all(v["gecti"] for v in R.values()) else "BAŞARISIZ TEST VAR")


if __name__ == "__main__":
    main()
