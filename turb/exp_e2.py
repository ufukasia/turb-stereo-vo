"""
E2 — Ön-uç karakterizasyonu ve korelasyon modeli doğrulaması (rapor §6/E2).

"Bu deney projenin bilimsel kalbidir." Bozulma sentetik olduğu için TEMİZ görüntü
elimizdedir; dolayısıyla her özniteliğin GERÇEK yer değiştirmesi bilinir ve teorik
model ampirik olarak sınanabilir.

Ölçülenler (rapor §7.3):
    1. Öznitelik konum hatası dağılımı (ölçülen - GT eğim alanı)
    2. Ampirik açısal korelasyon fonksiyonu rho_hat(dtheta)
    3. Dedektör tekrarlanabilirliği
    4. Eşleştirme skoru / iç-değer oranı
    5. Betimleyici mesafe dağılımı

Kritik ayrım — DAİRESELLİK YOK: rho_hat, simülatörün korelasyon MODELİNDEN değil,
ALIKED+LightGlue'nun bozuk görüntüde ölçtüğü GERÇEK öznitelik kaymalarından
hesaplanır. Karşılaştırılan iki model:
    (M1) iso53    : exp[-(dtheta/theta_c)^{5/3}]        (rapor §4.2 mühendislik modeli)
    (M2) physical : Kolmogorov açıklık-ortalamalı tilt tensörü (turb/correlation.py)

    python -m turb.exp_e2 --sequences MH01 MH03 --D_over_r0 1 3 5 8 --n_frames 40
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import torch

from . import correlation as _corr
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


# =============================================================================
def forward_displacement(sim: TurbulenceSimulator, frame: int, cam: int,
                         u: np.ndarray, v: np.ndarray, iters: int = 3):
    """TEMİZ (u,v) konumundaki sahne noktasının bozuk görüntüdeki konumu.

    Simülatör geri-yönlü harita kullanır: out(x) = in(x - d(x)). Dolayısıyla temiz
    konum u_c, bozuk görüntüde u_d = u_c + d(u_d) noktasında görünür; bu örtük
    denklem sabit-nokta yinelemesiyle çözülür (alan düzgün olduğundan 2-3 adım yeter).
    """
    ud, vd = np.array(u, dtype=np.float64), np.array(v, dtype=np.float64)
    for _ in range(iters):
        du, dv = sim.displacement_at(frame, cam, ud, vd)
        ud = u + du
        vd = v + dv
    return ud - u, vd - v


def _extract(matcher, img):
    """Ön-uçtaki dedektörle keypoint + descriptor çıkar."""
    t = torch.from_numpy(img / 255.).float()[None, None].to(matcher.device)
    h, w = t.shape[-2:]
    data = {"image": t, "image_size": t.new_tensor([[w, h]])}
    with torch.no_grad():
        f = matcher.aliked(data)
    return (f["keypoints"][0].cpu().numpy(),
            f["descriptors"][0].cpu().numpy(),
            {k: v for k, v in f.items()})


def _match(matcher, feats_a, feats_b, shape):
    """LightGlue ile eşleştir (ön-uçtaki stereo eşleştiricinin aynısı)."""
    h, w = shape
    d0 = {k: (v if torch.is_tensor(v) else torch.as_tensor(v)) for k, v in feats_a.items()}
    d1 = {k: (v if torch.is_tensor(v) else torch.as_tensor(v)) for k, v in feats_b.items()}
    for d in (d0, d1):
        d["image_size"] = torch.tensor([[w, h]], dtype=torch.float32, device=matcher.device)
    with torch.no_grad():
        out = matcher.lightglue_stereo({"image0": d0, "image1": d1})
    m = out["matches0"][0].cpu().numpy()
    sc = out["matching_scores0"][0].cpu().numpy()
    return m, sc


# =============================================================================
def run(sequences, base_path, D_over_r0_list, theta_c_list, seeds, n_frames,
        frame_stride, blur, outdir, max_kp=1024, det_thr=None):
    import config
    from stereo_front_end import StereoMatcher

    os.makedirs(outdir, exist_ok=True)
    cfg = config.load()
    fcfg = cfg["frontend"]
    rows = []

    for seq in sequences:
        matcher = StereoMatcher(base_path=base_path, sequence=seq,
                                matcher_type="lightglue", filter_method=fcfg["filter"],
                                y_axis_threshold=fcfg["y_axis_threshold"],
                                extractor=fcfg["extractor"], max_track_length=0,
                                max_keypoints=int(max_kp),
                                detection_threshold=det_thr or fcfg["detection_threshold"])
        # SAYISAL sırala (bkz. stereo_front_end.process_sequence): sözlüksel sıralama
        # karışık-basamaklı damgalarda (TartanAir) kare sırasını bozar; sabit-hane
        # damgalarda (EuRoC) sonuç birebir aynıdır.
        ts_all = [str(t) for t in sorted({int(t) for t in matcher.cam0_df["timestamp"]})]
        import cv2
        for dr in D_over_r0_list:
            for tc in theta_c_list:
                for seed in seeds:
                    p = TurbParams(enabled=True, seed=int(seed), D_over_r0=float(dr),
                                   theta_c_over_fov=float(tc), tau0_over_dt=1.0,
                                   blur=blur)
                    sim = None
                    t0 = time.time()
                    acc = {"err": [], "gt_mag": [], "meas_mag": [], "rep": [], "ms": [],
                           "desc": [], "n_clean": [], "n_dist": []}
                    pair_sep, pair_prod_l, pair_prod_t = [], [], []
                    var_acc = []
                    for fi in range(0, n_frames * frame_stride, frame_stride):
                        if fi >= len(ts_all):
                            break
                        ts = ts_all[fi]
                        fn = matcher.cam0_df.loc[
                            matcher.cam0_df["timestamp"] == ts, "filename"].values[0]
                        img = cv2.imread(str(matcher.base_path / seq / "mav0/cam0/data" / fn),
                                         cv2.IMREAD_GRAYSCALE)
                        if img is None:
                            continue
                        clean = cv2.remap(img, matcher.map1x, matcher.map1y, cv2.INTER_LINEAR)
                        if sim is None:
                            sim = TurbulenceSimulator(p, image_size=(clean.shape[1],
                                                                     clean.shape[0]))
                        dist = sim.apply(clean, fi, cam=0)

                        kp_c, ds_c, f_c = _extract(matcher, clean)
                        kp_d, ds_d, f_d = _extract(matcher, dist)
                        if len(kp_c) < 10 or len(kp_d) < 10:
                            continue
                        acc["n_clean"].append(len(kp_c))
                        acc["n_dist"].append(len(kp_d))

                        # GT yer değiştirme (temiz konumlarda, ileri yön)
                        gdu, gdv = forward_displacement(sim, fi, 0, kp_c[:, 0], kp_c[:, 1])

                        m, sc = _match(matcher, f_c, f_d, clean.shape)
                        ok = m >= 0
                        idx_c = np.nonzero(ok)[0]
                        idx_d = m[ok]
                        if idx_c.size < 8:
                            continue
                        meas = kp_d[idx_d] - kp_c[idx_c]
                        gt = np.stack([gdu[idx_c], gdv[idx_c]], axis=1)
                        err = meas - gt
                        # aşırı uç eşleşmeleri (yanlış eşleşme) ayıkla: |err| > 20 px
                        good = np.linalg.norm(err, axis=1) < 20.0
                        acc["err"].append(err[good])
                        acc["gt_mag"].append(np.linalg.norm(gt[good], axis=1))
                        acc["meas_mag"].append(np.linalg.norm(meas[good], axis=1))
                        acc["ms"].append(float(idx_c.size) / max(len(kp_c), 1))
                        acc["desc"].append(float(np.mean(sc[ok])))

                        # tekrarlanabilirlik: GT ile öngörülen konumda dedektör var mı (<=2 px)
                        pred = kp_c + np.stack([gdu, gdv], axis=1)
                        from scipy.spatial import cKDTree
                        tree = cKDTree(kp_d)
                        d2, _ = tree.query(pred, k=1)
                        acc["rep"].append(float(np.mean(d2 <= 2.0)))

                        # --- ampirik korelasyon: ÖLÇÜLEN yer değiştirme alanından ---
                        mu = meas[good]
                        pc = kp_c[idx_c][good]
                        if len(mu) >= 12:
                            # KARE-BASINA ORTALAMA CIKARILMAZ: temiz ve bozuk kare AYNI
                            # ana aittir, dolayisiyla olculen yer degistirmenin TEK
                            # kaynagi turbulanstir (ego-hareket yok). Ortalamayi cikarmak
                            # alanin ORTAK-MOD (izoplanatik) bilesenini — yani tam da
                            # olcmek istedigimiz uzun menzilli korelasyonu — siler ve
                            # rho_hat'in kuyrugunu yapay olarak bastirirdi.
                            var_acc.append(np.mean(np.sum(mu ** 2, axis=1)) / 2.0)
                            n = len(mu)
                            ia, ib = np.triu_indices(n, k=1)
                            if ia.size > 20000:
                                sel = np.random.default_rng(0).choice(ia.size, 20000,
                                                                      replace=False)
                                ia, ib = ia[sel], ib[sel]
                            d = pc[ia] - pc[ib]
                            r = np.hypot(d[:, 0], d[:, 1])
                            good_r = r > 1e-6
                            ia, ib, d, r = ia[good_r], ib[good_r], d[good_r], r[good_r]
                            ux, uy = d[:, 0] / r, d[:, 1] / r
                            # boyuna / enine bileşenlere ayır
                            la = mu[ia, 0] * ux + mu[ia, 1] * uy
                            lb = mu[ib, 0] * ux + mu[ib, 1] * uy
                            ta = -mu[ia, 0] * uy + mu[ia, 1] * ux
                            tb = -mu[ib, 0] * uy + mu[ib, 1] * ux
                            pair_sep.append(r)
                            pair_prod_l.append(la * lb)
                            pair_prod_t.append(ta * tb)

                    if not acc["err"] or not pair_sep:
                        continue
                    err = np.concatenate(acc["err"])
                    sep = np.concatenate(pair_sep)
                    pl = np.concatenate(pair_prod_l)
                    pt = np.concatenate(pair_prod_t)
                    var = float(np.mean(var_acc))

                    # Kutulanmis ampirik korelasyon. Kutu araligi theta_c ile
                    # OLCEKLENIR: sabit 0-500 px araligi, theta_c = 37 px iken
                    # kuyrugu, theta_c = 750 px iken yalnizca basi orneklerdi ve
                    # kosullar arasi R^2 karsilastirilamaz hale gelirdi.
                    _tc_px = sim.theta_c * p.px_per_rad
                    _rmax = float(min(6.0 * _tc_px, 700.0))
                    edges = np.linspace(0, _rmax, 26)
                    ctr = 0.5 * (edges[1:] + edges[:-1])
                    bi = np.digitize(sep, edges) - 1
                    okb = (bi >= 0) & (bi < len(ctr))
                    cnt = np.bincount(bi[okb], minlength=len(ctr))
                    sl = np.bincount(bi[okb], weights=pl[okb], minlength=len(ctr))
                    st = np.bincount(bi[okb], weights=pt[okb], minlength=len(ctr))
                    with np.errstate(invalid="ignore", divide="ignore"):
                        rho_l = sl / np.maximum(cnt, 1) / var
                        rho_t = st / np.maximum(cnt, 1) / var
                    m_ok = cnt >= 50

                    theta_c_px = sim.theta_c * p.px_per_rad
                    # --- GÜRÜLTÜ SEYRELTMESİ (nugget) DÜZELTMESİ ---
                    # Ölçülen yer değiştirme = türbülans eğimi + İKİ bağımsız dedektör
                    # konum hatası (temiz ve bozuk kare). Dedektör gürültüsü uzaysal
                    # olarak KORELASYONSUZ olduğundan ampirik rho_hat'i sabit bir
                    # A = sigma_tilt^2 / (sigma_tilt^2 + 2 sigma_det^2) çarpanıyla
                    # SEYRELTİR (jeoistatistikte "nugget" etkisi). Model karşılaştırması
                    # bu yüzden GENLİK SERBEST yapılır: rho_hat ~ A * rho_model.
                    # A'nın kendisi de bilgi taşır: sqrt(A * var) türbülansın ölçülen
                    # tilt genliğidir ve teorik sigma ile bağımsız olarak kıyaslanır.
                    # R^2 ve genlik uydurmasi KUTU CIFT SAYISI ile agirliklandirilir:
                    # uzak ayrim kutularinda cift sayisi hizla duser ve agirliksiz
                    # olcut, birkac gurultulu kutunun etkisinde kalir.
                    _w = cnt[m_ok].astype(float)
                    _w = _w / _w.sum() if _w.sum() > 0 else _w

                    def _fit_amp(rho_m):
                        y = rho_l[m_ok]
                        rm = rho_m[m_ok]
                        num = float(np.sum(_w * y * rm))
                        den = float(np.sum(_w * rm ** 2))
                        A = num / den if den > 0 else np.nan
                        pred = A * rm
                        ybar = float(np.sum(_w * y))
                        sst = float(np.sum(_w * (y - ybar) ** 2))
                        r2 = (1.0 - float(np.sum(_w * (y - pred) ** 2)) / sst
                              if sst > 1e-9 else np.nan)
                        return A, r2

                    mod_phys_l = _corr.rho_long(ctr / theta_c_px)
                    mod_phys_t = _corr.rho_tran(ctr / theta_c_px)
                    mod_iso = _corr.rho_iso53(ctr, theta_c_px)
                    A_phys, r2_phys = _fit_amp(mod_phys_l)
                    A_iso, r2_iso = _fit_amp(mod_iso)

                    # serbest (A, theta_c) uydurması — biçim sınaması
                    from scipy.optimize import curve_fit

                    def _fit_free(fun, p0, bounds):
                        try:
                            sig = 1.0 / np.sqrt(np.maximum(cnt[m_ok], 1.0))
                            popt, _ = curve_fit(fun, ctr[m_ok], rho_l[m_ok], p0=p0,
                                                sigma=sig, bounds=bounds, maxfev=20000)
                            pred = fun(ctr[m_ok], *popt)
                            y = rho_l[m_ok]
                            ybar = float(np.sum(_w * y))
                            sst = float(np.sum(_w * (y - ybar) ** 2))
                            return popt, (1.0 - float(np.sum(_w * (y - pred) ** 2)) / sst
                                          if sst > 1e-9 else np.nan)
                        except Exception:
                            return [np.nan] * len(p0), np.nan

                    pf_phys, r2_phys_fit = _fit_free(
                        lambda x, A, tc_: A * _corr.rho_long(x / tc_),
                        [max(A_phys, 1e-3), theta_c_px], ([1e-4, 1.0], [2.0, 5e3]))
                    pf_iso, r2_iso_fit = _fit_free(
                        lambda x, A, tc_: A * np.exp(-(x / tc_) ** (5 / 3)),
                        [max(A_iso, 1e-3), theta_c_px], ([1e-4, 1.0], [2.0, 5e3]))
                    pf_gen, r2_free = _fit_free(
                        lambda x, A, tc_, e: A * np.exp(-(x / tc_) ** e),
                        [max(A_iso, 1e-3), theta_c_px, 1.0],
                        ([1e-4, 1.0, 0.2], [2.0, 5e3, 3.0]))
                    tc53, r2_53fit = float(pf_iso[1]), r2_iso_fit
                    tcf, pf = float(pf_gen[1]), float(pf_gen[2])
                    sigma_tilt_implied = float(np.sqrt(max(A_phys, 0.0) * var))

                    rows.append({
                        "sequence": seq, "D_over_r0": float(dr), "theta_c_over_fov": float(tc),
                        "seed": int(seed), "blur": blur,
                        "sigma_px_theory": float(sim.sigma_px_theory),
                        "theta_c_px": float(theta_c_px),
                        "n_frames": len(acc["err"]),
                        "meas_sigma_px": float(np.sqrt(var)),
                        "pos_err_rms_px": float(np.sqrt(np.mean(np.sum(err ** 2, axis=1)) / 2)),
                        "pos_err_median_px": float(np.median(np.linalg.norm(err, axis=1))),
                        "gt_disp_rms_px": float(np.sqrt(np.mean(
                            np.concatenate(acc["gt_mag"]) ** 2))),
                        "repeatability": float(np.mean(acc["rep"])),
                        "matching_score": float(np.mean(acc["ms"])),
                        "desc_score": float(np.mean(acc["desc"])),
                        "n_kp_clean": float(np.mean(acc["n_clean"])),
                        "n_kp_dist": float(np.mean(acc["n_dist"])),
                        "sep_px": ctr, "rho_long_emp": rho_l, "rho_tran_emp": rho_t,
                        "rho_long_model_phys": mod_phys_l, "rho_tran_model_phys": mod_phys_t,
                        "rho_model_iso53": mod_iso, "bin_counts": cnt,
                        "R2_physical": r2_phys, "R2_iso53": r2_iso,
                        "A_physical": A_phys, "A_iso53": A_iso,
                        "sigma_tilt_implied_px": sigma_tilt_implied,
                        "R2_physical_fitted": r2_phys_fit,
                        "theta_c_px_physical_fit": float(pf_phys[1]),
                        "R2_iso53_fitted": r2_53fit, "theta_c_px_iso53_fit": tc53,
                        "R2_free": r2_free, "free_exponent": pf, "theta_c_px_free": tcf,
                        "elapsed_s": time.time() - t0,
                    })
                    print(f"[E2] {seq} dr={dr} tc={tc} s={seed}: "
                          f"sig_ölç={rows[-1]['meas_sigma_px']:.2f} "
                          f"sig_tilt={sigma_tilt_implied:.2f} (teorik {sim.sigma_px_theory:.2f}) "
                          f"A={A_phys:.2f} hata={rows[-1]['pos_err_rms_px']:.2f} "
                          f"tekrar={rows[-1]['repeatability']:.2f} "
                          f"R2fiz={r2_phys:.3f} R2iso={r2_iso:.3f} "
                          f"R2fizfit={r2_phys_fit:.3f} R2isofit={r2_iso_fit:.3f} "
                          f"üs={pf:.2f}(R2={r2_free:.3f}) ({rows[-1]['elapsed_s']:.0f}s)",
                          flush=True)
        del matcher
        torch.cuda.empty_cache()

    with open(os.path.join(outdir, "e2_results.json"), "w") as f:
        json.dump(_json_safe(rows), f, indent=1)
    print(f"\n✓ {os.path.join(outdir, 'e2_results.json')} ({len(rows)} hücre)")
    return rows


def main():
    ap = argparse.ArgumentParser(description="E2 ön-uç karakterizasyonu")
    ap.add_argument("--sequences", nargs="+", default=["MH01", "MH03", "MH05"])
    ap.add_argument("--base_path", default="/data/sdata/")
    ap.add_argument("--D_over_r0", nargs="+", type=float, default=[1, 2, 3, 5, 8, 12])
    ap.add_argument("--theta_c", nargs="+", type=float, default=[0.05, 0.15, 0.40, 1.0])
    ap.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    ap.add_argument("--n_frames", type=int, default=30)
    ap.add_argument("--frame_stride", type=int, default=20)
    ap.add_argument("--blur", default="speckle_global")
    ap.add_argument("--max_kp", type=int, default=1024)
    ap.add_argument("--outdir", default="evidence/E2")
    a = ap.parse_args()
    run(a.sequences, a.base_path, a.D_over_r0, a.theta_c, a.seeds, a.n_frames,
        a.frame_stride, a.blur, a.outdir, a.max_kp)


if __name__ == "__main__":
    main()
