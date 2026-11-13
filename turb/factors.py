"""
Yama-gruplu ölçüm faktörü (rapor §4.4'teki PatchGroupedProjectionFactor).

Tek bir kareye ait, AYNI izoplanatik yamadaki N gözlem tek bir faktörde toplanır ve
tam (köşegen-olmayan) kovaryansla verilir:

    graph.add( PatchGroupedProjectionFactor(
                   patch.measurements,
                   noiseModel.Gaussian.Covariance(Sigma_patch),   # M x M
                   X(i), patch.landmarkKeys ) )

GTSAM Python'da hazır böyle bir faktör olmadığından `gtsam.CustomFactor` üzerine
kurulur. `noiseModel.Gaussian.Covariance` karekök bilgi matrisini kendi hesaplar ve
beyazlatmayı otomatik yapar — filtre tabanlı sistemlerde (MSCKF) gereken elle
Cholesky beyazlatma adımı burada gerekmez.

PERFORMANS: CustomFactor'ın hata/Jakobiyen fonksiyonu her yeniden-doğrusallaştırmada
Python'a döner. Bu yüzden yama içi TÜM gözlemler tek çağrıda VEKTÖRİZE hesaplanır ve
Jakobiyenler analitik olarak yazılır (gözlem başına gtsam nesnesi kurulmaz).
Doğruluk, `verify_jacobians()` ile GTSAM'in kendi faktörlerine karşı sınanır.

GTSAM konvansiyonu (doğrulanmıştır):
    X(f) = T_WC (kameranın dünyadaki pozu), retraction T <- T * Exp(xi),
    xi = [rot(3), trans(3)]  ->  d(p_c)/d(xi) = [ [p_c]_x , -I ]
"""
from __future__ import annotations

import numpy as np
import gtsam
from gtsam.symbol_shorthand import X, L


def _skew(v):
    return np.array([[0.0, -v[2], v[1]],
                     [v[2], 0.0, -v[0]],
                     [-v[1], v[0], 0.0]])


def make_patch_factor(frame_key: int, lm_keys, measurements: np.ndarray,
                      is_stereo: np.ndarray, Sigma: np.ndarray,
                      fx: float, fy: float, cx: float, cy: float, baseline: float,
                      robust: bool = False, huber_k: float = 1.345,
                      huber_k_mono: float | None = None,
                      huber_k_stereo: float | None = None):
    """Yama-gruplu faktör üret.

    Args:
        frame_key: X(f) anahtarı
        lm_keys: [L(j)] anahtar listesi (N adet)
        measurements: (M,) ölçüm vektörü; stereo gözlem (uL,uR,v), mono (u,v)
        is_stereo: (N,) bool
        Sigma: (M,M) tam ölçüm kovaryansı
        robust: Huber m-kestiricisi ile sar. VARSAYILAN False'tur ama ARKA-UÇ
            BUNU AÇIKÇA GEÇMEK ZORUNDADIR (bkz. aşağıdaki not).

    ROBUST ÇEKİRDEK ASİMETRİSİ (2026-08-21'de ölçülerek bulundu — KRİTİK).
    Arka-uçtaki bireysel-faktör yolları (A0/A0C/A1/A2) gürültü modellerini
    `--turb_robust` bayrağına göre Huber ile sarıyordu; bu yama yolu ise
    SARMIYORDU. Yani skaler şişirme kollarında bir aykırı-değer çekirdeği
    varken yapılı kollarda yoktu ve A1 ile A5 "yalnız kovaryans yapısında
    farklı" DEĞİLDİ. m-kestiricisi büyük artıkları bastırarak kestiricinin
    ETKİN kovaryansını büyütür; NEES tam da beyan edilen ile etkin kovaryansı
    karşılaştırdığı için bu, C3'ün ana karşılaştırmasını doğrudan kirletir.
    Ölçülen etki: A1 NEES 48, A5 NEES 193 (yapı 3.8 kat "kötü" görünüyordu),
    oysa kapalı-form türetim A5'in 1.2-1.5 kat İYİ olmasını gerektiriyor.

    GRANÜLERLİK FARKI (makalede beyan edilir): bireysel faktörde Huber gözlem
    başına, yama faktöründe BLOK başına uygulanır (GTSAM m-kestiriciyi
    beyazlatılmış artığın normuna uygular). Bunlar birebir aynı değildir;
    tam eşleşmiş karşılaştırma için `--turb_robust 0` ablasyonu koşulur.

    GÖZLEM BAŞINA HUBER (pyeco_turb, 2026-11-07; `huber_k_mono`/`huber_k_stereo`).
    Yukarıdaki granülerlik farkını giderir: her gözlem i için kendi marjinal
    bloğuyla beyazlatılmış norm x_i = sqrt(e_i^T Sigma_ii^-1 e_i) hesaplanır ve
    GTSAM'in Huber ağırlığı w_i = 1 (x_i <= k) | k/x_i (x_i > k) alınır; mono için
    k = huber_k_mono (2 boyutlu norm), stereo için k = huber_k_stereo (3 boyutlu
    norm) — yani bireysel GenericProjection/GenericStereo faktörlerindeki
    Robust(Huber) ile AYNI eşik ve AYNI granülerlik. Ağırlık, gözlemin artığını ve
    Jakobiyen satırlarını sqrt(w_i) ile ölçekleyerek uygulanır; bu, kovaryansı
    Sigma_w = W^-1/2 Sigma W^-1/2 ile değiştirmeye denktir (korelasyon yapısı
    korunur). Maliyet değerlendirmesinde (Jakobiyen istenmeyen çağrı) artık
    sqrt(2 rho(x_i))/x_i ile ölçeklenir. Köşegen Sigma'da hem doğrusallaştırılmış
    sistem hem maliyet GTSAM Robust(Huber) ile birebir aynıdır; korelasyonlu yamada
    eşiği aşan gözlemler için maliyet, gözlemler arası çapraz terimler yüzünden
    yaklaşıktır (eşik altındaki gözlemlerde tam Gauss).
    `robust=True` (blok başına Huber) ile birlikte kullanılmaz.
    """
    n = len(lm_keys)
    is_stereo = np.asarray(is_stereo, dtype=bool)
    z = np.asarray(measurements, dtype=np.float64).ravel()
    offs = np.zeros(n, dtype=np.int64)
    dim = 0
    for i in range(n):
        offs[i] = dim
        dim += 3 if is_stereo[i] else 2
    assert z.size == dim and Sigma.shape == (dim, dim)

    per_obs = huber_k_mono is not None or huber_k_stereo is not None
    if per_obs and robust:
        raise ValueError("blok-başına ve gözlem-başına Huber birlikte kullanılamaz")
    if per_obs:
        _km = float(huber_k_mono) if huber_k_mono else np.inf
        _ks = float(huber_k_stereo) if huber_k_stereo else np.inf
        _kk = np.array([_ks if is_stereo[i] else _km for i in range(n)])
        _Sinv = [np.linalg.inv(Sigma[offs[i]:offs[i] + (3 if is_stereo[i] else 2),
                                     offs[i]:offs[i] + (3 if is_stereo[i] else 2)])
                 for i in range(n)]
        _sl = [slice(int(offs[i]), int(offs[i]) + (3 if is_stereo[i] else 2)) for i in range(n)]
        # blok-köşegen marjinal bilgi matrisi ve bileşen -> gözlem eşlemi (vektörize norm için)
        _Sblk = np.zeros((dim, dim))
        _obs_of = np.empty(dim, dtype=np.int64)
        for i in range(n):
            _Sblk[_sl[i], _sl[i]] = _Sinv[i]
            _obs_of[_sl[i]] = i

    model = gtsam.noiseModel.Gaussian.Covariance(Sigma)
    if robust:
        model = gtsam.noiseModel.Robust.Create(
            gtsam.noiseModel.mEstimator.Huber.Create(float(huber_k)), model)
    keys = gtsam.KeyVector([frame_key] + list(lm_keys))
    st_idx = np.nonzero(is_stereo)[0]
    mo_idx = np.nonzero(~is_stereo)[0]

    def error_func(this, values, H):
        P = values.atPose3(gtsam.Symbol(frame_key).key() if False else frame_key)
        R = P.rotation().matrix()
        t = P.translation()
        pw = np.empty((n, 3))
        for i, k in enumerate(lm_keys):
            pw[i] = values.atPoint3(k)
        pc = (pw - t) @ R                       # = R^T (pw - t)  (satır bazlı)
        zc = np.maximum(pc[:, 2], 1e-6)
        u_pred = fx * pc[:, 0] / zc + cx
        v_pred = fy * pc[:, 1] / zc + cy
        uR_pred = fx * (pc[:, 0] - baseline) / zc + cx

        e = np.empty(dim)
        e[offs[mo_idx]] = u_pred[mo_idx] - z[offs[mo_idx]]
        e[offs[mo_idx] + 1] = v_pred[mo_idx] - z[offs[mo_idx] + 1]
        if st_idx.size:
            e[offs[st_idx]] = u_pred[st_idx] - z[offs[st_idx]]
            e[offs[st_idx] + 1] = uR_pred[st_idx] - z[offs[st_idx] + 1]
            e[offs[st_idx] + 2] = v_pred[st_idx] - z[offs[st_idx] + 2]

        if H is not None:
            # VEKTÖRİZE Jakobiyen: yama içindeki tüm gözlemler tek seferde.
            # d(pi)/d(p_c) satırları; d(p_c)/d(xi) = [ [p_c]_x , -I ] ve a^T [p]_x = (a x p)^T.
            x_, y_ = pc[:, 0], pc[:, 1]
            iz = 1.0 / zc
            iz2 = iz * iz
            Ju = np.zeros((n, 3)); Jv = np.zeros((n, 3)); Jr = np.zeros((n, 3))
            Ju[:, 0] = fx * iz; Ju[:, 2] = -fx * x_ * iz2
            Jv[:, 1] = fy * iz; Jv[:, 2] = -fy * y_ * iz2
            Jr[:, 0] = fx * iz; Jr[:, 2] = -fx * (x_ - baseline) * iz2

            def _P(J):
                # [J x p_c , -J]  (np.cross'tan hızlı, aynı sonuç)
                out = np.empty((n, 6))
                out[:, 0] = J[:, 1] * pc[:, 2] - J[:, 2] * pc[:, 1]
                out[:, 1] = J[:, 2] * pc[:, 0] - J[:, 0] * pc[:, 2]
                out[:, 2] = J[:, 0] * pc[:, 1] - J[:, 1] * pc[:, 0]
                out[:, 3:] = -J
                return out
            Pu, Pv, Pr = _P(Ju), _P(Jv), _P(Jr)
            Lu = Ju @ R.T                              # (n,3)
            Lv = Jv @ R.T
            Lr = Jr @ R.T

            Hp = np.empty((dim, 6))
            HL = np.zeros((n, dim, 3))
            if st_idx.size:
                o = offs[st_idx]
                Hp[o] = Pu[st_idx]; Hp[o + 1] = Pr[st_idx]; Hp[o + 2] = Pv[st_idx]
                HL[st_idx, o] = Lu[st_idx]; HL[st_idx, o + 1] = Lr[st_idx]; HL[st_idx, o + 2] = Lv[st_idx]
            if mo_idx.size:
                o = offs[mo_idx]
                Hp[o] = Pu[mo_idx]; Hp[o + 1] = Pv[mo_idx]
                HL[mo_idx, o] = Lu[mo_idx]; HL[mo_idx, o + 1] = Lv[mo_idx]
            H[0] = Hp
            for i in range(n):
                H[i + 1] = HL[i]
        if per_obs:
            # gözlem başına Huber. Doğrusallaştırmada (H verilmiş) IRLS ağırlığı sqrt(w_i)
            # (GTSAM Robust::WhitenSystem ile aynı); yalnız maliyet istendiğinde (H None)
            # ölçek sqrt(2 rho(x_i))/x_i, yani 0.5 x'^2 = rho(x_i) (GTSAM Robust::loss ile aynı).
            xs = np.sqrt(np.maximum(np.bincount(_obs_of, weights=e * (_Sblk @ e), minlength=n), 0.0))
            sw = np.ones(dim)
            big = xs > _kk
            if np.any(big):
                if H is not None:
                    f_ = np.sqrt(_kk[big] / xs[big])
                else:
                    f_ = np.sqrt(2.0 * (_kk[big] * xs[big] - 0.5 * _kk[big] ** 2)) / xs[big]
                fo = np.ones(n); fo[big] = f_
                sw = fo[_obs_of]
            if np.any(sw != 1.0):
                e = e * sw
                if H is not None:
                    for j in range(len(H)):
                        H[j] = H[j] * sw[:, None]
        return e

    return gtsam.CustomFactor(model, keys, error_func), error_func


# =============================================================================
def verify_jacobians(seed: int = 0, tol: float = 1e-5) -> dict:
    """Analitik Jakobiyenleri GTSAM'in kendi faktörlerine karşı doğrula.

    Mono için GenericProjectionFactorCal3_S2, stereo için GenericStereoFactor3D
    referans alınır; ayrıca sayısal türevle de karşılaştırılır.
    """
    rng = np.random.default_rng(seed)
    fx, fy, cx, cy, b = 458.0, 457.0, 367.0, 248.0, 0.11
    K = gtsam.Cal3_S2(fx, fy, 0.0, cx, cy)
    Ks = gtsam.Cal3_S2Stereo(fx, fy, 0.0, cx, cy, b)
    pose = gtsam.Pose3(gtsam.Rot3.Expmap(rng.normal(0, 0.3, 3)),
                       gtsam.Point3(*rng.normal(0, 1.0, 3)))
    # Noktalar kamera önünde üretilir (cheirality garantisi)
    pts = [np.asarray(pose.transformFrom(
        gtsam.Point3(*(rng.normal(0, 1.5, 2).tolist() + [rng.uniform(5.0, 12.0)]))),
        dtype=float) for _ in range(3)]
    def _mkvals(pose_, pts_):
        v = gtsam.Values()
        v.insert(X(0), pose_)
        for i, p in enumerate(pts_):
            v.insert(L(i), p)
        return v

    vals = _mkvals(pose, pts)

    is_st = np.array([True, False, True])
    meas = []
    for i, p in enumerate(pts):
        cam = gtsam.PinholeCameraCal3_S2(pose, K)
        uv = cam.project(p)
        if is_st[i]:
            scam = gtsam.StereoCamera(pose, Ks)
            sp = scam.project(p)
            meas += [sp.uL(), sp.uR(), sp.v()]
        else:
            meas += [uv[0], uv[1]]
    meas = np.array(meas) + rng.normal(0, 0.3, len(meas))
    dim = len(meas)
    Sigma = np.eye(dim) * 0.25
    f, efun = make_patch_factor(X(0), [L(0), L(1), L(2)], meas, is_st, Sigma,
                                fx, fy, cx, cy, b)
    f.linearize(vals)   # doğrusallaştırmanın hata vermediğini de sınar

    def raw_error(values):
        return np.asarray(efun(None, values, None))

    e0 = raw_error(vals)
    eps = 1e-6
    num_pose = np.zeros((dim, 6))
    for k in range(6):
        d = np.zeros(6); d[k] = eps
        v2 = _mkvals(pose.retract(d), pts)
        num_pose[:, k] = (raw_error(v2) - e0) / eps
    Ha = [np.zeros((dim, 6))] + [np.zeros((dim, 3)) for _ in range(3)]
    efun(None, vals, Ha)
    err_pose = float(np.max(np.abs(num_pose - Ha[0])))

    err_pt = 0.0
    for li in range(3):
        num_pt = np.zeros((dim, 3))
        for k in range(3):
            d = np.zeros(3); d[k] = eps
            pts2 = [np.array(q, dtype=float) for q in pts]
            pts2[li] = pts2[li] + d
            v2 = _mkvals(pose, pts2)
            num_pt[:, k] = (raw_error(v2) - e0) / eps
        err_pt = max(err_pt, float(np.max(np.abs(num_pt - Ha[li + 1]))))

    # GTSAM referans faktörleriyle hata değeri karşılaştırması
    ref = gtsam.GenericStereoFactor3D(gtsam.StereoPoint2(meas[0], meas[1], meas[2]),
                                      gtsam.noiseModel.Isotropic.Sigma(3, 0.5),
                                      X(0), L(0), Ks)
    ref_e = ref.unwhitenedError(vals)
    err_val = float(np.max(np.abs(np.asarray(ref_e) - e0[:3])))

    return {"jac_pose_max_err": err_pose, "jac_point_max_err": err_pt,
            "error_vs_gtsam_stereo": err_val,
            "ok": bool(err_pose < 1e-3 and err_pt < 1e-3 and err_val < 1e-9)}
