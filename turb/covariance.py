"""
İzoplanatik yama gruplaması ve türbülans-farkında ölçüm kovaryansı (C3).

Ölçüm modeli (rapor §4.1):
    z_ij = pi(T_i, l_j) + d(theta_j, t_i) + n_ij
d, sıfır ortalamalı ama UZAYSAL OLARAK KORELASYONLU türbülans eğim alanıdır.
Yerleşik VIO/SLAM arka-uçları Cov(d) = 0 ve n ~ N(0, sigma^2 I) varsayar; bu ihlal
kestiriciyi aşırı-güvene ve tutarsızlığa sokar.

Bu modül, bir karedeki gözlemleri izoplanatik yamalara ayırır ve her yama için
2N x 2N (mono) / 3N x 3N (stereo dahil) ölçüm kovaryansını üretir:

    Sigma = Sigma_taban  +  sigma_px^2 * C_turb

  Sigma_taban : hattın mevcut gürültü modeli (mono px_noise, stereo stereo_sigma)
  C_turb      : turb.correlation'dan gelen korelasyon yapısı
                - "physical": tam Kolmogorov tensörü (boyuna != enine, C_xy != 0)
                - "iso53"   : rapor §4.2 mühendislik modeli  rho (x) I_2

STEREO GÖZLEMLER: (uL, uR, v) üçlüsünde uL ve v cam0'ın, uR ise cam1'in eğim
alanından etkilenir. İki kameranın ayak izleri stereo taban çizgisi B kadar
ayrıldığından kameralar arası korelasyon rho(B/D)'dir — bu, modelin doğal bir
parçasıdır ve B >> r0 iken uR neredeyse bağımsız bir bozulma taşır.
"""
from __future__ import annotations

import numpy as np

from . import correlation as _corr

# E2 ile ölçülen ampirik korelasyon parametreleri (bkz. turb_correlation_blocks).
# Medyanlar: genlik 0.68-0.83 -> 0.75 ; üs 0.32-0.74 (D/r0=3-8 aralığında ~0.40).
EMP_A = 0.75
EMP_P = 0.40


# =============================================================================
def group_isoplanatic(uv: np.ndarray, theta_c_px: float, max_per_patch: int = 12,
                      grid_scale: float = 1.0):
    """Gözlemleri izoplanatik yamalara ayır.

    Yamalar, korelasyon açısı theta_c'nin piksel karşılığında bir ızgara ile
    tanımlanır (basit, belirlenimci ve O(N)). Aynı hücredeki öznitelikler güçlü
    korelasyonludur; hücre sınırındaki zayıf korelasyon ihmal edilir — bu, tam
    kovaryansın BLOK-KÖŞEGEN yaklaşımıdır ve kestirimi tutucu (conservative) kılar.

    Args:
        uv: (N,2) piksel konumları
        theta_c_px: korelasyon açısının piksel karşılığı
        max_per_patch: bir yamadaki en çok öznitelik (matris boyutunu sınırlar)
    Dönüş: yama listesi — her biri gözlem indeksleri dizisi
    """
    if len(uv) == 0:
        return []
    cell = max(theta_c_px * grid_scale, 8.0)
    keys = np.floor(uv / cell).astype(np.int64)
    order = np.lexsort((keys[:, 0], keys[:, 1]))
    patches = []
    start = 0
    ks = keys[order]
    for i in range(1, len(order) + 1):
        if i == len(order) or not np.array_equal(ks[i], ks[start]):
            idx = order[start:i]
            for c in range(0, len(idx), max_per_patch):
                patches.append(idx[c:c + max_per_patch])
            start = i
    return patches


def select_patch_aware(uv: np.ndarray, scores: np.ndarray | None, theta_c_px: float,
                       max_per_patch: int = 5, min_total_frac: float = 0.70,
                       sigma_px: float = 0.0, sigma_min_px: float = 0.25):
    """A4: izoplanatik-yama farkındalıklı öznitelik seçimi.

    Aynı izoplanatik yamadaki öznitelikler neredeyse AYNI bozulmayı yaşar; bir
    yamadan çok sayıda gözlem almak bağımsız bilgi eklemeden kestiriciyi o yamanın
    ortak hatasına bağlar. Yama başına en iyi K gözlem tutulur.

    İKİ EMNİYET KAPISI (ölçümle eklendi — bunlarsız kol ıraksıyordu):

    1. **Gözlem tabanı.** Seçim, kare başına gözlem sayısını `min_total_frac`
       oranının altına düşüremez. Aksi halde graf yeterince kısıtlanmaz ve
       BatchFixedLagSmoother indeterminant sisteme girer (ölçüldü: K=5'te
       8 hücrenin 5'i "IRAKSAK config" ile düştü). K, oran sağlanana kadar artırılır.

    2. **Türbülans eşiği.** Bozulma yokken (sigma_px < sigma_min_px) seçim
       UYGULANMAZ. Türbülans yokken yama kavramı anlamsızdır; gözlem atmak yalnız
       bilgi kaybıdır (ölçüldü: temiz veride A4 ATE = 4.70 m, A0 = 0.11 m).

    Dönüş: tutulacak gözlem indeksleri (bool maske).
    """
    n = len(uv)
    if n == 0:
        return np.zeros(0, dtype=bool)
    if sigma_px < sigma_min_px:
        return np.ones(n, dtype=bool)          # kapı 2: bozulma yok -> seçim yok

    patches = group_isoplanatic(uv, theta_c_px, max_per_patch=10 ** 6)
    if scores is not None:
        sc = np.asarray(scores)
        patches = [idx[np.argsort(-sc[idx])] for idx in patches]

    need = int(np.ceil(min_total_frac * n))
    K = int(max_per_patch)
    while True:                                 # kapı 1: gözlem tabanı
        keep = np.zeros(n, dtype=bool)
        for idx in patches:
            keep[idx[:K]] = True
        if keep.sum() >= need or K >= n:
            return keep
        K += 1


# =============================================================================
def turb_correlation_blocks(theta_x: np.ndarray, theta_y: np.ndarray,
                            z_layer: float, D: float, model: str,
                            theta_c: float | None = None,
                            baseline_m: float = 0.0):
    """Öznitelikler arası 2x2 korelasyon blokları (N,N,2,2).

    baseline_m: iki ayak izi arasına EKLENEN doğrusal ötelenme [m] (stereo taban
    çizgisi). Sıfırdan farklıysa bloklar cam0<->cam1 ayrımını temsil eder.
    TÜM modeller bu ötelenmeyi görmek ZORUNDADIR: ampirik/iso53 modeller açısal
    ayrımla parametrelidir, bu yüzden ötelenme z_layer ile açıya çevrilir
    (dth += baseline_m / z_layer). Aksi halde aynı özniteliğin uL ve uR'si TAM
    KORELASYONLU sayılır, disparitenin türbülans belirsizliği SIFIR olur ve
    derinlik aşırı güvenli kestirilir (ölçüldü: A5 NEES 194 vs skaler A1 34).
    """
    tx = np.asarray(theta_x, dtype=np.float64)
    ty = np.asarray(theta_y, dtype=np.float64)
    dx = z_layer * (tx[:, None] - tx[None, :]) + float(baseline_m)
    dy = z_layer * (ty[:, None] - ty[None, :])
    if model == "empirical":
        # ÖN-UÇ SEVİYESİNDE KALİBRE EDİLMİŞ KORELASYON (E2 ile ölçüldü).
        #
        # Açıklık-ortalamalı tilt tensörü FAZIN doğru fiziğidir, ama ön-uç fazı
        # değil EŞLEŞTİRİLMİŞ ÖZNİTELİK KAYMASINI ölçer. E2 ikisinin AYNI OLMADIĞINI
        # gösterdi (ampirik rho_hat, dairesel olmayan biçimde ALIKED+LightGlue'nun
        # gerçek kaymalarından hesaplanır):
        #     D/r0 :   3      5      8     12
        #     R^2 fizik tensörü : 0.08  0.51  0.46  0.68
        #     R^2 serbest üs    : 0.86  0.92  0.94  0.94
        #     uydurulan üs      : 0.32  0.41  0.42  0.74   (tilt tensörü çok daha dik)
        #     genlik A          : 0.68  0.83  0.71  0.67
        # Yani gözlenen korelasyon hem DAHA ZAYIF (A<1) hem DAHA YAVAŞ SÖNÜYOR.
        # Nedeni fizik hatası değil ÖLÇME ZİNCİRİ: dedektör sonlu bir pencerede
        # ağırlık merkezi bulur, kısa-poz PSF her özniteliği bağımsız kaydırır ve
        # eşleştirme kendi hatasını ekler. Bu terimleri modele koymadan tilt
        # tensörünü doğrudan ölçüm kovaryansı saymak, DİFERANSİYEL yönleri aşırı
        # güvenli yapar (ölçüldü: A5 NEES 194 vs skaler A1 34).
        if theta_c is None:
            raise ValueError("empirical theta_c ister")
        dth = np.hypot(dx, dy) / max(z_layer, 1e-12)   # taban çizgisi DAHİL açısal ayrım
        R = EMP_A * np.exp(-(np.maximum(dth, 0.0) / max(theta_c, 1e-12)) ** EMP_P)
        if baseline_m == 0.0:
            np.fill_diagonal(R, 1.0)      # aynı öznitelik, aynı kamera: tam korelasyon
        Cxx = R; Cyy = R; Cxy = np.zeros_like(R)
    elif model == "iso53":
        if theta_c is None:
            raise ValueError("iso53 theta_c ister")
        R = _corr.rho_iso53(np.hypot(dx, dy) / max(z_layer, 1e-12), theta_c)
        Cxx = R; Cyy = R; Cxy = np.zeros_like(R)
    else:
        Cxx, Cxy, Cyy = _corr.tensor_blocks(dx, dy, D)
    return Cxx, Cxy, Cyy


def patch_covariance(obs_theta: np.ndarray, is_stereo: np.ndarray,
                     sigma_px: float, base_sigma_mono: float, base_sigma_stereo: float,
                     z_layer: float, D: float, stereo_baseline: float,
                     model: str = "physical", theta_c: float | None = None,
                     jitter: float = 1e-9, shrinkage: float = 0.05,
                     sigma_ind_px: float = 0.0) -> np.ndarray:
    """Bir yama için tam ölçüm kovaryansı.

    Gözlem sırası: her gözlem için stereo ise (uL, uR, v), mono ise (u, v).
    Args:
        obs_theta: (N,2) açısal konumlar [rad]
        is_stereo: (N,) bool
    Dönüş: (M, M) simetrik pozitif tanımlı kovaryans; M = sum(3 if stereo else 2)
    """
    tx = obs_theta[:, 0]
    ty = obs_theta[:, 1]
    N = len(tx)
    Cxx, Cxy, Cyy = turb_correlation_blocks(tx, ty, z_layer, D, model, theta_c)

    # Kameralar arası korelasyon: ayak izleri stereo taban çizgisi kadar ötelenir.
    Sxx, Sxy, _ = turb_correlation_blocks(tx, ty, z_layer, D, model, theta_c,
                                          baseline_m=stereo_baseline)

    # blok yerleşimi
    offs = np.zeros(N, dtype=np.int64)
    dim = 0
    for i in range(N):
        offs[i] = dim
        dim += 3 if is_stereo[i] else 2
    S = np.zeros((dim, dim), dtype=np.float64)
    s2 = sigma_px ** 2

    for i in range(N):
        oi = offs[i]
        si = is_stereo[i]
        for j in range(N):
            oj = offs[j]
            sj = is_stereo[j]
            # cam0-cam0 (u,v) blokları
            uu = s2 * Cxx[i, j]
            uv_ = s2 * Cxy[i, j]
            vv = s2 * Cyy[i, j]
            iu, iv = oi, oi + (2 if si else 1)
            ju, jv = oj, oj + (2 if sj else 1)
            S[iu, ju] += uu
            S[iu, jv] += uv_
            S[iv, ju] += uv_
            S[iv, jv] += vv
            if si:
                ir = oi + 1                     # uR (cam1)
                # uR(i) - u(j): cam1 <-> cam0 ayrımı
                S[ir, ju] += s2 * Sxx[i, j]
                S[ir, jv] += s2 * Sxy[i, j]
                S[iu, oj + 1] += s2 * Sxx[j, i] if sj else 0.0
                if sj:
                    S[ir, oj + 1] += s2 * Cxx[i, j]   # cam1 <-> cam1 (aynı ayrım)
                    S[iv, oj + 1] += s2 * Sxy[j, i]
    S = 0.5 * (S + S.T)

    # KÜÇÜLTME (shrinkage): korelasyon yapısı 1'e çok yaklaştığında ortak-mod yönü
    # neredeyse kısıtsız kalır ve doğrusal sistem belirsizleşir (ampirik olarak
    # gözlendi: A5 kolunda IndeterminantLinearSystem). Sigma_turb <- (1-l) Sigma_turb
    # + l diag(Sigma_turb) düzenlemesi modeli TUTUCU yönde bozar (korelasyonu bir miktar
    # azaltır) ve iyi koşullu tutar. l=0 tam model, l=1 köşegen.
    if shrinkage > 0:
        d = np.diag(np.diag(S))
        S = (1.0 - shrinkage) * S + shrinkage * d

    # taban gürültü + TÜRBÜLANSIN BAĞIMSIZ BİLEŞENİ
    #
    # KRİTİK (2026-08-19'da ölçümle bulundu): türbülans öznitelik konum hatasına İKİ
    # bileşen katar ve model uzun süre yalnız birincisini taşıyordu:
    #   (1) KORELASYONLU eğim (tilt) alanı -> sigma_px^2 * C   [zaten vardı]
    #   (2) BAĞIMSIZ bileşen: kısa-poz PSF'in her özniteliğin ağırlık merkezini
    #       birbirinden BAĞIMSIZ kaydırması + eşleştirme hatasının bozulmayla
    #       büyümesi -> sigma_ind^2 * I
    # E2 bunu ZATEN ÖLÇMÜŞTÜ (`pos_err_rms_px` = ölçülen kaymanın GT eğim alanıyla
    # AÇIKLANAMAYAN kısmı) ama kovaryansa hiç beslenmiyordu:
    #     D/r0=3 : tilt 0.967 px, bağımsız 1.406 px  (oran 1.46)
    #     D/r0=5 : tilt 1.480 px, bağımsız 2.596 px  (oran 1.75)
    #     D/r0=8 : tilt 2.189 px, bağımsız 3.893 px  (oran 1.78)
    #     D/r0=12: tilt 3.069 px, bağımsız 5.573 px  (oran 1.82)
    # Yani bağımsız bileşen korelasyonludan BÜYÜKTÜR, ama model onun yerine temiz-veri
    # taban gürültüsünü (0.5 px) kullanıyordu — D/r0=8'de ~8 kat eksik.
    # SONUCU: Sigma'nın DİFERANSİYEL yönleri aşırı-güvenli oluyor, yama faktörü
    # öznitelikler arası farkları olduğundan çok daha kesin sanıyor ve NEES patlıyor.
    # Ölçüldü: aynı ölçekle skaler şişirme (A1) NEES 32 verirken yapısal kovaryans
    # (A5) 254 veriyordu — yapı suçlu değildi, EKSİK BAĞIMSIZ TERİM suçluydu.
    ind2 = float(sigma_ind_px) ** 2
    for i in range(N):
        oi = offs[i]
        if is_stereo[i]:
            S[oi, oi] += base_sigma_stereo ** 2 + ind2
            S[oi + 1, oi + 1] += base_sigma_stereo ** 2 + ind2
            S[oi + 2, oi + 2] += base_sigma_stereo ** 2 + ind2
        else:
            S[oi, oi] += base_sigma_mono ** 2 + ind2
            S[oi + 1, oi + 1] += base_sigma_mono ** 2 + ind2
    S[np.diag_indices_from(S)] += jitter
    return S


def scalar_inflation(base_sigma: float, sigma_px: float) -> float:
    """A1: skaler kovaryans şişirme — sqrt(sigma_taban^2 + sigma_px^2)."""
    return float(np.sqrt(base_sigma ** 2 + sigma_px ** 2))


# =============================================================================
def sigma_px_from_r0(D_over_r0: float, lam: float, D: float, px_per_rad: float) -> float:
    """Teorik tek eksen tilt jitteri [px] (A1/A3/A5 için ölçek parametresi)."""
    from .optics import sigma_tilt_px
    return sigma_tilt_px(D_over_r0, lam, D, px_per_rad)
