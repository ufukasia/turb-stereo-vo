"""
Deney tanımları — E0/E1/E4/E5 hücre kümeleri (rapor §6).

    python -m turb.experiments --exp E1a --jobs 12 --gpu_jobs 3
    python -m turb.experiments --exp E4  --jobs 10 --gpu_jobs 2

Hücre kümeleri rapor §5.3'teki denetim parametrelerini birebir izler. Temiz
(D/r0 = 0) hücreler TEK tohumla koşulur: hat belirlenimcidir, türbülans yoksa
gerçekleme varyansı da yoktur (rapor §6.0). Bu, eşleştirilmiş tasarımın referans
kolunu bedava verir.
"""
from __future__ import annotations

import argparse
import os

from .runner import Cell, run_experiment, run_experiment_pipelined

EUROC = "/data/sdata/"
TUMVI = "/data/sdata_tumvi_asl/"
TARTAN = "/data/sdata_tartanair/"

EUROC_SEQS = ["MH01", "MH02", "MH03", "MH04", "MH05"]
TUMVI_SEQS = ["slides1", "slides2", "slides3"]
# TartanAir: uzun yörüngeler (>= 600 kare) — pencere protokolü buna göre kısaltılır
TARTAN_SEQS = ["westerndesert_Easy_P004", "seasonsforest_Easy_P004"]

# Rapor §5.3
DR_LEVELS = [1, 2, 3, 5, 8, 12]
TAU_LEVELS = [0.1, 0.5, 1.0, 2.0, 5.0]
TC_LEVELS = [0.05, 0.15, 0.40, 1.0]

# TartanAir profili: 24 m/s'ye varan hareket ve uzak sahne -> EuRoC için Tabu ile
# ayarlanmış derinlik/sıçrama kapıları geçersiz. Bu ayarlar ölçümle seçildi
# (bkz. evidence/tartanair_profile.md).
TARTAN_CFG = (("frontend", "max_landmark_depth", 60.0),
              ("frontend", "replenish_threshold", 90),
              ("frontend", "target_feature_count", 250),
              ("frontend", "max_keypoints", 600))
TARTAN_ARGV = ("--max_init_depth", "60.0", "--disparity_bias", "0.0",
               "--jump_max_trans", "0.0", "--pnp_min_inliers", "8")

# --- KITTI (arac, dis mekan, gercek kamera) ---------------------------------
# TRUBA'da kosar; yerelde base_path yoksa hucre uretilmez.
KITTI = "/arf/scratch/uasil/data/sdata_kitti/"
# Pencere protokolu window_len=300 + son pencere 700 -> >=1000 kare gerekir.
# kitti04 (271) ve kitti03 (801) DISLANDI: kisa sekansta pencere protokolu
# dejenere olur (bkz. sanity I1-dejenere-pencere).
KITTI_SEQS = ["kitti00", "kitti02", "kitti05", "kitti08", "kitti09", "kitti10",
              "kitti01", "kitti06", "kitti07"]

# KITTI profili — HEPSI OLCULDU (evidence/KITTI_KOK_NEDEN.md, faktoriyel kitti07):
#   kare basina hareket 0.707 m (EuRoC 0.022 m, 32x)  -> jump kapisi 0.25 m ile
#     karelerin %83'unde GERCEK hareket "sicrama" sayiliyordu.
#   taze stereo derinlik medyani 16.9 m (EuRoC 4.7 m) -> 19 m esigi gozlemlerin
#     %45'ini atiyordu.
# Faktoriyel: J tek basina %10.55, D tek basina %11.30, JD %4.47, JMD %3.09
# (taban %9.90). Yani parametreler AYRISMIYOR, birlikte verilmeli.
# NOT: SHA (TRUBA) bitince bu degerler outputs/sha_kitti/best_config.yaml ile
# guncellenecek; asagidakiler faktoriyel olcumun en iyisidir.
KITTI_CFG = (("frontend", "max_landmark_depth", 80.0),)
KITTI_ARGV = ("--max_init_depth", "35.9", "--jump_max_trans", "3.24",
              "--motion_trans_sigma", "1.62")


def _euroc(seq, **kw):
    return Cell(dataset="euroc", base_path=EUROC, sequence=seq, **kw)


def _kitti(seq, **kw):
    """KITTI hucresi. TartanAir kalibiyla ayni; farklar OLCULEN degerlerden gelir."""
    kw.setdefault("windows", (100, 400, 700))
    kw.setdefault("window_len", 300)
    kw.setdefault("frame_dt", 0.1)          # 10 Hz (times.txt)
    kw.setdefault("image_width", 1226)      # 00-03 icin 1241/1242; turbulans olcegi icin fark onemsiz
    kw.setdefault("stereo_baseline", 0.5372)
    kw.setdefault("cfg_overrides", KITTI_CFG)
    kw.setdefault("backend_argv", KITTI_ARGV)
    return Cell(dataset="kitti", base_path=KITTI, sequence=seq, **kw)


def _tartan(seq, **kw):
    kw.setdefault("windows", (100, 300, 500))
    kw.setdefault("window_len", 200)
    kw.setdefault("frame_dt", 0.1)
    kw.setdefault("image_width", 640)
    kw.setdefault("stereo_baseline", 0.25)
    kw.setdefault("cfg_overrides", TARTAN_CFG)
    kw.setdefault("backend_argv", TARTAN_ARGV)
    return Cell(dataset="tartanair", base_path=TARTAN, sequence=seq, **kw)


# =============================================================================
def cells_E0():
    """E0 — temel çizgi doğrulaması (türbülans kapalı)."""
    return [_euroc(s, D_over_r0=0.0, seed=0, arms=("A0",), nees=True, cov_every=10)
            for s in EUROC_SEQS]


# E4/E3/E8 tarafindan yeniden kullanilacak on-uc kosumlari (dump onbellegi):
# ayni (sekans, dr, tohum) hucresi bu deneylerde tekrar gerekiyor.
SHARED_LEVELS = (3, 5, 8)
SHARED_SEEDS = 6


def _shared(dr, k) -> bool:
    return dr in SHARED_LEVELS and k < SHARED_SEEDS


def cells_E1a(seeds=20, tilt_only_seeds=8):
    """E1-a — şiddet ekseni (+ kanal ayrıştırması: yalnız-tilt kolu)."""
    cells = [_euroc(s, D_over_r0=0.0, seed=0, arms=("A0",), cache_dump=True)
             for s in EUROC_SEQS]
    for s in EUROC_SEQS:
        for dr in DR_LEVELS:
            for k in range(seeds):
                cells.append(_euroc(s, D_over_r0=dr, seed=k, arms=("A0",),
                                    cache_dump=_shared(dr, k)))
    # Kanal ayrıştırması: geometrik (tilt) etkiyi radyometrik (bulanıklık) etkiden
    # ayırmak için aynı hücreler yalnız-tilt olarak da koşulur. C3'ün hedefi tilt
    # kanalıdır; bu kol, kazancın hangi kanaldan geldiğini gösterir.
    for s in EUROC_SEQS:
        for dr in (3, 5, 8):
            for k in range(tilt_only_seeds):
                cells.append(_euroc(s, D_over_r0=dr, seed=k, blur="none",
                                    arms=("A0",), tag="tiltonly"))
    return cells


def cells_E1b(seeds=10):
    """E1-b — zamansal rejim ekseni (D/r0 = 5 sabit)."""
    return [_euroc(s, D_over_r0=5.0, tau0_over_dt=t, seed=k, arms=("A0",))
            for s in EUROC_SEQS for t in TAU_LEVELS for k in range(seeds)]


def cells_E1c(seeds=10):
    """E1-c — uzaysal korelasyon ekseni (D/r0 = 5 sabit)."""
    return [_euroc(s, D_over_r0=5.0, theta_c_over_fov=t, seed=k, arms=("A0",))
            for s in EUROC_SEQS for t in TC_LEVELS for k in range(seeds)]


def cells_E4(seeds=10, levels=(3, 5, 8)):
    """E4/E5 — ablasyon (7 kol) + tutarlılık (NEES). Loop closure KAPALI.

    KOLLAR TEK BİR DENEYDE (2026-08-17 kararı): önceden A0C ve A4 ayrı bir deneyde
    (E4extra) koşuluyordu; bu, eşleştirilmiş testlerin ortak hücre KESİŞİMİNİ ikiye
    bölüyor ve `turb.sanity` I3 kuralını ihlal ediyordu (kollar farklı hücre
    kümelerinde ölçülürse etki büyüklükleri karşılaştırılamaz). Artık yedi kolun
    hepsi aynı hücrede, aynı türbülans gerçeklemesi üzerinde koşar.

    LOOP CLOSURE KAPALI (kasıtlı): bkz. Cell.loop_closure gerekçesi. Kapanışın
    katkıyı maskeleyip maskelemediği AYRI bir eksende (E9) sınanır — hakemin
    "gerçek bir SLAM bunu zaten düzeltirdi" itirazına kaçmadan cevap verilir.
    """
    arms = ("A0", "A0C", "A1", "A2", "A3", "A4", "A5")
    cells = [_euroc(s, D_over_r0=0.0, seed=0, arms=arms, nees=True, cov_every=10,
                    cache_dump=True) for s in EUROC_SEQS]
    for s in EUROC_SEQS:
        for dr in levels:
            for k in range(seeds):
                cells.append(_euroc(s, D_over_r0=dr, seed=k, arms=arms,
                                    nees=True, cov_every=10, cache_dump=True))
    return cells


def cells_E4_model(seeds=6, levels=(5,)):
    """Korelasyon MODELİ ablasyonu: fizik tensörü vs rapor §4.2 iso53 modeli."""
    out = []
    for s in EUROC_SEQS:
        for dr in levels:
            for k in range(seeds):
                out.append(_euroc(s, D_over_r0=dr, seed=k, arms=("A3",),
                                  corr_model="iso53", nees=True, cov_every=10,
                                  tag="iso53", cache_dump=True))
    return out


def cells_E9(seeds=6, levels=(5,)):
    """E9 — DIŞ GEÇERLİLİK: katkı, tam bir SLAM sisteminde (loop closure + PGO + GBA)
    hayatta kalıyor mu?

    NEDEN AYRI BİR EKSEN: hakemin en meşru itirazı şudur — "temel çizginiz loop
    closure'sız bir görsel odometridir; gerçek bir SLAM sistemi bu sürüklenmeyi zaten
    düzeltirdi, dolayısıyla kazancınız yapaydır." Bu itiraza iki yanlış cevap vardır:
    (a) loop closure'ı hiç konuşmamak (kaçınma), (b) ana deneyi loop closure AÇIK
    koşmak (mekanizmayı maskeleme + NEES/ATE'nin farklı kestiricilere düşmesi).
    Doğru cevap üçüncüsüdür: ayrı, KÜÇÜK ve denetimli bir eksen.

    Ölçülen soru: loop closure türbülans hasarını ne kadar maskeliyor ve maskeledikten
    SONRA türbülans-farkında kovaryansın kazancı sürüyor mu? Kapanışın kendisi de bir
    sonuçtur: kapanışlar türbülans altında seyrekleşir/yanlışlaşırsa bu raporlanır.

    Boyut: ana deneyin ~%12'si (3 sekans × 1 şiddet × 4 kol × 6 tohum).
    """
    arms = ("A0", "A1", "A3", "A5")
    seqs = ["MH01", "MH03", "MH05"]
    cells = [_euroc(s, D_over_r0=0.0, seed=0, arms=arms, nees=True, cov_every=10,
                    loop_closure=True, tag="loop", cache_dump=True) for s in seqs]
    for s in seqs:
        for dr in levels:
            for k in range(seeds):
                cells.append(_euroc(s, D_over_r0=float(dr), seed=k, arms=arms,
                                    nees=True, cov_every=10, loop_closure=True,
                                    tag="loop", cache_dump=True))
    return cells


def cells_tartan(seeds=8, levels=(0, 3, 5, 8)):
    """TartanAir — dış mekân / fiziksel olarak makul rejim (dış geçerlilik)."""
    cells = []
    for s in TARTAN_SEQS:
        for dr in levels:
            for k in (range(seeds) if dr > 0 else [0]):
                cells.append(_tartan(s, D_over_r0=float(dr), seed=k, arms=("A0",)))
    return cells


def cells_tartan_arms(seeds=5, levels=(5,)):
    arms = ("A0", "A1", "A3", "A5")
    return [_tartan(s, D_over_r0=float(dr), seed=k, arms=arms, nees=True, cov_every=10)
            for s in TARTAN_SEQS for dr in levels for k in range(seeds)]


def cells_kitti(seeds=8, levels=(0, 3, 5, 8)):
    """KITTI — GERCEK dis mekan bozulma haritasi (A0 tek kol).

    TartanAir'in (sentetik) yerine GERCEK kamera/arac verisiyle ayni egriyi cikarir.
    Makalenin dis gecerlilik iddiasi bunun uzerine kurulur.
    """
    return [_kitti(s, D_over_r0=float(dr), seed=k, arms=("A0",))
            for s in KITTI_SEQS for dr in levels
            for k in (range(seeds) if dr > 0 else [0])]


def cells_kitti_arms(seeds=5, levels=(5,)):
    """KITTI — kol karsilastirmasi (A0/A1/A3/A5), NEES ile.

    D/r0=5 secildi: bozulma haritasinda mimari rejim degisiminin oldugu nokta
    (bkz. evidence/FIRSATCI_STEREO.md). Kollarin AYNI hucrelerde eslestirilmis
    kosmasi sarttir (sanity I2-hayatta-kalma-yanliligi).
    """
    arms = ("A0", "A1", "A3", "A5")
    return [_kitti(s, D_over_r0=float(dr), seed=k, arms=arms, nees=True, cov_every=10)
            for s in KITTI_SEQS for dr in levels for k in range(seeds)]


def cells_E3(seeds=5, levels=(3, 5, 8),
             methods=("none", "temporal_median", "lucky", "elastic"),
             seqs=("MH01", "MH03", "MH05")):
    """E3 — restorasyon boru hattı (C2: PSNR ≠ ATE).

    Her hücre hem ALGISAL kaliteyi (PSNR/SSIM, temiz kareye göre) hem de
    NAVİGASYON hatasını (ATE/RPE) aynı koşumda kaydeder; C2'nin ana analizi bu iki
    eksenin saçılımıdır. Referans kolları: "none" (bozuk -> VIO) ve temiz (üst sınır).
    """
    cells = [_euroc(s, D_over_r0=0.0, seed=0, arms=("A0",), tag="clean") for s in seqs]
    for s in seqs:
        for dr in levels:
            for meth in methods:
                for k in range(seeds):
                    cells.append(_euroc(s, D_over_r0=float(dr), seed=k, arms=("A0",),
                                        restoration=meth, tag=f"res_{meth}",
                                        cache_dump=(meth == "none")))
    return cells



def cells_E3win(seeds=3, levels=(3, 5), windows=(3, 5, 7),
                methods=("temporal_median", "elastic"),
                seqs=("MH01", "MH03", "MH05")):
    """E3-win — restorasyon PENCERE BOYU süpürmesi (C2 kök nedeni).

    NEDEN (evidence/C2_KOK_NEDEN.md): E3'te üç restorasyon yöntemi de hem PSNR'ı
    hem ATE'yi kötüleştirdi. Kök neden ölçüldü: klasik türbülans giderme STATİK
    kamera varsayar, oysa EuRoC'ta 7 karelik (0.35 s) füzyon penceresinde kamera
    16-34 cm ötelenip 2-3.5° dönüyor — stereo taban çizgisi (11 cm) BUNDAN KÜÇÜK.
    Kontrol deneyi (`turb/_c2_probe.py`) kanıtladı: aynı yöntemler TÜRBÜLANS HİÇ
    YOKKEN de PSNR'ı 9-13 dB'ye düşürüyor. Yani hasar hareketten.

    Bu süpürme şunu sınar: pencereyi kısaltmak kurtarır mı?
        7 kare = 0.35 s -> ~16-34 cm öteleme
        5 kare = 0.25 s -> ~11-24 cm
        3 kare = 0.15 s -> ~ 7-14 cm
    3 karede de çöküyorsa sonuç "HİÇBİR pencere boyu işe yaramıyor" olur ve
    "zamansal füzyon çevik platformla mimari olarak bağdaşmıyor" tezini
    ÖLÇÜMLE kurar. Kurtarıyorsa eşik bulunur — o da bir sonuçtur.

    `lucky` DIŞARIDA: yapısal kusuru pencere boyundan bağımsız (seçilen kare
    MEVCUT kare olmayabilir -> zaman ekseninde kayma). Pencereyi kısaltmak onu
    düzeltmez, yalnız seçim havuzunu daraltır.
    """
    cells = []
    for s in seqs:
        for dr in levels:
            # referans: restorasyonsuz (pencereden bağımsız, bir kez)
            for k in range(seeds):
                cells.append(_euroc(s, D_over_r0=float(dr), seed=k, arms=("A0",),
                                    restoration="none", tag="win_none",
                                    cache_dump=True))
            for w in windows:
                for meth in methods:
                    for k in range(seeds):
                        cells.append(_euroc(
                            s, D_over_r0=float(dr), seed=k, arms=("A0",),
                            restoration=meth, tag=f"win{w}_{meth}",
                            cfg_overrides=(("turbulence", "restoration_window", int(w)),),
                            cache_dump=False))
    return cells


def cells_E8(seeds=6, levels=(0, 3, 5, 8)):
    """E8-a — genellenebilirlik: bozulma ON-UC secimine ozgu mu?

    Rapor E8-a, turbulansin BUTUN sistemleri bozdugunu gostermeyi ister ("bu bir hat
    kusuru degil"). Burada ayni bozulma, hattin ON-UC yapilandirmasi degistirilerek
    kosulur:
        (i)  ALIKED + LightGlue      (ana hat)
        (ii) SuperPoint + LightGlue  (farkli dedektor/betimleyici)
        (iii) ALIKED + klasik BF     (ogrenilmis eslestirici yok)
    Ayni D/r0 ekseninde ucunun de bozulmasi, etkinin dedektor/eslestirici secimine
    degil goruntu olusumuna ait oldugunu gosterir. (Tam bagimsiz sistemler —
    ORB-SLAM3 / VINS / OpenVINS — icin turb.make_dataset ile diske yazilan
    Turb-VIO sekanslari kullanilir; bkz. turb/make_dataset.py.)
    """
    variants = {
        "aliked_lg": (),
        "superpoint_lg": (("frontend", "extractor", "superpoint"),),
        "aliked_bf": (("frontend", "matcher", "classical_bf"),),
    }
    seqs = ["MH01", "MH03", "MH05"]
    out = []
    for name, ov in variants.items():
        for sq in seqs:
            for dr in levels:
                for k in (range(seeds) if dr > 0 else [0]):
                    out.append(_euroc(sq, D_over_r0=float(dr), seed=k, arms=("A0",),
                                      cfg_overrides=ov, tag=f"fe_{name}",
                                      cache_dump=(name == "aliked_lg")))
    return out


# =============================================================================
#  FAZ 2 (ON_KAYIT_FAZ2.md §1) — hareket-seçimli pencereler, fiziksel PSF
# =============================================================================
# Veri kökü: bu makinede EUROC sabiti (yukarıda) olmayan bir yolu gösteriyor; Faz 2 hücreleri
# $SLAM_DATA_ROOT'u (yoksa yerel diski) kullanır. Eski deney tanımları DEĞİŞMEDİ.
EUROC_FAZ2 = os.environ.get("SLAM_DATA_ROOT", "/data/sdata")
FAZ2_SEQS = ["MH01", "MH02", "MH03", "MH04", "MH05",
             "V1_01_easy", "V1_02_medium", "V2_01_easy", "V2_02_medium"]
FAZ2_DR = (0, 4, 8, 12, 16, 20, 24)
FAZ2_ARMS = ("A0", "A0C", "A1", "A3", "AG", "A0C+A3")
FAZ2_A5_SEQS = ("MH02", "MH04")           # A5 (oracle yama, yavaş) yalnız bu alt kümede
# A0C / A0C+A3 taban çarpanı k: V-serisi temiz koşumlarında örneklem dışı kalibre edilir
# (ON_TEST_FAZ2_EUROC.md T5). config.yaml'daki 16 ESKİ hattın değeridir; burada CLI ile ezilir.
# 2026-11-08 ön-test T5 (4060 Ti, cov_every=10, gürültüsüz): k=12 -> 2.54.
# 2026-11-08 ön-test R5 (YENİDEN, uzun-koşum yapılandırması: RTX 3090, cov_every=2, sigma=1 DN, tohum 0;
# V1_01/V1_02/V2_01/V2_02 D/r0=0, 12 pencere, gerçek A0C koşumları): k=12 -> pencere medyanlarının medyanı
# 2.38 (k=8: 4.40, k=16: 1.80). Seçilen k = 12 (değişmedi).
FAZ2_K_A0C = 12.0


# Kararlar (2026-11-08, ON_TEST_FAZ2_EUROC.md §4.10/§4.11/§5; kullanıcı onayı koordinatör üzerinden):
#   1) her hücre (D/r0=0 dahil) sigma=1 DN sensör gürültüsü, gürültü tohumu = hücre tohumu, tohum 0-4;
#      eşleme (sekans, tohum, pencere); 2) cov_every=2; 3) pencere listesi pencereler_faz2_kapsam1m.json;
#   5) bütün ızgara TEK kartta (RTX 3090) — bkz. turb/faz2_kos.py.
FAZ2_NOISE_DN = 1.0
FAZ2_COV_EVERY = 2
FAZ2_WINDOWS_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                 "evidence", "pencereler_faz2_kapsam1m.json")


def cells_FAZ2(seeds=5, k_a0c=None, windows_file=None):
    """Faz 2 ana ızgarası: 9 sekans x 7 D/r0 x `seeds` tohum; D/r0=0 da tohumlu (sensör gürültüsü)."""
    from .runner import windows_for
    k = k_a0c if k_a0c is not None else FAZ2_K_A0C
    if k is None:
        raise SystemExit("FAZ2_K_A0C tanımsız: önce A0C k kalibrasyonu (ON_TEST T5)")
    wf = windows_file or FAZ2_WINDOWS_FILE
    argv = ("--turb_cal_factor", str(float(k)))
    cells = []
    for s in FAZ2_SEQS:
        w = windows_for(s, wf)
        arms = FAZ2_ARMS + (("A5",) if s in FAZ2_A5_SEQS else ())
        for dr in FAZ2_DR:
            for sd in range(seeds):
                cells.append(Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence=s, seed=sd,
                                  D_over_r0=float(dr), arms=arms, windows=w, window_len=400,
                                  nees=True, cov_every=FAZ2_COV_EVERY, psf_legacy_scale=False,
                                  sensor_noise_dn=FAZ2_NOISE_DN, noise_seed=None,
                                  backend_argv=argv, tag="faz2", cache_dump=False))
    return cells


# =============================================================================
#  FAZ 2 küçük A0 ızgaraları (ON_KAYIT §1; ON_TEST_FAZ2_EUROC.md §7)
# =============================================================================
# Ortak kurallar (v0.3): sigma=1 DN, gürültü tohumu = hücre tohumu, tohum 0-4, kapsam1m pencereleri,
# tek kart RTX 3090, eşli (sekans, tohum, pencere). Etiket "faz2": böylece ana ızgarayla AYNI
# yapılandırmadaki hücrelerin anahtarı ana ızgaranınkiyle aynıdır (yeniden koşulmaz, eşlenir).

FAKT_KOSULLAR = {            # ad -> (tilt, blur)
    "yok": (False, "none"),
    "tilt": (True, "none"),
    "blur": (False, "speckle_global"),
    "ikisi": (True, "speckle_global"),
}


def cells_FAZ2_faktoriyel(seeds=5, levels=(8, 16), seqs=None, kosullar=("tilt", "blur")):
    """H2 tilt x blur faktöriyeli, yalnız A0.

    VARSAYILAN yalnız "tilt" ve "blur" koşullarını üretir: "ikisi" ana ızgaranın A0 hücresiyle
    (aynı anahtar, aynı yapılandırma) ve "yok" ana ızgaranın D/r0=0 hücresiyle özdeştir (ön-test §7:
    D/r0>0'da kanal yokken kanca çıktısı D/r0=0 gürültülü çıktıyla bayt-aynı). Dört koşulun hepsi
    istenirse kosullar=tuple(FAKT_KOSULLAR).
    """
    cells = []
    for s in (seqs or FAZ2_SEQS):
        w = __import__("turb.runner", fromlist=["windows_for"]).windows_for(s, FAZ2_WINDOWS_FILE)
        for dr in levels:
            for kos in kosullar:
                tilt, blur = FAKT_KOSULLAR[kos]
                for sd in range(seeds):
                    cells.append(Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence=s, seed=sd,
                                      D_over_r0=float(dr), tilt=tilt, blur=blur, arms=("A0",),
                                      windows=w, window_len=400, nees=True, cov_every=FAZ2_COV_EVERY,
                                      psf_legacy_scale=False, sensor_noise_dn=FAZ2_NOISE_DN,
                                      backend_argv=("--turb_cal_factor", str(float(FAZ2_K_A0C))),
                                      tag="faz2", cache_dump=False))
    return cells


def butce_overrides(kp: int, replenish_ratio: float):
    """max_keypoints / replenish_ratio_kp değişikliğinin TÜRETİLMİŞ değerleriyle birlikte cfg_overrides'ı.

    config.load() oranları yükleme anında uygular (config._scale_with_budget); runner'ın cfg_overrides'ı
    ondan SONRA yazıldığı için yalnız max_keypoints'i değiştirmek T / N_rep / stereo_uzak_min'i ESKİ
    bütçede bırakırdı. Burada aynı formüller config.yaml oranlarıyla yeniden hesaplanır. kp > butce_kare_esik
    (128) olduğundan anahtar kare penceresi değişmez (kp <= 128 desteklenmez).
    """
    import config as _c
    import yaml
    with open(_c.DEFAULT_PATH) as f:
        raw = yaml.safe_load(f)
    fr, br = raw["frontend"], raw["backend"]
    if kp <= int(br.get("butce_kare_esik", 0) or 0):
        raise ValueError("kp <= butce_kare_esik: kare-pencere ayarları da değişir; desteklenmiyor")
    return (("frontend", "max_keypoints", int(kp)),
            ("frontend", "replenish_ratio_kp", float(replenish_ratio)),
            ("frontend", "target_feature_count", round(float(fr["target_ratio_kp"]) * kp)),
            ("frontend", "replenish_threshold", round(float(replenish_ratio) * kp)),
            ("backend", "stereo_uzak_min", round(float(br["stereo_uzak_min_ratio_kp"]) * kp)))


H6_KP = (256, 512)
H6_REPL = (0.2, 0.4)


def cells_FAZ2_h6(kirilma, seeds=5, seqs=None, varsayilan_dahil=False):
    """H6 nedensel müdahale: D/r0 in {kırılma-4, kırılma, kırılma+4} x replenish {0.2,0.4} x kp {256,512}, A0.

    (0.2, 256) = PyEco varsayılanı = ana ızgaranın A0 hücresi; varsayılan olarak üretilmez (eşlenir).
    Etiket her (r, kp) için ayrıdır (Cell.key cfg_overrides'ı içermez; çakışmayı önler).
    """
    levels = [float(kirilma) + d for d in (-4, 0, 4)]
    if min(levels) < 0:
        raise ValueError("kırılma - 4 < 0")
    cells = []
    for s in (seqs or FAZ2_SEQS):
        w = __import__("turb.runner", fromlist=["windows_for"]).windows_for(s, FAZ2_WINDOWS_FILE)
        for r in H6_REPL:
            for kp in H6_KP:
                if (r, kp) == (0.2, 256) and not varsayilan_dahil:
                    continue
                ov = butce_overrides(kp, r)
                for dr in levels:
                    for sd in range(seeds):
                        cells.append(Cell(dataset="euroc", base_path=EUROC_FAZ2, sequence=s, seed=sd,
                                          D_over_r0=dr, arms=("A0",), windows=w, window_len=400,
                                          nees=True, cov_every=FAZ2_COV_EVERY, psf_legacy_scale=False,
                                          sensor_noise_dn=FAZ2_NOISE_DN, cfg_overrides=ov,
                                          backend_argv=("--turb_cal_factor", str(float(FAZ2_K_A0C))),
                                          tag=("faz2" if (r, kp) == (0.2, 256) else f"faz2_h6_r{r:g}_kp{kp}"),
                                          cache_dump=False))
    return cells


REGISTRY = {
    "FAZ2_FAKT": cells_FAZ2_faktoriyel,
    "FAZ2": cells_FAZ2,
    "E0": cells_E0,
    "E1a": cells_E1a,
    "E1b": cells_E1b,
    "E1c": cells_E1c,
    "E3": cells_E3,
    "E3win": cells_E3win,
    "E4": cells_E4,
    "E4model": cells_E4_model,
    "E8": cells_E8,
    "E9": cells_E9,
    "TARTAN": cells_tartan,
    "TARTANarms": cells_tartan_arms,
    "KITTI": cells_kitti,            # gercek dis mekan bozulma haritasi (A0)
    "KITTIarms": cells_kitti_arms,   # kol karsilastirmasi + NEES
}


def main():
    ap = argparse.ArgumentParser(description="Türbülans-VIO deney koşucusu")
    ap.add_argument("--exp", required=True, choices=sorted(REGISTRY))
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--gpu_jobs", type=int, default=2)
    ap.add_argument("--chunk", type=int, default=24)
    ap.add_argument("--legacy", action="store_true",
                    help="eski tek-havuz (GPU+CPU aynı işçide) koşum yolu")
    ap.add_argument("--seeds", type=int, default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--no_resume", action="store_true")
    # KOL SUZGECI: yalnizca belirtilen kollari kos. DIMM duzeltmesi (2026-08-23)
    # SADECE cevrimici kestirim kullanan kollari etkiler; `current_sigma_px`
    # A1 ve A5 icin oracle'i DOGRUDAN dondurur ve DIMM'e hic dokunmaz, A0 ise
    # turb arka-ucunu hic kullanmaz. Dolayisiyla E4'un 7 kolunu bastan kosmak
    # yerine 3'unu (A2, A3, A4) kosmak yeterlidir: ~20 sa yerine ~8.5 sa.
    ap.add_argument("--only_arms", default="",
                    help="virgulle ayrilmis kol listesi, or. A2,A3,A4")
    a = ap.parse_args()

    fn = REGISTRY[a.exp]
    cells = fn(a.seeds) if a.seeds is not None else fn()
    if a.only_arms:
        keep = tuple(x.strip() for x in a.only_arms.split(",") if x.strip())
        cells = [c.__class__(**{**c.__dict__,
                                "arms": tuple(x for x in c.arms if x in keep)})
                 for c in cells]
        cells = [c for c in cells if c.arms]
        print(f"[{a.exp}] kol suzgeci {keep} -> {len(cells)} hucre")
    out = a.out or f"evidence/{a.exp}/records.jsonl"
    print(f"[{a.exp}] {len(cells)} hücre -> {out}")
    if a.legacy:
        run_experiment(cells, out, jobs=a.jobs, gpu_jobs=a.gpu_jobs,
                       resume=not a.no_resume)
    else:
        # İKİ-FAZLI HAT (varsayılan): GPU ön-uçları ayrı havuzda seri/sıralı,
        # CPU arka-uçları ayrı havuzda paralel; yığınlar örtüşür.
        run_experiment_pipelined(cells, out, cpu_jobs=a.jobs, gpu_jobs=a.gpu_jobs,
                                 chunk=a.chunk, resume=not a.no_resume)


if __name__ == "__main__":
    main()
