"""
TAHLİL KATMANI — sonuçların fizik ve literatürle tutarlılığının otomatik denetimi.

Amaç (kullanıcı direktifi, 2026-08-15): "makale yanlış sonuçlar ile literatüre
meydan okumasın". Bir deney bittiğinde sonucu HEMEN sınayan, hatayı haftalar sonra
makale yazarken değil o anda yakalayan bir denetçi. Kuyruk (run_queue_v3.sh) her
deneyden sonra bunu çağırır; KIRMIZI bulgu varsa çıkış kodu 2'dir.

Denetimin FELSEFESİ — üç ayrı soru sorulur ve karıştırılmaz:

  (S) SİMÜLATÖR DOĞRU MU?   Ölçülen fizik kapalı-form literatür değerine uyuyor mu?
                            (Kolmogorov/Fried/Tyler sabitleri — tolerans dar.)
  (D) DAVRANIŞ MAKUL MÜ?    Bozulma-hata haritası niteliksel olarak literatürün
                            beklediği yönde mi? (monotonluk, kanal sırası, temiz
                            temel çizginin yayımlanmış EuRoC aralığında olması.)
  (I) İDDİA DÜRÜST MÜ?      İstatistik, sonucu olduğundan güçlü göstermiyor mu?
                            (hayatta-kalma yanlılığı, dejenere pencere, tutarlılığın
                            doğruluk feda edilerek satın alınması, literatürde
                            görülmemiş büyüklükte etki.)

(S) ihlali = YAZILIM HATASI. (D) ihlali = ya yazılım hatası ya gerçek bir keşif —
her hâlde açıklanmadan yayımlanamaz. (I) ihlali = MAKALE HATASI; sonuç doğru olsa
bile sunumu yanıltıcıdır.

Kullanım:
    python -m turb.sanity --exp E4
    python -m turb.sanity --all           # bütün deneyler + fizik denetimi
    python -m turb.sanity --physics       # yalnız simülatör fiziği (hızlı)

Çıkış kodu: 0 = temiz, 1 = SARI uyarı var, 2 = KIRMIZI bulgu var.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

EVID = "evidence"
ALERT = os.path.join(EVID, "SANITY_ALERT.md")

OK, WARN, RED = "OK", "SARI", "KIRMIZI"
ACK = "KABUL"          # bulgu ele alındı; gerekçesiyle birlikte görünür kalır

# Ele alınmış bulguların kaydı. Faz 2.5 kapısı (makale üretim skilli,
# references/sonuc-tahlili.md) bir kırmızı bulgunun üç yoldan biriyle kapanmasını
# ister: (a) düzeltilir, (b) ölçüm artefaktı gösterilip KURAL düzeltilir, (c)
# gerekçelendirilerek makalede sınırlılık olarak yazılacağına karar verilir.
# Bu dosya (c) ve (a)'nın kaydıdır: {kural_kimligi: {"nasil": ..., "gerekce": ...}}.
#
# GÜVENLİK: kabul, bulguyu SİLMEZ. Gerekçesiyle raporda kalır ve rapor genel durumu
# en iyi ihtimalle SARI olur — böylece "kapatıp unutma" mümkün değildir.
ACK_FILE = os.path.join(EVID, "sanity_ack.json")


# =============================================================================
@dataclass
class Finding:
    rule: str                 # kural kimliği, ör. "S1-tilt-katsayisi"
    group: str                # S | D | I
    level: str                # OK | SARI | KIRMIZI
    claim: str                # neyin sınandığı
    measured: str             # ölçülen
    expected: str             # literatür/fizik beklentisi
    reference: str            # kaynak
    action: str = ""          # KIRMIZI/SARI ise ne yapılmalı

    ack: dict | None = None   # {"nasil": ..., "gerekce": ...}

    def md(self) -> str:
        icon = {OK: "✅", WARN: "⚠️", RED: "🛑", ACK: "🟩"}[self.level]
        s = (f"### {icon} `{self.rule}` — {self.claim}\n"
             f"- **Ölçülen:** {self.measured}\n"
             f"- **Beklenen:** {self.expected}\n"
             f"- **Kaynak:** {self.reference}\n")
        if self.action:
            s += f"- **Yapılacak:** {self.action}\n"
        if self.ack:
            s += (f"- **Ele alındı ({self.ack.get('nasil','?')}):** "
                  f"{self.ack.get('gerekce','')}\n")
        return s


@dataclass
class Report:
    findings: list = field(default_factory=list)

    acks: dict = field(default_factory=dict)

    def add(self, *a, **kw):
        f = Finding(*a, **kw)
        # Kabul edilmiş bulgu: seviyesi KABUL'e düşer ama metinde kalır.
        if f.level in (RED, WARN) and f.rule in self.acks:
            f.ack = self.acks[f.rule]
            f.level = ACK
        self.findings.append(f)

    @property
    def worst(self) -> str:
        lv = [f.level for f in self.findings]
        if RED in lv:
            return RED
        if WARN in lv or ACK in lv:
            return WARN
        return OK

    def exit_code(self) -> int:
        return {OK: 0, WARN: 1, RED: 2}[self.worst]


def load_acks() -> dict:
    if not os.path.exists(ACK_FILE):
        return {}
    try:
        d = json.load(open(ACK_FILE))
        return {k: v for k, v in d.items() if isinstance(v, dict)}
    except (json.JSONDecodeError, OSError):
        return {}


# =============================================================================
# Literatür sabitleri — her biri kaynağıyla birlikte (denetimin çıpası)
# =============================================================================
from .metrics import NEES_BAND as _NEES_BAND

LIT = {
    # Kolmogorov faz yapı fonksiyonu sabiti
    "D_phi_const": (6.88, "Fried 1965; Roddier 1981 Prog. Opt. 19"),
    # G-tilt (gradyan eğimi) tek-eksen varyans katsayısı:
    #   sigma_a^2 = C (D/r0)^{5/3} (lam/D)^2 ,  C_G = 0.170
    "C_gtilt": (0.170, "Tyler 1994 JOSA A 11:358; Lofdahl 2010 (K_G=0.340/2)"),
    "C_ztilt": (0.1816, "Noll 1976 JOSA 66:207 (Zernike tilt)"),
    # Eğim varyansının D/r0 üsteli
    "tilt_exponent": (5.0 / 3.0, "Kolmogorov spektrumu, f^{-11/3}"),
    # Yayımlanmış stereo-VIO ATE aralığı, EuRoC MH (temiz veri, tam sekans):
    # VINS-Fusion 0.08-0.30 m, ORB-SLAM3 0.03-0.09 m, OKVIS 0.16-0.40 m
    "euroc_ate_band": ((0.01, 0.60), "Campos 2021 ORB-SLAM3 T.5; Qin 2018 VINS-Mono T.2; "
                                     "Delmerico 2018 ICRA benchmark"),
    # Tutarlı bir kestiricinin NEES'i, 3 SD icin chi2 kabul bandi (%99, N~20 pencere)
    # TEK KAYNAK: turb/metrics.py::NEES_BAND (dof=3, iki yanli %95).
    "nees_band": (_NEES_BAND, "Bar-Shalom 2001 §5.4 (NEES chi-kare tutarlilik testi), "
                              "tek ornek icin iki yanli %95 bolgesi, dof=3"),
    "coverage_nominal": (0.95, "%95 guven elipsoidi tanimi"),
    # Korelasyonlu ölçüm gürültüsünün düzeltilmesinden literatürde bildirilen
    # tutarlılık kazancı büyüklüğü (kat olarak). Bunun ÇOK ötesi şüphelidir.
    "consistency_gain_plausible": (50.0, "Julier&Uhlmann 2001 (kovaryans kesisimi); "
                                         "Huang 2010 FEJ; Barrau&Bonnabel 2017 IEKF — "
                                         "bildirilen NEES iyilesmeleri tipik 2-20 kat"),
    # YAPI kazancinin TEORIK ust siniri: A1 (skaler sisirme) ile A5 (yapili) AYNI
    # kosegene sahiptir (turb/backend.py ölçek/nugget ayrimi), dolayisiyla fark
    # YALNIZ kosegen-disi yapidan gelir. Kapali-formda olculdu (evidence/
    # E2_ETKI_BUYUKLUGU.md): yamada K=3-8 gozlemle kazanc 1.11-1.83 kat. Zamansal
    # korelasyon etkin ornek sayisini daha da dusurdugu icin gercek hatta birkac
    # kat gorulebilir; ama 10 katin ustu bu mekanizmayla ACIKLANAMAZ.
    "structure_gain_bound": (10.0, "kapali-form turetim, evidence/E2_ETKI_BUYUKLUGU.md "
                                   "(K=3-8, sigma_det/sigma_tilt=1.75): 1.11-1.83 kat"),
}


def _band(x, lo, hi):
    return (x is not None) and np.isfinite(x) and (lo <= x <= hi)


# =============================================================================
# (S) SİMÜLATÖR FİZİĞİ
# =============================================================================
def check_physics(rep: Report, quick: bool = True):
    """Simülatörün ürettiği istatistikleri kapalı-form literatüre karşı sınar.

    Bu denetim, deney sonuçlarından BAĞIMSIZDIR: doğrudan turb.screen/fields
    çağrılır. Böylece "sonuç güzel çıktı ama simülatör bozuktu" senaryosu kapanır.
    """
    from .optics import sigma_tilt_rad
    from .screen import phase_screen, structure_function_shift
    from .fields import TiltField
    from .simulator import TurbParams
    from . import correlation as corr

    P = TurbParams()                      # HATTIN GERÇEK AYARLARI ile sınanır
    D, lam = P.aperture_D, P.wavelength
    ppd = P.screen_px_per_D
    delta = D / ppd
    N = 1024 if quick else P.screen_N     # hızlı kipte küçük ekran (ölçek değişmez)
    rng = np.random.default_rng(7)

    def _tilt_field(dr0, cal):
        return TiltField(N, delta, r0=D / dr0, D=D, lam=lam, rng=rng,
                         n_subharmonics=P.subharmonics, calibration=cal)

    # --- S1: hattın ÜRETTİĞİ eğim varyansı kapalı-forma uyuyor mu? ------------
    # Ölçüt, "ham ekran" değil hattın gerçekten kullandığı (kalibrasyonlu) alandır:
    # E1/E4 sonuçlarını belirleyen budur.
    n_real = 4 if quick else 16
    dr0 = 5.0
    cal = P.tilt_calibration
    meas = float(np.mean([_tilt_field(dr0, cal).tilt_std() for _ in range(n_real)]))
    theo = sigma_tilt_rad(dr0, lam, D, kind="gtilt")
    ratio = meas / theo if theo > 0 else np.nan
    tol = 0.10
    rep.add("S1-tilt-varyansi", "S", OK if abs(ratio - 1) < tol else RED,
            "Hattın ürettiği açıklık-ortalamalı G-tilt varyansı kapalı-form değere uyuyor mu?",
            f"sigma_hat/sigma_teori = {ratio:.3f}  (D/r0={dr0}, N={N}, ppd={ppd}, "
            f"kalibrasyon={cal:g}, n={n_real})",
            f"1.00 ± {tol:.2f};  sigma_a^2 = {LIT['C_gtilt'][0]} (D/r0)^(5/3) (lam/D)^2",
            LIT["C_gtilt"][1],
            "" if abs(ratio - 1) < tol else
            "turb/fields.py TiltField transfer fonksiyonunu ve turb/screen.py spektral "
            "normalizasyonunu denetle — E1/E4'ün BÜTÜN sonuçları geçersizdir.")

    # --- S1b: kalibrasyon sabiti savunulabilir mi? ----------------------------
    # Kalibrasyon bir "uydurma katsayı" DEĞİLDİR: kaba ızgaranın (delta = D/ppd)
    # ayrıklaştırma açığını kapatır ve bu açık S3'te bağımsız olarak ölçülür.
    # Sabit 1.25'i aşarsa artık düzeltme değil maskeleme olur.
    raw = float(np.mean([_tilt_field(dr0, 1.0).tilt_std() for _ in range(n_real)])) / theo
    okc = 0.80 <= cal <= 1.25 and abs(raw * cal - 1) < 0.12
    rep.add("S1b-kalibrasyon-mesruiyeti", "S", OK if okc else RED,
            "Eğim kalibrasyon sabiti, ölçülen ayrıklaştırma açığıyla açıklanıyor mu?",
            f"ham oran = {raw:.3f};  1/ham = {1/max(raw,1e-9):.3f};  "
            f"kullanılan kalibrasyon = {cal:g};  çarpım = {raw*cal:.3f}",
            "kalibrasyon in [0.80, 1.25] ve ham x kalibrasyon = 1.00 ± 0.12",
            "Schmidt 2010 'Numerical Simulation of Optical Wave Propagation' §9 "
            "(ayrık faz ekranlarında düşük/yüksek frekans açığı ve düzeltmesi)",
            "" if okc else
            "Kalibrasyon açığı açıklamıyor: ya ekran çözünürlüğü (screen_px_per_D) "
            "artırılmalı ya da katsayı makalede gerekçelendirilemez.")

    # --- S2: (D/r0)^{5/3} ölçekleme üsteli ------------------------------------
    drs = np.array([2.0, 4.0, 8.0])
    sig2 = [float(np.mean([_tilt_field(d, cal).tilt_std()
                           for _ in range(max(2, n_real // 2))])) for d in drs]
    slope = float(np.polyfit(np.log(drs), np.log(np.array(sig2)), 1)[0])
    exp_meas = 2 * slope                       # varyans üsteli
    ok = abs(exp_meas - 5 / 3) < 0.20
    rep.add("S2-olcekleme-usteli", "S", OK if ok else RED,
            "Eğim varyansı (D/r0)^{5/3} ile ölçekleniyor mu?",
            f"ölçülen üstel = {exp_meas:.3f}  (D/r0 = 2, 4, 8)",
            "5/3 = 1.667 ± 0.20", LIT["tilt_exponent"][1],
            "" if ok else "Kolmogorov spektrum üsteli (-11/3) veya r0 ölçeklemesi hatalı.")

    # --- S3: faz yapı fonksiyonu ve AÇIĞIN KANALLARA YAYILIMI -----------------
    # Ayrık FFT faz ekranları Kolmogorov gücünü sistematik olarak EKSİK üretir
    # (Schmidt 2010 §9; subharmonikler açığı azaltır, kapatmaz). Burada sorulan
    # "açık var mı" değil — vardır ve bilinir — "açık İDDİA ETTİĞİMİZ niceliklere
    # ne kadar geçiyor" sorusudur:
    #   * EĞİM kanalı: açık S1b'de ölçülür ve kalibrasyonla kapatılır -> etkisiz.
    #   * BULANIKLIK kanalı: kalibre EDİLMEZ -> bulanıklık bu oran kadar iyimserdir.
    # Bu yüzden sonuç, bulanıklık kanalının sınırlılığı olarak raporlanır.
    r0 = D / 4.0
    n_sf = 3 if quick else 8

    def _sf(dl, n):
        acc, rr = None, None
        for _ in range(n_sf):
            phz = phase_screen(n, dl, r0, rng, n_subharmonics=P.subharmonics)
            rr, Dm = structure_function_shift(phz, dl, range(1, min(40, n // 8)))
            acc = Dm if acc is None else acc + Dm
        Dm = acc / n_sf
        # geçerli aralık: ızgaranın çözebildiği en küçük ölçeğin üstü, dış-ölçek
        # etkisinin altı. Kaba ızgarada alt sınır örnekleme aralığı belirler.
        sel = (rr >= max(2.0 * dl, 0.15 * r0)) & (rr < 4.0 * r0)
        if sel.sum() < 3:
            return float("nan")
        return float(np.median(Dm[sel] / (6.88 * (rr[sel] / r0) ** (5 / 3))))

    try:
        # ekran GENİŞLİĞİ iki ızgarada eşit tutulur (64 D) — açık yalnız
        # örnekleme aralığına bağlı kalsın diye
        ext = 64.0 * D
        fine = _sf(D / P.fine_px_per_D, int(round(ext / (D / P.fine_px_per_D))))
        coarse = _sf(delta, int(round(ext / delta)))
        gap = 1.0 - fine
        lvl = OK if abs(gap) < 0.10 else (RED if abs(gap) > 0.25 else WARN)
        rep.add("S3-yapi-fonksiyonu", "S", lvl,
                "Faz yapı fonksiyonu D_phi(r) = 6.88 (r/r0)^{5/3} — ekran güç açığı",
                f"ince ızgara (D/{P.fine_px_per_D}): {fine:.3f} · "
                f"kaba ızgara (D/{ppd}): {coarse:.3f}  (0.15 r0 < r < 1.5 r0, "
                f"ekran genişliği {ext/D:g} D)",
                "1.00 ± 0.10 (ayrık ekranlarda %10'a kadar açık literatürde beklenir; "
                "> %25 ise normalizasyon hatası)",
                LIT["D_phi_const"][1] + "; Schmidt 2010 §9 (ayrık ekran güç açığı)",
                "" if lvl == OK else
                "Açık %25'i aşıyorsa screen.py spektral normalizasyonu hatalıdır.")
        # açığın bulanıklık kanalına yayılımı — makalede SINIRLILIK olarak yazılacak
        if np.isfinite(fine) and abs(gap) > 0.05:
            rep.add("S3c-bulaniklik-kanali-iyimser", "S", WARN,
                    "Ekran güç açığı bulanıklık kanalına doğrudan geçiyor "
                    "(eğim kanalı kalibreli, bulanıklık değil)",
                    f"faz varyansı teoriden %{100*gap:.1f} eksik -> kısa-poz PSF "
                    f"genişliği ~%{100*(1-(fine)**0.6):.1f} dar (D_phi ~ (D/r0)^{{5/3}})",
                    "sıfır açık (ulaşılamaz) — bu yüzden BİLDİRİLİR, düzeltilmez",
                    "Roggemann & Welsh 1996 §3.5 (kısa-poz PSF ve r0 bağımlılığı)",
                    "Makalenin sınırlılıklar bölümüne: 'bulanıklık kanalı beyan edilen "
                    "D/r0'a göre bir miktar İYİMSERDİR; bildirilen bozulma bu yüzden "
                    "gerçek atmosferin ALT sınırıdır' — bu, iddiayı zayıflatmaz, "
                    "tutuculaştırır.")
    except Exception as e:                                   # pragma: no cover
        rep.add("S3-yapi-fonksiyonu", "S", WARN, "Faz yapı fonksiyonu",
                f"ölçülemedi: {type(e).__name__}: {e}", "—", LIT["D_phi_const"][1])

    # --- S4: korelasyon tensörü — BAĞIMSIZ kapalı-forma karşı ------------------
    # Sarazin & Roddier 1990 (DIMM) diferansiyel görüntü hareketi varyansları:
    #   sigma^2_boyuna = 2 lam^2 r0^{-5/3} (0.179 D^{-1/3} - 0.0968 s^{-1/3})
    #   sigma^2_enine  = 2 lam^2 r0^{-5/3} (0.179 D^{-1/3} - 0.1450 s^{-1/3})
    # Buradan normalize korelasyon (s >> D):
    #   rho_boyuna(s) = (0.0968/0.179) (D/s)^{1/3} = 0.541 (D/s)^{1/3}
    #   rho_enine(s)  = (0.1450/0.179) (D/s)^{1/3} = 0.810 (D/s)^{1/3}
    # DİKKAT: enine korelasyon boyunadan BÜYÜKTÜR (oran 0.0968/0.145 = 2/3).
    # Bu, gradyan kovaryansının C_boy = (1/2)D''(s), C_en = (1/2)D'(s)/s
    # yapısından gelir ve sezgiye aykırıdır — kural buna göre yazılmıştır.
    sD = np.array([1.0, 2.0, 4.0, 6.0, 8.0])
    rl = np.asarray(corr.rho_long(sD), dtype=float)
    rt = np.asarray(corr.rho_tran(sD), dtype=float)
    dl = 0.541 * sD ** (-1 / 3.0)
    dt = 0.810 * sD ** (-1 / 3.0)
    err = float(np.max(np.abs(np.concatenate([rl / dl, rt / dt]) - 1.0)))
    r0_ok = abs(float(corr.rho_long(np.array([0.0]))[0]) - 1) < 1e-6 and \
            abs(float(corr.rho_tran(np.array([0.0]))[0]) - 1) < 1e-6
    mono = not (np.any(np.diff(rl) > 1e-3) or np.any(np.diff(rt) > 1e-3))
    lvl = OK if (err < 0.15 and r0_ok and mono) else (RED if err > 0.30 or not r0_ok
                                                     or not mono else WARN)
    rep.add("S4-korelasyon-tensoru", "S", lvl,
            "Eğim korelasyon tensörü BAĞIMSIZ DIMM kapalı-formuyla uyuşuyor mu?",
            "s/D=" + ", ".join(f"{s:g}:({a:.3f}|{b:.3f})" for s, a, b in zip(sD, rl, rt)) +
            f"  — en büyük bağıl sapma {100*err:.1f}%; rho(0)=1 {r0_ok}; tekdüze {mono}",
            "rho_boy = 0.541 (D/s)^{1/3}, rho_en = 0.810 (D/s)^{1/3}, sapma < %15",
            "Sarazin & Roddier 1990 A&A 227:294 (DIMM); Fried 1975 JOSA 65:1427",
            "" if lvl == OK else "turb/correlation.py tensor_blocks denetlenmeli.")

    # --- S5: enine/boyuna oranı 2/3'e gidiyor mu? -----------------------------
    ratio_lt = float(np.median(rl / rt))
    ok5 = abs(ratio_lt - 2 / 3) < 0.10
    rep.add("S5-anizotropi-orani", "S", OK if ok5 else WARN,
            "Eğim korelasyonunun anizotropisi (boyuna/enine) doğru mu? "
            "Bu oran, köşegen-olmayan Sigma'nın VAROLUŞ GEREKÇESİDİR",
            f"rho_boyuna/rho_enine medyanı = {ratio_lt:.3f}",
            "2/3 = 0.667 ± 0.10 (Kolmogorov, s >> D)",
            "Sarazin & Roddier 1990: 0.0968/0.1450 = 2/3",
            "" if ok5 else
            "Oran 1'e yakınsa tensör izotroptur ve C3 katkısının anizotropi "
            "argümanı çöker; 0'a yakınsa kovaryans tekilleşir.")

    # --- S6/S7: bulanıklık çekirdeğinin GÖRÜNTÜ PİKSELİNDEKİ ölçeği -------------
    check_psf_scale(rep, quick=quick)


# =============================================================================
# (S6/S7) PSF ÖLÇEĞİ — filtrelemede kullanılan çekirdeğin piksel ölçeği
# =============================================================================
#: S7 kapalı formu (kısa-poz, tilt çıkarılmış ortalama OTF):
#:   OTF_SE(nu) = T0(nu) * exp{ -3.44 (D/r0)^{5/3} nu^{5/3} [1 - alpha nu^{1/3}] }
#:   T0(nu) = (2/pi) [acos(nu) - nu sqrt(1 - nu^2)],  nu = lam f / D in [0, 1]
#: alpha = 1 (yakın alan; yalnız faz bozulması — simülatörün ince ekranı faz-tek).
#: Kaynak: D. L. Fried, "Optical resolution through a randomly inhomogeneous medium
#: for very long and very short exposures", JOSA 56(10):1372-1379, 1966 (Denk. 7.x
#: kısa-poz MTF); aynı bağıntı Goodman, Statistical Optics (2. bs., 2015) §8.7 ve
#: Roggemann & Welsh, Imaging Through Turbulence (1996) §3.4'te. PSF, OTF'nin
#: Hankel dönüşümüdür: PSF(t) ~ int_0^1 OTF(nu) J0(2 pi nu t) nu dnu, t [lam/D].
FRIED_SE_ALPHA = 1.0


def fried_psf_profile(D_over_r0: float, short: bool = True,
                      alpha: float = FRIED_SE_ALPHA, t_max: float = 40.0):
    """Fried 1966 ortalama OTF'sinin Hankel dönüşümü: (t [lam/D], PSF(t)).

    Normalizasyon analitiktir: 2 pi int_0^inf PSF(t) t dt = OTF(0) = 1."""
    from scipy.special import j0
    nu = np.linspace(0.0, 1.0, 4001)
    T0 = (2.0 / np.pi) * (np.arccos(nu) - nu * np.sqrt(np.clip(1.0 - nu ** 2, 0, None)))
    a = 3.44 * (float(D_over_r0) ** (5.0 / 3.0)) * nu ** (5.0 / 3.0)
    if short:
        a = a * (1.0 - alpha * nu ** (1.0 / 3.0))
    otf = T0 * np.exp(-a)
    t = np.linspace(0.0, t_max, 8001)
    integrand = otf[None, :] * j0(2.0 * np.pi * nu[None, :] * t[:, None]) * nu[None, :]
    prof = 2.0 * np.pi * np.trapz(integrand, nu, axis=1)
    return t, prof


def fried_psf_fwhm_lamD(D_over_r0: float, short: bool = True,
                        alpha: float = FRIED_SE_ALPHA) -> float:
    """Fried 1966 ortalama PSF'inin FWHM'i [lam/D]."""
    t, prof = fried_psf_profile(D_over_r0, short, alpha)
    half = 0.5 * prof[0]
    i = int(np.argmax(prof < half))
    if i <= 0:
        return float("nan")
    t1, t2, p1, p2 = t[i - 1], t[i], prof[i - 1], prof[i]
    return float(2.0 * (t1 + (half - p1) * (t2 - t1) / (p2 - p1)))


def fried_psf_ee50_lamD(D_over_r0: float, short: bool = True,
                        alpha: float = FRIED_SE_ALPHA) -> float:
    """Fried 1966 ortalama PSF'inin %50 çevrelenmiş-enerji ÇAPI [lam/D]."""
    t, prof = fried_psf_profile(D_over_r0, short, alpha)
    ee = np.concatenate([[0.0], np.cumsum(0.5 * (prof[1:] * t[1:] + prof[:-1] * t[:-1])
                                          * np.diff(t))]) * 2.0 * np.pi
    return float(2.0 * np.interp(0.5, ee, t))


def _psf_fwhm(K: np.ndarray, d: float, upsample: int = 1) -> float:
    """Azimut-ortalamalı profilden FWHM [d birimi x d]. upsample>1: kübik büyütme
    (piksel ızgarasında örneklenmiş dar çekirdekler için)."""
    if upsample > 1:
        from scipy.ndimage import zoom
        K = np.clip(zoom(K, upsample, order=3), 0.0, None)
        d = d / upsample
    n0, n1 = K.shape
    yy, xx = np.mgrid[0:n0, 0:n1]
    cy, cx = (K * yy).sum() / K.sum(), (K * xx).sum() / K.sum()
    r = np.hypot(yy - cy, xx - cx).ravel()
    v = K.ravel()
    nb = int(np.ceil(r.max())) + 1
    b = np.floor(r + 0.5).astype(int)
    prof = np.bincount(b, weights=v, minlength=nb) / np.maximum(
        np.bincount(b, minlength=nb), 1)
    peak = float(K.max())
    half = 0.5 * peak
    rc = np.arange(nb, dtype=float)
    prof[0] = peak
    i = int(np.argmax(prof < half))
    if i <= 0:
        return float("nan")
    r1, r2, p1, p2 = rc[i - 1], rc[i], prof[i - 1], prof[i]
    return float(2.0 * (r1 + (half - p1) * (r2 - r1) / (p2 - p1)) * d)


def _psf_ee50(K: np.ndarray, d: float) -> float:
    """%50 çevrelenmiş-enerji çapı [d birimi x d] (ağırlık merkezi etrafında)."""
    n0, n1 = K.shape
    yy, xx = np.mgrid[0:n0, 0:n1]
    w = K / K.sum()
    cy, cx = (w * yy).sum(), (w * xx).sum()
    r = np.hypot(yy - cy, xx - cx).ravel()
    o = np.argsort(r)
    c = np.cumsum(w.ravel()[o])
    return float(2.0 * np.interp(0.5, c, r[o]) * d)


def _psf_rms(K: np.ndarray, d: float) -> float:
    """Eksen başına RMS genişlik sqrt(<x^2 + y^2>/2) [d birimi x d]."""
    n0, n1 = K.shape
    yy, xx = np.mgrid[0:n0, 0:n1]
    w = K / K.sum()
    cy, cx = (w * yy).sum(), (w * xx).sum()
    return float(np.sqrt(0.5 * (w * ((yy - cy) ** 2 + (xx - cx) ** 2)).sum()) * d)


def psf_width_table(levels=(0, 1, 3, 5, 8, 12), n_seeds: int = 2, n_side: int = 4,
                    oversample: int = 8) -> list:
    """Hattın gerçek ayarlarıyla (TurbParams varsayılanları) eski/yeni PSF genişlikleri.

    Her D/r0 için n_seeds x n_side^2 açıklık konumunda (3D aralıklı) kısa-poz PSF
    üretilir; ağırlık merkezine göre hizalı oldukları için ortalamaları Fried'in
    tilt çıkarılmış ortalama PSF'ine karşılık gelir. Eski kip: uygulanan çekirdek
    (bir örnek = bir piksel). Yeni kip: FWHM 1/oversample piksel ızgarasında, RMS
    filtrelemede kullanılan (oversample=1) çekirdekte ölçülür."""
    from .simulator import TurbParams
    from .fields import SpeckleBlurField
    from .optics import short_exposure_blur_sigma_rad
    P = TurbParams()
    D, lam, ppr = P.aperture_D, P.wavelength, P.px_per_rad
    fdelta = D / P.fine_px_per_D
    lamD_px = lam / D * ppr
    rows = []
    for dr in levels:
        acc = {"leg": None, "new1": None, "new8": None}
        n = 0
        for k in range(n_seeds):
            r0 = D / max(float(dr), 1.0)
            psf_px = min(int(P.psf_max_px), _odd_i(int(6 * max(1.0, float(dr))) + 5))
            mk = lambda **kw: SpeckleBlurField(P.fine_N, fdelta, r0, D, lam,  # noqa: E731
                                               np.random.default_rng(1000 + k),
                                               pupil_px=P.fine_px_per_D,
                                               psf_px=psf_px, **kw)
            f_leg = mk()
            f_new = mk(px_per_rad=ppr, max_px=P.psf_max_px)
            if float(dr) == 0.0:
                f_leg.phz[:] = 0.0
                f_new.phz[:] = 0.0
            for i in range(n_side):
                for j in range(n_side):
                    X, Y = (i - n_side / 2) * 3 * D, (j - n_side / 2) * 3 * D
                    ks = {"leg": f_leg.psf_at(X, Y), "new1": f_new.psf_at(X, Y),
                          "new8": f_new.psf_at(X, Y, oversample=oversample)}
                    for kk, v in ks.items():
                        acc[kk] = v if acc[kk] is None else acc[kk] + v
                    n += 1
        m = {kk: v / n for kk, v in acc.items()}
        th_se = fried_psf_fwhm_lamD(dr, short=True) * lamD_px
        th_le = fried_psf_fwhm_lamD(dr, short=False) * lamD_px
        rows.append({
            "D_over_r0": float(dr), "n_psf": n,
            "legacy_fwhm_px": _psf_fwhm(m["leg"], 1.0, upsample=8),
            "legacy_ee50_px": _psf_ee50(m["leg"], 1.0),
            "new_ee50_px": _psf_ee50(m["new8"], 1.0 / oversample),
            "fried_se_ee50_px": fried_psf_ee50_lamD(dr, short=True) * lamD_px,
            "legacy_rms_px": _psf_rms(m["leg"], 1.0),
            "legacy_kernel_px": int(m["leg"].shape[0]),
            "new_fwhm_px": _psf_fwhm(m["new8"], 1.0 / oversample),
            "new_rms_px": _psf_rms(m["new1"], 1.0),
            "new_rms_os_px": _psf_rms(m["new8"], 1.0 / oversample),
            "new_kernel_px": int(m["new1"].shape[0]),
            "new_kernel_sum": float(m["new1"].sum()),
            "fried_se_fwhm_px": th_se, "fried_le_fwhm_px": th_le,
            "gauss_se_fwhm_px": (2.3548 * short_exposure_blur_sigma_rad(r0, lam, D) * ppr
                                 if dr > 0 else float("nan")),
            "lamD_px": lamD_px, "sample_px": f_new.sample_px,
        })
    return rows


def _odd_i(n: int) -> int:
    n = max(int(n), 1)
    return n if n % 2 == 1 else n + 1


def check_psf_scale(rep: "Report", quick: bool = True, rows: list | None = None):
    """S6 kırınım FWHM'i ve S7 kısa-poz genişliğinin D/r0 ölçeklemesi."""
    from .simulator import TurbParams
    P = TurbParams()
    if rows is None:
        rows = psf_width_table(n_seeds=1 if quick else 3, n_side=4 if quick else 6)
    r0w = [r for r in rows if r["D_over_r0"] == 0.0]
    lamD = rows[0]["lamD_px"]
    if r0w:
        fw = r0w[0]["new_fwhm_px"]
        ratio = fw / lamD
        sums_ok = all(abs(r["new_kernel_sum"] - 1.0) < 1e-5 for r in rows)
        ok6 = abs(ratio - 1.0) <= 0.15 and sums_ok
        leg_ratio = r0w[0]["legacy_fwhm_px"] / lamD
        rep.add("S6-kirinim-fwhm", "S", OK if ok6 else RED,
                "Kırınım sınırlı çekirdeğin FWHM'i görüntü pikselinde lam/D*f/p mi?",
                f"FWHM = {fw:.3f} px; lam/D*f/p = {lamD:.3f} px; oran = {ratio:.3f} "
                f"(Airy 1.029); çekirdek toplamı 1 = {sums_ok}; "
                f"psf_legacy_scale={P.psf_legacy_scale}; eski kip oranı = {leg_ratio:.2f}",
                "oran 1.00 ± 0.15 ve çekirdek toplamı 1",
                "Born & Wolf, Principles of Optics §8.5.2 (Airy FWHM = 1.029 lam/D)",
                "" if ok6 else
                "PSF örnek boyu piksele çevrilmiyor (turb/fields.py SpeckleBlurField "
                "_resampled / simulator psf_legacy_scale). Bulanıklık kanalı ölçeksiz.")
    # S7 genişlik ölçütü: %50 çevrelenmiş-enerji çapı (EE50). Ortalama kısa-poz PSF
    # D/r0 >= 5'te kırınım çekirdeği + hale biçimindedir; FWHM bu yapıda çekirdeğin
    # yarı-maksimumu geçip geçmemesine bağlı olarak sıçrar ve kararsızdır. EE50 enerji
    # temellidir; FWHM oranı S7b'de ayrıca bildirilir.
    bad, parts, fparts, fbad = [], [], [], []
    for r in rows:
        if r["D_over_r0"] <= 0:
            continue
        q = r["new_ee50_px"] / r["fried_se_ee50_px"]
        parts.append(f"{r['D_over_r0']:g}:{r['new_ee50_px']:.2f}/{r['fried_se_ee50_px']:.2f}"
                     f"={q:.2f}")
        if not (abs(q - 1.0) <= 0.20):
            bad.append(r["D_over_r0"])
        qf = r["new_fwhm_px"] / r["fried_se_fwhm_px"]
        fparts.append(f"{r['D_over_r0']:g}:{qf:.2f}")
        if not (abs(qf - 1.0) <= 0.20):
            fbad.append(r["D_over_r0"])
    ok7 = not bad and bool(parts)
    rep.add("S7-kisa-poz-olcekleme", "S", OK if ok7 else RED,
            "Tilt çıkarılmış kısa-poz PSF genişliği (EE50 çapı) D/r0 ile Fried kısa-poz "
            "bağıntısına göre ölçekleniyor mu?",
            "D/r0: EE50 simülasyon/teori [px] = " + "; ".join(parts),
            "her D/r0'da oran 1.00 ± 0.20 (alpha = 1, yakın alan)",
            "Fried 1966 JOSA 56:1372 (kısa-poz OTF); Goodman Statistical Optics §8.7",
            "" if ok7 else
            f"Sapan D/r0: {bad}. İnce ekran örneklemesi (fine_px_per_D) veya FFT alanı "
            "(lam/delta) kanatları kırpıyor olabilir; turb/fields.py denetlenmeli.")
    rep.add("S7b-kisa-poz-fwhm", "S", OK if not fbad else WARN,
            "Aynı karşılaştırma FWHM ile (bilgi amaçlı; çekirdek+hale yapısında kararsız)",
            "D/r0: FWHM simülasyon/teori = " + "; ".join(fparts),
            "1.00 ± 0.20",
            "Fried 1966 JOSA 56:1372",
            "" if not fbad else
            f"FWHM oranı {fbad} düzeylerinde bandın dışında: simülasyonun ortalama "
            "kısa-poz PSF çekirdeği Fried'in kısa-poz yaklaşımından dar (EE50 uyumlu).")
    return rows


# =============================================================================
# Kayıt yükleme yardımcıları
# =============================================================================
def load(exp: str) -> list:
    p = os.path.join(EVID, exp, "records.jsonl")
    if not os.path.exists(p):
        return []
    out = []
    for line in open(p):
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return out


def _cellkey(r) -> tuple:
    c = r["cell"]
    return (r["key"], r["window"])


MIN_PATH_M = 1.0      # bir pencerenin bilgi taşıması için asgari GT yol uzunluğu


# =============================================================================
# (I) İSTATİSTİKSEL DÜRÜSTLÜK — deneyden bağımsız, her kayıt kümesine uygulanır
# =============================================================================
def check_integrity(rep: Report, recs: list, exp: str):
    if not recs:
        return

    # --- I0: KOD KÖKENİ — kayıtlar tek bir yazılım sürümünden mi geliyor? --------
    # Bu kural en başta gelir: köken karışmışsa aşağıdaki HİÇBİR bulgu anlamlı değildir.
    fps = collections.Counter(r.get("code_fp") or "damgasız" for r in recs)
    try:
        from .runner import code_fingerprint
        cur = code_fingerprint()
    except BaseException:
        cur = None
    stale = {k: v for k, v in fps.items() if cur and k not in (cur, "damgasız")}
    n_unstamped = fps.get("damgasız", 0)
    known = {k for k in fps if k != "damgasız"}
    # DAMGASIZ = "kökeni bilinmiyor", "farklı sürüm" DEĞİLDİR. Damgalama özelliği
    # eklenmeden önce üretilmiş kayıtlar böyledir; bunları tek başına KIRMIZI saymak
    # yanlış alarm üretir. KIRMIZI yalnız iki AYRI bilinen sürüm karışmışsa ya da
    # bilinen damga güncel koddan farklıysa verilir.
    if len(known) > 1 or stale:
        lvl = RED
    elif n_unstamped:
        lvl = WARN
    else:
        lvl = OK
    rep.add("I0-kod-kokeni", "I", lvl,
            "Bütün kayıtlar AYNI yazılım sürümüyle mi üretildi?",
            f"parmak izi dağılımı: {dict(fps)}"
            + (f"; güncel kod = {cur}" if cur else ""),
            "tek parmak izi ve güncel kodla aynı",
            "Tekrarlanabilirlik: Sandve 2013 PLoS CB 'Ten Simple Rules for "
            "Reproducible Computational Research' (kural 6: hangi sürümün hangi "
            "sonucu ürettiği kayıt altına alınır)",
            "" if lvl == OK else
            "İki farklı SLAM sürümünün sonuçları AYNI tabloda birleşiyor. Sürümler "
            "arası fark, yöntemin etkisi sanılır — bu kayıt kümesi kullanılamaz. "
            "Eski damgalı kayıtları arşive taşıyıp deneyi güncel kodla yeniden koş.")

    # --- I1: dejenere (hareketsiz) pencereler ---------------------------------
    ok_recs = [r for r in recs if r.get("m_valid")]
    if ok_recs:
        pl = np.array([r.get("m_path_length_gt") or 0.0 for r in ok_recs])
        frac = float(np.mean(pl < MIN_PATH_M))
        by_w = collections.Counter(
            (r["cell"]["sequence"], r["window"]) for r, p in zip(ok_recs, pl)
            if p < MIN_PATH_M)
        rep.add("I1-dejenere-pencere", "I",
                OK if frac < 0.02 else (RED if frac > 0.10 else WARN),
                "Pencerelerin GT yol uzunluğu anlamlı mı? (hareketsiz pencere ATE'yi seyreltir)",
                f"kayıtların %{100*frac:.1f}'i < {MIN_PATH_M} m yol; "
                f"suçlular: {dict(by_w) if by_w else '—'}",
                f"tüm pencereler > {MIN_PATH_M} m",
                "Zhang & Scaramuzza 2018 IROS (yörünge değerlendirme protokolü): "
                "ATE yalnız yeterli hareket içeren bölütlerde anlamlıdır",
                "" if frac < 0.02 else
                "Bu pencereler analizden ÇIKARILMALI (turb/analyze.py'de yol filtresi) "
                "ve makalede pencere seçim ölçütü açıkça yazılmalı.")

    # --- I2: kol başına başarısızlık oranı / hayatta-kalma yanlılığı -----------
    tot = collections.Counter(r["arm"] for r in recs)
    fail = collections.Counter(r["arm"] for r in recs if r.get("status") != "ok")
    rates = {a: fail[a] / tot[a] for a in tot if tot[a]}
    if rates:
        worst = max(rates, key=lambda a: rates[a])
        spread = max(rates.values()) - min(rates.values())
        lvl = OK if spread < 0.05 else (RED if spread > 0.15 else WARN)
        rep.add("I2-hayatta-kalma-yanliligi", "I", lvl,
                "Kollar aynı hücrelerde mi hayatta kalıyor? (ıraksayan hücreler düşerse "
                "kalan ortalama YANLI olur)",
                "; ".join(f"{a}: %{100*rates[a]:.0f} ({fail[a]}/{tot[a]})"
                          for a in sorted(rates)),
                "kollar arası başarısızlık farkı < %5",
                "Little & Rubin 2019 (eksik veri mekanizması: MNAR); "
                "Sturm 2012 IROS (SLAM değerlendirmesinde başarısız koşumların raporlanması)",
                "" if lvl == OK else
                f"'{worst}' kolu diğerlerinden sık ıraksıyor. (a) ıraksama oranı ATE "
                "tablosunun YANINDA birincil sonuç olarak raporlanmalı; (b) eşleştirilmiş "
                "testler yalnız BÜTÜN kolların başarılı olduğu hücre KESİŞİMİNDE "
                "koşulmalı; (c) aksi hâlde iddia geçersizdir.")

    # --- I3: eşleştirme dengesizliği ------------------------------------------
    cells_by_arm = collections.defaultdict(set)
    for r in recs:
        if r.get("status") == "ok" and r.get("m_valid"):
            cells_by_arm[r["arm"]].add(_cellkey(r))
    if len(cells_by_arm) > 1:
        common = set.intersection(*cells_by_arm.values())
        sizes = {a: len(s) for a, s in cells_by_arm.items()}
        mn = min(sizes.values())
        lost = 1 - len(common) / max(max(sizes.values()), 1)
        rep.add("I3-eslestirme-dengesi", "I",
                OK if lost < 0.05 else (RED if lost > 0.30 else WARN),
                "Eşleştirilmiş test için ortak hücre kümesi yeterli mi?",
                f"kol başına hücre: {sizes}; ortak kesişim: {len(common)} "
                f"(en büyük kolun %{100*(1-lost):.0f}'i)",
                "kesişim, en büyük kolun >%95'i",
                "Demsar 2006 JMLR (eşleştirilmiş karşılaştırma protokolü)",
                "" if lost < 0.05 else
                "Eşleştirilmiş istatistikler SADECE kesişim üzerinde koşulmalı; "
                "kol başına farklı n ile bildirilen etki büyüklükleri karşılaştırılamaz.")

    # --- I4: NEES örnek sayısı asimetrisi -------------------------------------
    nrec = [r for r in recs if r.get("n_valid")]
    if nrec:
        n_by_arm = collections.Counter(r["arm"] for r in nrec)
        if len(n_by_arm) > 1:
            mn, mx = min(n_by_arm.values()), max(n_by_arm.values())
            lvl = OK if mn >= 0.9 * mx else (RED if mn < 0.5 * mx else WARN)
            rep.add("I4-nees-orneklem", "I", lvl,
                    "NEES pencereleri kollar arasında eşit mi?",
                    f"{dict(n_by_arm)}",
                    "kollar arasında ±%10",
                    "Bar-Shalom 2001 §5.4: NEES testi ancak AYNI koşullardaki "
                    "bağımsız koşumlar üzerinde geçerlidir",
                    "" if lvl == OK else
                    "Az örnekli kolun pencereleri neden düşüyor? (cov dump yok / "
                    "GT örtüşmesi yok / ıraksama). Seçilim yanlılığı riski — "
                    "NEES da ortak hücre kesişiminde hesaplanmalı.")


# =============================================================================
# (D) DAVRANIŞ — literatürün beklediği niteliksel örüntüler
# =============================================================================
def _agg(recs, keyfn, valfn):
    d = collections.defaultdict(list)
    for r in recs:
        v = valfn(r)
        if v is not None and np.isfinite(v):
            d[keyfn(r)].append(v)
    return d


def check_behaviour(rep: Report, recs: list, exp: str):
    good = [r for r in recs if r.get("status") == "ok" and r.get("m_valid")
            and (r.get("m_path_length_gt") or 0) >= MIN_PATH_M]
    if not good:
        return

    # --- D1: temiz temel çizgi yayımlanmış EuRoC aralığında mı? ---------------
    clean = [r for r in good if r["cell"]["D_over_r0"] == 0 and r["arm"] in ("A0",)
             and r["cell"]["dataset"] == "euroc"]
    if clean:
        ate = np.array([r["m_ate_se3"] for r in clean])
        med = float(np.median(ate))
        lo, hi = LIT["euroc_ate_band"][0]
        ok = _band(med, lo, hi)
        rep.add("D1-temiz-temel-cizgi", "D", OK if ok else RED,
                "Türbülanssız temel çizgi, yayımlanmış stereo-VIO başarımıyla uyumlu mu? "
                "(uyumsuzsa 'zayıf temel çizgiyi düzelttiler' eleştirisi haklı olur)",
                f"MH temiz ATE medyanı = {med:.3f} m  (n={len(ate)} pencere, "
                f"%25-%75: {np.percentile(ate,25):.3f}-{np.percentile(ate,75):.3f})",
                f"{lo}-{hi} m arası (pencere-bazlı olduğu için alt sınır gevşek)",
                LIT["euroc_ate_band"][1],
                "" if ok else
                "Temel çizgi literatürün dışında. Yüksekse hat ayarı bozuk (iddia "
                "geçersiz); çok düşükse pencereler dejenere (bkz. I1).")

    # --- D2: bozulma-hata haritası tekdüze mi? --------------------------------
    a0 = [r for r in good if r["arm"] == "A0"]
    by_dr = _agg(a0, lambda r: r["cell"]["D_over_r0"], lambda r: r.get("m_ate_se3"))
    if len(by_dr) >= 3:
        xs = sorted(by_dr)
        ys = [float(np.median(by_dr[x])) for x in xs]
        # Spearman (bağımsız uygulama — scipy'siz)
        rx = np.argsort(np.argsort(xs)); ry = np.argsort(np.argsort(ys))
        rho = float(np.corrcoef(rx, ry)[0, 1])
        rep.add("D2-monotonluk", "D", OK if rho > 0.7 else (RED if rho < 0.3 else WARN),
                "ATE, türbülans şiddetiyle (D/r0) tekdüze artıyor mu?",
                f"Spearman rho = {rho:.2f};  " +
                ", ".join(f"D/r0={x:g}:{y:.2f}m" for x, y in zip(xs, ys)),
                "rho > 0.7 (tekdüze artan)",
                "Fried 1966 (kısa poz görüntüleme); Roggemann & Welsh 1996 §3 — "
                "görüntü bozulması D/r0 ile tekdüze artar",
                "" if rho > 0.7 else
                "Tekdüze olmayan harita ya simülatör hatası ya da ön-uçta doyum "
                "(yüksek D/r0'da zaten hiç eşleşme kalmıyor) demektir; ayırt edilmeli.")

    # --- D3: türbülans gerçekten zarar veriyor mu? ----------------------------
    if len(by_dr) >= 2 and 0.0 in by_dr:
        hi_k = max(by_dr)
        if hi_k > 0:
            r_hi = float(np.median(by_dr[hi_k])); r_0 = float(np.median(by_dr[0.0]))
            gain = r_hi / max(r_0, 1e-9)
            rep.add("D3-bozulma-etkisi", "D", OK if gain > 1.5 else RED,
                    "En yüksek şiddette bozulma ölçülebilir mi?",
                    f"ATE(D/r0={hi_k:g}) / ATE(temiz) = {gain:.2f}",
                    "> 1.5",
                    "Roggemann & Welsh 1996; Zhu & Milanfar 2013 TPAMI (türbülanslı "
                    "görüntülerde öznitelik konum hatası)",
                    "" if gain > 1.5 else
                    "Bozulma hattı etkilemiyor: enjeksiyon kancası (stereo_front_end.py) "
                    "gerçekten çalışıyor mu? E0 bit-özdeşlik testini tersine sına.")

    # --- D4: ön-uç bozulmayla zayıflıyor mu? ----------------------------------
    fkeys = ("n_matches", "n_tracks", "mean_matches", "n_features")
    fs = None
    for k in fkeys:
        d = _agg(a0, lambda r: r["cell"]["D_over_r0"],
                 lambda r: (r.get("frontend_stats") or {}).get(k))
        if len(d) >= 3:
            fs = (k, d); break
    if fs:
        k, d = fs
        xs = sorted(d); ys = [float(np.median(d[x])) for x in xs]
        drop = ys[-1] / max(ys[0], 1e-9)
        rep.add("D4-onuc-zayiflamasi", "D", OK if drop < 0.95 else WARN,
                "Türbülans ön-uçta ölçülebilir bir zayıflama yaratıyor mu?",
                f"{k}: temizde {ys[0]:.0f} -> D/r0={xs[-1]:g}'de {ys[-1]:.0f} "
                f"(oran {drop:.2f})",
                "< 0.95 (eşleşme sayısı azalmalı)",
                "Zhu & Milanfar 2013; Anantrasirichai 2013 (türbülansta öznitelik "
                "kararlılığının düşmesi)")


# =============================================================================
# (I/D) KOL KARŞILAŞTIRMASI — asıl iddianın dürüstlük denetimi
# =============================================================================
def check_arms(rep: Report, recs: list, exp: str):
    good = [r for r in recs if r.get("status") == "ok" and r.get("m_valid")
            and (r.get("m_path_length_gt") or 0) >= MIN_PATH_M]
    arms = sorted({r["arm"] for r in good})
    if len(arms) < 2:
        return

    # ortak hücre kesişimi (I3'ün zorunlu kıldığı doğru taban)
    per = collections.defaultdict(dict)
    for r in good:
        per[r["arm"]][_cellkey(r)] = r
    common = set.intersection(*[set(per[a]) for a in arms])
    if len(common) < 5:
        rep.add("A0-yetersiz-kesisim", "I", WARN,
                "Kollar arası ortak hücre sayısı", f"{len(common)}", ">= 5",
                "eşleştirilmiş tasarım gereği")
        return
    ck = sorted(common)
    ate = {a: np.array([per[a][c].get("m_ate_se3", np.nan) for c in ck]) for a in arms}
    nees = {a: np.array([per[a][c].get("n_nees_pos_median", np.nan) for c in ck])
            for a in arms}
    cov = {a: np.array([per[a][c].get("n_coverage_pos", np.nan) for c in ck])
           for a in arms}

    base = "A0" if "A0" in arms else arms[0]

    # --- I5: tutarlılık doğruluk pahasına mı satın alındı? --------------------
    for a in arms:
        if a == base:
            continue
        if np.all(~np.isfinite(nees[a])) or np.all(~np.isfinite(nees[base])):
            continue
        m = np.isfinite(nees[a]) & np.isfinite(nees[base]) & \
            np.isfinite(ate[a]) & np.isfinite(ate[base])
        if m.sum() < 5:
            continue
        nr = float(np.median(nees[a][m]) / max(np.median(nees[base][m]), 1e-12))
        da = float(np.median(ate[a][m] - ate[base][m]))
        rel = da / max(float(np.median(ate[base][m])), 1e-9)
        if nr < 0.5 and rel > 0.15:
            rep.add(f"I5-tutarlilik-bedeli-{a}", "I", WARN,
                    f"'{a}' kolunun tutarlılık kazancı doğruluk feda edilerek mi alındı?",
                    f"NEES x{nr:.3f} ama ATE %{100*rel:+.0f} ({da:+.3f} m)",
                    "NEES iyileşirken ATE bozulmamalı (veya bozulma açıkça bildirilmeli)",
                    "Bar-Shalom 2001 §5.4; Julier & Uhlmann 2001 — kovaryansı "
                    "büyütmek tutarlılığı HER ZAMAN 'satın alabilir'; katkı ancak "
                    "doğruluk korunuyorsa anlamlıdır",
                    "Makalede bu ödünleşim gizlenmeden, A1 (skaler şişirme) kolu "
                    "ile YAN YANA raporlanmalı: A1 aynı bedeli ödüyorsa yapı "
                    "gerçekten kazandırıyor demektir.")

    # --- I6: literatürde görülmemiş etki büyüklüğü ----------------------------
    lim = LIT["consistency_gain_plausible"][0]
    for a in arms:
        if a == base:
            continue
        m = np.isfinite(nees[a]) & np.isfinite(nees[base])
        if m.sum() < 5:
            continue
        nr = float(np.median(nees[base][m]) / max(np.median(nees[a][m]), 1e-12))
        if nr > lim:
            rep.add(f"I6-asiri-etki-{a}", "I", RED,
                    f"'{a}' kolunun bildirdiği kazanç literatürdeki büyüklüklerin "
                    "çok üstünde — kabul edilmeden önce açıklanmalı",
                    f"NEES iyileşmesi {nr:.0f} kat (temel çizgi {base})",
                    f"< {lim:.0f} kat; tipik bildirilen aralık 2-20 kat",
                    LIT["consistency_gain_plausible"][1],
                    "İki olasılık var ve AYIRT EDİLMELİ: (a) temel çizgi patolojik "
                    "biçimde aşırı-güvenli (NEES~1000) — bu durumda kazanç gerçek ama "
                    "iddia 'kalibre temel çizgiye (A0C) göre' yeniden ifade edilmeli; "
                    "(b) NEES hesabında hata/ölçek karışıklığı var. A0C kolu bunun "
                    "için koşuluyor; sonuç gelene dek bu oran makalede kullanılmamalı.")

    # --- D5: kapsama oranı aşırı-tutucu mu? ----------------------------------
    for a in arms:
        c = cov[a][np.isfinite(cov[a])]
        if len(c) < 5:
            continue
        mc = float(np.median(c))
        if mc > 0.995:
            rep.add(f"D5-asiri-tutucu-{a}", "D", WARN,
                    f"'{a}' kolu aşırı-tutucu mu? (kovaryans şişirilmiş olabilir)",
                    f"kapsama = {mc:.3f}", "~0.95 (nominal)",
                    LIT["coverage_nominal"][1],
                    "Kapsama 1.0'a dayanmışsa kestirici bilgi kaybediyordur; "
                    "ATE ile birlikte değerlendirilmeli.")

    # --- I8: YAPI kazanci teorik siniri (C3'un ana iddiasi) --------------------
    # A1 ve A5 AYNI kosegeni tasir; fark yalnizca kosegen-disi yapidir. Bu yuzden
    # kazanc iki YONDE de sinanabilir: teorik olarak >= 1 olmali (yapi hicbir zaman
    # zarar vermez) ve teorik ust siniri asmamali.
    if "A1" in nees and "A5" in nees:
        m = np.isfinite(nees["A1"]) & np.isfinite(nees["A5"])
        if m.sum() >= 5:
            g = float(np.median(nees["A1"][m]) / max(np.median(nees["A5"][m]), 1e-12))
            bnd = LIT["structure_gain_bound"][0]
            if g < 1.0:
                lvl, note = RED, (
                    "Yapi, skaler sisirmeden KOTU cikti. Kapali-form turetim bunun "
                    "IMKANSIZ oldugunu soyluyor (ayni kosegen + korelasyon => daha az "
                    "bilgi => daha buyuk posterior => daha kucuk NEES). Bu bir BULGU "
                    "degil, bir KUSURDUR; kol karsilastirmasi yapilmadan once "
                    "turb/backend.py kovaryans yolu denetlenmelidir. (Yasanmis: "
                    "corr_model 'physical'e sabitlenmis, nugget iki kez sayilmis ve "
                    "stereo capraz-kamera korelasyonu 1.0 alinmisti.)")
            elif g > bnd:
                lvl, note = RED, (
                    "Kazanc teorik ust sinirin ustunde. Kosegenlerin GERCEKTEN esit "
                    "oldugu dogrulanmali (turb.preflight P5 ve backend'deki "
                    "olcek/nugget ayrimi); esit degilse fark yapidan degil OLCEKTEN "
                    "geliyordur ve iddia yanlis kurulmustur.")
            else:
                lvl, note = OK, ""
            rep.add("I8-yapi-kazanci", "I", lvl,
                    "Yapili kovaryansin skaler sisirmeye gore kazanci teorik "
                    "aralikta mi? (C3'un ANA iddiasi)",
                    f"A1/A5 NEES orani = {g:.2f} kat",
                    f"[1.0, {bnd:.0f}] kat",
                    LIT["structure_gain_bound"][1], note)

    # --- I7: temiz veride tutarlılık taban çizgisi ----------------------------
    cl = [r for r in good if r["cell"]["D_over_r0"] == 0 and r["arm"] == base
          and np.isfinite(r.get("n_nees_pos_median") or np.nan)]
    if cl:
        v = float(np.median([r["n_nees_pos_median"] for r in cl]))
        lo, hi = LIT["nees_band"][0]
        rep.add("I7-temiz-nees", "I", OK if _band(v, lo, hi) else WARN,
                "Temel çizgi TÜRBÜLANSSIZ veride tutarlı mı? (değilse gözlenen "
                "tutarsızlık türbülansa değil hattın kendisine aittir)",
                f"{base} temiz NEES medyanı = {v:.1f}",
                f"[{lo}, {hi}] bandı içinde",
                LIT["nees_band"][1],
                "" if _band(v, lo, hi) else
                "Bu bir BULGU olarak yazılmalı ('hat ATE için ayarlanmış, tutarlılık "
                "için değil') ve A0C kalibre kolu ana karşılaştırma tabanı yapılmalı. "
                "Aksi hâlde 'kötü temel çizgiyi düzelttiniz' eleştirisi geçerlidir.")


# =============================================================================
def run_experiment_checks(exp: str, rep: Report):
    recs = load(exp)
    if not recs:
        rep.add(f"{exp}-veri", "I", WARN, f"{exp} kayıtları", "records.jsonl yok/boş",
                "kayıt bekleniyor", "—")
        return
    check_integrity(rep, recs, exp)
    check_behaviour(rep, recs, exp)
    check_arms(rep, recs, exp)


def write_report(rep: Report, title: str, path: str = ALERT):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    import datetime
    lines = [f"# Tahlil raporu — {title}",
             f"Üretim: {datetime.datetime.now():%Y-%m-%d %H:%M} · "
             f"genel durum: **{rep.worst}**", ""]
    n_ack = sum(1 for f in rep.findings if f.level == ACK)
    if n_ack:
        lines += [f"> {n_ack} bulgu **ele alınmış** olarak kayıtlı "
                  f"(`{ACK_FILE}`); gerekçeleriyle aşağıda duruyorlar — "
                  f"kapatılmaları onları görünmez yapmaz.", ""]
    for grp, name in (("S", "S — Simülatör fiziği (yazılım doğruluğu)"),
                      ("D", "D — Davranış (literatürle niteliksel uyum)"),
                      ("I", "I — İddia dürüstlüğü (istatistiksel geçerlilik)")):
        fs = [f for f in rep.findings if f.group == grp]
        if not fs:
            continue
        lines.append(f"## {name}\n")
        for f in sorted(fs, key=lambda f: {RED: 0, WARN: 1, ACK: 2, OK: 3}[f.level]):
            lines.append(f.md())
    txt = "\n".join(lines)
    with open(path, "w") as fh:
        fh.write(txt)
    return txt


def main():
    ap = argparse.ArgumentParser(description="Sonuç tahlili — fizik/literatür denetimi")
    ap.add_argument("--exp", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--physics", action="store_true", help="yalnız simülatör fiziği")
    ap.add_argument("--full", action="store_true", help="fizikte daha çok gerçekleme")
    ap.add_argument("--out", default=None)
    ap.add_argument("--psf_table", action="store_true",
                    help="yalnız S6/S7 ve eski/yeni PSF genişlik tablosu")
    a = ap.parse_args()

    if a.psf_table:
        rep = Report(acks=load_acks())
        rows = psf_width_table(n_seeds=3 if a.full else 2, n_side=6 if a.full else 4)
        check_psf_scale(rep, rows=rows)
        print("| D/r0 | eski FWHM | eski RMS | eski EE50 | eski çekirdek | yeni FWHM | "
              "yeni RMS | yeni EE50 | yeni çekirdek | Fried SE FWHM | Fried SE EE50 | "
              "Fried LE FWHM |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in rows:
            print(f"| {r['D_over_r0']:g} | {r['legacy_fwhm_px']:.2f} | {r['legacy_rms_px']:.2f} | "
                  f"{r['legacy_ee50_px']:.2f} | {r['legacy_kernel_px']} | "
                  f"{r['new_fwhm_px']:.2f} | {r['new_rms_px']:.2f} | {r['new_ee50_px']:.2f} | "
                  f"{r['new_kernel_px']} | {r['fried_se_fwhm_px']:.2f} | "
                  f"{r['fried_se_ee50_px']:.2f} | {r['fried_le_fwhm_px']:.2f} |")
        print(json.dumps(rows, indent=1))
        txt = write_report(rep, "PSF ölçeği", a.out or os.path.join(EVID, "sanity_psf.md"))
        print(txt)
        sys.exit(rep.exit_code())

    rep = Report(acks=load_acks())
    title = a.exp or ("tüm deneyler" if a.all else "fizik")
    if a.physics or a.all:
        check_physics(rep, quick=not a.full)
    if a.all:
        for e in sorted(os.listdir(EVID)):
            if os.path.exists(os.path.join(EVID, e, "records.jsonl")):
                run_experiment_checks(e, rep)
    elif a.exp:
        run_experiment_checks(a.exp, rep)

    out = a.out or (ALERT if rep.worst != OK else
                    os.path.join(EVID, f"sanity_{title.replace(' ', '_')}.md"))
    txt = write_report(rep, title, out)
    print(txt)
    print(f"\n>>> genel: {rep.worst}  ->  {out}")
    sys.exit(rep.exit_code())


if __name__ == "__main__":
    main()
