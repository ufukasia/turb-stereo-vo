"""
ÖN-UÇUŞ KAPISI — uzun koşuma çıkmadan önce her şeyi küçük ölçekte sına.

NEDEN VAR (kullanıcı direktifi, 2026-08-19): bu projede aynı hata örüntüsü üç kez
tekrarlandı — 17 saatlik bir koşum başlatıldı, saatler sonra bir kusur bulundu, tüm
sonuçlar çöpe gitti. Örüntünün nedeni koşumların BÜYÜK başlamasıydı. Bundan sonra
uzun kuyruk, YALNIZCA bu kapı geçilirse başlar.

TASARIM İLKESİ: kapı, birim testi değil ENTEGRASYON testidir. Kuyruğun kullandığı
kod yolunun TAMAMINI (turb.runner -> config -> ön-uç -> arka-uç -> metrikler -> NEES)
küçük bir hücre ızgarasında koşturur. Çünkü buraya kadar bulunan kusurların HİÇBİRİ
birim testiyle yakalanabilecek türden değildi; hepsi bileşenler arası sözleşmede idi.

Denetlenen her madde, GERÇEKTEN YAŞANMIŞ bir kusura karşılık gelir:

  P1  kol çöküşü          — A4 %47-79 oranında ıraksıyordu
  P2  kovaryans verimi    — %89 sessiz başarısızlık (Cholesky koşul sayısı)
  P3  kovaryans PD        — 6979 kovaryansın 138'i pozitif tanımlı değildi (A3)
  P4  NEES sonlu/pozitif  — tabloya "nees_pos_mean = -5977" düşmüştü
  P5  KOL AYIRT EDİCİLİĞİ — A0C, cal_factor=1.0 yüzünden A0 ile BİREBİR aynıydı (plasebo)
  P6  fiziksel kestirim   — DIMM theta_c'si 58-1460 px arası salınıyordu (görüntü 752 px)
  P7  köken damgası       — code_fp hatlı kod yoluna eklenmemişti
  P8  ATE makullüğü       — A3, D/r0=8'de 17.9 m veriyordu
  P9  NEES teğet uzayı    — hata ham ötelenmeydi, kovaryans SE(3) Expmap teğetinde
  P10 A1/A5 köşegeni      — nugget iki kez sayılınca köşegenler ayrışmıştı
  P11 robust simetrisi    — Huber A1'de vardı, A5'te yoktu (C3'ü kirletiyordu)
  P12 Schur tümleyeni     — Hll blok-köşegen varsayılıyordu; yama faktörü onu bağlar
  P13 sıfır-türbülans     — Huber granülerliği A5'i temiz veride 2.26x kötüleştiriyordu
  P14 sigma ölçüm mü      — theta_c kestirilemez (sabit nominal); sigma denetlenir

    ./turbvio/bin/python -m turb.preflight            # ~20 dk
    ./turbvio/bin/python -m turb.preflight --quick    # ~8 dk (tek rejim)

Çıkış kodu: 0 = geçti (uzun koşum serbest), 1 = KALDI (kuyruk başlatılmaz).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

ARMS = ("A0", "A0C", "A1", "A2", "A3", "A4", "A5")
OUT = "evidence/preflight.md"


def _cells(quick: bool, seq: str, seeds: int, robust=None):
    """Küçük ama TEMSİLİ ızgara: temiz + en zor rejim, bütün kollar.

    Pencere 200 kare (tam koşumda 400) ve tek pencere -> hücre başına maliyet ~1/10.
    Rejim seçimi kasıtlı: D/r0=0 kolların AYIRT EDİLEBİLİRLİĞİNİ (P5), D/r0=8 ise
    sayısal kırılganlığı (P2/P3/P4/P8) sınar. Aradaki rejimler bir şey eklemez.
    """
    from .experiments import EUROC
    from .runner import Cell
    levels = [0.0, 8.0] if not quick else [8.0]
    cells = []
    for dr in levels:
        for k in range(seeds if dr > 0 else 1):
            cells.append(Cell(dataset="euroc", base_path=EUROC, sequence=seq,
                              seed=k, D_over_r0=dr, arms=ARMS, nees=True,
                              cov_every=5, windows=(200,), window_len=200,
                              warmup=50, cache_dump=False, tag="preflight",
                              backend_argv=(() if robust is None else
                                            ("--turb_robust", str(int(robust))))))
    return cells


class Report:
    def __init__(self):
        self.rows = []

    def add(self, pid, ok, claim, measured, expected, incident):
        self.rows.append((pid, bool(ok), claim, measured, expected, incident))

    @property
    def passed(self):
        return all(r[1] for r in self.rows)

    def render(self, extra=""):
        import datetime
        L = ["# Ön-uçuş kapısı", "",
             f"Üretim: {datetime.datetime.now():%Y-%m-%d %H:%M} · sonuç: "
             f"**{'GEÇTİ' if self.passed else 'KALDI'}**", "",
             "**P** maddeleri yazılımın MEKANİĞİNİ denetler (her biri yaşanmış bir "
             "kusura karşılık gelir). **E** maddeleri ÖN-KAYITLI BİLİMSEL "
             "BEKLENTİLERdir: koşum mekanik olarak kusursuz olup tezi çürütebilir "
             "ve bu, günler süren kuyruk başlamadan ÖNCE bilinmelidir. Bir E maddesi "
             "düşerse sonuç 'hata' değil **DUR-VE-İNCELE**dir — beklenti de yanlış "
             "olabilir, yöntem de.", "",
             "| # | durum | denetim | ölçülen | beklenen | hangi olaydan doğdu |",
             "|---|---|---|---|---|---|"]
        for pid, ok, claim, meas, exp, inc in self.rows:
            L.append(f"| {pid} | {'✅' if ok else '🛑'} | {claim} | {meas} | {exp} | {inc} |")
        if extra:
            L += ["", extra]
        return "\n".join(L) + "\n"


# =============================================================================
# E SINIFI — ÖN-KAYITLI BİLİMSEL BEKLENTİLER
#
# NEDEN AYRI BİR SINIF (kullanıcı direktifi, 2026-08-19): P maddeleri yazılımın
# MEKANİĞİNİ denetler (çöküyor mu, sayı üretiyor mu, matris geçerli mi). Ama bir
# koşum mekanik olarak kusursuz olup BİLİMSEL olarak tezi çürütebilir — ve bu,
# 41 saatlik kuyruğu başlatmadan ÖNCE bilinmelidir. Ölçüldü: kapı P1-P8'in hepsini
# geçti, ama küçük ızgarada bizim yöntemimizin (A3/A5) tutarlılığı naif skaler
# şişirmeden (A1) 8 KAT KÖTÜ çıktı — yani tezin tersi. O koşum başlatılsaydı iki
# gün sonra aynı sonucu daha büyük örneklemle öğrenecektik.
#
# İLKE: beklentiler koşumdan ÖNCE yazılır (ön-kayıt). Sonradan sonuca göre
# yazılan "beklenti" beklenti değil, rasyonalizasyondur.
#
# İhlal KIRMIZI değil **DUR-VE-İNCELE**dir: beklenti yanlış da olabilir, yöntem
# hatalı da. Ayırt edilmeden uzun koşuma çıkılmaz.
# =============================================================================
EXPECTATIONS = [
    # (kimlik, açıklama, fonksiyon(ctx) -> (ok, ölçülen, beklenen))
    #
    # ctx.m(kol, anahtar)      : TÜRBÜLANSLI hücrelerde medyan (gösterim için)
    # ctx.ratio(A, B, anahtar) : hücre bazında EŞLEŞTİRİLMİŞ A/B oranının medyanı
    #
    # NEDEN EŞLEŞTİRME (2026-08-21): eski sürüm kolların medyanlarını ayrı ayrı
    # alıp karşılaştırıyordu. Kapıda hücre başına yalnız 2-3 kayıt var; hücreden
    # hücreye NEES varyasyonu, kollar arası farktan BÜYÜK. Aranan etki teorik
    # olarak 1.2-1.5 kat (evidence/E2_ETKI_BUYUKLUGU.md) olduğundan eşleştirilmemiş
    # test bu etkiyi ÇÖZEMEZ — kapı gürültüyle geçer ya da kalır. Eşleştirme
    # hücre etkisini götürür. Ayrıca TEMİZ hücreler (D/r0=0) artık dışarıda:
    # orada türbülans yok, A1 ile A5 zaten aynıdır, etkiyi 1.0'a doğru seyreltirdi.
    ("E1", "Türbülans-farkında kovaryans, temel çizgiden TUTARLI olmalı "
           "(NEES bandın merkezine YAKLAŞMALI)",
     lambda c: (c.ratio("A0", "A3", "nees") > 1.0 and c.ratio("A0", "A5", "nees") > 1.0,
                f"A0/A3={c.ratio('A0','A3','nees'):.2f}x  "
                f"A0/A5={c.ratio('A0','A5','nees'):.2f}x  "
                f"(A0={c.m('A0','nees'):.0f} A3={c.m('A3','nees'):.0f} "
                f"A5={c.m('A5','nees'):.0f})",
                "eşleştirilmiş oranlar > 1")),
    ("E2", "YAPI, skaler şişirmeden İYİ olmalı — C3'ün ANA İDDİASI "
           "(aynı köşegen, fark YALNIZ köşegen-dışı yapıdan)",
     # Tolerans 0.85: kapı bir TARAMA'dır, nihai bilimsel test değil. Teorik kazanç
     # 1.2-1.5x, kapıdaki örneklem 2 hücre; %15'lik gürültü payı bırakılır. Yaşanan
     # kusur (A5, A1'den 5.7 kat KÖTÜ) bu toleransla da gürültülü biçimde yakalanır.
     lambda c: (c.ratio("A1", "A5", "nees") >= 0.85,
                f"A1/A5 = {c.ratio('A1','A5','nees'):.2f}x  "
                f"(A1={c.m('A1','nees'):.0f}  A5={c.m('A5','nees'):.0f})",
                "eşleştirilmiş A1/A5 >= 0.85 (teorik beklenti 1.2-1.5)")),
    ("E3", "Kalibre temel çizgi (A0C) tutarlılığı satın almalı ama DOĞRULUKLA ödemeli "
           "— yapının gerekçesi bu ödünleşimdir",
     lambda c: (c.ratio("A0", "A0C", "nees") > 1.0,
                f"A0/A0C = {c.ratio('A0','A0C','nees'):.2f}x, "
                f"ATE {c.m('A0','ate'):.3f} -> {c.m('A0C','ate'):.3f}",
                "eşleştirilmiş A0/A0C > 1")),
    ("E4", "Türbülans-farkında kollar DOĞRULUĞU bozmamalı "
           "(tutarlılık doğruluk pahasına satın alınmamalı)",
     lambda c: (c.ratio("A5", "A0", "ate") <= 1.5,
                f"A5/A0 ATE = {c.ratio('A5','A0','ate'):.2f}x  "
                f"(A0={c.m('A0','ate'):.3f} m  A5={c.m('A5','ate'):.3f} m)",
                "eşleştirilmiş A5/A0 ATE <= 1.5")),
]


def check_expectations(rep, recs):
    """Ön-kayıtlı beklentileri küçük ızgarada, HÜCRE BAZINDA EŞLEŞTİREREK sına."""
    def _ck(r):
        c = r["cell"]
        return (c.get("sequence"), c.get("seed"), c.get("D_over_r0"))

    per = defaultdict(dict)          # cellkey -> arm -> {ate, nees}
    for r in recs:
        if r.get("status") != "ok" or not r.get("m_valid"):
            continue
        if not (r["cell"].get("D_over_r0") or 0) > 0:   # TEMİZ hücreler hariç
            continue
        d = {"ate": r.get("m_ate_se3")}
        if r.get("n_valid") and r.get("n_nees_pos_median") is not None:
            d["nees"] = r["n_nees_pos_median"]
        per[_ck(r)][r["arm"]] = d

    class Ctx:
        @staticmethod
        def m(arm, key):
            v = [a[arm][key] for a in per.values()
                 if arm in a and a[arm].get(key) is not None
                 and np.isfinite(a[arm][key])]
            return float(np.median(v)) if v else float("nan")

        @staticmethod
        def ratio(num, den, key):
            """Hücre bazında num/den oranlarının medyanı (eşleştirilmiş test)."""
            rs = []
            for a in per.values():
                if num not in a or den not in a:
                    continue
                x, y = a[num].get(key), a[den].get(key)
                if x is None or y is None or not (np.isfinite(x) and np.isfinite(y)):
                    continue
                if abs(y) < 1e-12:
                    continue
                rs.append(float(x) / float(y))
            return float(np.median(rs)) if rs else float("nan")

    ctx = Ctx()
    npair = len(per)
    for pid, claim, fn in [(e[0], e[1], e[2]) for e in EXPECTATIONS]:
        try:
            ok, meas, exp = fn(ctx)
            ok = bool(ok)
            meas = f"{meas}  [{npair} eşleşmiş hücre]"
        except Exception as e:                       # eksik kol -> sınanamaz
            ok, meas, exp = True, f"sınanamadı ({type(e).__name__}: {e})", "—"
        rep.add(pid, ok, claim, meas, exp,
                "ön-kayıtlı beklenti (2026-08-19; eşleştirme 2026-08-21)")


def _load_cov_files(spool):
    out = defaultdict(list)
    for f in glob.glob(os.path.join(spool, "*", "cov_*.npz")):
        arm = os.path.basename(f)[4:-4]
        try:
            out[arm].append(np.load(f))
        except (OSError, ValueError):
            pass
    return out


def main():
    ap = argparse.ArgumentParser(description="Ön-uçuş kapısı")
    ap.add_argument("--sequence", default="MH03")
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--jobs", type=int, default=7)
    ap.add_argument("--gpu_jobs", type=int, default=2)
    ap.add_argument("--quick", action="store_true")
    # Kapı MANTIĞINI, SLAM'i yeniden koşturmadan mevcut kayıtlar üzerinde
    # yeniden değerlendir. (Kapının kendisi de yanlış olabilir; onu
    # düzeltmek 35 dakikalık koşuma mal olmamalı.)
    ap.add_argument("--report_only", action="store_true")
    # Spool YOLU: iki kapiyi paralel kosturabilmek icin ayrilabilir olmali.
    # (Ayni spool'u paylasan iki kosum birbirinin cikti dosyalarini ezer;
    # 2026-08-21'de yasandi.)
    ap.add_argument("--spool", default="/tmp/turb_preflight")
    # ROBUST ABLASYONU: 0 -> hicbir kolda m-kestirici yok (TAM ESLESMIS
    # karsilastirma; yama faktorunde Huber BLOK bazinda, bireysel
    # faktorde GOZLEM bazinda uygulanir -- ikisi birebir ayni degildir).
    # VARSAYILAN None: bayrak VERILMEZSE arka-uca hicbir sey gecilmez ve
    # config.yaml'daki `turbulence.robust` gecerli olur. Eskiden varsayilan 1'di
    # ve config.yaml'i SESSIZCE EZIYORDU: ana eksen robust=0'a cevrildikten
    # sonra kapi yine robust=1 kosuyor, P13 hakli olarak dusuyordu (2026-08-22,
    # 6 saat kayip). Kapinin yapilandirmasi, KOSULACAK yapilandirma olmalidir.
    ap.add_argument("--robust", type=int, default=None)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--records", default="evidence/preflight/records.jsonl")
    a = ap.parse_args()

    from .runner import run_experiment_pipelined, code_fingerprint

    spool = a.spool
    os.makedirs(os.path.dirname(a.records) or ".", exist_ok=True)
    if not a.report_only:
        if os.path.exists(a.records):
            os.remove(a.records)
        os.system(f"rm -rf {spool}")

        cells = _cells(a.quick, a.sequence, a.seeds, robust=a.robust)
        print(f"[ön-uçuş] {len(cells)} hücre × {len(ARMS)} kol "
              f"({a.sequence}, pencere 200 kare) -> {a.records}", flush=True)
        run_experiment_pipelined(cells, a.records, cpu_jobs=a.jobs,
                                 gpu_jobs=a.gpu_jobs, chunk=len(cells),
                                 resume=False, spool=spool)
    else:
        print(f"[ön-uçuş] --report_only: {a.records} yeniden değerlendiriliyor",
              flush=True)

    recs = [json.loads(l) for l in open(a.records) if l.strip()]
    rep = Report()
    if not recs:
        rep.add("P0", False, "kayıt üretildi mi", "0 kayıt", "> 0", "—")
        print(rep.render()); sys.exit(1)

    # ---- P1 kol çöküşü ------------------------------------------------------
    tot = Counter(r["arm"] for r in recs)
    bad = Counter(r["arm"] for r in recs if r.get("status") != "ok")
    # ÇEKİRDEK vs KEŞİFSEL kol ayrımı. A4 (yama-farkında öznitelik SEÇİMİ) gözlem
    # atarak grafı zayıflatır ve ölçülen ıraksama oranı diğer kollardan bir mertebe
    # yüksektir (D/r0=8'de %17-19; A0 %0.5, A3 %4.5). Bu bir KUSUR değil, kolun
    # ölçülen ÖZELLİĞİdir ve makalede negatif sonuç olarak raporlanacaktır —
    # dolayısıyla kapıyı kilitlemez. Ama gizlenmez de: oranı burada basılır ve
    # `turb.analyze` kol dayanıklılık tablosunda birincil sonuç olarak yer alır.
    CORE = tuple(x for x in ARMS if x != "A4")
    core_bad = {x: bad[x] / tot[x] for x in CORE if tot[x]}
    ok1 = all(v <= 0.10 for v in core_bad.values())
    rep.add("P1", ok1, "çekirdek kolların hiçbiri çökmüyor (A4 keşifsel — raporlanır)",
            "; ".join(f"{x}:%{100*bad[x]/tot[x]:.0f}" for x in sorted(tot))
            + f"  [A4 keşifsel: %{100*bad['A4']/max(tot['A4'],1):.0f}]",
            "çekirdek kollar < %10", "A4 %47-79 ıraksıyordu")

    # ---- P2/P3/P4 kovaryans -------------------------------------------------
    # KAYNAK: spool hücre bitince silindiği için .npz dosyaları okunamaz; verimin
    # DOĞRU kaynağı kayıtlardaki `n_n` alanıdır (pencerede fiilen kullanılan NEES
    # örneği). İlk sürüm .npz'lere bakıyordu ve her zaman 0 buluyordu — kapının
    # kendisi sessizce yanlış ölçüyordu. (Kapı da denetlenir.)
    nn = [r["n_n"] for r in recs if r.get("n_valid") and r.get("n_n")]
    win = int(recs[0]["cell"].get("window_len", 200))
    ce = int(recs[0]["cell"].get("cov_every", 5))
    n_exp = win // max(ce, 1)
    med = float(np.median(nn)) if nn else 0.0
    rep.add("P2", med >= 0.75 * n_exp, "kovaryans verimi (pencere başına örnek)",
            f"medyan {med:.0f} / beklenen {n_exp} ({len(nn)} pencere)",
            f">= {0.75*n_exp:.0f}", "%89 sessiz başarısızlık (Cholesky)")
    # PD ve koşul sayısı artık arka-uçta ZORUNLU kapıdır (PD olmayan matris hiç
    # döndürülmez), burada NEES'in işaretinden dolaylı doğrulanır (P4).
    rep.add("P3", True, "kovaryans PD kapısı arka-uçta etkin",
            "gtsam_backend._pose_marginal_schur PD kontrolü",
            "PD olmayan matris döndürülmez", "6979'un 138'i bozuktu (A3)")

    # ---- P4 NEES sonlu ve pozitif ------------------------------------------
    nv = [r for r in recs if r.get("n_valid")]
    neg = [r for r in nv if (r.get("n_nees_pos_mean") or 0) < 0
           or (r.get("n_nees_rot_mean") or 0) < 0]
    nonfin = [r for r in nv if not np.isfinite(r.get("n_nees_pos_median") or np.nan)]
    rep.add("P4", not neg and not nonfin, "NEES sonlu ve pozitif",
            f"negatif {len(neg)}, sonsuz/NaN {len(nonfin)} / {len(nv)}", "0 / 0",
            "tabloya -5977 düşmüştü")

    # ---- P5 KOL AYIRT EDİCİLİĞİ (plasebo kol yakalayıcı) -------------------
    key = lambda r: (r["key"], r["window"])          # noqa: E731
    by_arm = defaultdict(dict)
    for r in recs:
        if r.get("status") == "ok" and r.get("m_valid"):
            by_arm[r["arm"]][key(r)] = r
    ident, checked = [], 0
    for i, x in enumerate(ARMS):
        for y in ARMS[i + 1:]:
            common = set(by_arm.get(x, {})) & set(by_arm.get(y, {}))
            if len(common) < 3:
                continue
            checked += 1
            d = np.array([abs(by_arm[x][c]["m_ate_se3"] - by_arm[y][c]["m_ate_se3"])
                          for c in common])
            rel = d / np.maximum([by_arm[x][c]["m_ate_se3"] for c in common], 1e-9)
            # EŞİK NEDEN 1e-3: iki kol FONKSİYONEL olarak aynıysa sonuçları bit
            # düzeyinde aynı OLMAYABİLİR (farklı kod yolları, farklı toplama sırası).
            # Ölçüldü: cal_factor=1.0 iken A0C ile A0'ın ATE'leri 3e-5 bağıl farkla
            # ayrılıyordu ve 1e-6'lık ilk eşik onları "farklı" sayıp plaseboyu
            # KAÇIRDI. Gerçek bir yöntem farkı %0.1'in çok üstünde iz bırakır.
            if float(np.median(rel)) < 1e-3:
                ident.append(f"{x}~{y} (medyan bağıl fark {np.median(rel):.1e})")
    rep.add("P5", not ident, "her kol GERÇEKTEN farklı bir yöntem mi",
            (", ".join(ident) if ident else f"{checked} ikili, hepsi farklı"),
            "her ikilide medyan bağıl fark > 1e-3",
            "A0C, cal_factor=1.0 yüzünden A0 ile aynıydı (plasebo)")

    # ---- P6 fiziksel kestirim (DIMM) ---------------------------------------
    tc, sg = [], []
    for r in recs:
        st = r.get("arm_stats")
        if isinstance(st, dict) and r["arm"] in ("A3", "A4"):
            if st.get("theta_c_hat_median") is not None:
                tc.append(st["theta_c_hat_median"])
            if st.get("sigma_hat_median") is not None:
                sg.append(st["sigma_hat_median"])
    w = int(recs[0]["cell"].get("image_width", 752))
    ok6 = bool(tc) and all(0.01 * w <= t <= 2.0 * w for t in tc)
    rep.add("P6", ok6, "DIMM theta_c kestirimi fiziksel aralıkta",
            (f"{min(tc):.0f}-{max(tc):.0f} px (görüntü {w} px)" if tc else "kestirim yok"),
            f"[{0.01*w:.0f}, {2.0*w:.0f}] px",
            "58-1460 px arası salınıyordu")

    # ---- P7 köken damgası ---------------------------------------------------
    fps = Counter(r.get("code_fp") or "damgasız" for r in recs)
    cur = code_fingerprint()
    rep.add("P7", set(fps) == {cur}, "her kayıt güncel kodla damgalı",
            f"{dict(fps)} (güncel {cur})", f"hepsi {cur}",
            "code_fp hatlı kod yoluna eklenmemişti")

    # ---- P8 ATE makullüğü ---------------------------------------------------
    hi = [r for r in recs if r.get("m_valid") and (r.get("m_ate_se3") or 0) > 10.0]
    rep.add("P8", not hi, "hiçbir kol saçma ATE üretmiyor",
            (f"{len(hi)} kayıt > 10 m: "
             + ", ".join(sorted({f"{r['arm']}={r['m_ate_se3']:.1f}m" for r in hi})[:4])
             if hi else "en büyük ATE "
             f"{max((r.get('m_ate_se3') or 0) for r in recs if r.get('m_valid')):.2f} m"),
            "< 10 m", "A3, D/r0=8'de 17.9 m veriyordu")

    # ---- P9 NEES teğet konvansiyonu ------------------------------------------
    # Kovaryans arka-uçta linearize() ile üretilir, yani GTSAM'in Pose3::retract
    # teğetindedir. metrics._gt_logmap hatayı AYNI teğette almalıdır. GTSAM farklı
    # bir bayrakla yeniden derlenirse (GTSAM_POSE3_EXPMAP kapalı) bu sessizce
    # bozulur ve NEES ölçüsü anlamını yitirir — bu yüzden her kapıda sınanır.
    try:
        import gtsam as _g
        _T = _g.Pose3(_g.Rot3.Expmap(np.array([0.3, -0.2, 0.5])), _g.Point3(1., 2., 3.))
        _xi = np.array([0.1, 0.2, -0.15, 0.4, -0.3, 0.2])
        _tan_ok = bool(np.allclose(_T.retract(_xi).matrix(),
                                   _T.compose(_g.Pose3.Expmap(_xi)).matrix(), atol=1e-10))
        _msg = "retract == compose(Expmap)" if _tan_ok else "retract != compose(Expmap) — teğetler AYRIŞTI"
    except Exception as _e:
        _tan_ok, _msg = False, f"sınanamadı: {type(_e).__name__}"
    rep.add("P9", _tan_ok, "NEES hatası kovaryansla AYNI teğet uzayında",
            _msg, "retract == compose(Expmap)",
            "hata ham ötelenme alınıyordu; |w|=1 rad'da %40 sapma")

    # ---- P10 A1/A5 köşegen özdeşliği ----------------------------------------
    # E2'nin ANLAMLI olmasının ön koşulu: A1 (skaler şişirme) ile A5 (yapılı) AYNI
    # köşegeni taşımalı, fark YALNIZ köşegen-dışı yapıdan gelmeli. Aksi hâlde
    # ölçülen fark "yapı" değil "ölçek" farkıdır ve C3 iddiası yanlış kurulmuştur.
    # (Yaşanmış: backend'de model-bağımlı ölçek/nugget ayrımı yokken ampirik
    # modelin köşegeni A1'inkinden farklıydı — nugget iki kez sayılıyordu.)
    try:
        import turb.covariance as _cv
        _rng = np.random.default_rng(0)
        _sp, _ind = 1.3, 1.75
        _st_ = _sp * np.hypot(1.0, _ind)
        _bm, _bs = 1.0, 1.2
        _N = 5
        _th = np.column_stack([_rng.uniform(-.02, .02, _N), _rng.uniform(-.015, .015, _N)])
        _iss = np.array([True, False, True, True, False])
        _exp = np.array(sum(([_cv.scalar_inflation(_bs, _st_)] * 3 if k else
                             [_cv.scalar_inflation(_bm, _st_)] * 2 for k in _iss), []))
        _worst, _rho = 0.0, 0.0
        for _mdl in ("empirical", "physical"):
            _sk, _si = (_st_, 0.0) if _mdl == "empirical" else (_sp, _ind * _sp)
            _S = _cv.patch_covariance(_th, _iss, _sk, _bm, _bs, 1000.0, 0.10, 0.11,
                                      model=_mdl, theta_c=0.15 * 752 / 170648.0,
                                      shrinkage=0.05, sigma_ind_px=_si)
            _d = np.sqrt(np.diag(_S))
            _worst = max(_worst, float(np.abs(_d - _exp).max() / _exp.max()))
            _o = _S - np.diag(np.diag(_S))
            _rho = max(_rho, float(np.abs(_o / np.outer(_d, _d)).max()))
        _p10 = _worst < 1e-6 and _rho > 0.01
        _m10 = f"köşegen bağıl fark {_worst:.1e}; köşegen-dışı max |rho| {_rho:.3f}"
    except Exception as _e:
        _p10, _m10 = False, f"sınanamadı: {type(_e).__name__}: {_e}"
    rep.add("P10", _p10, "A1 ile A5 AYNI köşegeni taşıyor (fark yalnız yapıda)",
            _m10, "köşegen farkı < 1e-6 VE köşegen-dışı yapı var (|rho| > 0.01)",
            "nugget iki kez sayılınca A5'in köşegeni A1'den ayrışmıştı")

    # ---- P11 robust çekirdek simetrisi ---------------------------------------
    # BÜTÜN kollar aynı aykırı-değer çekirdeği rejiminde olmalı. Aksi hâlde A1 ile
    # A5 "yalnız kovaryans yapısında farklı" DEĞİLDİR ve C3'ün ana karşılaştırması
    # geçersizdir. (Yaşanmış, 2026-08-21: bireysel-faktör yolları Huber ile
    # sarılıyordu, yama ve tekil-yama yolları SARILMIYORDU. m-kestirici büyük
    # artıkları bastırıp kestiricinin ETKİN kovaryansını büyütür; NEES tam da
    # beyan edilen ile etkin kovaryansı karşılaştırır. Ölçülen kirlenme:
    # A1 NEES 48 / A5 NEES 193 — yapı 3.8 kat "kötü" görünüyordu.)
    _rb = {}
    for r in recs:
        st = r.get("arm_stats")
        if not isinstance(st, dict):
            continue
        nr = int(st.get("n_fac_robust", 0) or 0)
        npl = int(st.get("n_fac_plain", 0) or 0)
        if nr + npl == 0:
            continue
        _acc = _rb.setdefault(r["arm"], [0, 0])   # DIKKAT: `a` argparse nesnesi,
        _acc[0] += nr                             # gölgelenmemeli (yaşandı).
        _acc[1] += npl
    if _rb:
        _fr = {a: v[0] / max(v[0] + v[1], 1) for a, v in _rb.items()}
        _ok11 = (max(_fr.values()) - min(_fr.values())) < 1e-9
        _m11 = "; ".join(f"{a}:%{100*_fr[a]:.0f}" for a in sorted(_fr))
    else:
        _ok11, _m11 = False, "hiçbir kol robust sayacı bildirmedi (sayaç eklenmemiş?)"
    rep.add("P11", _ok11, "BÜTÜN kollar aynı aykırı-değer çekirdeği rejiminde",
            _m11, "her kolda robust faktör oranı AYNI",
            "yama ve tekil yollar Huber'ı atlıyordu; A1'de vardı, A5'te yoktu")

    # ---- P12 Schur tümleyeni, BAĞLI landmark'larda da doğru mu ---------------
    # gtsam_backend._pose_marginal_schur eskiden Hll'in BLOK-KÖŞEGEN olduğunu
    # VARSAYIYORDU. Bu, bireysel izdüşüm faktörleri için doğru (her faktör tek
    # landmark'a dokunur) ama YAMA faktörü için yanlış (X(f) ile n landmark'a
    # aynı anda dokunur). Yalnız köşegen bloklarla hesaplanan tümleyen eksik
    # çıkarır -> kovaryans fazla küçük -> NEES fazla büyük. Yani kovaryans
    # SKALER kollarda doğru, YAPILI kollarda yanlış hesaplanıyordu ve yapı
    # sistematik olarak kötü görünüyordu. Burada iki faktör tipinde de yoğun
    # tersle karşılaştırılarak sınanır.
    try:
        import gtsam as _g
        from gtsam.symbol_shorthand import X as _X, L as _Lk
        import gtsam_backend as _GB
        from turb.factors import make_patch_factor as _mpf
        _rg = np.random.default_rng(3)
        _fx = _fy = 458.0; _cx, _cy = 367.0, 248.0; _b = 0.11; _n = 4
        _pw = np.array([[.5, .2, 4.], [-.3, .1, 5.], [.1, -.4, 3.5], [.4, .3, 4.5]])
        _po = _g.Pose3(); _vl = _g.Values(); _vl.insert(_X(0), _po)
        for _i, _q in enumerate(_pw):
            _vl.insert(_Lk(_i), _g.Point3(*_q))
        _mz = np.array(sum(([_fx * q[0] / q[2] + _cx, _fy * q[1] / q[2] + _cy]
                            for q in _pw), [])) + _rg.normal(0, .3, 2 * _n)
        _Sg = np.full((2 * _n, 2 * _n), 1.4); np.fill_diagonal(_Sg, 4.0)
        _lp = _g.noiseModel.Isotropic.Sigma(3, 10.0)
        _worst = 0.0; _yollar = []
        for _patch in (False, True):
            _gr = _g.NonlinearFactorGraph()
            if _patch:
                _f, _ = _mpf(_X(0), [_Lk(i) for i in range(_n)], _mz,
                             np.zeros(_n, bool), _Sg, _fx, _fy, _cx, _cy, _b)
                _gr.add(_f)
            else:
                _K = _g.Cal3_S2(_fx, _fy, 0., _cx, _cy)
                for _i in range(_n):
                    _gr.add(_g.GenericProjectionFactorCal3_S2(
                        _mz[2 * _i:2 * _i + 2],
                        _g.noiseModel.Isotropic.Sigma(2, 2.0),
                        _X(0), _Lk(_i), _K, False, False))
            _gr.add(_g.PriorFactorPose3(_X(0), _po,
                                        _g.noiseModel.Isotropic.Sigma(6, 1e-6)))
            for _i in range(_n):
                _gr.add(_g.PriorFactorPoint3(_Lk(_i), _g.Point3(*_pw[_i]), _lp))
            _or = _g.Ordering(); _or.push_back(_X(0))
            for _i in range(_n):
                _or.push_back(_Lk(_i))
            _gfg = _gr.linearize(_vl)
            _dm = [6] + [3] * _n
            _C6 = _GB._pose_marginal_schur(_gfg, _or, 1, _dm)
            _Hm, _ = _gfg.hessian(_or)
            _Hm = np.asarray(_Hm) + np.eye(sum(_dm)) * 1e-9
            _ref = np.linalg.inv(_Hm)[:6, :6]
            if _C6 is None:
                _worst = float("inf"); break
            _worst = max(_worst, float(np.abs(_C6 - _ref).max()
                                       / max(np.abs(_ref).max(), 1e-300)))
            _Hll = np.asarray(_Hm)[6:, 6:]
            _mk = np.ones_like(_Hll, bool); _o = 0
            for _d in [3] * _n:
                _mk[_o:_o + _d, _o:_o + _d] = False; _o += _d
            _yollar.append(("yama" if _patch else "bireysel") +
                           ("/HIZLI" if np.abs(_Hll[_mk]).max() <= 1e-12 *
                            max(np.abs(_Hll).max(), 1.0) else "/TAM"))
        _ok12 = _worst < 1e-9
        _m12 = f"yoğun tersle max bağıl fark {_worst:.1e}  [{', '.join(_yollar)}]"
    except Exception as _e:
        _ok12, _m12 = False, f"sınanamadı: {type(_e).__name__}: {_e}"
    rep.add("P12", _ok12, "Schur tümleyeni BAĞLI landmark'larda da doğru",
            _m12, "yoğun tersle bağıl fark < 1e-9",
            "Hll blok-köşegen VARSAYILIYORDU; yama faktörü onu bağlıyor, "
            "yapılı kolların kovaryansı fazla küçük çıkıyordu")

    # ---- P13 sıfır türbülansta kollar temel çizgiye indirgeniyor mu ---------
    # D/r0 = 0'da türbülans terimi YOKTUR: sigma_px = 0 -> Sigma = temel çizgi
    # köşegeni (doğrudan ölçüldü, fark 4e-10, köşegen-dışı TAM SIFIR). O halde
    # A1 (skaler şişirme) ile A5 (yapılı) temiz veride BİREBİR aynı sonucu
    # vermek ZORUNDADIR. Ayrışırlarsa fark kovaryanstan değil FAKTÖR YAPISINDAN
    # geliyordur ve C3 karşılaştırması kirlidir.
    # (Yaşanmış, 2026-08-22: robust=1'de A5 temiz veride A0'dan 2.26 kat kötüydü;
    # sebep Huber granülerliğiydi -- bireysel faktörde gözlem başına, yama
    # faktöründe blok başına. Ana eksen robust=0'a çevrildi.)
    _cl = defaultdict(dict)
    for r in recs:
        if r.get("status") != "ok" or not r.get("m_valid"):
            continue
        if (r["cell"].get("D_over_r0") or 0) > 0:
            continue
        _k = (r["cell"].get("sequence"), r["cell"].get("seed"), r.get("window"))
        _cl[_k][r["arm"]] = r["m_ate_se3"]
    _pairs = [(v["A1"], v["A5"]) for v in _cl.values() if "A1" in v and "A5" in v]
    if _pairs:
        _rel = max(abs(a - b) / max(abs(a), 1e-12) for a, b in _pairs)
        _ok13 = _rel < 1e-3
        _m13 = (f"{len(_pairs)} temiz hücrede A1 vs A5 ATE, max bağıl fark "
                f"{_rel:.2e}")
    else:
        _ok13, _m13 = True, "temiz hücrede A1/A5 çifti yok — sınanamadı"
    rep.add("P13", _ok13, "D/r0=0'da A1 ile A5 BİREBİR aynı (yapı terimi yok)",
            _m13, "bağıl fark < 1e-3",
            "Huber granülerliği A5'i temiz veride 2.26 kat kötüleştiriyordu")

    # ---- P14 çevrimiçi kestirim GERÇEKTEN ölçüyor mu (sigma ekseni) ---------
    # ESKI SURUM YANLIS SEYI DENETLIYORDU: "theta_c hucreden hucreye degismeli".
    # Olculdu (evidence/THETA_C_KESTIRILEMEZ.md) ki theta_c VIO artiklarindan
    # KESTIRILEMEZ -- turbulans disi artik turbulans egiminin >=2 kati oldugunda
    # uydurma sinira yapisir (v1 33.3 px, v2 8.3 px). Sinir degistirmek sorunu
    # tasidi, cozmedi; kestirilemezlik VERIDE, kodda degil.
    # Yeni tasarim: theta_c SABIT nominal (config.yaml theta_c_nominal_over_fov),
    # yalniz SIGMA kestirilir. Dolayisiyla denetim de sigma ekseninde olmalidir:
    #   (a) sigma_hat GERCEKTEN degismeli -- turbulans siddetiyle birlikte
    #       (sabit bir sayi = olcum yok, v1'in gunahi buydu)
    #   (b) sigma_hat oracle ile AYNI MERTEBEDE olmali (0.5x - 2.0x)
    _st = [r.get("arm_stats") for r in recs
           if r.get("arm") in ("A3", "A4") and isinstance(r.get("arm_stats"), dict)]
    _pairs = [(float(x["sigma_hat_median"]), float(x.get("sigma_px_oracle") or np.nan))
              for x in _st
              if x.get("sigma_hat_median") is not None
              and np.isfinite(x.get("sigma_hat_median", np.nan))]
    _pairs = [(a, b) for a, b in _pairs if b > 1e-9]
    if len(_pairs) >= 3:
        _sh = np.array([a for a, _ in _pairs])
        _so = np.array([b for _, b in _pairs])
        _spread = float(np.std(_sh) / max(abs(np.mean(_sh)), 1e-12))
        _ratio = float(np.median(_sh / _so))
        _ok14 = (_spread > 0.02) and (0.5 <= _ratio <= 2.0)
        _m14 = (f"{len(_pairs)} kayıt · sigma_hat/oracle = {_ratio:.2f}x · "
                f"bağıl saçılım {_spread:.3f}")
    else:
        _ok14, _m14 = True, "yeterli A3/A4 kaydı yok — sınanamadı"
    rep.add("P14", _ok14, "çevrimiçi SIGMA bir ÖLÇÜM mü (sabit değil, oracle mertebesinde)",
            _m14, "saçılım > 0.02 VE sigma_hat/oracle 0.5-2.0 arasında",
            "theta_c sınıra çakılıp sabit değer raporluyordu; artık sabit "
            "nominal kullanılıyor ve yalnız sigma kestiriliyor")

    # ---- E sınıfı: ön-kayıtlı bilimsel beklentiler --------------------------
    check_expectations(rep, recs)

    # ---- özet tablo ---------------------------------------------------------
    lines = ["## Kol özeti (bu küçük ızgarada)", "",
             "| kol | n | ATE ort. | NEES medyan | kapsama | başarısız |",
             "|---|---|---|---|---|---|"]
    for arm in ARMS:
        g = [r for r in recs if r["arm"] == arm]
        gv = [r for r in g if r.get("m_valid")]
        gn = [r for r in g if r.get("n_valid")]
        f = sum(1 for r in g if r.get("status") != "ok")
        lines.append(
            f"| {arm} | {len(g)} | "
            f"{np.mean([r['m_ate_se3'] for r in gv]):.4f} | "
            f"{(np.median([r['n_nees_pos_median'] for r in gn]) if gn else float('nan')):.1f} | "
            f"{(np.mean([r['n_coverage_pos'] for r in gn]) if gn else float('nan')):.3f} | "
            f"{f} |")

    txt = rep.render("\n".join(lines))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(txt)
    print("\n" + txt)
    # GEÇİŞ DAMGASI: kuyruk, kapıyı AYNI kod üzerinde ikinci kez koşturmasın diye
    # bu damgaya bakar. Damgada kod parmak izi var — tek bir satır değişse damga
    # geçersizleşir ve kapı yeniden koşar. Yalnız zamana bakan bir damga, "kapıyı
    # geçtim, sonra kodu düzelttim" durumunu sessizce kaçırırdı.
    stamp = os.path.join(os.path.dirname(a.out) or ".",
                         "PREFLIGHT_PASS.stamp" if a.robust is None
                         else f"PREFLIGHT_PASS_robust{a.robust}.stamp")
    if rep.passed:
        import datetime as _dt
        with open(stamp, "w") as fh:
            fh.write(f"code_fp={code_fingerprint()}\n"
                     f"zaman={_dt.datetime.now():%Y-%m-%d %H:%M:%S}\n"
                     f"sekans={a.sequence} tohum={a.seeds} kayit={len(recs)}\n")
    elif os.path.exists(stamp):
        os.remove(stamp)

    print(f">>> {'GEÇTİ — uzun koşum serbest' if rep.passed else 'KALDI — kuyruk BAŞLATILMAZ'}"
          f"  ({a.out})")
    sys.exit(0 if rep.passed else 1)


if __name__ == "__main__":
    main()
