"""Marjinalizasyonsuz kayan pencere demet ayarlamasi (ORB-SLAM3 LocalBundleAdjustment karsiligi).

BatchFixedLagSmoother pencereden cikan pozu marjinalize eder; pozun gordugu ve pencerede kalan butun noktalar
yogun bir marjinal faktorle birbirine baglanir, her guncellemede marjinalizasyon icin ayrica dogrusallastirma ve
eleme yapilir. Burada cikan poz marjinalize edilmez, dusurulur:
  * pencerenin en eski pozu o anki kestiriminde sabitlenir (ayar serbestligi),
  * dusen pozun pencerede kalan noktalara olan gozlemleri, poz sabit kabul edilerek noktaya 3x3 bir oncel olarak
    eklenir. ORB-SLAM3 yerel BA'si yerel haritayi goren eski anahtar kareleri sabit tutar; bu onun dogrusallastirilmis
    karsiligidir ve noktalar arasinda yogun baglanti kurmaz.
Dogrusal sistem noktalar once (Schur tumleyeni sirasi) elenir; LM sonumu kucuk baslatilmaz (zayif kisitli noktalarda
tekil sistem nedeniyle bosa eleme olmasin diye).

Arayuz, gtsam_backend'in kullandigi kadariyla BatchFixedLagSmoother ile aynidir: update(graf, degerler, zaman),
calculateEstimate(), calculateEstimatePose3/Point3.
"""
import os
import time
import numpy as np
import gtsam
from gtsam import NonlinearFactorGraph, Values


class ZamanHaritasi(dict):
    """FixedLagSmootherKeyTimestampMap yerine: anahtar -> zaman damgasi."""

    def insert(self, kv):
        self[kv[0]] = float(kv[1])


class KayanPencere:
    def __init__(self, pencere, lm_params, sabit_ek=0, sabit_n=1, sabit_sigma=1e-4, nokta_once=True,
                 dusen_oncel=True, kf_aralik=0, yakin=10):
        self.win = float(pencere)
        self.lmp = lm_params
        self.ek = float(sabit_ek)          # pencerenin disinda sabit tutulan gozlemci pozlarin derinligi [kare]
        self.sabit_n = int(sabit_n)        # pencere icinde sabitlenen en eski poz sayisi (ek=0 iken ayar serbestligi)
        self.nokta_once = bool(nokta_once)
        self.dusen_oncel = bool(dusen_oncel)
        # anahtar kare seyreltmesi: son `yakin` kare tam kalir; daha eskilerden yalniz her `kf_aralik`'inci kare
        # pencere sonuna kadar tutulur, aradakiler dusurulur (gozlemleri noktalara oncel olarak aktarilir)
        self.kf_aralik = int(kf_aralik)
        self.yakin = float(yakin)
        self.sabit_gurultu = gtsam.noiseModel.Isotropic.Sigma(6, float(sabit_sigma))
        self.theta = Values()              # pencere + sabit bolgedeki degiskenlerin son kestirimi
        self.ts = {}                       # theta'daki anahtar -> zaman damgasi
        self.poz = set()                   # theta'daki poz anahtarlari
        self.kova = {}                     # en kucuk poz zaman damgasi -> faktorler
        self.oncel_H, self.oncel_g = {}, {}   # nokta -> dusen gozlemlerin bilgi matrisi ve bilgi-agirlikli konumu
        self.oncel_f = {}                  # nokta -> guncel oncel faktoru
        self.t = float('-inf')
        self.son = Values()                # son cozumun pencere ici kismi
        self.sure = {'kur': 0.0, 'lm': 0.0, 'it': 0, 'n': 0, 'fak': 0, 'deg': 0}

    @staticmethod
    def _poz_mu(k):
        return gtsam.Symbol(k).chr() == ord('x')

    def update(self, graph, values, ts_map):
        # BatchFixedLagSmoother gibi: faktorler ve degerler once eklenir, cozum hata verse de durumda kalir
        for k in values.keys():
            # gozlemi kalmadigi icin kestirimde gorunmeyen nokta arka uc tarafindan yeniden eklenebilir
            if self.theta.exists(k):
                self.theta.erase(k)        # update(k, Point3) Python'da Vector olarak saklaniyor
            if gtsam.Symbol(k).chr() == ord('x'):
                self.theta.insert(k, values.atPose3(k))
            else:
                self.theta.insert(k, values.atPoint3(k))
        for k, v in ts_map.items():
            self.ts[k] = max(v, self.ts.get(k, v))
            if self._poz_mu(k):
                self.poz.add(k)
        if ts_map:
            self.t = max(self.t, max(ts_map.values()))
        for k in values.keys():
            if k not in self.ts:
                self.ts[k] = self.t
                if self._poz_mu(k):
                    self.poz.add(k)
        for i in range(graph.size()):
            fa = graph.at(i)
            if fa is None:
                continue
            # kova = (faktorun dusecegi an, en eski pozu): faktor pozlarindan biri dusunce duser
            pz = [self.ts[k] for k in fa.keys() if k in self.poz]
            b = (min(self._dz(p_) for p_ in pz), min(pz)) if pz else (self.t + self.win + self.ek, self.t)
            g = self.kova.get(b)
            if g is None:
                g = self.kova[b] = NonlinearFactorGraph()
            g.add(fa)
        self._coz()

    def _dusur(self, kg, alt):
        """Dusen kovadaki gozlemleri, pozlar sabit kabul edilerek pencerede kalan noktalara 3x3 oncel olarak aktar."""
        for i in range(kg.size()):
            fa = kg.at(i)
            ks = list(fa.keys())
            hedef = [k for k in ks if k not in self.poz and self.ts.get(k, -np.inf) >= alt]
            if not hedef or not all(self.theta.exists(k) for k in ks):
                continue
            try:
                A, b = fa.linearize(self.theta).jacobian()
            except Exception:
                continue
            ofs, c = {}, 0
            for k in ks:
                d = 6 if k in self.poz else 3
                ofs[k] = (c, d)
                c += d
            for k in hedef:
                o, _ = ofs[k]
                Al = A[:, o:o + 3]
                H = Al.T @ Al
                p = np.asarray(self.theta.atPoint3(k))
                self.oncel_H[k] = self.oncel_H.get(k, 0) + H
                self.oncel_g[k] = self.oncel_g.get(k, 0) + H @ p + Al.T @ b
                self.oncel_f.pop(k, None)

    def _dz(self, p):
        """Zaman damgasi p olan pozun dusecegi an: anahtar kare pencere (+ sabit bolge) sonunda, ara kare `yakin`
        kare sonra."""
        if self.kf_aralik > 1 and int(round(p)) % self.kf_aralik != 0:
            return p + self.yakin
        return p + self.win + self.ek

    def _coz(self):
        _t0 = time.time()
        alt = self.t - self.win            # zaman damgasi alt'tan kucuk anahtar pencere disi
        sinir = alt - self.ek
        for b in [b for b in self.kova if b[0] < self.t]:
            if self.dusen_oncel:
                self._dusur(self.kova[b], alt)
            del self.kova[b]
        for k in [k for k, v in self.ts.items()
                  if (self._dz(v) < self.t if k in self.poz else v < alt)]:
            self.theta.erase(k)
            del self.ts[k]
            self.poz.discard(k)
            self.oncel_H.pop(k, None); self.oncel_g.pop(k, None); self.oncel_f.pop(k, None)

        g = NonlinearFactorGraph()
        for b in sorted(self.kova):
            kg = self.kova[b]
            if b[1] >= alt:
                g.push_back(kg)
                continue
            # sabit bolge: yalniz pencerede kalan noktalara olan gozlemler
            for i in range(kg.size()):
                fa = kg.at(i)
                if all(k in self.ts and (k in self.poz or self.ts[k] >= alt) for k in fa.keys()):
                    g.add(fa)
        _faktorlu = set(g.keyVector())
        for k, H in self.oncel_H.items():
            if k not in _faktorlu:
                continue                   # gozlemi kalmamis nokta cozulmez (theta'da durur)
            fo = self.oncel_f.get(k)
            if fo is None:
                Hr = H + 1e-4 * np.eye(3)
                gr = self.oncel_g[k] + 1e-4 * np.asarray(self.theta.atPoint3(k))
                fo = self.oncel_f[k] = gtsam.PriorFactorPoint3(
                    k, np.linalg.solve(Hr, gr), gtsam.noiseModel.Gaussian.Information(Hr))
            g.add(fo)
        anahtar = g.keyVector()
        v = Values()
        pozlar, noktalar = [], []
        for k in anahtar:
            if k in self.poz:
                v.insert(k, self.theta.atPose3(k))
                pozlar.append((self.ts[k], k))
            else:
                v.insert(k, self.theta.atPoint3(k))
                noktalar.append(k)
        pozlar.sort()
        sabit = [k for t_, k in pozlar if t_ < alt]
        if not sabit:
            sabit = [k for t_, k in pozlar[:self.sabit_n]]
        for k in sabit:
            g.add(gtsam.PriorFactorPose3(k, self.theta.atPose3(k), self.sabit_gurultu))
        lmp = self.lmp
        if self.nokta_once:
            o = gtsam.Ordering()
            for k in noktalar:
                o.push_back(k)
            for _, k in pozlar:
                o.push_back(k)
            lmp.setOrdering(o)            # parametre nesnesi yalniz bu duzlestiricide kullanilir
        _t1 = time.time()
        self.k = getattr(self, 'k', 0) + 1
        if str(self.k) == os.environ.get('PYECO_KAYAN_DOKUM_KARE'):
            # tani: pencere grafini dosyaya dok (cozucu karsilastirmasi icin)
            with open(os.environ['PYECO_KAYAN_DOKUM'], 'w') as _d:
                _d.write(g.serialize() + '\n@@@\n' + v.serialize())
        _opt = gtsam.LevenbergMarquardtOptimizer(g, v, lmp)
        r = _opt.optimize()
        _t2 = time.time()
        self.sure['kur'] += _t1 - _t0; self.sure['lm'] += _t2 - _t1; self.sure['it'] += _opt.iterations()
        self.sure['n'] += 1; self.sure['fak'] += g.size(); self.sure['deg'] += v.size()
        self.theta.update(r)
        son = Values()
        for k in anahtar:
            if self.ts[k] >= alt:
                if k in self.poz:
                    son.insert(k, r.atPose3(k))
                else:
                    son.insert(k, r.atPoint3(k))
        self.son = son

    def calculateEstimate(self):
        return self.son

    def calculateEstimatePose3(self, k):
        return self.theta.atPose3(k)

    def calculateEstimatePoint3(self, k):
        return self.theta.atPoint3(k)
