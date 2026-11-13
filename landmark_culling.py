"""
Landmark ayıklama: kare-yerel kapıların aksine ısrarcı kötü landmark'lar arka uçta kara-listeye
alınır; kötü lid'e bir daha factor eklenmez, PnP/motion-only/SBP'de kullanılmaz ve fixed-lag
penceresinden marjinalize olur. gtsam_backend ve RealtimeBackend aynı sınıfı kullanır.

Kurallar (hepsi varsayılan kapalı): cull_found_ratio (found/visible oranı), cull_chi2_strikes ve
cull_pnp_strikes (ısrarcı outlier), cull_allow_reinit (kara-listedeki lid taze stereo ile yeniden doğabilir mi).
"""
from collections import Counter

import numpy as np


class LandmarkCuller:
    """Landmark başına istatistik + kara-liste; tüm kurallar kapalıyken `active=False`."""

    def __init__(self, args, fx, fy, cx, cy, img_w=None, img_h=None):
        g = lambda k, d: getattr(args, k, d) if getattr(args, k, None) is not None else d  # noqa: E731
        self.found_ratio = float(g('cull_found_ratio', 0.0) or 0.0)
        self.found_min_visible = int(g('cull_found_min_visible', 8) or 8)
        self.chi2_strikes = int(g('cull_chi2_strikes', 0) or 0)
        self.pnp_strikes = int(g('cull_pnp_strikes', 0) or 0)
        self.allow_reinit = int(g('cull_allow_reinit', 1)) > 0
        self.fx, self.fy, self.cx, self.cy = float(fx), float(fy), float(cx), float(cy)
        # Görüntü sınırı verilmezse principal point'ten simetrik
        self.img_w = float(img_w) if img_w else 2.0 * self.cx
        self.img_h = float(img_h) if img_h else 2.0 * self.cy
        self.stats = {}          # lid -> dict(created, n_obs, n_visible, n_found, n_chi2, n_pnp)
        self.bad = set()         # kara-liste
        self.audit = Counter()
        self.active = any([self.found_ratio > 0, self.chi2_strikes > 0, self.pnp_strikes > 0])

    # ------------------------------------------------------------------ durum
    def describe(self):
        on = []
        if self.found_ratio > 0:
            on.append(f"found_ratio<{self.found_ratio} (n_vis≥{self.found_min_visible})")
        if self.chi2_strikes > 0:
            on.append(f"chi2_strikes={self.chi2_strikes}")
        if self.pnp_strikes > 0:
            on.append(f"pnp_strikes={self.pnp_strikes}")
        return ", ".join(on) + f" | reinit={'açık' if self.allow_reinit else 'kapalı'}"

    def is_bad(self, lid):
        return lid in self.bad

    def can_reinit(self, lid):
        """Kara-listedeki lid taze stereo ile yeniden yaratılabilir mi?"""
        return self.allow_reinit

    # ------------------------------------------------------------------ olaylar
    def on_init(self, lid, f):
        """Landmark grafa (yeniden) girdi → istatistik sıfırla, kara-listeden çıkar."""
        if not self.active:
            return
        self.stats[lid] = dict(created=f, n_obs=0, n_visible=0, n_found=0, n_chi2=0, n_pnp=0)
        if lid in self.bad:
            self.bad.discard(lid)
            self.audit['cull_reinit'] += 1

    def on_observation(self, lid):
        if not self.active:
            return
        s = self.stats.get(lid)
        if s is not None:
            s['n_obs'] += 1

    def on_chi2_outliers(self, lids):
        if not self.active or not lids:
            return
        for lid in lids:
            s = self.stats.get(lid)
            if s is not None:
                s['n_chi2'] += 1
                self.audit['strike_chi2'] += 1

    def on_pnp_outliers(self, lids):
        if not self.active or not lids:
            return
        for lid in lids:
            s = self.stats.get(lid)
            if s is not None:
                s['n_pnp'] += 1
                self.audit['strike_pnp'] += 1

    def visibility_pass(self, T_WC, lm_world, active_set, obs_ids):
        """Aktif landmark'ları güncel pozla yansıt: görüntü içindeyse visible++, gözlemlendiyse found++."""
        if self.found_ratio <= 0:
            return
        T_CW = np.linalg.inv(T_WC)
        for lid in active_set:
            if lid in self.bad:
                continue
            s = self.stats.get(lid)
            pw = lm_world.get(lid)
            if s is None or pw is None:
                continue
            pc = T_CW @ np.array([pw[0], pw[1], pw[2], 1.0])
            if pc[2] <= 1e-3:
                continue
            u = self.fx * pc[0] / pc[2] + self.cx
            v = self.fy * pc[1] / pc[2] + self.cy
            if 0.0 <= u < self.img_w and 0.0 <= v < self.img_h:
                s['n_visible'] += 1
                if lid in obs_ids:
                    s['n_found'] += 1

    # ------------------------------------------------------------------ karar
    def evaluate(self, f, active_set):
        """Kuralları uygula; bu karede kara-listeye alınan lid kümesini döndür.
        Çağıran: lm_world / active_set / lm_desc / pending_mv'den siler."""
        if not self.active:
            return set()
        culled = set()
        for lid in list(active_set):
            if lid in self.bad:
                continue
            s = self.stats.get(lid)
            if s is None:
                continue
            why = None
            if self.chi2_strikes > 0 and s['n_chi2'] >= self.chi2_strikes:
                why = 'chi2_strikes'
            elif self.pnp_strikes > 0 and s['n_pnp'] >= self.pnp_strikes:
                why = 'pnp_strikes'
            elif (self.found_ratio > 0 and s['n_visible'] >= self.found_min_visible
                  and s['n_found'] / float(s['n_visible']) < self.found_ratio):
                why = 'found_ratio'
            if why is not None:
                culled.add(lid)
                self.bad.add(lid)
                self.audit['culled_' + why] += 1
        if culled:
            self.audit['culled_total'] += len(culled)
        return culled

    def apply_cull(self, culled, lm_world, *sets, dicts=()):
        """Kara-listeye alınanları çağıranın defterlerinden sil (lm_world + verilen küme/sözlükler)."""
        for lid in culled:
            lm_world.pop(lid, None)
            for st in sets:
                st.discard(lid)
            for d in dicts:
                d.pop(lid, None)

    def on_window(self, resident_lms):
        """Pencereden marjinalize olan landmark'ların istatistiğini bırak; reinit açıksa
        kara-listeyi pencere-içiyle sınırla."""
        for lid in list(self.stats):
            if lid not in resident_lms:
                self.stats.pop(lid, None)
        if self.allow_reinit:
            self.bad &= set(resident_lms)

    def forget(self, lid):
        """Landmark pencereden düştü → istatistiğini bırak."""
        self.stats.pop(lid, None)

    def report(self):
        return dict(self.audit, bad_now=len(self.bad), tracked_stats=len(self.stats))
