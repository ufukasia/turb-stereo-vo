"""
E0 — türbülans kancalarının BİT DÜZEYİNDE etkisiz olduğunun kanıtı (yeni SLAM).

Kabul ölçütü (rapor §6/E0): türbülans kapalıyken hattın çıktısı, kancalar EKLENMEDEN
önceki kodun çıktısıyla özdeş olmalıdır. Bu, sonraki bütün farkların türbülanstan
geldiğinin garantisidir — ve SLAM her güncellendiğinde YENİDEN kanıtlanmalıdır.

İki bağımsız karşılaştırma yapılır:

  (1) ÖN-UÇ dump'ı — kancalı ve kancasız kod aynı sekansı işler, üretilen landmark
      dump'ları ALAN ALAN karşılaştırılır (float dizileri bit düzeyinde eşit mi).
  (2) ARKA-UÇ yörüngesi — aynı dump üzerinde kancalı ve kancasız backend koşulur,
      ATE ve poz dizisi karşılaştırılır.

Kancasız kod, bu depoya kopyalanmadan ÖNCEKİ orijinal SLAM'in bir kopyasıdır
(varsayılan: /tmp/e0/orijinal, `--pristine` ile değiştirilebilir).

    ./turbvio/bin/python -m turb.e0_bit_identity --frames 120 --sequence MH01
"""
from __future__ import annotations

import argparse
import os
import pickle
import subprocess
import sys

import numpy as np

HOOKED = "."


FRONTEND_SNIPPET = r"""
import sys, os
sys.path.insert(0, {root!r})
os.chdir({root!r})
import cv2; cv2.setNumThreads(2)
import config
from stereo_front_end import StereoMatcher
cfg = config.load()
cfg['dataset']['sequence'] = {seq!r}
cfg['frontend']['visualize'] = False
cfg['frontend']['save_visualizations'] = False
m = StereoMatcher(base_path=cfg['dataset']['base_path'], sequence={seq!r},
                  matcher_type=cfg['frontend']['matcher'],
                  filter_method=cfg['frontend']['filter'],
                  y_axis_threshold=cfg['frontend']['y_axis_threshold'],
                  extractor=cfg['frontend']['extractor'],
                  max_keypoints=cfg['frontend']['max_keypoints'],
                  detection_threshold=cfg['frontend']['detection_threshold'])
config.apply_frontend(m, cfg)
assert getattr(m, 'turb', None) is None, 'türbülans AÇIK — E0 geçersiz'
m.process_sequence(max_frames={n}, visualize=False, save_visualizations=False,
                   export_landmarks_path={out!r})
print('FE_OK')
"""

BACKEND_SNIPPET = r"""
import sys, os
sys.path.insert(0, {root!r})
os.chdir({root!r})
import cv2; cv2.setNumThreads(2)
import config
from gtsam_backend import main as bmain
cfg = config.load()
cfg['dataset']['sequence'] = {seq!r}
a = config.backend_args(cfg, argv=['--dump', {dump!r}, '--out', {out!r},
                                   '--max_frames', str({n})])
ate = bmain(a)
print('ATE_REPR=' + repr(float(ate)))
"""


def _run(py, code, tag):
    p = subprocess.run([py, "-c", code], capture_output=True, text=True)
    if p.returncode != 0:
        print(f"[{tag}] BAŞARISIZ:\n{p.stdout[-3000:]}\n{p.stderr[-3000:]}")
        sys.exit(2)
    return p.stdout


def _cmp_dumps(a_path, b_path):
    """İki landmark dump'ını alan alan karşılaştır; farkları döndür."""
    with open(a_path, "rb") as f:
        A = pickle.load(f)
    with open(b_path, "rb") as f:
        B = pickle.load(f)
    diffs = []
    fa, fb = A["frames"], B["frames"]
    if len(fa) != len(fb):
        return [f"kare sayısı {len(fa)} != {len(fb)}"], 0
    n_val = 0
    for i, (ra, rb) in enumerate(zip(fa, fb)):
        if set(ra) != set(rb):
            diffs.append(f"kare {i}: anahtar kümesi farklı")
            continue
        for k in ra:
            va, vb = ra[k], rb[k]
            if isinstance(va, np.ndarray):
                if va.shape != vb.shape:
                    diffs.append(f"kare {i} '{k}': şekil {va.shape} != {vb.shape}")
                elif va.dtype.kind == "f":
                    if not np.array_equal(va, vb, equal_nan=True):
                        m = float(np.nanmax(np.abs(va - vb))) if va.size else 0.0
                        diffs.append(f"kare {i} '{k}': max|fark| = {m:.3e}")
                    else:
                        n_val += va.size
                elif not np.array_equal(va, vb):
                    diffs.append(f"kare {i} '{k}': tamsayı/bool dizisi farklı")
                else:
                    n_val += va.size
            elif va != vb:
                diffs.append(f"kare {i} '{k}': {va!r} != {vb!r}")
            else:
                n_val += 1
    return diffs, n_val


def main():
    ap = argparse.ArgumentParser(description="E0 bit-düzeyi özdeşlik (yeni SLAM)")
    ap.add_argument("--pristine", default="/tmp/e0/orijinal",
                    help="kancasız orijinal SLAM kopyasının kökü")
    ap.add_argument("--sequence", default="MH01")
    ap.add_argument("--frames", type=int, default=120)
    ap.add_argument("--python", default=os.path.join(HOOKED, "turbvio/bin/python"))
    ap.add_argument("--out", default="evidence/E0_bit_identity.md")
    a = ap.parse_args()

    if not os.path.isdir(a.pristine):
        sys.exit(f"kancasız kopya yok: {a.pristine}")
    work = "/tmp/e0"
    os.makedirs(work, exist_ok=True)
    d_pri = os.path.join(work, "lm_pristine.pkl")
    d_hok = os.path.join(work, "lm_hooked.pkl")

    print("[1/4] kancasız ön-uç ...", flush=True)
    _run(a.python, FRONTEND_SNIPPET.format(root=a.pristine, seq=a.sequence,
                                           n=a.frames, out=d_pri), "kancasız-FE")
    print("[2/4] kancalı ön-uç ...", flush=True)
    _run(a.python, FRONTEND_SNIPPET.format(root=HOOKED, seq=a.sequence,
                                           n=a.frames, out=d_hok), "kancalı-FE")

    diffs, n_val = _cmp_dumps(d_pri, d_hok)
    fe_ok = not diffs
    print(f"   ön-uç: {'ÖZDEŞ' if fe_ok else 'FARKLI'} "
          f"({n_val} değer karşılaştırıldı, {len(diffs)} fark)")

    print("[3/4] kancasız arka-uç ...", flush=True)
    o1 = _run(a.python, BACKEND_SNIPPET.format(root=a.pristine, seq=a.sequence,
                                               dump=d_pri, n=a.frames,
                                               out=os.path.join(work, "t_pri.csv")),
              "kancasız-BE")
    print("[4/4] kancalı arka-uç ...", flush=True)
    o2 = _run(a.python, BACKEND_SNIPPET.format(root=HOOKED, seq=a.sequence,
                                               dump=d_pri, n=a.frames,
                                               out=os.path.join(work, "t_hok.csv")),
              "kancalı-BE")

    def _ate(o):
        for line in o.splitlines():
            if line.startswith("ATE_REPR="):
                return float(line.split("=", 1)[1])
        return float("nan")

    ate_p, ate_h = _ate(o1), _ate(o2)
    be_ok = (ate_p == ate_h)
    print(f"   arka-uç ATE: kancasız {ate_p!r} · kancalı {ate_h!r} -> "
          f"{'ÖZDEŞ' if be_ok else 'FARKLI'}")

    # yörünge dosyaları da bire bir mi?
    import filecmp
    traj_ok = filecmp.cmp(os.path.join(work, "t_pri.csv"),
                          os.path.join(work, "t_hok.csv"), shallow=False)

    ok = fe_ok and be_ok and traj_ok
    import datetime
    md = [
        "# E0 — Temel çizgi bit-düzeyi özdeşlik doğrulaması (YENİ SLAM)", "",
        f"Üretim: {datetime.datetime.now():%Y-%m-%d %H:%M} · "
        f"sekans **{a.sequence}**, {a.frames} kare · sonuç: "
        f"**{'GEÇTİ' if ok else 'KALDI'}**", "",
        "**Kabul ölçütü (rapor §6/E0):** türbülans kapalıyken hattın çıktısı, türbülans",
        "kodu eklenmeden ÖNCEKİ kodun çıktısıyla özdeş olmalıdır. SLAM her",
        "güncellendiğinde bu YENİDEN kanıtlanır — 2026-08-17'de arka-uç köklü",
        "değiştiği için kancalar yeniden uygulandı ve test yeniden koşuldu.", "",
        "## Ölçüm", "",
        "| Karşılaştırma | Sonuç |",
        "|---|---|",
        f"| Ön-uç landmark dump'ı (alan alan, {n_val} değer) | "
        f"{'**ÖZDEŞ**' if fe_ok else f'FARKLI ({len(diffs)} fark)'} |",
        f"| Arka-uç ATE (16 anlamlı basamak) | `{ate_p!r}` vs `{ate_h!r}` → "
        f"{'**ÖZDEŞ**' if be_ok else 'FARKLI'} |",
        f"| Yörünge CSV (bayt bayt) | {'**ÖZDEŞ**' if traj_ok else 'FARKLI'} |", "",
    ]
    if diffs:
        md += ["### Ön-uç farkları (ilk 20)", ""] + [f"- {d}" for d in diffs[:20]] + [""]
    md += [
        "## Neden özdeş", "",
        "Enjeksiyon noktalarının tamamı tek bir `None` kapısının arkasındadır ve",
        "kapalıyken hiçbir dizi kopyalanmaz, hiçbir hesap yapılmaz:", "",
        "* `stereo_front_end.StereoMatcher.turb is None` → rektifikasyondan sonraki",
        "  warp çağrıları hiç yürütülmez (`process_frame_pair` ve",
        "  `process_frame_temporal_only`).",
        "* `gtsam_backend` içinde `turb_builder is None` (kol A0) → gözlemler eskisi",
        "  gibi doğrudan `GenericStereoFactor3D` / `GenericProjectionFactorCal3_S2`",
        "  olarak eklenir; `search_by_projection` yolu da aynı kapıya bağlıdır.",
        "* `--cov_out` boş → pencere-göreli marjinal kovaryans hesabı hiç çağrılmaz.",
        "* `config.yaml`'daki `turbulence:` bölümü `enabled: false` ile gelir ve",
        "  `injector_from_cfg` `None` döndürür.", "",
        "Kancasız referans, orijinal SLAM deposunun değiştirilmemiş bir kopyasıdır",
        f"(`{a.pristine}`); yalnız `config.yaml`'ın `turbulence:` bölümü çıkarılmıştır,",
        "diğer bütün parametreler (loop closure/PGO/GBA dahil) AYNIDIR.",
    ]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"\n{'✓ GEÇTİ' if ok else '✗ KALDI'} -> {a.out}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
