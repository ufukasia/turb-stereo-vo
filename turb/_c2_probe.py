"""C2 KESIN SINAVI: restorasyon TURBULANSTAN mi, HAREKETTEN mi bozuluyor?

Kontrol: ayni restorasyon yontemleri TEMIZ (turbulanssiz) karelere uygulanir.
  - PSNR temiz veride de cokuyorsa  -> sebep KAMERA HAREKETI (yontem alan disi)
  - yalniz turbulansli veride cokuyorsa -> sebep turbulans etkilesimi

    ./turbvio/bin/python -m turb._c2_probe
"""
import glob, os
import numpy as np
import cv2


def main():
    from .restoration import Restorer
    from .metrics import psnr, ssim
    base = "/data/sdata"
    cand = sorted(glob.glob(os.path.join(base, "*/mav0/cam0/data/*.png")))
    if not cand:
        cand = sorted(glob.glob(os.path.join(base, "*/*/mav0/cam0/data/*.png")))
    if not cand:
        print("goruntu bulunamadi"); return
    seq = cand[0].split("/mav0/")[0].split("/")[-1]
    imgs = [cv2.imread(f, cv2.IMREAD_GRAYSCALE) for f in cand[200:320]]
    imgs = [im for im in imgs if im is not None]
    print(f"sekans {seq} · {len(imgs)} TEMIZ kare (turbulans YOK)\n")
    print(f"  {'yontem':>16} {'PSNR(temiz,restore)':>21} {'SSIM':>7}")
    print(f"  {'none (referans)':>16} {'inf (birebir)':>21} {'1.000':>7}")
    for m in ("temporal_median", "lucky", "elastic"):
        r = Restorer(m, window=7, sharpen=0.6)
        ps, ss = [], []
        for im in imgs:
            out = r(im)
            if out.shape != im.shape:
                continue
            ps.append(psnr(im, out)); ss.append(ssim(im, out))
        ps = [p for p in ps if np.isfinite(p)]
        print(f"  {m:>16} {np.median(ps):>21.2f} {np.median(ss):>7.3f}")
    print("\nTURBULANSLI veride olculen (E3, D/r0=3): none 23.32 · elastic 18.03")
    print("                                          lucky 17.94 · median 17.04")


if __name__ == "__main__":
    main()
