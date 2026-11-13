"""Rejim hesabi: nadir irtifa vs egik/siyiran gorus hatti icin tilt ve stereo derinlik.
Cn2(h) = C0 (h/2m)^(-4/3) (Wyngaard 1971, gunduz yuzey tabakasi), h<1 m sabit.
Kuresel dalga r0 = [0.423 k^2 int Cn2(h(z)) (z/L)^(5/3) dz]^(-3/5), z nesneden olculur.
Tilt (tek eksen) sigma = sqrt(0.170 (D/r0)^(5/3)) lambda/D  [rad] x f [px/rad]."""
import numpy as np
lam=550e-9;k=2*np.pi/lam
def cn2(h,C0): h=np.maximum(h,1.0); return C0*(h/2.0)**(-4/3)
def r0_path(hfun,L,C0,n=4000):
    z=np.linspace(1e-3,L,n)
    return (0.423*k*k*np.trapz(cn2(hfun(z),C0)*(z/L)**(5/3),z))**(-3/5)
def tilt_px(r0,D,f): return np.sqrt(0.170*(D/r0)**(5/3))*(lam/D)*f
cams={'drone(f731,D1.3mm)':(731,1.26e-3),'zoom(f5000,D10mm)':(5000,10e-3)}
for C0 in [1e-13,1e-12]:
  print(f"\n=== Cn2(2 m) = {C0:.0e} ===")
  print("NADIR: irtifa H")
  for H in [15,30,60,120,300]:
    r0=r0_path(lambda z:z,H,C0)
    s=" | ".join(f"{n}: {tilt_px(r0,D,f):.4f} px" for n,(f,D) in cams.items())
    disp=731*0.4/H
    print(f"  H={H:4d} m  r0={r0*1e3:8.1f} mm  {s} | disparite(B=0.4) {disp:5.1f} px, dZ/Z={0.25/disp*100:4.1f}%")
  print("EGIK/SIYIRAN: kamera yuksekligi Hc, zemin noktasi yatay mesafe L")
  for Hc in [3,5,10]:
    for L in [50,100,300,1000]:
      r0=r0_path(lambda z,Hc=Hc,L=L:Hc*z/L,L,C0)
      s=" | ".join(f"{n}: {tilt_px(r0,D,f):.3f} px" for n,(f,D) in cams.items())
      print(f"  Hc={Hc:3d} m L={L:5d} m  r0={r0*1e3:8.1f} mm  {s}")
