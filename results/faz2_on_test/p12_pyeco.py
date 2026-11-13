"""P12'nin PyEco imzasıyla (gfg, keys, dims, n_pose_keys) yeniden koşulması — preflight'a dokunmadan."""
import sys, numpy as np
sys.path.insert(0, "/home/masa24/Desktop/makaleler/turbulans vio/pyeco_turb")
import gtsam
from gtsam.symbol_shorthand import X, L
import gtsam_backend as GB
from turb.factors import make_patch_factor
rg = np.random.default_rng(3)
fx = fy = 458.0; cx, cy = 367.0, 248.0; b = 0.11; n = 4
pw = np.array([[.5, .2, 4.], [-.3, .1, 5.], [.1, -.4, 3.5], [.4, .3, 4.5]])
po = gtsam.Pose3(); vl = gtsam.Values(); vl.insert(X(0), po)
for i, q in enumerate(pw):
    vl.insert(L(i), gtsam.Point3(*q))
mz = np.array(sum(([fx*q[0]/q[2]+cx, fy*q[1]/q[2]+cy] for q in pw), [])) + rg.normal(0, .3, 2*n)
Sg = np.full((2*n, 2*n), 1.4); np.fill_diagonal(Sg, 4.0)
lp = gtsam.noiseModel.Isotropic.Sigma(3, 10.0)
worst = 0.0
for patch in (False, True):
    gr = gtsam.NonlinearFactorGraph()
    if patch:
        f, _ = make_patch_factor(X(0), [L(i) for i in range(n)], mz, np.zeros(n, bool), Sg, fx, fy, cx, cy, b)
        gr.add(f)
    else:
        K = gtsam.Cal3_S2(fx, fy, 0., cx, cy)
        for i in range(n):
            gr.add(gtsam.GenericProjectionFactorCal3_S2(mz[2*i:2*i+2], gtsam.noiseModel.Isotropic.Sigma(2, 2.0),
                                                        X(0), L(i), K, False, False))
    gr.add(gtsam.PriorFactorPose3(X(0), po, gtsam.noiseModel.Isotropic.Sigma(6, 1e-6)))
    for i in range(n):
        gr.add(gtsam.PriorFactorPoint3(L(i), gtsam.Point3(*pw[i]), lp))
    keys = [X(0)] + [L(i) for i in range(n)]
    dims = [6] + [3]*n
    gfg = gr.linearize(vl)
    C6 = GB._pose_marginal_schur(gfg, keys, dims, 1)
    o = gtsam.Ordering()
    for k in keys:
        o.push_back(k)
    H, _ = gfg.hessian(o)
    ref = np.linalg.inv(np.asarray(H) + np.eye(sum(dims))*1e-9)[:6, :6]
    d = float(np.abs(C6 - ref).max() / np.abs(ref).max())
    worst = max(worst, d)
    Hll = np.asarray(H)[6:, 6:]
    off = Hll.copy()
    for j in range(n):
        off[3*j:3*j+3, 3*j:3*j+3] = 0
    print(("yama" if patch else "bireysel"), f"bağıl fark {d:.2e}", "Hll bağlı" if np.abs(off).max() > 1e-12 else "Hll blok-köşegen")
print("P12-PyEco", "GEÇTİ" if worst < 1e-6 else "KALDI", f"{worst:.2e}")
