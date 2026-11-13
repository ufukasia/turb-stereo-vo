"""SE(3) (Umeyama, ölçeksiz) hizalı ATE — eski ve yeni hattı aynı ölçütle kıyaslamak için."""
import sys, numpy as np, pandas as pd
def ate(est_csv, gt_csv):
    e=pd.read_csv(est_csv); e.columns=[c.strip() for c in e.columns]
    g=np.loadtxt(gt_csv,delimiter=',',comments='#')
    te=e.iloc[:,0].values.astype(np.int64); E=e.iloc[:,1:4].values
    gi={int(t):i for i,t in enumerate(g[:,0].astype(np.int64))}
    idx=[(i,gi[t]) for i,t in enumerate(te) if t in gi]
    A=E[[a for a,_ in idx]]; B=g[[b for _,b in idx],1:4]
    ma,mb=A.mean(0),B.mean(0); H=(A-ma).T@(B-mb); U,S,Vt=np.linalg.svd(H)
    D=np.eye(3); D[2,2]=np.sign(np.linalg.det(Vt.T@U.T)); R=Vt.T@D@U.T
    err=np.linalg.norm((A-ma)@R.T+mb-B,axis=1); yol=np.linalg.norm(np.diff(B,axis=0),axis=1).sum()
    return np.sqrt((err**2).mean()), yol, len(idx)
for est,gt in zip(sys.argv[1::2],sys.argv[2::2]):
    a,y,n=ate(est,gt); print(f"{est.split('/')[-2]}/{est.split('/')[-1]}: ATE_SE3={a:.3f} m  yol={y:.1f} m  ATE/yol={100*a/y:.2f}%  n={n}")
