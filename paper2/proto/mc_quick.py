import os
import sys as _s; _s.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sys

from wvar import sim
import ext
from check_cpp import var_ipcw
rng=np.random.default_rng(2026)
for alpha in (0.0, 1.0):
  for cmax in (6.0, 2.5):
    est=[]; sen=[]; sei=[]
    for r in range(400):
        x,k,Z=sim(300,alpha=alpha,cmax=cmax,rng=rng)
        try: f=ext.subdist_fit(x,k,Z,av=alpha)
        except AssertionError: continue
        if not f['converged']: continue
        V=var_ipcw(x,k,Z,f['theta'],f['var_model'],av=alpha)
        est.append(f['theta']); sen.append(np.sqrt(np.diag(f['var']))); sei.append(np.sqrt(np.diag(V)))
    est=np.array(est); sen=np.array(sen); sei=np.array(sei)
    sd=est.std(0,ddof=1)
    print(f"alpha={alpha} cmax={cmax} reps={len(est)} empSD={np.round(sd,4)} naiveSE/SD={np.round(sen.mean(0)/sd,3)} ipcwSE/SD={np.round(sei.mean(0)/sd,3)}")
