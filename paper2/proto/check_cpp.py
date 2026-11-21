import os
import numpy as np, sys, ctypes as C
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wvar import W, sim
import ext
L=ext.L; dp=ext.dp; ip=ext.ip
def var_ipcw(t,c,X,th,vm,af=0,av=0.0,base='gompertz',knots=None):
    t,c,X=ext.prep(t,c,X); n,p=X.shape; kn=np.ascontiguousarray(knots if knots is not None else np.zeros(1),float); nk=0 if knots is None else len(knots)
    P=len(th); v=np.zeros((P,P)); vm=np.ascontiguousarray(vm,float); th=np.ascontiguousarray(th,float)
    rc=L.sc_subdist_var_ipcw(n,p,dp(t),ip(c),dp(X),1,af,C.c_double(av),ext.BASE[base],nk,dp(kn),dp(th),dp(vm),dp(v)); assert rc==0,rc
    return v
if __name__ == '__main__':
  rng=np.random.default_rng(11)
  for alpha,af in [(0.0,0),(1.0,0),(0.5,1)]:
      for ties in (False, True):
          x,k,Z=sim(300,alpha=alpha,rng=rng)
          if ties: x=np.round(x,1)
          f=ext.subdist_fit(x,k,Z,af=af,av=alpha)
          th=f['theta']; Vc=var_ipcw(x,k,Z,th,f['var_model'],af=af,av=alpha)
          m=W(x,k,Z,alpha=alpha,alpha_free=bool(af))
          S=m.score_matrix(th); Ps=m.psi(th); Hi=f['var_model']
          Vp=Hi@((S+Ps).T@(S+Ps))@Hi; Vn=Hi@(S.T@S)@Hi
          print(f"alpha={alpha} free={af} ties={ties}: rel.diff C++ vs Python (psi) = {np.abs(Vc-Vp).max()/np.abs(Vp).max():.2e};"
                f" naive C++ vs Python = {np.abs(f['var']-Vn).max()/np.abs(Vn).max():.2e}")
          print("   SE naive", np.round(np.sqrt(np.diag(f['var'])),4), " SE ipcw", np.round(np.sqrt(np.diag(Vc)),4))
