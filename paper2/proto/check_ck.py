import os
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ext
from ck import *
rng=np.random.default_rng(4)
# ---- data from a K-cause constrained model (Gompertz PH for cause 1, logit shares, Weibull G_k)
def simK(n,K,beta=0.6,rho=-0.8,ltau=np.log(0.5),gam=None,dlt=None,rng=rng):
    z=rng.normal(size=n)
    E=np.exp(beta*z); tau=np.exp(ltau)
    p1=1-np.exp(E*tau/rho)                  # F1(inf) = 1 - exp(-E u(inf)), u(inf) = -tau/rho
    l=np.zeros((n,K-1))
    for k in range(3,K+1): l[:,k-2]=gam[k-3]+dlt[k-3]*z
    pi=np.exp(l); pi/=pi.sum(1,keepdims=True)
    u=rng.uniform(size=n); cause=np.ones(n,int)
    rest=u>=p1
    cum=p1[:,None]+(1-p1)[:,None]*np.cumsum(pi,1)
    for i in np.where(rest)[0]: cause[i]=2+np.searchsorted(cum[i],u[i])
    cause=np.minimum(cause,K)
    T=np.empty(n)
    # cause 1 time given cause 1: invert F1(t)/p1
    i1=cause==1; v=rng.uniform(size=i1.sum())*p1[i1]
    T[i1]=np.log1p(rho*(-np.log1p(-v))/(E[i1]*tau))/rho
    scales=[1.0,2.0,0.7]; shapes=[1.5,0.8,2.0]
    for k in range(2,K+1):
        ik=cause==k; T[ik]=scales[k-2]*rng.weibull(shapes[k-2],ik.sum())
    Cc=rng.uniform(0,6,n); x=np.minimum(T,Cc); d=np.where(T<=Cc,cause,0)
    return np.maximum(x,1e-3),d,z[:,None]
if __name__ == "__main__":
  # 1. K=2 equals the two-cause model
  x,d,Z=simK(400,2,gam=[],dlt=[])
  for fam in ('logistic','weibull','lognormal'):
      a=ext.con_fit(x,d,Z,fam=fam); b=ck_fit(x,d,Z,2,[fam])
      print(f"K=2 {fam:9s}: |theta diff|={np.abs(a['theta']-b['theta']).max():.1e} |ll diff|={abs(a['loglik']-b['loglik']):.1e} |var diff|={np.abs(a['var']-b['var']).max():.1e}")
  # 2. derivatives against finite differences, K=3 and 4
  for K,fam,af in [(3,['weibull','loglogistic'],0),(4,['weibull','lognormal','logistic'],0),(3,['weibull','weibull'],1)]:
      x,d,Z=simK(500,K,gam=[-0.3,0.4][:K-2],dlt=[0.5,-0.4][:K-2])
      P=P_of('gompertz',af,1,K)
      th=np.r_[[-0.7,np.log(0.4)],([0.3] if af else []),[0.5]]
      th=np.r_[th,np.tile([0.1,0.2],K-1),np.tile([0.1,0.2],K-2)]+rng.normal(0,0.05,P)
      ll,g,H=ck_eval(x,d,Z,K,fam,th,af=af,av=0.3 if af else 0.0)
      h=1e-6; gn=np.zeros(P); Hn=np.zeros((P,P))
      for a in range(P):
          e=np.zeros(P); e[a]=h
          lp,gp,_=ck_eval(x,d,Z,K,fam,th+e,af=af); lm,gm,_=ck_eval(x,d,Z,K,fam,th-e,af=af)
          gn[a]=(lp-lm)/(2*h); Hn[:,a]=-(gp-gm)/(2*h)
      print(f"K={K} af={af}: grad rel err {np.abs(g-gn).max()/np.abs(g).max():.1e}  info rel err {np.abs(H-Hn).max()/np.abs(H).max():.1e}  events {np.bincount(d)}")
