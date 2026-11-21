import os
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ck import *
from check_ck import simK
rng=np.random.default_rng(9)
# independent log-likelihood: Gompertz PH cause 1, Weibull G_k: G = 1 - exp(-exp(e^b (a + log t)))
def ll_ind(th,x,d,z,K):
    rho,lt,beta=th[:3]; tau=np.exp(lt); E=np.exp(beta*z)
    u=tau*np.expm1(rho*x)/rho; F1=1-np.exp(-E*u); f1=E*tau*np.exp(rho*x)*np.exp(-E*u)
    p1=1-np.exp(E*tau/rho); q=1-p1
    G=[];g=[]
    for k in range(2,K+1):
        a,b=th[3+2*(k-2)],th[4+2*(k-2)]; lw=np.exp(b)*(a+np.log(x)); w=np.exp(lw)
        G.append(1-np.exp(-w)); g.append(np.exp(b)/x*w*np.exp(-w))
    l=np.zeros((len(x),K-1)); o=3+2*(K-1)
    for k in range(3,K+1): l[:,k-2]=th[o+2*(k-3)]+th[o+2*(k-3)+1]*z
    pi=np.exp(l); pi/=pi.sum(1,keepdims=True)
    S=1-F1-q*sum(pi[:,k-2]*G[k-2] for k in range(2,K+1))
    out=np.where(d==1,np.log(f1),0.0)
    for k in range(2,K+1): out=np.where(d==k,np.log(q*pi[:,k-2]*g[k-2]),out)
    out=np.where(d==0,np.log(S),out)
    return out.sum()
K=3; x,d,Z=simK(600,K,gam=[-0.3],dlt=[0.5],rng=rng)
th=np.array([-0.8,np.log(0.5),0.6, 0.0,np.log(1.5), -np.log(2.0),np.log(0.8), -0.3,0.5])
print("loglik C++", ck_eval(x,d,Z,K,['weibull','weibull'],th)[0], " independent", ll_ind(th,x,d,Z[:,0],K))
# predictions: sum of the K CIFs at infinity is 1 for every z
f=ck_fit(x,d,Z,K,['weibull','weibull'])
zz=np.linspace(-3,3,7)[:,None]
tot=sum(ck_predict(f['theta'],zz,[np.inf],k,K,['weibull','weibull'])[:,0] for k in range(1,K+1))
print("converged",f['converged'],"sum of long-term incidences:",np.round(tot,12))
# recovery: 300 replications, n = 800
true=th.copy(); est=[]; se=[]; cov=[]
for r in range(300):
    x,d,Z=simK(800,K,gam=[-0.3],dlt=[0.5],rng=rng)
    try: f=ck_fit(x,d,Z,K,['weibull','weibull'])
    except AssertionError: continue
    if not f['converged']: continue
    s=np.sqrt(np.diag(f['var'])); est.append(f['theta']); se.append(s); cov.append(np.abs(f['theta']-true)<=1.96*s)
est=np.array(est); se=np.array(se); cov=np.array(cov)
print("reps",len(est))
print("bias      ",np.round(est.mean(0)-true,3))
print("SE/SD     ",np.round(se.mean(0)/est.std(0,ddof=1),3))
print("coverage  ",np.round(cov.mean(0),3))
