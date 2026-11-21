import os
import sys as _s; _s.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sys, ctypes as C

import ext
L=ext.L; dp=ext.dp; ip=ext.ip; Info=ext.Info
def P_of(base,af,p,K): return ext.NB[base]+af+p + 2*(K-1) + (K-2)*(1+p)
def ck_fit(t,c,X,K,fam,af=0,av=0.0,base='gompertz',init=None,max_iter=500,tol=1e-10):
    t,c,X=ext.prep(t,c,X); n,p=X.shape; P=P_of(base,af,p,K); th=np.zeros(P); v=np.zeros((P,P)); info=Info()
    fa=np.ascontiguousarray([ext.COND[f] for f in fam],np.int32)
    ini=None if init is None else dp(np.ascontiguousarray(init,float))
    rc=L.sc_constrainedK_fit(n,p,dp(t),ip(c),dp(X),K,af,C.c_double(av),ext.BASE[base],ip(fa),ini,max_iter,C.c_double(tol),dp(th),dp(v),C.byref(info)); assert rc==0,rc
    return dict(theta=th,var=v,loglik=info.loglik,converged=bool(info.converged),flags=info.flags,iterations=info.iterations,grad=info.grad_maxabs)
def ck_eval(t,c,X,K,fam,th,af=0,av=0.0,base='gompertz'):
    t,c,X=ext.prep(t,c,X); n,p=X.shape; th=np.ascontiguousarray(th,float); P=len(th); ll=C.c_double(0); g=np.zeros(P); H=np.zeros((P,P))
    fa=np.ascontiguousarray([ext.COND[f] for f in fam],np.int32)
    rc=L.sc_constrainedK_eval(n,p,dp(t),ip(c),dp(X),K,af,C.c_double(av),ext.BASE[base],ip(fa),dp(th),C.byref(ll),dp(g),dp(H)); assert rc==0,rc
    return ll.value,g,H
def ck_predict(th,Xn,times,k,K,fam,af=0,av=0.0,base='gompertz',var=None,average=0):
    Xn=np.ascontiguousarray(Xn,float); nn,p=Xn.shape; times=np.ascontiguousarray(times,float); th=np.ascontiguousarray(th,float)
    fa=np.ascontiguousarray([ext.COND[f] for f in fam],np.int32)
    rows=1 if average else nn
    out=np.zeros((rows,len(times))); se=np.zeros((rows,len(times)))
    rc=L.sc_constrainedK_predict(p,K,af,C.c_double(av),ext.BASE[base],ip(fa),dp(th),None if var is None else dp(np.ascontiguousarray(var,float)),k,average,nn,dp(Xn),len(times),dp(times),dp(out),None if var is None else dp(se)); assert rc==0,rc
    return (out,se) if var is not None else out
