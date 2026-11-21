import ctypes as C, numpy as np
import os
L = C.CDLL(os.environ.get('SURVX_LIB', '/tmp/sx2.so'))
class Info(C.Structure):
    _fields_ = [("loglik_init", C.c_double), ("loglik", C.c_double), ("grad_maxabs", C.c_double), ("iterations", C.c_int), ("converged", C.c_int), ("flags", C.c_int), ("n_events", C.c_int)]
dp = lambda a: a.ctypes.data_as(C.POINTER(C.c_double)); ip = lambda a: a.ctypes.data_as(C.POINTER(C.c_int))
BASE = dict(gompertz=0, gompertz3=1, weibull=2, mix_weibull=3, mix_loglogistic=4, mix_lognormal=5, logistic3=6, spline=7)
NB = dict(gompertz=2, gompertz3=3, weibull=2, mix_weibull=3, mix_loglogistic=3, mix_lognormal=3, logistic3=3)
COND = dict(logistic=0, weibull=1, loglogistic=2, lognormal=3)
def prep(t, c, X):
    t = np.ascontiguousarray(t, float); c = np.ascontiguousarray(c, np.int32); X = np.ascontiguousarray(X, float).reshape(len(t), -1); return t, c, X
def P1(base, af, p, kn=None): return (max(len(kn), 2) if base == 'spline' else NB[base]) + af + p
def subdist_fit(t, c, X, interest=1, af=0, av=0.0, base='gompertz', knots=None, grid=0, init=None, max_iter=500, tol=1e-10):
    t, c, X = prep(t, c, X); n, p = X.shape; kn = np.ascontiguousarray(knots if knots is not None else np.zeros(1), float); nk = 0 if knots is None else len(knots)
    P = P1(base, af, p, knots); th = np.zeros(P); vr = np.zeros((P, P)); vm = np.zeros((P, P)); info = Info(); ns = C.c_int(0)
    ini = None if init is None else dp(np.ascontiguousarray(init, float))
    rc = L.sc_subdist_fit(n, p, dp(t), ip(c), dp(X), interest, af, C.c_double(av), BASE[base], nk, dp(kn), grid, ini, max_iter, C.c_double(tol), dp(th), dp(vr), dp(vm), C.byref(info), C.byref(ns))
    assert rc == 0, rc
    return dict(theta=th, var=vr, var_model=vm, loglik=info.loglik, converged=bool(info.converged), flags=info.flags, iterations=info.iterations, grad=info.grad_maxabs, n_support=ns.value)
def subdist_eval(t, c, X, th, interest=1, af=0, av=0.0, base='gompertz', knots=None, grid=0):
    t, c, X = prep(t, c, X); n, p = X.shape; kn = np.ascontiguousarray(knots if knots is not None else np.zeros(1), float); nk = 0 if knots is None else len(knots)
    th = np.ascontiguousarray(th, float); P = len(th); ll = C.c_double(0); g = np.zeros(P); H = np.zeros((P, P))
    rc = L.sc_subdist_eval(n, p, dp(t), ip(c), dp(X), interest, af, C.c_double(av), BASE[base], nk, dp(kn), grid, dp(th), C.byref(ll), dp(g), dp(H)); assert rc == 0, rc
    return ll.value, g, H
def con_fit(t, c, X, af=0, av=0.0, base='gompertz', fam='logistic', init=None, max_iter=500, tol=1e-10):
    t, c, X = prep(t, c, X); n, p = X.shape; P = P1(base, af, p) + 2; th = np.zeros(P); v = np.zeros((P, P)); info = Info()
    ini = None if init is None else dp(np.ascontiguousarray(init, float))
    rc = L.sc_constrained_fit(n, p, dp(t), ip(c), dp(X), af, C.c_double(av), BASE[base], COND[fam], ini, max_iter, C.c_double(tol), dp(th), dp(v), C.byref(info)); assert rc == 0, rc
    return dict(theta=th, var=v, loglik=info.loglik, converged=bool(info.converged), flags=info.flags, iterations=info.iterations, grad=info.grad_maxabs)
def con_eval(t, c, X, th, af=0, av=0.0, base='gompertz', fam='logistic'):
    t, c, X = prep(t, c, X); n, p = X.shape; th = np.ascontiguousarray(th, float); P = len(th); ll = C.c_double(0); g = np.zeros(P); H = np.zeros((P, P))
    rc = L.sc_constrained_eval(n, p, dp(t), ip(c), dp(X), af, C.c_double(av), BASE[base], COND[fam], dp(th), C.byref(ll), dp(g), dp(H)); assert rc == 0, rc
    return ll.value, g, H
def con_predict(th, Xn, times, k, af=0, av=0.0, base='gompertz', fam='logistic', var=None):
    Xn = np.ascontiguousarray(Xn, float); nn, p = Xn.shape; times = np.ascontiguousarray(times, float); th = np.ascontiguousarray(th, float)
    out = np.zeros((nn, len(times))); se = np.zeros((nn, len(times)))
    rc = L.sc_constrained_predict(p, af, C.c_double(av), BASE[base], COND[fam], dp(th), None if var is None else dp(np.ascontiguousarray(var, float)), k, nn, dp(Xn), len(times), dp(times), dp(out), None if var is None else dp(se)); assert rc == 0, rc
    return (out, se) if var is not None else out
def direct_predict1(th, Xn, times, af=0, av=0.0, base='gompertz', knots=None, var=None):
    Xn = np.ascontiguousarray(Xn, float); nn, p = Xn.shape; times = np.ascontiguousarray(times, float); th = np.ascontiguousarray(th, float)
    kn = np.ascontiguousarray(knots if knots is not None else np.zeros(1), float); nk = np.array([0 if knots is None else len(knots)], np.int32)
    out = np.zeros((nn, len(times))); se = np.zeros((nn, len(times)))
    rc = L.sc_directcif_predict(p, 1, ip(np.array([af], np.int32)), dp(np.array([av], float)), ip(np.array([BASE[base]], np.int32)), ip(nk), dp(kn), dp(th), None if var is None else dp(np.ascontiguousarray(var, float)), 1, nn, dp(Xn), len(times), dp(times), dp(out), None if var is None else dp(se)); assert rc == 0, rc
    return (out, se) if var is not None else out
