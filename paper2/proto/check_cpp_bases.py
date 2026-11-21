import os
# psi for every baseline: Python formula with F(t; z) taken from the C++ prediction routine,
# versus the C++ meat_ipcw. Also the naive sandwich, as a check of the per-subject scores.
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wvar import W, sim
import ext
from check_cpp import var_ipcw
class WG(W):
    def __init__(s, x, k, Z, base, af, av, knots):
        super().__init__(x, k, Z); s.base, s.af_, s.av, s.kn = base, af, av, knots
    def F(s, th, Zr, t):            # Zr: rows, t: vector -> matrix rows x len(t)
        return ext.direct_predict1(th, Zr, np.atleast_1d(t), af=s.af_, av=s.av, base=s.base, knots=s.kn)
    def contrib(s, th, v=None, Gs=None):
        n = len(s.x)
        if Gs is None: Gs = s.km(np.ones(n) if v is None else v)[0]
        pts, mass = s.support(Gs); out = np.zeros(n)
        for i in range(n):
            z = s.Z[i:i+1]; t = s.x[i]
            if s.k[i] == 1:
                h = 1e-6 * max(t, 1e-3); Fp = s.F(th, z, [t + h, t - h])[0]; out[i] = np.log((Fp[0] - Fp[1]) / (2 * h))
            elif s.k[i] == 0:
                if t > 0: out[i] = np.log1p(-s.F(th, z, [t])[0, 0])
            else:
                sel = pts >= t; Fv = s.F(th, z, pts[sel])[0]
                out[i] = (mass[sel] * np.log1p(-Fv)).sum() / mass[sel].sum()
        return out
    def psi(s, th, h=1e-6):
        n = len(s.x); Gs, R, D = s.km(np.ones(n)); pts, mass = s.support(Gs)
        m = len(s.ct); M = len(pts); P = len(th); i2 = np.where(s.k == 2)[0]
        def LSm(t_): return np.log1p(-s.F(t_, s.Z[i2], pts))
        B = np.zeros((len(i2), M, P))
        for a in range(P):
            e = np.zeros(P); e[a] = h; B[:, :, a] = (LSm(th + e) - LSm(th - e)) / (2 * h)
        Q = np.zeros((m, P))
        for idx, kk in enumerate(i2):
            f = np.searchsorted(s.ct, s.x[kk], 'left')
            if f >= m: continue
            Gprev = 1.0 if f == 0 else Gs[f - 1]; rho = Gs[f:m] / Gprev; b = B[idx]
            bnext = np.vstack([b[1:M], np.zeros((1, P))])
            if M == m: bnext[m - 1] = 0.0
            Q[f:m] += np.cumsum((rho[:, None] * (bnext[f:m] - b[f:m]))[::-1], 0)[::-1]
        Psi = np.zeros((n, P))
        for r in range(m):
            c = s.ct[r]; inR = (s.x > c) | ((s.x == c) & (s.k == 0)); cens = (s.x == c) & (s.k == 0)
            if R[r] - D[r] <= 0: continue
            Psi += np.outer(inR * D[r] / (R[r] * (R[r] - D[r])) - cens / (R[r] - D[r]), Q[r])
        return Psi
rng = np.random.default_rng(3)
cases = [('gompertz', 0, 0.0), ('gompertz', 1, 0.5), ('weibull', 0, 1.0), ('gompertz3', 0, 0.0),
         ('logistic3', 0, 0.0), ('mix_weibull', 0, 0.0), ('mix_loglogistic', 1, 0.5), ('mix_lognormal', 0, 1.0), ('spline', 0, 0.0), ('spline', 1, 0.5)]
for base, af, av in cases:
    x, k, Z = sim(160, rng=rng); x = np.maximum(np.round(x, 2), 0.01)
    kn = None
    if base == 'spline':
        lt = np.log(np.sort(x[k == 1])); kn = np.quantile(lt, [0, 0.5, 1.0])
    try:
        f = ext.subdist_fit(x, k, Z, af=af, av=av, base=base, knots=kn)
    except AssertionError as e:
        print(base, af, 'fit failed', e); continue
    th = f['theta']; Hi = f['var_model']
    if not np.all(np.isfinite(Hi)): print(base, af, 'no information matrix', f['converged'], f['flags']); continue
    m = WG(x, k, Z, base, af, av, kn)
    S = m.score_matrix(th, h=1e-5); Ps = m.psi(th)
    Vn = Hi @ (S.T @ S) @ Hi; Vp = Hi @ ((S + Ps).T @ (S + Ps)) @ Hi
    Vc = var_ipcw(x, k, Z, th, Hi, af=af, av=av, base=base, knots=kn)
    rd = lambda A, B: np.abs(A - B).max() / np.abs(B).max()
    print(f"{base:16s} af={af} conv={f['converged']}  naive: {rd(f['var'], Vn):.1e}   ipcw: {rd(Vc, Vp):.1e}   "
          f"SE ratio ipcw/naive: {np.round(np.sqrt(np.diag(Vc)/np.diag(f['var'])),3)}")
