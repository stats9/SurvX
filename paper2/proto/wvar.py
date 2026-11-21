# Reference implementation of the weighted (IPCW) direct-CIF estimator and of its variance
# with and without the term for the estimated censoring distribution.
# Exact check: the influence of each subject from the implicit derivative of theta-hat in case
# weights (the KM of censoring re-estimated with the weights), compared with the analytic psi.
import numpy as np
from scipy.optimize import minimize

def u_gomp(t, rho, ltau):
    tau = np.exp(ltau)
    return tau * np.expm1(rho * t) / rho

def logS(t, eta, th, alpha):            # log{1 - F(t; z)}
    u = u_gomp(t, th[0], th[1]); a = alpha
    A = np.exp(eta) * u if a == 0 else np.log1p(a * np.exp(eta) * u) / a
    return -A

def logf(t, eta, th, alpha):            # log dF/dt
    rho, ltau = th[0], th[1]; tau = np.exp(ltau)
    u = u_gomp(t, rho, ltau); du = tau * np.exp(rho * t)
    E = np.exp(eta)
    if alpha == 0: return eta + np.log(du) - E * u
    return eta + np.log(du) - (1 / alpha + 1) * np.log1p(alpha * E * u)

class W:
    """kind: 0 cens, 1 interest, 2 competing. theta = (rho, log tau, [alpha], beta...)"""
    def __init__(s, x, kind, Z, alpha=0.0, alpha_free=False):
        s.x, s.k, s.Z = np.asarray(x, float), np.asarray(kind), np.atleast_2d(Z).reshape(len(x), -1)
        s.alpha, s.af = alpha, alpha_free
        s.ct = np.unique(s.x[s.k == 0]); s.tmax = s.x.max()
    def split(s, th):
        b = 2; al = s.alpha
        if s.af: al = th[2]; b = 3
        return th[:2], al, th[b:]
    def km(s, v):
        # weighted KM of censoring; failures precede censorings at ties
        x, k = s.x, s.k; G = 1.0; Gs = []; R = []; D = []
        for c in s.ct:
            Rj = v[(x > c) | ((x == c) & (k == 0))].sum(); dj = v[(x == c) & (k == 0)].sum()
            R.append(Rj); D.append(dj); G = G * (1 - dj / Rj); Gs.append(G)
        return np.array(Gs), np.array(R), np.array(D)
    def support(s, Gs):
        # support points and masses: censoring times, plus the remaining mass at tmax
        prev = np.r_[1.0, Gs[:-1]]; mass = prev - Gs; pts = s.ct.copy()
        if Gs[-1] > 1e-14: pts = np.r_[pts, s.tmax]; mass = np.r_[mass, Gs[-1]]
        return pts, mass
    def contrib(s, th, v=None, Gs=None):
        n = len(s.x)
        if Gs is None: Gs = s.km(np.ones(n) if v is None else v)[0]
        pts, mass = s.support(Gs)
        bl, al, be = s.split(th); eta = s.Z @ be
        out = np.zeros(n)
        i1 = s.k == 1; out[i1] = logf(s.x[i1], eta[i1], bl, al)
        i0 = (s.k == 0) & (s.x > 0); out[i0] = logS(s.x[i0], eta[i0], bl, al)
        i2 = np.where(s.k == 2)[0]
        LS = logS(pts[None, :], eta[i2, None], bl, al)              # n2 x m
        W_ = np.where(pts[None, :] >= s.x[i2, None], mass[None, :], 0.0)
        out[i2] = (W_ * LS).sum(1) / W_.sum(1)
        return out
    def fit(s, th0):
        f = lambda th: -s.contrib(th).sum()
        r = minimize(f, th0, method='BFGS', options={'gtol': 1e-10, 'maxiter': 2000})
        th = r.x
        for _ in range(5):                                           # polish with Newton
            g = s.score_matrix(th).sum(0); H = s.hess(th)
            th = th + np.linalg.solve(H, g)
        return th
    def score_matrix(s, th, h=1e-6, **kw):
        P = len(th); S = np.zeros((len(s.x), P))
        for a in range(P):
            e = np.zeros(P); e[a] = h
            S[:, a] = (s.contrib(th + e, **kw) - s.contrib(th - e, **kw)) / (2 * h)
        return S
    def hess(s, th, h=1e-5):                                         # minus second derivative
        P = len(th); H = np.zeros((P, P))
        for a in range(P):
            e = np.zeros(P); e[a] = h
            H[:, a] = -(s.score_matrix(th + e).sum(0) - s.score_matrix(th - e).sum(0)) / (2 * h)
        return (H + H.T) / 2
    # --- analytic psi (the formula to be implemented in C++) -----------------------------
    def psi(s, th, h=1e-6):
        n = len(s.x); Gs, R, D = s.km(np.ones(n)); pts, mass = s.support(Gs)
        m = len(s.ct); M = len(pts); P = len(th)
        bl, al, be = s.split(th)
        i2 = np.where(s.k == 2)[0]
        # b_kr = d log S(pts_r; z_k) / d theta
        def LSm(t_): b_, a_, e_ = s.split(t_); return logS(pts[None, :], (s.Z[i2] @ e_)[:, None], b_, a_)
        B = np.zeros((len(i2), M, P))
        for a in range(P):
            e = np.zeros(P); e[a] = h; B[:, :, a] = (LSm(th + e) - LSm(th - e)) / (2 * h)
        Q = np.zeros((m, P))                         # Q_r' = sum_k sum_{r>=r'} rho_kr (b_k,r+1 - b_kr)
        for idx, kk in enumerate(i2):
            f = np.searchsorted(s.ct, s.x[kk], 'left')       # first censoring index >= T_k
            if f >= m: continue
            Gprev = 1.0 if f == 0 else Gs[f - 1]
            rho = Gs[f:m] / Gprev                             # r = f..m-1
            b = B[idx]                                        # M x P, index r over pts
            bnext = np.vstack([b[1:M], np.zeros((1, P))])     # b_{r+1}; zero if no tmax point
            if M == m: bnext[m - 1] = 0.0                     # G_last = 0: no term
            term = rho[:, None] * (bnext[f:m] - b[f:m])
            Q[f:m] += np.cumsum(term[::-1], 0)[::-1]
        # dh_r/dv_i = -1(i cens at c_r)/(R-D) + 1(i in R_r) D/(R(R-D))
        Psi = np.zeros((n, P))
        for r in range(m):
            c = s.ct[r]
            inR = (s.x > c) | ((s.x == c) & (s.k == 0)); cens = (s.x == c) & (s.k == 0)
            if R[r] - D[r] <= 0: continue
            Psi += np.outer(inR * D[r] / (R[r] * (R[r] - D[r])) - cens / (R[r] - D[r]), Q[r])
        return Psi
    # --- exact influence via implicit derivative in case weights ------------------------------
    def dU_dv(s, th, eps=1e-5):
        n = len(s.x); P = len(th); out = np.zeros((n, P))
        for i in range(n):
            vp = np.ones(n); vp[i] += eps; vm = np.ones(n); vm[i] -= eps
            Up = (vp[:, None] * s.score_matrix(th, v=vp)).sum(0)
            Um = (vm[:, None] * s.score_matrix(th, v=vm)).sum(0)
            out[i] = (Up - Um) / (2 * eps)
        return out

def sim(n, beta=0.5, rho=-0.6, ltau=np.log(0.4), alpha=0.0, cmax=6.0, rng=None):
    rng = rng or np.random.default_rng()
    z = rng.binomial(1, 0.5, n).astype(float)
    eta = beta * z
    th = np.array([rho, ltau])
    p1 = 1 - np.exp(logS(np.inf * np.ones(n), eta, th, alpha)) if rho < 0 else np.ones(n)
    p1 = -np.expm1(logS(np.full(n, 1e6), eta, th, alpha))
    is1 = rng.uniform(size=n) < p1
    T = np.empty(n)
    # cause-1 times by inversion of F/F(inf), on a fine grid
    grid = np.linspace(0, 40, 40001)
    for g in (0.0, 1.0):
        sel = is1 & (z == g)
        Fg = -np.expm1(logS(grid, beta * g, th, alpha)); Fg = Fg / Fg[-1]
        T[sel] = np.interp(rng.uniform(size=sel.sum()), Fg, grid)
    T[~is1] = rng.exponential(1.5, (~is1).sum())
    C = rng.uniform(0, cmax, n)
    x = np.minimum(T, C); kind = np.where(T <= C, np.where(is1, 1, 2), 0)
    return x, kind, z[:, None]
