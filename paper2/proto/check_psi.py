import os
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wvar import W, sim
rng = np.random.default_rng(7)
for alpha, af in [(0.0, False), (1.0, False), (0.5, True)]:
    x, k, Z = sim(150, alpha=alpha, rng=rng)
    x = np.round(x, 2)                                  # create ties
    m = W(x, k, Z, alpha=alpha, alpha_free=af)
    th0 = np.array([-0.5, np.log(0.3)] + ([0.4] if af else []) + [0.3])
    th = m.fit(th0)
    S = m.score_matrix(th); H = m.hess(th)
    Psi = m.psi(th); dU = m.dU_dv(th)
    print(f"alpha={alpha} free={af} n_events={np.bincount(k)} theta={np.round(th,4)}")
    print("  max |dU/dv - (s+psi)| =", np.abs(dU - (S + Psi)).max(), " max|psi| =", np.abs(Psi).max())
    Hi = np.linalg.inv(H)
    se_naive = np.sqrt(np.diag(Hi @ (S.T @ S) @ Hi)); se_psi = np.sqrt(np.diag(Hi @ ((S+Psi).T @ (S+Psi)) @ Hi))
    print("  SE naive", np.round(se_naive, 4), " SE with psi", np.round(se_psi, 4))
