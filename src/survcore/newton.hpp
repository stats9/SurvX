// survcore: Newton-Raphson maximiser with step halving for concave
// (partial / pseudo) log-likelihoods, with honest convergence diagnostics.
// Convergence rule follows the relative log-likelihood change used for the
// Cox model in Therneau & Grambsch (2000), Modeling Survival Data, Sec. 3.
#pragma once
#include <cmath>
#include <vector>
#include "linalg.hpp"

namespace survcore {

struct Eval {
  double ll = 0.0;
  std::vector<double> g;   // gradient
  Mat H;                   // NEGATIVE Hessian (observed information)
};

enum FitFlag : int {
  FIT_OK = 0,
  FIT_MAX_ITER = 1,            // iteration limit reached
  FIT_INFO_NOT_PD = 2,         // information matrix not positive definite
  FIT_NONFINITE = 4,           // log-likelihood became non-finite
  FIT_POSSIBLE_MONOTONE = 8    // coefficient may be diverging (separation)
};

struct FitResult {
  std::vector<double> beta;
  std::vector<double> last_step;   // last accepted Newton step
  Mat var;                 // inverse information; empty (n=0) if unavailable
  double ll_init = 0.0, ll = 0.0;
  double grad_maxabs = 0.0;
  int iterations = 0;
  int flags = FIT_OK;
  bool converged = false;
};

template <class Model>
FitResult newton_fit(const Model& m, std::vector<double> beta, int max_iter,
                     double tol) {
  FitResult r;
  const std::size_t p = beta.size();
  Eval ev;
  m.eval(beta, ev, true);
  r.ll_init = ev.ll;
  if (!std::isfinite(ev.ll)) { r.flags |= FIT_NONFINITE; r.beta = beta; r.ll = ev.ll; return r; }
  std::vector<double> step(p), trial(p);
  Eval evt;
  for (int it = 1; it <= max_iter; ++it) {
    r.iterations = it;
    Mat L = ev.H;
    if (!cholesky(L)) { r.flags |= FIT_INFO_NOT_PD; break; }
    step = ev.g;
    chol_solve(L, step);
    double scale = 1.0;
    bool accepted = false;
    for (int h = 0; h < 40; ++h) {
      for (std::size_t j = 0; j < p; ++j) trial[j] = beta[j] + scale * step[j];
      m.eval(trial, evt, true);
      if (std::isfinite(evt.ll) && evt.ll >= ev.ll - 1e-13 * (std::fabs(ev.ll) + 1.0)) { accepted = true; break; }
      scale *= 0.5;
    }
    if (!accepted) { r.flags |= FIT_NONFINITE; break; }
    const double ll_old = ev.ll;
    r.last_step.assign(p, 0.0);
    for (std::size_t j = 0; j < p; ++j) r.last_step[j] = scale * step[j];
    beta = trial;
    std::swap(ev, evt);
    if (std::fabs(ev.ll - ll_old) <= tol * (std::fabs(ev.ll) + tol)) { r.converged = true; break; }
  }
  if (!r.converged && r.flags == FIT_OK) r.flags |= FIT_MAX_ITER;
  r.beta = beta;
  r.ll = ev.ll;
  for (double gj : ev.g) r.grad_maxabs = std::fmax(r.grad_maxabs, std::fabs(gj));
  Mat L = ev.H;
  if (cholesky(L)) r.var = chol_inverse(L);
  else r.flags |= FIT_INFO_NOT_PD;
  return r;
}

}  // namespace survcore
