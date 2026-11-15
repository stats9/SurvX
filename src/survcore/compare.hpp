// survcore: side-by-side comparison of regression models for the cumulative
// incidence of one cause, by apparent and K-fold cross-validated prediction
// performance (Brier score, IPA, time-dependent AUC, calibration error).
//
// Models (all fitted by this library): any number of direct models, each defined by a
// baseline family (the same for every cause), a link (0 proportional hazards, 1
// proportional odds, 2 generalized odds rate with alpha estimated), for spline
// baselines the degrees of freedom, and the estimator (0 joint likelihood of all causes,
// 1 weighted likelihood for the cause of interest, 2 constrained model of Shi, Cheng & Jeong
// 2013 with conditional family `cond`; see directcif_ext.hpp); followed by
//   Fine-Gray                               (Fine & Gray 1999)
//   cause-specific Cox, product-limit CIF   (Prentice et al. 1978; Ozenne et al. 2017)
//
// The joint and the constrained direct models have a full likelihood of the same data, so
// their log-likelihood and AIC are comparable with each other; the weighted direct model
// has a pseudo-likelihood (no AIC is reported). Fine-Gray and cause-specific Cox use partial likelihoods;
// across all models only the prediction metrics are comparable.
// Cross-validation: Stone (1974) JRSS-B 36:111-147; predictions of the held-out folds
// are pooled and scored once (Gerds & Schumacher 2007, Biometrics 63:1283-1287).
#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include "cox.hpp"
#include "directcif.hpp"
#include "directcif_ext.hpp"
#include "finegray.hpp"
#include "metrics.hpp"

namespace survcore {

struct DirectChoice { int baseline = BASE_GOMPERTZ; int link = 0; int spline_df = 3; int estimation = 0; int cond = COND_LOGISTIC; };
constexpr int COMPARE_METRICS = 4;      // brier, ipa, auc, calibration error

struct ModelPrediction {
  std::vector<double> risk;             // n_test x nt
  double loglik = std::nan(""), npar = std::nan(""), aic = std::nan("");
  int flags = 0;
  bool converged = false;
};

inline double step_value(const std::vector<double>& x, const std::vector<double>& y, double at) {
  const std::size_t k = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), at) - x.begin());
  return k == 0 ? 0.0 : y[k - 1];
}

// Fits one model on the training data and predicts the risk of `interest` (1..K) for the
// test rows at each time. id < n_direct: direct model `choices[id]`; id == n_direct:
// Fine-Gray; id == n_direct + 1: cause-specific Cox.
inline ModelPrediction fit_predict(int id, const std::vector<DirectChoice>& choices, std::size_t n, std::size_t p, int K,
                                   const double* time, const int* cause, const double* X, int interest, std::size_t nte,
                                   const double* Xte, const std::vector<double>& times) {
  const std::size_t nt = times.size();
  const int nd = static_cast<int>(choices.size());
  ModelPrediction out;
  out.risk.assign(nte * nt, std::nan(""));
  if (id < nd && choices[id].estimation != 0) {
    // one directly modelled cause: weighted likelihood (1) or constrained model for two causes (2)
    const DirectChoice& ch = choices[id];
    DirectSpec s;
    s.K = 1; s.p = p;
    s.alpha_free.assign(1, ch.link == 2 ? 1 : 0);
    s.alpha_value.assign(1, ch.link == 1 ? 1.0 : (ch.link == 2 ? 0.5 : 0.0));
    s.baseline.assign(1, ch.baseline);
    s.knots.assign(1, std::vector<double>());
    std::vector<int> c12(n);
    for (std::size_t i = 0; i < n; ++i) c12[i] = cause[i] == 0 ? 0 : (cause[i] == interest ? 1 : 2);
    if (ch.estimation == 1) {
      if (ch.baseline == BASE_SPLINE) s.knots[0] = direct_default_knots(n, time, c12.data(), 1, ch.spline_df);
      SubdistCIFModel m(n, time, c12.data(), X, s, 1, 0);
      const SubdistFit f = subdist_fit(m, n, time, c12.data(), nullptr, 500, 1e-10);
      out.loglik = f.fit.ll; out.npar = static_cast<double>(s.nparam());
      out.flags = f.fit.flags; out.converged = f.fit.converged;
      if (!std::isfinite(f.fit.ll)) return out;
      for (std::size_t i = 0; i < nte; ++i)
        for (std::size_t t = 0; t < nt; ++t) out.risk[i * nt + t] = m.cif(f.fit.beta, 0, Xte + i * p, times[t]);
    } else {
      if (K != 2 || ch.baseline == BASE_WEIBULL || ch.baseline == BASE_SPLINE) { out.flags = FIT_NONFINITE; return out; }
      ConstrainedCIFModel m(n, time, c12.data(), X, s, ch.cond);
      const FitResult f = constrained_fit(m, n, time, c12.data(), nullptr, 500, 1e-10);
      out.loglik = f.ll; out.npar = static_cast<double>(m.nparam()); out.aic = -2.0 * f.ll + 2.0 * out.npar;
      out.flags = f.flags; out.converged = f.converged;
      if (!std::isfinite(f.ll)) return out;
      for (std::size_t i = 0; i < nte; ++i)
        for (std::size_t t = 0; t < nt; ++t) out.risk[i * nt + t] = m.cif12(f.beta, 1, Xte + i * p, times[t]);
    }
  } else if (id < nd) {
    const DirectChoice& ch = choices[id];
    DirectSpec s;
    s.K = K; s.p = p;
    s.alpha_free.assign(K, ch.link == 2 ? 1 : 0);
    s.alpha_value.assign(K, ch.link == 1 ? 1.0 : (ch.link == 2 ? 0.5 : 0.0));
    s.baseline.assign(K, ch.baseline);
    s.knots.assign(K, std::vector<double>());
    if (ch.baseline == BASE_SPLINE) for (int k = 0; k < K; ++k) s.knots[k] = direct_default_knots(n, time, cause, k + 1, ch.spline_df);
    DirectCIFModel m(n, time, cause, X, s);
    const DirectFit f = direct_fit(m, n, time, cause, nullptr, false, 500, 1e-10);
    out.loglik = f.fit.ll; out.npar = static_cast<double>(s.nparam()); out.aic = -2.0 * f.fit.ll + 2.0 * out.npar;
    out.flags = f.fit.flags; out.converged = f.fit.converged;
    if (!std::isfinite(f.fit.ll)) return out;
    for (std::size_t i = 0; i < nte; ++i)
      for (std::size_t t = 0; t < nt; ++t) out.risk[i * nt + t] = m.cif(f.fit.beta, interest - 1, Xte + i * p, times[t]);
  } else if (id == nd) {
    FineGrayModel m(n, p, time, cause, X, interest);
    if (m.n_events() == 0) { out.flags = FIT_NONFINITE; return out; }
    const FitResult f = newton_fit(m, std::vector<double>(p, 0.0), 50, 1e-10);
    out.loglik = f.ll; out.npar = static_cast<double>(p); out.flags = f.flags; out.converged = f.converged;
    std::vector<double> H;
    m.baseline(f.beta, nullptr, H);
    const std::vector<double>& et = m.event_times();
    for (std::size_t i = 0; i < nte; ++i) {
      double lp = 0.0;
      for (std::size_t j = 0; j < p; ++j) lp += Xte[i * p + j] * f.beta[j];
      for (std::size_t t = 0; t < nt; ++t) out.risk[i * nt + t] = -std::expm1(-step_value(et, H, times[t]) * std::exp(lp));
    }
  } else {
    // K cause-specific Cox models (Efron ties), product-limit CIF on the union of jump times.
    std::vector<std::vector<double>> bt(K), bh(K), beta(K);
    out.npar = static_cast<double>(p * static_cast<std::size_t>(K)); out.converged = true; out.loglik = 0.0;
    std::vector<int> ev(n);
    for (int k = 0; k < K; ++k) {
      for (std::size_t i = 0; i < n; ++i) ev[i] = cause[i] == k + 1;
      CoxModel m(n, p, time, ev.data(), X, nullptr, EFRON);
      if (m.n_events() == 0) { beta[k].assign(p, 0.0); continue; }
      const FitResult f = newton_fit(m, std::vector<double>(p, 0.0), 50, 1e-10);
      beta[k] = f.beta; out.flags |= f.flags; out.converged = out.converged && f.converged; out.loglik += f.ll;
      m.baseline(f.beta, nullptr, bt[k], bh[k]);
    }
    std::vector<double> grid;
    for (int k = 0; k < K; ++k) grid.insert(grid.end(), bt[k].begin(), bt[k].end());
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
    const std::size_t M = grid.size(), Ks = static_cast<std::size_t>(K);
    std::vector<double> inc(M * Ks, 0.0), F(M), r(Ks);
    for (std::size_t k = 0; k < Ks; ++k) {
      double prev = 0.0;
      for (std::size_t j = 0; j < bt[k].size(); ++j) {
        const std::size_t g = static_cast<std::size_t>(std::lower_bound(grid.begin(), grid.end(), bt[k][j]) - grid.begin());
        inc[g * Ks + k] = bh[k][j] - prev; prev = bh[k][j];
      }
    }
    for (std::size_t i = 0; i < nte; ++i) {
      for (std::size_t k = 0; k < Ks; ++k) {
        double lp = 0.0;
        for (std::size_t j = 0; j < p; ++j) lp += Xte[i * p + j] * beta[k][j];
        r[k] = std::exp(lp);
      }
      double S = 1.0, Fk = 0.0;
      for (std::size_t g = 0; g < M; ++g) {
        double tot = 0.0;
        for (std::size_t k = 0; k < Ks; ++k) tot += inc[g * Ks + k] * r[k];
        Fk += S * inc[g * Ks + static_cast<std::size_t>(interest - 1)] * r[static_cast<std::size_t>(interest - 1)];
        S *= std::fmax(0.0, 1.0 - tot);
        F[g] = Fk;
      }
      for (std::size_t t = 0; t < nt; ++t) out.risk[i * nt + t] = step_value(grid, F, times[t]);
    }
  }
  return out;
}

struct CompareResult {
  // models = direct choices, then Fine-Gray, then cause-specific Cox
  // metrics[((model * nt + t) * COMPARE_METRICS + metric) * 2 + (0 apparent | 1 cross-validated)]
  std::vector<double> metrics;
  // fit[model * 5 + {0 loglik, 1 npar, 2 aic, 3 flags, 4 number of folds whose fit did not converge}]
  std::vector<double> fit;
};

inline std::uint64_t splitmix64(std::uint64_t& s) {       // Steele, Lea & Flood (2014), OOPSLA
  std::uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

// Fold of every subject: a seeded Fisher-Yates shuffle, identical in every language binding.
inline std::vector<int> cv_folds(std::size_t n, int folds, std::uint64_t seed) {
  std::vector<std::size_t> perm(n);
  for (std::size_t i = 0; i < n; ++i) perm[i] = i;
  std::uint64_t s = seed;
  for (std::size_t i = n; i > 1; --i) { const std::size_t j = static_cast<std::size_t>(splitmix64(s) % i); std::swap(perm[i - 1], perm[j]); }
  std::vector<int> f(n);
  for (std::size_t q = 0; q < n; ++q) f[perm[q]] = static_cast<int>(q % static_cast<std::size_t>(folds));
  return f;
}

inline CompareResult compare_models(const std::vector<DirectChoice>& choices, std::size_t n, std::size_t p, int K,
                                    const double* time, const int* cause, const double* X, int interest,
                                    const std::vector<double>& times, int folds, std::uint64_t seed, int calib_groups = 10) {
  const std::size_t nt = times.size();
  const int COMPARE_MODELS = static_cast<int>(choices.size()) + 2;
  CompareResult res;
  res.metrics.assign(COMPARE_MODELS * nt * COMPARE_METRICS * 2, std::nan(""));
  res.fit.assign(COMPARE_MODELS * 5, std::nan(""));
  const std::vector<double> tv(time, time + n);
  const std::vector<int> cv(cause, cause + n);
  auto score = [&](int model, const std::vector<double>& risk, int which) {
    std::vector<double> pr(n);
    for (std::size_t t = 0; t < nt; ++t) {
      bool ok = true;
      for (std::size_t i = 0; i < n; ++i) { pr[i] = risk[i * nt + t]; ok = ok && std::isfinite(pr[i]); }
      if (!ok) continue;
      const BrierResult b = brier_score(tv, cv, interest, pr, times[t]);
      double* m = &res.metrics[((static_cast<std::size_t>(model) * nt + t) * COMPARE_METRICS) * 2];
      m[0 * 2 + which] = b.brier; m[1 * 2 + which] = b.ipa;
      m[2 * 2 + which] = auc_t(tv, cv, interest, pr, times[t]);
      m[3 * 2 + which] = calibration(tv, cv, interest, pr, times[t], calib_groups).error;
    }
  };
  for (int id = 0; id < COMPARE_MODELS; ++id) {
    const ModelPrediction a = fit_predict(id, choices, n, p, K, time, cause, X, interest, n, X, times);
    res.fit[id * 5 + 0] = a.loglik; res.fit[id * 5 + 1] = a.npar; res.fit[id * 5 + 2] = a.aic;
    res.fit[id * 5 + 3] = a.flags; res.fit[id * 5 + 4] = 0.0;
    score(id, a.risk, 0);
  }
  if (folds >= 2) {
    const std::vector<int> fold = cv_folds(n, folds, seed);
    std::vector<std::vector<double>> pooled(COMPARE_MODELS, std::vector<double>(n * nt, std::nan("")));
    for (int f = 0; f < folds; ++f) {
      std::vector<double> ttr, Xtr, Xte; std::vector<int> ctr; std::vector<std::size_t> te;
      for (std::size_t i = 0; i < n; ++i) {
        if (fold[i] == f) { te.push_back(i); Xte.insert(Xte.end(), X + i * p, X + (i + 1) * p); }
        else { ttr.push_back(time[i]); ctr.push_back(cause[i]); Xtr.insert(Xtr.end(), X + i * p, X + (i + 1) * p); }
      }
      for (int id = 0; id < COMPARE_MODELS; ++id) {
        const ModelPrediction a = fit_predict(id, choices, ttr.size(), p, K, ttr.data(), ctr.data(), Xtr.data(), interest,
                                              te.size(), Xte.data(), times);
        if (!a.converged) res.fit[id * 5 + 4] += 1.0;
        for (std::size_t q = 0; q < te.size(); ++q)
          for (std::size_t t = 0; t < nt; ++t) pooled[id][te[q] * nt + t] = a.risk[q * nt + t];
      }
    }
    for (int id = 0; id < COMPARE_MODELS; ++id) score(id, pooled[id], 1);
  }
  return res;
}

}  // namespace survcore
