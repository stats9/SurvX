// survcore: the numerical parts of adjusted risks (adjusted_risk() in R): a multinomial logistic
// model of the treatment, stabilised inverse probability weights with optional trimming, balance
// (standardised mean differences) and weight diagnostics, the weighted Aalen-Johansen estimator,
// bootstrap summaries, and the seeded resampling used by every bootstrap of the package.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "compare.hpp"
#include "linalg.hpp"
#include "summaries.hpp"

namespace survcore {

// ---------------------------------------------------------------------------------------------
// multinomial logistic regression of a (codes 0..K-1, level 0 the reference) on X (n x p,
// row-major, with its intercept column) by Newton-Raphson with step halving; returns the n x K
// fitted probabilities (row-major). Iterations as in the R implementation it replaces: at most 100,
// ridge 1e-8 on the information, a step is halved until the log-likelihood does not decrease by
// more than 1e-10 (or the step is below 1e-8 of a full step), stop when the change of the
// log-likelihood is below 1e-10 (|l| + 1).
// ---------------------------------------------------------------------------------------------
inline void mlogit_probs(std::size_t n, std::size_t p, std::size_t K, const double* X, const std::vector<double>& b,
                         std::vector<double>& P) {
  P.assign(n * K, 0.0);
  std::vector<double> e(K);
  for (std::size_t i = 0; i < n; ++i) {
    const double* x = X + i * p;
    e[0] = 0.0;
    for (std::size_t k = 1; k < K; ++k) {
      double v = 0.0;
      for (std::size_t r = 0; r < p; ++r) v += x[r] * b[(k - 1) * p + r];
      e[k] = v;
    }
    const double mx = *std::max_element(e.begin(), e.end());
    double sum = 0.0;
    for (std::size_t k = 0; k < K; ++k) { e[k] = std::exp(e[k] - mx); sum += e[k]; }
    for (std::size_t k = 0; k < K; ++k) P[i * K + k] = e[k] / sum;
  }
}

inline double mlogit_ll(std::size_t n, std::size_t K, const int* a, const std::vector<double>& P) {
  double ll = 0.0;
  for (std::size_t i = 0; i < n; ++i) ll += std::log(std::fmax(P[i * K + static_cast<std::size_t>(a[i])], 1e-300));
  return ll;
}

inline int mlogit_fit(std::size_t n, std::size_t p, std::size_t K, const double* X, const int* a, std::vector<double>& P,
                      int& iterations) {
  const std::size_t q = p * (K - 1);
  std::vector<double> b(q, 0.0), bn(q), g(q);
  mlogit_probs(n, p, K, X, b, P);
  double l0 = mlogit_ll(n, K, a, P), ln = l0;
  std::vector<double> Pn;
  iterations = 0;
  for (int it = 0; it < 100; ++it) {
    iterations = it + 1;
    std::fill(g.begin(), g.end(), 0.0);
    Mat H(q);
    for (std::size_t i = 0; i < n; ++i) {
      const double* x = X + i * p;
      for (std::size_t j = 1; j < K; ++j) {
        const double y = static_cast<std::size_t>(a[i]) == j ? 1.0 : 0.0;
        const double rj = y - P[i * K + j];
        for (std::size_t r = 0; r < p; ++r) g[(j - 1) * p + r] += x[r] * rj;
        for (std::size_t k = 1; k < K; ++k) {
          const double w = P[i * K + j] * ((j == k ? 1.0 : 0.0) - P[i * K + k]);
          if (w == 0.0) continue;
          for (std::size_t r = 0; r < p; ++r)
            for (std::size_t u = 0; u < p; ++u) H((j - 1) * p + r, (k - 1) * p + u) += x[r] * w * x[u];
        }
      }
    }
    for (std::size_t r = 0; r < q; ++r) H(r, r) += 1e-8;
    std::vector<double> step = g;
    if (!cholesky(H)) break;                       // the information is not positive definite: stop
    chol_solve(H, step);
    double s = 1.0;
    for (;;) {
      for (std::size_t r = 0; r < q; ++r) bn[r] = b[r] + step[r] * s;
      mlogit_probs(n, p, K, X, bn, Pn);
      ln = mlogit_ll(n, K, a, Pn);
      if (std::isfinite(ln) && ln >= l0 - 1e-10) break;
      s /= 2.0;
      if (s < 1e-8) break;
    }
    b = bn; P = Pn;
    if (std::fabs(ln - l0) < 1e-10 * (std::fabs(l0) + 1.0)) break;
    l0 = ln;
  }
  return 0;
}

// stabilised weights P(A = a) / P(A = a | x), trimmed to the [trim, 1 - trim] quantiles (type 7)
inline void ipw_weights(std::size_t n, std::size_t K, const double* P, const int* a, double trim, std::vector<double>& w) {
  std::vector<double> marg(K, 0.0);
  for (std::size_t i = 0; i < n; ++i) marg[static_cast<std::size_t>(a[i])] += 1.0;
  for (double& m : marg) m /= static_cast<double>(n);
  w.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) w[i] = marg[static_cast<std::size_t>(a[i])] / P[i * K + static_cast<std::size_t>(a[i])];
  if (trim > 0.0) {
    std::vector<double> s;
    for (double v : w) if (!std::isnan(v)) s.push_back(v);
    std::sort(s.begin(), s.end());
    const double lo = quantile7_sorted(s, trim), hi = quantile7_sorted(s, 1.0 - trim);
    for (double& v : w) v = std::fmin(std::fmax(v, lo), hi);
  }
}

// ---------------------------------------------------------------------------------------------
// weighted Aalen-Johansen estimate of the cumulative incidence of cause k at times tt (one event
// type: 1 - weighted Kaplan-Meier); risk set {entry < s <= time} (entry may be null)
// ---------------------------------------------------------------------------------------------
inline void weighted_cif(std::size_t n, const double* time, const int* status, const double* entry, const double* w, int k,
                         std::size_t nt, const double* tt, double* out) {
  std::vector<double> et;
  for (std::size_t i = 0; i < n; ++i) if (status[i] > 0) et.push_back(time[i]);
  std::sort(et.begin(), et.end());
  et.erase(std::unique(et.begin(), et.end()), et.end());
  if (et.empty()) { for (std::size_t j = 0; j < nt; ++j) out[j] = 0.0; return; }
  const std::size_t m = et.size();
  std::vector<double> Y(m, 0.0), dall(m, 0.0), dk(m, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    // at risk at every event time s with entry < s <= time
    const std::size_t hi = static_cast<std::size_t>(std::upper_bound(et.begin(), et.end(), time[i]) - et.begin());
    const std::size_t lo = entry ? static_cast<std::size_t>(std::upper_bound(et.begin(), et.end(), entry[i]) - et.begin()) : 0;
    for (std::size_t j = lo; j < hi; ++j) Y[j] += w[i];
    if (status[i] > 0) {
      const std::size_t j = static_cast<std::size_t>(std::lower_bound(et.begin(), et.end(), time[i]) - et.begin());
      dall[j] += w[i];
      if (status[i] == k) dk[j] += w[i];
    }
  }
  std::vector<double> F(m);
  double S = 1.0, Fc = 0.0;
  for (std::size_t j = 0; j < m; ++j) {
    const bool ok = Y[j] > 0.0;
    const double ha = ok ? dall[j] / Y[j] : 0.0, hk = ok ? dk[j] / Y[j] : 0.0;
    Fc += S * hk;
    S *= 1.0 - ha;
    F[j] = Fc;
  }
  for (std::size_t j = 0; j < nt; ++j) out[j] = step_at(et.data(), F.data(), m, tt[j], 0.0);
}

// ---------------------------------------------------------------------------------------------
// bootstrap summaries of a column of draws: standard deviation (n - 1) and the type-7 quantiles
// a, 1 - a. finite_only: non-finite draws are dropped from both; otherwise a non-finite draw
// makes the standard deviation NaN and is dropped from the quantiles only (NA values).
// ---------------------------------------------------------------------------------------------
inline void boot_column(std::size_t nb, const double* v, std::size_t stride, bool finite_only, double a, double& sd,
                        double& lo, double& hi) {
  std::vector<double> x, q;
  bool any_nan = false;
  for (std::size_t b = 0; b < nb; ++b) {
    const double y = v[b * stride];
    if (finite_only) { if (std::isfinite(y)) x.push_back(y); continue; }
    if (std::isnan(y)) any_nan = true;
    x.push_back(y);
  }
  for (double y : x) if (!std::isnan(y)) q.push_back(y);
  sd = kNaN;
  if (!any_nan && x.size() >= 2) {
    double m = 0.0;
    for (double y : x) m += y;
    m /= static_cast<double>(x.size());
    double ss = 0.0;
    for (double y : x) ss += (y - m) * (y - m);
    sd = std::sqrt(ss / static_cast<double>(x.size() - 1));
  }
  std::sort(q.begin(), q.end());
  lo = quantile7_sorted(q, a);
  hi = quantile7_sorted(q, 1.0 - a);
}

// ---------------------------------------------------------------------------------------------
// resampling and folds, from the index source of the core (compare.hpp: splitmix64 from the seed, or
// the generator of the front end); with R's generator row b is sample.int(n, n, replace = TRUE) - 1
// ---------------------------------------------------------------------------------------------
// B samples of size n drawn with replacement from 0..n-1 (row b of out)
inline void resample_indices(std::uint64_t seed, std::size_t n, std::size_t B, int* out) {
  std::uint64_t s = seed;
  for (std::size_t b = 0; b < B; ++b)
    for (std::size_t i = 0; i < n; ++i) out[b * n + i] = static_cast<int>(draw_index(s, n));
}

// folds 1..folds balanced within strata, strata in order of first appearance; with R's generator the
// folds of stratum st are sample(rep_len(1:folds, m)), as x[strata == st] <- ... in R
inline void stratified_folds(std::size_t n, const int* strata, int folds, std::uint64_t seed, int* out) {
  std::uint64_t s = seed;
  std::vector<int> seen;
  for (std::size_t i = 0; i < n; ++i)
    if (std::find(seen.begin(), seen.end(), strata[i]) == seen.end()) seen.push_back(strata[i]);
  for (int st : seen) {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < n; ++i) if (strata[i] == st) idx.push_back(i);
    const std::vector<std::size_t> perm = random_permutation(s, idx.size());
    for (std::size_t q = 0; q < idx.size(); ++q) out[idx[q]] = static_cast<int>(perm[q] % static_cast<std::size_t>(folds)) + 1;
  }
}

}  // namespace survcore
