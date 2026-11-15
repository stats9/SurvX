// survcore: K-sample tests.
//   Mantel (1966); Peto & Peto (1972); Harrington & Fleming (1982) Biometrika
//   69:553-566 -- log-rank test and the G^rho family, weight S(t-)^rho with S the
//   pooled Kaplan-Meier estimate.
//   Gray (1988) Ann Stat 16:1141-1154 -- K-sample test for the cumulative
//   incidence of a competing risk: scores (2.6) with weights K_k = L R_k,
//   L = [G_1^0(t-)]^rho from (3.2), variance (2.10) as implemented by the author
//   in cmprsk (see gray_test below).
// Stratified versions add score and variance contributions over strata, as
// described by Gray (1988) after Theorem 1.
#pragma once
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>
#include "linalg.hpp"
#include "special.hpp"

namespace survcore {

struct KSampleResult {
  std::vector<double> score;   // length K (log-rank: observed - expected)
  std::vector<double> expected;// log-rank only
  std::vector<double> observed;// log-rank only
  Mat var;                     // K x K covariance of the scores
  double statistic = NAN, p_value = NAN;
  int df = 0;
};

inline double chisq_upper(double x, double df) { return chisq_sf(x, df); }

// Indices of each stratum, sorted by time (strata == nullptr: one stratum).
inline std::vector<std::vector<std::size_t>> strata_orders(const std::vector<double>& time, const int* strata) {
  std::map<int, std::vector<std::size_t>> by;
  for (std::size_t i = 0; i < time.size(); ++i) by[strata ? strata[i] : 0].push_back(i);
  std::vector<std::vector<std::size_t>> out;
  for (auto& kv : by) {
    std::stable_sort(kv.second.begin(), kv.second.end(), [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
    out.push_back(std::move(kv.second));
  }
  return out;
}

inline void finish_ksample(KSampleResult& r, int K) {
  // Quadratic form on the first K-1 scores (they sum to zero).
  const std::size_t q = static_cast<std::size_t>(K - 1);
  Mat V(q);
  for (std::size_t a = 0; a < q; ++a) for (std::size_t b = 0; b < q; ++b) V(a, b) = r.var(a, b);
  std::vector<double> z(r.score.begin(), r.score.begin() + q), u = z;
  r.df = K - 1;
  if (q == 0 || !cholesky(V)) return;
  chol_solve(V, u);
  double x = 0.0;
  for (std::size_t a = 0; a < q; ++a) x += z[a] * u[a];
  r.statistic = x;
  r.p_value = chisq_sf(x, static_cast<double>(q));
}

// group: 0..K-1. status: 0 censored, nonzero event.
inline KSampleResult logrank_test(const std::vector<double>& time, const std::vector<int>& status,
                                  const std::vector<int>& group, int K, const int* strata, double rho) {
  KSampleResult r;
  r.score.assign(K, 0.0); r.expected.assign(K, 0.0); r.observed.assign(K, 0.0); r.var = Mat(K);
  for (const std::vector<std::size_t>& ord : strata_orders(time, strata)) {
  const std::size_t n = ord.size();
  std::vector<double> Y(K, 0.0), d(K), leave(K);
  for (std::size_t i = 0; i < n; ++i) Y[group[ord[i]]] += 1.0;
  double Yall = static_cast<double>(n), S = 1.0;
  for (std::size_t i = 0; i < n;) {
    const double t = time[ord[i]];
    std::fill(d.begin(), d.end(), 0.0); std::fill(leave.begin(), leave.end(), 0.0);
    double dall = 0.0;
    std::size_t j = i;
    for (; j < n && time[ord[j]] == t; ++j) {
      const std::size_t k = ord[j];
      leave[group[k]] += 1.0;
      if (status[k] != 0) { d[group[k]] += 1.0; dall += 1.0; }
    }
    if (dall > 0.0) {
      const double w = rho == 0.0 ? 1.0 : std::pow(S, rho);   // S(t-)
      const double vfac = Yall > 1.0 ? dall * (Yall - dall) / (Yall * Yall * (Yall - 1.0)) : 0.0;
      for (int a = 0; a < K; ++a) {
        const double e = dall * Y[a] / Yall;
        r.observed[a] += w * d[a]; r.expected[a] += w * e;
        for (int b = 0; b < K; ++b)
          r.var(a, b) += w * w * vfac * Y[a] * ((a == b ? Yall : 0.0) - Y[b]);
      }
      S *= 1.0 - dall / Yall;
    }
    for (int a = 0; a < K; ++a) Y[a] -= leave[a];
    Yall -= static_cast<double>(j - i);
    i = j;
  }
  }
  for (int a = 0; a < K; ++a) r.score[a] = r.observed[a] - r.expected[a];
  finish_ksample(r, K);
  return r;
}

// cause: 0 censored, `interest` = type 1, any other value = type 2 (competing).
//
// Scores: Gray (1988) eq. (2.6) with K_k = L R_k and L = [1 - F_1^0(t-)]^rho, eq. (3.2).
// Variance: eq. (2.10), evaluated exactly as in the author's own implementation
// (cmprsk, src/crstm.f, subroutine crst; Copyright (C) 2000 Robert Gray, GPL), of
// which this function is a port. Relative to a literal reading of the plug-in
// estimates listed after Theorem 1, that implementation
//   * uses d_kr = L H_k (1{k=r} - H_r / H.), i.e. K_k / G_1^0 with G_1k replaced by
//     its null value G_1^0 (H_r = Y_r / S_r(t-));
//   * evaluates G_1^0(t) / S_r(t) in b_jkr with right-continuous estimates;
//   * applies finite-population corrections for tied failures.
inline KSampleResult gray_test(const std::vector<double>& time, const std::vector<int>& cause,
                               const std::vector<int>& group, int K, const int* strata, int interest, double rho) {
  const std::size_t Ks = static_cast<std::size_t>(K);
  KSampleResult res;
  res.score.assign(K, 0.0); res.var = Mat(Ks);
  for (const std::vector<std::size_t>& ord : strata_orders(time, strata)) {
    const std::size_t n = ord.size();
    std::vector<double> rs(K, 0.0), f1m(K, 0.0), f1(K, 0.0), skmm(K, 1.0), skm(K, 1.0), v3(K, 0.0);
    std::vector<double> c(Ks * Ks, 0.0), a(Ks * Ks, 0.0), v2(Ks * Ks, 0.0), v(Ks * Ks, 0.0);
    std::vector<double> d0(K), d1(K), d2(K);
    for (std::size_t i = 0; i < n; ++i) rs[group[ord[i]]] += 1.0;
    double fm = 0.0, f = 0.0;
    for (std::size_t i = 0; i < n;) {
      const double t = time[ord[i]];
      std::fill(d0.begin(), d0.end(), 0.0); std::fill(d1.begin(), d1.end(), 0.0); std::fill(d2.begin(), d2.end(), 0.0);
      double nd1 = 0.0, nd2 = 0.0;
      std::size_t j = i;
      for (; j < n && time[ord[j]] == t; ++j) {
        const std::size_t q = ord[j];
        if (cause[q] == 0) d0[group[q]] += 1.0;
        else if (cause[q] == interest) { d1[group[q]] += 1.0; nd1 += 1.0; }
        else { d2[group[q]] += 1.0; nd2 += 1.0; }
      }
      if (nd1 + nd2 > 0.0) {
        double tr = 0.0, tq = 0.0;                         // sum_r H_r and sum_r R_r
        for (int g = 0; g < K; ++g) {
          if (rs[g] <= 0.0) continue;
          skm[g] = skmm[g] * (rs[g] - d1[g] - d2[g]) / rs[g];
          f1[g] = f1m[g] + skmm[g] * d1[g] / rs[g];
          tr += rs[g] / skmm[g];
          tq += rs[g] * (1.0 - f1m[g]) / skmm[g];
        }
        f = fm + nd1 / tr;
        const double fb = rho == 0.0 ? 1.0 : std::pow(1.0 - fm, rho);
        std::fill(a.begin(), a.end(), 0.0);
        for (std::size_t g = 0; g < Ks; ++g) {
          if (rs[g] <= 0.0) continue;
          const double t1 = rs[g] / skmm[g];
          for (std::size_t h = 0; h < Ks; ++h) {
            if (rs[h] <= 0.0) continue;
            a[g * Ks + h] = g == h ? fb * t1 * (1.0 - t1 / tr) : -fb * t1 * rs[h] / (skmm[h] * tr);
            c[g * Ks + h] += a[g * Ks + h] * nd1 / (tr * (1.0 - fm));
          }
        }
        for (int g = 0; g < K; ++g)
          if (rs[g] > 0.0) res.score[g] += fb * (d1[g] - nd1 * rs[g] * (1.0 - f1m[g]) / (skmm[g] * tq));
        if (nd1 > 0.0)
          for (std::size_t k = 0; k < Ks; ++k) {
            if (rs[k] <= 0.0) continue;
            const double t4 = skm[k] > 0.0 ? 1.0 - (1.0 - f) / skm[k] : 1.0;
            const double t5 = nd1 > 1.0 ? 1.0 - (nd1 - 1.0) / (tr * skmm[k] - 1.0) : 1.0;
            const double t3 = t5 * skmm[k] * nd1 / (tr * rs[k]);
            v3[k] += t4 * t4 * t3;
            for (std::size_t g = 0; g < Ks; ++g) {
              const double t1 = a[g * Ks + k] - t4 * c[g * Ks + k];
              v2[g * Ks + k] += t1 * t4 * t3;
              for (std::size_t h = 0; h <= g; ++h) v[g * Ks + h] += t1 * (a[h * Ks + k] - t4 * c[h * Ks + k]) * t3;
            }
          }
        if (nd2 > 0.0)
          for (std::size_t k = 0; k < Ks; ++k) {
            if (skm[k] <= 0.0 || d2[k] <= 0.0) continue;
            const double t4 = (1.0 - f) / skm[k];
            const double t5 = d2[k] > 1.0 ? 1.0 - (d2[k] - 1.0) / (rs[k] - 1.0) : 1.0;
            const double t3 = t5 * skmm[k] * skmm[k] * d2[k] / (rs[k] * rs[k]);
            v3[k] += t4 * t4 * t3;
            for (std::size_t g = 0; g < Ks; ++g) {
              const double t1 = t4 * c[g * Ks + k];
              v2[g * Ks + k] -= t1 * t4 * t3;
              for (std::size_t h = 0; h <= g; ++h) v[g * Ks + h] += t1 * t4 * c[h * Ks + k] * t3;
            }
          }
      }
      for (int g = 0; g < K; ++g) rs[g] -= d0[g] + d1[g] + d2[g];
      fm = f; f1m = f1; skmm = skm;
      i = j;
    }
    for (std::size_t g = 0; g < Ks; ++g)
      for (std::size_t h = 0; h <= g; ++h) {
        double x = v[g * Ks + h];
        for (std::size_t k = 0; k < Ks; ++k)
          x += c[g * Ks + k] * c[h * Ks + k] * v3[k] + c[g * Ks + k] * v2[h * Ks + k] + c[h * Ks + k] * v2[g * Ks + k];
        res.var(g, h) += x;
        if (h != g) res.var(h, g) += x;
      }
  }
  finish_ksample(res, K);
  return res;
}

}  // namespace survcore
