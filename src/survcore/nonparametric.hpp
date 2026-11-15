// survcore: nonparametric estimators.
//   Kaplan & Meier (1958) JASA 53:457-481  -- product-limit estimator
//   Greenwood (1926)                        -- variance formula
//   Nelson (1972) Technometrics 14:945-966; Aalen (1978) Ann Statist 6:701-726
//                                           -- cumulative hazard and its variance
//   Irwin (1949) J Hygiene 47:188-189; Klein & Moeschberger (2003), Sec. 4.5
//                                           -- restricted mean survival time
//   Aalen & Johansen (1978) Scand J Stat 5:141-150 -- cumulative incidence
//   Kalbfleisch & Prentice (2002), Sec. 8.2 -- CIF as sum of S(t-) dLambda_k
//   Marubini & Valsecchi (1995), Analysing Survival Data from Clinical Trials
//   and Observational Studies, Sec. 10.3 -- delta-method variance of the CIF
#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace survcore {

struct KMResult {
  std::vector<double> time, n_risk, n_event, surv, std_err;  // at distinct event times
  std::vector<double> cumhaz, cumhaz_se;                     // Nelson-Aalen
};

// status: 0 = censored, nonzero = event. weights may be empty (all 1).
inline KMResult kaplan_meier(const std::vector<double>& time,
                             const std::vector<int>& status,
                             const std::vector<double>& weights) {
  const std::size_t n = time.size();
  std::vector<std::size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
  double at_risk = 0.0;
  for (std::size_t i = 0; i < n; ++i) at_risk += weights.empty() ? 1.0 : weights[i];
  KMResult r;
  double s = 1.0, gw = 0.0;  // gw: Greenwood sum d / (Y (Y - d))
  double na = 0.0, nav = 0.0;
  for (std::size_t i = 0; i < n;) {
    const double t = time[ord[i]];
    double d = 0.0, leave = 0.0;
    std::size_t j = i;
    for (; j < n && time[ord[j]] == t; ++j) {
      const double w = weights.empty() ? 1.0 : weights[ord[j]];
      leave += w;
      if (status[ord[j]] != 0) d += w;
    }
    if (d > 0.0) {
      s *= 1.0 - d / at_risk;
      gw += (at_risk > d) ? d / (at_risk * (at_risk - d)) : INFINITY;
      r.time.push_back(t);
      r.n_risk.push_back(at_risk);
      r.n_event.push_back(d);
      r.surv.push_back(s);
      // Greenwood's formula is undefined once S(t) = 0: report NaN rather than 0.
      r.std_err.push_back(s > 0.0 ? s * std::sqrt(gw) : std::nan(""));
      na += d / at_risk; nav += d / (at_risk * at_risk);
      r.cumhaz.push_back(na); r.cumhaz_se.push_back(std::sqrt(nav));
    }
    at_risk -= leave;
    i = j;
  }
  return r;
}

// Restricted mean survival time on [0, tau] from a Kaplan-Meier fit, with the
// Greenwood-type standard error  sum_j A_j^2 d_j / (Y_j (Y_j - d_j)),  A_j = int_{t_j}^{tau} S.
inline void rmst(const KMResult& km, double tau, double& mean, double& se) {
  const std::size_t m = km.time.size();
  std::size_t last = 0;
  while (last < m && km.time[last] <= tau) ++last;      // event times <= tau
  mean = 0.0;
  double prev_t = 0.0, prev_s = 1.0;
  for (std::size_t j = 0; j < last; ++j) { mean += prev_s * (km.time[j] - prev_t); prev_t = km.time[j]; prev_s = km.surv[j]; }
  mean += prev_s * (tau - prev_t);
  double var = 0.0, area = 0.0, right = tau;            // area = int_{t_j}^{tau} S
  for (std::size_t j = last; j-- > 0;) {
    area += km.surv[j] * (right - km.time[j]);
    right = km.time[j];
    const double Y = km.n_risk[j], d = km.n_event[j];
    if (Y > d) var += area * area * d / (Y * (Y - d));
  }
  se = std::sqrt(var);
}

struct AJResult {
  std::vector<double> time;   // distinct times with at least one event of any cause
  std::vector<double> cif;    // row-major, time.size() x K
  std::vector<double> surv;   // overall event-free survival
  std::vector<double> std_err;// row-major, time.size() x K (delta method)
};

// cause: 0 = censored, 1..K = cause of first event.
inline AJResult aalen_johansen(const std::vector<double>& time,
                               const std::vector<int>& cause, int K) {
  const std::size_t n = time.size();
  std::vector<std::size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
  AJResult r;
  double at_risk = static_cast<double>(n), s = 1.0;
  std::vector<double> F(K, 0.0), d(K);
  // running sums for the variance, per cause: B, C, D, E, Fq; A is common
  double A = 0.0;
  std::vector<double> vB(K, 0.0), vC(K, 0.0), vD(K, 0.0), vE(K, 0.0), vF(K, 0.0);
  for (std::size_t i = 0; i < n;) {
    const double t = time[ord[i]];
    std::fill(d.begin(), d.end(), 0.0);
    double dall = 0.0;
    std::size_t j = i;
    for (; j < n && time[ord[j]] == t; ++j) {
      const int c = cause[ord[j]];
      if (c >= 1 && c <= K) { d[c - 1] += 1.0; dall += 1.0; }
    }
    if (dall > 0.0) {
      for (int k = 0; k < K; ++k) F[k] += s * d[k] / at_risk;   // S(t-) dLambda_k(t)
      const double Y = at_risk, gj = (Y > dall) ? dall / (Y * (Y - dall)) : 0.0;
      A += gj;
      for (int k = 0; k < K; ++k) {
        const double u = s * d[k] / (Y * Y);                    // s is still S(t-)
        vB[k] += F[k] * gj;
        vC[k] += F[k] * F[k] * gj;
        vD[k] += s * s * d[k] * (Y - d[k]) / (Y * Y * Y);
        vE[k] += u;
        vF[k] += F[k] * u;
        const double v = F[k] * F[k] * A - 2.0 * F[k] * vB[k] + vC[k] + vD[k] - 2.0 * F[k] * vE[k] + 2.0 * vF[k];
        r.std_err.push_back(v > 0.0 ? std::sqrt(v) : 0.0);
      }
      s *= 1.0 - dall / at_risk;
      r.time.push_back(t);
      r.surv.push_back(s);
      r.cif.insert(r.cif.end(), F.begin(), F.end());
    }
    at_risk -= static_cast<double>(j - i);
    i = j;
  }
  return r;
}

}  // namespace survcore
