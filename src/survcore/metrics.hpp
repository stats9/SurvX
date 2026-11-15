// survcore: prediction performance for a predicted risk pi_i = F_1(t | z_i) of the
// event of interest by time t, with right censoring and competing risks.
//
//   Graf et al. (1999) Stat Med 18:2529-2545; Gerds & Schumacher (2006) Biom J
//   48:1029-1040; Schoop et al. (2011) Biometrics 67:39-49 -- IPCW Brier score.
//   Kattan & Gerds (2018) Diagn Progn Res 2:7 -- index of prediction accuracy (IPA).
//   Blanche, Dartigues & Jacqmin-Gadda (2013) Stat Med 32:5381-5397 -- time-dependent
//   AUC with competing risks; controls are all subjects without the event of interest
//   by t (event-free, or failed from a competing cause).
//   Gerds, Andersen & Kattan (2014) Stat Med 33:3191-3203 -- calibration in risk groups.
//
// Inverse probability of censoring weights (marginal Kaplan-Meier G of the censoring
// time, events preceding censorings at tied times):
//   subject with an event at X_i <= t : 1 / G(X_i-);  subject with X_i > t : 1 / G(t);
//   subject censored before t        : 0.
#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "nonparametric.hpp"

namespace survcore {

// IPCW weight of every subject at horizon t (see header). cause: 0 = censored.
inline std::vector<double> ipcw_weights(const std::vector<double>& time, const std::vector<int>& cause, double t) {
  const std::size_t n = time.size();
  std::vector<std::size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
  std::vector<double> w(n, 0.0);
  double G = 1.0, Gt = 1.0, at_risk = static_cast<double>(n);
  bool passed = false;
  for (std::size_t i = 0; i < n;) {
    const double u = time[ord[i]];
    if (u > t && !passed) { Gt = G; passed = true; }
    double ev = 0.0, cens = 0.0;
    std::size_t j = i;
    for (; j < n && time[ord[j]] == u; ++j) (cause[ord[j]] == 0 ? cens : ev) += 1.0;
    for (std::size_t q = i; q < j; ++q) {
      const std::size_t k = ord[q];
      if (u <= t && cause[k] != 0) w[k] = 1.0 / G;                 // G(X_i-)
    }
    if (cens > 0.0) G *= 1.0 - cens / (at_risk - ev);              // events leave before censorings
    at_risk -= static_cast<double>(j - i);
    i = j;
  }
  if (!passed) Gt = G;                                             // G(t) with t beyond the last time
  for (std::size_t k = 0; k < n; ++k) if (time[k] > t) w[k] = Gt > 0.0 ? 1.0 / Gt : 0.0;
  return w;
}

// Aalen-Johansen CIF of `interest` at time t (all other non-zero causes pooled).
inline double aj_at(const std::vector<double>& time, const std::vector<int>& cause, int interest, double t) {
  std::vector<int> c(cause.size());
  for (std::size_t i = 0; i < cause.size(); ++i) c[i] = cause[i] == 0 ? 0 : (cause[i] == interest ? 1 : 2);
  const AJResult aj = aalen_johansen(time, c, 2);
  double F = 0.0;
  for (std::size_t j = 0; j < aj.time.size() && aj.time[j] <= t; ++j) F = aj.cif[j * 2];
  return F;
}

struct BrierResult { double brier, brier_null, ipa; };

inline BrierResult brier_score(const std::vector<double>& time, const std::vector<int>& cause, int interest,
                               const std::vector<double>& pred, double t) {
  const std::size_t n = time.size();
  const std::vector<double> w = ipcw_weights(time, cause, t);
  const double f0 = aj_at(time, cause, interest, t);
  double bs = 0.0, b0 = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double y = (time[i] <= t && cause[i] == interest) ? 1.0 : 0.0;
    bs += w[i] * (y - pred[i]) * (y - pred[i]);
    b0 += w[i] * (y - f0) * (y - f0);
  }
  BrierResult r;
  r.brier = bs / static_cast<double>(n);
  r.brier_null = b0 / static_cast<double>(n);
  r.ipa = r.brier_null > 0.0 ? 1.0 - r.brier / r.brier_null : std::nan("");
  return r;
}

// Weighted probability that a case has a larger predicted risk than a control (ties count 1/2).
inline double auc_t(const std::vector<double>& time, const std::vector<int>& cause, int interest,
                    const std::vector<double>& pred, double t) {
  const std::size_t n = time.size();
  const std::vector<double> w = ipcw_weights(time, cause, t);
  std::vector<std::size_t> ord;
  for (std::size_t i = 0; i < n; ++i) if (w[i] > 0.0) ord.push_back(i);
  std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return pred[a] < pred[b]; });
  double wcase = 0.0, wctrl = 0.0, below = 0.0, num = 0.0;
  for (std::size_t i = 0; i < ord.size();) {
    std::size_t j = i;
    double tc = 0.0, tk = 0.0;                                      // case / control weight in this tie block
    for (; j < ord.size() && pred[ord[j]] == pred[ord[i]]; ++j) {
      const std::size_t k = ord[j];
      if (time[k] <= t && cause[k] == interest) tc += w[k]; else tk += w[k];
    }
    num += tc * (below + 0.5 * tk);
    below += tk; wcase += tc; wctrl += tk;
    i = j;
  }
  return (wcase > 0.0 && wctrl > 0.0) ? num / (wcase * wctrl) : std::nan("");
}

struct CalibrationResult {
  std::vector<double> mean_pred, observed, size;   // per risk group, ordered by predicted risk
  double error = 0.0;                              // sum_g (n_g / n) |mean_pred_g - observed_g|
};

// Groups of (nearly) equal size by rank of the predicted risk; observed = Aalen-Johansen within the group.
inline CalibrationResult calibration(const std::vector<double>& time, const std::vector<int>& cause, int interest,
                                     const std::vector<double>& pred, double t, int groups) {
  const std::size_t n = time.size();
  std::vector<std::size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return pred[a] < pred[b]; });
  CalibrationResult r;
  for (int g = 0; g < groups; ++g) {
    const std::size_t lo = n * static_cast<std::size_t>(g) / groups, hi = n * static_cast<std::size_t>(g + 1) / groups;
    if (hi <= lo) continue;
    std::vector<double> tt; std::vector<int> cc;
    double mp = 0.0;
    for (std::size_t q = lo; q < hi; ++q) { tt.push_back(time[ord[q]]); cc.push_back(cause[ord[q]]); mp += pred[ord[q]]; }
    mp /= static_cast<double>(hi - lo);
    const double obs = aj_at(tt, cc, interest, t);
    r.mean_pred.push_back(mp); r.observed.push_back(obs); r.size.push_back(static_cast<double>(hi - lo));
    r.error += static_cast<double>(hi - lo) / static_cast<double>(n) * std::fabs(mp - obs);
  }
  return r;
}

}  // namespace survcore
