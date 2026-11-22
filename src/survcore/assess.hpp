// survcore: further measures of the performance of predicted risks, and of risk scores.
//
//   Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA (1982) JAMA 247:2543-2546, and Harrell,
//     Lee & Mark (1996) Stat Med 15:361-387 -- concordance index C.
//   Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011) Stat Med 30:1105-1117 -- C with
//     inverse probability of censoring weights 1/G(T_i-)^2, restricted to event times <= tau.
//   Wolbers M, Koller MT, Witteman JCM, Steyerberg EW (2009) Epidemiology 20:555-561 -- C for one
//     cause with competing risks: a subject who failed from a competing cause stays comparable
//     (as a control) with every case.
//   Therneau TM, Atkinson EJ (2024) Concordance, vignette of the survival package -- the
//     infinitesimal jackknife variance used here (var = sum of squared dfbetas).
//   Blanche P, Dartigues JF, Jacqmin-Gadda H (2013) Stat Med 32:5381-5397 -- influence functions
//     of the IPCW AUC with known censoring weights (the "conservative" variance of
//     riskRegression::Score); Brier score: Gerds & Schumacher (2006) Biom J 48:1029-1040.
//   Austin PC, Harrell FE, van Klaveren D (2020) Stat Med 39:2714-2742 -- smooth calibration
//     curve and the integrated calibration index ICI, E50, E90: a Cox model (or, with competing
//     risks, a Fine-Gray model) for the outcome with a restricted cubic spline (3 knots at the 10th,
//     50th and 90th percentiles) of the complementary log-log of the predicted risk.
//   Vickers AJ, Elkin EB (2006) Med Decis Making 26:565-574; Vickers AJ, Cronin AM, Elkin EB,
//     Gonen M (2008) BMC Med Inform Decis Mak 8:53 -- decision curve analysis, net benefit with
//     censored data (Kaplan-Meier or Aalen-Johansen among subjects above the threshold).
#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "cox.hpp"
#include "finegray.hpp"
#include "metrics.hpp"
#include "newton.hpp"
#include "nonparametric.hpp"
#include "summaries.hpp"
#include "timeeffects.hpp"

namespace survcore {

namespace assess_detail {

// Fenwick (binary indexed) tree of weights over 0..m-1
class Fenwick {
 public:
  explicit Fenwick(std::size_t m) : t_(m + 1, 0.0) {}
  void add(std::size_t i, double v) { for (++i; i < t_.size(); i += i & (~i + 1)) t_[i] += v; }
  double prefix(std::size_t i) const {                    // sum of the weights of 0..i-1
    double s = 0.0;
    for (; i > 0; i -= i & (~i + 1)) s += t_[i];
    return s;
  }
  double total() const { return prefix(t_.size() - 1); }
 private:
  std::vector<double> t_;
};

// value of a right-continuous step function with jumps at increasing x (0 before the first)
inline double step0(const std::vector<double>& x, const std::vector<double>& y, double at) {
  const std::size_t k = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), at) - x.begin());
  return k == 0 ? 0.0 : y[k - 1];
}

// rank (0-based, ties share the rank) of every value among the distinct values; m = number of ranks
inline std::vector<std::size_t> dense_ranks(const double* v, std::size_t n, std::size_t& m) {
  std::vector<double> u(v, v + n);
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  m = u.size();
  std::vector<std::size_t> r(n);
  for (std::size_t i = 0; i < n; ++i) r[i] = static_cast<std::size_t>(std::lower_bound(u.begin(), u.end(), v[i]) - u.begin());
  return r;
}

// G(t-) of the censoring distribution at every subject's own time (Kaplan-Meier of the censoring
// times; at tied times events come before censorings), as in ipcw_weights()
inline std::vector<double> g_minus(const double* time, const int* cause, std::size_t n) {
  std::vector<std::size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
  std::vector<double> g(n, 1.0);
  double G = 1.0, at_risk = static_cast<double>(n);
  for (std::size_t i = 0; i < n;) {
    const double u = time[ord[i]];
    double ev = 0.0, cens = 0.0;
    std::size_t j = i;
    for (; j < n && time[ord[j]] == u; ++j) (cause[ord[j]] == 0 ? cens : ev) += 1.0;
    for (std::size_t q = i; q < j; ++q) g[ord[q]] = G;
    if (cens > 0.0) G *= 1.0 - cens / (at_risk - ev);
    at_risk -= static_cast<double>(j - i);
    i = j;
  }
  return g;
}

inline double sd_over_sqrt_n(const std::vector<double>& x) {   // sd(x) / sqrt(n), sd with n - 1
  const std::size_t n = x.size();
  if (n < 2) return std::nan("");
  double m = 0.0;
  for (double v : x) m += v;
  m /= static_cast<double>(n);
  double s = 0.0;
  for (double v : x) s += (v - m) * (v - m);
  return std::sqrt(s / static_cast<double>(n - 1)) / std::sqrt(static_cast<double>(n));
}

}  // namespace assess_detail

// ------------------------------------------------------------------------------------------------
// concordance index
// ------------------------------------------------------------------------------------------------
// A case is a subject with an event of `interest` at a time <= tau. Its controls are the subjects
// still event-free after it (time larger, or the same time and censored) and, with competing
// risks, every subject who failed from another cause (Wolbers et al. 2009). Events of interest at
// the same time are not compared. A pair is concordant when the case has the higher risk score;
// tied scores count 1/2. Each pair has the weight of its case: 1 (Harrell) or 1/G(T_i-)^2 (Uno).
// C = sum of weight x (1, 1/2 or 0) / sum of weights; standard error from the infinitesimal
// jackknife: dfbeta_u = (N_u - C D_u) / D, with N_u and D_u the weighted concordance and the
// weight of the pairs that contain subject u (in either role); var = sum dfbeta_u^2.
struct ConcordanceResult {
  double c = std::nan(""), se = std::nan("");
  double concordant = 0.0, discordant = 0.0, tied = 0.0;   // weighted counts of the pairs
  double pairs = 0.0;                                       // unweighted number of comparable pairs
};

inline ConcordanceResult concordance_index(std::size_t n, const double* time, const int* cause, int interest,
                                           const double* risk, double tau, bool uno, double* dfbeta = nullptr) {
  using namespace assess_detail;
  ConcordanceResult res;
  std::size_t R = 0;
  const std::vector<std::size_t> rk = dense_ranks(risk, n, R);
  const std::vector<double> gm = g_minus(time, cause, n);
  std::vector<double> a(n, 0.0);                           // case weights
  for (std::size_t i = 0; i < n; ++i)
    if (cause[i] == interest && time[i] <= tau) a[i] = uno ? (gm[i] > 0.0 ? 1.0 / (gm[i] * gm[i]) : 0.0) : 1.0;
  auto competing = [&](std::size_t i) { return cause[i] != 0 && cause[i] != interest; };
  std::vector<std::size_t> asc(n);
  std::iota(asc.begin(), asc.end(), 0);
  std::stable_sort(asc.begin(), asc.end(), [&](std::size_t x, std::size_t y) { return time[x] < time[y]; });
  std::vector<double> Nu(n, 0.0), Du(n, 0.0);              // per subject: weighted concordance and weight of its pairs
  double cnt_lo = 0.0, cnt_eq = 0.0, cnt_hi = 0.0, npairs = 0.0;
  // the controls of every case: subjects at risk after it (descending pass) ...
  {
    Fenwick F(R), Fc(R);                                    // weights 1; Fc counts each subject once (for the pair count)
    std::size_t j = n;
    while (j > 0) {
      std::size_t i = j;
      const double u = time[asc[j - 1]];
      while (i > 0 && time[asc[i - 1]] == u) --i;          // block asc[i..j-1] has time u
      for (std::size_t q = i; q < j; ++q) if (cause[asc[q]] == 0) F.add(rk[asc[q]], 1.0);
      for (std::size_t q = i; q < j; ++q) {
        const std::size_t k = asc[q];
        if (a[k] == 0.0) continue;
        const double lo = F.prefix(rk[k]), eq = F.prefix(rk[k] + 1) - lo, tot = F.total(), hi = tot - lo - eq;
        Nu[k] += a[k] * (lo + 0.5 * eq); Du[k] += a[k] * tot;
        cnt_lo += a[k] * lo; cnt_eq += a[k] * eq; cnt_hi += a[k] * hi; npairs += tot;
      }
      for (std::size_t q = i; q < j; ++q) if (cause[asc[q]] != 0) F.add(rk[asc[q]], 1.0);
      j = i;
    }
  }
  // ... and the subjects who failed from a competing cause at or before it (ascending pass)
  {
    Fenwick F(R);
    std::size_t i = 0;
    while (i < n) {
      std::size_t j = i;
      const double u = time[asc[i]];
      while (j < n && time[asc[j]] == u) ++j;
      for (std::size_t q = i; q < j; ++q) if (competing(asc[q])) F.add(rk[asc[q]], 1.0);
      for (std::size_t q = i; q < j; ++q) {
        const std::size_t k = asc[q];
        if (a[k] == 0.0) continue;
        const double lo = F.prefix(rk[k]), eq = F.prefix(rk[k] + 1) - lo, tot = F.total(), hi = tot - lo - eq;
        Nu[k] += a[k] * (lo + 0.5 * eq); Du[k] += a[k] * tot;
        cnt_lo += a[k] * lo; cnt_eq += a[k] * eq; cnt_hi += a[k] * hi; npairs += tot;
      }
      i = j;
    }
  }
  // every subject as a control: the weighted cases it is compared with (ascending pass)
  {
    Fenwick C(R);
    std::size_t i = 0;
    auto as_control = [&](std::size_t k) {
      const double tot = C.total(), le = C.prefix(rk[k] + 1), eq = le - C.prefix(rk[k]), gt = tot - le;
      Nu[k] += gt + 0.5 * eq; Du[k] += tot;
    };
    while (i < n) {
      std::size_t j = i;
      const double u = time[asc[i]];
      while (j < n && time[asc[j]] == u) ++j;
      for (std::size_t q = i; q < j; ++q) if (cause[asc[q]] == interest) as_control(asc[q]);   // cases strictly before u
      for (std::size_t q = i; q < j; ++q) if (a[asc[q]] > 0.0) C.add(rk[asc[q]], a[asc[q]]);
      for (std::size_t q = i; q < j; ++q) if (cause[asc[q]] == 0) as_control(asc[q]);          // cases up to u
      i = j;
    }
    for (std::size_t k = 0; k < n; ++k) if (competing(k)) as_control(k);                       // every case
  }
  const double den = cnt_lo + cnt_eq + cnt_hi;
  res.concordant = cnt_lo; res.discordant = cnt_hi; res.tied = cnt_eq; res.pairs = npairs;
  if (!(den > 0.0)) return res;
  res.c = (cnt_lo + 0.5 * cnt_eq) / den;
  double v = 0.0;
  for (std::size_t k = 0; k < n; ++k) {
    const double d = (Nu[k] - res.c * Du[k]) / den;
    if (dfbeta) dfbeta[k] = d;
    v += d * d;
  }
  res.se = std::sqrt(v);
  return res;
}

// ------------------------------------------------------------------------------------------------
// influence functions of the IPCW AUC and Brier score at t with the censoring weights treated as
// known (conservative; riskRegression::Score(..., conservative = TRUE)); se = sd(IF) / sqrt(n)
// ------------------------------------------------------------------------------------------------
inline void score_influence(const std::vector<double>& time, const std::vector<int>& cause, int interest,
                            const std::vector<double>& pred, double t, std::vector<double>& if_auc, std::vector<double>& if_brier,
                            double& auc, double& brier) {
  const std::size_t n = time.size();
  const std::vector<double> w = ipcw_weights(time, cause, t);
  std::vector<double> cw(n, 0.0), vw(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    if (w[i] <= 0.0) continue;
    if (time[i] <= t && cause[i] == interest) cw[i] = w[i]; else vw[i] = w[i];
  }
  const double nn = static_cast<double>(n);
  double sc = 0.0, sv = 0.0;
  for (std::size_t i = 0; i < n; ++i) { sc += cw[i]; sv += vw[i]; }
  // A_i = (1/n) sum_j vw_j s(r_i, r_j), B_i = (1/n) sum_j cw_j s(r_j, r_i); s = 1 if first larger, 1/2 if tied
  std::vector<std::size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](std::size_t x, std::size_t y) { return pred[x] < pred[y]; });
  std::vector<double> A(n), B(n);
  double below_v = 0.0, below_c = 0.0;
  for (std::size_t i = 0; i < n;) {
    std::size_t j = i;
    double bv = 0.0, bc = 0.0;
    for (; j < n && pred[ord[j]] == pred[ord[i]]; ++j) { bv += vw[ord[j]]; bc += cw[ord[j]]; }
    for (std::size_t q = i; q < j; ++q) {
      A[ord[q]] = (below_v + 0.5 * bv) / nn;
      B[ord[q]] = ((sc - below_c - bc) + 0.5 * bc) / nn;
    }
    below_v += bv; below_c += bc;
    i = j;
  }
  const double Wbar = sc / nn, Vbar = sv / nn, mu = Wbar * Vbar;
  double nu = 0.0;
  for (std::size_t i = 0; i < n; ++i) nu += cw[i] * A[i];
  nu /= nn;
  if_auc.assign(n, std::nan(""));
  auc = mu > 0.0 ? nu / mu : std::nan("");
  if (mu > 0.0)
    for (std::size_t i = 0; i < n; ++i)
      if_auc[i] = ((cw[i] * A[i] + vw[i] * B[i]) * mu - nu * (cw[i] * Vbar + vw[i] * Wbar)) / (mu * mu);
  brier = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double y = (time[i] <= t && cause[i] == interest) ? 1.0 : 0.0;
    brier += w[i] * (y - pred[i]) * (y - pred[i]);
  }
  brier /= nn;
  if_brier.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double y = (time[i] <= t && cause[i] == interest) ? 1.0 : 0.0;
    if_brier[i] = w[i] * (y - pred[i]) * (y - pred[i]) - brier;
  }
}

// ------------------------------------------------------------------------------------------------
// smooth calibration at t (Austin, Harrell & van Klaveren 2020)
// ------------------------------------------------------------------------------------------------
// The predicted risks are clamped to [1e-10, 1 - 1e-10] before the complementary log-log. observed
// receives the calibrated risk of every subject. Returns false (and NaN measures) when the three
// knots are not distinct or the model cannot be fitted.
struct SmoothCalibration { double ici = std::nan(""), e50 = std::nan(""), e90 = std::nan(""); bool ok = false; };

inline SmoothCalibration smooth_calibration(const std::vector<double>& time, const std::vector<int>& cause, int interest,
                                            const std::vector<double>& pred, double t, std::vector<double>* observed = nullptr) {
  using namespace assess_detail;
  SmoothCalibration out;
  const std::size_t n = time.size();
  if (observed) observed->assign(n, std::nan(""));
  if (n < 5) return out;
  std::vector<double> x(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double r = std::fmin(std::fmax(pred[i], 1e-10), 1.0 - 1e-10);
    x[i] = std::log(-std::log1p(-r));
  }
  std::vector<double> xs(x);
  std::sort(xs.begin(), xs.end());
  const double k1 = quantile7_sorted(xs, 0.1), k2 = quantile7_sorted(xs, 0.5), k3 = quantile7_sorted(xs, 0.9);
  if (!(k1 < k2 && k2 < k3)) return out;
  const NaturalSpline ns = ns_setup(std::vector<double>{k2}, k1, k3);
  const std::size_t p = 2;
  std::vector<double> X(n * p);
  for (std::size_t i = 0; i < n; ++i) {
    const std::vector<double> b = ns_row(ns, x[i]);
    if (b.size() != p) return out;
    X[i * p] = b[0]; X[i * p + 1] = b[1];
  }
  bool comp = false; std::size_t events = 0;
  for (std::size_t i = 0; i < n; ++i) { comp = comp || (cause[i] != 0 && cause[i] != interest); events += cause[i] == interest; }
  if (events < 3) return out;
  std::vector<double> beta, et, H;
  if (!comp) {
    std::vector<int> st(n);
    for (std::size_t i = 0; i < n; ++i) st[i] = cause[i] == interest;
    CoxModel m(n, p, time.data(), st.data(), X.data(), nullptr, EFRON);
    const FitResult f = newton_fit(m, std::vector<double>(p, 0.0), 50, 1e-10);
    if (!std::isfinite(f.ll) || !f.converged) return out;
    beta = f.beta;
    m.baseline(beta, nullptr, et, H);
  } else {
    FineGrayModel m(n, p, time.data(), cause.data(), X.data(), interest);
    if (m.n_events() == 0) return out;
    const FitResult f = newton_fit(m, std::vector<double>(p, 0.0), 50, 1e-10);
    if (!std::isfinite(f.ll) || !f.converged) return out;
    beta = f.beta;
    m.baseline(beta, nullptr, H);
    et = m.event_times();
  }
  const double H0 = step0(et, H, t);
  std::vector<double> d(n);
  double s = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double o = -std::expm1(-H0 * std::exp(X[i * p] * beta[0] + X[i * p + 1] * beta[1]));
    if (observed) (*observed)[i] = o;
    d[i] = std::fabs(o - pred[i]);
    s += d[i];
  }
  std::sort(d.begin(), d.end());
  out.ici = s / static_cast<double>(n);
  out.e50 = quantile7_sorted(d, 0.5);
  out.e90 = quantile7_sorted(d, 0.9);
  out.ok = true;
  return out;
}

// observed (Kaplan-Meier or Aalen-Johansen) over expected (mean predicted) risk at t, with the
// interval exp(log(O/E) -/+ z se(O)/O) from the standard error of the non-parametric estimate
inline void observed_expected(const std::vector<double>& time, const std::vector<int>& cause, int interest,
                              const std::vector<double>& pred, double t, double z, double& oe, double& lower, double& upper) {
  const std::size_t n = time.size();
  std::vector<int> c(n);
  for (std::size_t i = 0; i < n; ++i) c[i] = cause[i] == 0 ? 0 : (cause[i] == interest ? 1 : 2);
  const AJResult aj = aalen_johansen(time, c, 2);
  double O = 0.0, se = std::nan("");
  for (std::size_t j = 0; j < aj.time.size() && aj.time[j] <= t; ++j) { O = aj.cif[j * 2]; se = aj.std_err[j * 2]; }
  double E = 0.0;
  for (double v : pred) E += v;
  E /= static_cast<double>(n);
  oe = E > 0.0 ? O / E : std::nan("");
  lower = upper = std::nan("");
  if (std::isfinite(oe) && O > 0.0 && std::isfinite(se)) {
    lower = oe * std::exp(-z * se / O);
    upper = oe * std::exp(z * se / O);
  }
}

// ------------------------------------------------------------------------------------------------
// decision curve: net benefit of treating the subjects whose risk is at least the threshold p,
// TP/n - FP/n x p/(1-p), with TP/n = pi F_+(t) and FP/n = pi (1 - F_+(t)), pi the share above the
// threshold and F_+ the Kaplan-Meier (Aalen-Johansen) risk among them; "treat all": F(t) - (1 -
// F(t)) p/(1-p)
// ------------------------------------------------------------------------------------------------
inline void net_benefit(const std::vector<double>& time, const std::vector<int>& cause, int interest,
                        const std::vector<double>& pred, double t, const std::vector<double>& thr,
                        std::vector<double>& nb, std::vector<double>& nb_all) {
  const std::size_t n = time.size(), m = thr.size();
  nb.assign(m, std::nan("")); nb_all.assign(m, std::nan(""));
  const double F = aj_at(time, cause, interest, t);
  for (std::size_t k = 0; k < m; ++k) {
    const double p = thr[k];
    if (!(p >= 0.0 && p < 1.0)) continue;
    const double odds = p / (1.0 - p);
    nb_all[k] = F - (1.0 - F) * odds;
    std::vector<double> ts; std::vector<int> cs;
    for (std::size_t i = 0; i < n; ++i) if (pred[i] >= p) { ts.push_back(time[i]); cs.push_back(cause[i]); }
    if (ts.empty()) { nb[k] = 0.0; continue; }
    const double pi = static_cast<double>(ts.size()) / static_cast<double>(n);
    const double Fp = aj_at(ts, cs, interest, t);
    nb[k] = pi * Fp - pi * (1.0 - Fp) * odds;
  }
}

// ------------------------------------------------------------------------------------------------
// all measures of M models at nt times, and the paired comparison of every model with the one of
// lowest Brier score at that time. risk[(m * nt + j) * n + i] (one n x nt column-major
// matrix per model). out[((m * nt + j) * kAssessCols) + c]:
//   0 brier, 1 brier_null, 2 ipa, 3 auc, 4 calibration error (risk groups), 5 se(brier),
//   6 se(auc), 7 O/E, 8 lower, 9 upper, 10 ICI, 11 E50, 12 E90,
//   13 brier - brier of the reference, 14 its se, 15 its p-value,
//   16 auc - auc of the reference, 17 its se, 18 its p-value.
// ref[j]: the reference model (0-based) at time j, -1 if none.
// ------------------------------------------------------------------------------------------------
constexpr int kAssessCols = 19;

inline void assess_models(std::size_t n, const double* time, const int* cause, int interest, std::size_t nt, const double* times,
                          std::size_t M, const double* risk, int groups, double z, double* out, int* ref) {
  using namespace assess_detail;
  const std::vector<double> tv(time, time + n);
  const std::vector<int> cv(cause, cause + n);
  const double nan = std::nan("");
  for (std::size_t k = 0; k < M * nt * kAssessCols; ++k) out[k] = nan;
  std::vector<std::vector<double>> ifa(M), ifb(M);
  for (std::size_t j = 0; j < nt; ++j) {
    const double t = times[j];
    std::vector<char> okm(M, 0);
    for (std::size_t m = 0; m < M; ++m) {
      std::vector<double> pr(n);
      bool ok = true;
      for (std::size_t i = 0; i < n; ++i) { pr[i] = risk[(m * nt + j) * n + i]; ok = ok && std::isfinite(pr[i]); }
      if (!ok) { ifa[m].clear(); ifb[m].clear(); continue; }
      okm[m] = 1;
      double* o = out + (m * nt + j) * kAssessCols;
      const BrierResult b = brier_score(tv, cv, interest, pr, t);
      o[0] = b.brier; o[1] = b.brier_null; o[2] = b.ipa;
      o[3] = auc_t(tv, cv, interest, pr, t);
      o[4] = calibration(tv, cv, interest, pr, t, groups).error;
      double a_if = 0.0, b_if = 0.0;
      score_influence(tv, cv, interest, pr, t, ifa[m], ifb[m], a_if, b_if);
      o[5] = sd_over_sqrt_n(ifb[m]);
      o[6] = std::isfinite(o[3]) ? sd_over_sqrt_n(ifa[m]) : nan;
      observed_expected(tv, cv, interest, pr, t, z, o[7], o[8], o[9]);
      const SmoothCalibration sc = smooth_calibration(tv, cv, interest, pr, t);
      o[10] = sc.ici; o[11] = sc.e50; o[12] = sc.e90;
    }
    int best = -1;
    for (std::size_t m = 0; m < M; ++m)
      if (okm[m] && std::isfinite(out[(m * nt + j) * kAssessCols]) &&
          (best < 0 || out[(m * nt + j) * kAssessCols] < out[(static_cast<std::size_t>(best) * nt + j) * kAssessCols]))
        best = static_cast<int>(m);
    ref[j] = best;
    if (best < 0) continue;
    const std::size_t r = static_cast<std::size_t>(best);
    for (std::size_t m = 0; m < M; ++m) {
      if (!okm[m]) continue;
      double* o = out + (m * nt + j) * kAssessCols;
      const double* orf = out + (r * nt + j) * kAssessCols;
      o[13] = o[0] - orf[0]; o[16] = o[3] - orf[3];
      if (m == r) continue;
      std::vector<double> d(n);
      for (std::size_t i = 0; i < n; ++i) d[i] = ifb[m][i] - ifb[r][i];
      o[14] = sd_over_sqrt_n(d);
      o[15] = o[14] > 0.0 ? 2.0 * pnorm_std(-std::fabs(o[13]) / o[14]) : nan;
      if (std::isfinite(o[3]) && std::isfinite(orf[3])) {
        for (std::size_t i = 0; i < n; ++i) d[i] = ifa[m][i] - ifa[r][i];
        o[17] = sd_over_sqrt_n(d);
        o[18] = o[17] > 0.0 ? 2.0 * pnorm_std(-std::fabs(o[16]) / o[17]) : nan;
      }
    }
  }
}

}  // namespace survcore
