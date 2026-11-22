// survcore: score test of proportional hazards of a Cox model (Grambsch & Therneau 1994), per
// covariate column, per term (groups of columns) and globally, with the scaled Schoenfeld residuals
// for the plot (ph_test() and time_effects() in R; computed as survival::cox.zph does).
//
// For every event (status 1), ordered by stratum and time: the risk-set weighted mean xbar(t) and
// covariance V(t) of the covariates with weights w_j exp(x_j' beta) over {same stratum, entry < t
// <= time}; Schoenfeld residuals s = x - xbar(t); a transform g of the event times centred over the
// events; U = sum s g, A = sum V g^2 - (sum V g)' (sum V)^-1 (sum V g); chi-square U_j^2 / A_jj per
// column, U_J' A_JJ^-1 U_J per term J and U' A^-1 U globally. Without delayed entry and strata the
// risk-set sums are reverse cumulative sums in time order; otherwise they are summed directly.
#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "linalg.hpp"
#include "summaries.hpp"

namespace survcore {

enum PhTransform : int { PHT_IDENTITY = 0, PHT_LOG = 1, PHT_RANK = 2, PHT_KM = 3, PHT_KM_WEIGHTED = 4 };

struct RiskSets {
  std::vector<std::size_t> rows;          // event rows, by stratum then time
  std::vector<double> time, xbar, V;      // xbar: nd x p, V: nd x p x p
};

inline RiskSets cox_risksets(std::size_t n, std::size_t p, const double* time, const int* status, const double* X,
                             const double* beta, const double* w, const double* entry, const int* strata) {
  RiskSets rs;
  std::vector<double> lp(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) for (std::size_t r = 0; r < p; ++r) lp[i] += X[i * p + r] * beta[r];
  std::vector<std::size_t> ev;
  for (std::size_t i = 0; i < n; ++i) if (status[i] == 1) ev.push_back(i);
  if (!entry && !strata) {
    std::vector<double> rr(n);
    for (std::size_t i = 0; i < n; ++i) rr[i] = (w ? w[i] : 1.0) * std::exp(lp[i]);
    std::vector<std::size_t> o(n);
    std::iota(o.begin(), o.end(), 0);
    std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
    std::stable_sort(ev.begin(), ev.end(), [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
    // reverse cumulative sums in time order: S0, S1 (p), S2 (p x p, lower triangle)
    std::vector<double> S0(n), S1(n * p), S2(n * p * p);
    double s0 = 0.0;
    std::vector<double> s1(p, 0.0), s2(p * p, 0.0);
    for (std::size_t q = n; q-- > 0;) {
      const std::size_t i = o[q];
      s0 += rr[i];
      for (std::size_t a = 0; a < p; ++a) {
        s1[a] += X[i * p + a] * rr[i];
        for (std::size_t b = 0; b <= a; ++b) s2[a * p + b] += X[i * p + a] * X[i * p + b] * rr[i];
      }
      S0[q] = s0;
      for (std::size_t a = 0; a < p; ++a) { S1[q * p + a] = s1[a]; for (std::size_t b = 0; b <= a; ++b) S2[(q * p + a) * p + b] = s2[a * p + b]; }
    }
    std::vector<double> ts(n);
    for (std::size_t q = 0; q < n; ++q) ts[q] = time[o[q]];
    for (std::size_t i : ev) {
      const std::size_t f = static_cast<std::size_t>(std::lower_bound(ts.begin(), ts.end(), time[i]) - ts.begin());
      rs.rows.push_back(i); rs.time.push_back(time[i]);
      const double d0 = S0[f];
      std::vector<double> xb(p);
      for (std::size_t a = 0; a < p; ++a) { xb[a] = S1[f * p + a] / d0; rs.xbar.push_back(xb[a]); }
      for (std::size_t a = 0; a < p; ++a)
        for (std::size_t b = 0; b < p; ++b) {
          const std::size_t hi = std::max(a, b), lo = std::min(a, b);
          rs.V.push_back(S2[(f * p + hi) * p + lo] / d0 - xb[a] * xb[b]);
        }
    }
    return rs;
  }
  double mx = -INFINITY;
  for (double v : lp) mx = std::fmax(mx, v);
  std::vector<double> rr(n);
  for (std::size_t i = 0; i < n; ++i) rr[i] = (w ? w[i] : 1.0) * std::exp(lp[i] - mx);
  std::stable_sort(ev.begin(), ev.end(), [&](std::size_t a, std::size_t b) {
    const int sa = strata ? strata[a] : 0, sb = strata ? strata[b] : 0;
    return sa < sb || (sa == sb && time[a] < time[b]);
  });
  for (std::size_t k : ev) {
    const double t = time[k];
    const int sk = strata ? strata[k] : 0;
    double s0 = 0.0;
    std::vector<double> m(p, 0.0), c(p * p, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
      if ((strata ? strata[j] : 0) != sk || !((entry ? entry[j] : -INFINITY) < t) || !(time[j] >= t)) continue;
      s0 += rr[j];
      for (std::size_t a = 0; a < p; ++a) m[a] += X[j * p + a] * rr[j];
    }
    for (std::size_t a = 0; a < p; ++a) m[a] /= s0;
    for (std::size_t j = 0; j < n; ++j) {
      if ((strata ? strata[j] : 0) != sk || !((entry ? entry[j] : -INFINITY) < t) || !(time[j] >= t)) continue;
      const double sr = std::sqrt(rr[j]);
      for (std::size_t a = 0; a < p; ++a)
        for (std::size_t b = 0; b < p; ++b) c[a * p + b] += (X[j * p + a] * sr) * (X[j * p + b] * sr);
    }
    rs.rows.push_back(k); rs.time.push_back(t);
    for (std::size_t a = 0; a < p; ++a) rs.xbar.push_back(m[a]);
    for (std::size_t a = 0; a < p; ++a) for (std::size_t b = 0; b < p; ++b) rs.V.push_back(c[a * p + b] / s0 - m[a] * m[b]);
  }
  return rs;
}

// the transform g of the event times (rows: the event rows, tev their times)
inline std::vector<double> ph_transform(int transform, std::size_t n, const double* time, const int* status, const double* entry,
                                        const double* w, const std::vector<std::size_t>& rows, const std::vector<double>& tev) {
  std::vector<double> g(tev.size());
  if (transform == PHT_IDENTITY) return tev;
  if (transform == PHT_LOG) { for (std::size_t q = 0; q < tev.size(); ++q) g[q] = std::log(tev[q]); return g; }
  if (transform == PHT_RANK) {                         // average ranks of all observed times
    std::vector<std::size_t> o(n);
    std::iota(o.begin(), o.end(), 0);
    std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return time[a] < time[b]; });
    std::vector<double> rk(n);
    for (std::size_t a = 0; a < n;) {
      std::size_t b = a;
      while (b < n && time[o[b]] == time[o[a]]) ++b;
      const double avg = (static_cast<double>(a + 1) + static_cast<double>(b)) / 2.0;
      for (std::size_t c = a; c < b; ++c) rk[o[c]] = avg;
      a = b;
    }
    for (std::size_t q = 0; q < rows.size(); ++q) g[q] = rk[rows[q]];
    return g;
  }
  // 1 - Kaplan-Meier (weighted with w for PHT_KM_WEIGHTED) just before each event time
  std::vector<double> et;
  for (std::size_t i = 0; i < n; ++i) if (status[i] == 1) et.push_back(time[i]);
  std::sort(et.begin(), et.end());
  et.erase(std::unique(et.begin(), et.end()), et.end());
  const bool wt = transform == PHT_KM_WEIGHTED && w;
  std::vector<double> S(et.size());
  double s = 1.0;
  for (std::size_t j = 0; j < et.size(); ++j) {
    double d = 0.0, Y = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double wi = wt ? w[i] : 1.0;
      if (status[i] == 1 && time[i] == et[j]) d += wi;
      if (time[i] >= et[j] && (!entry || entry[i] < et[j])) Y += wi;
    }
    s *= 1.0 - d / Y;
    S[j] = s;
  }
  for (std::size_t q = 0; q < tev.size(); ++q) {
    const std::size_t j = static_cast<std::size_t>(std::lower_bound(et.begin(), et.end(), tev[q]) - et.begin());
    g[q] = 1.0 - (j == 0 ? 1.0 : S[j - 1]);
  }
  return g;
}

struct PhTest {
  std::vector<double> chi_col, chi_term, df_term, g0, scaled, time;
  std::vector<std::size_t> rows;
  double chi_glob = kNaN;
  bool ok = true;
};

// term_of: term index (0..nterm-1) of every column
inline PhTest ph_test(std::size_t n, std::size_t p, const double* time, const int* status, const double* X, const double* beta,
                      const double* vcov, const double* w, const double* entry, const int* strata, int transform,
                      const int* term_of, std::size_t nterm) {
  PhTest out;
  const RiskSets rs = cox_risksets(n, p, time, status, X, beta, w, entry, strata);
  const std::size_t nd = rs.rows.size();
  if (nd == 0) { out.ok = false; return out; }
  out.rows = rs.rows; out.time = rs.time;
  out.g0 = ph_transform(transform, n, time, status, entry, w, rs.rows, rs.time);
  double gm = 0.0;
  for (double v : out.g0) gm += v;
  gm /= static_cast<double>(nd);
  std::vector<double> sch(nd * p), U(p, 0.0), Ibb(p * p, 0.0), Ibg(p * p, 0.0), Igg(p * p, 0.0);
  for (std::size_t q = 0; q < nd; ++q) {
    const double g = out.g0[q] - gm;
    for (std::size_t a = 0; a < p; ++a) {
      sch[q * p + a] = X[rs.rows[q] * p + a] - rs.xbar[q * p + a];
      U[a] += sch[q * p + a] * g;
    }
    for (std::size_t a = 0; a < p * p; ++a) {
      const double v = rs.V[q * p * p + a];
      Ibb[a] += v; Ibg[a] += v * g; Igg[a] += v * g * g;
    }
  }
  // A = Igg - Ibg' Ibb^-1 Ibg
  std::vector<double> Z = Ibg;
  if (!lu_solve(p, Ibb, Z, p)) { out.ok = false; return out; }
  std::vector<double> A(p * p);
  for (std::size_t a = 0; a < p; ++a)
    for (std::size_t b = 0; b < p; ++b) {
      double v = 0.0;
      for (std::size_t c = 0; c < p; ++c) v += Ibg[c * p + a] * Z[c * p + b];
      A[a * p + b] = Igg[a * p + b] - v;
    }
  out.chi_col.resize(p);
  for (std::size_t a = 0; a < p; ++a) out.chi_col[a] = U[a] * U[a] / A[a * p + a];
  auto quad = [&](const std::vector<std::size_t>& J) {
    const std::size_t m = J.size();
    std::vector<double> AJ(m * m), u(m);
    for (std::size_t a = 0; a < m; ++a) { u[a] = U[J[a]]; for (std::size_t b = 0; b < m; ++b) AJ[a * m + b] = A[J[a] * p + J[b]]; }
    std::vector<double> x = u;
    if (!lu_solve(m, AJ, x, 1)) return kNaN;
    double s = 0.0;
    for (std::size_t a = 0; a < m; ++a) s += u[a] * x[a];
    return s;
  };
  out.chi_term.assign(nterm, kNaN); out.df_term.assign(nterm, 0.0);
  for (std::size_t t = 0; t < nterm; ++t) {
    std::vector<std::size_t> J;
    for (std::size_t a = 0; a < p; ++a) if (term_of[a] == static_cast<int>(t)) J.push_back(a);
    out.df_term[t] = static_cast<double>(J.size());
    if (!J.empty()) out.chi_term[t] = quad(J);
  }
  std::vector<std::size_t> all(p);
  std::iota(all.begin(), all.end(), 0);
  out.chi_glob = quad(all);
  out.scaled.assign(nd * p, 0.0);
  for (std::size_t q = 0; q < nd; ++q)
    for (std::size_t a = 0; a < p; ++a) {
      double v = 0.0;
      for (std::size_t c = 0; c < p; ++c) v += sch[q * p + c] * vcov[c * p + a];
      out.scaled[q * p + a] = v * static_cast<double>(nd) + beta[a];
    }
  return out;
}

}  // namespace survcore
