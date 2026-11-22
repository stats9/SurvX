// Data expansions and small matrix summaries that the front ends used to compute themselves:
//   fg_expand  - the weighted counting-process data of the Fine-Gray model with delayed entry
//                (Geskus 2011): subjects with a competing event stay in the risk set with weights
//                G(t) H(t) / {G(T-) H(T-)} from the censoring (G) and truncation (H) distributions;
//   time_fix   - merging of times that differ only by rounding error (as survival's timefix);
//   sandwich   - the robust variance V (sum_c U_c U_c') V, with scores summed within clusters.
#ifndef SURVCORE_EXPAND_HPP
#define SURVCORE_EXPAND_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace survcore {

struct FGRows {
  std::vector<int> id, status;          // id 0-based
  std::vector<double> start, stop, w;
};

inline std::size_t count_lt(const std::vector<double>& s, double a) {
  return static_cast<std::size_t>(std::lower_bound(s.begin(), s.end(), a) - s.begin());
}
inline std::size_t count_le(const std::vector<double>& s, double a) {
  return static_cast<std::size_t>(std::upper_bound(s.begin(), s.end(), a) - s.begin());
}

// time, status (0 censored, k > 0 cause k), entry (may be null: 0); rows with weight 0 or empty
// intervals are left out; the original rows come first, then the extra rows of each subject with a
// competing event, in subject order
inline FGRows fg_expand(std::size_t n, const double* time, const int* status, const double* entry, int cause) {
  std::vector<double> se(n), st(time, time + n);
  for (std::size_t i = 0; i < n; ++i) se[i] = entry ? entry[i] : 0.0;
  std::vector<double> ent = se;
  std::sort(se.begin(), se.end()); std::sort(st.begin(), st.end());
  std::vector<double> ct, lt;
  for (std::size_t i = 0; i < n; ++i) { if (status[i] == 0) ct.push_back(time[i]); if (ent[i] > 0) lt.push_back(ent[i]); }
  std::vector<double> ctall = ct, ltall = lt;
  std::sort(ct.begin(), ct.end()); ct.erase(std::unique(ct.begin(), ct.end()), ct.end());
  std::sort(lt.begin(), lt.end()); lt.erase(std::unique(lt.begin(), lt.end()), lt.end());
  std::sort(ctall.begin(), ctall.end()); std::sort(ltall.begin(), ltall.end());
  // censoring survival G at and after each censoring time: cumG[k] = prod_{j < k} (1 - d_j / Y_j)
  std::vector<double> cumG(ct.size() + 1, 1.0);
  for (std::size_t k = 0; k < ct.size(); ++k) {
    const double Y = static_cast<double>(count_lt(se, ct[k])) - static_cast<double>(count_lt(st, ct[k]));
    const double d = static_cast<double>(count_le(ctall, ct[k]) - count_lt(ctall, ct[k]));
    cumG[k + 1] = cumG[k] * (1.0 - d / Y);
  }
  // truncation: Hrev[k] = prod_{j >= k} (1 - e_j / R_j), Hrev[m] = 1
  std::vector<double> Hrev(lt.size() + 1, 1.0);
  for (std::size_t k = lt.size(); k-- > 0;) {
    const double R = static_cast<double>(count_le(se, lt[k])) - static_cast<double>(count_le(st, lt[k]));
    const double e = static_cast<double>(count_le(ltall, lt[k]) - count_lt(ltall, lt[k]));
    Hrev[k] = (1.0 - e / R) * Hrev[k + 1];
  }
  auto G_right = [&](double a) { return cumG[count_le(ct, a)]; };
  auto G_left = [&](double a) { return cumG[count_lt(ct, a)]; };
  auto H_after = [&](double a) { return Hrev[count_le(lt, a)]; };
  auto H_left = [&](double a) { return Hrev[count_lt(lt, a)]; };
  std::vector<double> jumps = ct;
  jumps.insert(jumps.end(), lt.begin(), lt.end());
  std::sort(jumps.begin(), jumps.end()); jumps.erase(std::unique(jumps.begin(), jumps.end()), jumps.end());
  double tmax = 0.0;
  for (std::size_t i = 0; i < n; ++i) if (status[i] == cause) tmax = std::max(tmax, time[i]);
  FGRows out;
  auto push = [&](int id, double a, double b, int s, double w) {
    if (!(w > 0) || !(b > a)) return;
    out.id.push_back(id); out.start.push_back(a); out.stop.push_back(b); out.status.push_back(s); out.w.push_back(w);
  };
  for (std::size_t i = 0; i < n; ++i) push(static_cast<int>(i), ent[i], time[i], status[i] == cause, 1.0);
  for (std::size_t i = 0; i < n; ++i) {
    if (status[i] == 0 || status[i] == cause || !(time[i] < tmax)) continue;
    const double Ti = time[i];
    std::vector<double> cuts(1, Ti);
    for (std::size_t j = count_le(jumps, Ti); j < jumps.size() && jumps[j] < tmax; ++j) cuts.push_back(jumps[j]);
    const double den = G_left(Ti) * H_left(Ti);
    for (std::size_t c = 0; c < cuts.size(); ++c) {
      const double b = c + 1 < cuts.size() ? cuts[c + 1] : tmax;
      const double w = den > 0 ? G_right(cuts[c]) * H_after(cuts[c]) / den : 0.0;
      push(static_cast<int>(i), cuts[c], b, 0, w);
    }
  }
  return out;
}

// R's mean(): long-double sum, then one correction pass
inline double r_mean(const std::vector<double>& x) {
  long double s = 0;
  for (double v : x) s += v;
  s /= static_cast<long double>(x.size());
  long double t = 0;
  for (double v : x) t += v - s;
  return static_cast<double>(s + t / static_cast<long double>(x.size()));
}

// v: all time values (non-finite ones are left alone); returns true when some values were merged
inline bool time_fix(std::size_t n, double* v, double tol) {
  std::vector<double> u;
  for (std::size_t i = 0; i < n; ++i) if (std::isfinite(v[i])) u.push_back(v[i]);
  std::sort(u.begin(), u.end()); u.erase(std::unique(u.begin(), u.end()), u.end());
  if (u.size() < 2) return false;
  std::vector<double> au(u.size());
  for (std::size_t k = 0; k < u.size(); ++k) au[k] = std::fabs(u[k]);
  const double ma = r_mean(au);
  std::vector<double> cuts(1, u[0]);
  bool any = false;
  for (std::size_t k = 0; k + 1 < u.size(); ++k) {
    const double du = u[k + 1] - u[k];
    const bool tied = du <= tol || du / ma <= tol;
    any = any || tied;
    if (!tied) cuts.push_back(u[k + 1]);
  }
  if (!any) return false;
  for (std::size_t i = 0; i < n; ++i)
    if (std::isfinite(v[i])) v[i] = cuts[count_le(cuts, v[i]) - 1];
  return true;
}

// out = V (sum_c U_c U_c') V; U n x p row-major, cl 0-based cluster codes 0..nc-1 (null: each row)
inline void sandwich(std::size_t n, std::size_t p, const double* U, const int* cl, std::size_t nc, const double* V, double* out) {
  const std::size_t m = cl ? nc : n;
  std::vector<double> S(m * p, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t c = cl ? static_cast<std::size_t>(cl[i]) : i;
    for (std::size_t j = 0; j < p; ++j) S[c * p + j] += U[i * p + j];
  }
  std::vector<double> M(p * p, 0.0), VM(p * p, 0.0);
  for (std::size_t j = 0; j < p; ++j)
    for (std::size_t k = 0; k < p; ++k) {
      double s = 0; for (std::size_t c = 0; c < m; ++c) s += S[c * p + j] * S[c * p + k];
      M[j * p + k] = s;
    }
  for (std::size_t i = 0; i < p; ++i)
    for (std::size_t k = 0; k < p; ++k) {
      double s = 0; for (std::size_t j = 0; j < p; ++j) s += V[i * p + j] * M[j * p + k];
      VM[i * p + k] = s;
    }
  for (std::size_t i = 0; i < p; ++i)
    for (std::size_t k = 0; k < p; ++k) {
      double s = 0; for (std::size_t j = 0; j < p; ++j) s += VM[i * p + j] * V[j * p + k];
      out[i * p + k] = s;
    }
}


// sample covariance (divisor n - 1) of the rows of X (nb x p, row-major) without a NaN, as
// stats::cov(X[stats::complete.cases(X), ]): means and sums of products accumulated in long double
// with the correction pass of R; *nok receives the number of rows used (NaN output when below 2)
inline void cov_complete(std::size_t nb, std::size_t p, const double* X, double* out, std::size_t& nok) {
  std::vector<std::size_t> rows;
  for (std::size_t b = 0; b < nb; ++b) {
    bool ok = true;
    for (std::size_t j = 0; j < p; ++j) if (std::isnan(X[b * p + j])) { ok = false; break; }
    if (ok) rows.push_back(b);
  }
  nok = rows.size();
  if (nok < 2) { for (std::size_t k = 0; k < p * p; ++k) out[k] = std::nan(""); return; }
  const long double nn = static_cast<long double>(nok);
  std::vector<long double> mu(p);
  for (std::size_t j = 0; j < p; ++j) {
    long double s = 0; for (std::size_t b : rows) s += X[b * p + j];
    s /= nn;
    long double t = 0; for (std::size_t b : rows) t += X[b * p + j] - s;
    mu[j] = s + t / nn;
  }
  for (std::size_t j = 0; j < p; ++j)
    for (std::size_t k = 0; k <= j; ++k) {
      long double s = 0;
      for (std::size_t b : rows) s += (X[b * p + j] - mu[j]) * (X[b * p + k] - mu[k]);
      out[j * p + k] = out[k * p + j] = static_cast<double>(s / (nn - 1));
    }
}

}  // namespace survcore

#endif
