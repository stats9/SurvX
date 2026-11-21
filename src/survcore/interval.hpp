// survcore: non-parametric maximum likelihood estimation (NPMLE) of a survival function from
// interval-censored data.
//
// Data: subject i is known to have the event in (L_i, R_i]; L_i = R_i is an exact time, R_i = Inf
// right censoring at L_i, L_i = 0 left censoring at R_i.
//
// Support (Turnbull 1976, JRSS B 38:290-295; Peto 1973, Appl Statist 22:86-91): the NPMLE puts all
// its mass on the innermost intervals, the intervals (q, p] whose left end q is a left end point and
// whose right end p is a right end point of the data with no other end point in between. For the
// ordering of tied end points, an exact time t is the closed set [t, t]: its left end sorts just
// before t, a closed right end at v before an open left end at v (an interval (v, .] does not contain
// v). Every subject covers a contiguous block a_i..b_i of innermost intervals.
//
// Estimation: the self-consistency (EM) algorithm of Turnbull (1976),
//   p_j <- p_j / W sum_i w_i 1{a_i <= j <= b_i} / P_i,   P_i = sum_{j = a_i}^{b_i} p_j,
// in O(n + m) per step with prefix sums, accelerated by SQUAREM (Varadhan & Roland 2008, Scand J
// Statist 35:335-353, scheme S3) with a monotonicity safeguard. Convergence is checked on the
// Kuhn-Tucker conditions of the concave problem (Gentleman & Geyer 1994, Biometrika 81:618-623):
// d_j = sum_i w_i 1{a_i <= j <= b_i} / P_i <= W for every j, with equality where p_j > 0; the reported
// value max_j d_j / W - 1 is zero at the NPMLE.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace survcore {

struct ICKey {
  double v; int o;   // o: 0 left end of an exact time, 1 closed right end, 2 open left end
  bool operator<(const ICKey& b) const { return v < b.v || (v == b.v && o < b.o); }
  bool operator<=(const ICKey& b) const { return !(b < *this); }
};

struct TurnbullSupport {
  std::vector<double> lo, hi;     // innermost intervals (lo, hi]; lo = hi: a point mass at an exact time
  std::vector<int> a, b;          // block of innermost intervals covered by subject i (0-based)
  std::vector<ICKey> klo, khi;
};

inline ICKey ic_left_key(double L, double R) { return L == R ? ICKey{L, 0} : ICKey{L, 2}; }
inline ICKey ic_right_key(double R) { return ICKey{R, 1}; }

inline TurnbullSupport turnbull_support(int n, const double* L, const double* R) {
  struct E { ICKey k; int left; };
  std::vector<E> e; e.reserve(2 * static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) { e.push_back({ic_left_key(L[i], R[i]), 1}); e.push_back({ic_right_key(R[i]), 0}); }
  std::sort(e.begin(), e.end(), [](const E& x, const E& y) {
    if (x.k < y.k) return true;
    if (y.k < x.k) return false;
    return x.left > y.left;
  });
  TurnbullSupport s;
  for (std::size_t m = 0; m + 1 < e.size(); ++m) {
    if (e[m].left && !e[m + 1].left) {
      s.klo.push_back(e[m].k); s.khi.push_back(e[m + 1].k);
      s.lo.push_back(e[m].k.v); s.hi.push_back(e[m + 1].k.v);
    }
  }
  const int J = static_cast<int>(s.lo.size());
  s.a.resize(n); s.b.resize(n);
  for (int i = 0; i < n; ++i) {
    const ICKey kl = ic_left_key(L[i], R[i]), kr = ic_right_key(R[i]);
    // first innermost interval whose left end is at or after the subject's left end
    int lo = 0, hi = J;
    while (lo < hi) { const int mid = (lo + hi) / 2; if (s.klo[mid] < kl) lo = mid + 1; else hi = mid; }
    s.a[i] = lo;
    // last innermost interval whose right end is at or before the subject's right end
    lo = 0; hi = J;
    while (lo < hi) { const int mid = (lo + hi) / 2; if (s.khi[mid] <= kr) lo = mid + 1; else hi = mid; }
    s.b[i] = lo - 1;
  }
  return s;
}

struct TurnbullResult {
  std::vector<double> p;
  double loglik = -INFINITY, kkt = INFINITY;
  int iterations = 0;
  bool converged = false;
};

namespace ic_detail {

inline double ic_loglik(const std::vector<double>& p, const std::vector<int>& a, const std::vector<int>& b,
                        const double* w, std::vector<double>& cum) {
  const std::size_t J = p.size();
  cum.assign(J + 1, 0.0);
  for (std::size_t j = 0; j < J; ++j) cum[j + 1] = cum[j] + p[j];
  double ll = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (w && w[i] == 0.0) continue;
    const double P = cum[b[i] + 1] - cum[a[i]];
    if (!(P > 0.0)) return -INFINITY;
    ll += (w ? w[i] : 1.0) * std::log(P);
  }
  return ll;
}

// d_j = sum_i w_i 1{a_i <= j <= b_i} / P_i
inline void ic_gradient(const std::vector<double>& p, const std::vector<int>& a, const std::vector<int>& b,
                        const double* w, std::vector<double>& cum, std::vector<double>& d) {
  const std::size_t J = p.size();
  cum.assign(J + 1, 0.0);
  for (std::size_t j = 0; j < J; ++j) cum[j + 1] = cum[j] + p[j];
  std::vector<double> diff(J + 1, 0.0);
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (w && w[i] == 0.0) continue;
    const double P = cum[b[i] + 1] - cum[a[i]];
    const double v = (w ? w[i] : 1.0) / std::fmax(P, 1e-300);
    diff[a[i]] += v; diff[b[i] + 1] -= v;
  }
  d.assign(J, 0.0);
  double run = 0.0;
  for (std::size_t j = 0; j < J; ++j) { run += diff[j]; d[j] = run; }
}

inline void ic_em(std::vector<double>& p, const std::vector<int>& a, const std::vector<int>& b, const double* w, double W,
                  std::vector<double>& cum, std::vector<double>& d) {
  ic_gradient(p, a, b, w, cum, d);
  double s = 0.0;
  for (std::size_t j = 0; j < p.size(); ++j) { p[j] *= d[j] / W; s += p[j]; }
  if (s > 0.0) for (double& v : p) v /= s;
}

}  // namespace ic_detail

// NPMLE on the innermost intervals of `s`; w (may be null) are case weights (bootstrap counts)
inline TurnbullResult turnbull_fit(const TurnbullSupport& s, const double* w, int max_iter, double tol) {
  using namespace ic_detail;
  const std::size_t J = s.lo.size(), n = s.a.size();
  TurnbullResult r;
  double W = 0.0;
  for (std::size_t i = 0; i < n; ++i) W += w ? w[i] : 1.0;
  if (!J || !(W > 0.0)) return r;
  // start: equal mass on the innermost intervals covered by at least one (weighted) subject
  std::vector<double> p(J, 0.0), cum, d;
  {
    std::vector<int> cov(J + 1, 0);
    for (std::size_t i = 0; i < n; ++i) if (!w || w[i] > 0.0) { cov[s.a[i]] += 1; cov[s.b[i] + 1] -= 1; }
    int run = 0; double tot = 0.0;
    for (std::size_t j = 0; j < J; ++j) { run += cov[j]; p[j] = run > 0 ? 1.0 : 0.0; tot += p[j]; }
    for (double& v : p) v /= tot;
  }
  double ll = ic_loglik(p, s.a, s.b, w, cum);
  std::vector<double> p1, p2, q(J), rr(J), vv(J);
  for (int it = 1; it <= max_iter; ++it) {
    r.iterations = it;
    p1 = p; ic_em(p1, s.a, s.b, w, W, cum, d);
    p2 = p1; ic_em(p2, s.a, s.b, w, W, cum, d);
    double nr = 0.0, nv = 0.0;
    for (std::size_t j = 0; j < J; ++j) { rr[j] = p1[j] - p[j]; vv[j] = p2[j] - p1[j] - rr[j]; nr += rr[j] * rr[j]; nv += vv[j] * vv[j]; }
    std::vector<double> next = p2;
    if (nv > 0.0 && nr > 0.0) {
      double alpha = -std::sqrt(nr / nv);
      if (alpha > -1.0) alpha = -1.0;
      double tot = 0.0;
      for (std::size_t j = 0; j < J; ++j) { q[j] = std::fmax(p[j] - 2.0 * alpha * rr[j] + alpha * alpha * vv[j], 0.0); tot += q[j]; }
      if (tot > 0.0) {
        for (double& v : q) v /= tot;
        ic_em(q, s.a, s.b, w, W, cum, d);             // stabilising EM step
        const double lq = ic_loglik(q, s.a, s.b, w, cum), l2 = ic_loglik(p2, s.a, s.b, w, cum);
        if (std::isfinite(lq) && lq >= l2) next = q;
      }
    }
    const double lnew = ic_loglik(next, s.a, s.b, w, cum);
    p.swap(next);
    const double dll = lnew - ll;
    ll = lnew;
    if (it % 5 == 0 || dll < tol) {
      ic_gradient(p, s.a, s.b, w, cum, d);
      double mx = 0.0;
      for (std::size_t j = 0; j < J; ++j) mx = std::fmax(mx, d[j] / W);
      r.kkt = mx - 1.0;
      if (r.kkt < 1e-7 && std::fabs(dll) < tol) { r.converged = true; break; }
    }
  }
  // masses that the algorithm drives to zero only asymptotically
  {
    std::vector<double> pz = p;
    double tot = 0.0;
    for (double& v : pz) { if (v < 1e-12) v = 0.0; tot += v; }
    for (double& v : pz) v /= tot;
    if (std::isfinite(ic_loglik(pz, s.a, s.b, w, cum))) p.swap(pz);
  }
  r.loglik = ic_loglik(p, s.a, s.b, w, cum);
  ic_gradient(p, s.a, s.b, w, cum, d);
  double mx = 0.0;
  for (std::size_t j = 0; j < J; ++j) mx = std::fmax(mx, d[j] / W);
  r.kkt = mx - 1.0;
  r.p = p;
  return r;
}

}  // namespace survcore
