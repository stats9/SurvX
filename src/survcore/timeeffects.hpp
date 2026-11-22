// survcore: the numerical parts of time-varying effects (time_effects() in R): splitting of
// follow-up at event times, the built-in time functions (log, linear, steps, natural cubic spline
// in log time, as splines::ns), the effect beta + gamma' B(t) with its delta-method interval, and
// survival and cumulative incidence from hazard increments.
#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

#include "linalg.hpp"
#include "summaries.hpp"

namespace survcore {

// ---------------------------------------------------------------------------------------------
// splitting (start, stop] at the event times et (sorted) strictly inside the interval
// ---------------------------------------------------------------------------------------------
inline std::size_t split_lo(const std::vector<double>& et, double s) {          // #{et <= s}
  return static_cast<std::size_t>(std::upper_bound(et.begin(), et.end(), s) - et.begin());
}
inline std::size_t split_hi(const std::vector<double>& et, double s) {          // #{et < s}
  return static_cast<std::size_t>(std::lower_bound(et.begin(), et.end(), s) - et.begin());
}
inline double split_size(std::size_t n, const double* start, const double* stop, const std::vector<double>& et) {
  double tot = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double d = static_cast<double>(split_hi(et, stop[i])) - static_cast<double>(split_lo(et, start[i]));
    tot += std::fmax(d, 0.0) + 1.0;
  }
  return tot;
}
// rows of the split data: source row (0-based), start, stop, status
inline void split_rows(std::size_t n, const double* start, const double* stop, const int* status, const std::vector<double>& et,
                       std::vector<int>& idx, std::vector<double>& s0, std::vector<double>& s1, std::vector<int>& st) {
  idx.clear(); s0.clear(); s1.clear(); st.clear();
  const std::size_t m = et.size();
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t lo = split_lo(et, start[i]), hi = split_hi(et, stop[i]);
    const std::size_t rows = (hi > lo ? hi - lo : 0) + 1;
    for (std::size_t r = 1; r <= rows; ++r) {
      const std::size_t ci = lo + r;                                  // 1-based index of the cut after this piece
      const bool last = r == rows;
      idx.push_back(static_cast<int>(i));
      s0.push_back(r == 1 ? start[i] : et[std::max<std::size_t>(ci - 1, 1) - 1]);
      s1.push_back(last ? stop[i] : et[std::min(ci, m) - 1]);
      st.push_back(last ? status[i] : 0);
    }
  }
}

// ---------------------------------------------------------------------------------------------
// B-splines (Cox-de Boor) and natural cubic splines as splines::ns
// ---------------------------------------------------------------------------------------------
// value (deriv = 0) or derivative of every B-spline of order k on the knots t at x; x equal to the
// right boundary knot is taken in the last non-empty interval (left continuity, as splineDesign)
inline std::vector<double> bspline_row(const std::vector<double>& t, int k, double x, int deriv) {
  const std::size_t nk = t.size(), nb = nk - static_cast<std::size_t>(k);
  // interval mu: t[mu] <= x < t[mu + 1], or the last interval with t[mu] < t[mu + 1] at the boundary
  std::size_t mu = static_cast<std::size_t>(k - 1);
  const double right = t[nb];
  if (x >= right) {
    mu = nb - 1;
    while (mu > 0 && !(t[mu] < t[mu + 1])) --mu;
  } else {
    for (std::size_t j = static_cast<std::size_t>(k - 1); j < nb; ++j) if (t[j] <= x && x < t[j + 1]) { mu = j; break; }
  }
  // order-1 basis: indicator of interval mu; raise the order, differentiating at the last `deriv` steps
  std::vector<double> B(nk - 1, 0.0);
  B[mu] = 1.0;
  for (int ord = 2; ord <= k; ++ord) {
    std::vector<double> N(nk - static_cast<std::size_t>(ord), 0.0);
    const bool diff = ord > k - deriv;
    for (std::size_t i = 0; i < N.size(); ++i) {
      const double d1 = t[i + ord - 1] - t[i], d2 = t[i + ord] - t[i + 1];
      double a = 0.0, b = 0.0;
      if (diff) {
        if (d1 > 0.0) a = (ord - 1) * B[i] / d1;
        if (d2 > 0.0) b = (ord - 1) * B[i + 1] / d2;
        N[i] = a - b;
      } else {
        if (d1 > 0.0) a = (x - t[i]) / d1 * B[i];
        if (d2 > 0.0) b = (t[i + ord] - x) / d2 * B[i + 1];
        N[i] = a + b;
      }
    }
    B.swap(N);
  }
  return B;                                         // nk - k values
}

// natural cubic spline basis without intercept at x (as splines::ns(x, knots, Boundary.knots)):
// B-splines on the augmented knots, the first column dropped, projected on the null space of the
// second derivatives at the boundary knots by the Householder QR of LINPACK (dqrdc2); linear
// beyond the boundary knots
struct NaturalSpline {
  std::vector<double> aknots;
  std::vector<double> u1, u2;      // Householder vectors (full length nb) of the two constraint columns
  std::size_t nb = 0;              // B-splines after dropping the first
};

inline NaturalSpline ns_setup(const std::vector<double>& knots, double bl, double br) {
  NaturalSpline s;
  s.aknots = {bl, bl, bl, bl};
  for (double v : knots) s.aknots.push_back(v);
  s.aknots.insert(s.aknots.end(), {br, br, br, br});
  std::sort(s.aknots.begin(), s.aknots.end());
  const std::size_t nall = s.aknots.size() - 4;
  s.nb = nall - 1;
  // constraint: second derivatives at the boundary knots (2 x nall), first column dropped; QR of its transpose (nb x 2)
  std::vector<double> c0 = bspline_row(s.aknots, 4, bl, 2), c1 = bspline_row(s.aknots, 4, br, 2);
  std::vector<std::vector<double>> A = {std::vector<double>(c0.begin() + 1, c0.end()), std::vector<double>(c1.begin() + 1, c1.end())};
  std::vector<double>* U[2] = {&s.u1, &s.u2};
  for (std::size_t l = 0; l < 2; ++l) {
    std::vector<double>& x = A[l];
    double nrm = 0.0;
    for (std::size_t i = l; i < s.nb; ++i) nrm += x[i] * x[i];
    nrm = std::sqrt(nrm);
    std::vector<double> u(s.nb, 0.0);
    if (nrm != 0.0) {
      if (x[l] != 0.0) nrm = std::copysign(nrm, x[l]);
      const double sc = 1.0 / nrm;                       // dscal by 1 / nrmxl, as LINPACK
      for (std::size_t i = l; i < s.nb; ++i) u[i] = x[i] * sc;
      u[l] += 1.0;
      // apply to the remaining column
      for (std::size_t j = l + 1; j < 2; ++j) {
        double tt = 0.0;
        for (std::size_t i = l; i < s.nb; ++i) tt -= u[i] * A[j][i];
        tt /= u[l];
        for (std::size_t i = l; i < s.nb; ++i) A[j][i] += tt * u[i];
      }
    }
    *U[l] = u;
  }
  return s;
}

// Q' y for the two Householder transformations (LINPACK dqrsl), then rows 3.. (the projected basis)
inline void ns_project(const NaturalSpline& s, std::vector<double>& y) {
  const std::vector<double>* U[2] = {&s.u1, &s.u2};
  for (std::size_t l = 0; l < 2; ++l) {
    const std::vector<double>& u = *U[l];
    if (u[l] == 0.0) continue;
    double tt = 0.0;
    for (std::size_t i = l; i < s.nb; ++i) tt -= u[i] * y[i];
    tt /= u[l];
    for (std::size_t i = l; i < s.nb; ++i) y[i] += tt * u[i];
  }
}

inline std::vector<double> ns_row(const NaturalSpline& s, double x) {
  const double bl = s.aknots.front(), br = s.aknots.back();
  std::vector<double> b(s.nb + 1, 0.0);
  if (x < bl || x > br) {
    const double piv = x < bl ? bl : br;
    const std::vector<double> v0 = bspline_row(s.aknots, 4, piv, 0), v1 = bspline_row(s.aknots, 4, piv, 1);
    for (std::size_t j = 0; j < b.size(); ++j) b[j] = v0[j] + (x - piv) * v1[j];
  } else b = bspline_row(s.aknots, 4, x, 0);
  std::vector<double> y(b.begin() + 1, b.end());
  ns_project(s, y);
  return std::vector<double>(y.begin() + 2, y.end());
}

// ---------------------------------------------------------------------------------------------
// built-in time functions B(t): 1 log(max(t, tmin / 2) / tm), 2 t - tm, 3 steps 1{t > s_j},
// 4 natural spline in log(max(t, tmin / 2)) with knots and boundary knots on the log scale
// ---------------------------------------------------------------------------------------------
enum TeForm : int { TEF_LOG = 1, TEF_LINEAR = 2, TEF_STEP = 3, TEF_SPLINE = 4 };

inline std::size_t te_basis_cols(int form, std::size_t nbreak, std::size_t nknot) {
  return form == TEF_STEP ? nbreak : form == TEF_SPLINE ? nknot + 1 : 1;
}

inline void te_basis(int form, std::size_t nt, const double* t, double tmin, double tm, const double* breaks, std::size_t nbreak,
                     const double* knots, std::size_t nknot, double bl, double br, double* out) {
  const std::size_t q = te_basis_cols(form, nbreak, nknot);
  NaturalSpline ns;
  if (form == TEF_SPLINE) ns = ns_setup(std::vector<double>(knots, knots + nknot), bl, br);
  for (std::size_t i = 0; i < nt; ++i) {
    const double tp = std::fmax(t[i], tmin / 2.0);
    double* o = out + i * q;
    switch (form) {
      case TEF_LOG: o[0] = std::log(tp / tm); break;
      case TEF_LINEAR: o[0] = t[i] - tm; break;
      case TEF_STEP: for (std::size_t j = 0; j < nbreak; ++j) o[j] = t[i] > breaks[j] ? 1.0 : 0.0; break;
      default: { const std::vector<double> r = ns_row(ns, std::log(tp)); for (std::size_t j = 0; j < q; ++j) o[j] = r[j]; }
    }
  }
}

// ---------------------------------------------------------------------------------------------
// effect of one column at times: log HR = b_j + gamma' B(t) (B: nt x q, may be q = 0), delta-method
// standard error from V (P x P) over the indices ix = (j, tv...)
// ---------------------------------------------------------------------------------------------
inline void te_effect(std::size_t nt, std::size_t q, const double* B, std::size_t P, const double* b, const double* V,
                      const int* ix, double zq, double* est, double* se, double* hr, double* lower, double* upper, double* p) {
  // ix[0] < 0: no constant coefficient (the effect is the time function part alone)
  const std::size_t c0 = ix[0] < 0 ? 1 : 0, m = q + 1;
  std::vector<double> a(m);
  for (std::size_t i = 0; i < nt; ++i) {
    a[0] = 1.0;
    for (std::size_t c = 0; c < q; ++c) a[c + 1] = B[i * q + c];
    double e = 0.0;
    for (std::size_t c = c0; c < m; ++c) e += a[c] * b[ix[c]];
    double v = 0.0;
    for (std::size_t r = c0; r < m; ++r) {
      double av = 0.0;
      for (std::size_t c = c0; c < m; ++c) av += a[c] * V[static_cast<std::size_t>(ix[c]) * P + static_cast<std::size_t>(ix[r])];
      v += av * a[r];
    }
    const double s = std::sqrt(std::fmax(v, 0.0));
    est[i] = e; se[i] = s; hr[i] = std::exp(e);
    lower[i] = std::exp(e - zq * s); upper[i] = std::exp(e + zq * s);
    p[i] = 2.0 * pnorm_std(-std::fabs(e / s));
  }
}

// ---------------------------------------------------------------------------------------------
// survival or cumulative incidence from hazard increments dH (K causes, each nr x m on the event
// times tj), evaluated at times: mode 0 survival exp(-H) of one hazard, 1 cumulative incidence
// 1 - exp(-H) of one subdistribution hazard, 2 cumulative incidence of cause `cause` (0-based)
// from K cause-specific hazards with S(t-) by product limit (product_limit) or exp(-sum H)
// ---------------------------------------------------------------------------------------------
inline void cum_predict(int mode, std::size_t nr, std::size_t m, const double* tj, std::size_t K, const double* dH, int cause,
                        bool product_limit, std::size_t nt, const double* times, double* out) {
  std::vector<std::size_t> kk(nt);
  for (std::size_t q = 0; q < nt; ++q) kk[q] = static_cast<std::size_t>(std::upper_bound(tj, tj + m, times[q]) - tj);
  std::vector<double> path(m);
  for (std::size_t i = 0; i < nr; ++i) {
    if (mode == 0 || mode == 1) {
      double H = 0.0;
      for (std::size_t j = 0; j < m; ++j) { H += dH[i * m + j]; path[j] = H; }
      for (std::size_t q = 0; q < nt; ++q) {
        const double h = kk[q] > 0 ? path[kk[q] - 1] : 0.0;
        out[i * nt + q] = mode == 0 ? std::exp(-h) : 1.0 - std::exp(-h);
      }
      if (mode == 0) for (std::size_t q = 0; q < nt; ++q) if (kk[q] == 0) out[i * nt + q] = 1.0;
      continue;
    }
    double S = 1.0, Ht = 0.0, F = 0.0;
    for (std::size_t j = 0; j < m; ++j) {
      double tot = 0.0;
      for (std::size_t c = 0; c < K; ++c) tot += dH[(c * nr + i) * m + j];
      const double Sm = S;                                   // S just before tj
      F += Sm * dH[(static_cast<std::size_t>(cause) * nr + i) * m + j];
      path[j] = F;
      if (product_limit) S *= std::fmax(1.0 - tot, 0.0);
      else { Ht += tot; S = std::exp(-Ht); }
    }
    for (std::size_t q = 0; q < nt; ++q) out[i * nt + q] = kk[q] > 0 ? path[kk[q] - 1] : 0.0;
  }
}


// ---------------------------------------------------------------------------------------------
// smoother of the residual plots: least squares fit of y on an intercept and ns(x, df) (knots at
// the type-7 quantiles j / df of x, boundary knots the range of x), as stats::lm(y ~ splines::ns(x, df)),
// with the fitted line and pointwise limits fit +- z se at the grid. Householder QR; false when the
// design is rank deficient or there are fewer than df + 1 points.
// ---------------------------------------------------------------------------------------------
inline bool ns_smooth(std::size_t n, const double* x, const double* y, int df, std::size_t ng, const double* grid, double z,
                      double* fit, double* lower, double* upper) {
  const std::size_t q = static_cast<std::size_t>(df) + 1;
  if (df < 1 || n < q) return false;
  std::vector<double> xs(x, x + n);
  std::sort(xs.begin(), xs.end());
  std::vector<double> kn;
  for (int j = 1; j < df; ++j) kn.push_back(quantile7_sorted(xs, static_cast<double>(j) / df));
  const double bl = xs.front(), br = xs.back();
  if (!(bl < br)) return false;
  const NaturalSpline sp = ns_setup(kn, bl, br);
  auto row = [&](double v) { std::vector<double> r(1, 1.0); const std::vector<double> b = ns_row(sp, v); r.insert(r.end(), b.begin(), b.end()); return r; };
  std::vector<double> A(n * q), b(y, y + n);                       // A row-major n x q
  for (std::size_t i = 0; i < n; ++i) { const std::vector<double> r = row(x[i]); for (std::size_t j = 0; j < q; ++j) A[i * q + j] = r[j]; }
  double scale = 0.0;
  for (double v : A) scale = std::fmax(scale, std::fabs(v));
  for (std::size_t j = 0; j < q; ++j) {                            // Householder on column j
    double nrm = 0.0;
    for (std::size_t i = j; i < n; ++i) nrm += A[i * q + j] * A[i * q + j];
    nrm = std::sqrt(nrm);
    if (!(nrm > 1e-10 * scale * std::sqrt(static_cast<double>(n)))) return false;
    const double alpha = A[j * q + j] > 0 ? -nrm : nrm;
    std::vector<double> u(n - j);
    for (std::size_t i = j; i < n; ++i) u[i - j] = A[i * q + j];
    u[0] -= alpha;
    double uu = 0.0; for (double v : u) uu += v * v;
    if (uu > 0) {
      for (std::size_t c = j; c < q; ++c) {
        double t = 0.0; for (std::size_t i = j; i < n; ++i) t += u[i - j] * A[i * q + c];
        t = 2.0 * t / uu; for (std::size_t i = j; i < n; ++i) A[i * q + c] -= t * u[i - j];
      }
      double t = 0.0; for (std::size_t i = j; i < n; ++i) t += u[i - j] * b[i];
      t = 2.0 * t / uu; for (std::size_t i = j; i < n; ++i) b[i] -= t * u[i - j];
    }
  }
  std::vector<double> beta(q);                                     // back substitution R beta = (Q'y)[0..q)
  for (std::size_t j = q; j-- > 0;) {
    double s = b[j]; for (std::size_t c = j + 1; c < q; ++c) s -= A[j * q + c] * beta[c];
    beta[j] = s / A[j * q + j];
  }
  double rss = 0.0; for (std::size_t i = q; i < n; ++i) rss += b[i] * b[i];
  const double sigma = n > q ? std::sqrt(rss / static_cast<double>(n - q)) : std::nan("");
  for (std::size_t g = 0; g < ng; ++g) {
    const std::vector<double> r = row(grid[g]);
    double f = 0.0; for (std::size_t j = 0; j < q; ++j) f += r[j] * beta[j];
    std::vector<double> v(q);                                      // R' v = r, se = sigma |v|
    for (std::size_t j = 0; j < q; ++j) { double s = r[j]; for (std::size_t c = 0; c < j; ++c) s -= A[c * q + j] * v[c]; v[j] = s / A[j * q + j]; }
    double vv = 0.0; for (double e : v) vv += e * e;
    const double se = sigma * std::sqrt(vv);
    fit[g] = f; lower[g] = f - z * se; upper[g] = f + z * se;
  }
  return true;
}

}  // namespace survcore
