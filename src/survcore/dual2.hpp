// survcore: second-order forward-mode automatic differentiation in five variables.
// Used for the per-cause "local" variables (up to three baseline quantities, alpha, eta) of the direct
// cumulative incidence models, so that gradients and Hessians are exact to rounding
// error (no finite differences). Forward-mode AD with truncated Taylor arithmetic:
//   Griewank & Walther (2008), Evaluating Derivatives, 2nd ed., Chapter 13.
#pragma once
#include <cmath>

namespace survcore {

struct D2 {
  static constexpr int N = 5;
  double v = 0.0;
  double g[N] = {0, 0, 0, 0, 0};
  double h[N][N] = {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
  D2() = default;
  explicit D2(double c) : v(c) {}
  static D2 variable(int i, double x) { D2 r(x); r.g[i] = 1.0; return r; }
};

inline D2 operator+(const D2& a, const D2& b) {
  D2 r(a.v + b.v);
  for (int i = 0; i < D2::N; ++i) { r.g[i] = a.g[i] + b.g[i]; for (int j = 0; j < D2::N; ++j) r.h[i][j] = a.h[i][j] + b.h[i][j]; }
  return r;
}
inline D2 operator-(const D2& a, const D2& b) {
  D2 r(a.v - b.v);
  for (int i = 0; i < D2::N; ++i) { r.g[i] = a.g[i] - b.g[i]; for (int j = 0; j < D2::N; ++j) r.h[i][j] = a.h[i][j] - b.h[i][j]; }
  return r;
}
inline D2 operator*(const D2& a, const D2& b) {
  D2 r(a.v * b.v);
  for (int i = 0; i < D2::N; ++i) {
    r.g[i] = a.g[i] * b.v + a.v * b.g[i];
    for (int j = 0; j < D2::N; ++j) r.h[i][j] = a.h[i][j] * b.v + a.g[i] * b.g[j] + a.g[j] * b.g[i] + a.v * b.h[i][j];
  }
  return r;
}
inline D2 operator*(double c, const D2& a) {
  D2 r(c * a.v);
  for (int i = 0; i < D2::N; ++i) { r.g[i] = c * a.g[i]; for (int j = 0; j < D2::N; ++j) r.h[i][j] = c * a.h[i][j]; }
  return r;
}
inline D2 operator+(const D2& a, double c) { D2 r = a; r.v += c; return r; }
// Chain rule for a scalar function with value f, first derivative f1 and second derivative f2 at a.v.
inline D2 chain(const D2& a, double f, double f1, double f2) {
  D2 r(f);
  for (int i = 0; i < D2::N; ++i) { r.g[i] = f1 * a.g[i]; for (int j = 0; j < D2::N; ++j) r.h[i][j] = f1 * a.h[i][j] + f2 * a.g[i] * a.g[j]; }
  return r;
}
inline D2 exp(const D2& a) { const double e = std::exp(a.v); return chain(a, e, e, e); }

// phi(w) = (e^w - 1)/w with phi(0) = 1, and its first two derivatives. Power series
// for |w| < 1 (no cancellation), closed forms otherwise.
inline void phi_fun(double w, double& f, double& f1, double& f2) {
  if (std::fabs(w) < 1.0) {
    // phi = sum_{n>=0} w^n/(n+1)!, phi' = sum_{n>=1} n w^(n-1)/(n+1)!, phi'' = sum_{n>=2} n(n-1) w^(n-2)/(n+1)!
    double pw[30];
    pw[0] = 1.0;
    for (int n = 1; n < 30; ++n) pw[n] = pw[n - 1] * w;
    f = 0.0; f1 = 0.0; f2 = 0.0;
    double c = 1.0;                                    // 1/(n+1)!
    for (int n = 0; n < 30; ++n) {
      if (n > 0) c /= (n + 1);
      f += c * pw[n];
      if (n >= 1) f1 += n * c * pw[n - 1];
      if (n >= 2) f2 += n * (n - 1) * c * pw[n - 2];
    }
  } else {
    const double e = std::exp(w);
    f = std::expm1(w) / w;
    f1 = (e * (w - 1.0) + 1.0) / (w * w);
    f2 = (e * (w * w - 2.0 * w + 2.0) - 2.0) / (w * w * w);
  }
}
// psi(v) = log(1+v)/v with psi(0) = 1, and its first two derivatives (v > -1).
inline void psi_fun(double v, double& f, double& f1, double& f2) {
  if (std::fabs(v) < 0.1) {
    // psi = sum (-v)^n/(n+1), psi' = sum_{n>=1} (-1)^n n v^(n-1)/(n+1), psi'' = sum_{n>=2} (-1)^n n(n-1) v^(n-2)/(n+1)
    double pw[40];
    pw[0] = 1.0;
    for (int n = 1; n < 40; ++n) pw[n] = pw[n - 1] * v;
    f = 0.0; f1 = 0.0; f2 = 0.0;
    double sgn = 1.0;
    for (int n = 0; n < 40; ++n) {
      f += sgn * pw[n] / (n + 1);
      if (n >= 1) f1 += sgn * n * pw[n - 1] / (n + 1);
      if (n >= 2) f2 += sgn * n * (n - 1) * pw[n - 2] / (n + 1);
      sgn = -sgn;
    }
  } else {
    const double l = std::log1p(v), q = 1.0 + v;
    f = l / v;
    f1 = (v / q - l) / (v * v);
    f2 = -1.0 / (v * q * q) - 2.0 / (v * v * q) + 2.0 * l / (v * v * v);
  }
}

// ---- one set of formulas for plain values (double) and for derivatives (D2) ----
inline double val(double a) { return a; }
inline double val(const D2& a) { return a.v; }
inline double unary(double, double f, double, double) { return f; }
inline D2 unary(const D2& a, double f, double f1, double f2) { return chain(a, f, f1, f2); }
template <class T> struct MakeVar;
template <> struct MakeVar<double> { static double make(int, double x) { return x; } };
template <> struct MakeVar<D2> { static D2 make(int i, double x) { return D2::variable(i, x); } };

template <class T> T t_exp(const T& a) { const double e = std::exp(val(a)); return unary(a, e, e, e); }
template <class T> T t_log(const T& a) { const double x = val(a); return unary(a, std::log(x), 1.0 / x, -1.0 / (x * x)); }
template <class T> T t_expm1(const T& a) { const double e = std::exp(val(a)); return unary(a, std::expm1(val(a)), e, e); }
template <class T> T t_neglog1m(const T& a) {            // -log(1 - a)
  const double q = 1.0 - val(a);
  return unary(a, -std::log1p(-val(a)), 1.0 / q, 1.0 / (q * q));
}
template <class T> T t_phi(const T& a) { double f, f1, f2; phi_fun(val(a), f, f1, f2); return unary(a, f, f1, f2); }
template <class T> T t_psi(const T& a) { double f, f1, f2; psi_fun(val(a), f, f1, f2); return unary(a, f, f1, f2); }
template <class T> T t_sigmoid(const T& a) {             // 1 / (1 + e^-a), stable for both signs
  const double x = val(a), s = x >= 0.0 ? 1.0 / (1.0 + std::exp(-x)) : std::exp(x) / (1.0 + std::exp(x));
  return unary(a, s, s * (1.0 - s), s * (1.0 - s) * (1.0 - 2.0 * s));
}
template <class T> T t_softplus(const T& a) {            // log(1 + e^a)
  const double x = val(a), s = x >= 0.0 ? 1.0 / (1.0 + std::exp(-x)) : std::exp(x) / (1.0 + std::exp(x));
  return unary(a, x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x)), s, s * (1.0 - s));
}
template <class T> T t_normcdf(const T& a) {             // standard normal distribution function
  const double z = val(a), d = std::exp(-0.5 * z * z) * 0.3989422804014327;
  return unary(a, 0.5 * std::erfc(-z * 0.7071067811865476), d, -z * d);
}
template <class T> T t_recip_neg(const T& a) {           // -1 / a
  const double x = val(a);
  return unary(a, -1.0 / x, 1.0 / (x * x), -2.0 / (x * x * x));
}

}  // namespace survcore
