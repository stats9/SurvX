// survcore: Fine-Gray proportional subdistribution hazards model.
//   Fine & Gray (1999) JASA 94:496-509, Section 4: weighted score function with
//   inverse-probability-of-censoring weights
//        w_i(t) = r_i(t) G(t) / G(X_i ^ t),
//   where G is the Kaplan-Meier estimate of the censoring distribution.
//   Subjects with a competing event stay in the risk set with weight
//   G(t-) / G(X_i-); subjects still under observation have weight 1.
//
// Algorithm: forward-backward scan, O(n p^2) per evaluation after one sort, in
// the spirit of Kawaguchi et al. (2021, R Journal 12(2), fastcmprsk). All
// running sums are accumulated by addition only.
//
// Conventions (to be cross-checked against cmprsk::crr on tied data):
//  * the censoring Kaplan-Meier is evaluated left-continuously, G(t-);
//  * ties among events of interest use the Breslow form.
//
// Robust variance: Fine & Gray (1999), Section 4, Sigma-hat built from
// eta_i + psi_i, where eta_i is the weighted score residual and psi_i corrects
// for the estimation of the censoring distribution. Implemented in
// robust_meat() with cumulative sums, O(n p) time.
#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "newton.hpp"

namespace survcore {

class FineGrayModel {
 public:
  // cause: 0 = censored, `interest` = event of interest, anything else = competing.
  FineGrayModel(std::size_t n, std::size_t p, const double* time, const int* cause,
                const double* X, int interest)
      : n_(n), p_(p), time_(time), X_(X), ord_(n), kind_(n), gminus_(n), jle_(n) {
    std::iota(ord_.begin(), ord_.end(), 0);
    std::stable_sort(ord_.begin(), ord_.end(),
                     [&](std::size_t a, std::size_t b) { return time_[a] < time_[b]; });
    for (std::size_t i = 0; i < n; ++i)
      kind_[i] = cause[i] == 0 ? 0 : (cause[i] == interest ? 1 : 2);
    // Left-continuous censoring Kaplan-Meier G(t-) and the distinct event times.
    double G = 1.0, at_risk = static_cast<double>(n);
    for (std::size_t i = 0; i < n;) {
      const double t = time_[ord_[i]];
      double cens = 0.0, d = 0.0;
      std::size_t j = i;
      for (; j < n && time_[ord_[j]] == t; ++j) {
        const std::size_t k = ord_[j];
        gminus_[k] = G;
        if (kind_[k] == 0) cens += 1.0;
        if (kind_[k] == 1) d += 1.0;
      }
      if (d > 0.0) { et_.push_back(t); ed_.push_back(d); eg_.push_back(G); }
      for (std::size_t q = i; q < j; ++q) jle_[ord_[q]] = et_.size();  // # event times <= X_i
      if (cens > 0.0) G *= 1.0 - cens / at_risk;
      at_risk -= static_cast<double>(j - i);
      i = j;
    }
  }

  std::size_t n_events() const {
    double d = 0.0;
    for (double x : ed_) d += x;
    return static_cast<std::size_t>(d);
  }
  const std::vector<double>& event_times() const { return et_; }

  void eval(const std::vector<double>& beta, Eval& out, bool need_hess) const {
    const std::size_t p = p_, m = et_.size();
    out.ll = 0.0;
    out.g.assign(p, 0.0);
    if (out.H.n != p) out.H = Mat(p);
    out.H.zero();
    std::vector<double> e(n_);
    double c = -INFINITY;
    for (std::size_t i = 0; i < n_; ++i) {
      const double* x = X_ + i * p;
      double s = 0.0;
      for (std::size_t j = 0; j < p; ++j) s += x[j] * beta[j];
      e[i] = s;
      c = std::fmax(c, s);
    }
    for (std::size_t i = 0; i < n_; ++i) {
      if (kind_[i] == 1) {
        out.ll += e[i];
        const double* x = X_ + i * p;
        for (std::size_t a = 0; a < p; ++a) out.g[a] += x[a];
      }
      e[i] = std::exp(e[i] - c);
    }
    std::vector<double> S0, S1;
    risk_sums(e, S0, S1);
    // Log-likelihood, gradient and the outer-product part of the information.
    std::vector<double> P(m + 1, 0.0), Q(m + 1, 0.0), a1(p);
    for (std::size_t jj = 0; jj < m; ++jj) {
      const double d = ed_[jj];
      out.ll -= d * (c + std::log(S0[jj]));
      for (std::size_t a = 0; a < p; ++a) {
        a1[a] = S1[jj * p + a] / S0[jj];
        out.g[a] -= d * a1[a];
      }
      if (need_hess)
        for (std::size_t a = 0; a < p; ++a)
          for (std::size_t b = 0; b <= a; ++b) out.H(a, b) -= d * a1[a] * a1[b];
      P[jj + 1] = P[jj] + d / S0[jj];                       // sum over t_j <= x
    }
    if (need_hess) {
      for (std::size_t jj = m; jj-- > 0;) Q[jj] = Q[jj + 1] + ed_[jj] * eg_[jj] / S0[jj];  // sum over t_j > x
      for (std::size_t i = 0; i < n_; ++i) {
        const std::size_t jl = jle_[i];
        double B = P[jl];
        if (kind_[i] == 2) B += Q[jl] / gminus_[i];
        const double eb = e[i] * B;
        if (eb == 0.0) continue;
        const double* x = X_ + i * p;
        for (std::size_t a = 0; a < p; ++a) {
          const double v = eb * x[a];
          for (std::size_t b = 0; b <= a; ++b) out.H(a, b) += v * x[b];
        }
      }
      for (std::size_t a = 0; a < p; ++a)
        for (std::size_t b = a + 1; b < p; ++b) out.H(a, b) = out.H(b, a);
    }
  }

 private:
  // S0_j and S1_j of the weighted risk set at every distinct event time, given
  // e_i = exp(eta_i - c). Backward scan for subjects still under observation,
  // forward scan for competing-event subjects (IPCW weight G(t_j-) / G(X_i-)).
  void risk_sums(const std::vector<double>& e, std::vector<double>& S0, std::vector<double>& S1) const {
    const std::size_t p = p_, m = et_.size();
    S0.assign(m, 0.0);
    S1.assign(m * p, 0.0);
    std::vector<double> run1(p, 0.0);
    double run0 = 0.0;
    {
      std::size_t i = n_;
      for (std::size_t jj = m; jj-- > 0;) {
        while (i > 0 && time_[ord_[i - 1]] >= et_[jj]) {
          const std::size_t k = ord_[--i];
          const double* x = X_ + k * p;
          run0 += e[k];
          for (std::size_t a = 0; a < p; ++a) run1[a] += e[k] * x[a];
        }
        S0[jj] = run0;
        for (std::size_t a = 0; a < p; ++a) S1[jj * p + a] = run1[a];
      }
    }
    // Forward scan: add G(t_j-) * sum over competing subjects with X_i < t_j of e_i / G(X_i-).
    {
      run0 = 0.0;
      std::fill(run1.begin(), run1.end(), 0.0);
      std::size_t i = 0;
      for (std::size_t jj = 0; jj < m; ++jj) {
        while (i < n_ && time_[ord_[i]] < et_[jj]) {
          const std::size_t k = ord_[i++];
          if (kind_[k] != 2) continue;
          const double v = e[k] / gminus_[k];
          const double* x = X_ + k * p;
          run0 += v;
          for (std::size_t a = 0; a < p; ++a) run1[a] += v * x[a];
        }
        S0[jj] += eg_[jj] * run0;
        for (std::size_t a = 0; a < p; ++a) S1[jj * p + a] += eg_[jj] * run1[a];
      }
    }
  }

  std::vector<double> scaled_risk(const std::vector<double>& beta, const double* zref, double& c) const {
    std::vector<double> e(n_);
    c = -INFINITY;
    for (std::size_t i = 0; i < n_; ++i) {
      const double* x = X_ + i * p_;
      double v = 0.0;
      for (std::size_t j = 0; j < p_; ++j) v += (x[j] - (zref ? zref[j] : 0.0)) * beta[j];
      e[i] = v;
      c = std::fmax(c, v);
    }
    for (std::size_t i = 0; i < n_; ++i) e[i] = std::exp(e[i] - c);
    return e;
  }

 public:
  // Breslow-type baseline cumulative subdistribution hazard at covariate value
  // zref (null = zero vector): Fine & Gray (1999), Section 5.
  void baseline(const std::vector<double>& beta, const double* zref, std::vector<double>& cumhaz) const {
    double c;
    std::vector<double> e = scaled_risk(beta, zref, c), S0, S1;
    risk_sums(e, S0, S1);
    cumhaz.resize(et_.size());
    double H = 0.0;
    for (std::size_t j = 0; j < et_.size(); ++j) { H += ed_[j] / S0[j] * std::exp(-c); cumhaz[j] = H; }
  }

  // "Meat" of the sandwich: sum_i (eta_i + psi_i)(eta_i + psi_i)^T.
  // The robust variance is  I^{-1} * meat * I^{-1}  with I the information from eval().
  void robust_meat(const std::vector<double>& beta, Mat& meat) const {
    const std::size_t p = p_, m = et_.size(), n = n_;
    double c;
    std::vector<double> e = scaled_risk(beta, nullptr, c), S0, S1;
    risk_sums(e, S0, S1);
    // prefix sums over event times t_j <= x and suffix sums over t_j > x
    std::vector<double> P0(m + 1, 0.0), P1((m + 1) * p, 0.0), Q0(m + 1, 0.0), Q1((m + 1) * p, 0.0);
    for (std::size_t j = 0; j < m; ++j) {
      const double f = ed_[j] / S0[j];
      P0[j + 1] = P0[j] + f;
      for (std::size_t a = 0; a < p; ++a) P1[(j + 1) * p + a] = P1[j * p + a] + f * S1[j * p + a] / S0[j];
    }
    for (std::size_t j = m; j-- > 0;) {
      const double f = ed_[j] * eg_[j] / S0[j];
      Q0[j] = Q0[j + 1] + f;
      for (std::size_t a = 0; a < p; ++a) Q1[j * p + a] = Q1[(j + 1) * p + a] + f * S1[j * p + a] / S0[j];
    }
    std::vector<double> res(n * p, 0.0);      // eta_i, then eta_i + psi_i
    for (std::size_t i = 0; i < n; ++i) {
      const double* x = X_ + i * p;
      const std::size_t jl = jle_[i];
      double* r = &res[i * p];
      if (kind_[i] == 1)                       // event of interest at t_{jl-1} = X_i
        for (std::size_t a = 0; a < p; ++a) r[a] += x[a] - S1[(jl - 1) * p + a] / S0[jl - 1];
      for (std::size_t a = 0; a < p; ++a) r[a] -= e[i] * (x[a] * P0[jl] - P1[jl * p + a]);
      if (kind_[i] == 2) {
        const double v = e[i] / gminus_[i];
        for (std::size_t a = 0; a < p; ++a) r[a] -= v * (x[a] * Q0[jl] - Q1[jl * p + a]);
      }
    }
    // psi_i: forward scan over distinct observed times; h(u) is non-zero only at censoring times.
    std::vector<double> A1(p, 0.0), Hc(p, 0.0), h(p);
    double A0 = 0.0, at_risk = static_cast<double>(n);
    for (std::size_t i = 0; i < n;) {
      const double t = time_[ord_[i]];
      double cens = 0.0;
      std::size_t j = i;
      for (; j < n && time_[ord_[j]] == t; ++j) {
        const std::size_t k = ord_[j];
        if (kind_[k] == 0) cens += 1.0;
        if (kind_[k] == 2) {                   // competing subjects with X_i <= u
          const double v = e[k] / gminus_[k];
          const double* x = X_ + k * p;
          A0 += v;
          for (std::size_t a = 0; a < p; ++a) A1[a] += v * x[a];
        }
      }
      if (cens > 0.0) {
        const std::size_t jl = jle_[ord_[i]];  // number of event times <= u
        for (std::size_t a = 0; a < p; ++a) {
          h[a] = (A1[a] * Q0[jl] - A0 * Q1[jl * p + a]) / at_risk;
          Hc[a] += h[a] * cens / at_risk;
        }
      }
      for (std::size_t q = i; q < j; ++q) {
        const std::size_t k = ord_[q];
        double* r = &res[k * p];
        for (std::size_t a = 0; a < p; ++a) r[a] += ((kind_[k] == 0 && cens > 0.0) ? h[a] : 0.0) - Hc[a];
      }
      at_risk -= static_cast<double>(j - i);
      i = j;
    }
    meat = Mat(p);
    for (std::size_t i = 0; i < n; ++i) {
      const double* r = &res[i * p];
      for (std::size_t a = 0; a < p; ++a)
        for (std::size_t b = 0; b <= a; ++b) meat(a, b) += r[a] * r[b];
    }
    for (std::size_t a = 0; a < p; ++a)
      for (std::size_t b = a + 1; b < p; ++b) meat(a, b) = meat(b, a);
  }

 private:
  std::size_t n_, p_;
  const double* time_;
  const double* X_;
  std::vector<std::size_t> ord_;
  std::vector<int> kind_;            // 0 censored, 1 interest, 2 competing
  std::vector<double> gminus_;       // G(X_i-)
  std::vector<std::size_t> jle_;     // number of distinct event times <= X_i
  std::vector<double> et_, ed_, eg_; // distinct event times, counts, G(t_j-)
};

}  // namespace survcore
