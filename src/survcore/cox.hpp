// survcore: Cox proportional hazards partial likelihood (right-censored data).
//   Cox (1972) JRSS-B 34:187-220; Cox (1975) Biometrika 62:269-276
//   Breslow (1974) Biometrics 30:89-99      -- ties = BRESLOW
//   Efron (1977) JASA 72:557-565            -- ties = EFRON
// A cause-specific Cox model (Prentice et al. 1978, Biometrics 34:541-554) is
// this model with status_i = 1{cause_i == k}: other causes count as censored.
//
// Numerical notes:
//  * risk-set sums are accumulated by addition only, scanning from the largest
//    time downwards (no subtraction, hence no cancellation);
//  * the linear predictor is shifted by its maximum before exponentiation, so
//    exp() cannot overflow; the log-likelihood is corrected exactly.
// Cost of one evaluation: O(n p^2) after a single O(n log n) sort.
#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "newton.hpp"

namespace survcore {

enum Ties : int { BRESLOW = 0, EFRON = 1 };

class CoxModel {
 public:
  // X is row-major n x p. weights may be null (all 1).
  CoxModel(std::size_t n, std::size_t p, const double* time, const int* status,
           const double* X, const double* weights, Ties ties)
      : n_(n), p_(p), time_(time), status_(status), X_(X), w_(weights), ties_(ties), ord_(n) {
    std::iota(ord_.begin(), ord_.end(), 0);
    std::stable_sort(ord_.begin(), ord_.end(),
                     [&](std::size_t a, std::size_t b) { return time_[a] > time_[b]; });
  }

  std::size_t n_events() const {
    std::size_t d = 0;
    for (std::size_t i = 0; i < n_; ++i) d += status_[i] != 0;
    return d;
  }

  void eval(const std::vector<double>& beta, Eval& out, bool need_hess) const {
    const std::size_t p = p_;
    out.ll = 0.0;
    out.g.assign(p, 0.0);
    if (out.H.n != p) out.H = Mat(p);
    out.H.zero();
    std::vector<double> eta(n_);
    double c = -INFINITY;
    for (std::size_t i = 0; i < n_; ++i) {
      const double* x = X_ + i * p;
      double e = 0.0;
      for (std::size_t j = 0; j < p; ++j) e += x[j] * beta[j];
      eta[i] = e;
      c = std::fmax(c, e);
    }
    double S0 = 0.0;
    std::vector<double> S1(p, 0.0), D1(p), a1(p);
    Mat S2(p), D2(p);
    for (std::size_t i = 0; i < n_;) {
      const double t = time_[ord_[i]];
      double D0 = 0.0, dw = 0.0;
      std::size_t d = 0;
      std::fill(D1.begin(), D1.end(), 0.0);
      if (need_hess) D2.zero();
      std::size_t j = i;
      for (; j < n_ && time_[ord_[j]] == t; ++j) {
        const std::size_t k = ord_[j];
        const double* x = X_ + k * p;
        const double w = w_ ? w_[k] : 1.0;
        const double r = w * std::exp(eta[k] - c);
        S0 += r;
        for (std::size_t a = 0; a < p; ++a) S1[a] += r * x[a];
        if (need_hess)
          for (std::size_t a = 0; a < p; ++a) {
            const double rx = r * x[a];
            for (std::size_t b = 0; b <= a; ++b) S2(a, b) += rx * x[b];
          }
        if (status_[k] != 0) {
          ++d;
          dw += w;
          out.ll += w * eta[k];
          for (std::size_t a = 0; a < p; ++a) out.g[a] += w * x[a];
          if (ties_ == EFRON) {
            D0 += r;
            for (std::size_t a = 0; a < p; ++a) D1[a] += r * x[a];
            if (need_hess)
              for (std::size_t a = 0; a < p; ++a) {
                const double rx = r * x[a];
                for (std::size_t b = 0; b <= a; ++b) D2(a, b) += rx * x[b];
              }
          }
        }
      }
      if (d > 0) {
        const std::size_t terms = (ties_ == EFRON) ? d : 1;
        const double wt = (ties_ == EFRON) ? dw / static_cast<double>(d) : dw;
        for (std::size_t q = 0; q < terms; ++q) {
          const double f = (ties_ == EFRON) ? static_cast<double>(q) / static_cast<double>(d) : 0.0;
          const double a0 = S0 - f * D0;
          out.ll -= wt * (c + std::log(a0));
          for (std::size_t a = 0; a < p; ++a) {
            a1[a] = (S1[a] - f * D1[a]) / a0;
            out.g[a] -= wt * a1[a];
          }
          if (need_hess)
            for (std::size_t a = 0; a < p; ++a)
              for (std::size_t b = 0; b <= a; ++b)
                out.H(a, b) += wt * ((S2(a, b) - f * D2(a, b)) / a0 - a1[a] * a1[b]);
        }
      }
      i = j;
    }
    if (need_hess)
      for (std::size_t a = 0; a < p; ++a)
        for (std::size_t b = a + 1; b < p; ++b) out.H(a, b) = out.H(b, a);
  }

  // Baseline cumulative hazard at covariate value zref (null = zero vector), at
  // the distinct event times in increasing order.
  //   Breslow (1972, discussion of Cox) / Breslow (1974): dLambda0 = d / S0
  //   Efron-type increment for ties = EFRON, as in Therneau & Grambsch (2000), Sec. 3.3.
  void baseline(const std::vector<double>& beta, const double* zref,
                std::vector<double>& times, std::vector<double>& cumhaz) const {
    const std::size_t p = p_;
    std::vector<double> eta(n_);
    double c = -INFINITY;
    for (std::size_t i = 0; i < n_; ++i) {
      const double* x = X_ + i * p;
      double e = 0.0;
      for (std::size_t j = 0; j < p; ++j) e += (x[j] - (zref ? zref[j] : 0.0)) * beta[j];
      eta[i] = e;
      c = std::fmax(c, e);
    }
    std::vector<double> t_desc, inc_desc;
    double S0 = 0.0;
    for (std::size_t i = 0; i < n_;) {
      const double t = time_[ord_[i]];
      double D0 = 0.0, dw = 0.0;
      std::size_t d = 0, j = i;
      for (; j < n_ && time_[ord_[j]] == t; ++j) {
        const std::size_t k = ord_[j];
        const double w = w_ ? w_[k] : 1.0, r = w * std::exp(eta[k] - c);
        S0 += r;
        if (status_[k] != 0) { ++d; dw += w; D0 += r; }
      }
      if (d > 0) {
        double inc = 0.0;
        if (ties_ == EFRON)
          for (std::size_t q = 0; q < d; ++q)
            inc += (dw / static_cast<double>(d)) / (S0 - static_cast<double>(q) / static_cast<double>(d) * D0);
        else
          inc = dw / S0;
        t_desc.push_back(t);
        inc_desc.push_back(inc * std::exp(-c));
      }
      i = j;
    }
    times.assign(t_desc.rbegin(), t_desc.rend());
    cumhaz.resize(times.size());
    double H = 0.0;
    for (std::size_t q = 0; q < times.size(); ++q) { H += inc_desc[inc_desc.size() - 1 - q]; cumhaz[q] = H; }
  }

 private:
  std::size_t n_, p_;
  const double* time_;
  const int* status_;
  const double* X_;
  const double* w_;
  Ties ties_;
  std::vector<std::size_t> ord_;
};

}  // namespace survcore
