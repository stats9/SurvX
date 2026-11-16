// survcore: Cox proportional hazards partial likelihood.
//   Cox (1972) JRSS-B 34:187-220; Cox (1975) Biometrika 62:269-276
//   Breslow (1974) Biometrics 30:89-99      -- ties = BRESLOW
//   Efron (1977) JASA 72:557-565            -- ties = EFRON
// Extensions (all optional):
//   * strata: a separate baseline hazard in every stratum, common coefficients
//     (Kalbfleisch & Prentice 2002, Sec. 4.4);
//   * (entry, time] intervals: delayed entry and time-varying covariates in the
//     counting process form of Andersen & Gill (1982, Ann Stat 10:1100-1120);
//     subject i is at risk at t when entry_i < t <= time_i;
//   * score residuals for the robust (sandwich) variance of Lin & Wei (1989,
//     JASA 84:1074-1078), summed within clusters by the caller.
// A cause-specific Cox model (Prentice et al. 1978, Biometrics 34:541-554) is
// this model with status_i = 1{cause_i == k}: other causes count as censored.
//
// Numerical notes:
//  * risk-set sums are accumulated scanning from the largest time downwards;
//    without delayed entry this uses additions only (no cancellation). With
//    entry times, subjects leave the risk set by subtraction, as in the
//    counting-process algorithm of survival::agreg.fit;
//  * the linear predictor is shifted by its maximum before exponentiation, so
//    exp() cannot overflow; the log-likelihood is corrected exactly.
// Cost of one evaluation: O(n p^2) after an O(n log n) sort.
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
  // X is row-major n x p. weights, entry and strata may be null.
  CoxModel(std::size_t n, std::size_t p, const double* time, const int* status,
           const double* X, const double* weights, Ties ties,
           const double* entry = nullptr, const int* strata = nullptr)
      : n_(n), p_(p), time_(time), status_(status), X_(X), w_(weights), ties_(ties),
        entry_(entry), strata_(strata), ord_(n) {
    std::iota(ord_.begin(), ord_.end(), 0);
    std::stable_sort(ord_.begin(), ord_.end(), [&](std::size_t a, std::size_t b) {
      if (st(a) != st(b)) return st(a) < st(b);
      return time_[a] > time_[b];
    });
    if (entry_) {
      ord_entry_.resize(n);
      std::iota(ord_entry_.begin(), ord_entry_.end(), 0);
      std::stable_sort(ord_entry_.begin(), ord_entry_.end(), [&](std::size_t a, std::size_t b) {
        if (st(a) != st(b)) return st(a) < st(b);
        return entry_[a] > entry_[b];
      });
    }
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
    std::vector<double> eta;
    const double c = linear_predictor(beta, nullptr, eta);
    std::vector<double> a1(p);
    scan(eta, c, need_hess, [&](const Sums& s, std::size_t, double, const std::vector<std::size_t>& dead) {
      double dw = 0.0;
      for (std::size_t k : dead) {
        const double w = wt(k);
        dw += w;
        out.ll += w * eta[k];
        const double* x = X_ + k * p;
        for (std::size_t a = 0; a < p; ++a) out.g[a] += w * x[a];
      }
      const std::size_t d = dead.size();
      const std::size_t terms = (ties_ == EFRON) ? d : 1;
      const double wq = (ties_ == EFRON) ? dw / static_cast<double>(d) : dw;
      for (std::size_t q = 0; q < terms; ++q) {
        const double f = (ties_ == EFRON) ? static_cast<double>(q) / static_cast<double>(d) : 0.0;
        const double a0 = s.S0 - f * s.D0;
        out.ll -= wq * (c + std::log(a0));
        for (std::size_t a = 0; a < p; ++a) {
          a1[a] = (s.S1[a] - f * s.D1[a]) / a0;
          out.g[a] -= wq * a1[a];
        }
        if (need_hess)
          for (std::size_t a = 0; a < p; ++a)
            for (std::size_t b = 0; b <= a; ++b)
              out.H(a, b) += wq * ((s.S2(a, b) - f * s.D2(a, b)) / a0 - a1[a] * a1[b]);
      }
    });
    if (need_hess)
      for (std::size_t a = 0; a < p; ++a)
        for (std::size_t b = a + 1; b < p; ++b) out.H(a, b) = out.H(b, a);
  }

  // Baseline cumulative hazard at covariate value zref (null = zero vector), at the distinct
  // event times of every stratum, ordered by stratum and then by time.
  //   Breslow (1972, discussion of Cox) / Breslow (1974): dLambda0 = d / S0
  //   Efron-type increment for ties = EFRON, as in Therneau & Grambsch (2000), Sec. 3.3.
  void baseline(const std::vector<double>& beta, const double* zref, std::vector<double>& times,
                std::vector<double>& cumhaz, std::vector<int>* strata_out = nullptr) const {
    std::vector<double> eta;
    const double c = linear_predictor(beta, zref, eta);
    std::vector<double> t_desc, inc_desc;
    std::vector<int> s_desc;
    scan(eta, c, false, [&](const Sums& s, std::size_t stratum, double t, const std::vector<std::size_t>& dead) {
      double dw = 0.0;
      for (std::size_t k : dead) dw += wt(k);
      const std::size_t d = dead.size();
      double inc = 0.0;
      if (ties_ == EFRON)
        for (std::size_t q = 0; q < d; ++q)
          inc += (dw / static_cast<double>(d)) / (s.S0 - static_cast<double>(q) / static_cast<double>(d) * s.D0);
      else
        inc = dw / s.S0;
      t_desc.push_back(t);
      inc_desc.push_back(inc * std::exp(-c));
      s_desc.push_back(static_cast<int>(stratum));
    });
    // the scan visits strata in increasing order and times in decreasing order within a stratum
    times.clear(); cumhaz.clear();
    if (strata_out) strata_out->clear();
    std::size_t i = 0;
    while (i < t_desc.size()) {
      std::size_t j = i;
      while (j < t_desc.size() && s_desc[j] == s_desc[i]) ++j;
      double H = 0.0;
      for (std::size_t q = j; q-- > i;) {
        H += inc_desc[q];
        times.push_back(t_desc[q]);
        cumhaz.push_back(H);
        if (strata_out) strata_out->push_back(s_desc[q]);
      }
      i = j;
    }
  }

  // Score residuals U_i (n x p, row-major): the contribution of row i to the score at beta,
  // U_i = sum over event times t in (entry_i, time_i] of (x_i - xbar(t)) dM_i(t), in the Efron or
  // Breslow form (Therneau & Grambsch 2000, Sec. 7.2). Their sum over all rows is the score.
  void score_residuals(const std::vector<double>& beta, std::vector<double>& U) const {
    const std::size_t p = p_;
    std::vector<double> eta;
    const double c = linear_predictor(beta, nullptr, eta);
    struct EvTime { std::size_t stratum; double t, A, AD; std::vector<double> B, BD, abar; };
    std::vector<EvTime> ev;
    std::vector<double> a1(p);
    scan(eta, c, false, [&](const Sums& s, std::size_t stratum, double t, const std::vector<std::size_t>& dead) {
      double dw = 0.0;
      for (std::size_t k : dead) dw += wt(k);
      const std::size_t d = dead.size();
      const std::size_t terms = (ties_ == EFRON) ? d : 1;
      const double wq = (ties_ == EFRON) ? dw / static_cast<double>(d) : dw;
      EvTime e{stratum, t, 0.0, 0.0, std::vector<double>(p, 0.0), std::vector<double>(p, 0.0), std::vector<double>(p, 0.0)};
      for (std::size_t q = 0; q < terms; ++q) {
        const double f = (ties_ == EFRON) ? static_cast<double>(q) / static_cast<double>(d) : 0.0;
        const double a0 = s.S0 - f * s.D0;
        e.A += wq / a0;
        e.AD += wq * (1.0 - f) / a0;
        for (std::size_t a = 0; a < p; ++a) {
          a1[a] = (s.S1[a] - f * s.D1[a]) / a0;
          e.B[a] += wq * a1[a] / a0;
          e.BD[a] += wq * (1.0 - f) * a1[a] / a0;
          e.abar[a] += a1[a] / static_cast<double>(terms);
        }
      }
      ev.push_back(std::move(e));
    });
    // event times of each stratum in increasing order, with prefix sums of A and B
    std::size_t nstr = 0;
    for (const auto& e : ev) nstr = std::max(nstr, e.stratum + 1);
    std::vector<std::vector<std::size_t>> idx(nstr);
    for (std::size_t q = ev.size(); q-- > 0;) idx[ev[q].stratum].push_back(q);   // ascending in time
    std::vector<std::vector<double>> cA(nstr), cB(nstr);                         // cB: (m+1) x p
    for (std::size_t s = 0; s < nstr; ++s) {
      const std::size_t m = idx[s].size();
      cA[s].assign(m + 1, 0.0); cB[s].assign((m + 1) * p, 0.0);
      for (std::size_t q = 0; q < m; ++q) {
        const EvTime& e = ev[idx[s][q]];
        cA[s][q + 1] = cA[s][q] + e.A;
        for (std::size_t a = 0; a < p; ++a) cB[s][(q + 1) * p + a] = cB[s][q * p + a] + e.B[a];
      }
    }
    U.assign(n_ * p, 0.0);
    for (std::size_t i = 0; i < n_; ++i) {
      const std::size_t s = st(i);
      if (s >= nstr) continue;
      const auto& ix = idx[s];
      auto tval = [&](std::size_t q) { return ev[ix[q]].t; };
      // event times t with entry_i < t <= time_i: positions [lo, hi)
      std::size_t lo = 0, hi = ix.size();
      { std::size_t a = 0, b = ix.size(); while (a < b) { std::size_t mid = (a + b) / 2; if (tval(mid) <= time_[i]) a = mid + 1; else b = mid; } hi = a; }
      if (entry_) { std::size_t a = 0, b = hi; while (a < b) { std::size_t mid = (a + b) / 2; if (tval(mid) <= entry_[i]) a = mid + 1; else b = mid; } lo = a; }
      const double r = wt(i) * std::exp(eta[i] - c);
      const double* x = X_ + i * p;
      double* u = U.data() + i * p;
      const double sA = cA[s][hi] - cA[s][lo];
      for (std::size_t a = 0; a < p; ++a) u[a] = -r * (x[a] * sA - (cB[s][hi * p + a] - cB[s][lo * p + a]));
      if (status_[i] != 0 && hi > lo && tval(hi - 1) == time_[i]) {
        const EvTime& e = ev[ix[hi - 1]];
        for (std::size_t a = 0; a < p; ++a)
          u[a] += -r * (x[a] * (e.AD - e.A) - (e.BD[a] - e.B[a])) + wt(i) * (x[a] - e.abar[a]);
      }
    }
  }

 private:
  struct Sums {
    double S0 = 0.0, D0 = 0.0;
    std::vector<double> S1, D1;
    Mat S2, D2;
  };

  std::size_t st(std::size_t i) const { return strata_ ? static_cast<std::size_t>(strata_[i]) : 0; }
  double wt(std::size_t i) const { return w_ ? w_[i] : 1.0; }

  double linear_predictor(const std::vector<double>& beta, const double* zref, std::vector<double>& eta) const {
    eta.assign(n_, 0.0);
    double c = -INFINITY;
    for (std::size_t i = 0; i < n_; ++i) {
      const double* x = X_ + i * p_;
      double e = 0.0;
      for (std::size_t j = 0; j < p_; ++j) e += (x[j] - (zref ? zref[j] : 0.0)) * beta[j];
      eta[i] = e;
      c = std::fmax(c, e);
    }
    return c;
  }

  void add(Sums& s, std::size_t k, const std::vector<double>& eta, double c, double sign, bool hess) const {
    const double* x = X_ + k * p_;
    const double r = sign * wt(k) * std::exp(eta[k] - c);
    s.S0 += r;
    for (std::size_t a = 0; a < p_; ++a) s.S1[a] += r * x[a];
    if (hess)
      for (std::size_t a = 0; a < p_; ++a) {
        const double rx = r * x[a];
        for (std::size_t b = 0; b <= a; ++b) s.S2(a, b) += rx * x[b];
      }
  }

  // Visit every event time of every stratum (strata in increasing order, times decreasing),
  // with the risk-set sums of that time and the rows with an event there.
  template <class F>
  void scan(const std::vector<double>& eta, double c, bool hess, F&& visit) const {
    const std::size_t p = p_;
    Sums s;
    s.S1.assign(p, 0.0); s.D1.assign(p, 0.0);
    if (hess) { s.S2 = Mat(p); s.D2 = Mat(p); }
    std::vector<std::size_t> dead;
    std::size_t i = 0, e = 0;
    while (i < n_) {
      const std::size_t stratum = st(ord_[i]);
      s.S0 = 0.0; std::fill(s.S1.begin(), s.S1.end(), 0.0);
      if (hess) s.S2.zero();
      while (entry_ && e < n_ && st(ord_entry_[e]) < stratum) ++e;
      while (i < n_ && st(ord_[i]) == stratum) {
        const double t = time_[ord_[i]];
        s.D0 = 0.0; std::fill(s.D1.begin(), s.D1.end(), 0.0);
        if (hess) s.D2.zero();
        dead.clear();
        std::size_t j = i;
        for (; j < n_ && st(ord_[j]) == stratum && time_[ord_[j]] == t; ++j) {
          const std::size_t k = ord_[j];
          add(s, k, eta, c, 1.0, hess);
          if (status_[k] != 0) {
            dead.push_back(k);
            if (ties_ == EFRON) {
              const double* x = X_ + k * p;
              const double r = wt(k) * std::exp(eta[k] - c);
              s.D0 += r;
              for (std::size_t a = 0; a < p; ++a) s.D1[a] += r * x[a];
              if (hess)
                for (std::size_t a = 0; a < p; ++a) {
                  const double rx = r * x[a];
                  for (std::size_t b = 0; b <= a; ++b) s.D2(a, b) += rx * x[b];
                }
            }
          }
        }
        // rows whose interval starts at or after t are not at risk at t
        if (entry_)
          while (e < n_ && st(ord_entry_[e]) == stratum && entry_[ord_entry_[e]] >= t) {
            add(s, ord_entry_[e], eta, c, -1.0, hess);
            ++e;
          }
        if (!dead.empty()) visit(s, stratum, t, dead);
        i = j;
      }
    }
  }

  std::size_t n_, p_;
  const double* time_;
  const int* status_;
  const double* X_;
  const double* w_;
  Ties ties_;
  const double* entry_;
  const int* strata_;
  std::vector<std::size_t> ord_, ord_entry_;
};

}  // namespace survcore
