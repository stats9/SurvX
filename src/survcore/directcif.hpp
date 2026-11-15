// survcore: direct parametric models for the cumulative incidence function (CIF).
//
// Likelihood (all causes jointly; censored subjects contribute 1 - sum_k F_k):
//   Jeong & Fine (2006), Appl Statist 55:187-200, eq. (11);
//   Jeong & Fine (2007), Biostatistics 8:184-196, eq. (4.1)-(4.2).
// Link: the generalized odds rate (GOR) transformation, Jeong & Fine (2007) eq. (3.5),
//   after Dabrowska & Doksum (1988), JASA 83:744-749. alpha -> 0 is the proportional
//   (subdistribution) hazards model, alpha = 1 the proportional odds model.
//
// Two ways a baseline enters, both reducing to  F = 1 - {1 + alpha e^eta u}^(-1/alpha):
//
//  (A) baselines given as a cumulative subdistribution hazard u(t), as in Jeong & Fine
//      (2007) eq. (3.9):  F_k(t; z) = 1 - {1 + alpha e^eta u_k(t)}^(-1/alpha).
//        GOMPERTZ   u = tau (e^(rho t) - 1)/rho              Jeong & Fine (2006) eq. (6); (2007) eq. (3.7)
//        GOMPERTZ3  u = tau/(rho zeta) (e^(zeta e^(rho t)) - e^zeta)
//                                                             Haile, Jeong, Chen & Cheng (2016), J Appl Stat
//                                                             43:2239-2253, eq. (3)-(4); zeta = 0 is GOMPERTZ
//        WEIBULL    u = (lambda t)^kappa                      the Weibull form used for comparison in
//                                                             Jeong & Fine (2006), Sec. 2.1 and 5
//        SPLINE     log u = s(log t; gamma), a restricted cubic spline
//                                                             Royston & Parmar (2002), Stat Med 21:2175-2197,
//                                                             Sec. 2.2; Mozumder, Rutherford & Lambert (2018),
//                                                             Stat Med 37:82-97, eq. (12)-(19). With alpha = 0
//                                                             this is their log-cumulative SDH model (16)-(18),
//                                                             with alpha = 1 their logit model (19).
//  (B) baselines given as a baseline CIF F_0(t) = F(t; z = 0), with the regression
//      g{F(t; z)} = g{F_0(t)} + z'beta of Shi, Cheng & Jeong (2013), Biom J 55:82-96,
//      eq. (2) and (4). This is form (A) with u = [(1 - F_0)^(-alpha) - 1]/alpha.
//        LOGISTIC3  F_0 = p [e^(b(t-c)) - e^(-bc)] / [1 + e^(b(t-c))]
//                                                             Cheng (2009), Int J Biostat 5(1):29, eq. (2.1);
//                                                             Shi et al. (2013) eq. (3)
//        MIX_WEIBULL      F_0 = p {1 - exp(-(lambda t)^kappa)}  Cheng (2009) eq. (2.4), the "conditional
//                                                             Weibull" CIF: P(cause) x P(T <= t | cause), the
//                                                             decomposition of Larson & Dinse (1985), Appl
//                                                             Statist 34:201-211
//        MIX_LOGLOGISTIC  F_0 = p (lambda t)^kappa / {1 + (lambda t)^kappa}   same construction; no equation
//        MIX_LOGNORMAL    F_0 = p Phi{(log t - mu)/sigma}                       in the cited papers
//      Note: this is not the regression model of Larson & Dinse (1985), who put a logistic
//      model on P(cause) and a piecewise-exponential model on the conditional hazards.
//
// Parameters per cause, in this order: baseline parameters (unconstrained scale, see
// base_names()), [alpha_k if estimated], beta_k (p values).
//
// Numerical formulation, identical to the formulas above but free of cancellation:
//     v = alpha e^eta u,  A = e^eta u psi(v),  psi(v) = log(1+v)/v,  F = -expm1(-A),
//     log f = eta + log u' - A - log1p(v);     for (B): u = L phi(alpha L), L = -log(1 - F_0).
// Gradient and Hessian are exact (second-order forward-mode AD on at most five local
// variables per cause, then the chain rule for the linear forms z'beta and the spline).
//
// The models do not force sum_k F_k <= 1 (Shi, Cheng & Jeong 2013). Parameter values for
// which a censored subject would have 1 - sum_k F_k(x_i) <= 0 have likelihood zero and are
// rejected, never clamped.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "dual2.hpp"
#include "newton.hpp"
#include "nonparametric.hpp"

namespace survcore {

enum Baseline : int {
  BASE_GOMPERTZ = 0, BASE_GOMPERTZ3 = 1, BASE_WEIBULL = 2, BASE_MIX_WEIBULL = 3,
  BASE_MIX_LOGLOGISTIC = 4, BASE_MIX_LOGNORMAL = 5, BASE_LOGISTIC3 = 6, BASE_SPLINE = 7
};
constexpr int BASE_COUNT = 8;

inline bool base_needs_positive_time(int b) {
  return b == BASE_WEIBULL || b == BASE_MIX_WEIBULL || b == BASE_MIX_LOGLOGISTIC || b == BASE_MIX_LOGNORMAL || b == BASE_SPLINE;
}

struct DirectSpec {
  int K = 2;                                   // number of causes
  std::size_t p = 0;                           // covariates per cause
  std::vector<int> alpha_free;                 // per cause: 1 = alpha estimated
  std::vector<double> alpha_value;             // per cause: fixed value or start value
  std::vector<int> baseline;                   // per cause: Baseline
  std::vector<std::vector<double>> knots;      // per cause (SPLINE): knots in log time, boundary knots first and last;
                                               // empty = 1 d.f. (log u linear in log t)
  std::size_t nbase(int k) const {
    switch (baseline[k]) {
      case BASE_GOMPERTZ: case BASE_WEIBULL: return 2;
      case BASE_SPLINE: return knots[k].size() >= 2 ? knots[k].size() : 2;
      default: return 3;
    }
  }
  std::size_t block(int k) const { return nbase(k) + (alpha_free[k] ? 1 : 0) + p; }
  std::size_t offset(int k) const { std::size_t o = 0; for (int j = 0; j < k; ++j) o += block(j); return o; }
  std::size_t nparam() const { return offset(K); }
};

// Restricted cubic spline basis in x = log t and its derivative with respect to x:
// Royston & Parmar (2002), Sec. 2.2: s = g0 + g1 x + sum_j g_(j+1) v_j(x),
// v_j = (x-k_j)+^3 - l_j (x-kmin)+^3 - (1-l_j)(x-kmax)+^3, l_j = (kmax-k_j)/(kmax-kmin).
// With knots, the basis is evaluated on the rescaled variable (x - kmin)/(kmax - kmin), so
// that all basis functions are of order one. This spans exactly the same functions s(x) as
// the formula above (an affine change of variable maps one basis onto the other); only the
// numerical values of gamma differ, and the information matrix is far better conditioned.
// (Royston & Parmar's software orthogonalises the basis for the same reason.)
inline void rcs_basis(const std::vector<double>& kn, double x, double* z, double* dz) {
  z[0] = 1.0; dz[0] = 0.0;
  const std::size_t M = kn.size();
  if (M < 3) { z[1] = x; dz[1] = 1.0; return; }
  const double kmin = kn[0], R = kn[M - 1] - kn[0], xs = (x - kmin) / R;
  z[1] = xs; dz[1] = 1.0 / R;
  auto c3 = [](double a) { return a > 0.0 ? a * a * a : 0.0; };
  auto c2 = [](double a) { return a > 0.0 ? 3.0 * a * a : 0.0; };
  for (std::size_t j = 1; j + 1 < M; ++j) {
    const double kj = (kn[j] - kmin) / R, l = 1.0 - kj;       // scaled knots: kmin -> 0, kmax -> 1
    z[j + 1] = c3(xs - kj) - l * c3(xs) - (1.0 - l) * c3(xs - 1.0);
    dz[j + 1] = (c2(xs - kj) - l * c2(xs) - (1.0 - l) * c2(xs - 1.0)) / R;
  }
}

// Baseline of form (A)/(B) at time t (t may be +infinity). q[0..2] are the local baseline
// variables (for SPLINE: q[0] = s(log t), q[1] = ds/dlog t). Returns 0 = ok, 1 = outside
// the support, 2 = u is infinite (F = 1). logdu (log du/dt) is filled only for finite t.
template <class T>
int base_u(int type, double t, const T* q, const T& alpha, T& u, T* logdu) {
  const bool inf = std::isinf(t);
  const double lt = inf ? 0.0 : std::log(t);
  switch (type) {
    case BASE_GOMPERTZ: {
      if (inf) { if (val(q[0]) >= 0.0) return 2; u = t_exp(q[1]) * t_recip_neg(q[0]); return 0; }
      const T w = t * q[0];
      u = t * (t_exp(q[1]) * t_phi(w));
      if (logdu) *logdu = q[1] + w;
      return 0;
    }
    case BASE_GOMPERTZ3: {                       // u = tau e^zeta t phi(rho t) phi(zeta m), m = e^(rho t) - 1
      if (inf) {
        const double r = val(q[0]), ze = val(q[2]);
        if (r < 0.0) { u = t_exp(q[1] + q[2]) * t_phi((-1.0) * q[2]) * t_recip_neg(q[0]); return 0; }
        if (r > 0.0 && ze < 0.0) { u = t_exp(q[1] + q[2]) * (t_recip_neg(q[0]) * t_recip_neg(q[2])); return 0; }
        return 2;
      }
      const T w = t * q[0], m = t_expm1(w);
      u = t * (t_exp(q[1] + q[2]) * t_phi(w) * t_phi(q[2] * m));
      if (logdu) *logdu = q[1] + w + q[2] * t_exp(w);
      if (!std::isfinite(val(u))) return 1;
      return 0;
    }
    case BASE_WEIBULL: {                         // q = (log lambda, log kappa)
      if (inf) return 2;
      const T lw = t_exp(q[1]) * (q[0] + lt);
      u = t_exp(lw);
      if (logdu) *logdu = q[1] + lw + (-lt);
      return 0;
    }
    case BASE_SPLINE: {
      if (inf) return 2;
      if (!(val(q[1]) > 0.0)) return 1;          // the log cumulative hazard must increase
      u = t_exp(q[0]);
      if (logdu) *logdu = q[0] + t_log(q[1]) + (-lt);
      return 0;
    }
    default: break;
  }
  // form (B): L = -log(1 - F0), log L' = log f0 + L; then u = L phi(alpha L), log u' = log L' + alpha L
  T L = T(0.0), logdL = T(0.0);
  if (type == BASE_LOGISTIC3) {                  // q = (log b, c, logit p); F0 = p sigma(b(t-c)) (1 - e^(-bt))
    if (inf) { L = t_softplus(q[2]); }
    else {
      const T b = t_exp(q[0]), x = b * ((-1.0) * q[1] + t), bt = t * b;
      const T F0 = t_sigmoid(q[2]) * t_sigmoid(x) * ((-1.0) * t_expm1((-1.0) * bt));
      if (!(val(F0) < 1.0)) return 1;
      L = t_neglog1m(F0);
      // log f0 = log p + log b + x + log(1 + e^(-bc)) - 2 log(1 + e^x)
      logdL = (-1.0) * t_softplus((-1.0) * q[2]) + q[0] + x + t_softplus((-1.0) * (b * q[1])) + (-2.0) * t_softplus(x) + L;
    }
  } else {                                       // mixtures: q = (logit p, log lambda | mu, log kappa | log sigma)
    if (inf) { L = t_softplus(q[0]); }
    else {
      T G, logg;
      if (type == BASE_MIX_LOGNORMAL) {
        const T z = t_exp((-1.0) * q[2]) * ((-1.0) * q[1] + lt);
        G = t_normcdf(z);
        logg = (-0.5) * (z * z) + (-1.0) * q[2] + (-0.9189385332046727 - lt);
      } else {
        const T lw = t_exp(q[2]) * (q[1] + lt);
        if (type == BASE_MIX_WEIBULL) {
          const T w = t_exp(lw);
          G = (-1.0) * t_expm1((-1.0) * w);
          logg = q[2] + lw + (-1.0) * w + (-lt);
        } else {
          G = t_sigmoid(lw);
          logg = q[2] + lw + (-2.0) * t_softplus(lw) + (-lt);
        }
      }
      const T F0 = t_sigmoid(q[0]) * G;
      if (!(val(F0) < 1.0)) return 1;
      L = t_neglog1m(F0);
      logdL = (-1.0) * t_softplus((-1.0) * q[0]) + logg + L;
    }
  }
  const T aL = alpha * L;
  u = L * t_phi(aL);
  if (logdu) *logdu = logdL + aL;
  if (!std::isfinite(val(u))) return 1;
  return 0;
}

// F and log f from the baseline, alpha and eta. Returns false outside the support.
// At t = infinity only F is produced; F = 1 when the baseline is proper.
template <class T>
bool direct_Ff(int type, double t, const T* q, const T& alpha, const T& eta, T* F, T* logf) {
  T u, logdu;
  const int rc = base_u(type, t, q, alpha, u, std::isinf(t) ? nullptr : &logdu);
  if (rc == 1) return false;
  if (rc == 2) { if (F) *F = T(1.0); return true; }
  const T eu = t_exp(eta) * u, v = alpha * eu;
  if (!(val(v) > -1.0) || !std::isfinite(val(v))) return false;
  const T A = eu * t_psi(v);
  if (F) { const double e = std::exp(-val(A)); *F = unary(A, -std::expm1(-val(A)), e, -e); }
  if (logf && !std::isinf(t)) *logf = eta + logdu + (-1.0) * A + (-1.0) * unary(v, std::log1p(val(v)), 1.0 / (1.0 + val(v)), -1.0 / ((1.0 + val(v)) * (1.0 + val(v))));
  return true;
}

class DirectCIFModel {
 public:
  // cause: 0 = censored, 1..K. X row-major n x p (may be null when p = 0).
  DirectCIFModel(std::size_t n, const double* time, const int* cause, const double* X, const DirectSpec& spec)
      : n_(n), time_(time), cause_(cause), X_(X), s_(spec), P_(spec.nparam()) {
    for (int k = 0; k < s_.K; ++k) { off_.push_back(s_.offset(k)); nb_.push_back(s_.nbase(k)); }
    tmax_ = 0.0;
    for (std::size_t i = 0; i < n; ++i) tmax_ = std::fmax(tmax_, time[i]);
    zrow_.resize(s_.K); dzrow_.resize(s_.K);
    for (int k = 0; k < s_.K; ++k) {
      if (s_.baseline[k] != BASE_SPLINE) continue;
      zrow_[k].resize(n * nb_[k]); dzrow_[k].resize(n * nb_[k]);
      for (std::size_t i = 0; i < n; ++i) rcs_basis(s_.knots[k], std::log(time[i]), &zrow_[k][i * nb_[k]], &dzrow_[k][i * nb_[k]]);
    }
  }
  const DirectSpec& spec() const { return s_; }
  std::size_t nparam() const { return P_; }
  double tmax() const { return tmax_; }

  // CIF of cause k (0-based) for covariate row x at time t (t may be +infinity). 1 outside the support.
  double cif(const std::vector<double>& th, int k, const double* x, double t) const {
    if (t <= 0.0) return 0.0;
    double q[3], F = 1.0;
    std::vector<double> z, dz;
    locals<double>(th, k, t, nullptr, nullptr, q, z, dz);
    const double al = alpha_of(th, k), et = eta_of(th, k, x);
    if (!direct_Ff<double>(s_.baseline[k], t, q, al, et, &F, nullptr)) return 1.0;
    return F;
  }
  double total(const std::vector<double>& th, std::size_t i, double t) const {
    double s = 0.0;
    for (int k = 0; k < s_.K; ++k) s += cif(th, k, X_ + i * s_.p, t);
    return s;
  }
  bool feasible_in_range(const std::vector<double>& th) const {
    for (std::size_t i = 0; i < n_; ++i) if (!(total(th, i, tmax_) < 1.0)) return false;
    return true;
  }

  // Log-likelihood, gradient and observed information (out.H = minus the Hessian).
  // out.ll = -infinity when the parameter value is outside the support of the model.
  void eval(const std::vector<double>& th, Eval& out, bool need_hess) const {
    const std::size_t P = P_;
    const int K = s_.K;
    out.ll = 0.0;
    out.g.assign(P, 0.0);
    if (out.H.n != P) out.H = Mat(P);
    out.H.zero();
    std::vector<D2> F(K);
    std::vector<double> gs(P), z, dz;
    std::vector<Term> terms;
    D2 lf, q[3];
    for (std::size_t i = 0; i < n_; ++i) {
      const double* x = X_ + i * s_.p;
      const double t = time_[i];
      const int c = cause_[i];
      if (c >= 1) {
        const int k = c - 1;
        const double* zr = zrow_[k].empty() ? nullptr : &zrow_[k][i * nb_[k]];
        const double* dr = dzrow_[k].empty() ? nullptr : &dzrow_[k][i * nb_[k]];
        locals<D2>(th, k, t, zr, dr, q, z, dz);
        if (!direct_Ff<D2>(s_.baseline[k], t, q, D2::variable(3, alpha_of(th, k)), D2::variable(4, eta_of(th, k, x)), nullptr, &lf) ||
            !std::isfinite(lf.v)) { bad(out); return; }
        out.ll += lf.v;
        make_terms(k, x, zr ? zr : z.data(), dr ? dr : dz.data(), terms);
        add_terms(terms, lf, 1.0, out, need_hess, nullptr);
      } else {
        double S = 1.0;
        std::fill(gs.begin(), gs.end(), 0.0);
        bool ok = true;
        for (int k = 0; k < K && ok; ++k) {
          const double* zr = zrow_[k].empty() ? nullptr : &zrow_[k][i * nb_[k]];
          const double* dr = dzrow_[k].empty() ? nullptr : &dzrow_[k][i * nb_[k]];
          locals<D2>(th, k, t, zr, dr, q, z, dz);
          ok = direct_Ff<D2>(s_.baseline[k], t, q, D2::variable(3, alpha_of(th, k)), D2::variable(4, eta_of(th, k, x)), &F[k], nullptr);
          S -= F[k].v;
        }
        if (!ok || !(S > 0.0)) { bad(out); return; }
        out.ll += std::log(S);
        for (int k = 0; k < K; ++k) {
          const double* zr = zrow_[k].empty() ? nullptr : &zrow_[k][i * nb_[k]];
          const double* dr = dzrow_[k].empty() ? nullptr : &dzrow_[k][i * nb_[k]];
          make_terms(k, x, zr ? zr : z.data(), dr ? dr : dz.data(), terms);
          add_terms(terms, F[k], -1.0 / S, out, need_hess, gs.data());
        }
        if (need_hess)                                 // -(grad sum F)(grad sum F)^T / S^2
          for (std::size_t a = 0; a < P; ++a) {
            if (gs[a] == 0.0) continue;
            const double ga = gs[a] / (S * S);
            for (std::size_t b = 0; b <= a; ++b) out.H(a, b) -= ga * gs[b];
          }
      }
    }
    if (!std::isfinite(out.ll)) { bad(out); return; }
    if (need_hess)
      for (std::size_t a = 0; a < P; ++a)
        for (std::size_t b = 0; b <= a; ++b) { out.H(a, b) = -out.H(a, b); out.H(b, a) = out.H(a, b); }
  }

  // Gradient of F_k(t; x) with respect to all parameters (delta method). Returns F.
  double cif_grad(const std::vector<double>& th, int k, const double* x, double t, std::vector<double>& grad) const {
    grad.assign(P_, 0.0);
    if (t <= 0.0) return 0.0;
    D2 q[3], F;
    std::vector<double> z, dz;
    locals<D2>(th, k, t, nullptr, nullptr, q, z, dz);
    if (!direct_Ff<D2>(s_.baseline[k], t, q, D2::variable(3, alpha_of(th, k)), D2::variable(4, eta_of(th, k, x)), &F, nullptr)) return 1.0;
    std::vector<Term> terms;
    make_terms(k, x, z.data(), dz.data(), terms);
    for (const Term& e : terms) grad[e.idx] += F.g[e.loc] * e.coef;
    return F.v;
  }

 protected:
  struct Term { std::size_t idx; int loc; double coef; };   // local variable `loc` depends on parameter idx with weight coef
  static void bad(Eval& out) { out.ll = -std::numeric_limits<double>::infinity(); }

  double alpha_of(const std::vector<double>& th, int k) const { return s_.alpha_free[k] ? th[off_[k] + nb_[k]] : s_.alpha_value[k]; }
  double eta_of(const std::vector<double>& th, int k, const double* x) const {
    const double* b = &th[off_[k] + nb_[k] + (s_.alpha_free[k] ? 1 : 0)];
    double e = 0.0;
    for (std::size_t j = 0; j < s_.p; ++j) e += x[j] * b[j];
    return e;
  }
  // Local baseline variables of cause k at time t. For a spline, z and dz receive the basis
  // rows (taken from zr / dr when precomputed), and q[0] = s, q[1] = ds/dlog t.
  template <class T>
  void locals(const std::vector<double>& th, int k, double t, const double* zr, const double* dr, T* q,
              std::vector<double>& z, std::vector<double>& dz) const {
    const std::size_t o = off_[k], nb = nb_[k];
    if (s_.baseline[k] == BASE_SPLINE) {
      if (!zr) {
        z.assign(nb, 0.0); dz.assign(nb, 0.0);
        if (std::isfinite(t)) rcs_basis(s_.knots[k], std::log(t), z.data(), dz.data());
        zr = z.data(); dr = dz.data();
      }
      double sv = 0.0, dv = 0.0;
      for (std::size_t j = 0; j < nb; ++j) { sv += th[o + j] * zr[j]; dv += th[o + j] * dr[j]; }
      q[0] = MakeVar<T>::make(0, sv); q[1] = MakeVar<T>::make(1, dv); q[2] = T(0.0);
    } else {
      for (std::size_t j = 0; j < 3; ++j) q[j] = j < nb ? MakeVar<T>::make(static_cast<int>(j), th[o + j]) : T(0.0);
    }
  }
  void make_terms(int k, const double* x, const double* z, const double* dz, std::vector<Term>& terms) const {
    terms.clear();
    const std::size_t o = off_[k], nb = nb_[k];
    if (s_.baseline[k] == BASE_SPLINE) {
      for (std::size_t j = 0; j < nb; ++j) { terms.push_back({o + j, 0, z[j]}); if (dz[j] != 0.0) terms.push_back({o + j, 1, dz[j]}); }
    } else {
      for (std::size_t j = 0; j < nb; ++j) terms.push_back({o + j, static_cast<int>(j), 1.0});
    }
    std::size_t q = o + nb;
    if (s_.alpha_free[k]) terms.push_back({q++, 3, 1.0});
    for (std::size_t j = 0; j < s_.p; ++j) terms.push_back({q + j, 4, x[j]});
  }
  // Adds w * (gradient, Hessian) of the local function d to the totals (lower triangle of H).
  // gsum (optional) accumulates the unweighted parameter gradient.
  void add_terms(const std::vector<Term>& terms, const D2& d, double w, Eval& out, bool need_hess, double* gsum) const {
    for (const Term& e : terms) {
      const double g = d.g[e.loc] * e.coef;
      out.g[e.idx] += w * g;
      if (gsum) gsum[e.idx] += g;
    }
    if (!need_hess) return;
    for (const Term& e : terms) {
      const double we = w * e.coef;
      for (const Term& f : terms) if (f.idx <= e.idx) out.H(e.idx, f.idx) += we * d.h[e.loc][f.loc] * f.coef;
    }
  }

  std::size_t n_;
  const double* time_;
  const int* cause_;
  const double* X_;
  DirectSpec s_;
  std::size_t P_;
  std::vector<std::size_t> off_, nb_;
  std::vector<std::vector<double>> zrow_, dzrow_;
  double tmax_;
};

// ---------------------------------------------------------------------------
// Default knots of a spline baseline: boundary knots at the smallest and largest
// uncensored log time of the cause, interior knots at equally spaced centiles of
// the uncensored log times (Royston & Parmar 2002, Sec. 2.2 and Table I; Mozumder
// et al. 2018, Sec. 2.5). df = 1 gives no knots (Weibull / log-logistic).
// ---------------------------------------------------------------------------
inline std::vector<double> direct_default_knots(std::size_t n, const double* time, const int* cause, int k1, int df) {
  std::vector<double> x;
  for (std::size_t i = 0; i < n; ++i) if (cause[i] == k1 && time[i] > 0.0) x.push_back(std::log(time[i]));
  std::sort(x.begin(), x.end());
  std::vector<double> kn;
  if (df < 2 || x.size() < 2 || !(x.back() > x.front())) return kn;
  for (int j = 0; j <= df; ++j) {
    const double pos = static_cast<double>(j) / df * static_cast<double>(x.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(pos);
    const double fr = pos - static_cast<double>(lo);
    const double v = lo + 1 < x.size() ? x[lo] + fr * (x[lo + 1] - x[lo]) : x.back();
    if (kn.empty() || v > kn.back()) kn.push_back(v);
  }
  if (kn.size() < 3) kn.clear();
  return kn;
}

// ---------------------------------------------------------------------------
// Fitting: Levenberg-Marquardt damped Newton on the exact Hessian, in stages
// (baseline only; + regression coefficients; + alpha), started from the
// Aalen-Johansen estimate of each cause. Marquardt (1963), SIAM J 11:431-441;
// Nocedal & Wright (2006), Numerical Optimization, Sec. 10.3.
// ---------------------------------------------------------------------------
enum DirectFlag : int {
  DIRECT_SUM_EXCEEDS_IN_RANGE = 16,    // some subject has sum_k F_k(tmax) > 1
  DIRECT_SUM_EXCEEDS_AT_INFINITY = 32, // some subject has sum_k F_k(infinity) > 1
  DIRECT_CONSTRAINT_ACTIVE = 64        // constrained fit stopped on the boundary
};

struct DirectDiagnostics {
  double n_exceed_range = 0, max_sum_range = 0, n_exceed_inf = 0, max_sum_inf = 0;
  double first_crossing_time = std::numeric_limits<double>::quiet_NaN();
};

template <class Model>
bool lm_stage(const Model& m, std::vector<double>& th, const std::vector<char>& active, int max_iter,
              double tol, bool constrain, Eval& ev, int& iters, bool& on_boundary, bool& stalled) {
  const std::size_t P = th.size();
  std::vector<std::size_t> idx;
  for (std::size_t j = 0; j < P; ++j) if (active[j]) idx.push_back(j);
  const std::size_t q = idx.size();
  m.eval(th, ev, true);
  if (!std::isfinite(ev.ll)) return false;
  double lambda = 1e-3;
  std::vector<double> step(q), trial(P);
  Eval et;
  on_boundary = false; stalled = false;
  int flat = 0;                        // consecutive iterations that move the parameters but not the likelihood
  for (int it = 0; it < max_iter; ++it) {
    ++iters;
    bool accepted = false, blocked_by_constraint = false;
    double pred = 0.0, smax = 0.0;
    for (int tries = 0; tries < 60 && !accepted; ++tries) {
      // Marquardt (1963) scaling: the damped system is solved for D H D with D = diag(|H_aa|)^(-1/2),
      // so that neither the damping nor the positive-definiteness test depends on the units of
      // time or of the covariates.
      Mat M(q);
      std::vector<double> dsc(q);
      for (std::size_t a = 0; a < q; ++a) {
        const double haa = std::fabs(ev.H(idx[a], idx[a]));
        dsc[a] = 1.0 / std::sqrt(haa > 0.0 && std::isfinite(haa) ? haa : 1.0);
      }
      for (std::size_t a = 0; a < q; ++a) {
        for (std::size_t b = 0; b < q; ++b) M(a, b) = ev.H(idx[a], idx[b]) * dsc[a] * dsc[b];
        M(a, a) += lambda;
      }
      if (!cholesky(M)) { lambda *= 10.0; if (lambda > 1e12) break; continue; }
      for (std::size_t a = 0; a < q; ++a) step[a] = ev.g[idx[a]] * dsc[a];
      chol_solve(M, step);
      for (std::size_t a = 0; a < q; ++a) step[a] *= dsc[a];
      trial = th;
      pred = 0.0; smax = 0.0;
      for (std::size_t a = 0; a < q; ++a) {
        trial[idx[a]] += step[a];
        pred += step[a] * ev.g[idx[a]];
        smax = std::fmax(smax, std::fabs(step[a]) / (1.0 + std::fabs(th[idx[a]])));
      }
      if (constrain && !m.feasible_in_range(trial)) { blocked_by_constraint = true; lambda *= 10.0; continue; }
      m.eval(trial, et, true);
      if (std::isfinite(et.ll) && et.ll >= ev.ll - 1e-12 * (std::fabs(ev.ll) + 1.0)) accepted = true;
      else lambda *= 10.0;
      if (lambda > 1e12) break;
    }
    if (!accepted) { on_boundary = blocked_by_constraint; return blocked_by_constraint; }
    const double gain = et.ll - ev.ll;
    th = trial;
    std::swap(ev, et);
    lambda = std::fmax(lambda / 10.0, 1e-12);
    if (gain <= tol * (std::fabs(ev.ll) + 1.0) && smax <= 1e-6 && lambda <= 1e-6) return true;
    if (0.5 * pred <= 1e-13 * (std::fabs(ev.ll) + 1.0) && smax <= 1e-8) return true;
    // A flat ridge: the likelihood has stopped increasing while a parameter keeps drifting
    // (it is not identified by the data). Stop and report instead of iterating to the limit.
    // An exploding score: the iterates approach a singularity of the likelihood (with alpha < -1 the
    // density is unbounded where 1 + alpha e^eta u(t) reaches 0 at an event time). No maximum exists.
    double gmax = 0.0;
    for (std::size_t a = 0; a < q; ++a) gmax = std::fmax(gmax, std::fabs(ev.g[idx[a]]));
    if (!(gmax <= 1e10 * (std::fabs(ev.ll) + 1.0))) { stalled = true; return false; }
    flat = gain <= 1e-9 * (std::fabs(ev.ll) + 1.0) ? flat + 1 : 0;
    if (flat >= 15) { stalled = true; return false; }
  }
  return false;
}

struct DirectFit {
  FitResult fit;                       // beta = full parameter vector
  DirectDiagnostics diag;
};

namespace detail {
// least squares y ~ a + b x; returns false when x has no spread
inline bool ols(const std::vector<double>& x, const std::vector<double>& y, double& a, double& b) {
  const std::size_t m = x.size();
  if (m < 2) return false;
  double mx = 0, my = 0, sxx = 0, sxy = 0;
  for (std::size_t i = 0; i < m; ++i) { mx += x[i]; my += y[i]; }
  mx /= static_cast<double>(m); my /= static_cast<double>(m);
  for (std::size_t i = 0; i < m; ++i) { sxx += (x[i] - mx) * (x[i] - mx); sxy += (x[i] - mx) * (y[i] - my); }
  if (!(sxx > 0.0)) return false;
  b = sxy / sxx; a = my - b * mx;
  return std::isfinite(a) && std::isfinite(b);
}
inline double qnorm(double p) {        // inverse normal cdf by bisection (start values only)
  double lo = -10.0, hi = 10.0;
  for (int i = 0; i < 80; ++i) { const double mid = 0.5 * (lo + hi); if (0.5 * std::erfc(-mid * 0.7071067811865476) < p) lo = mid; else hi = mid; }
  return 0.5 * (lo + hi);
}
}  // namespace detail

// Index (within the baseline block) of the parameter that is lowered to shrink the CIF.
inline std::size_t base_level_index(int type) {
  switch (type) {
    case BASE_GOMPERTZ: case BASE_GOMPERTZ3: return 1;       // log tau
    case BASE_LOGISTIC3: return 2;                           // logit p
    default: return 0;                                       // log lambda, logit p, gamma_0
  }
}

// Start values: each baseline is matched to the Aalen-Johansen estimate of its cause.
// The least-squares starts for the Weibull-type forms follow Cheng (2009), Sec. 2.4.
inline std::vector<double> direct_start(std::size_t n, const double* time, const int* cause, const DirectSpec& s) {
  std::vector<double> t(time, time + n);
  std::vector<int> c(cause, cause + n);
  const AJResult aj = aalen_johansen(t, c, s.K);
  const std::size_t m = aj.time.size();
  std::vector<double> th(s.nparam(), 0.0);
  for (int k = 0; k < s.K; ++k) {
    const std::size_t o = s.offset(k);
    const int type = s.baseline[k];
    std::vector<double> at, af;                               // points where the CIF of cause k jumps
    double prev = 0.0;
    for (std::size_t j = 0; j < m; ++j) { const double F = aj.cif[j * s.K + k]; if (F > prev && aj.time[j] > 0.0) { at.push_back(aj.time[j]); af.push_back(F); prev = F; } }
    const double pmax = af.empty() ? 1e-3 : af.back();
    double pk = std::fmin(std::fmax(pmax, 1e-4), 0.98);
    const double q0 = std::fmin(pk * 1.05, 0.99);             // plateau start, slightly above the last estimate
    double th_time = at.empty() ? 1.0 : at.back();            // first time the CIF reaches half its final value
    for (std::size_t j = 0; j < at.size(); ++j) if (af[j] >= 0.5 * pk) { th_time = at[j]; break; }
    std::vector<double> lx, ly;
    double a = 0.0, b = 1.0;
    if (type == BASE_GOMPERTZ || type == BASE_GOMPERTZ3) {
      const double cc = -std::log1p(-pk), r = std::log1p(-0.5 * pk) / std::log1p(-pk);
      const double rho = std::log1p(-r) / th_time;
      th[o] = rho; th[o + 1] = std::log(-rho * cc);
    } else if (type == BASE_WEIBULL) {                        // log(-log(1 - F)) = kappa (log lambda + log t)
      for (std::size_t j = 0; j < at.size(); ++j) { lx.push_back(std::log(at[j])); ly.push_back(std::log(-std::log1p(-std::fmin(af[j], 0.999)))); }
      if (!detail::ols(lx, ly, a, b) || !(b > 0.05)) { b = 1.0; a = std::log(-std::log1p(-0.5 * pk)) - std::log(th_time); }
      th[o] = a / b; th[o + 1] = std::log(b);
    } else if (type == BASE_SPLINE) {
      const std::size_t nb = s.nbase(k);
      for (std::size_t j = 0; j < at.size(); ++j) { lx.push_back(std::log(at[j])); ly.push_back(std::log(-std::log1p(-std::fmin(af[j], 0.999)))); }
      if (!detail::ols(lx, ly, a, b) || !(b > 0.05)) { b = 1.0; a = std::log(-std::log1p(-0.5 * pk)) - std::log(th_time); }
      th[o] = a; th[o + 1] = b;                               // linear start; refined below when it stays monotone
      if (nb > 2 && lx.size() > 2 * nb) {
        Mat A(nb);
        std::vector<double> rhs(nb, 0.0), z(nb), dz(nb);
        for (std::size_t j = 0; j < lx.size(); ++j) {
          rcs_basis(s.knots[k], lx[j], z.data(), dz.data());
          for (std::size_t u = 0; u < nb; ++u) { rhs[u] += z[u] * ly[j]; for (std::size_t v = 0; v < nb; ++v) A(u, v) += z[u] * z[v]; }
        }
        for (std::size_t u = 0; u < nb; ++u) A(u, u) *= 1.0 + 1e-8;
        if (cholesky(A)) {
          chol_solve(A, rhs);
          bool mono = true;
          for (std::size_t i = 0; i < n && mono; ++i) {
            if (!(time[i] > 0.0)) continue;
            rcs_basis(s.knots[k], std::log(time[i]), z.data(), dz.data());
            double d = 0.0;
            for (std::size_t u = 0; u < nb; ++u) d += rhs[u] * dz[u];
            mono = d > 1e-3;
          }
          if (mono) for (std::size_t u = 0; u < nb; ++u) th[o + u] = rhs[u];
        }
      }
    } else if (type == BASE_LOGISTIC3) {                      // logit(F / q0) ~ b (t - c)
      for (std::size_t j = 0; j < at.size(); ++j) { const double r = std::fmin(std::fmax(af[j] / q0, 1e-4), 0.999); lx.push_back(at[j]); ly.push_back(std::log(r / (1.0 - r))); }
      if (!detail::ols(lx, ly, a, b) || !(b > 0.0)) { b = 2.0 / std::fmax(th_time, 1e-8); a = -b * th_time; }
      th[o] = std::log(b); th[o + 1] = -a / b; th[o + 2] = std::log(q0 / (1.0 - q0));
    } else {                                                  // mixtures: transform of F / q0 linear in log t
      for (std::size_t j = 0; j < at.size(); ++j) {
        const double r = std::fmin(std::fmax(af[j] / q0, 1e-4), 0.999);
        lx.push_back(std::log(at[j]));
        ly.push_back(type == BASE_MIX_WEIBULL ? std::log(-std::log1p(-r)) : (type == BASE_MIX_LOGLOGISTIC ? std::log(r / (1.0 - r)) : detail::qnorm(r)));
      }
      if (!detail::ols(lx, ly, a, b) || !(b > 0.05)) { b = 1.0; a = -std::log(th_time); }
      th[o] = std::log(q0 / (1.0 - q0));
      if (type == BASE_MIX_LOGNORMAL) { th[o + 1] = -a / b; th[o + 2] = std::log(1.0 / b); }   // z = (log t - mu)/sigma
      else { th[o + 1] = a / b; th[o + 2] = std::log(b); }                                       // kappa (log lambda + log t)
    }
    if (s.alpha_free[k]) th[o + s.nbase(k)] = s.alpha_value[k];
  }
  return th;
}

inline DirectDiagnostics direct_diagnose(const DirectCIFModel& m, std::size_t n, const std::vector<double>& th) {
  DirectDiagnostics d;
  const double inf = std::numeric_limits<double>::infinity(), tmax = m.tmax();
  double first = inf;
  for (std::size_t i = 0; i < n; ++i) {
    const double sr = m.total(th, i, tmax), si = m.total(th, i, inf);
    d.max_sum_range = std::fmax(d.max_sum_range, sr);
    d.max_sum_inf = std::fmax(d.max_sum_inf, si);
    if (si > 1.0) d.n_exceed_inf += 1;
    if (sr > 1.0) {
      d.n_exceed_range += 1;
      double lo = 0.0, hi = tmax;                      // the sum is increasing in t: bisection for the crossing
      for (int b = 0; b < 60; ++b) { const double mid = 0.5 * (lo + hi); if (m.total(th, i, mid) > 1.0) hi = mid; else lo = mid; }
      first = std::fmin(first, hi);
    }
  }
  if (first < inf) d.first_crossing_time = first;
  return d;
}

// Staged fit shared by all direct models. role[j]: 0 = baseline parameter, 1 = regression
// coefficient, 2 = alpha. Stage 1 fits the baseline parameters, stage 2 adds the regression
// coefficients, the last stage all parameters. shrink(th) moves a start value with a non-finite
// likelihood towards smaller CIFs.
template <class Model, class Shrink>
FitResult staged_fit(const Model& m, std::vector<double> th, const std::vector<int>& role, bool have_init,
                     bool constrain, int max_iter, double tol, Shrink shrink) {
  const std::size_t P = th.size();
  FitResult r;
  Eval ev;
  m.eval(th, ev, false);
  for (int it = 0; it < 80 && (!std::isfinite(ev.ll) || (constrain && !m.feasible_in_range(th))); ++it) {
    shrink(th);
    m.eval(th, ev, false);
  }
  r.ll_init = ev.ll;
  r.beta = th;
  if (!std::isfinite(ev.ll)) { r.flags |= FIT_NONFINITE; return r; }
  std::vector<char> act(P, 0);
  bool ok = true, boundary = false, stalled = false;
  const bool any_beta = std::any_of(role.begin(), role.end(), [](int a) { return a == 1; });
  const bool any_alpha = std::any_of(role.begin(), role.end(), [](int a) { return a == 2; });
  const int stages = have_init ? 1 : (1 + (any_beta ? 1 : 0) + (any_alpha ? 1 : 0));
  for (int st = 0; st < stages; ++st) {
    const bool last = st == stages - 1;
    for (std::size_t j = 0; j < P; ++j) act[j] = last || role[j] == 0 || (role[j] == 1 && st >= 1);
    ok = lm_stage(m, th, act, max_iter, last ? tol : std::fmax(tol, 1e-8), constrain, ev, r.iterations, boundary, stalled);
  }
  r.beta = th;
  r.ll = ev.ll;
  r.converged = ok && !boundary;
  if (!ok) r.flags |= stalled ? FIT_POSSIBLE_MONOTONE : FIT_MAX_ITER;   // 8: a parameter is drifting on a flat likelihood
  if (boundary) r.flags |= DIRECT_CONSTRAINT_ACTIVE;
  for (double gj : ev.g) r.grad_maxabs = std::fmax(r.grad_maxabs, std::fabs(gj));
  {                                    // variance = H^-1, factorised on the correlation scale (unit-free PD test)
    Mat L = ev.H;
    std::vector<double> dsc(P, 1.0);
    bool okd = true;
    for (std::size_t a = 0; a < P; ++a) { const double h = ev.H(a, a); if (h > 0.0 && std::isfinite(h)) dsc[a] = 1.0 / std::sqrt(h); else okd = false; }
    for (std::size_t a = 0; a < P; ++a) for (std::size_t b = 0; b < P; ++b) L(a, b) = ev.H(a, b) * dsc[a] * dsc[b];
    if (okd && cholesky(L)) {
      r.var = chol_inverse(L);
      for (std::size_t a = 0; a < P; ++a) for (std::size_t b = 0; b < P; ++b) r.var(a, b) *= dsc[a] * dsc[b];
    } else r.flags |= FIT_INFO_NOT_PD;
  }
  if (r.flags & FIT_INFO_NOT_PD) r.converged = false;   // a maximum needs a positive definite information matrix
  return r;
}

// role of every parameter of a DirectSpec (see staged_fit)
inline std::vector<int> direct_roles(const DirectSpec& s) {
  std::vector<int> role(s.nparam(), 0);
  for (int k = 0; k < s.K; ++k) {
    const std::size_t o = s.offset(k), nb = s.nbase(k), na = s.alpha_free[k] ? 1 : 0;
    if (na) role[o + nb] = 2;
    for (std::size_t j = 0; j < s.p; ++j) role[o + nb + na + j] = 1;
  }
  return role;
}

// constrain = false: the models exactly as published; the sum of the fitted CIFs is only
// reported. constrain = true: steps leaving the region where sum_k F_k(tmax; z_i) < 1 for
// every subject are rejected.
inline DirectFit direct_fit(const DirectCIFModel& m, std::size_t n, const double* time, const int* cause,
                            const std::vector<double>* init, bool constrain, int max_iter, double tol) {
  const DirectSpec& s = m.spec();
  DirectFit out;
  out.fit = staged_fit(m, init ? *init : direct_start(n, time, cause, s), direct_roles(s), init != nullptr, constrain, max_iter, tol,
                       [&s](std::vector<double>& th) { for (int k = 0; k < s.K; ++k) th[s.offset(k) + base_level_index(s.baseline[k])] -= 0.5; });
  if (out.fit.flags & FIT_NONFINITE) return out;
  out.diag = direct_diagnose(m, n, out.fit.beta);
  if (out.diag.n_exceed_range > 0) out.fit.flags |= DIRECT_SUM_EXCEEDS_IN_RANGE;
  if (out.diag.n_exceed_inf > 0) out.fit.flags |= DIRECT_SUM_EXCEEDS_AT_INFINITY;
  return out;
}

}  // namespace survcore
