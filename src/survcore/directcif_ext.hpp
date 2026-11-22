// survcore: two further ways of fitting the direct parametric model for a cumulative incidence
// function (same parametrisation of F_1(t; z) as directcif.hpp):
//
//  SubdistCIFModel      the cause of interest alone, by the inverse-probability-of-censoring
//                       weighted likelihood of the subdistribution ("parametric Fine-Gray").
//  ConstrainedCIFModel  two causes with asymptotes that add up to one,
//                       F_2(t; z) = {1 - F_1(infinity; z)} G(t): Shi, Cheng & Jeong (2013), eq. (7).
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "directcif.hpp"

namespace survcore {

// ---------------------------------------------------------------------------
// 1. Weighted likelihood for the cause of interest
//
// With complete censoring information (the censoring time C_i of every subject known) the
// subdistribution of cause 1 has the ordinary censored-data likelihood (Fine & Gray 1999,
// Sec. 3.1): a subject failing from cause 1 at T contributes f_1(T; z), every other subject
// 1 - F_1(C; z). When C_i is unknown for the subjects who failed from a competing cause at T_i,
// it is replaced by its estimated conditional distribution given C_i >= T_i, obtained from the
// Kaplan-Meier estimate G of the censoring distribution (Fine & Gray 1999, Sec. 3.2; Geskus
// 2011, Biometrics 67:39-49; parametric use: Lambert, Wilkes & Crowther 2017, Stat Med
// 36:1429-1446):
//
//   l(theta) = sum_{cause 1} log f_1(T_i; z_i) + sum_{censored} log{1 - F_1(T_i; z_i)}
//            + sum_{competing} sum_{c_j >= T_i} w_ij log{1 - F_1(c_j; z_i)},
//   w_ij = {G(c_j-) - G(c_j)} / G(T_i-),
//
// where c_j runs over the distinct censoring times and the mass G leaves beyond the last
// censoring time is placed at the largest observed time. This is the time-integrated form of
// the weighted likelihood with weights G(t)/G(T_i) (integration by parts of the weighted
// cumulative subdistribution hazard). Ties: failures precede censorings.
//
// The variance is the sandwich H^-1 (sum_i s_i s_i') H^-1 with one score s_i per subject; the
// sampling variability of G is ignored, as in the cited weighted implementations.
// ---------------------------------------------------------------------------
class SubdistCIFModel : public DirectCIFModel {
 public:
  // cause: 0 = censored, `interest` = cause of interest, anything else = competing event.
  // spec must have K = 1. grid > 0: the censoring distribution is coarsened to at most `grid`
  // support points (groups of consecutive censoring times replaced by their mass-weighted mean).
  SubdistCIFModel(std::size_t n, const double* time, const int* cause, const double* X, const DirectSpec& spec,
                  int interest, int grid)
      : DirectCIFModel(n, time, cause, X, spec), kind_(n, 0), first_(n, 0), denom_(n, 1.0) {
    std::vector<std::size_t> ord(n);
    for (std::size_t i = 0; i < n; ++i) { ord[i] = i; kind_[i] = cause[i] == 0 ? 0 : (cause[i] == interest ? 1 : 2); }
    std::sort(ord.begin(), ord.end(), [time](std::size_t a, std::size_t b) { return time[a] < time[b]; });
    // Kaplan-Meier estimate of the censoring distribution, failures before censorings at ties
    std::vector<double> ct, mass;
    double G = 1.0;
    for (std::size_t a = 0; a < n;) {
      std::size_t b = a;
      double dc = 0.0, de = 0.0;
      while (b < n && time[ord[b]] == time[ord[a]]) { if (kind_[ord[b]] == 0) dc += 1.0; else de += 1.0; ++b; }
      if (dc > 0.0) {
        const double atrisk = static_cast<double>(n - a) - de;
        const double Gn = G * (1.0 - dc / atrisk);
        ct.push_back(time[ord[a]]); mass.push_back(G - Gn);
        kc_.push_back(time[ord[a]]); kR_.push_back(atrisk); kD_.push_back(dc); kG_.push_back(Gn);
        G = Gn;
      }
      a = b;
    }
    if (grid > 0 && ct.size() > static_cast<std::size_t>(grid)) {
      std::vector<double> gt, gm;
      const std::size_t m = ct.size(), g = static_cast<std::size_t>(grid);
      for (std::size_t q = 0; q < g; ++q) {
        const std::size_t lo = q * m / g, hi = (q + 1) * m / g;
        double w = 0.0, tw = 0.0;
        for (std::size_t j = lo; j < hi; ++j) { w += mass[j]; tw += mass[j] * ct[j]; }
        if (w > 0.0) { gt.push_back(tw / w); gm.push_back(w); }
      }
      ct.swap(gt); mass.swap(gm);
    }
    if (G > 1e-14) { ct.push_back(tmax_); mass.push_back(G); }      // mass beyond the last censoring time
    grid_ = grid > 0 && ct.size() != kc_.size();
    ct_ = ct; mass_ = mass;
    const std::size_t m = ct_.size();
    std::vector<double> suffix(m + 1, 0.0);
    for (std::size_t j = m; j-- > 0;) suffix[j] = suffix[j + 1] + mass_[j];
    for (std::size_t i = 0; i < n; ++i) {
      if (kind_[i] != 2) continue;
      first_[i] = static_cast<std::size_t>(std::lower_bound(ct_.begin(), ct_.end(), time[i]) - ct_.begin());
      denom_[i] = suffix[first_[i]];
    }
  }

  std::size_t n_support() const { return ct_.size(); }

  void eval(const std::vector<double>& th, Eval& out, bool need_hess) const {
    const std::size_t P = P_;
    out.ll = 0.0;
    out.g.assign(P, 0.0);
    if (out.H.n != P) out.H = Mat(P);
    out.H.zero();
    prepare(th);
    for (std::size_t i = 0; i < n_; ++i)
      if (!contribution(th, i, out, need_hess)) { bad(out); return; }
    if (!std::isfinite(out.ll)) { bad(out); return; }
    if (need_hess)
      for (std::size_t a = 0; a < P; ++a)
        for (std::size_t b = 0; b <= a; ++b) { out.H(a, b) = -out.H(a, b); out.H(b, a) = out.H(a, b); }
  }

  // sum_i s_i s_i' with s_i the score of subject i. Returns false outside the support.
  bool meat(const std::vector<double>& th, Mat& M) const {
    const std::size_t P = P_;
    M = Mat(P);
    M.zero();
    Eval e;
    e.H = Mat(P);
    prepare(th);
    for (std::size_t i = 0; i < n_; ++i) {
      e.ll = 0.0;
      e.g.assign(P, 0.0);
      if (!contribution(th, i, e, false)) return false;
      for (std::size_t a = 0; a < P; ++a) for (std::size_t b = 0; b < P; ++b) M(a, b) += e.g[a] * e.g[b];
    }
    return true;
  }

  // sum_i (s_i + psi_i)(s_i + psi_i)': the score of subject i plus the effect of subject i on the
  // estimated censoring distribution, i.e. the exact derivative of the total score with respect
  // to a case weight of subject i with the Kaplan-Meier estimate of G recomputed with that
  // weight (infinitesimal jackknife). With r indexing the censoring times c_r, R_r the number at
  // risk for censoring and D_r the number censored at c_r, h_r = log(1 - D_r / R_r):
  //   psi_i = sum_r Q_r dh_r/dv_i,  dh_r/dv_i = 1(i in R_r) D_r / {R_r (R_r - D_r)} - 1(i censored at c_r) / (R_r - D_r),
  //   Q_r   = sum_{k competing, T_k <= c_r} sum_{r' >= r} {G(c_r') / G(T_k-)} (b_k,r'+1 - b_k,r'),
  // b_k,r = d log{1 - F(c_r; z_k)} / d theta (b at the point beyond the last censoring time: the
  // largest observed time, which carries the remaining mass). This is the discrete form of the
  // Fine & Gray (1999, Sec. 3.2) correction term. Available for the exact censoring
  // distribution (grid = 0); returns false otherwise or outside the support.
  bool meat_ipcw(const std::vector<double>& th, Mat& M) const {
    const std::size_t P = P_;
    M = Mat(P);
    M.zero();
    if (grid_) return false;
    const std::size_t mc = kc_.size(), mt = ct_.size();         // mt = mc + 1 if the tmax point exists
    const bool spline = s_.baseline[0] == BASE_SPLINE;
    const std::size_t o = off_[0], nb = nb_[0];
    const bool af = s_.alpha_free[0] != 0;
    const std::size_t ia = o + nb, b0 = o + nb + (af ? 1 : 0), p = s_.p;
    const double a = alpha_of(th, 0);
    // quantities at the support points that do not depend on the subject
    std::vector<D2> U(mt, D2(0.0));
    std::vector<double> sv(mt, 0.0), zc;
    if (spline) {
      zc.assign(mt * nb, 0.0);
      std::vector<double> dz(nb);
      for (std::size_t j = 0; j < mt; ++j) {
        rcs_basis(s_.knots[0], std::log(ct_[j]), &zc[j * nb], dz.data());
        double v = 0.0;
        for (std::size_t k = 0; k < nb; ++k) v += th[o + k] * zc[j * nb + k];
        sv[j] = v;
      }
    } else {
      D2 q[3];
      std::vector<double> z, dz;
      locals<D2>(th, 0, 1.0, nullptr, nullptr, q, z, dz);
      const D2 al = D2::variable(3, a);
      for (std::size_t j = 0; j < mt; ++j) if (base_u<D2>(s_.baseline[0], ct_[j], q, al, U[j], nullptr) != 0 || !std::isfinite(U[j].v)) return false;
    }
    // G(c_r) after the censoring at c_r, r = 0..mc-1
    std::vector<double> Q(mc * P, 0.0), b(P), bn(P), acc(P);
    auto grad_logS = [&](std::size_t j, const double* x, double E, double* g) -> bool {
      std::fill(g, g + P, 0.0);
      double hd[6];
      if (spline) {
        const double w = E * std::exp(sv[j]);
        if (!(a * w > -1.0) || !std::isfinite(w)) return false;
        h_derivs(w, a, hd);
        for (std::size_t k = 0; k < nb; ++k) g[o + k] = -hd[1] * w * zc[j * nb + k];
        if (af) g[ia] = -hd[3];
        for (std::size_t r = 0; r < p; ++r) g[b0 + r] = -hd[1] * w * x[r];
      } else {
        const D2& u = U[j];
        const double w = E * u.v;
        if (!(a * w > -1.0) || !std::isfinite(w)) return false;
        h_derivs(w, a, hd);
        for (std::size_t k = 0; k < nb; ++k) g[o + k] = -hd[1] * E * u.g[k];
        if (af) g[ia] = -(hd[1] * E * u.g[3] + hd[3]);
        for (std::size_t r = 0; r < p; ++r) g[b0 + r] = -hd[1] * w * x[r];
      }
      return true;
    };
    for (std::size_t i = 0; i < n_; ++i) {
      if (kind_[i] != 2) continue;
      const std::size_t f = static_cast<std::size_t>(std::lower_bound(kc_.begin(), kc_.end(), time_[i]) - kc_.begin());
      if (f >= mc) continue;
      const double Gprev = f == 0 ? 1.0 : kG_[f - 1];
      if (!(Gprev > 0.0)) continue;
      const double* x = X_ + i * p;
      const double E = std::exp(eta_of(th, 0, x));
      // b at the point after the last censoring time (tmax), or zero if G reaches 0
      if (mt > mc) { if (!grad_logS(mc, x, E, bn.data())) return false; } else std::fill(bn.begin(), bn.end(), 0.0);
      std::fill(acc.begin(), acc.end(), 0.0);
      for (std::size_t r = mc; r-- > f;) {
        if (!grad_logS(r, x, E, b.data())) return false;
        const double rho = kG_[r] / Gprev;
        double* Qr = &Q[r * P];
        for (std::size_t c = 0; c < P; ++c) { acc[c] += rho * (bn[c] - b[c]); Qr[c] += acc[c]; }
        b.swap(bn);
      }
    }
    // prefix sums of Q_r D_r / {R_r (R_r - D_r)} over the censoring times
    std::vector<double> pre((mc + 1) * P, 0.0);
    for (std::size_t r = 0; r < mc; ++r) {
      const double RD = kR_[r] - kD_[r];
      const double c1 = RD > 0.0 ? kD_[r] / (kR_[r] * RD) : 0.0;
      for (std::size_t c = 0; c < P; ++c) pre[(r + 1) * P + c] = pre[r * P + c] + c1 * Q[r * P + c];
    }
    Eval e;
    e.H = Mat(P);
    prepare(th);
    std::vector<double> sc(P);
    for (std::size_t i = 0; i < n_; ++i) {
      e.ll = 0.0;
      e.g.assign(P, 0.0);
      if (!contribution(th, i, e, false)) return false;
      const std::size_t r = static_cast<std::size_t>(std::lower_bound(kc_.begin(), kc_.end(), time_[i]) - kc_.begin());
      for (std::size_t c = 0; c < P; ++c) sc[c] = e.g[c] + pre[r * P + c];      // c_r < T_i
      if (kind_[i] == 0 && r < mc && kc_[r] == time_[i]) {                       // censored at c_r: in R_r and censored
        const double RD = kR_[r] - kD_[r];
        if (RD > 0.0) {
          const double c1 = kD_[r] / (kR_[r] * RD) - 1.0 / RD;
          for (std::size_t c = 0; c < P; ++c) sc[c] += c1 * Q[r * P + c];
        }
      }
      for (std::size_t aa = 0; aa < P; ++aa) for (std::size_t bb = 0; bb < P; ++bb) M(aa, bb) += sc[aa] * sc[bb];
    }
    return true;
  }

 private:
  // log{1 - F(t; x)} as a function of the local variables
  bool log_surv(const std::vector<double>& th, const double* x, double t, const double* zr, const double* dr,
                std::vector<double>& z, std::vector<double>& dz, D2& ls) const {
    D2 q[3], F;
    locals<D2>(th, 0, t, zr, dr, q, z, dz);
    if (!direct_Ff<D2>(s_.baseline[0], t, q, D2::variable(3, alpha_of(th, 0)), D2::variable(4, eta_of(th, 0, x)), &F, nullptr)) return false;
    const double S = 1.0 - F.v;
    if (!(S > 0.0)) return false;
    ls = chain(F, std::log1p(-F.v), -1.0 / S, -1.0 / (S * S));
    return true;
  }
  bool contribution(const std::vector<double>& th, std::size_t i, Eval& out, bool need_hess) const {
    const double* x = X_ + i * s_.p;
    const double t = time_[i];
    const bool spline = s_.baseline[0] == BASE_SPLINE;
    const double* zr = spline ? &zrow_[0][i * nb_[0]] : nullptr;
    const double* dr = spline ? &dzrow_[0][i * nb_[0]] : nullptr;
    std::vector<double> z, dz;
    std::vector<Term> terms;
    D2 d;
    if (kind_[i] == 1) {
      D2 q[3];
      locals<D2>(th, 0, t, zr, dr, q, z, dz);
      if (!direct_Ff<D2>(s_.baseline[0], t, q, D2::variable(3, alpha_of(th, 0)), D2::variable(4, eta_of(th, 0, x)), nullptr, &d) ||
          !std::isfinite(d.v)) return false;
      out.ll += d.v;
      make_terms(0, x, spline ? zr : z.data(), spline ? dr : dz.data(), terms);
      add_terms(terms, d, 1.0, out, need_hess, nullptr);
      return true;
    }
    if (kind_[i] == 0) {
      if (!(t > 0.0)) return true;                       // censored at time 0: no information
      if (!log_surv(th, x, t, zr, dr, z, dz, d)) return false;
      out.ll += d.v;
      make_terms(0, x, spline ? zr : z.data(), spline ? dr : dz.data(), terms);
      add_terms(terms, d, 1.0, out, need_hess, nullptr);
      return true;
    }
    // competing event: expectation over the estimated censoring distribution beyond T_i
    const double den = denom_[i];
    if (!(den > 0.0)) return true;
    if (spline && !slow_) return spline_competing(th, i, x, den, out, need_hess);
    if (spline) {                                        // the local variables depend on the time point
      for (std::size_t j = first_[i]; j < ct_.size(); ++j) {
        if (!log_surv(th, x, ct_[j], nullptr, nullptr, z, dz, d)) return false;
        const double w = mass_[j] / den;
        out.ll += w * d.v;
        make_terms(0, x, z.data(), dz.data(), terms);
        add_terms(terms, d, w, out, need_hess, nullptr);
      }
      return true;
    }
    D2 acc(0.0);
    if (!slow_) {
      if (!competing_sum(th, i, x, den, acc)) return false;
    } else {
      for (std::size_t j = first_[i]; j < ct_.size(); ++j) {
        if (!log_surv(th, x, ct_[j], nullptr, nullptr, z, dz, d)) return false;
        acc = acc + (mass_[j] / den) * d;
      }
    }
    out.ll += acc.v;
    make_terms(0, x, z.data(), dz.data(), terms);
    add_terms(terms, acc, 1.0, out, need_hess, nullptr);
    return true;
  }

 public:
  // the original evaluation of every (subject, support point) pair through the full automatic
  // differentiation, kept to check the fast path
  void set_reference_path(bool on) { slow_ = on; }

 private:
  // Fast exact evaluation of sum_j w_ij log{1 - F(c_j; z_i)} for a non-spline baseline.
  // log{1 - F(c; z)} = -A with A = h(u(c), alpha, eta) = log1p(alpha e^eta u) / alpha (= e^eta u at
  // alpha = 0). u(c_j) and its derivatives in the baseline parameters and alpha (local variables
  // 0-3) depend on c_j only, so they are computed once per evaluation (prepare); the pair (i, j)
  // then needs the closed-form derivatives of h in (u, alpha, eta) only. With alpha fixed at 0 (PH
  // link) A = e^eta u is linear in u and the sum over j is a suffix sum: O(1) per subject.
  void prepare(const std::vector<double>& th) const {
    if (slow_) return;
    const std::size_t m = ct_.size();
    if (s_.baseline[0] == BASE_SPLINE) { prepare_spline(th); return; }
    U_.assign(m, D2(0.0)); ok_.assign(m + 1, 1);
    D2 q[3];
    std::vector<double> z, dz;
    locals<D2>(th, 0, 1.0, nullptr, nullptr, q, z, dz);            // the baseline variables do not depend on t
    const D2 al = D2::variable(3, alpha_of(th, 0));
    for (std::size_t j = 0; j < m; ++j) {
      D2 u;
      const int rc = base_u<D2>(s_.baseline[0], ct_[j], q, al, u, nullptr);
      if (rc != 0 || !std::isfinite(u.v)) { ok_[j] = 0; U_[j] = D2(0.0); } else U_[j] = u;
    }
    for (std::size_t j = m; j-- > 0;) ok_[j] = ok_[j] && ok_[j + 1];
    ph_ = !s_.alpha_free[0] && alpha_of(th, 0) == 0.0;
    if (ph_) {
      suf_.assign(m + 1, D2(0.0));
      for (std::size_t j = m; j-- > 0;) suf_[j] = suf_[j + 1] + mass_[j] * U_[j];
    }
  }

  bool competing_sum(const std::vector<double>& th, std::size_t i, const double* x, double den, D2& acc) const {
    const std::size_t f = first_[i];
    if (!ok_[f]) return false;
    const double eta = eta_of(th, 0, x), E = std::exp(eta);
    acc = D2(0.0);
    if (ph_) {                                              // -A = -e^eta sum_j w_j u_j / den
      const D2& S = suf_[f];
      const double c = -E / den;
      acc.v = c * S.v;
      for (int k = 0; k < 4; ++k) {
        acc.g[k] = c * S.g[k];
        acc.h[k][4] = acc.h[4][k] = c * S.g[k];
        for (int l = 0; l < 4; ++l) acc.h[k][l] = c * S.h[k][l];
      }
      acc.g[4] = c * S.v; acc.h[4][4] = c * S.v;
      return std::isfinite(acc.v);
    }
    const double a = alpha_of(th, 0);
    const bool af = s_.alpha_free[0] != 0;                  // local variable 3 (alpha) is live only if estimated
    const int nl = af ? 4 : 3;
    // series coefficients in v for |v| < 0.05 (13 terms: truncation error below 0.05^13 / 14 ~ 1e-18)
    static const struct Series {
      double h[13], a[13], aa[13];
      Series() {
        for (int k = 0; k < 13; ++k) {
          const double sg = (k % 2 == 0) ? 1.0 : -1.0;
          h[k] = sg / (k + 1.0);                                       // log1p(v) / v
          a[k] = -sg * (k + 1.0) / (k + 2.0);                          // {v/(1+v) - log1p(v)} / v^2
          aa[k] = sg * (k + 1.0) * (k + 2.0) / (k + 3.0);              // {-v^2/(1+v)^2 - 2(v/(1+v) - log1p(v))} / v^3
        }
      }
    } ser;
    double G[5] = {0, 0, 0, 0, 0}, Hs[5][5] = {};
    for (std::size_t j = f; j < ct_.size(); ++j) {
      const D2& u = U_[j];
      const double w = E * u.v, v = a * w;
      if (!(v > -1.0) || !std::isfinite(v)) return false;
      const double op = 1.0 + v, iop = 1.0 / op;
      double h, ha, haa;
      if (std::fabs(v) < 0.05) {
        double sh = ser.h[12], sa = ser.a[12], saa = ser.aa[12];
        for (int k = 11; k >= 0; --k) { sh = sh * v + ser.h[k]; sa = sa * v + ser.a[k]; saa = saa * v + ser.aa[k]; }
        h = w * sh; ha = w * w * sa; haa = w * w * w * saa;
      } else {
        const double l = std::log1p(v), gg = v * iop - l, ia = 1.0 / a;
        h = l * ia; ha = gg * ia * ia; haa = (-(v * v) * iop * iop - 2.0 * gg) * ia * ia * ia;
      }
      const double hw = iop, hww = -a * iop * iop, haw = -w * iop * iop;
      // derivatives of h in (u, alpha, eta) from those in (w = e^eta u, alpha)
      const double c = -mass_[j] / den;                     // log S = -A
      const double hu = c * E * hw, huu = c * E * E * hww, hue = c * (E * hw + E * w * hww), hua = c * E * haw;
      acc.v += c * h;
      for (int k = 0; k < nl; ++k) {
        G[k] += hu * u.g[k];
        for (int l = 0; l <= k; ++l) Hs[k][l] += huu * u.g[k] * u.g[l] + hu * u.h[k][l];
        Hs[4][k] += hue * u.g[k];
      }
      G[4] += c * w * hw;
      Hs[4][4] += c * (w * hw + w * w * hww);
      if (af) {
        G[3] += c * ha;
        for (int k = 0; k < 4; ++k) Hs[3 > k ? 3 : k][3 > k ? k : 3] += hua * u.g[k];
        Hs[3][3] += hua * u.g[3] + c * haa;
        Hs[4][3] += c * w * haw;
      }
    }
    for (int k = 0; k < 5; ++k) {
      acc.g[k] = G[k];
      for (int l = 0; l <= k; ++l) acc.h[k][l] = acc.h[l][k] = Hs[k][l];
    }
    return std::isfinite(acc.v);
  }

  // Spline baseline: log u(c_j) = s_j = gamma'z_j, so log{1 - F} = -h(w, alpha) with
  // w = e^(eta + s_j); s_j and eta enter through w only. The gamma derivatives are linear forms
  // in the basis row z_j, accumulated per subject (O(nb^2) per pair, no automatic
  // differentiation); with alpha fixed at 0, h = w and the sums over j are suffix sums.
  void prepare_spline(const std::vector<double>& th) const {
    const std::size_t m = ct_.size(), nb = nb_[0], o = off_[0];
    if (zc_.size() != m * nb) {                           // basis rows at the support points (fixed)
      zc_.assign(m * nb, 0.0); std::vector<double> dz(nb);
      for (std::size_t j = 0; j < m; ++j) rcs_basis(s_.knots[0], std::log(ct_[j]), &zc_[j * nb], dz.data());
      dsc_.assign(m * nb, 0.0);
      for (std::size_t j = 0; j < m; ++j) rcs_basis(s_.knots[0], std::log(ct_[j]), dz.data(), &dsc_[j * nb]);
    }
    sv_.assign(m, 0.0); ok_.assign(m + 1, 1);
    for (std::size_t j = 0; j < m; ++j) {
      double sv = 0.0, dv = 0.0;
      for (std::size_t k = 0; k < nb; ++k) { sv += th[o + k] * zc_[j * nb + k]; dv += th[o + k] * dsc_[j * nb + k]; }
      sv_[j] = sv;
      if (!(dv > 0.0) || !std::isfinite(std::exp(sv))) ok_[j] = 0;   // as base_u: the log cumulative hazard must increase
    }
    for (std::size_t j = m; j-- > 0;) ok_[j] = ok_[j] && ok_[j + 1];
    ph_ = !s_.alpha_free[0] && alpha_of(th, 0) == 0.0;
    if (ph_) {                                            // suffix sums of mass_j e^s_j (1, z_j, z_j z_j')
      const std::size_t L = 1 + nb + nb * nb;
      ssuf_.assign((m + 1) * L, 0.0);
      for (std::size_t j = m; j-- > 0;) {
        const double* nx = &ssuf_[(j + 1) * L];
        double* cu = &ssuf_[j * L];
        const double e = mass_[j] * std::exp(sv_[j]);
        const double* z = &zc_[j * nb];
        cu[0] = nx[0] + e;
        for (std::size_t k = 0; k < nb; ++k) {
          cu[1 + k] = nx[1 + k] + e * z[k];
          for (std::size_t l = 0; l < nb; ++l) cu[1 + nb + k * nb + l] = nx[1 + nb + k * nb + l] + e * z[k] * z[l];
        }
      }
    }
  }

  bool spline_competing(const std::vector<double>& th, std::size_t i, const double* x, double den, Eval& out,
                        bool need_hess) const {
    const std::size_t f = first_[i], m = ct_.size(), nb = nb_[0], o = off_[0];
    if (!ok_[f]) return false;
    const bool af = s_.alpha_free[0] != 0;
    const double eta = eta_of(th, 0, x), E = std::exp(eta);
    // sums over j of c_j times: h, h_s (= h_eta), h_ss, h_sa, h_a, h_aa, and the z-weighted versions
    double S0 = 0.0, Ss = 0.0, Sss = 0.0, Ssa = 0.0, Sa = 0.0, Saa = 0.0;
    std::vector<double>& Gz = wz1_; std::vector<double>& Hzz = wzz_; std::vector<double>& Zss = wz2_; std::vector<double>& Zsa = wz3_;
    Gz.assign(nb, 0.0); Zss.assign(nb, 0.0); Zsa.assign(nb, 0.0); Hzz.assign(nb * nb, 0.0);
    if (ph_) {
      const std::size_t L = 1 + nb + nb * nb;
      const double* su = &ssuf_[f * L];
      const double c = -E / den;
      S0 = Ss = Sss = c * su[0];
      for (std::size_t k = 0; k < nb; ++k) { Gz[k] = Zss[k] = c * su[1 + k]; for (std::size_t l = 0; l < nb; ++l) Hzz[k * nb + l] = c * su[1 + nb + k * nb + l]; }
    } else {
      const double a = alpha_of(th, 0);
      for (std::size_t j = f; j < m; ++j) {
        const double w = E * std::exp(sv_[j]), v = a * w;
        if (!(v > -1.0) || !std::isfinite(v)) return false;
        double hd[6];
        h_derivs(w, a, hd);                                // h, hw, hww, ha, haw, haa
        const double c = -mass_[j] / den;
        const double hs = hd[1] * w, hss = hd[2] * w * w + hd[1] * w, hsa = hd[4] * w;
        S0 += c * hd[0]; Ss += c * hs; Sss += c * hss; Ssa += c * hsa; Sa += c * hd[3]; Saa += c * hd[5];
        const double* z = &zc_[j * nb];
        for (std::size_t k = 0; k < nb; ++k) {
          Gz[k] += c * hs * z[k]; Zss[k] += c * hss * z[k]; Zsa[k] += c * hsa * z[k];
          if (need_hess) for (std::size_t l = 0; l <= k; ++l) Hzz[k * nb + l] += c * hss * z[k] * z[l];
        }
      }
    }
    if (!std::isfinite(S0)) return false;
    out.ll += S0;
    const std::size_t ia = o + nb, b0 = o + nb + (af ? 1 : 0), p = s_.p;
    for (std::size_t k = 0; k < nb; ++k) out.g[o + k] += Gz[k];
    if (af) out.g[ia] += Sa;
    for (std::size_t r = 0; r < p; ++r) out.g[b0 + r] += Ss * x[r];
    if (!need_hess) return true;
    for (std::size_t k = 0; k < nb; ++k) for (std::size_t l = 0; l <= k; ++l) out.H(o + k, o + l) += Hzz[k * nb + l];
    if (af) {
      for (std::size_t k = 0; k < nb; ++k) out.H(ia, o + k) += Zsa[k];
      out.H(ia, ia) += Saa;
    }
    for (std::size_t r = 0; r < p; ++r) {
      for (std::size_t k = 0; k < nb; ++k) out.H(b0 + r, o + k) += Zss[k] * x[r];
      if (af) out.H(b0 + r, ia) += Ssa * x[r];
      for (std::size_t q = 0; q <= r; ++q) out.H(b0 + r, b0 + q) += Sss * x[r] * x[q];
    }
    return true;
  }

  // h(w, alpha) = log1p(alpha w) / alpha and its derivatives: d = (h, h_w, h_ww, h_a, h_wa, h_aa).
  // A series in v = alpha w is used for |v| < 0.05 (13 terms; truncation below 1e-18).
  static void h_derivs(double w, double a, double* d) {
    static const struct Series {
      double h[13], a[13], aa[13];
      Series() {
        for (int k = 0; k < 13; ++k) {
          const double sg = (k % 2 == 0) ? 1.0 : -1.0;
          h[k] = sg / (k + 1.0); a[k] = -sg * (k + 1.0) / (k + 2.0); aa[k] = sg * (k + 1.0) * (k + 2.0) / (k + 3.0);
        }
      }
    } ser;
    const double v = a * w, iop = 1.0 / (1.0 + v);
    if (std::fabs(v) < 0.05) {
      double sh = ser.h[12], sa = ser.a[12], saa = ser.aa[12];
      for (int k = 11; k >= 0; --k) { sh = sh * v + ser.h[k]; sa = sa * v + ser.a[k]; saa = saa * v + ser.aa[k]; }
      d[0] = w * sh; d[3] = w * w * sa; d[5] = w * w * w * saa;
    } else {
      const double l = std::log1p(v), gg = v * iop - l, ia = 1.0 / a;
      d[0] = l * ia; d[3] = gg * ia * ia; d[5] = (-(v * v) * iop * iop - 2.0 * gg) * ia * ia * ia;
    }
    d[1] = iop; d[2] = -a * iop * iop; d[4] = -w * iop * iop;
  }

  bool slow_ = false;
  mutable std::vector<double> zc_, dsc_, sv_, ssuf_, wz1_, wz2_, wz3_, wzz_;
  mutable std::vector<D2> U_, suf_;
  mutable std::vector<char> ok_;
  mutable bool ph_ = false;

  std::vector<int> kind_;                 // 0 censored, 1 cause of interest, 2 competing event
  std::vector<std::size_t> first_;        // competing events: first support point >= T_i
  std::vector<double> denom_;             // competing events: G(T_i-)
  std::vector<double> ct_, mass_;         // support and masses of the estimated censoring distribution
  std::vector<double> kc_, kR_, kD_, kG_; // exact Kaplan-Meier of censoring: times, at risk, censored, G after
  bool grid_ = false;                     // the support was coarsened
};

struct SubdistFit {
  FitResult fit;                          // fit.var = sandwich variance
  Mat var_model;                          // inverse of the information of the weighted likelihood
};

// cause01: 1 = cause of interest, 2 = competing event, 0 = censored (used for the start values).
inline SubdistFit subdist_fit(const SubdistCIFModel& m, std::size_t n, const double* time, const int* cause12,
                              const std::vector<double>* init, int max_iter, double tol) {
  const DirectSpec& s = m.spec();
  const std::size_t P = s.nparam();
  std::vector<double> th;
  if (init) th = *init;
  else {                                  // block of cause 1 of the joint start values
    DirectSpec s2 = s;
    s2.K = 2;
    s2.alpha_free.assign(2, s.alpha_free[0]); s2.alpha_value.assign(2, s.alpha_value[0]);
    s2.baseline.assign(2, s.baseline[0]); s2.knots.assign(2, s.knots[0]);
    const std::vector<double> th2 = direct_start(n, time, cause12, s2);
    th.assign(th2.begin(), th2.begin() + static_cast<std::ptrdiff_t>(P));
  }
  SubdistFit out;
  out.fit = staged_fit(m, th, direct_roles(s), init != nullptr, false, max_iter, tol,
                       [&s](std::vector<double>& v) { v[base_level_index(s.baseline[0])] -= 0.5; });
  out.var_model = out.fit.var;
  if (out.fit.var.n == P) {
    Mat M;
    if (m.meat(out.fit.beta, M)) {
      Mat V(P);
      const Mat& A = out.var_model;
      std::vector<double> AM(P * P, 0.0);
      for (std::size_t a = 0; a < P; ++a) for (std::size_t b = 0; b < P; ++b) { double v = 0.0; for (std::size_t c = 0; c < P; ++c) v += A(a, c) * M(c, b); AM[a * P + b] = v; }
      for (std::size_t a = 0; a < P; ++a) for (std::size_t b = 0; b < P; ++b) { double v = 0.0; for (std::size_t c = 0; c < P; ++c) v += AM[a * P + c] * A(c, b); V(a, b) = v; }
      out.fit.var = V;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// 2. Constrained model for two causes: Shi, Cheng & Jeong (2013), Biom J 55:82-96.
//
//   F_1(t; z) : any direct model with a plateau below one (directcif.hpp),
//   F_2(t; z) = p_2(z) G(t),   p_2(z) = 1 - F_1(infinity; z)                         (their eq. 7)
//
// so that F_1(infinity; z) + F_2(infinity; z) = 1 for every z (their eq. 1). G is the
// distribution function of the time to a cause 2 event given that it occurs; it has no
// covariates. The paper uses the modified logistic form, G(t) = [e^{b(t-c)} - e^{-bc}] /
// [1 + e^{b(t-c)}]; the Weibull, log-logistic and log-normal forms are offered in addition.
// Likelihood (their eq. 5): f_1 for cause 1 events, p_2(z) g(t) for cause 2 events and
// 1 - F_1(t; z) - p_2(z) G(t) for censored observations.
// Parameters: block of cause 1 as in DirectSpec, then the two parameters of G:
//   logistic (log b, c); Weibull, log-logistic (log lambda, log kappa); log-normal (mu, log sigma).
// ---------------------------------------------------------------------------
enum CondFamily : int { COND_LOGISTIC = 0, COND_WEIBULL = 1, COND_LOGLOGISTIC = 2, COND_LOGNORMAL = 3 };
constexpr int COND_COUNT = 4;

// G(t) and log g(t) with a, b the two parameters. t = infinity: G = 1.
template <class T>
void cond_Gg(int family, double t, const T& a, const T& b, T& G, T* logg) {
  if (std::isinf(t)) { G = T(1.0); return; }
  if (family == COND_LOGISTIC) {
    const T rate = t_exp(a), x = rate * ((-1.0) * b + t), bt = t * rate;
    G = t_sigmoid(x) * ((-1.0) * t_expm1((-1.0) * bt));
    if (logg) *logg = a + x + t_softplus((-1.0) * (rate * b)) + (-2.0) * t_softplus(x);
    return;
  }
  const double lt = std::log(t);
  if (family == COND_LOGNORMAL) {
    const T z = t_exp((-1.0) * b) * ((-1.0) * a + lt);
    G = t_normcdf(z);
    if (logg) *logg = (-0.5) * (z * z) + (-1.0) * b + (-0.9189385332046727 - lt);
    return;
  }
  const T lw = t_exp(b) * (a + lt);
  if (family == COND_WEIBULL) {
    const T w = t_exp(lw);
    G = (-1.0) * t_expm1((-1.0) * w);
    if (logg) *logg = b + lw + (-1.0) * w + (-lt);
  } else {
    G = t_sigmoid(lw);
    if (logg) *logg = b + lw + (-2.0) * t_softplus(lw) + (-lt);
  }
}

class ConstrainedCIFModel : public DirectCIFModel {
 public:
  // cause: 0 = censored, 1 = modelled directly, 2 = competing cause. spec must have K = 1 and a
  // baseline with a plateau (not WEIBULL, not SPLINE).
  ConstrainedCIFModel(std::size_t n, const double* time, const int* cause, const double* X, const DirectSpec& spec, int family)
      : DirectCIFModel(n, time, cause, X, spec), family_(family), P1_(spec.nparam()) {}

  std::size_t nparam() const { return P1_ + 2; }
  int family() const { return family_; }

  // CIF of cause k (1 or 2) at time t (may be infinite) for covariate row x
  double cif12(const std::vector<double>& th, int k, const double* x, double t) const {
    if (t <= 0.0) return 0.0;
    if (k == 1) return cif(th, 0, x, t);
    const double p2 = 1.0 - cif(th, 0, x, std::numeric_limits<double>::infinity());
    double G = 1.0;
    cond_Gg<double>(family_, t, th[P1_], th[P1_ + 1], G, nullptr);
    return p2 * G;
  }
  // gradient of the CIF of cause k with respect to all parameters; returns the CIF
  double cif12_grad(const std::vector<double>& th, int k, const double* x, double t, std::vector<double>& grad) const {
    std::vector<double> g1;
    grad.assign(P1_ + 2, 0.0);
    if (t <= 0.0) return 0.0;
    if (k == 1) {
      const double F = cif_grad(th, 0, x, t, g1);
      for (std::size_t j = 0; j < P1_; ++j) grad[j] = g1[j];
      return F;
    }
    const double Finf = cif_grad(th, 0, x, std::numeric_limits<double>::infinity(), g1);
    D2 G;
    cond_Gg<D2>(family_, t, D2::variable(0, th[P1_]), D2::variable(1, th[P1_ + 1]), G, nullptr);
    for (std::size_t j = 0; j < P1_; ++j) grad[j] = -g1[j] * G.v;
    grad[P1_] = (1.0 - Finf) * G.g[0];
    grad[P1_ + 1] = (1.0 - Finf) * G.g[1];
    return (1.0 - Finf) * G.v;
  }

  void eval(const std::vector<double>& th, Eval& out, bool need_hess) const {
    const std::size_t P = P1_ + 2;
    const double inf = std::numeric_limits<double>::infinity();
    out.ll = 0.0;
    out.g.assign(P, 0.0);
    if (out.H.n != P) out.H = Mat(P);
    out.H.zero();
    std::vector<double> z, dz, gs(P), gB(P);
    std::vector<Term> t1, t2;
    t2.push_back({P1_, 0, 1.0});
    t2.push_back({P1_ + 1, 1, 1.0});
    const D2 a = D2::variable(0, th[P1_]), b = D2::variable(1, th[P1_ + 1]);
    D2 q[3], lf, Ft, Finf, G, logg;
    for (std::size_t i = 0; i < n_; ++i) {
      const double* x = X_ + i * s_.p;
      const double t = time_[i];
      const int c = cause_[i];
      const D2 al = D2::variable(3, alpha_of(th, 0)), et = D2::variable(4, eta_of(th, 0, x));
      locals<D2>(th, 0, t, nullptr, nullptr, q, z, dz);
      make_terms(0, x, z.data(), dz.data(), t1);
      if (c == 1) {
        if (!direct_Ff<D2>(s_.baseline[0], t, q, al, et, nullptr, &lf) || !std::isfinite(lf.v)) { bad(out); return; }
        out.ll += lf.v;
        add_terms(t1, lf, 1.0, out, need_hess, nullptr);
        continue;
      }
      if (!direct_Ff<D2>(s_.baseline[0], inf, q, al, et, &Finf, nullptr)) { bad(out); return; }
      const double p2 = 1.0 - Finf.v;
      if (!(p2 > 0.0)) { bad(out); return; }
      if (c == 2) {                                       // log p_2(z) + log g(t)
        cond_Gg<D2>(family_, t, a, b, G, &logg);
        if (!std::isfinite(logg.v)) { bad(out); return; }
        const D2 lp = chain(Finf, std::log(p2), -1.0 / p2, -1.0 / (p2 * p2));
        out.ll += lp.v + logg.v;
        add_terms(t1, lp, 1.0, out, need_hess, nullptr);
        add_terms(t2, logg, 1.0, out, need_hess, nullptr);
        continue;
      }
      if (!(t > 0.0)) continue;                           // censored at time 0
      // censored: S = 1 - F_1(t) - p_2 G(t) = 1 - C1 - C2 with C1 = F_1(t) - G F_1(inf) + G, C2 = p_2 G
      if (!direct_Ff<D2>(s_.baseline[0], t, q, al, et, &Ft, nullptr)) { bad(out); return; }
      cond_Gg<D2>(family_, t, a, b, G, nullptr);
      const double S = 1.0 - Ft.v - p2 * G.v;
      if (!(S > 0.0)) { bad(out); return; }
      out.ll += std::log(S);
      std::fill(gs.begin(), gs.end(), 0.0);
      const D2 C1 = Ft + (-G.v) * Finf;                   // derivatives with respect to the cause 1 block (G held fixed)
      const D2 C2 = p2 * G;                               // derivatives with respect to the parameters of G (p_2 held fixed)
      add_terms(t1, C1, -1.0 / S, out, need_hess, gs.data());
      add_terms(t2, C2, -1.0 / S, out, need_hess, gs.data());
      if (need_hess) {
        // mixed second derivatives of S: -(d p_2 / d theta_1)(d G / d theta_2)' = +(d F_1(inf))(d G)'
        std::fill(gB.begin(), gB.end(), 0.0);
        for (const Term& e : t1) gB[e.idx] += Finf.g[e.loc] * e.coef;
        for (std::size_t j = 0; j < P1_; ++j) {
          if (gB[j] == 0.0) continue;
          out.H(P1_, j) += gB[j] * G.g[0] / S;
          out.H(P1_ + 1, j) += gB[j] * G.g[1] / S;
        }
        for (std::size_t r = 0; r < P; ++r) {             // -(grad S)(grad S)' / S^2
          if (gs[r] == 0.0) continue;
          const double gr = gs[r] / (S * S);
          for (std::size_t cc = 0; cc <= r; ++cc) out.H(r, cc) -= gr * gs[cc];
        }
      }
    }
    if (!std::isfinite(out.ll)) { bad(out); return; }
    if (need_hess)
      for (std::size_t r = 0; r < P; ++r)
        for (std::size_t cc = 0; cc <= r; ++cc) { out.H(r, cc) = -out.H(r, cc); out.H(cc, r) = out.H(r, cc); }
  }

 private:
  int family_;
  std::size_t P1_;
};

// Start values: cause 1 block from the joint start values; G from the observed cause 2 event times.
inline std::vector<double> constrained_start(std::size_t n, const double* time, const int* cause, const DirectSpec& s, int family) {
  DirectSpec s2 = s;
  s2.K = 2;
  s2.alpha_free.assign(2, s.alpha_free[0]); s2.alpha_value.assign(2, s.alpha_value[0]);
  s2.baseline.assign(2, s.baseline[0]); s2.knots.assign(2, s.knots[0]);
  const std::vector<double> th2 = direct_start(n, time, cause, s2);
  const std::size_t P1 = s.nparam();
  std::vector<double> th(th2.begin(), th2.begin() + static_cast<std::ptrdiff_t>(P1));
  std::vector<double> t2;
  for (std::size_t i = 0; i < n; ++i) if (cause[i] == 2 && time[i] > 0.0) t2.push_back(time[i]);
  double a = 0.0, b = 0.0;
  if (!t2.empty()) {
    std::sort(t2.begin(), t2.end());
    const double med = t2[t2.size() / 2];
    double ml = 0.0, sl = 0.0;
    for (double v : t2) ml += std::log(v);
    ml /= static_cast<double>(t2.size());
    for (double v : t2) sl += (std::log(v) - ml) * (std::log(v) - ml);
    sl = std::sqrt(sl / static_cast<double>(t2.size()));
    if (!(sl > 0.05)) sl = 0.5;
    switch (family) {
      case COND_LOGISTIC: a = std::log(2.0 / med); b = 0.5 * med; break;     // rate 2/median, centre below the median
      case COND_WEIBULL: a = -std::log(med) + std::log(std::log(2.0)) * sl / 1.2825; b = std::log(1.2825 / sl); break;
      case COND_LOGLOGISTIC: a = -ml; b = std::log(1.8138 / sl); break;       // sd of a logistic variable = pi / sqrt(3)
      default: a = ml; b = std::log(sl); break;
    }
  }
  th.push_back(a); th.push_back(b);
  return th;
}

inline FitResult constrained_fit(const ConstrainedCIFModel& m, std::size_t n, const double* time, const int* cause,
                                 const std::vector<double>* init, int max_iter, double tol) {
  const DirectSpec& s = m.spec();
  std::vector<int> role = direct_roles(s);
  role.push_back(0); role.push_back(0);
  return staged_fit(m, init ? *init : constrained_start(n, time, cause, s, m.family()), role, init != nullptr, false, max_iter, tol,
                    [&s](std::vector<double>& v) { v[base_level_index(s.baseline[0])] -= 0.5; });
}

// ---------------------------------------------------------------------------
// 3. Constrained model for K >= 2 causes (extension of Shi, Cheng & Jeong 2013).
//
//   F_1(t; z) : direct model with a plateau p_1(z) = F_1(infinity; z) < 1 (directcif.hpp),
//   F_k(t; z) = {1 - p_1(z)} pi_k(z) G_k(t),   k = 2, ..., K,
//   pi_k(z)   = exp(l_k) / sum_{j=2..K} exp(l_j),  l_2 = 0,  l_k = gamma_k + delta_k' z (k >= 3),
//
// so that sum_k F_k(infinity; z) = 1 for every z: the cause of interest keeps the direct
// regression model of Shi et al. (effects on its cumulative incidence), and the remaining
// long-run risk is shared among the competing causes by a multinomial logit, each with its own
// distribution G_k of the time to the event given that it occurs (no covariates, as in Shi et
// al.). With K = 2, pi_2 = 1 and the model is that of ConstrainedCIFModel.
// Parameters: block of cause 1 (DirectSpec), then (a_k, b_k) of G_k for k = 2..K, then
// (gamma_k, delta_k) for k = 3..K.
// Likelihood: f_1(t) for cause 1 events; {1 - p_1} pi_k g_k(t) for cause k events;
// S = 1 - F_1(t) - {1 - p_1} sum_k pi_k G_k(t) for censored observations.
//
// Covariates of the time distributions (optional, CondCovariates). G_k may depend on covariates
// through an accelerated failure time model with, optionally, a regression for its shape:
//   G_k(t | w, v) = G_k0(t e^{-u}; shape shifted by s),  u = eta_k' w,  s = zeta_k' v,
// so that e^{eta} is a time ratio (> 1: events of type k come later). In the two parameters of
// cond_Gg this is, exactly,
//   Weibull, log-logistic (a = log lambda, b = log kappa):  a - u,  b + s
//   log-normal            (a = mu,         b = log sigma):  a + u,  b + s
//   modified logistic     (a = log rate,   b = centre c):   a - u + s,  b e^u
// (for the logistic form s multiplies the rate with the centre fixed, i.e. it changes the
// dimensionless steepness rate * c). Every G_k(. | w, v) is a distribution function, so the
// incidences still add up to one for every subject. The proportion of the competing causes and
// their shares are not affected by w and v; only the timing of the events is.
// Parameters (appended, so that the layout without these covariates is unchanged): eta_k for
// k = 2..K (the columns wcol[k-2] of W), then zeta_k for k = 2..K (the columns vcol[k-2] of V).
// Derivatives: the local variables of G_k are (a_k, b_k, u, s), slots 0..3 of D2, so that
// gradients and Hessians stay exact.
// ---------------------------------------------------------------------------
struct CondCovariates {
  const double* W = nullptr;              // n x nw, row-major: covariates of the location (time ratio)
  std::size_t nw = 0;
  const double* V = nullptr;              // n x nv, row-major: covariates of the shape
  std::size_t nv = 0;
  std::vector<std::vector<std::size_t>> wcol, vcol;   // per competing cause k = 2..K (index k - 2): columns used
  bool empty() const {
    for (const auto& c : wcol) if (!c.empty()) return false;
    for (const auto& c : vcol) if (!c.empty()) return false;
    return true;
  }
};

class ConstrainedKModel : public DirectCIFModel {
 public:
  // cause: 0 = censored, 1 = modelled directly, 2..K = competing causes; family[k-2] for G_k.
  ConstrainedKModel(std::size_t n, const double* time, const int* cause, const double* X, const DirectSpec& spec,
                    int K, const std::vector<int>& family, const CondCovariates& cov = CondCovariates())
      : DirectCIFModel(n, time, cause, X, spec), K_(K), fam_(family), P1_(spec.nparam()), cov_(cov) {
    const std::size_t Kc = static_cast<std::size_t>(K_ - 1);
    PG_ = P1_ + 2 * Kc;
    PT_ = PG_ + static_cast<std::size_t>(K_ - 2) * (1 + spec.p);
    cov_.wcol.resize(Kc);
    cov_.vcol.resize(Kc);
    has_cov_ = !cov_.empty();
    epos_.assign(Kc, 0); zpos_.assign(Kc, 0);
    std::size_t o = PT_;
    for (std::size_t k = 0; k < Kc; ++k) { epos_[k] = o; o += cov_.wcol[k].size(); }
    for (std::size_t k = 0; k < Kc; ++k) { zpos_[k] = o; o += cov_.vcol[k].size(); }
    PA_ = o;
  }
  std::size_t nparam() const { return PA_; }
  std::size_t nparam_base() const { return PT_; }          // without the covariates of the G_k
  bool has_cond_covariates() const { return has_cov_; }
  const double* data_X() const { return X_; }
  const CondCovariates& cond_covariates() const { return cov_; }
  int K() const { return K_; }
  const std::vector<int>& families() const { return fam_; }
  std::size_t gpos(int k) const { return P1_ + 2 * static_cast<std::size_t>(k - 2); }               // a_k, b_k
  std::size_t lpos(int k) const { return PG_ + static_cast<std::size_t>(k - 3) * (1 + s_.p); }     // gamma_k, delta_k
  std::size_t epos(int k) const { return epos_[static_cast<std::size_t>(k - 2)]; }                // eta_k
  std::size_t zpos(int k) const { return zpos_[static_cast<std::size_t>(k - 2)]; }                // zeta_k
  bool feasible_in_range(const std::vector<double>&) const { return true; }

  // G_k(t) (and log g_k(t) when logg is given) for the covariate rows w (length nw) and v (length
  // nv) of one subject, as D2 in the local variables (a_k, b_k, u, s); terms (optional) receive
  // the map from the local variables to the parameters. Without covariates for G_k this is
  // exactly cond_Gg in (a_k, b_k).
  void cond_local(const std::vector<double>& th, int k, double t, const double* w, const double* v, D2& G, D2* logg,
                  std::vector<Term>* terms) const {
    const std::size_t o = gpos(k), kk = static_cast<std::size_t>(k - 2);
    const std::vector<std::size_t>& wc = cov_.wcol[kk];
    const std::vector<std::size_t>& vc = cov_.vcol[kk];
    D2 a = D2::variable(0, th[o]), b = D2::variable(1, th[o + 1]);
    if (!wc.empty() || !vc.empty()) {
      double uv = 0.0, sv = 0.0;
      for (std::size_t j = 0; j < wc.size(); ++j) uv += th[epos_[kk] + j] * w[wc[j]];
      for (std::size_t j = 0; j < vc.size(); ++j) sv += th[zpos_[kk] + j] * v[vc[j]];
      const D2 u = D2::variable(2, uv), s = D2::variable(3, sv);
      switch (fam_[kk]) {
        case COND_LOGISTIC: a = a - u + s; b = b * exp(u); break;
        case COND_LOGNORMAL: a = a + u; b = b + s; break;
        default: a = a - u; b = b + s; break;
      }
    }
    cond_Gg<D2>(fam_[kk], t, a, b, G, logg);
    if (terms) {
      terms->clear();
      terms->push_back({o, 0, 1.0});
      terms->push_back({o + 1, 1, 1.0});
      for (std::size_t j = 0; j < wc.size(); ++j) terms->push_back({epos_[kk] + j, 2, w[wc[j]]});
      for (std::size_t j = 0; j < vc.size(); ++j) terms->push_back({zpos_[kk] + j, 3, v[vc[j]]});
    }
  }
  // covariate rows of subject i (null without covariates)
  const double* wrow(std::size_t i) const { return cov_.W ? cov_.W + i * cov_.nw : nullptr; }
  const double* vrow(std::size_t i) const { return cov_.V ? cov_.V + i * cov_.nv : nullptr; }

  // shares pi_k, k = 2..K (index k - 2)
  void shares(const std::vector<double>& th, const double* x, std::vector<double>& pi) const {
    pi.assign(static_cast<std::size_t>(K_ - 1), 0.0);
    std::vector<double> l(static_cast<std::size_t>(K_ - 1), 0.0);
    double mx = 0.0;
    for (int k = 3; k <= K_; ++k) {
      const std::size_t o = lpos(k);
      double v = th[o];
      for (std::size_t r = 0; r < s_.p; ++r) v += th[o + 1 + r] * x[r];
      l[k - 2] = v; mx = std::fmax(mx, v);
    }
    double sum = 0.0;
    for (std::size_t j = 0; j < l.size(); ++j) { pi[j] = std::exp(l[j] - mx); sum += pi[j]; }
    for (double& v : pi) v /= sum;
  }

  // CIF of cause k at time t for the covariate row x (and the rows w, v of the covariates of the
  // G_k, needed only when the model has them).
  double cifk(const std::vector<double>& th, int k, const double* x, double t, const double* w = nullptr,
              const double* v = nullptr) const {
    if (t <= 0.0) return 0.0;
    if (k == 1) return cif(th, 0, x, t);
    const double q = 1.0 - cif(th, 0, x, std::numeric_limits<double>::infinity());
    std::vector<double> pi;
    shares(th, x, pi);
    double G = 1.0;
    if (has_cov_) { D2 Gd; cond_local(th, k, t, w, v, Gd, nullptr, nullptr); G = Gd.v; }
    else cond_Gg<double>(fam_[k - 2], t, th[gpos(k)], th[gpos(k) + 1], G, nullptr);
    return q * pi[k - 2] * G;
  }
  double cifk_grad(const std::vector<double>& th, int k, const double* x, double t, std::vector<double>& grad,
                   const double* w = nullptr, const double* v = nullptr) const {
    grad.assign(PA_, 0.0);
    if (t <= 0.0) return 0.0;
    std::vector<double> g1;
    if (k == 1) {
      const double F = cif_grad(th, 0, x, t, g1);
      for (std::size_t j = 0; j < P1_; ++j) grad[j] = g1[j];
      return F;
    }
    const double Finf = cif_grad(th, 0, x, std::numeric_limits<double>::infinity(), g1);
    const double q = 1.0 - Finf;
    std::vector<double> pi;
    shares(th, x, pi);
    D2 G;
    std::vector<Term> tg;
    cond_local(th, k, t, w, v, G, nullptr, &tg);
    const double pk = pi[k - 2];
    for (std::size_t j = 0; j < P1_; ++j) grad[j] = -g1[j] * pk * G.v;
    for (const Term& e : tg) grad[e.idx] += q * pk * G.g[e.loc] * e.coef;
    for (int j = 3; j <= K_; ++j) {                          // d pi_k / d l_j = pi_k (1{k=j} - pi_j)
      const double dl = q * G.v * pk * ((k == j ? 1.0 : 0.0) - pi[j - 2]);
      const std::size_t oj = lpos(j);
      grad[oj] += dl;
      for (std::size_t r = 0; r < s_.p; ++r) grad[oj + 1 + r] += dl * x[r];
    }
    return q * pk * G.v;
  }

  void eval(const std::vector<double>& th, Eval& out, bool need_hess) const {
    const std::size_t P = PA_, p = s_.p, Kc = static_cast<std::size_t>(K_ - 1);
    const double inf = std::numeric_limits<double>::infinity();
    out.ll = 0.0;
    out.g.assign(P, 0.0);
    if (out.H.n != P) out.H = Mat(P);
    out.H.zero();
    std::vector<double> z, dz, gS(P), pi, HS(need_hess ? P * P : 0);
    std::vector<Term> t1;
    std::vector<std::vector<Term>> tg(Kc);                   // local variables of G_k -> parameters
    D2 q[3], lf, Ft, Finf;
    std::vector<D2> G(Kc), logg(Kc);
    // d l_j / d theta: l_j = gamma_j + delta_j' x  (j = 3..K)
    auto lidx = [&](int j, std::size_t r) { return lpos(j) + r; };   // r = 0: gamma, r >= 1: delta_{r-1}
    for (std::size_t i = 0; i < n_; ++i) {
      const double* x = X_ + i * p;
      const double* wi = wrow(i);
      const double* vi = vrow(i);
      const double t = time_[i];
      const int c = cause_[i];
      const D2 al = D2::variable(3, alpha_of(th, 0)), et = D2::variable(4, eta_of(th, 0, x));
      locals<D2>(th, 0, t, nullptr, nullptr, q, z, dz);
      make_terms(0, x, z.data(), dz.data(), t1);
      auto lcoef = [&](std::size_t r) { return r == 0 ? 1.0 : x[r - 1]; };
      if (c == 1) {
        if (!direct_Ff<D2>(s_.baseline[0], t, q, al, et, nullptr, &lf) || !std::isfinite(lf.v)) { bad(out); return; }
        out.ll += lf.v;
        add_terms(t1, lf, 1.0, out, need_hess, nullptr);
        continue;
      }
      if (c == 0 && !(t > 0.0)) continue;                  // censored at time 0
      if (!direct_Ff<D2>(s_.baseline[0], inf, q, al, et, &Finf, nullptr)) { bad(out); return; }
      const double qv = 1.0 - Finf.v;
      if (!(qv > 0.0)) { bad(out); return; }
      shares(th, x, pi);
      if (c >= 2) {                                         // log(1 - p_1) + log pi_c + log g_c(t)
        const int k = c;
        cond_local(th, k, t, wi, vi, G[0], &logg[0], &tg[0]);
        if (!std::isfinite(logg[0].v) || !(pi[k - 2] > 0.0)) { bad(out); return; }
        const D2 lq = chain(Finf, std::log(qv), -1.0 / qv, -1.0 / (qv * qv));
        out.ll += lq.v + std::log(pi[k - 2]) + logg[0].v;
        add_terms(t1, lq, 1.0, out, need_hess, nullptr);
        add_terms(tg[0], logg[0], 1.0, out, need_hess, nullptr);
        for (int j = 3; j <= K_; ++j) {                     // d log pi_k / d l_j = 1{k=j} - pi_j
          const double dj = (k == j ? 1.0 : 0.0) - pi[j - 2];
          for (std::size_t r = 0; r <= p; ++r) out.g[lidx(j, r)] += dj * lcoef(r);
          if (!need_hess) continue;
          for (int m = 3; m <= j; ++m) {                    // d2 / dl_j dl_m = -pi_j (1{j=m} - pi_m)
            const double hjm = -pi[j - 2] * ((j == m ? 1.0 : 0.0) - pi[m - 2]);
            for (std::size_t r = 0; r <= p; ++r)
              for (std::size_t u = 0; u <= p; ++u) {
                const std::size_t a = lidx(j, r), b = lidx(m, u);
                if (b <= a) out.H(a, b) += hjm * lcoef(r) * lcoef(u);
              }
          }
        }
        continue;
      }
      // censored: S = 1 - F_1(t) - q M,  M = sum_k pi_k G_k(t)
      if (!direct_Ff<D2>(s_.baseline[0], t, q, al, et, &Ft, nullptr)) { bad(out); return; }
      double M = 0.0;
      for (int k = 2; k <= K_; ++k) {
        cond_local(th, k, t, wi, vi, G[k - 2], nullptr, &tg[k - 2]);
        M += pi[k - 2] * G[k - 2].v;
      }
      const double S = 1.0 - Ft.v - qv * M;
      if (!(S > 0.0)) { bad(out); return; }
      out.ll += std::log(S);
      // gradient of S
      std::fill(gS.begin(), gS.end(), 0.0);
      for (const Term& e : t1) gS[e.idx] += (-Ft.g[e.loc] + M * Finf.g[e.loc]) * e.coef;
      std::vector<double> dMl(Kc, 0.0);                      // dM / dl_j (index j - 2)
      for (int j = 3; j <= K_; ++j) {
        double s = 0.0;
        for (int k = 2; k <= K_; ++k) s += G[k - 2].v * pi[k - 2] * ((k == j ? 1.0 : 0.0) - pi[j - 2]);
        dMl[j - 2] = s;
      }
      for (int k = 2; k <= K_; ++k)
        for (const Term& e : tg[k - 2]) gS[e.idx] += -qv * pi[k - 2] * G[k - 2].g[e.loc] * e.coef;
      for (int j = 3; j <= K_; ++j) for (std::size_t r = 0; r <= p; ++r) gS[lidx(j, r)] += -qv * dMl[j - 2] * lcoef(r);
      for (std::size_t a = 0; a < P; ++a) out.g[a] += gS[a] / S;
      if (!need_hess) continue;
      std::fill(HS.begin(), HS.end(), 0.0);
      auto H = [&](std::size_t a, std::size_t b) -> double& { return a >= b ? HS[a * P + b] : HS[b * P + a]; };
      // block 1: -H_Ft + M H_Finf
      for (const Term& e : t1) for (const Term& f : t1)
        if (f.idx <= e.idx) H(e.idx, f.idx) += e.coef * f.coef * (-Ft.h[e.loc][f.loc] + M * Finf.h[e.loc][f.loc]);
      // dFinf as a parameter vector
      std::vector<double> dFi(P1_, 0.0);
      for (const Term& e : t1) dFi[e.idx] += Finf.g[e.loc] * e.coef;
      for (int k = 2; k <= K_; ++k) {
        const double pk = pi[k - 2];
        const D2& Gk = G[k - 2];
        const std::vector<Term>& tk = tg[k - 2];
        // G_k block: -q pi_k H_Gk (parameters of G_k and of its covariates; lower triangle)
        for (const Term& e : tk)
          for (const Term& f : tk)
            if (f.idx <= e.idx) H(e.idx, f.idx) += -qv * pk * Gk.h[e.loc][f.loc] * e.coef * f.coef;
        // block 1 x G_k: dFinf pi_k dG_k
        for (const Term& e : tk)
          for (std::size_t j = 0; j < P1_; ++j) H(e.idx, j) += dFi[j] * pk * Gk.g[e.loc] * e.coef;
        // G_k x l_j: -q dG_k dpi_k/dl_j
        for (int j = 3; j <= K_; ++j) {
          const double dp = pk * ((k == j ? 1.0 : 0.0) - pi[j - 2]);
          for (std::size_t r = 0; r <= p; ++r)
            for (const Term& e : tk) H(lidx(j, r), e.idx) += -qv * Gk.g[e.loc] * e.coef * dp * lcoef(r);
        }
      }
      for (int j = 3; j <= K_; ++j) {
        // block 1 x l_j: dFinf dM/dl_j
        for (std::size_t r = 0; r <= p; ++r) for (std::size_t a = 0; a < P1_; ++a) H(lidx(j, r), a) += dFi[a] * dMl[j - 2] * lcoef(r);
        // l_j x l_m: -q sum_k G_k d2 pi_k / dl_j dl_m
        for (int m = 3; m <= j; ++m) {
          double s = 0.0;
          for (int k = 2; k <= K_; ++k) {
            const double ej = (k == j ? 1.0 : 0.0) - pi[j - 2], em = (k == m ? 1.0 : 0.0) - pi[m - 2];
            s += G[k - 2].v * pi[k - 2] * (ej * em - pi[j - 2] * ((j == m ? 1.0 : 0.0) - pi[m - 2]));
          }
          for (std::size_t r = 0; r <= p; ++r)
            for (std::size_t u = 0; u <= p; ++u) {
              const std::size_t a = lidx(j, r), b = lidx(m, u);
              if (b <= a) H(a, b) += -qv * s * lcoef(r) * lcoef(u);
              else if (j != m) H(b, a) += -qv * s * lcoef(r) * lcoef(u);
            }
        }
      }
      // log S: H_S / S - gS gS' / S^2  (lower triangle)
      for (std::size_t a = 0; a < P; ++a)
        for (std::size_t b = 0; b <= a; ++b) out.H(a, b) += HS[a * P + b] / S - gS[a] * gS[b] / (S * S);
    }
    if (!std::isfinite(out.ll)) { bad(out); return; }
    if (need_hess)
      for (std::size_t r = 0; r < P; ++r)
        for (std::size_t cc = 0; cc <= r; ++cc) { out.H(r, cc) = -out.H(r, cc); out.H(cc, r) = out.H(r, cc); }
  }

 private:
  int K_;
  std::vector<int> fam_;
  std::size_t P1_, PG_ = 0, PT_ = 0, PA_ = 0;
  CondCovariates cov_;
  bool has_cov_ = false;
  std::vector<std::size_t> epos_, zpos_;     // first parameter of eta_k, zeta_k (index k - 2)
};

inline std::vector<int> constrainedK_roles(const DirectSpec& s, int K) {
  std::vector<int> role = direct_roles(s);
  for (int k = 2; k <= K; ++k) { role.push_back(0); role.push_back(0); }
  for (int k = 3; k <= K; ++k) { role.push_back(0); for (std::size_t r = 0; r < s.p; ++r) role.push_back(1); }
  return role;
}

// Start values: cause 1 block and every G_k as in constrained_start (cause k recoded to 2); the
// logits from the observed shares of the competing events.
inline std::vector<double> constrainedK_start(std::size_t n, const double* time, const int* cause, const DirectSpec& s,
                                              int K, const std::vector<int>& family) {
  std::vector<int> c2(cause, cause + n);
  for (int& v : c2) if (v >= 2) v = 2;
  std::vector<double> th = constrained_start(n, time, c2.data(), s, family[0]);
  th.resize(s.nparam());
  std::vector<double> cnt(static_cast<std::size_t>(K + 1), 0.0);
  for (std::size_t i = 0; i < n; ++i) if (cause[i] >= 1 && cause[i] <= K) cnt[static_cast<std::size_t>(cause[i])] += 1.0;
  for (int k = 2; k <= K; ++k) {
    std::vector<int> ck(n);
    for (std::size_t i = 0; i < n; ++i) ck[i] = cause[i] == 1 ? 1 : (cause[i] == k ? 2 : 0);
    const std::vector<double> tk = constrained_start(n, time, ck.data(), s, family[k - 2]);
    th.push_back(tk[tk.size() - 2]); th.push_back(tk[tk.size() - 1]);
  }
  for (int k = 3; k <= K; ++k) {
    th.push_back(std::log(std::fmax(cnt[k], 0.5) / std::fmax(cnt[2], 0.5)));
    for (std::size_t r = 0; r < s.p; ++r) th.push_back(0.0);
  }
  return th;
}

// Without covariates of the G_k: the staged fit. With them (and no start values) the fit is
// nested: first the model without them (exactly the fit above), then all parameters from that
// maximum with the new coefficients at zero. Every accepted step of the final stage increases
// the likelihood, so the fit is never worse than the model without these covariates.
inline FitResult constrainedK_fit(const ConstrainedKModel& m, std::size_t n, const double* time, const int* cause,
                                  const std::vector<int>& family, const std::vector<double>* init, int max_iter, double tol) {
  const DirectSpec& s = m.spec();
  auto shrink = [&s](std::vector<double>& v) { v[base_level_index(s.baseline[0])] -= 0.5; };
  const std::vector<int> role0 = constrainedK_roles(s, m.K());
  if (!m.has_cond_covariates()) {
    return staged_fit(m, init ? *init : constrainedK_start(n, time, cause, s, m.K(), family), role0,
                      init != nullptr, false, max_iter, tol, shrink);
  }
  std::vector<int> role = role0;
  role.resize(m.nparam(), 1);
  if (init) return staged_fit(m, *init, role, true, false, max_iter, tol, shrink);
  const ConstrainedKModel m0(n, time, cause, m.data_X(), s, m.K(), family);
  const FitResult f0 = staged_fit(m0, constrainedK_start(n, time, cause, s, m.K(), family), role0, false, false, max_iter, tol, shrink);
  if (f0.flags & FIT_NONFINITE) { FitResult r = f0; r.beta.resize(m.nparam(), 0.0); return r; }
  std::vector<double> th = f0.beta;
  th.resize(m.nparam(), 0.0);
  FitResult f = staged_fit(m, th, role, true, false, max_iter, tol, shrink);
  f.ll_init = f0.ll_init;
  f.iterations += f0.iterations;
  return f;
}

// The parameters of G_k fitted with centred covariates (W - mw, V - mv) expressed for the
// covariates themselves, and the Jacobian of that map (P x P, row-major). The map is the identity
// except for (a_k, b_k): with c = eta_k' mw and d = zeta_k' mv,
//   Weibull, log-logistic: a = a^c + c,  b = b^c - d;   log-normal: a = a^c - c,  b = b^c - d;
//   modified logistic:     a = a^c + c - d,  b = b^c e^{-c}.
// The likelihood is the same function of either set, so the maximum and the inverse information
// transform exactly (the Hessian term with the gradient vanishes at a stationary point).
inline void cond_uncentre(const ConstrainedKModel& m, const std::vector<double>& mw, const std::vector<double>& mv,
                          std::vector<double>& th, std::vector<double>& J) {
  const std::size_t P = m.nparam();
  J.assign(P * P, 0.0);
  for (std::size_t a = 0; a < P; ++a) J[a * P + a] = 1.0;
  const CondCovariates& cv = m.cond_covariates();
  const std::vector<double> th0 = th;
  for (int k = 2; k <= m.K(); ++k) {
    const std::size_t kk = static_cast<std::size_t>(k - 2), o = m.gpos(k), e = m.epos(k), z = m.zpos(k);
    const std::vector<std::size_t>& wc = cv.wcol[kk];
    const std::vector<std::size_t>& vc = cv.vcol[kk];
    if (wc.empty() && vc.empty()) continue;
    double c = 0.0, d = 0.0;
    for (std::size_t j = 0; j < wc.size(); ++j) c += th0[e + j] * mw[wc[j]];
    for (std::size_t j = 0; j < vc.size(); ++j) d += th0[z + j] * mv[vc[j]];
    const int fam = m.families()[kk];
    if (fam == COND_LOGISTIC) {
      const double ec = std::exp(-c);
      th[o] = th0[o] + c - d;
      th[o + 1] = th0[o + 1] * ec;
      for (std::size_t j = 0; j < wc.size(); ++j) { J[o * P + e + j] = mw[wc[j]]; J[(o + 1) * P + e + j] = -th0[o + 1] * ec * mw[wc[j]]; }
      for (std::size_t j = 0; j < vc.size(); ++j) J[o * P + z + j] = -mv[vc[j]];
      J[(o + 1) * P + o + 1] = ec;
    } else {
      const double sg = fam == COND_LOGNORMAL ? -1.0 : 1.0;
      th[o] = th0[o] + sg * c;
      th[o + 1] = th0[o + 1] - d;
      for (std::size_t j = 0; j < wc.size(); ++j) J[o * P + e + j] = sg * mw[wc[j]];
      for (std::size_t j = 0; j < vc.size(); ++j) J[(o + 1) * P + z + j] = -mv[vc[j]];
    }
  }
}

}  // namespace survcore
