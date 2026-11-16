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

}  // namespace survcore
