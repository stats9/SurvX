// Memory and undefined-behaviour checks of the whole C API of survcore (every sc_* entry point),
// on random data of many sizes with ties, zero times, one covariate column or none, and invalid
// input. Run under AddressSanitizer + UBSan and under valgrind. Output buffers are sized as the
// R wrappers size them.
#include "survcore/survcore.h"
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using V = std::vector<double>;
using I = std::vector<int>;
static std::mt19937_64 rng(11);
static double U() { return std::uniform_real_distribution<double>(0, 1)(rng); }
static double N() { return std::normal_distribution<double>(0, 1)(rng); }
static int calls = 0, nonzero = 0;
static void tally(int rc) { ++calls; nonzero += rc != 0; }

struct Data { int n, p, K; V time, entry, X; I status, cause, group, strata, cluster; };

static Data make(int n, int p, int K, bool ties) {
  Data d; d.n = n; d.p = p; d.K = K;
  d.time.resize(n); d.entry.resize(n); d.X.resize(static_cast<size_t>(n) * (p ? p : 1)); d.status.resize(n); d.cause.resize(n);
  d.group.resize(n); d.strata.resize(n); d.cluster.resize(n);
  for (int i = 0; i < n; ++i) {
    double lp = 0;
    for (int j = 0; j < p; ++j) { const double x = j == 0 ? (U() < 0.5) : N(); d.X[static_cast<size_t>(i) * p + j] = x; lp += 0.3 * x; }
    double best = 1e300; int cb = 0;
    for (int k = 1; k <= K; ++k) { const double e = -std::log(U()) / (0.3 * k * std::exp(lp)); if (e < best) { best = e; cb = k; } }
    const double c = 4.0 * U();
    double t = std::fmin(best, c) + 1e-3;
    if (ties) t = std::ceil(t * 5) / 5;
    d.time[i] = t; d.cause[i] = best < c ? cb : 0; d.status[i] = d.cause[i] > 0;
    d.entry[i] = U() < 0.3 ? 0.5 * t * U() : 0.0;
    d.group[i] = i % 3; d.strata[i] = i % 2; d.cluster[i] = i / 4;
  }
  for (int k = 0; k <= K && k < n; ++k) { d.cause[k] = k; d.status[k] = k > 0; }   // every cause present
  return d;
}

int main() {
  sc_fit_info info;
  #ifndef REPS
#define REPS 36
#endif
  for (int rep = 0; rep < REPS; ++rep) {
    const int n = 30 + (rep * 37) % 260, p = rep % 3, K = 1 + rep % 3;
    Data d = make(n, p, K, rep % 4 == 0);
    const double* X = p ? d.X.data() : nullptr;
    const int pp = p ? p : 1;
    // ---- non-parametric
    { V ot(n), nr(n), ne(n), s(n), se(n); int m = 0;
      tally(sc_kaplan_meier(n, d.time.data(), d.status.data(), nullptr, ot.data(), nr.data(), ne.data(), s.data(), se.data(), &m));
      V w(n, 1.5); tally(sc_kaplan_meier(n, d.time.data(), d.status.data(), w.data(), ot.data(), nr.data(), ne.data(), s.data(), se.data(), &m));
      V ch(n); tally(sc_nelson_aalen(n, d.time.data(), d.status.data(), nullptr, ot.data(), ch.data(), se.data(), &m));
      double mean, sem; tally(sc_rmst(n, d.time.data(), d.status.data(), 2.0, &mean, &sem));
      V sv(n), cf(static_cast<size_t>(n) * K), cse(static_cast<size_t>(n) * K);
      tally(sc_aalen_johansen(n, d.time.data(), d.cause.data(), K, ot.data(), sv.data(), cf.data(), cse.data(), &m));
      double st, pv; int df; V obs(3), ex(3), var(9);
      tally(sc_logrank(n, d.time.data(), d.status.data(), d.group.data(), 3, nullptr, 0.0, &st, &df, &pv, obs.data(), ex.data(), var.data()));
      tally(sc_logrank(n, d.time.data(), d.status.data(), d.group.data(), 3, d.strata.data(), 1.0, &st, &df, &pv, obs.data(), ex.data(), var.data()));
      V sc(3);
      tally(sc_gray_test(n, d.time.data(), d.cause.data(), d.group.data(), 3, nullptr, 1, 0.0, &st, &df, &pv, sc.data(), var.data()));
      tally(sc_gray_test(n, d.time.data(), d.cause.data(), d.group.data(), 3, d.strata.data(), 1, 0.5, &st, &df, &pv, sc.data(), var.data()));
      (void)sc_chisq_upper(3.2, 2.0);
    }
    // ---- Cox and Fine-Gray
    if (p) {
      V b(p, 0.0), v(static_cast<size_t>(p) * p), g(p), H(static_cast<size_t>(p) * p), ot(n), ch(n);
      int m = 0; double ll;
      for (int ties = 0; ties < 2; ++ties) {
        std::fill(b.begin(), b.end(), 0.0);
        tally(sc_coxph(n, p, d.time.data(), d.status.data(), X, nullptr, ties, 50, 1e-9, b.data(), v.data(), &info));
        tally(sc_coxph_eval(n, p, d.time.data(), d.status.data(), X, nullptr, ties, b.data(), &ll, g.data(), H.data()));
        tally(sc_coxph_baseline(n, p, d.time.data(), d.status.data(), X, nullptr, ties, b.data(), nullptr, ot.data(), ch.data(), &m));
        V b2(p, 0.0); I os(n);
        tally(sc_coxph_ext(n, p, d.entry.data(), d.time.data(), d.status.data(), d.strata.data(), X, nullptr, ties, 50, 1e-9, b2.data(), v.data(), &info));
        tally(sc_coxph_ext_baseline(n, p, d.entry.data(), d.time.data(), d.status.data(), d.strata.data(), X, nullptr, ties, b2.data(), nullptr,
                                    ot.data(), ch.data(), os.data(), &m));
        V mart(n), xbar(static_cast<size_t>(n) * p), Us(static_cast<size_t>(n) * p);
        tally(sc_coxph_resid(n, p, d.entry.data(), d.time.data(), d.status.data(), d.strata.data(), X, nullptr, ties, b2.data(), mart.data(), xbar.data()));
        tally(sc_coxph_score(n, p, nullptr, d.time.data(), d.status.data(), nullptr, X, nullptr, ties, b2.data(), Us.data()));
      }
      tally(sc_coxph_baseline(n, p, d.time.data(), d.status.data(), X, nullptr, 1, b.data(), nullptr, ot.data(), ch.data(), &m));
      V lp(5, 0.2), out(static_cast<size_t>(5) * m);
      tally(sc_predict_from_cumhaz(m, ch.data(), 5, lp.data(), 0, out.data()));
      tally(sc_predict_from_cumhaz(m, ch.data(), 5, lp.data(), 1, out.data()));
      V bf(p, 0.0), vm(static_cast<size_t>(p) * p), vr(static_cast<size_t>(p) * p);
      tally(sc_finegray(n, p, d.time.data(), d.cause.data(), X, 1, 50, 1e-9, bf.data(), vm.data(), vr.data(), &info));
      tally(sc_finegray_eval(n, p, d.time.data(), d.cause.data(), X, 1, bf.data(), &ll, g.data(), H.data()));
      tally(sc_finegray_baseline(n, p, d.time.data(), d.cause.data(), X, 1, bf.data(), nullptr, ot.data(), ch.data(), &m));
      V betas(static_cast<size_t>(K) * p, 0.1); const int nn = 3;
      V cif(static_cast<size_t>(nn) * n * K);
      for (int pl = 0; pl < 2; ++pl)
        tally(sc_csc_cif(n, p, d.time.data(), d.cause.data(), X, K, betas.data(), pl, nn, X, pl, ot.data(), cif.data(), &m));
      V tm = {0.5, 1.0, 2.0};
      for (int model = 0; model < 2; ++model)
        for (int avg = 0; avg < 2; ++avg) {
          const int cells = avg ? 3 : nn * 3; V e(cells), s(cells), lo(cells), hi(cells); int ok = 0;
          tally(sc_cif_bootstrap(model, n, p, d.time.data(), d.cause.data(), X, K, 1, 1, 1, 1, nn, X, avg, 3, tm.data(), 5, 7, 0.95,
                                 e.data(), s.data(), lo.data(), hi.data(), &ok));
        }
    }
    // ---- prediction metrics
    { V pred(n); for (int i = 0; i < n; ++i) pred[i] = U();
      double br, bn, ipa, auc, err; int m = 0; V mp(10), ob(10), sz(10);
      tally(sc_brier(n, d.time.data(), d.cause.data(), 1, pred.data(), 1.0, &br, &bn, &ipa));
      tally(sc_auc(n, d.time.data(), d.cause.data(), 1, pred.data(), 1.0, &auc));
      tally(sc_calibration(n, d.time.data(), d.cause.data(), 1, pred.data(), 1.0, 10, mp.data(), ob.data(), sz.data(), &err, &m));
      I fold(n); tally(sc_cv_folds(n, 5, 3, fold.data()));
    }
    // ---- direct CIF models (joint, interval-censored, weighted, constrained)
    for (int base = 0; base < 8; ++base) {
      I af(K, rep % 2), bs(K, base), nk(K, 0); V av(K, rep % 2 ? 0.5 : 0.0), kn;
      if (base == 7) for (int k = 1; k <= K; ++k) { V kk(4); int m = 0; if (sc_directcif_default_knots(n, d.time.data(), d.cause.data(), k, 3, kk.data(), &m)) goto next_base; nk[k - 1] = m; kn.insert(kn.end(), kk.begin(), kk.begin() + m); }
      {
        const int P = sc_directcif_nparam(p, K, af.data(), bs.data(), nk.data());
        V th(P), var(static_cast<size_t>(P) * P), dg(5), g(P), H(static_cast<size_t>(P) * P); double ll;
        const double* knp = kn.empty() ? nullptr : kn.data();
        tally(sc_directcif_fit(n, p, K, d.time.data(), d.cause.data(), X, af.data(), av.data(), bs.data(), nk.data(), knp, rep % 3 == 0, nullptr, 100, 1e-8,
                               th.data(), var.data(), &info, dg.data()));
        tally(sc_directcif_eval(n, p, K, d.time.data(), d.cause.data(), X, af.data(), av.data(), bs.data(), nk.data(), knp, th.data(), &ll, g.data(), H.data()));
        V tm = {0.0, 0.7, 2.0, INFINITY}; const int nn = 4; V out(nn * 4), se(nn * 4), om(4), sm(4);
        for (int k = 1; k <= K; ++k) {
          tally(sc_directcif_predict(p, K, af.data(), av.data(), bs.data(), nk.data(), knp, th.data(), var.data(), k, nn, X, 4, tm.data(), out.data(), se.data()));
          tally(sc_directcif_predict_mean(p, K, af.data(), av.data(), bs.data(), nk.data(), knp, th.data(), var.data(), k, nn, X, 4, tm.data(), om.data(), sm.data()));
        }
        // interval-censored version: every event in (0.8 t, t]
        V t1(n), t2(n);
        for (int i = 0; i < n; ++i) { t1[i] = d.cause[i] ? 0.8 * d.time[i] : d.time[i]; t2[i] = d.time[i]; if (i % 7 == 0) t1[i] = t2[i]; }
        tally(sc_directcif_fit_ic(n, p, K, t1.data(), t2.data(), d.cause.data(), X, af.data(), av.data(), bs.data(), nk.data(), knp, 0, th.data(), 50, 1e-8,
                                  th.data(), var.data(), &info, dg.data()));
        tally(sc_directcif_eval_ic(n, p, K, t1.data(), t2.data(), d.cause.data(), X, af.data(), av.data(), bs.data(), nk.data(), knp, th.data(), &ll, g.data(), H.data()));
        // weighted, one cause
        const int nkb = base == 7 ? nk[0] : 0; const int P1 = sc_directcif_nparam(p, 1, af.data(), bs.data(), nk.data());
        V tw(P1), vr(static_cast<size_t>(P1) * P1), vmm(static_cast<size_t>(P1) * P1), vip(static_cast<size_t>(P1) * P1), gw(P1), Hw(static_cast<size_t>(P1) * P1); int ns = 0;
        for (int grid = 0; grid <= 10; grid += 10) {
          tally(sc_subdist_fit(n, p, d.time.data(), d.cause.data(), X, 1, af[0], av[0], base, nkb, knp, grid, nullptr, 100, 1e-8, tw.data(), vr.data(), vmm.data(), &info, &ns));
          tally(sc_subdist_eval(n, p, d.time.data(), d.cause.data(), X, 1, af[0], av[0], base, nkb, knp, grid, tw.data(), &ll, gw.data(), Hw.data()));
        }
        tally(sc_subdist_var_ipcw(n, p, d.time.data(), d.cause.data(), X, 1, af[0], av[0], base, nkb, knp, tw.data(), vmm.data(), vip.data()));
        // constrained (needs a baseline with a plateau and two or more causes)
        if (K >= 2 && base != 2 && base != 7)
          for (int fam = 0; fam < 4; ++fam) {
            const int Pc = P1 + 2; V tc(Pc), vc(static_cast<size_t>(Pc) * Pc), gc(Pc), Hc(static_cast<size_t>(Pc) * Pc);
            I c2(d.cause); for (int& c : c2) if (c > 2) c = 2;
            tally(sc_constrained_fit(n, p, d.time.data(), c2.data(), X, af[0], av[0], base, fam, nullptr, 100, 1e-8, tc.data(), vc.data(), &info));
            tally(sc_constrained_eval(n, p, d.time.data(), c2.data(), X, af[0], av[0], base, fam, tc.data(), &ll, gc.data(), Hc.data()));
            for (int k = 1; k <= 2; ++k) {
              tally(sc_constrained_predict(p, af[0], av[0], base, fam, tc.data(), vc.data(), k, nn, X, 4, tm.data(), out.data(), se.data()));
              tally(sc_constrained_predict_mean(p, af[0], av[0], base, fam, tc.data(), vc.data(), k, nn, X, 4, tm.data(), om.data(), sm.data()));
            }
            I famK(K - 1, fam); const int PK = P1 + 2 * (K - 1) + (K - 2) * (1 + p);
            V tk(PK), vk(static_cast<size_t>(PK) * PK), gk(PK), Hk(static_cast<size_t>(PK) * PK);
            tally(sc_constrainedK_fit(n, p, d.time.data(), d.cause.data(), X, K, af[0], av[0], base, famK.data(), nullptr, 100, 1e-8, tk.data(), vk.data(), &info));
            tally(sc_constrainedK_eval(n, p, d.time.data(), d.cause.data(), X, K, af[0], av[0], base, famK.data(), tk.data(), &ll, gk.data(), Hk.data()));
            for (int k = 1; k <= K; ++k) for (int a = 0; a < 2; ++a)
              tally(sc_constrainedK_predict(p, K, af[0], av[0], base, famK.data(), tk.data(), vk.data(), k, a, nn, X, 4, tm.data(), a ? om.data() : out.data(), a ? sm.data() : se.data()));
          }
      }
      next_base:;
    }
    // ---- model comparison (direct models, Fine-Gray, cause-specific Cox), cross-validated
    if (p && K >= 2 && rep % 3 == 0) {
      I bs = {0, 3, 7}, lk = {0, 1, 0}, sdf = {3, 3, 2}, est = {0, 1, 2}, cnd = {0, 1, 1};
      V tm = {0.5, 1.5}; const int M = 3 + 2;
      V met(static_cast<size_t>(M) * 2 * 4 * 2), fit(static_cast<size_t>(M) * 5);
      I c2(d.cause); for (int& c : c2) if (c > 2) c = 2;
      tally(sc_compare_models(n, p, 2, d.time.data(), c2.data(), X, 1, 2, tm.data(), 3, 5, 3, bs.data(), lk.data(), sdf.data(), met.data(), fit.data()));
      I bs2 = {0, 0, 3}; I est2 = {0, 1, 2};
      tally(sc_compare_models2(n, p, 2, d.time.data(), c2.data(), X, 1, 2, tm.data(), 3, 5, 3, bs2.data(), lk.data(), sdf.data(), est2.data(), cnd.data(),
                               met.data(), fit.data()));
    }
    // ---- parametric models, interval censoring, cure, frailty, cause-specific parametric
    {
      const int dists[] = {0, 1, 2, 3, 4, 5, 6, 7};
      for (int di : dists) for (int aft = 0; aft < 2; ++aft) {
        if (aft && (di == 2 || di == 6 || di == 7)) continue;
        V kn; if (di == 6) kn = {std::log(0.05), std::log(0.8), std::log(4.5)}; if (di == 7) kn = {0.7, 1.6};
        const int nk = static_cast<int>(kn.size()); const double* knp = nk ? kn.data() : nullptr;
        V th(64, 0.0), var(64 * 64), out(12), se(12);
        tally(sc_par_fit(n, p, di, aft, nk, knp, rep % 2 ? d.entry.data() : nullptr, d.time.data(), d.status.data(), X, nullptr, 100, 1e-8, th.data(), var.data(), &info));
        V tm = {0.5, 1.0, INFINITY};
        for (int type = 0; type < 5; ++type) {
          V tt = type == 3 ? V{0.25, 0.5, 0.75} : (type == 4 ? V{1.0, 2.0, 3.0} : tm);
          if (type == 1 || type == 2) tt[2] = 3.0;
          for (int avg = 0; avg < 2; ++avg)
            tally(sc_par_predict(di, aft, nk, knp, p, th.data(), var.data(), 4, X, 3, tt.data(), type, avg, out.data(), se.data()));
        }
        { V res(n); for (int rt = 0; rt < 3; ++rt)
            tally(sc_par_residuals(di, aft, nk, knp, p, th.data(), n, X, d.time.data(), rep % 2 ? d.entry.data() : nullptr, d.status.data(), rt, res.data())); }
        // interval-censored: types 0..3, start from the exact fit
        V t1(n), t2(n); I ty(n);
        for (int i = 0; i < n; ++i) { ty[i] = d.status[i] ? (i % 5 == 0 ? 3 : (i % 3 == 0 ? 1 : 2)) : 0; t1[i] = ty[i] == 2 ? 0.7 * d.time[i] : d.time[i]; t2[i] = d.time[i]; }
        V th2(th), var2(64 * 64);
        tally(sc_par_fit_ic(n, p, di, aft, nk, knp, t1.data(), t2.data(), ty.data(), X, th.data(), 100, 1e-8, th2.data(), var2.data(), &info));
        // without a start: the closed-form start on the midpoints
        tally(sc_par_fit_ic(n, p, di, aft, nk, knp, t1.data(), t2.data(), ty.data(), X, nullptr, 100, 1e-8, th2.data(), var2.data(), &info));
        // shared gamma frailty
        const int C = (n + 3) / 4; V th3(65, 0.0), var3(65 * 65), sh(C), rt(C);
        V ini(th.begin(), th.begin() + 64); ini.push_back(std::log(0.5));
        tally(sc_par_frailty(n, p, di, aft, nk, knp, d.time.data(), d.status.data(), X, d.cluster.data(), C, ini.data(), 100, 1e-8, th3.data(), var3.data(), &info, sh.data(), rt.data()));
        // cure
        for (int mix = 0; mix < 2; ++mix) {
          V tc(64, 0.0), vc(64 * 64);
          tally(sc_cure_fit(n, p, p, di, aft, nk, knp, mix, nullptr, d.time.data(), d.status.data(), X, X, nullptr, 100, 1e-8, tc.data(), vc.data(), &info));
          for (int type = 0; type < 2; ++type)
            tally(sc_cure_predict(di, aft, nk, knp, p, p, mix, tc.data(), vc.data(), 4, X, X, 3, tm.data(), type, 0, out.data(), se.data()));
          V res(n); for (int rt = 0; rt < 3; ++rt)
            tally(sc_cure_residuals(di, aft, nk, knp, p, p, mix, tc.data(), n, X, X, d.time.data(), nullptr, d.status.data(), rt, res.data()));
        }
      }
      if (K >= 2) {
        I ds = {1, 3, 4}; ds.resize(K, 1); I af(K, 0), nks(K, 0);
        V th(static_cast<size_t>(K) * (2 + p), 0.1), var(th.size() * th.size(), 0.0);
        for (size_t a = 0; a < th.size(); ++a) var[a * th.size() + a] = 0.01;
        V tm = {0.5, 1.0, 2.0}, out(12), se(12);
        for (int k = 0; k <= K; ++k) for (int avg = 0; avg < 2; ++avg)
          tally(sc_parcsc_predict(K, ds.data(), af.data(), nks.data(), nullptr, p, th.data(), var.data(), k, 4, X, 3, tm.data(), 4.0, avg, out.data(), se.data()));
      }
      // Cox with frailty
      if (p) {
        const int C = (n + 3) / 4; V est(p + 1), var(static_cast<size_t>(p + 1) * (p + 1)), tk(n), hk(n), sh(C), rt(C); int Kout = 0;
        V ini(p + 1, 0.0); ini[p] = std::log(0.5);
        tally(sc_cox_frailty(n, p, d.time.data(), d.status.data(), X, d.cluster.data(), C, ini.data(), 100, 1e-8, est.data(), var.data(), &info, &Kout, tk.data(), hk.data(), sh.data(), rt.data()));
      }
      // Turnbull
      { V L(n), R(n), lo(n), hi(n), mass(n), idd(2); I a(n), b(n), ii(2); int J = 0;
        for (int i = 0; i < n; ++i) { L[i] = d.status[i] ? 0.6 * d.time[i] : d.time[i]; R[i] = d.status[i] ? d.time[i] : INFINITY; if (i % 9 == 0) L[i] = 0; if (i % 11 == 0) L[i] = R[i] = d.time[i]; }
        tally(sc_turnbull(n, L.data(), R.data(), nullptr, 2000, 1e-8, &J, lo.data(), hi.data(), mass.data(), a.data(), b.data(), idd.data(), ii.data()));
        V w(n, 2.0);
        tally(sc_turnbull(n, L.data(), R.data(), w.data(), 2000, 1e-8, &J, lo.data(), hi.data(), mass.data(), a.data(), b.data(), idd.data(), ii.data()));
      }
      // kernel hazard
      { V ot(n), nr(n), ne(n), s(n), se(n); int m = 0;
        sc_kaplan_meier(n, d.time.data(), d.status.data(), nullptr, ot.data(), nr.data(), ne.data(), s.data(), se.data(), &m);
        V et, dd, Y; for (int j = 0; j < m; ++j) if (ne[j] > 0) { et.push_back(ot[j]); dd.push_back(ne[j]); Y.push_back(nr[j]); }
        V x = {0.1, 0.5, 1.0, 2.0, 3.0}, h(5), hs(5); double bw;
        for (int bnd = 0; bnd < 2; ++bnd) {
          tally(sc_kernel_hazard(static_cast<int>(et.size()), et.data(), dd.data(), Y.data(), -1.0, bnd, 5, x.data(), h.data(), hs.data(), &bw));
          tally(sc_kernel_hazard(static_cast<int>(et.size()), et.data(), dd.data(), Y.data(), 0.5, bnd, 5, x.data(), h.data(), hs.data(), &bw));
        }
      }
    }
  }
  // ---- invalid input: refused, never read or written out of bounds
  int refused = 0, tried = 0;
  { const int n = 10; V t(n, 1.0), o(n); I s(n, 1); int m;
    V tneg(t); tneg[3] = -1; V tnan(t); tnan[2] = NAN;
    ++tried; refused += sc_kaplan_meier(0, t.data(), s.data(), nullptr, o.data(), o.data(), o.data(), o.data(), o.data(), &m) != 0;
    ++tried; refused += sc_kaplan_meier(n, nullptr, s.data(), nullptr, o.data(), o.data(), o.data(), o.data(), o.data(), &m) != 0;
    ++tried; refused += sc_kaplan_meier(n, tneg.data(), s.data(), nullptr, o.data(), o.data(), o.data(), o.data(), o.data(), &m) != 0;
    ++tried; refused += sc_kaplan_meier(n, tnan.data(), s.data(), nullptr, o.data(), o.data(), o.data(), o.data(), o.data(), &m) != 0;
    V b(1), v(1); ++tried; refused += sc_coxph(n, 1, tnan.data(), s.data(), t.data(), nullptr, 0, 10, 1e-8, b.data(), v.data(), &info) != 0;
    I c(n, 0); ++tried; refused += sc_finegray(n, 1, t.data(), c.data(), t.data(), 1, 10, 1e-8, b.data(), v.data(), v.data(), &info) != 0;
    V L(n, 1.0), R(n, 0.5), lo(n), hi(n), ms(n), dd(2); I a(n), bb(n), ii(2); int J;
    ++tried; refused += sc_turnbull(n, L.data(), R.data(), nullptr, 100, 1e-8, &J, lo.data(), hi.data(), ms.data(), a.data(), bb.data(), dd.data(), ii.data()) != 0;
    I cl(n, 0); V e(2), vv(4), tk(n), hk(n), sh(1), rt(1); int Ko;
    ++tried; refused += sc_cox_frailty(n, 1, t.data(), s.data(), t.data(), cl.data(), 5, nullptr, 10, 1e-8, e.data(), vv.data(), &info, &Ko, tk.data(), hk.data(), sh.data(), rt.data()) != 0;
    I bad = {9}; I af = {0}; V av = {0}; I nk = {0}; V th(8), var(64), dg(5);
    ++tried; refused += sc_directcif_fit(n, 0, 1, t.data(), s.data(), nullptr, af.data(), av.data(), bad.data(), nk.data(), nullptr, 0, nullptr, 10, 1e-8, th.data(), var.data(), &info, dg.data()) != 0;
  }
  std::printf("calls %d (rc != 0: %d), invalid inputs refused %d of %d\n", calls, nonzero, refused, tried);
  return 0;
}
