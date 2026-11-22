// Memory and undefined-behaviour checks of the summary entry points of the core
// (src/capi_support.cpp): distributions, quantiles, intervals, frailty, adjusted risks, Turnbull
// summaries, prediction measures, the test of proportional hazards and time-varying effects, on
// random data with ties, empty groups and invalid input. Build as memcheck_capi.cpp (see README).
#include "survcore/survcore.h"
#include <cmath>
#include <cstdio>
#include <random>
#include <algorithm>
#include <vector>

using V = std::vector<double>;
using I = std::vector<int>;
static std::mt19937_64 rng(5);
static double U() { return std::uniform_real_distribution<double>(0, 1)(rng); }
static double N() { return std::normal_distribution<double>(0, 1)(rng); }
static int calls = 0, nonzero = 0;
static void tally(int rc) { ++calls; nonzero += rc != 0; }

#ifndef REPS
#define REPS 40
#endif

int main() {
  for (int rep = 0; rep < REPS; ++rep) {
    const int n = 5 + (rep * 29) % 220, p = 1 + rep % 3, K = 2 + rep % 2;
    V time(n), entry(n), X(static_cast<size_t>(n) * p), w(n), risk(static_cast<size_t>(n) * 4);
    I status(n), cause(n), g(n), strata(n), cl(n);
    for (int i = 0; i < n; ++i) {
      time[i] = std::round((-std::log(U()) * 3 + 0.05) * 10) / 10 + 0.1;
      entry[i] = U() < 0.3 ? time[i] * U() : 0.0;
      for (int j = 0; j < p; ++j) X[static_cast<size_t>(i) * p + j] = N();
      w[i] = 0.3 + 2 * U(); cause[i] = static_cast<int>(U() * (K + 1)); status[i] = cause[i] > 0;
      g[i] = i % 3; strata[i] = i % 2; cl[i] = i / 4;
      for (int j = 0; j < 4; ++j) risk[static_cast<size_t>(i) * 4 + j] = U();
    }
    status[0] = 1; cause[0] = 1;
    // distributions, quantiles, intervals, steps
    (void)sc_qnorm(U()); (void)sc_pnorm(N()); (void)sc_pchisq_upper(3 * U(), 1 + rep % 5); (void)sc_qgamma(U(), 0.01 + 5 * U(), 1 + U());
    { V pr = {0, 0.1, 0.5, 0.9, 1, 1.5}, out(6); tally(sc_quantile(n, time.data(), 6, pr.data(), 7, out.data())); tally(sc_quantile(n, time.data(), 6, pr.data(), 1, out.data())); tally(sc_quantile(0, nullptr, 6, pr.data(), 7, out.data())); }
    { V F(n), se(n, 0.05), lo(n), hi(n); for (int i = 0; i < n; ++i) F[i] = U(); F[0] = 0; F[n - 1] = 1; tally(sc_cif_interval(n, F.data(), se.data(), 1.96, lo.data(), hi.data())); }
    { V t = {0, 1, 1, 2}, v = {0.1, 0.2, 0.3, 0.4}, s = {-1, 0.5, 1, 3, NAN}, out(5); tally(sc_step_at(4, t.data(), v.data(), 5, s.data(), 0, out.data())); }
    { V et(n), d(n), Y(n); int m = 0; tally(sc_event_counts(n, time.data(), cause.data(), 1, entry.data(), et.data(), d.data(), Y.data(), &m));
      V H(m), se(m), lo(m), hi(m); if (m) tally(sc_na_table(m, d.data(), Y.data(), 1.96, H.data(), se.data(), lo.data(), hi.data()));
      V S(m, 0.5), kl(m), ku(m); if (m) tally(sc_km_limits(m, S.data(), se.data(), 1.96, kl.data(), ku.data()));
      double md; if (m) tally(sc_first_below(m, et.data(), S.data(), 0.5, &md)); }
    // frailty
    { const int C = (n + 3) / 4; V sh(C, 2.0), rt(C, 1.5), sz(C), ev(C), mn(C), lo(C), hi(C);
      tally(sc_frailty_clusters(n, status.data(), cl.data(), C, sh.data(), rt.data(), rep % 4 == 0, 0.95, sz.data(), ev.data(), mn.data(), lo.data(), hi.data()));
      double st, pv; tally(sc_frailty_lrt(-10, -11, 0, &st, &pv));
      V bt(n), bh(n, 0.01), beta(p, 0.1), tt = {0.5, 2, 9}, out(static_cast<size_t>(n) * 3), S(out.size());
      for (int i = 0; i < n; ++i) bt[i] = time[i];
      tally(sc_breslow_cumhaz(n, bt.data(), bh.data(), n, p, X.data(), beta.data(), 3, tt.data(), out.data()));
      V shr(n, 2.0), rtr(n, 1.0);
      for (int mode = 0; mode < 3; ++mode) tally(sc_frailty_survival(n, 3, out.data(), 0.5, mode, shr.data(), rtr.data(), S.data())); }
    // adjusted risks
    if (n >= 12) {
      const int q = p + 1; V Xi(static_cast<size_t>(n) * q);
      for (int i = 0; i < n; ++i) { Xi[static_cast<size_t>(i) * q] = 1; for (int j = 0; j < p; ++j) Xi[static_cast<size_t>(i) * q + j + 1] = X[static_cast<size_t>(i) * p + j]; }
      I a(n); for (int i = 0; i < n; ++i) a[i] = i % K;
      V P(static_cast<size_t>(n) * K), wt(n); int it;
      tally(sc_mlogit(n, q, Xi.data(), a.data(), K, P.data(), &it));
      tally(sc_ipw_weights(n, K, P.data(), a.data(), rep % 2 ? 0.05 : 0.0, wt.data()));
      V bef(static_cast<size_t>(K - 1) * p), aft(bef.size()), sz(K), ess(K), mw(K), pmn(K), pmx(K);
      tally(sc_balance(n, p, X.data(), a.data(), K, wt.data(), P.data(), bef.data(), aft.data(), sz.data(), ess.data(), mw.data(), pmn.data(), pmx.data()));
      V tt = {0.5, 1, 3, 20}, out(4);
      tally(sc_weighted_cif(n, time.data(), cause.data(), entry.data(), wt.data(), 1, 4, tt.data(), out.data()));
      tally(sc_weighted_cif(n, time.data(), cause.data(), nullptr, wt.data(), 2, 4, tt.data(), out.data()));
      V dr(static_cast<size_t>(30) * 4), sd(4), lo(4), hi(4); for (auto& v : dr) v = N(); dr[5] = NAN; dr[9] = INFINITY;
      tally(sc_boot_summary(30, 4, dr.data(), 0, 0.95, sd.data(), lo.data(), hi.data()));
      tally(sc_boot_summary(30, 4, dr.data(), 1, 0.95, sd.data(), lo.data(), hi.data()));
      I rs(static_cast<size_t>(7) * n), fo(n); tally(sc_resample(rep, n, 7, rs.data())); tally(sc_stratified_folds(n, status.data(), 5, rep, fo.data()));
    }
    // Turnbull
    {
      V L(n), R(n);
      for (int i = 0; i < n; ++i) { L[i] = status[i] ? 0.6 * time[i] : time[i]; R[i] = status[i] ? time[i] : INFINITY; if (i % 7 == 0) L[i] = R[i] = time[i]; if (i % 11 == 1) L[i] = 0; }
      V lo(n), hi(n), mass(n), dd(2); I a(n), b(n), ii(2); int J = 0;
      if (sc_turnbull(n, L.data(), R.data(), nullptr, 5000, 1e-8, &J, lo.data(), hi.data(), mass.data(), a.data(), b.data(), dd.data(), ii.data()) == 0) {
        V t = {0, 0.5, 1, 2, 5, 50}, e(6), l6(6), u6(6);
        tally(sc_tb_surv(J, lo.data(), hi.data(), mass.data(), 6, t.data(), rep % 2, e.data(), l6.data(), u6.data()));
        V st(2 * (2 * J + 1)), ss(st.size()); int m = 0;
        tally(sc_tb_steps(J, lo.data(), hi.data(), mass.data(), st.data(), ss.data(), &m));
        double qv; tally(sc_tb_quantile(m, st.data(), ss.data(), 0.5, &qv));
        if (n >= 6) { const int G = 2 + rep % 2; I gg(n); for (int i = 0; i < n; ++i) gg[i] = i % G; double s1, p1; int df; V sz(G), Uv(G), sdv(G);
          tally(sc_tb_test(n, L.data(), R.data(), gg.data(), G, 5000, 1e-8, &s1, &df, &p1, sz.data(), Uv.data(), sdv.data())); }
        V grid = {0, 0.5, 1, 2, 4, 8}, bo(static_cast<size_t>(12) * 6);
        tally(sc_tb_bootstrap(n, L.data(), R.data(), 6, grid.data(), 12, rep, 2000, 1e-8, bo.data()));
        V tq = {0.7, 1, 9, -1}, lq(4), uq(4); tally(sc_grid_boot_limits(6, grid.data(), 12, bo.data(), 4, tq.data(), 0.9, 10, lq.data(), uq.data()));
        V cl6(6), cu6(6); tally(sc_column_limits(12, 6, bo.data(), 0.9, 10, cl6.data(), cu6.data()));
        double ml, mu; tally(sc_tb_median_limits(6, grid.data(), 12, bo.data(), 0.9, &ml, &mu));
        V go = {0.2, 3, 10}, ao(static_cast<size_t>(12) * 3); tally(sc_approx_rows(12, 6, grid.data(), bo.data(), 3, go.data(), ao.data()));
      }
    }
    // prediction measures and ranking
    { V tt = {1, 2, 3, 4}, b(4), b0(4), ip(4), au(4), ca(4); tally(sc_measures(n, time.data(), cause.data(), 1, 4, tt.data(), risk.data(), 10, b.data(), b0.data(), ip.data(), au.data(), ca.data()));
      double ib; tally(sc_integrate_grid(4, tt.data(), b.data(), &ib));
      V np(4); tally(sc_np_curve(n, time.data(), cause.data(), 1, 4, tt.data(), np.data())); tally(sc_np_curve(n, time.data(), status.data(), 1, 4, tt.data(), np.data()));
      V loss = {0.2, NAN, 0.19, 0.21}, fl(16); for (auto& v : fl) v = 0.2 + 0.01 * N(); I rk(4), eq(4); V df(4), se(4);
      tally(sc_rank_one_se(4, 4, loss.data(), fl.data(), rk.data(), df.data(), se.data(), eq.data())); }
    // proportional hazards
    if (n >= 10) {
      V beta(p, 0.1), vc(static_cast<size_t>(p) * p, 0.0); for (int j = 0; j < p; ++j) vc[static_cast<size_t>(j) * p + j] = 0.01;
      I term(p); for (int j = 0; j < p; ++j) term[j] = j % 2;
      const int nterm = p > 1 ? 2 : 1;
      V cc(p), pc(p), ct(nterm), dt(nterm), pt(nterm), et(n), gg(n), sc(static_cast<size_t>(n) * p); double cg, pg; int nd; I rows(n);
      for (int tr = 0; tr < 5; ++tr)
        tally(sc_ph_test(n, p, time.data(), status.data(), X.data(), beta.data(), vc.data(), rep % 2 ? w.data() : nullptr, rep % 3 ? nullptr : entry.data(),
                         rep % 4 ? nullptr : strata.data(), tr, term.data(), nterm, cc.data(), pc.data(), ct.data(), dt.data(), pt.data(), &cg, &pg, &nd,
                         rows.data(), et.data(), gg.data(), sc.data()));
      V pv = {0.01, NAN, 0.2, 0.03}, po(4); tally(sc_p_adjust_holm(4, pv.data(), po.data()));
    }
    // time-varying effects
    {
      V et; for (int i = 0; i < n; ++i) if (status[i]) et.push_back(time[i]);
      std::sort(et.begin(), et.end()); et.erase(std::unique(et.begin(), et.end()), et.end());
      double sz; tally(sc_te_split_size(n, entry.data(), time.data(), static_cast<int>(et.size()), et.data(), &sz));
      const int cap = static_cast<int>(sz); I idx(cap), st(cap); V s0(cap), s1(cap); int nr;
      tally(sc_te_split(n, entry.data(), time.data(), status.data(), static_cast<int>(et.size()), et.data(), cap, idx.data(), s0.data(), s1.data(), st.data(), &nr));
      V tt = {0.01, 0.5, 1, 3, 30}, kn = {0.2, 0.9}, br = {1.0, 2.5}, B(5 * 3);
      for (int form = 1; form <= 4; ++form)
        tally(sc_te_basis(form, 5, tt.data(), 0.1, 1.0, form == 3 ? 2 : 0, br.data(), form == 4 ? 2 : 0, kn.data(), -2.0, 2.0, B.data()));
      V b = {0.3, -0.2, 0.1, 0.05}, Vc(16, 0.0); for (int j = 0; j < 4; ++j) Vc[j * 4 + j] = 0.01; I ix = {0, 2, 3};
      V e(5), s(5), h(5), l(5), u(5), pp(5); tally(sc_te_effect(5, 2, B.data(), 4, b.data(), Vc.data(), ix.data(), 1.96, e.data(), s.data(), h.data(), l.data(), u.data(), pp.data()));
      V tj = {0.5, 1, 2}, dH(static_cast<size_t>(3) * 4 * 3, 0.02), out(4 * 5);
      for (int mode = 0; mode < 3; ++mode) tally(sc_cum_predict(mode, 4, 3, tj.data(), mode == 2 ? 3 : 1, dH.data(), 1, rep % 2, 5, tt.data(), out.data()));
    }
  }
  // data expansions, timefix, sandwich
  for (int rep = 0; rep < REPS; ++rep) {
    const int n = 1 + (rep * 13) % 90;
    V time(n), entry(n), Sc(static_cast<size_t>(n) * 2), Vv = {0.2, 0.01, 0.01, 0.3}, out(4); I status(n), cl(n);
    for (int i = 0; i < n; ++i) { time[i] = std::round((-std::log(U()) * 3 + 0.1) * 10) / 10 + 0.1; entry[i] = U() < 0.4 ? std::round(time[i] * U() * 10) / 10 : 0.0;
      if (entry[i] >= time[i]) entry[i] = 0; status[i] = static_cast<int>(U() * 4); cl[i] = i % 5; Sc[2 * i] = N(); Sc[2 * i + 1] = N(); }
    for (int cause = 1; cause <= 2; ++cause) {
      int nr = 0; tally(sc_fg_expand(n, time.data(), status.data(), rep % 2 ? entry.data() : nullptr, cause, 0, nullptr, nullptr, nullptr, nullptr, nullptr, &nr));
      if (nr) { I id(nr), st(nr); V a(nr), b(nr), w(nr); tally(sc_fg_expand(n, time.data(), status.data(), rep % 2 ? entry.data() : nullptr, cause, nr, id.data(), a.data(), b.data(), st.data(), w.data(), &nr)); }
    }
    V tv = time; tv.insert(tv.end(), entry.begin(), entry.end()); tv.push_back(INFINITY); tv.push_back(NAN); int ch;
    for (size_t i = 0; i + 1 < tv.size(); i += 3) tv[i + 1] = tv[i] + 1e-12;
    tally(sc_timefix(static_cast<int>(tv.size()), tv.data(), 1.5e-8, &ch));
    tally(sc_sandwich(n, 2, Sc.data(), cl.data(), 5, Vv.data(), out.data())); tally(sc_sandwich(n, 2, Sc.data(), nullptr, 0, Vv.data(), out.data()));
  }
  // smoother, covariance, elementwise formulas, and draws from an installed index source
  for (int rep = 0; rep < REPS; ++rep) {
    const int n = 3 + (rep * 17) % 120, df = 1 + rep % 5;
    V x(n), y(n), g = {0.0, 0.5, 1.0, 2.0, 2.5}, f(5), lo(5), hi(5);
    for (int i = 0; i < n; ++i) { x[i] = rep % 7 == 0 ? std::floor(3 * U()) : 3 * U(); y[i] = N(); }
    tally(sc_ns_smooth(n, x.data(), y.data(), df, 5, g.data(), 1.96, f.data(), lo.data(), hi.data()) == -4 ? 0 : sc_ns_smooth(n, x.data(), y.data(), df, 5, g.data(), 1.96, f.data(), lo.data(), hi.data()));
    V draws(static_cast<size_t>(n) * 3), cv(9); int nok;
    for (auto& v : draws) v = U() < 0.1 ? NAN : N();
    tally(sc_cov_complete(n, 3, draws.data(), cv.data(), &nok));
    V Xs(static_cast<size_t>(n) * 2), b = {0.3, -0.2}, H = {0, 0.1, 1, 50}, out(static_cast<size_t>(n) * 4);
    for (auto& v : Xs) v = N();
    tally(sc_subdist_cif(n, 2, Xs.data(), b.data(), 4, H.data(), out.data()));
    V e(n), se(n), l(n), u(n); for (int i = 0; i < n; ++i) { e[i] = N(); se[i] = i % 9 == 0 ? NAN : U(); }
    tally(sc_log_limits(n, e.data(), se.data(), 1.64, l.data(), u.data()));
    sc_set_index_source([](double m) { return std::floor(U() * m); });
    I rs(static_cast<size_t>(n) * 3), fo(n), st(n); for (int i = 0; i < n; ++i) st[i] = i % 3;
    tally(sc_resample(0, n, 3, rs.data())); tally(sc_stratified_folds(n, st.data(), 4, 0, fo.data()));
    int bad = 0; for (int v : rs) bad += v < 0 || v >= n; for (int v : fo) bad += v < 1 || v > 4;
    sc_set_index_source([](double m) { return m + 5.0; });             // out of range: clamped to n - 1
    tally(sc_resample(0, n, 2, rs.data())); for (int k = 0; k < 2 * n; ++k) bad += rs[k] != n - 1;
    sc_set_index_source(nullptr);
    if (bad) { std::printf("index source: %d indices out of range\n", bad); return 1; }
  }
  // toolkit: elementwise, reductions, means, running sums, predictions, tables, grids, increments
  for (int rep = 0; rep < REPS; ++rep) {
    const int n = rep % 9, nc = 1 + rep % 4;
    V x(n + 1), y(n + 1), o(n + 1), M(static_cast<size_t>(n) * nc + 1), oc(static_cast<size_t>(n) * nc + nc + n + 1);
    for (int i = 0; i <= n; ++i) { x[i] = i % 5 == 0 ? NAN : (rep % 2 ? N() : U()); y[i] = i % 7 == 3 ? 0.0 : U() + 0.1; }
    for (auto& v : M) v = U() < 0.1 ? NAN : U();
    for (int code = 1; code <= 23; ++code) tally(sc_elementwise(code, n, x.data(), y.data(), o.data()));
    double r; for (int code = 0; code <= 4; ++code) { tally(sc_reduce(code, n, x.data(), 1, &r)); tally(sc_reduce(code, n, x.data(), 0, &r)); }
    tally(sc_means(n, nc, M.data(), 0, 1, oc.data())); tally(sc_means(n, nc, M.data(), 1, 0, oc.data()));
    I g(n + 1); for (int i = 0; i <= n; ++i) g[i] = i % 3; V gm(3); if (n) tally(sc_group_means(n, g.data(), x.data(), 3, gm.data()));
    tally(sc_cumsum_rows(n, nc, M.data(), rep % 2 ? 1.0 : -1.0, oc.data()));
    V X(static_cast<size_t>(n) * 2 + 1), b = {0.2, -0.4}, H = {0, 0.3, 2}, Hr(static_cast<size_t>(n) * 3 + 1), pr(static_cast<size_t>(n) * 3 + 1);
    for (auto& v : X) v = N(); for (auto& v : Hr) v = U();
    tally(sc_linpred(n, 2, X.data(), b.data(), o.data()));
    for (int mode = 0; mode < 2; ++mode) { tally(sc_ph_predict(mode, n, 2, X.data(), b.data(), 3, H.data(), 0, pr.data())); tally(sc_ph_predict(mode, n, 2, X.data(), b.data(), 3, Hr.data(), 1, pr.data())); }
    V d(n + 1, 1.0), Yr(n + 1); for (int i = 0; i <= n; ++i) Yr[i] = n + 1 - i;
    for (int mode = 0; mode < 3; ++mode) tally(sc_np_curve_inc(n, d.data(), Yr.data(), mode, o.data()));
    V ll(n + 1), k(n + 1, 3.0), nob(n + 1, 50.0), a1(n + 1), a2(n + 1); for (auto& v : ll) v = -100 * U();
    tally(sc_info_criteria(n, ll.data(), k.data(), nob.data(), a1.data(), a2.data())); tally(sc_info_criteria(n, ll.data(), k.data(), nullptr, a1.data(), a2.data()));
    tally(sc_lr_test(n, ll.data(), a1.data(), k.data(), a1.data(), a2.data()));
    I ex(n + 1); for (int i = 0; i <= n; ++i) ex[i] = i % 2;
    V w1(n + 1), w2(n + 1), w3(n + 1), w4(n + 1), w5(n + 1), w6(n + 1), w7(n + 1);
    tally(sc_wald_table(n, x.data(), y.data(), 1.96, ex.data(), w1.data(), w2.data(), w3.data(), w4.data(), w5.data(), w6.data(), w7.data()));
    V Vm = {0.2, 0.01, 0.01, -0.1}, se(2); tally(sc_se_diag(2, Vm.data(), se.data()));
    V S(static_cast<size_t>(n) * 2 + 1), so(static_cast<size_t>(n) * 2 + 1); for (auto& v : S) v = N();
    tally(sc_scaled_schoenfeld(n, 2, S.data(), Vm.data(), b.data(), so.data()));
    I st(n + 1); for (int i = 0; i <= n; ++i) st[i] = i % 2; tally(sc_hazard_residuals(n, y.data(), st.data(), rep % 5, o.data()));
    int m; double from = rep % 3 ? 0.0 : -1.5, to = 0.5 + 40 * U();
    for (int kind = 0; kind < 3; ++kind) {
      if (sc_grid(kind, from, to, 0.3, 1 + rep % 7, 0, nullptr, &m) != 0) continue;
      V gr(m + 1); tally(sc_grid(kind, from, to, 0.3, 1 + rep % 7, m, gr.data(), &m));
    }
    V bt = {0.5, 1.0, 2.0, 0.7, 1.5}, bc = {0.1, 0.3, 0.4, 0.2, 0.6}, tj = {0.5, 0.7, 1.0, 1.5, 2.0}, LP(static_cast<size_t>(n) * 5 + 1), dH(static_cast<size_t>(n) * 5 + 1);
    I bs = {0, 0, 0, 1, 1}; for (auto& v : LP) v = N();
    tally(sc_hazard_increments(5, bt.data(), bc.data(), bs.data(), n, st.data(), 5, tj.data(), LP.data(), rep % 2, dH.data()));
  }
  // collinearity check, on designs with constant, duplicated and badly scaled columns
  for (int rep = 0; rep < REPS; ++rep) {
    const int n = rep % 30, p = 1 + rep % 5;
    V X(static_cast<size_t>(n) * p + 1); I dep(p);
    for (int i = 0; i < n; ++i) for (int j = 0; j < p; ++j) X[static_cast<size_t>(i) * p + j] = (j == 2 ? 2 * X[static_cast<size_t>(i) * p] : (j == 3 ? 1.0 : N() * (rep % 3 ? 1e6 : 1e-6)));
    int nd; tally(sc_collinear(n, p, X.data(), 1e-7, dep.data(), &nd));
  }
  // invalid input
  int refused = 0, tried = 0;
  { V x(3, 1.0), o(3); I s(3, 1);
    ++tried; refused += sc_quantile(3, x.data(), 3, x.data(), 5, o.data()) != 0;
    ++tried; refused += sc_step_at(3, nullptr, x.data(), 1, x.data(), 0, o.data()) != 0;
    V t = {2, 1}; ++tried; refused += sc_step_at(2, t.data(), x.data(), 1, x.data(), 0, o.data()) != 0;
    I bad = {0, 5, 1}; V P(9, 1.0 / 3); ++tried; refused += sc_ipw_weights(3, 3, P.data(), bad.data(), 0, o.data()) != 0;
    ++tried; refused += sc_te_basis(9, 1, x.data(), 1, 1, 0, nullptr, 0, nullptr, 0, 1, o.data()) != 0;
    I badcl = {0, 3, 1}; ++tried; refused += sc_sandwich(3, 1, x.data(), badcl.data(), 2, x.data(), o.data()) != 0;
    V neg = {1, 2, 3}, en = {0, 2.5, 0}; int nfe; ++tried; refused += sc_fg_expand(3, neg.data(), s.data(), en.data(), 1, 0, nullptr, nullptr, nullptr, nullptr, nullptr, &nfe) != 0;
    ++tried; refused += sc_cum_predict(2, 1, 1, x.data(), 1, x.data(), 3, 0, 1, x.data(), o.data()) != 0;
    I term = {7}; double d; int nd; I rows(3);
    ++tried; refused += sc_ph_test(3, 1, x.data(), s.data(), x.data(), x.data(), x.data(), nullptr, nullptr, nullptr, 0, term.data(), 1, o.data(), o.data(), o.data(), o.data(), o.data(), &d, &d, &nd, rows.data(), o.data(), o.data(), o.data()) != 0;
  }
  std::printf("calls %d (rc != 0: %d), invalid inputs refused %d of %d\n", calls, nonzero, refused, tried);
  return 0;
}
