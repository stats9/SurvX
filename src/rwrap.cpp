// SurvX: thin wrappers that adapt the C API of the core (survcore/survcore.h) to R's .C interface.
#include "rwrap.h"
#include "survcore/survcore.h"

#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <R.h>
#include <R_ext/Random.h>

namespace {
// every random draw of the core comes from R's generator while an RStream exists: the state is read
// from .Random.seed, R_unif_index gives the indices (as sample.int does), and the state is written back
double r_index(double n) { return R_unif_index(n); }
struct RStream {
  RStream() { GetRNGstate(); sc_set_index_source(r_index); }
  ~RStream() { sc_set_index_source(nullptr); PutRNGstate(); }
  RStream(const RStream&) = delete;
  RStream& operator=(const RStream&) = delete;
};
inline void put(const sc_fit_info& f, double* d, int* i) {
  d[0] = f.loglik_init; d[1] = f.loglik; d[2] = f.grad_maxabs;
  i[0] = f.iterations; i[1] = f.converged; i[2] = f.flags; i[3] = f.n_events;
}
}  // namespace

extern "C" {
void SX_km(int* n, double* time, int* status, double* out_time, double* out_nrisk, double* out_nevent, double* out_surv, double* out_se, int* m, int* rc) {
  *rc = sc_kaplan_meier(*n, time, status, nullptr, out_time, out_nrisk, out_nevent, out_surv, out_se, m);
}
void SX_rmst(int* n, double* time, int* status, double* tau, double* mean, double* se, int* rc) {
  *rc = sc_rmst(*n, time, status, *tau, mean, se);
}
void SX_aj(int* n, double* time, int* cause, int* K, double* out_time, double* out_surv, double* out_cif, double* out_se, int* m, int* rc) {
  *rc = sc_aalen_johansen(*n, time, cause, *K, out_time, out_surv, out_cif, out_se, m);
}
void SX_cox2(int* n, int* p, int* has_entry, double* entry, double* time, int* status, int* has_strata, int* strata, double* X, int* has_w, double* weights, int* ties, int* max_iter, double* tol, double* beta, double* var, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_coxph_ext(*n, *p, *has_entry ? entry : nullptr, time, status, *has_strata ? strata : nullptr, X, *has_w ? weights : nullptr,
                     *ties, *max_iter, *tol, beta, var, &f);
  put(f, info_d, info_i);
}
void SX_cox2_baseline(int* n, int* p, int* has_entry, double* entry, double* time, int* status, int* has_strata, int* strata, double* X, int* has_w, double* weights, int* ties, double* beta, double* out_time, double* out_cumhaz, int* out_strata, int* m, int* rc) {
  *rc = sc_coxph_ext_baseline(*n, *p, *has_entry ? entry : nullptr, time, status, *has_strata ? strata : nullptr, X, *has_w ? weights : nullptr,
                              *ties, beta, nullptr, out_time, out_cumhaz, out_strata, m);
}
void SX_cox_resid(int* n, int* p, int* has_entry, double* entry, double* time, int* status, int* has_strata, int* strata, double* X, int* has_w, double* weights, int* ties, double* beta, double* mart, double* xbar, int* rc) {
  *rc = sc_coxph_resid(*n, *p, *has_entry ? entry : nullptr, time, status, *has_strata ? strata : nullptr, X, *has_w ? weights : nullptr,
                       *ties, beta, mart, xbar);
}
void SX_cox_score(int* n, int* p, int* has_entry, double* entry, double* time, int* status, int* has_strata, int* strata, double* X, int* has_w, double* weights, int* ties, double* beta, double* U, int* rc) {
  *rc = sc_coxph_score(*n, *p, *has_entry ? entry : nullptr, time, status, *has_strata ? strata : nullptr, X, *has_w ? weights : nullptr,
                       *ties, beta, U);
}
void SX_fg(int* n, int* p, double* time, int* cause, double* X, int* interest, int* max_iter, double* tol, double* beta, double* var_model, double* var_robust, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_finegray(*n, *p, time, cause, X, *interest, *max_iter, *tol, beta, var_model, var_robust, &f);
  put(f, info_d, info_i);
}
void SX_fg_baseline(int* n, int* p, double* time, int* cause, double* X, int* interest, double* beta, double* out_time, double* out_cumhaz, int* m, int* rc) {
  *rc = sc_finegray_baseline(*n, *p, time, cause, X, *interest, beta, nullptr, out_time, out_cumhaz, m);
}
void SX_csc_cif(int* n, int* p, double* time, int* cause, double* X, int* K, double* betas, int* ties, int* n_new, double* Xnew, int* product_limit, double* out_time, double* out_cif, int* m, int* rc) {
  *rc = sc_csc_cif(*n, *p, time, cause, X, *K, betas, *ties, *n_new, Xnew, *product_limit, out_time, out_cif, m);
}
void SX_logrank(int* n, double* time, int* status, int* group, int* K, int* has_strata, int* strata, double* rho, double* stat, int* df, double* p_value, double* out_obs, double* out_exp, double* out_var, int* rc) {
  *rc = sc_logrank(*n, time, status, group, *K, *has_strata ? strata : nullptr, *rho, stat, df, p_value, out_obs, out_exp, out_var);
}
void SX_gray(int* n, double* time, int* cause, int* group, int* K, int* has_strata, int* strata, int* interest, double* rho, double* stat, int* df, double* p_value, double* out_score, double* out_var, int* rc) {
  *rc = sc_gray_test(*n, time, cause, group, *K, *has_strata ? strata : nullptr, *interest, *rho, stat, df, p_value, out_score, out_var);
}
void SX_direct_knots(int* n, double* time, int* cause, int* k, int* df, double* out_knots, int* n_knots, int* rc) {
  *rc = sc_directcif_default_knots(*n, time, cause, *k, *df, out_knots, n_knots);
}
void SX_direct_fit(int* n, int* p, int* K, double* time, int* cause, double* X, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, int* constraint, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, double* diag, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_directcif_fit(*n, *p, *K, time, cause, X, alpha_free, alpha_value, baseline, n_knots, knots, *constraint,
                         *has_init ? init : nullptr, *max_iter, *tol, theta, var, &f, diag);
  put(f, info_d, info_i);
}
void SX_direct_predict(int* p, int* K, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc) {
  *rc = sc_directcif_predict(*p, *K, alpha_free, alpha_value, baseline, n_knots, knots, theta, *has_var ? var : nullptr, *k, *n_new, Xnew, *nt, times,
                             out_cif, *has_var ? out_se : nullptr);
}
void SX_brier(int* n, double* time, int* cause, int* interest, double* pred, double* t, double* brier, double* brier_null, double* ipa, int* rc) {
  *rc = sc_brier(*n, time, cause, *interest, pred, *t, brier, brier_null, ipa);
}
void SX_auc(int* n, double* time, int* cause, int* interest, double* pred, double* t, double* auc, int* rc) {
  *rc = sc_auc(*n, time, cause, *interest, pred, *t, auc);
}
void SX_calibration(int* n, double* time, int* cause, int* interest, double* pred, double* t, int* groups, double* out_mean_pred, double* out_observed, double* out_size, double* error, int* m, int* rc) {
  *rc = sc_calibration(*n, time, cause, *interest, pred, *t, *groups, out_mean_pred, out_observed, out_size, error, m);
}
void SX_compare(int* n, int* p, int* K, double* time, int* cause, double* X, int* interest, int* nt, double* times, int* folds, int* n_direct, int* baseline, int* link, int* spline_df, int* estimation, int* cond, double* out_metrics, double* out_fit, int* rc) {
  RStream rs;
  *rc = sc_compare_models2(*n, *p, *K, time, cause, X, *interest, *nt, times, *folds, 0,
                           *n_direct, baseline, link, spline_df, estimation, cond, out_metrics, out_fit);
}
void SX_subdist_fit(int* n, int* p, double* time, int* cause, double* X, int* interest, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, int* grid, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var_robust, double* var_model, double* info_d, int* info_i, int* n_support, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_subdist_fit(*n, *p, time, cause, X, *interest, *alpha_free, *alpha_value, *baseline, *n_knots, knots, *grid,
                       *has_init ? init : nullptr, *max_iter, *tol, theta, var_robust, var_model, &f, n_support);
  put(f, info_d, info_i);
}
void SX_subdist_var_ipcw(int* n, int* p, double* time, int* cause, double* X, int* interest, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, double* theta, double* var_model, double* var, int* rc) {
  *rc = sc_subdist_var_ipcw(*n, *p, time, cause, X, *interest, *alpha_free, *alpha_value, *baseline, *n_knots, knots, theta, var_model, var);
}
void SX_constrained_fit(int* n, int* p, double* time, int* cause, double* X, int* alpha_free, double* alpha_value, int* baseline, int* family, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_constrained_fit(*n, *p, time, cause, X, *alpha_free, *alpha_value, *baseline, *family, *has_init ? init : nullptr, *max_iter, *tol, theta, var, &f);
  put(f, info_d, info_i);
}
void SX_constrained_predict(int* p, int* alpha_free, double* alpha_value, int* baseline, int* family, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc) {
  *rc = sc_constrained_predict(*p, *alpha_free, *alpha_value, *baseline, *family, theta, *has_var ? var : nullptr, *k, *n_new, Xnew, *nt, times,
                               out_cif, *has_var ? out_se : nullptr);
}
void SX_constrainedK_fit(int* n, int* p, double* time, int* cause, double* X, int* K, int* alpha_free, double* alpha_value, int* baseline, int* family, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_constrainedK_fit(*n, *p, time, cause, X, *K, *alpha_free, *alpha_value, *baseline, family, *has_init ? init : nullptr, *max_iter, *tol, theta, var, &f);
  put(f, info_d, info_i);
}
void SX_constrainedK_predict(int* p, int* K, int* alpha_free, double* alpha_value, int* baseline, int* family, double* theta, int* has_var, double* var, int* k, int* average, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc) {
  *rc = sc_constrainedK_predict(*p, *K, *alpha_free, *alpha_value, *baseline, family, theta, *has_var ? var : nullptr, *k, *average, *n_new, Xnew, *nt, times,
                                out_cif, *has_var ? out_se : nullptr);
}
// constrained model with covariates of the G_k: W = V = X, the flags (K - 1) x p (row-major) select
// the columns of X for the time ratio (wuse) and the shape (vuse) of each competing cause
void SX_constrainedK_fit_cov(int* n, int* p, double* time, int* cause, double* X, int* K, int* alpha_free, double* alpha_value, int* baseline, int* family, int* wuse, int* vuse, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_constrainedK_fit_cov(*n, *p, time, cause, X, *K, *alpha_free, *alpha_value, *baseline, family, *p, X, wuse, *p, X, vuse,
                                *has_init ? init : nullptr, *max_iter, *tol, theta, var, &f);
  put(f, info_d, info_i);
}
void SX_constrainedK_predict_cov(int* p, int* K, int* alpha_free, double* alpha_value, int* baseline, int* family, int* wuse, int* vuse, double* theta, int* has_var, double* var, int* k, int* average, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc) {
  *rc = sc_constrainedK_predict_cov(*p, *K, *alpha_free, *alpha_value, *baseline, family, *p, wuse, *p, vuse, theta, *has_var ? var : nullptr,
                                    *k, *average, *n_new, Xnew, Xnew, Xnew, *nt, times, out_cif, *has_var ? out_se : nullptr);
}
void SX_direct_predict_mean(int* p, int* K, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc) {
  *rc = sc_directcif_predict_mean(*p, *K, alpha_free, alpha_value, baseline, n_knots, knots, theta, *has_var ? var : nullptr, *k, *n_new, Xnew, *nt, times,
                                  out_cif, *has_var ? out_se : nullptr);
}
void SX_constrained_predict_mean(int* p, int* alpha_free, double* alpha_value, int* baseline, int* family, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc) {
  *rc = sc_constrained_predict_mean(*p, *alpha_free, *alpha_value, *baseline, *family, theta, *has_var ? var : nullptr, *k, *n_new, Xnew, *nt, times,
                                    out_cif, *has_var ? out_se : nullptr);
}
void SX_cif_bootstrap(int* model, int* n, int* p, double* time, int* cause, double* X, int* K, int* interest, int* ties, int* product_limit, int* adjust, int* n_new, double* Xnew, int* average, int* nt, double* times, int* B, double* level, double* out_est, double* out_se, double* out_lower, double* out_upper, int* n_ok, int* rc) {
  RStream rs;
  *rc = sc_cif_bootstrap(*model, *n, *p, time, cause, X, *K, *interest, *ties, *product_limit, *adjust, *n_new, Xnew, *average, *nt, times, *B,
                         0, *level, out_est, out_se, out_lower, out_upper, n_ok);
}
void SX_par_fit(int* n, int* p, int* dist, int* aft, int* nk, double* knots, int* has_entry, double* entry, double* time, int* status, double* X, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_par_fit(*n, *p, *dist, *aft, *nk, knots, *has_entry ? entry : nullptr, time, status, X, *has_init ? init : nullptr, *max_iter, *tol, theta, var, &f);
  put(f, info_d, info_i);
}
void SX_par_predict(int* dist, int* aft, int* nk, double* knots, int* p, double* theta, int* has_var, double* var, int* n_new, double* Xnew, int* nt, double* times, int* type, int* average, double* out, double* out_se, int* rc) {
  *rc = sc_par_predict(*dist, *aft, *nk, knots, *p, theta, *has_var ? var : nullptr, *n_new, Xnew, *nt, times, *type, *average, out, *has_var ? out_se : nullptr);
}
void SX_parcsc_predict(int* K, int* dist, int* aft, int* nk, double* knots, int* p, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* tmax, int* average, double* out, double* out_se, int* rc) {
  *rc = sc_parcsc_predict(*K, dist, aft, nk, knots, *p, theta, *has_var ? var : nullptr, *k, *n_new, Xnew, *nt, times, *tmax, *average, out, *has_var ? out_se : nullptr);
}
void SX_cure_fit(int* n, int* pl, int* pc, int* dist, int* aft, int* nk, double* knots, int* mixture, int* has_entry, double* entry, double* time, int* status, double* Xl, double* Xc, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_cure_fit(*n, *pl, *pc, *dist, *aft, *nk, knots, *mixture, *has_entry ? entry : nullptr, time, status, Xl, Xc,
                    *has_init ? init : nullptr, *max_iter, *tol, theta, var, &f);
  put(f, info_d, info_i);
}
void SX_cure_predict(int* dist, int* aft, int* nk, double* knots, int* pl, int* pc, int* mixture, double* theta, int* has_var, double* var, int* n_new, double* Xl, double* Xc, int* nt, double* times, int* type, int* average, double* out, double* out_se, int* rc) {
  *rc = sc_cure_predict(*dist, *aft, *nk, knots, *pl, *pc, *mixture, theta, *has_var ? var : nullptr, *n_new, Xl, Xc, *nt, times,
                        *type, *average, out, *has_var ? out_se : nullptr);
}
void SX_cv_folds(int* n, int* folds, int* out, int* rc) {
  RStream rs;
  *rc = sc_cv_folds(*n, *folds, 0, out);
}
void SX_kernel_hazard(int* m, double* et, double* d, double* Y, double* bw, int* boundary, int* nx, double* x, double* out_h, double* out_se, double* bw_used, int* rc) {
  *rc = sc_kernel_hazard(*m, et, d, Y, *bw, *boundary, *nx, x, out_h, out_se, bw_used);
}
void SX_par_fit_ic(int* n, int* p, int* dist, int* aft, int* nk, double* knots, double* time, double* time2, int* type, double* X, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_par_fit_ic(*n, *p, *dist, *aft, *nk, knots, time, time2, type, X, *has_init ? init : nullptr, *max_iter, *tol, theta, var, &f);
  put(f, info_d, info_i);
}
void SX_direct_fit_ic(int* n, int* p, int* K, double* time, double* time2, int* cause, double* X, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, int* constraint, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, double* diag, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_directcif_fit_ic(*n, *p, *K, time, time2, cause, X, alpha_free, alpha_value, baseline, n_knots, knots, *constraint,
                            *has_init ? init : nullptr, *max_iter, *tol, theta, var, &f, diag);
  put(f, info_d, info_i);
}
void SX_turnbull(int* n, double* L, double* R, int* has_w, double* w, int* max_iter, double* tol, int* n_int, double* lo, double* hi, double* mass, int* a, int* b, double* info_d, int* info_i, int* rc) {
  *rc = sc_turnbull(*n, L, R, *has_w ? w : nullptr, *max_iter, *tol, n_int, lo, hi, mass, a, b, info_d, info_i);
}
void SX_par_frailty(int* n, int* p, int* dist, int* aft, int* nk, double* knots, double* time, int* status, double* X, int* cluster, int* C, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, double* shape, double* rate, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_par_frailty(*n, *p, *dist, *aft, *nk, knots, time, status, X, cluster, *C, init, *max_iter, *tol, theta, var, &f, shape, rate);
  put(f, info_d, info_i);
}
void SX_cox_frailty(int* n, int* p, double* time, int* status, double* X, int* cluster, int* C, double* init, int* max_iter, double* tol, double* est, double* var, double* info_d, int* info_i, int* K, double* tk, double* hk, double* shape, double* rate, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_cox_frailty(*n, *p, time, status, X, cluster, *C, init, *max_iter, *tol, est, var, &f, K, tk, hk, shape, rate);
  put(f, info_d, info_i);
}

// ---------------- numerical summaries (capi_support.cpp); vector arguments of equal length ----------------
void SX_qnorm(int* n, double* p, double* out) { for (int i = 0; i < *n; ++i) out[i] = sc_qnorm(p[i]); }
void SX_pnorm(int* n, double* x, double* out) { for (int i = 0; i < *n; ++i) out[i] = sc_pnorm(x[i]); }
void SX_pchisq_upper(int* n, double* x, double* df, double* out) { for (int i = 0; i < *n; ++i) out[i] = sc_pchisq_upper(x[i], df[i]); }
void SX_qgamma(int* n, double* p, double* shape, double* rate, double* out) {
  for (int i = 0; i < *n; ++i) out[i] = sc_qgamma(p[i], shape[i], rate[i]);
}
void SX_quantile(int* n, double* x, int* np, double* p, int* type, double* out, int* rc) {
  *rc = sc_quantile(*n, x, *np, p, *type, out);
}
void SX_cif_interval(int* n, double* F, double* se, double* zq, double* lower, double* upper, int* rc) {
  *rc = sc_cif_interval(*n, F, se, *zq, lower, upper);
}
void SX_step_at(int* m, double* t, double* v, int* ns, double* s, double* before, double* out, int* rc) {
  *rc = sc_step_at(*m, t, v, *ns, s, *before, out);
}
void SX_frailty_lrt(double* ll, double* ll0, int* boundary, double* stat, double* p, int* rc) {
  *rc = sc_frailty_lrt(*ll, *ll0, *boundary, stat, p);
}
void SX_frailty_clusters(int* n, int* status, int* cluster, int* C, double* shape, double* rate, int* boundary, double* level,
                         double* size, double* events, double* mean, double* lower, double* upper, int* rc) {
  *rc = sc_frailty_clusters(*n, status, cluster, *C, shape, rate, *boundary, *level, size, events, mean, lower, upper);
}
void SX_breslow_cumhaz(int* nb, double* bh_time, double* bh_hazard, int* n, int* p, double* X, double* beta, int* nt, double* times,
                       double* out, int* rc) {
  *rc = sc_breslow_cumhaz(*nb, bh_time, bh_hazard, *n, *p, X, beta, *nt, times, out);
}
void SX_frailty_survival(int* n, int* nt, double* H, double* theta, int* mode, double* shape, double* rate, double* out, int* rc) {
  *rc = sc_frailty_survival(*n, *nt, H, *theta, *mode, shape, rate, out);
}
void SX_na_table(int* m, double* d, double* Y, double* zq, double* H, double* se, double* lower, double* upper, int* rc) {
  *rc = sc_na_table(*m, d, Y, *zq, H, se, lower, upper);
}

void SX_km_limits(int* m, double* S, double* se, double* zq, double* lower, double* upper, int* rc) {
  *rc = sc_km_limits(*m, S, se, *zq, lower, upper);
}
void SX_first_below(int* m, double* t, double* v, double* p, double* out, int* rc) { *rc = sc_first_below(*m, t, v, *p, out); }
void SX_event_counts(int* n, double* time, int* status, int* k, int* has_entry, double* entry, double* et, double* d, double* Y,
                     int* m, int* rc) {
  *rc = sc_event_counts(*n, time, status, *k, *has_entry ? entry : nullptr, et, d, Y, m);
}

void SX_mlogit(int* n, int* p, double* X, int* a, int* K, double* P, int* iterations, int* rc) {
  *rc = sc_mlogit(*n, *p, X, a, *K, P, iterations);
}
void SX_ipw_weights(int* n, int* K, double* P, int* a, double* trim, double* w, int* rc) { *rc = sc_ipw_weights(*n, *K, P, a, *trim, w); }
void SX_balance(int* n, int* q, double* Xc, int* a, int* K, double* w, double* P, double* before, double* after, double* size,
                double* ess, double* max_weight, double* ps_min, double* ps_max, int* rc) {
  *rc = sc_balance(*n, *q, Xc, a, *K, w, P, before, after, size, ess, max_weight, ps_min, ps_max);
}
void SX_weighted_cif(int* n, double* time, int* status, int* has_entry, double* entry, double* w, int* k, int* nt, double* tt,
                     double* out, int* rc) {
  *rc = sc_weighted_cif(*n, time, status, *has_entry ? entry : nullptr, w, *k, *nt, tt, out);
}
void SX_boot_summary(int* nb, int* m, double* draws, int* finite_only, double* level, double* sd, double* lower, double* upper,
                     int* rc) {
  *rc = sc_boot_summary(*nb, *m, draws, *finite_only, *level, sd, lower, upper);
}
void SX_resample(int* n, int* B, int* out, int* rc) {
  RStream rs;
  *rc = sc_resample(0, *n, *B, out);
}
void SX_stratified_folds(int* n, int* strata, int* folds, int* out, int* rc) {
  RStream rs;
  *rc = sc_stratified_folds(*n, strata, *folds, 0, out);
}

void SX_tb_surv(int* J, double* lo, double* hi, double* p, int* nt, double* t, int* left, double* est, double* lower,
                double* upper, int* rc) {
  *rc = sc_tb_surv(*J, lo, hi, p, *nt, t, *left, est, lower, upper);
}
void SX_tb_steps(int* J, double* lo, double* hi, double* p, double* time, double* surv, int* m, int* rc) {
  *rc = sc_tb_steps(*J, lo, hi, p, time, surv, m);
}
void SX_tb_quantile(int* m, double* time, double* surv, double* p, double* out, int* rc) { *rc = sc_tb_quantile(*m, time, surv, *p, out); }
void SX_tb_test(int* n, double* L, double* R, int* g, int* G, int* max_iter, double* tol, double* stat, int* df, double* p,
                double* size, double* U, double* sd, int* rc) {
  *rc = sc_tb_test(*n, L, R, g, *G, *max_iter, *tol, stat, df, p, size, U, sd);
}
void SX_tb_bootstrap(int* n, double* L, double* R, int* ng, double* grid, int* B, int* max_iter, double* tol,
                     double* out, int* rc) {
  RStream rs;
  *rc = sc_tb_bootstrap(*n, L, R, *ng, grid, *B, 0, *max_iter, *tol, out);
}
void SX_grid_boot_limits(int* ng, double* grid, int* B, double* draws, int* nt, double* times, double* level, int* min_n,
                         double* lower, double* upper, int* rc) {
  *rc = sc_grid_boot_limits(*ng, grid, *B, draws, *nt, times, *level, *min_n, lower, upper);
}
void SX_column_limits(int* B, int* m, double* draws, double* level, int* min_n, double* lower, double* upper, int* rc) {
  *rc = sc_column_limits(*B, *m, draws, *level, *min_n, lower, upper);
}
void SX_tb_median_limits(int* ng, double* grid, int* B, double* draws, double* level, double* lower, double* upper, int* rc) {
  *rc = sc_tb_median_limits(*ng, grid, *B, draws, *level, lower, upper);
}
void SX_approx_rows(int* nr, int* nx, double* x, double* Y, int* nout, double* xout, double* out, int* rc) {
  *rc = sc_approx_rows(*nr, *nx, x, Y, *nout, xout, out);
}

void SX_measures(int* n, double* time, int* cause, int* interest, int* nt, double* t, double* risk, int* groups, double* brier,
                 double* brier_null, double* ipa, double* auc, double* calib, int* rc) {
  *rc = sc_measures(*n, time, cause, *interest, *nt, t, risk, *groups, brier, brier_null, ipa, auc, calib);
}
void SX_integrate_grid(int* m, double* grid, double* y, double* out, int* rc) { *rc = sc_integrate_grid(*m, grid, y, out); }
void SX_np_curve(int* n, double* time, int* status, int* cause, int* nt, double* at, double* out, int* rc) {
  *rc = sc_np_curve(*n, time, status, *cause, *nt, at, out);
}

void SX_rank_one_se(int* M, int* folds, double* loss, double* fold_loss, int* rank, double* diff, double* se, int* equivalent,
                    int* rc) {
  *rc = sc_rank_one_se(*M, *folds, loss, fold_loss, rank, diff, se, equivalent);
}

void SX_ph_test(int* n, int* p, double* time, int* status, double* X, double* beta, double* vcov, int* has_w, double* w,
                int* has_entry, double* entry, int* has_strata, int* strata, int* transform, int* term_of, int* nterm, double* chi_col,
                double* p_col, double* chi_term, double* df_term, double* p_term, double* chi_glob, double* p_glob, int* nd, int* rows,
                double* ev_time, double* g, double* scaled, int* rc) {
  *rc = sc_ph_test(*n, *p, time, status, X, beta, vcov, *has_w ? w : nullptr, *has_entry ? entry : nullptr,
                   *has_strata ? strata : nullptr, *transform, term_of, *nterm, chi_col, p_col, chi_term, df_term, p_term, chi_glob,
                   p_glob, nd, rows, ev_time, g, scaled);
}

void SX_p_adjust_holm(int* n, double* p, double* out, int* rc) { *rc = sc_p_adjust_holm(*n, p, out); }

void SX_te_split_size(int* n, double* start, double* stop, int* m, double* et, double* out, int* rc) {
  *rc = sc_te_split_size(*n, start, stop, *m, et, out);
}
void SX_te_split(int* n, double* start, double* stop, int* status, int* m, double* et, int* cap, int* idx, double* s0, double* s1,
                 int* st, int* nrow, int* rc) {
  *rc = sc_te_split(*n, start, stop, status, *m, et, *cap, idx, s0, s1, st, nrow);
}
void SX_te_basis(int* form, int* nt, double* t, double* tmin, double* tm, int* nbreak, double* breaks, int* nknot, double* knots,
                 double* bl, double* br, double* out, int* rc) {
  *rc = sc_te_basis(*form, *nt, t, *tmin, *tm, *nbreak, breaks, *nknot, knots, *bl, *br, out);
}
void SX_te_effect(int* nt, int* q, double* B, int* P, double* b, double* V, int* ix, double* zq, double* est, double* se, double* hr,
                  double* lower, double* upper, double* p, int* rc) {
  *rc = sc_te_effect(*nt, *q, B, *P, b, V, ix, *zq, est, se, hr, lower, upper, p);
}
void SX_cum_predict(int* mode, int* nr, int* m, double* tj, int* K, double* dH, int* cause, int* product_limit, int* nt, double* times,
                    double* out, int* rc) {
  *rc = sc_cum_predict(*mode, *nr, *m, tj, *K, dH, *cause, *product_limit, *nt, times, out);
}


void SX_fg_expand(int* n, double* time, int* status, int* has_entry, double* entry, int* cause, int* cap, int* id, double* start,
                  double* stop, int* st, double* w, int* nrow, int* rc) {
  *rc = sc_fg_expand(*n, time, status, *has_entry ? entry : nullptr, *cause, *cap, id, start, stop, st, w, nrow);
}
void SX_timefix(int* n, double* v, double* tol, int* changed, int* rc) { *rc = sc_timefix(*n, v, *tol, changed); }
void SX_sandwich(int* n, int* p, double* U, int* has_cl, int* cl, int* nc, double* V, double* out, int* rc) {
  *rc = sc_sandwich(*n, *p, U, *has_cl ? cl : nullptr, *nc, V, out);
}


void SX_par_residuals(int* dist, int* aft, int* nk, double* knots, int* p, double* theta, int* n, double* X, double* time,
                      int* has_entry, double* entry, int* status, int* type, double* out, int* rc) {
  *rc = sc_par_residuals(*dist, *aft, *nk, knots, *p, theta, *n, X, time, *has_entry ? entry : nullptr, status, *type, out);
}
void SX_cure_residuals(int* dist, int* aft, int* nk, double* knots, int* pl, int* pc, int* mixture, double* theta, int* n,
                       double* Xl, double* Xc, double* time, int* has_entry, double* entry, int* status, int* type, double* out, int* rc) {
  *rc = sc_cure_residuals(*dist, *aft, *nk, knots, *pl, *pc, *mixture, theta, *n, Xl, Xc, time, *has_entry ? entry : nullptr, status,
                          *type, out);
}


void SX_cov_complete(int* nb, int* p, double* X, double* out, int* nok, int* rc) { *rc = sc_cov_complete(*nb, *p, X, out, nok); }


void SX_ns_smooth(int* n, double* x, double* y, int* df, int* ng, double* grid, double* z, double* fit, double* lower, double* upper,
                  int* rc) {
  *rc = sc_ns_smooth(*n, x, y, *df, *ng, grid, *z, fit, lower, upper);
}


void SX_subdist_cif(int* nr, int* p, double* X, double* beta, int* nt, double* H, double* out, int* rc) {
  *rc = sc_subdist_cif(*nr, *p, X, beta, *nt, H, out);
}
void SX_log_limits(int* n, double* est, double* se, double* z, double* lower, double* upper, int* rc) {
  *rc = sc_log_limits(*n, est, se, *z, lower, upper);
}


void SX_elementwise(int* code, int* n, double* x, int* has_y, double* y, double* out, int* rc) {
  *rc = sc_elementwise(*code, *n, x, *has_y ? y : nullptr, out);
}
void SX_reduce(int* code, int* n, double* x, int* na_rm, double* out, int* rc) { *rc = sc_reduce(*code, *n, x, *na_rm, out); }
void SX_means(int* nr, int* nc, double* M, int* by_row, int* na_rm, double* out, int* rc) { *rc = sc_means(*nr, *nc, M, *by_row, *na_rm, out); }
void SX_group_means(int* n, int* g, double* x, int* G, double* out, int* rc) { *rc = sc_group_means(*n, g, x, *G, out); }
void SX_cumsum_rows(int* nr, int* nc, double* M, double* cap, double* out, int* rc) { *rc = sc_cumsum_rows(*nr, *nc, M, *cap, out); }
void SX_linpred(int* nr, int* p, double* X, double* beta, double* out, int* rc) { *rc = sc_linpred(*nr, *p, X, beta, out); }
void SX_ph_predict(int* mode, int* nr, int* p, double* X, double* beta, int* nt, double* H, int* by_row, double* out, int* rc) {
  *rc = sc_ph_predict(*mode, *nr, *p, X, beta, *nt, H, *by_row, out);
}
void SX_np_curve_inc(int* m, double* d, double* Y, int* mode, double* out, int* rc) { *rc = sc_np_curve_inc(*m, d, Y, *mode, out); }
void SX_info_criteria(int* n, double* ll, double* k, int* has_nobs, double* nobs, double* aic, double* bic, int* rc) {
  *rc = sc_info_criteria(*n, ll, k, *has_nobs ? nobs : nullptr, aic, bic);
}
void SX_lr_test(int* n, double* ll1, double* ll0, double* df, double* stat, double* p, int* rc) { *rc = sc_lr_test(*n, ll1, ll0, df, stat, p); }
void SX_wald_table(int* n, double* est, double* se, double* z, int* expo, double* lo, double* hi, double* zs, double* p, double* eest,
                   double* elo, double* ehi, int* rc) {
  *rc = sc_wald_table(*n, est, se, *z, expo, lo, hi, zs, p, eest, elo, ehi);
}
void SX_se_diag(int* p, double* V, double* out, int* rc) { *rc = sc_se_diag(*p, V, out); }
void SX_scaled_schoenfeld(int* nd, int* p, double* S, double* V, double* beta, double* out, int* rc) {
  *rc = sc_scaled_schoenfeld(*nd, *p, S, V, beta, out);
}
void SX_hazard_residuals(int* n, double* H, int* status, int* type, double* out, int* rc) { *rc = sc_hazard_residuals(*n, H, status, *type, out); }


void SX_grid(int* kind, double* from, double* to, double* by, int* n, int* cap, double* out, int* m, int* rc) {
  *rc = sc_grid(*kind, *from, *to, *by, *n, *cap, out, m);
}


void SX_hazard_increments(int* nb, double* btime, double* bcum, int* bstr, int* nr, int* rstr, int* m, double* tj, double* LP,
                          int* lp_by_time, double* out, int* rc) {
  *rc = sc_hazard_increments(*nb, btime, bcum, bstr, *nr, rstr, *m, tj, LP, *lp_by_time, out);
}


void SX_collinear(int* n, int* p, double* X, double* tol, int* out, int* ndep, int* rc) { *rc = sc_collinear(*n, *p, X, *tol, out, ndep); }

}  // extern "C"
