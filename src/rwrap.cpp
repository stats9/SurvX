// SurvX: thin wrappers that adapt the C API of the core (survcore/survcore.h) to R's .C interface.
#include "rwrap.h"
#include "survcore/survcore.h"

namespace {
inline void put(const sc_fit_info& f, double* d, int* i) {
  d[0] = f.loglik_init; d[1] = f.loglik; d[2] = f.grad_maxabs;
  i[0] = f.iterations; i[1] = f.converged; i[2] = f.flags; i[3] = f.n_events;
}
}  // namespace

extern "C" {
void SX_km(int* n, double* time, int* status, double* out_time, double* out_nrisk, double* out_nevent, double* out_surv, double* out_se, int* m, int* rc) {
  *rc = sc_kaplan_meier(*n, time, status, nullptr, out_time, out_nrisk, out_nevent, out_surv, out_se, m);
}
void SX_na(int* n, double* time, int* status, double* out_time, double* out_cumhaz, double* out_se, int* m, int* rc) {
  *rc = sc_nelson_aalen(*n, time, status, nullptr, out_time, out_cumhaz, out_se, m);
}
void SX_rmst(int* n, double* time, int* status, double* tau, double* mean, double* se, int* rc) {
  *rc = sc_rmst(*n, time, status, *tau, mean, se);
}
void SX_aj(int* n, double* time, int* cause, int* K, double* out_time, double* out_surv, double* out_cif, double* out_se, int* m, int* rc) {
  *rc = sc_aalen_johansen(*n, time, cause, *K, out_time, out_surv, out_cif, out_se, m);
}
void SX_cox(int* n, int* p, double* time, int* status, double* X, int* has_w, double* weights, int* ties, int* max_iter, double* tol, double* beta, double* var, double* info_d, int* info_i, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_coxph(*n, *p, time, status, X, *has_w ? weights : nullptr, *ties, *max_iter, *tol, beta, var, &f);
  put(f, info_d, info_i);
}
void SX_cox_baseline(int* n, int* p, double* time, int* status, double* X, int* has_w, double* weights, int* ties, double* beta, double* out_time, double* out_cumhaz, int* m, int* rc) {
  *rc = sc_coxph_baseline(*n, *p, time, status, X, *has_w ? weights : nullptr, *ties, beta, nullptr, out_time, out_cumhaz, m);
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
void SX_direct_nparam(int* p, int* K, int* alpha_free, int* baseline, int* n_knots, int* out) {
  *out = sc_directcif_nparam(*p, *K, alpha_free, baseline, n_knots);
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
void SX_direct_eval(int* n, int* p, int* K, double* time, int* cause, double* X, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, double* theta, double* loglik, double* grad, double* info_matrix, int* rc) {
  *rc = sc_directcif_eval(*n, *p, *K, time, cause, X, alpha_free, alpha_value, baseline, n_knots, knots, theta, loglik, grad, info_matrix);
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
void SX_compare(int* n, int* p, int* K, double* time, int* cause, double* X, int* interest, int* nt, double* times, int* folds, int* seed, int* n_direct, int* baseline, int* link, int* spline_df, int* estimation, int* cond, double* out_metrics, double* out_fit, int* rc) {
  *rc = sc_compare_models2(*n, *p, *K, time, cause, X, *interest, *nt, times, *folds, static_cast<unsigned long long>(*seed < 0 ? -*seed : *seed),
                           *n_direct, baseline, link, spline_df, estimation, cond, out_metrics, out_fit);
}
void SX_subdist_fit(int* n, int* p, double* time, int* cause, double* X, int* interest, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, int* grid, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var_robust, double* var_model, double* info_d, int* info_i, int* n_support, int* rc) {
  sc_fit_info f = sc_fit_info();
  *rc = sc_subdist_fit(*n, *p, time, cause, X, *interest, *alpha_free, *alpha_value, *baseline, *n_knots, knots, *grid,
                       *has_init ? init : nullptr, *max_iter, *tol, theta, var_robust, var_model, &f, n_support);
  put(f, info_d, info_i);
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
void SX_direct_predict_mean(int* p, int* K, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc) {
  *rc = sc_directcif_predict_mean(*p, *K, alpha_free, alpha_value, baseline, n_knots, knots, theta, *has_var ? var : nullptr, *k, *n_new, Xnew, *nt, times,
                                  out_cif, *has_var ? out_se : nullptr);
}
void SX_constrained_predict_mean(int* p, int* alpha_free, double* alpha_value, int* baseline, int* family, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc) {
  *rc = sc_constrained_predict_mean(*p, *alpha_free, *alpha_value, *baseline, *family, theta, *has_var ? var : nullptr, *k, *n_new, Xnew, *nt, times,
                                    out_cif, *has_var ? out_se : nullptr);
}
void SX_cif_bootstrap(int* model, int* n, int* p, double* time, int* cause, double* X, int* K, int* interest, int* ties, int* product_limit, int* adjust, int* n_new, double* Xnew, int* average, int* nt, double* times, int* B, int* seed, double* level, double* out_est, double* out_se, double* out_lower, double* out_upper, int* n_ok, int* rc) {
  *rc = sc_cif_bootstrap(*model, *n, *p, time, cause, X, *K, *interest, *ties, *product_limit, *adjust, *n_new, Xnew, *average, *nt, times, *B,
                         static_cast<unsigned long long>(*seed < 0 ? -*seed : *seed), *level, out_est, out_se, out_lower, out_upper, n_ok);
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
void SX_cv_folds(int* n, int* folds, int* seed, int* out, int* rc) {
  *rc = sc_cv_folds(*n, *folds, static_cast<unsigned long long>(*seed < 0 ? -*seed : *seed), out);
}
void SX_kernel_hazard(int* m, double* et, double* d, double* Y, double* bw, int* boundary, int* nx, double* x, double* out_h, double* out_se, double* bw_used, int* rc) {
  *rc = sc_kernel_hazard(*m, et, d, Y, *bw, *boundary, *nx, x, out_h, out_se, bw_used);
}
}  // extern "C"
