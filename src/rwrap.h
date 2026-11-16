/* SurvX: entry points for R's .C interface. Every argument is a pointer; matrices are passed
 * row-major (the R code transposes). rc receives the return code of the C API. */
#pragma once
#ifdef __cplusplus
extern "C" {
#endif
void SX_km(int* n, double* time, int* status, double* out_time, double* out_nrisk, double* out_nevent, double* out_surv, double* out_se, int* m, int* rc);
void SX_na(int* n, double* time, int* status, double* out_time, double* out_cumhaz, double* out_se, int* m, int* rc);
void SX_rmst(int* n, double* time, int* status, double* tau, double* mean, double* se, int* rc);
void SX_aj(int* n, double* time, int* cause, int* K, double* out_time, double* out_surv, double* out_cif, double* out_se, int* m, int* rc);
void SX_cox(int* n, int* p, double* time, int* status, double* X, int* has_w, double* weights, int* ties, int* max_iter, double* tol, double* beta, double* var, double* info_d, int* info_i, int* rc);
void SX_cox_baseline(int* n, int* p, double* time, int* status, double* X, int* has_w, double* weights, int* ties, double* beta, double* out_time, double* out_cumhaz, int* m, int* rc);
void SX_fg(int* n, int* p, double* time, int* cause, double* X, int* interest, int* max_iter, double* tol, double* beta, double* var_model, double* var_robust, double* info_d, int* info_i, int* rc);
void SX_fg_baseline(int* n, int* p, double* time, int* cause, double* X, int* interest, double* beta, double* out_time, double* out_cumhaz, int* m, int* rc);
void SX_csc_cif(int* n, int* p, double* time, int* cause, double* X, int* K, double* betas, int* ties, int* n_new, double* Xnew, int* product_limit, double* out_time, double* out_cif, int* m, int* rc);
void SX_logrank(int* n, double* time, int* status, int* group, int* K, int* has_strata, int* strata, double* rho, double* stat, int* df, double* p_value, double* out_obs, double* out_exp, double* out_var, int* rc);
void SX_gray(int* n, double* time, int* cause, int* group, int* K, int* has_strata, int* strata, int* interest, double* rho, double* stat, int* df, double* p_value, double* out_score, double* out_var, int* rc);
void SX_direct_nparam(int* p, int* K, int* alpha_free, int* baseline, int* n_knots, int* out);
void SX_direct_knots(int* n, double* time, int* cause, int* k, int* df, double* out_knots, int* n_knots, int* rc);
void SX_direct_fit(int* n, int* p, int* K, double* time, int* cause, double* X, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, int* constraint, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, double* diag, int* rc);
void SX_direct_eval(int* n, int* p, int* K, double* time, int* cause, double* X, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, double* theta, double* loglik, double* grad, double* info_matrix, int* rc);
void SX_direct_predict(int* p, int* K, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc);
void SX_brier(int* n, double* time, int* cause, int* interest, double* pred, double* t, double* brier, double* brier_null, double* ipa, int* rc);
void SX_auc(int* n, double* time, int* cause, int* interest, double* pred, double* t, double* auc, int* rc);
void SX_calibration(int* n, double* time, int* cause, int* interest, double* pred, double* t, int* groups, double* out_mean_pred, double* out_observed, double* out_size, double* error, int* m, int* rc);
void SX_compare(int* n, int* p, int* K, double* time, int* cause, double* X, int* interest, int* nt, double* times, int* folds, int* seed, int* n_direct, int* baseline, int* link, int* spline_df, int* estimation, int* cond, double* out_metrics, double* out_fit, int* rc);
void SX_subdist_fit(int* n, int* p, double* time, int* cause, double* X, int* interest, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, int* grid, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var_robust, double* var_model, double* info_d, int* info_i, int* n_support, int* rc);
void SX_constrained_fit(int* n, int* p, double* time, int* cause, double* X, int* alpha_free, double* alpha_value, int* baseline, int* family, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, int* rc);
void SX_constrained_predict(int* p, int* alpha_free, double* alpha_value, int* baseline, int* family, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc);
void SX_direct_predict_mean(int* p, int* K, int* alpha_free, double* alpha_value, int* baseline, int* n_knots, double* knots, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc);
void SX_constrained_predict_mean(int* p, int* alpha_free, double* alpha_value, int* baseline, int* family, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* out_cif, double* out_se, int* rc);
void SX_cif_bootstrap(int* model, int* n, int* p, double* time, int* cause, double* X, int* K, int* interest, int* ties, int* product_limit, int* adjust, int* n_new, double* Xnew, int* average, int* nt, double* times, int* B, int* seed, double* level, double* out_est, double* out_se, double* out_lower, double* out_upper, int* n_ok, int* rc);
void SX_par_fit(int* n, int* p, int* dist, int* aft, int* has_entry, double* entry, double* time, int* status, double* X, int* has_init, double* init, int* max_iter, double* tol, double* theta, double* var, double* info_d, int* info_i, int* rc);
void SX_par_predict(int* dist, int* aft, int* p, double* theta, int* has_var, double* var, int* n_new, double* Xnew, int* nt, double* times, int* type, int* average, double* out, double* out_se, int* rc);
void SX_parcsc_predict(int* K, int* dist, int* aft, int* p, double* theta, int* has_var, double* var, int* k, int* n_new, double* Xnew, int* nt, double* times, double* tmax, int* average, double* out, double* out_se, int* rc);
void SX_kernel_hazard(int* m, double* et, double* d, double* Y, double* bw, int* boundary, int* nx, double* x, double* out_h, double* out_se, double* bw_used, int* rc);
#ifdef __cplusplus
}
#endif
