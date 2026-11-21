/* survcore: stable C API. Plain arrays in, plain arrays out.
 * Return value: 0 on success, negative on invalid input.
 *   -1 null pointer / bad dimension   -2 non-finite or negative input
 *   -3 no events of the requested type */
#pragma once
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double loglik_init, loglik, grad_maxabs;
  int iterations, converged, flags, n_events;
} sc_fit_info;

/* outputs must have room for n entries; *m receives the number used. */
int sc_kaplan_meier(int n, const double* time, const int* status, const double* weights,
                    double* out_time, double* out_nrisk, double* out_nevent,
                    double* out_surv, double* out_se, int* m);
/* Nelson-Aalen cumulative hazard and its standard error at the distinct event times. */
int sc_nelson_aalen(int n, const double* time, const int* status, const double* weights,
                    double* out_time, double* out_cumhaz, double* out_se, int* m);
/* Restricted mean survival time on [0, tau]. */
int sc_rmst(int n, const double* time, const int* status, double tau, double* mean, double* se);
/* out_cif and out_se (may be null): n x K row-major. */
int sc_aalen_johansen(int n, const double* time, const int* cause, int K,
                      double* out_time, double* out_surv, double* out_cif, double* out_se, int* m);

/* ties: 0 Breslow, 1 Efron. X row-major n x p. beta: in = start (or zeros), out = estimate.
 * var: p x p inverse information, filled with NaN when it is not positive definite. */
int sc_coxph(int n, int p, const double* time, const int* status, const double* X,
             const double* weights, int ties, int max_iter, double tol,
             double* beta, double* var, sc_fit_info* info);
/* Extended Cox model: entry (null: none) gives (entry, time] intervals for delayed entry and
 * time-varying covariates (Andersen & Gill 1982); strata (null: none) are codes 0, 1, ... with a
 * separate baseline hazard per stratum. Baseline output is ordered by stratum, then time;
 * out_strata receives the stratum of every entry. sc_coxph_score gives the n x p score residuals
 * (row-major) for the robust variance of Lin & Wei (1989). */
int sc_coxph_ext(int n, int p, const double* entry, const double* time, const int* status, const int* strata,
                 const double* X, const double* weights, int ties, int max_iter, double tol,
                 double* beta, double* var, sc_fit_info* info);
int sc_coxph_ext_baseline(int n, int p, const double* entry, const double* time, const int* status, const int* strata,
                          const double* X, const double* weights, int ties, const double* beta, const double* zref,
                          double* out_time, double* out_cumhaz, int* out_strata, int* m);
/* Martingale residuals (Efron or Breslow form, as survival) and, for rows with an event, the
 * risk-set covariate means of their Schoenfeld residuals (n x p row-major, NaN elsewhere). */
int sc_coxph_resid(int n, int p, const double* entry, const double* time, const int* status, const int* strata,
                   const double* X, const double* weights, int ties, const double* beta, double* mart, double* xbar);
int sc_coxph_score(int n, int p, const double* entry, const double* time, const int* status, const int* strata,
                   const double* X, const double* weights, int ties, const double* beta, double* U);
int sc_coxph_eval(int n, int p, const double* time, const int* status, const double* X,
                  const double* weights, int ties, const double* beta,
                  double* loglik, double* grad, double* info_matrix);

/* cause: 0 censored, `interest` the event of interest, other values competing.
 * var_model: inverse information. var_robust: sandwich variance of Fine & Gray (1999, Sec. 4).
 * Either may be null. */
int sc_finegray(int n, int p, const double* time, const int* cause, const double* X,
                int interest, int max_iter, double tol,
                double* beta, double* var_model, double* var_robust, sc_fit_info* info);
int sc_finegray_eval(int n, int p, const double* time, const int* cause, const double* X,
                     int interest, const double* beta,
                     double* loglik, double* grad, double* info_matrix);


/* Baseline cumulative (subdistribution) hazard at covariate value zref (null = zeros),
 * at the distinct event times. Outputs need room for n entries. */
int sc_coxph_baseline(int n, int p, const double* time, const int* status, const double* X,
                      const double* weights, int ties, const double* beta, const double* zref,
                      double* out_time, double* out_cumhaz, int* m);
int sc_finegray_baseline(int n, int p, const double* time, const int* cause, const double* X,
                         int interest, const double* beta, const double* zref,
                         double* out_time, double* out_cumhaz, int* m);
/* out: n_new x m row-major. mode 0: survival exp(-H e^lp); mode 1: CIF 1 - exp(-H e^lp). */
int sc_predict_from_cumhaz(int m, const double* cumhaz, int n_new, const double* lp, int mode, double* out);

/* Cumulative incidence from K cause-specific Cox models (causes coded 1..K).
 * betas: K x p row-major. product_limit 1: S(t) = prod(1 - sum_k dH_k); 0: exp(-sum_k H_k).
 * out_time needs n entries; out_cif needs n_new * n * K entries, laid out [new][time][cause]. */
int sc_csc_cif(int n, int p, const double* time, const int* cause, const double* X, int K,
               const double* betas, int ties, int n_new, const double* Xnew, int product_limit,
               double* out_time, double* out_cif, int* m);

/* K-sample tests. group codes 0..K-1; strata may be null.
 * out_score, out_var (K x K) and, for the log-rank test, out_obs / out_exp may be null. */
int sc_logrank(int n, const double* time, const int* status, const int* group, int K,
               const int* strata, double rho, double* stat, int* df, double* p_value,
               double* out_obs, double* out_exp, double* out_var);
int sc_gray_test(int n, const double* time, const int* cause, const int* group, int K,
                 const int* strata, int interest, double rho, double* stat, int* df,
                 double* p_value, double* out_score, double* out_var);
double sc_chisq_upper(double x, double df);

/* ---------------- Phase 2: direct parametric CIF models ----------------
 * cause: 0 censored, 1..K.
 * Link (Jeong & Fine 2007): alpha_free[k] = 1: alpha_k estimated (alpha_value[k] = start value);
 * alpha_free[k] = 0: alpha_k fixed at alpha_value[k] (0 = proportional hazards, 1 = proportional odds).
 * Baseline per cause (null = all Gompertz):
 *   0 Gompertz (Jeong & Fine 2006)            parameters rho, log tau
 *   1 three-parameter Gompertz (Haile et al. 2016)        rho, log tau, zeta
 *   2 Weibull                                  log lambda, log kappa
 *   3 mixture Weibull, 4 mixture log-logistic  logit p, log lambda, log kappa
 *   5 mixture log-normal                       logit p, mu, log sigma
 *   6 modified three-parameter logistic (Cheng 2009)      log b, c, logit p
 *   7 restricted cubic spline in log time (Royston & Parmar 2002; Mozumder et al. 2018)   gamma_0, gamma_1, ...
 * Spline knots: n_knots[k] knots of cause k (boundary knots first and last, in log time), all causes
 * concatenated in `knots`; n_knots[k] = 0 gives the 1-d.f. spline. n_knots / knots may be null.
 * Baselines 2-5 and 7 need strictly positive times.
 * Parameter vector theta, per cause in order: baseline parameters, [alpha_k if free], beta_k (p values).
 * constraint 0: the model as published; the sum of the fitted CIFs is only reported (default).
 * constraint 1: steps leaving { sum_k F_k(t_max; z_i) < 1 for all i } are rejected.
 * init: start values (length nparam) or null for the data-driven start.
 * diag (length 5, may be null): number of subjects with sum_k F_k(t_max) > 1, the largest such sum,
 * the same two at t = infinity, and the earliest time at which a subject's sum crosses 1 (NaN if none).
 * info->flags additionally uses 16 (sum exceeds 1 within the data range), 32 (exceeds 1 at infinity),
 * 64 (constrained fit stopped on the boundary). */
int sc_directcif_nparam(int p, int K, const int* alpha_free, const int* baseline, const int* n_knots);
/* Default knots of cause k (1..K) for a spline with df degrees of freedom: out_knots needs df + 1 entries. */
int sc_directcif_default_knots(int n, const double* time, const int* cause, int k, int df, double* out_knots, int* n_knots);
const char* sc_directcif_base_param_name(int baseline, int j);
int sc_directcif_fit(int n, int p, int K, const double* time, const int* cause, const double* X,
                     const int* alpha_free, const double* alpha_value, const int* baseline,
                     const int* n_knots, const double* knots, int constraint,
                     const double* init, int max_iter, double tol,
                     double* theta, double* var, sc_fit_info* info, double* diag);
int sc_directcif_eval(int n, int p, int K, const double* time, const int* cause, const double* X,
                      const int* alpha_free, const double* alpha_value, const int* baseline,
                      const int* n_knots, const double* knots, const double* theta,
                      double* loglik, double* grad, double* info_matrix);
/* CIF of cause `k` (1..K) for n_new covariate rows at nt times (a time may be +infinity: the
 * plateau; 1 for baselines without a plateau). out_cif, out_se: n_new x nt row-major; out_se
 * (delta method) needs var, both may be null. */
int sc_directcif_predict(int p, int K, const int* alpha_free, const double* alpha_value, const int* baseline,
                         const int* n_knots, const double* knots,
                         const double* theta, const double* var, int k, int n_new, const double* Xnew,
                         int nt, const double* times, double* out_cif, double* out_se);

/* ---------------- Phase 2: prediction performance at horizon t ----------------
 * pred[i] = predicted risk of `interest` by time t for subject i. */
int sc_brier(int n, const double* time, const int* cause, int interest, const double* pred, double t,
             double* brier, double* brier_null, double* ipa);
int sc_auc(int n, const double* time, const int* cause, int interest, const double* pred, double t, double* auc);
/* outputs need room for `groups` entries; *m receives the number of groups used. */
int sc_calibration(int n, const double* time, const int* cause, int interest, const double* pred, double t,
                   int groups, double* out_mean_pred, double* out_observed, double* out_size, double* error, int* m);

/* ---------------- Phase 2: model comparison ----------------
 * n_direct direct models, model j with baseline[j] (same family for every cause), link[j]
 * (0 proportional hazards, 1 proportional odds, 2 generalized odds rate with alpha estimated) and,
 * for spline baselines, spline_df[j] (null = 3); then Fine-Gray; then cause-specific Cox.
 * folds < 2: apparent performance only. With M = n_direct + 2:
 * out_metrics: M x nt x 4 x 2, [model][time][metric: brier, ipa, auc, calibration error][0 apparent, 1 cross-validated].
 * out_fit: M x 5, [model][loglik, npar, aic, flags, folds not converged]. */
int sc_compare_models(int n, int p, int K, const double* time, const int* cause, const double* X,
                      int interest, int nt, const double* times, int folds, unsigned long long seed,
                      int n_direct, const int* baseline, const int* link, const int* spline_df,
                      double* out_metrics, double* out_fit);
int sc_cv_folds(int n, int folds, unsigned long long seed, int* out_fold);

/* ---------------- Phase 2b: further estimators of the direct model ----------------
 * Weighted likelihood for one cause ("parametric Fine-Gray"): the direct model F(t; z) of
 * sc_directcif_* for the cause `interest` alone; subjects failing from other causes enter through
 * the Kaplan-Meier estimate of the censoring distribution (see directcif_ext.hpp).
 * cause: 0 = censored, `interest`, any other positive value = competing event.
 * Parameters as one cause block of sc_directcif_fit (baseline, [alpha], beta).
 * grid > 0 coarsens the censoring distribution to at most `grid` support points (0 = exact).
 * var_robust: sandwich variance (one score per subject); var_model: inverse information. */
int sc_subdist_fit(int n, int p, const double* time, const int* cause, const double* X, int interest,
                   int alpha_free, double alpha_value, int baseline, int n_knots, const double* knots, int grid,
                   const double* init, int max_iter, double tol,
                   double* theta, double* var_robust, double* var_model, sc_fit_info* info, int* n_support);
/* Sandwich variance of the weighted fit that includes the effect of estimating the censoring
 * distribution (exact influence of each subject, the Kaplan-Meier estimate recomputed):
 * var = var_model M var_model with M = sum_i (s_i + psi_i)(s_i + psi_i)'. Computed on the exact
 * censoring distribution. rc -4: outside the support. */
int sc_subdist_var_ipcw(int n, int p, const double* time, const int* cause, const double* X, int interest,
                        int alpha_free, double alpha_value, int baseline, int n_knots, const double* knots,
                        const double* theta, const double* var_model, double* var);
int sc_subdist_eval(int n, int p, const double* time, const int* cause, const double* X, int interest,
                    int alpha_free, double alpha_value, int baseline, int n_knots, const double* knots, int grid,
                    const double* theta, double* loglik, double* grad, double* info_matrix);

/* Constrained model for two causes, Shi, Cheng & Jeong (2013): cause 1 as in sc_directcif_* (a
 * baseline with a plateau: not Weibull, not spline), F_2(t; z) = {1 - F_1(inf; z)} G(t).
 * family (G): 0 modified logistic (log b, c), 1 Weibull, 2 log-logistic (log lambda, log kappa),
 * 3 log-normal (mu, log sigma). cause: 0, 1, 2. Parameters: cause 1 block, then the two of G. */
/* Constrained model for K >= 2 causes (extension of Shi, Cheng & Jeong 2013): cause 1 as in
 * sc_constrained_fit; F_k(t; z) = {1 - F_1(inf; z)} pi_k(z) G_k(t) for k = 2..K, with pi_k a
 * multinomial logit (cause 2 the reference) and G_k of family[k-2]. cause: 0 censored, 1..K.
 * Parameters: cause 1 block, (a_k, b_k) for k = 2..K, (gamma_k, delta_k) for k = 3..K.
 * K = 2 gives the model of sc_constrained_fit. predict: average = 1 gives the mean over the rows
 * of Xnew (out_cif and out_se of length nt). */
int sc_constrainedK_fit(int n, int p, const double* time, const int* cause, const double* X, int K,
                        int alpha_free, double alpha_value, int baseline, const int* family,
                        const double* init, int max_iter, double tol,
                        double* theta, double* var, sc_fit_info* info);
int sc_constrainedK_eval(int n, int p, const double* time, const int* cause, const double* X, int K,
                         int alpha_free, double alpha_value, int baseline, const int* family,
                         const double* theta, double* loglik, double* grad, double* info_matrix);
int sc_constrainedK_predict(int p, int K, int alpha_free, double alpha_value, int baseline, const int* family,
                            const double* theta, const double* var, int k, int average, int n_new, const double* Xnew,
                            int nt, const double* times, double* out_cif, double* out_se);
int sc_constrained_fit(int n, int p, const double* time, const int* cause, const double* X,
                       int alpha_free, double alpha_value, int baseline, int family,
                       const double* init, int max_iter, double tol,
                       double* theta, double* var, sc_fit_info* info);
int sc_constrained_eval(int n, int p, const double* time, const int* cause, const double* X,
                        int alpha_free, double alpha_value, int baseline, int family,
                        const double* theta, double* loglik, double* grad, double* info_matrix);
int sc_constrained_predict(int p, int alpha_free, double alpha_value, int baseline, int family,
                           const double* theta, const double* var, int k, int n_new, const double* Xnew,
                           int nt, const double* times, double* out_cif, double* out_se);

/* sc_compare_models with the estimator of every direct model: estimation[j] 0 joint, 1 weighted,
 * 2 constrained (two causes, baseline with a plateau) with conditional family cond[j] (see
 * sc_constrained_fit). Null estimation = all joint. AIC is reported for joint and constrained
 * models (full likelihoods of the same data), not for weighted models. */
int sc_compare_models2(int n, int p, int K, const double* time, const int* cause, const double* X,
                       int interest, int nt, const double* times, int folds, unsigned long long seed,
                       int n_direct, const int* baseline, const int* link, const int* spline_df,
                       const int* estimation, const int* cond, double* out_metrics, double* out_fit);

/* ---------------- curves for plotting ----------------
 * Mean over the rows of Xnew of the CIF of cause k, with the delta-method standard error of the mean
 * (gradient averaged over the rows). Arguments as in sc_directcif_predict / sc_constrained_predict;
 * outputs have nt entries. */
int sc_directcif_predict_mean(int p, int K, const int* alpha_free, const double* alpha_value, const int* baseline,
                              const int* n_knots, const double* knots, const double* theta, const double* var, int k,
                              int n_new, const double* Xnew, int nt, const double* times, double* out_cif, double* out_se);
int sc_constrained_predict_mean(int p, int alpha_free, double alpha_value, int baseline, int family, const double* theta,
                                const double* var, int k, int n_new, const double* Xnew, int nt, const double* times,
                                double* out_cif, double* out_se);
/* Nonparametric bootstrap (Efron & Tibshirani 1993) of the predicted CIF of cause `interest` for the Fine-Gray
 * model (model 0) or cause-specific Cox models (model 1; ties 1 Efron / 0 Breslow, product_limit as in
 * sc_csc_cif). adjust = 0: models without covariates (coefficients fixed at zero). average = 1: the mean
 * over the rows of Xnew (outputs nt entries), else n_new x nt. B resamples of the subjects drawn with
 * splitmix64(seed); out_est from the full data; out_se the bootstrap standard deviation; out_lower and
 * out_upper percentile limits at `level`; *n_ok resamples that could be fitted. Returns -4 when the
 * model cannot be fitted to the full data. */
int sc_cif_bootstrap(int model, int n, int p, const double* time, const int* cause, const double* X, int K, int interest,
                     int ties, int product_limit, int adjust, int n_new, const double* Xnew, int average, int nt,
                     const double* times, int B, unsigned long long seed, double level,
                     double* out_est, double* out_se, double* out_lower, double* out_upper, int* n_ok);

/* Parametric regression (survcore/parametric.hpp). dist: 0 exponential, 1 Weibull, 2 Gompertz,
 * 3 lognormal, 4 loglogistic, 6 Royston-Parmar spline (knots on the log-time scale, boundary knots
 * included), 7 piecewise exponential (knots = cut points); aft: 0 proportional hazards, 1 AFT. entry may be
 * null (no delayed entry); times must be positive and entry < time. theta: 1 + p (+1 shape)
 * parameters; init may be null. */
/* Parametric cure models: mixture (mixture = 1) or non-mixture (promotion time) with the cured
 * fraction expit(gamma'z); latency as in sc_par_fit. See survcore/cure.hpp. */
int sc_cure_fit(int n, int pl, int pc, int dist, int aft, int nk, const double* knots, int mixture,
                const double* entry, const double* time, const int* status, const double* Xl, const double* Xc,
                const double* init, int max_iter, double tol, double* theta, double* var, sc_fit_info* info);
int sc_cure_predict(int dist, int aft, int nk, const double* knots, int pl, int pc, int mixture, const double* theta,
                    const double* var, int n_new, const double* Xl, const double* Xc, int nt, const double* times,
                    int type, int average, double* out, double* out_se);
int sc_par_fit(int n, int p, int dist, int aft, int nk, const double* knots, const double* entry, const double* time,
               const int* status, const double* X, const double* init, int max_iter, double tol, double* theta,
               double* var, sc_fit_info* info);
/* type: 0 survival, 1 hazard, 2 cumulative hazard, 3 quantile (times = probabilities), 4 restricted
 * mean (times = tau). average = 1: mean over the rows of Xnew (nt outputs), else n_new x nt row-major.
 * out_se by the delta method; null var gives NaN. */
int sc_par_predict(int dist, int aft, int nk, const double* knots, int p, const double* theta, const double* var,
                   int n_new, const double* Xnew, int nt, const double* times, int type, int average,
                   double* out, double* out_se);
/* CIF of cause k (1..K) or overall survival (k = 0) from K parametric cause-specific hazards. */
int sc_parcsc_predict(int K, const int* dist, const int* aft, const int* nk, const double* knots, int p,
                      const double* theta, const double* var, int k, int n_new, const double* Xnew, int nt,
                      const double* times, double tmax, int average, double* out, double* out_se);
/* Kernel-smoothed hazard (survcore/hazard.hpp); bw <= 0: least-squares cross-validation. */
int sc_kernel_hazard(int m, const double* et, const double* d, const double* Y, double bw, int boundary,
                     int nx, const double* x, double* out_h, double* out_se, double* bw_used);

#ifdef __cplusplus
}
#endif
