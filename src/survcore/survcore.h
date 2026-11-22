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
/* Interval-censored data (survcore/directcif.hpp, DirectCIFModel::set_time2): an event of cause
 * cause_i in (time_i, time2_i]; time2_i = time_i an exact time, time_i = 0 left censoring; censored
 * subjects (cause 0) at time_i, time2 ignored. Other arguments as sc_directcif_fit. */
int sc_directcif_fit_ic(int n, int p, int K, const double* time, const double* time2, const int* cause, const double* X,
                        const int* alpha_free, const double* alpha_value, const int* baseline,
                        const int* n_knots, const double* knots, int constraint,
                        const double* init, int max_iter, double tol,
                        double* theta, double* var, sc_fit_info* info, double* diag);
int sc_directcif_eval_ic(int n, int p, int K, const double* time, const double* time2, const int* cause, const double* X,
                         const int* alpha_free, const double* alpha_value, const int* baseline,
                         const int* n_knots, const double* knots, const double* theta,
                         double* loglik, double* grad, double* info_matrix);
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
/* The same with covariates for the time distributions G_k: an accelerated failure time model,
 * G_k(t | w, v) = G_k0(t exp(-eta_k' w)) with the shape of G_k0 shifted by zeta_k' v (see
 * ConstrainedKModel in directcif_ext.hpp for the exact form per family). W: n x nw and V: n x nv,
 * row-major; wuse, vuse: (K - 1) x nw and (K - 1) x nv flags of 0/1, row k - 2 for cause k, the
 * columns used for that cause. Parameters: as above, then eta_k for k = 2..K (the flagged columns
 * in order), then zeta_k for k = 2..K. nw = nv = 0 gives sc_constrainedK_fit exactly. In
 * predict, Wnew and Vnew have n_new rows. */
int sc_constrainedK_fit_cov(int n, int p, const double* time, const int* cause, const double* X, int K,
                            int alpha_free, double alpha_value, int baseline, const int* family,
                            int nw, const double* W, const int* wuse, int nv, const double* V, const int* vuse,
                            const double* init, int max_iter, double tol,
                            double* theta, double* var, sc_fit_info* info);
int sc_constrainedK_eval_cov(int n, int p, const double* time, const int* cause, const double* X, int K,
                             int alpha_free, double alpha_value, int baseline, const int* family,
                             int nw, const double* W, const int* wuse, int nv, const double* V, const int* vuse,
                             const double* theta, double* loglik, double* grad, double* info_matrix);
int sc_constrainedK_predict_cov(int p, int K, int alpha_free, double alpha_value, int baseline, const int* family,
                                int nw, const int* wuse, int nv, const int* vuse,
                                const double* theta, const double* var, int k, int average, int n_new, const double* Xnew,
                                const double* Wnew, const double* Vnew,
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
/* Interval censoring (survcore/parametric.hpp, par_local_ic): type 0 right-censored at time, 1 exact
 * at time, 2 in (time, time2], 3 left-censored in (0, time2]. init is required (a start from the fit to
 * the midpoints). */
int sc_par_fit_ic(int n, int p, int dist, int aft, int nk, const double* knots, const double* time, const double* time2,
                  const int* type, const double* X, const double* init, int max_iter, double tol, double* theta,
                  double* var, sc_fit_info* info);
/* Shared gamma frailty (survcore/frailty.hpp). */
int sc_par_frailty(int n, int p, int dist, int aft, int nk, const double* knots, const double* time, const int* status,
                   const double* X, const int* cluster, int C, const double* init, int max_iter, double tol,
                   double* theta, double* var, sc_fit_info* info, double* shape, double* rate);
int sc_cox_frailty(int n, int p, const double* time, const int* status, const double* X, const int* cluster, int C,
                   const double* init, int max_iter, double tol, double* est, double* var, sc_fit_info* info,
                   int* K_out, double* tk, double* hk, double* shape, double* rate);
/* Turnbull NPMLE of a survival function from interval-censored data (survcore/interval.hpp). */
int sc_turnbull(int n, const double* L, const double* R, const double* w, int max_iter, double tol, int* n_int,
                double* lo, double* hi, double* mass, int* a, int* b, double* info_d, int* info_i);
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

/* ---------------- numerical summaries shared by the front ends (capi_support.cpp) ----------------
 * Distribution functions: standard normal (pnorm, qnorm), upper tail of the chi-square, gamma
 * quantile (shape, rate). sc_quantile: sample quantiles exactly as R's quantile() of type 1 or 7,
 * NA and NaN dropped. sc_cif_interval: limits of a cumulative incidence on the
 * complementary log-log scale (NaN where F is 0 or 1). sc_step_at: right-continuous step function
 * (jump times t increasing; `before` before the first jump). Frailty: half chi-square test of no
 * frailty; size, events and posterior mean and limits of every cluster (cluster codes 0..C-1);
 * cumulative hazards from Breslow jumps (out n x nt); survival in mode 0 marginal, 1 conditional,
 * 2 for the cluster of each row. sc_na_table: Nelson-Aalen estimate, standard error and limits. */
double sc_qnorm(double p);
double sc_pnorm(double x);
double sc_pchisq_upper(double x, double df);
double sc_qgamma(double p, double shape, double rate);
int sc_quantile(int n, const double* x, int np, const double* p, int type, double* out);
int sc_cif_interval(int n, const double* F, const double* se, double zq, double* lower, double* upper);
int sc_step_at(int m, const double* t, const double* v, int ns, const double* s, double before, double* out);
int sc_frailty_lrt(double ll, double ll0, int boundary, double* stat, double* p);
int sc_frailty_clusters(int n, const int* status, const int* cluster, int C, const double* shape, const double* rate,
                        int boundary, double level, double* size, double* events, double* mean, double* lower,
                        double* upper);
int sc_breslow_cumhaz(int nb, const double* bh_time, const double* bh_hazard, int n, int p, const double* X,
                      const double* beta, int nt, const double* times, double* out);
int sc_frailty_survival(int n, int nt, const double* H, double theta, int mode, const double* shape, const double* rate,
                        double* out);
int sc_na_table(int m, const double* d, const double* Y, double zq, double* H, double* se, double* lower,
                double* upper);
/* Kaplan-Meier limits on the log scale (NaN where S = 0); first time at which a curve is at or
 * below p (NaN if never); events and numbers at risk at the distinct event times of cause k. */
int sc_km_limits(int m, const double* S, const double* se, double zq, double* lower, double* upper);
int sc_first_below(int m, const double* t, const double* v, double p, double* out);
int sc_event_counts(int n, const double* time, const int* status, int k, const double* entry, double* et, double* d,
                    double* Y, int* m);
/* Adjusted risks (survcore/adjust.hpp). sc_mlogit: multinomial logistic model of a (codes 0..K-1)
 * on X (n x p, with an intercept column), P n x K row-major. sc_ipw_weights: stabilised weights,
 * trimmed at the trim and 1 - trim quantiles. sc_balance: standardised mean differences of the q
 * columns of Xc ((K - 1) x q, NaN where undefined) and weight diagnostics per level.
 * sc_weighted_cif: weighted Aalen-Johansen estimate of cause k (entry may be null).
 * sc_boot_summary: standard deviation and percentile limits of every column of nb x m draws.
 * sc_resample: B samples of 0..n-1 with replacement (B x n). sc_stratified_folds: folds 1..folds
 * balanced within strata. Both seeded (splitmix64), the same in every front end. */
int sc_mlogit(int n, int p, const double* X, const int* a, int K, double* P, int* iterations);
int sc_ipw_weights(int n, int K, const double* P, const int* a, double trim, double* w);
int sc_balance(int n, int q, const double* Xc, const int* a, int K, const double* w, const double* P, double* before,
               double* after, double* size, double* ess, double* max_weight, double* ps_min, double* ps_max);
int sc_weighted_cif(int n, const double* time, const int* status, const double* entry, const double* w, int k, int nt,
                    const double* tt, double* out);
int sc_boot_summary(int nb, int m, const double* draws, int finite_only, double level, double* sd, double* lower,
                    double* upper);
int sc_resample(unsigned long long seed, int n, int B, int* out);
int sc_stratified_folds(int n, const int* strata, int folds, unsigned long long seed, int* out);
/* Turnbull NPMLE summaries (survcore/interval_summary.hpp), from the innermost intervals and
 * masses of sc_turnbull. sc_tb_surv: P(T > t) (left = 1: P(T >= t)) with the identified bounds.
 * sc_tb_steps: the curve as a polyline (outputs 2 (2 J + 1)). sc_tb_quantile: first time the
 * curve drops to p. sc_tb_test: Finkelstein scores, permutational variance (group codes 0..G-1).
 * sc_tb_bootstrap: B x ng bootstrap curves (seeded). sc_grid_boot_limits, sc_column_limits,
 * sc_tb_median_limits: percentile limits from bootstrap curves. sc_approx_rows: rows interpolated
 * linearly at xout (as R's approx(rule = 1)). */
int sc_tb_surv(int J, const double* lo, const double* hi, const double* p, int nt, const double* t, int left, double* est,
               double* lower, double* upper);
int sc_tb_steps(int J, const double* lo, const double* hi, const double* p, double* time, double* surv, int* m);
int sc_tb_quantile(int m, const double* time, const double* surv, double p, double* out);
int sc_tb_test(int n, const double* L, const double* R, const int* g, int G, int max_iter, double tol, double* stat, int* df,
               double* p, double* size, double* U, double* sd);
int sc_tb_bootstrap(int n, const double* L, const double* R, int ng, const double* grid, int B, unsigned long long seed,
                    int max_iter, double tol, double* out);
int sc_grid_boot_limits(int ng, const double* grid, int B, const double* draws, int nt, const double* times, double level,
                        int min_n, double* lower, double* upper);
int sc_column_limits(int B, int m, const double* draws, double level, int min_n, double* lower, double* upper);
int sc_tb_median_limits(int ng, const double* grid, int B, const double* draws, double level, double* lower, double* upper);
int sc_approx_rows(int nr, int nx, const double* x, const double* Y, int nout, const double* xout, double* out);
/* Prediction measures of risks (n x nt) of cause `interest` at every time (sc_brier, sc_auc,
 * sc_calibration; NaN where undefined); trapezoidal mean of y over (0, grid]; non-parametric
 * probability of the event (1 - Kaplan-Meier) or cumulative incidence of `cause` at `at`. */
int sc_measures(int n, const double* time, const int* cause, int interest, int nt, const double* t, const double* risk,
                int groups, double* brier, double* brier_null, double* ipa, double* auc, double* calib);
int sc_integrate_grid(int m, const double* grid, const double* y, double* out);
int sc_np_curve(int n, const double* time, const int* status, int cause, int nt, const double* at, double* out);
/* Ranking of models on a cross-validated loss with the one-standard-error rule (rank 0 for models
 * without a finite loss); fold_loss M x folds. */
int sc_rank_one_se(int M, int folds, const double* loss, const double* fold_loss, int* rank, double* diff, double* se,
                   int* equivalent);
/* Score test of proportional hazards of a Cox model (Grambsch & Therneau 1994; survcore/phtest.hpp).
 * transform 0 identity, 1 log, 2 rank, 3 Kaplan-Meier, 4 weighted Kaplan-Meier of the event times. */
int sc_ph_test(int n, int p, const double* time, const int* status, const double* X, const double* beta, const double* vcov,
               const double* w, const double* entry, const int* strata, int transform, const int* term_of, int nterm,
               double* chi_col, double* p_col, double* chi_term, double* df_term, double* p_term, double* chi_glob,
               double* p_glob, int* nd, int* rows, double* ev_time, double* g, double* scaled);
/* Holm's adjustment of p-values (NaN left out), as R's p.adjust(method = "holm"). */
int sc_p_adjust_holm(int n, const double* p, double* out);
/* Time-varying effects (survcore/timeeffects.hpp). sc_te_split_size, sc_te_split: rows (start, stop]
 * cut at the sorted event times et strictly inside them (idx 0-based source rows). sc_te_basis:
 * built-in time functions, form 1 log(max(t, tmin/2)/tm), 2 t - tm, 3 steps 1{t > breaks_j},
 * 4 natural cubic spline in log(max(t, tmin/2)) (as splines::ns with knots and boundary knots on
 * the log scale); out nt x q. sc_te_effect: b[ix] combined with (1, B(t)), delta-method standard
 * error, hazard ratio, limits and p. sc_cum_predict: mode 0 survival, 1 cumulative incidence of one
 * hazard, 2 cumulative incidence of cause `cause` (0-based) from K cause-specific hazards, from
 * hazard increments dH (K x nr x m) on the event times tj. */
int sc_te_split_size(int n, const double* start, const double* stop, int m, const double* et, double* out);
int sc_te_split(int n, const double* start, const double* stop, const int* status, int m, const double* et, int cap, int* idx,
                double* s0, double* s1, int* st, int* nrow);
int sc_te_basis(int form, int nt, const double* t, double tmin, double tm, int nbreak, const double* breaks, int nknot,
                const double* knots, double bl, double br, double* out);
/* sc_te_effect: ix = (constant column, time-function columns), 0-based; ix[0] = -1 when there is no constant */
int sc_te_effect(int nt, int q, const double* B, int P, const double* b, const double* V, const int* ix, double zq, double* est,
                 double* se, double* hr, double* lower, double* upper, double* p);
int sc_cum_predict(int mode, int nr, int m, const double* tj, int K, const double* dH, int cause, int product_limit, int nt,
                   const double* times, double* out);

/* Fine-Gray counting-process data with delayed entry (Geskus 2011): call with cap = 0 for *nrow,
   then with arrays of that length; id 0-based; rows with weight 0 are left out */
int sc_fg_expand(int n, const double* time, const int* status, const double* entry, int cause, int cap, int* id,
                 double* start, double* stop, int* st, double* w, int* nrow);
/* merge time values closer than tol, absolutely or relative to the mean absolute time (in place) */
int sc_timefix(int n, double* v, double tol, int* changed);
/* robust (sandwich) variance V (sum_c U_c U_c') V; U n x p row-major; cl 0-based cluster codes or NULL */
int sc_sandwich(int n, int p, const double* U, const int* cl, int nc, const double* V, double* out);


/* Source of random indices. Every random draw of the core (bootstrap resamples, permutations,
   folds) is a uniform index in 0..n-1. With no source installed (the default) it comes from
   splitmix64 started at the seed argument of the call; fn, when installed, is called instead with
   n and must return an integer value in [0, n), and the seed arguments are then not used. The R
   front end installs R_unif_index around each call, so that set.seed() and RNGkind() govern the
   draws. Pass NULL to return to the default. Not thread-safe: install it on the calling thread. */
typedef double (*sc_index_fn)(double n);
void sc_set_index_source(sc_index_fn fn);


/* Residuals at the subjects' own times: type 0 Cox-Snell, 1 martingale, 2 deviance; entry NULL or
   delayed entry times (the cumulative hazard at entry is subtracted). Parametric and cure models. */
int sc_par_residuals(int dist, int aft, int nk, const double* knots, int p, const double* theta, int n, const double* X,
                     const double* time, const double* entry, const int* status, int type, double* out);
int sc_cure_residuals(int dist, int aft, int nk, const double* knots, int pl, int pc, int mixture, const double* theta, int n,
                      const double* Xl, const double* Xc, const double* time, const double* entry, const int* status, int type,
                      double* out);


/* sample covariance of the rows of X (nb x p, row-major) that have no NaN, as stats::cov of the
   complete cases; *nok: the number of rows used (NaN output below 2) */
int sc_cov_complete(int nb, int p, const double* X, double* out, int* nok);


/* least squares smoother y ~ 1 + ns(x, df), as stats::lm(y ~ splines::ns(x, df)), with the fit and
   pointwise limits fit +- z se at the grid; -4 when rank deficient or too few points */
int sc_ns_smooth(int n, const double* x, const double* y, int df, int ng, const double* grid, double z, double* fit, double* lower,
                 double* upper);


/* elementwise prediction formulas. sc_subdist_cif: Fine-Gray cumulative incidence
   1 - exp(-exp(x_i' beta) H_j), nr x nt row-major. sc_log_limits: limits est exp(-+ z se / est)
   of a positive quantity (hazard, cumulative hazard, quantile, RMST); NaN where undefined */
int sc_subdist_cif(int nr, int p, const double* X, const double* beta, int nt, const double* H, double* out);
int sc_log_limits(int n, const double* est, const double* se, double z, double* lower, double* upper);


/* Toolkit (survcore/toolkit.hpp): the small computations on results, so that no front end computes
   anything itself. sc_elementwise: out = f(x) or f(x, y) for code 1 exp, 2 log, 3 plogis, 4 logit,
   5 qnorm, 6 1 - x, 7 -log x, 8 x - y, 9 x / y, 10 log(x / y), 11 two-sided normal p of x / y,
   12 exp(x) sqrt(y), 13 interval midpoint (x left, y right), 14 -log(max(x, 1e-300)), 15 x y,
   16 sqrt x, 17 x + y, 18 qnorm(1 - x), 19 log((1 - x) / x), 20 exp(-x). sc_reduce: 0 sum, 1 mean,
   2 max |x|, 3 max, 4 min (as R; na_rm leaves out NaN). sc_means: column (or row) means.
   sc_group_means: means within groups. sc_cumsum_rows: running sums along rows, capped at cap >= 0.
   sc_linpred: X beta. sc_ph_predict: survival (mode 0) or cumulative incidence (1) from a cumulative
   baseline. sc_np_curve_inc: Nelson-Aalen (0), exp(-H) (1) or Kaplan-Meier (2) from counts.
   sc_info_criteria, sc_lr_test, sc_wald_table, sc_se_diag, sc_scaled_schoenfeld: as named. */
int sc_elementwise(int code, int n, const double* x, const double* y, double* out);
int sc_reduce(int code, int n, const double* x, int na_rm, double* out);
int sc_means(int nr, int nc, const double* M, int by_row, int na_rm, double* out);
int sc_group_means(int n, const int* g, const double* x, int G, double* out);
int sc_cumsum_rows(int nr, int nc, const double* M, double cap, double* out);
int sc_linpred(int nr, int p, const double* X, const double* beta, double* out);
int sc_ph_predict(int mode, int nr, int p, const double* X, const double* beta, int nt, const double* H, int by_row, double* out);
int sc_np_curve_inc(int m, const double* d, const double* Y, int mode, double* out);
int sc_info_criteria(int n, const double* ll, const double* k, const double* nobs, double* aic, double* bic);
int sc_lr_test(int n, const double* ll1, const double* ll0, const double* df, double* stat, double* p);
int sc_wald_table(int n, const double* est, const double* se, double z, const int* expo, double* lo, double* hi, double* zs,
                  double* p, double* eest, double* elo, double* ehi);
int sc_se_diag(int p, const double* V, double* out);
int sc_scaled_schoenfeld(int nd, int p, const double* S, const double* V, const double* beta, double* out);
/* residuals from cumulative hazards H (minus that at entry): type 0 Cox-Snell, 1 martingale, 2 deviance;
   with the martingale residuals given in place of H: 3 deviance, 4 Cox-Snell */
int sc_hazard_residuals(int n, const double* H, const int* status, int type, double* out);


/* grids as R: kind 0 seq(from, to, length.out = n), 1 seq(from, to, by), 2 pretty(c(from, to), n);
   call with cap = 0 for *m, then with an array of that length. sc_elementwise codes 21 |x - y|,
   22 signif(x, y), 23 round(x) (halves to even, as R) */
int sc_grid(int kind, double from, double to, double by, int n, int cap, double* out, int* m);


/* hazard increments dH0 exp(LP) of rows on a grid of times, from a stratified cumulative baseline */
int sc_hazard_increments(int nb, const double* btime, const double* bcum, const int* bstr, int nr, const int* rstr, int m,
                         const double* tj, const double* LP, int lp_by_time, double* out);

#ifdef __cplusplus
}
#endif
