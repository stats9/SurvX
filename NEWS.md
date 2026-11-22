# SurvX (development version)

## New methods

* **Interval-censored data** (the event is only known to lie between two visits):
  * New `Interval(left, right, cause)` outcome: interval-, left- and right-censored and exact times
    can be mixed in one data set.
  * New `turnbull()`: the Turnbull non-parametric maximum likelihood estimate of the survival curve
    (EM-ICM hybrid of the self-consistency and iterative convex minorant algorithms in compiled
    code, convergence checked on the Kuhn-Tucker conditions), by group, with bootstrap bands,
    medians with intervals, the range of survival values that the data leave unidentified within
    the innermost intervals, and a logrank-type test with Finkelstein scores and permutational
    variance. `print()`, `summary()`, `plot()`, `predict()`, `as.data.frame()`.
  * `parsurv()` fits every family (exponential, Weibull, Gompertz, log-normal, log-logistic,
    generalized gamma, Royston-Parmar spline, piecewise exponential; PH and AFT) to `Interval()`
    outcomes by the full likelihood with exact derivatives; `plot()` and `hazard_plot()` compare
    the fit with the Turnbull estimate.
  * `directcif()` (joint estimation) fits the direct parametric models for the cumulative
    incidence to interval-censored competing risks data with a known cause, using
    `F_k(R) - F_k(L)`.
* **Shared gamma frailty**: `+ frailty(id)` in the formula of `cox()` (Breslow baseline profiled
  out by EM, Newton on the profile likelihood of the coefficients and the frailty variance, whose
  curvature gives standard errors that account for the estimation of the variance) and of
  `parsurv()` (every family, marginal likelihood with exact derivatives). The output gives the
  conditional hazard ratios, the frailty variance with its interval, the likelihood ratio test of
  no frailty (half chi-square), Kendall's tau, the empirical Bayes frailty of every cluster, and
  `predict(fit, frailty = "marginal" | "conditional" | "cluster")`; `plot()` draws the cluster
  frailties or the marginal and conditional survival curves.
* **`time_effects()`**: for a `cox()`, `csc()` or `fg()` model, tests every term for a
  time-varying effect (Grambsch-Therneau score test, Holm-adjusted), compares forms of the change
  (constant, log-linear, a step at a change point chosen by profile likelihood, a natural spline
  in log time, or your own functions of time) by AIC or BIC and likelihood ratio tests, refits the
  model with the chosen forms on follow-up split at every event time, and reports the hazard ratio
  (cause-specific or subdistribution) at chosen times with a reading in plain words. `fix` imposes
  forms; `bootstrap` repeats the whole procedure for intervals and selection frequencies that
  include the choice. `plot()`, `predict()` (survival or cumulative incidence with the
  time-varying effects), `as.data.frame()` and an HTML `report()`.
* **`adjusted_risk()`**: adjusted absolute risks and causal contrasts between treatments, by
  standardisation of a fitted `cox()`, `fg()`, `csc()`, `directcif()`, `parsurv()`, `cure()` or
  `time_effects()` model (the g-formula) and/or by inverse probability of treatment weighting
  (weighted Kaplan-Meier or Aalen-Johansen with a logistic or multinomial treatment model), with
  risk differences and ratios at chosen times, curves, bootstrap intervals in which the models are
  refitted, balance (standardised mean differences) and weight diagnostics, a reading in plain
  words, `plot()` and an HTML `report()`.
* **`model_scan()`**: fits every model of the package that suits the data (one event type or
  competing risks; cure models when the Kaplan-Meier curve levels off) and compares them on
  criteria valid for the comparison: the cross-validated integrated Brier score, IPA, AUC and
  calibration error for all models, with the models within one standard error of the best
  grouped at rank 1; AIC only among the models with a full likelihood of the same data; and
  checks of every model (convergence, proportional hazards, incidences adding up to more than
  one, distance from the non-parametric estimate). `goal = "effect"` adds the effects of the
  covariates in every model side by side. `plot()` draws the ranking, the fit, the prediction
  error over time and the calibration; `report()` writes a self-contained HTML report with a
  methods paragraph and references.

## Improvements to existing functions

* `directcif(estimation = "weighted")`: new argument `variance`. The default, `"ipcw"`, is a
  sandwich variance that includes the effect of estimating the censoring distribution (the exact
  influence of each subject with the Kaplan-Meier estimate recomputed; the discrete form of the
  Fine-Gray correction term). `"robust"` gives the previous standard errors (censoring
  distribution treated as known), `"model"` the inverse information. The correction mostly
  shrinks the standard errors of the baseline parameters and of the predicted incidences; for
  the regression coefficients it is usually small. All three matrices are kept in the fit.
* `directcif(estimation = "constrained")` now takes any number of causes. With K >= 3 the
  cause of interest keeps the direct model of Shi, Cheng & Jeong (2013) and the remaining
  long-run risk is shared among the competing causes by a multinomial logit with covariates
  (`share_` coefficients, odds ratios of one competing cause against the first); the long-term
  incidences add up to one for every covariate value. `cond` may give one distribution per
  competing cause. K = 2 is unchanged.
* `ph_test()` and `model_scan()` point to `time_effects()` when a covariate fails the test of
  proportional hazards.
* Functions that need exact times (Kaplan-Meier, Cox, Fine-Gray, residuals, prediction metrics,
  `adjusted_risk()`) stop with a message that points to the functions for interval-censored data;
  `frailty()` in other functions, or with delayed entry, weights, `strata()`, `cluster()` or
  `tt()`, stops with a clear message.

## Documentation and validation

* New help pages `?Interval`, `?turnbull`, `?frailty`, `?time_effects`, `?adjusted_risk` and
  `?model_scan`; sections on interval-censored data in `?parsurv` and `?directcif`; `?cox` lists
  `frailty()`.
* Vignettes: new sections on time-varying effects, adjusted risks, interval-censored data and
  shared frailty (survival analysis), on `time_effects()`, adjusted cumulative incidence and
  interval-censored competing risks (competing risks), and updated tables of measures, validation
  and speed.
* The validation suite (`validation/` in the repository) now covers the new methods: accuracy
  against `survival::coxph(tt = )`, `riskRegression::ate`, `icenReg`, `interval`,
  `survival::survreg(interval2)`, `flexsurv`, `frailtyEM` and `parfm`; speed; and Monte Carlo
  scenarios for interval-censored Weibull regression and shared gamma frailty. Short check
  scripts `validation/check_time_effects.R`, `check_interval.R` and `check_frailty.R`.
* The fertility examples no longer use `age`, which is recorded at the end of follow-up and is
  therefore not a baseline covariate.

# SurvX 1.0.0

First release on CRAN.

* New data set `fertility`: time from marriage to a live birth or a stillbirth/abortion of 652
  women in Shiraz (Shayan, Ayatollahi and Zare 2011), a real competing risks example with many
  tied times.
* Six vignettes: getting started, survival analysis with one event, competing risks, parametric
  and cure models, interpreting and reporting results, and validation.
* New title and description; lifecycle stable.

* `cox()`: time-varying coefficients with `tt()` terms and the new argument `tt`, a function of
  time: `cox(Event(time, status) ~ x + tt(x), data, tt = function(t) log(t))` fits
  log HR(t) = b + g x log(t). Rows are split at the event times, which reproduces
  `survival::coxph(tt = )` (checked with Efron and Breslow ties and with strata); delayed entry and
  `cluster()` work as in other Cox fits.
  The print shows the hazard ratio at the quartiles of the event times, and the new function
  `tt_effect(fit, times)` gives it, with confidence intervals, at any times.
* `fg()`: `variance = "bootstrap"` (with `B` and `seed`) estimates the variance from bootstrap
  samples of the subjects, so that the estimation of the censoring and truncation weights is taken
  into account. With delayed entry the robust variance treats the Geskus weights as known and was
  5-11% too small in the Monte Carlo study; the print of such fits now says so and points to the
  bootstrap.
* Bootstrap intervals of `cif_plot()` and `cif_data()` no longer change the random number stream of
  the session.
* `Event()` gains `timefix = TRUE`: times (and entry times) that differ only by floating-point
  rounding are treated as tied, with the rule of `survival::aeqSurv()` used by `coxph()`,
  `survfit()` and `survdiff()`. Before, such near-ties (a few per 10000 computed times) were kept
  apart, and results could differ slightly from survival on large data sets.
  `Event(..., timefix = FALSE)` keeps the times as they are, as cmprsk and prodlim do.
  `logrank_test()` and `gray_test()` gain the same argument.
* All help pages rewritten: when to use the function and what to use instead,
  every argument with its default, the model, how to interpret and report the estimates, the
  assumptions and how to check them, messages and their remedies, the agreement with reference
  packages, and references grouped by what they are used for.
* New help page `?SurvX` with a table from research questions to functions.
* Compatibility with S7 1.0.0: the property holding the call of a fit is now `model_call` (S7 1.0.0
  does not keep a property named `call`, which broke `update()` and `getCall()`), and the package
  calls `S7_on_build()` when it is built, as S7 1.0.0 requires for methods of generics from other
  packages. `getCall(fit)` and `update(fit, ...)` work as before.

# SurvX 0.9.1

* The S3 class of `Event()` objects is now `"SurvXEvent"`: the former name `"Event"` was also used
  by timereg (loaded by riskRegression), whose `[` method then replaced SurvX's. `Event()` itself is
  unchanged.
* Cox residuals with Efron ties now follow the Efron form exactly, as in survival: a subject who
  fails at a tied time contributes the Efron-weighted hazard increments to its martingale residual,
  and Schoenfeld residuals use the risk-set mean averaged over the Efron terms (computed in the C++
  core). The martingale residuals of an Efron fit sum to zero. Before, the Breslow form was used for
  both, which differed from survival when event times were tied.
* `parsurv()`, `csc(baseline = )` and `cure()` report a fit with non-finite estimates or standard
  errors as not converged (with the warning in the printed output) instead of returning them silently.
* `directcif(estimation = "weighted")` is much faster, with identical results. The baseline
  quantities at the censoring support points are computed once per iteration instead of once per
  (subject, support point) pair, and the competing-event terms use closed-form derivatives. With the
  PH link (the default) the sums over the censoring distribution become suffix sums, so a fit is
  linear in the sample size: n = 10000 takes about 0.1 s instead of 31 s. Other links remain
  proportional to (competing events x support points) but are about 10 times faster. The new code
  also evaluates log{1 - F} directly, so it no longer rejects parameter values where F rounds to 1.
  Checked against the previous implementation (log-likelihood, score and information to 1e-9) for
  every baseline and link.
* `ph_test(transform = "rank")` now ranks each event time among all observed times, censored
  included, as `survival::cox.zph` does since version 3.0; before, the rank was taken among the event
  times only, and the test statistic differed slightly from `cox.zph`. The other transforms were
  already identical to `cox.zph`.
* `fg()`: the robust (Fine-Gray) variance now treats a censoring time tied with event times exactly
  as `cmprsk::crr` does: the censoring-correction term q(u) includes events of interest at u and
  excludes competing events at u. Estimates were already identical; standard errors could differ by
  up to a few per cent on data with many ties and are now identical to `crr` (checked against a
  line-by-line port of crr's variance routine).
* New validation suite in `validation/` (not part of the package build): accuracy against survival,
  cmprsk, riskRegression, flexsurv, flexsurvcure, survRM2 and timeROC; speed; a Monte Carlo study of
  bias, standard errors and coverage for every model; stability and invariance checks. One run
  writes all results to one text file and all figures to one PDF.

# SurvX 0.9.0

* New `cure()`: parametric cure models, **mixture** (Boag 1949; Berkson & Gage 1952; Farewell 1982)
  and **non-mixture / promotion time** (Yakovlev & Tsodikov 1996; Chen, Ibrahim & Sinha 1999), with
  a logistic model for the cured fraction (`cure = ~ covariates`) and any latency distribution of
  `parsurv()` (PH or AFT, splines, piecewise exponential, generalized gamma), with delayed entry.
  Exact derivatives in the C++ core; validated against independent code (log-likelihood and maximum
  to 1e-12) and, in the tests, against `flexsurvcure`. `predict()` gives the population survival,
  hazard, survival of the uncured and the cured fraction with delta-method intervals; `plot()`
  compares the fitted marginal survival with Kaplan-Meier and shows the cured fraction; the print
  reports odds ratios of being cured and latency effects with Reading lines.
* `fg()` accepts delayed entry (`Event(time, status, entry = )`) with the censoring and truncation
  weights of Geskus (2011), fitted as a weighted counting-process Cox model in the extended C++ core;
  agrees with `survival::finegray` + `coxph(cluster = id)`.
* `cif_plot()` and `cif_data()` now work for `csc()` with `strata()` or delayed entry (one curve per
  stratum at typical covariate values; bootstrap that resamples clusters when `cluster()` is given)
  and for `fg()` with delayed entry. The models are refitted in R for the bootstrap.

# SurvX 0.8.0

* `cox()` and `csc()` (Cox baseline) gain the three extensions that make them cover the uses of
  `survival::coxph`:
  * **strata**: `strata(center)` in the formula gives every stratum its own baseline hazard with
    common coefficients (Kalbfleisch & Prentice 2002, Sec. 4.4); `predict()` uses the stratum of
    each row of `newdata`.
  * **(start, stop] data**: `Event(stop, status, entry = start)` gives delayed entry and
    time-varying covariates in the counting-process form of Andersen & Gill (1982).
  * **robust variance**: `cluster(id)` in the formula, or `robust = TRUE`, gives the sandwich
    variance of Lin & Wei (1989) from score residuals in the Efron or Breslow form; the model-based
    variance stays in `x@vcov_model`.
  All three are implemented in the C++ core (one scan of the risk sets per evaluation) and checked
  against an independent brute-force implementation to 1e-12, and in the tests against
  `survival::coxph` (coefficients, model-based and robust variances, log-likelihood, predicted
  survival, martingale residuals) and `survival::cox.zph`.
* `residuals()` and `ph_test()` handle strata and (start, stop] data.
* The printed output states the strata, the interval form of the data and the kind of standard
  errors, and the Methods line cites the corresponding sources.
* `strata()` and `cluster()` are read by SurvX itself, so they need neither survival nor any other
  package, and they are refused with a clear message by the functions that do not support them
  (`fg()`, `directcif()`, parametric `csc()`).
* `cif_plot()` and `cif_data()` are not yet available for `csc()` fits with strata or delayed
  entry; `predict()` gives their cumulative incidence.

# SurvX 0.7.0

* `parsurv()` and `csc(baseline = )` gain `"gengamma"`, the generalized gamma distribution
  (Prentice 1974; the parametrisation of `flexsurv`), AFT, with delayed entry. It nests the
  log-normal (Q = 0), Weibull (Q = 1) and gamma (Q = sigma) distributions, so the estimate of Q
  helps to choose among them. The log density is computed without cancellation near Q = 0 and the
  survivor function switches from the incomplete gamma function to numerical integration of the
  density for |Q| < 0.01, so the likelihood is smooth across the log-normal limit. Validated against
  independent high-precision code and, in the tests, `flexsurv::flexsurvreg(dist = "gengamma")`.
* New `residuals()` methods for `cox()`, `csc()` and `parsurv()` fits: martingale, deviance and
  Cox-Snell residuals, and for Cox models Schoenfeld and scaled Schoenfeld residuals. Agree with
  `survival` (Breslow).
* New `ph_test()`: the score test of proportional hazards of Grambsch and Therneau (1994) for every
  covariate and globally, in the exact form of `survival::cox.zph` (>= 3.0), with time scales
  `"km"`, `"rank"`, `"identity"` and `"log"`; works for every cause of a cause-specific Cox model.
  The print gives a Reading line; `plot()` shows the scaled Schoenfeld residuals with a natural
  spline smooth and its confidence band (estimate of beta(t)).
* New `dist_plot()`: linearising plots of the Nelson-Aalen estimate for the exponential, Weibull,
  log-normal and log-logistic distributions, by group; with a fitted model, the Cox-Snell residual
  plot.
* Printed output of `cox()`, `fg()`, `csc()` and `directcif()` now has the same layout as the
  parametric models: Reading lines that state each effect with p < 0.05 in words on its own scale
  (rate for HR and csHR, cumulative incidence for sHR, odds of having had the event for OR), and a
  Fit line (log-likelihood, AIC, and for Cox models the likelihood ratio test against no
  covariates).
* Every result names its scientific sources. The printed output of every model, estimate, test and
  comparison ends with a "Methods:" line (for example "Cox (1972); Efron (1977); Breslow (1974)"),
  which depends on the options actually used (ties, estimation method, baseline, delayed entry).
  The new `sx_references(fit)` prints the full references behind a fitted object, and
  `sx_references("brier_score")` those of a function. Every help page has a References section;
  `predict()`, `estimates()` and the generic methods (AIC, BIC, confint) gained theirs.
* `ph_test()` returns an S7 object (`PHTest`, with `print`, `plot` and `as.data.frame`; all tests
  in `x@table`), like every other result of the package.
* `plot()` now works for parametric fits: survival curves for `parsurv()`, cumulative incidence
  curves for `csc(baseline = )`; `logLik()`, `AIC()` and `BIC()` work for cause-specific Cox models
  (the sum of the log partial likelihoods of the causes).
* New test file of all generic functions (print, summary, plot, coef, vcov, confint, logLik, AIC,
  BIC, nobs, formula, getCall, terms, model.matrix, update, predict, residuals, as.data.frame) on
  every class, run on all platforms including R-devel.
* `splines` (base R) is now imported, for the smooth of `plot(ph_test())`.

# SurvX 0.6.1

* `parsurv()` and `csc(baseline = )` gain two baselines: `"rp"`, the flexible parametric model of
  Royston and Parmar (restricted cubic spline on the log cumulative hazard, the basis of
  `flexsurv::flexsurvspline`; `df` or `knots`), and `"pwexp"`, the piecewise exponential model
  (`cuts`). Analytic gradient and Hessian; validated against independent code of the likelihoods and,
  in the tests, against `flexsurv::flexsurvspline` and a Poisson GLM on split data. The cumulative
  incidence of a piecewise exponential cause is integrated exactly up to its cut points.
* `compare_models()` gains `hazard`: parametric cause-specific hazard models, e.g.
  `hazard = list("weibull", "lognormal", c("rp", 4))`, are cross-validated with the same folds as the
  direct, Fine-Gray and Cox models, and their AIC is shown beside that of the joint and constrained
  direct models (all are full likelihoods of the same data).
* `cv_compare()` reports log-likelihood and AIC for parametric fits; internally `compare_models()` and
  `cv_compare()` share one cross-validation routine.
* (The generalized gamma distribution followed in 0.7.0.)

# SurvX 0.6.0

* New `parsurv()`: parametric survival regression with exponential, Weibull, Gompertz, log-normal and
  log-logistic distributions, as proportional hazards (PH) or accelerated failure time (AFT) models,
  with delayed entry. Exact (analytic) gradient and Hessian in C++; validated against
  `survival::survreg` and `flexsurv::flexsurvreg`. `predict()` gives survival, hazard, cumulative
  hazard, quantiles (e.g. the median) and restricted mean survival time, with delta-method standard
  errors and limits, for given covariate values or averaged over the sample.
* `csc()` gains `baseline`: with a parametric distribution every cause gets a parametric hazard model
  and the cumulative incidence is computed from the hazards by numerical integration, with
  delta-method bands; coherent by construction. Works in `cif_plot()`, `cif_data()`, `cv_compare()`,
  `AIC()` and the new `hazard_plot()`.
* New `hazard_np()`: Nelson-Aalen cumulative hazard and kernel-smoothed hazard (Epanechnikov kernel
  with the boundary kernels of Mueller and Wang; bandwidth by the least-squares cross-validation of
  Patil, taking the largest local minimum as recommended by Hall and Marron), by group and cause, with
  delayed entry. The smoothed hazard is shown only where at least max(10, 10% of the initial risk
  set) subjects are at risk.
* New `hazard_plot()`: fitted hazard, cumulative hazard or survival of parametric models with
  confidence bands and the non-parametric estimate as reference, to check the chosen distribution.
* `Event()` gains `entry` for delayed entry (left truncation). Functions that do not support it yet
  stop with a clear message.
* Output of the new models: effects on the ratio scale, a plain-language reading of every effect with
  p < 0.05 ("x (per one-unit increase): 12% lower rate of the event (HR 0.88, 95% CI ...)"), the
  baseline distribution on its natural scale in `summary()`, AIC, BIC and the likelihood ratio test.

# SurvX 0.5.4

* Test suite with testthat (`tests/testthat/`), replacing `tests/smoke.R`. It checks the estimators
  against survival, cmprsk and riskRegression when these are installed (Kaplan-Meier, Greenwood
  standard errors, log-rank and Gray tests, Aalen-Johansen, Cox, Fine-Gray, cause-specific Cox,
  Brier score and time-dependent AUC), the direct models against independent R code of the published
  likelihoods (joint, weighted and constrained), and every print, summary, predict and plot method.
  Run it with `devtools::test()` or as part of `R CMD check`.
* New vignette "Getting started with SurvX" (`vignette("SurvX")`): a complete competing risks
  analysis from the Aalen-Johansen estimator to cross-validated model comparison and the coherence check.
* New benchmark script `inst/benchmark/benchmark_SurvX.R`: accuracy (coefficients, standard
  errors, predicted CIFs, Gray test, Brier score, AUC), stability (convergence, errors and warnings
  with tied times, rare events and a nearly separating covariate) and speed (n = 500 to 50 000) of
  SurvX against cmprsk, survival and riskRegression on the same simulated data.
* The GitHub checks install pandoc, which the vignette needs.
* `R CMD check` notes removed: global variables of the ggplot2 layers declared, `utils` imported,
  an example line of `directcif()` shortened.

# SurvX 0.5.3

* `cif_plot()`: when several models are compared at one covariate profile (for example
  `type = "marginal", reference = TRUE`), the curves are labelled by the model names only; the long
  labels "model: average over the sample" no longer overflow the legend.
* Without grouping variables the single group is called "All patients" (was "all"); the top band of
  the stacked layout is labelled "no event yet" (was "event-free").

# SurvX 0.5.2

* New simulation study in `inst/simulations/`, replacing `direct_vs_finegray.R`. It has a shared engine (`sim_core.R`) and
  two studies: published models only (`sim_published_models.R`) and published models plus the SurvX
  extensions (`sim_survx_models.R`). The study covers eight scenarios, from proportional subdistribution
  hazards to a strongly waning effect, plus two scenarios unfavourable to the direct models. It uses
  paired Monte Carlo comparisons with Fine-Gray, parallel computing with per-replication seeds, and
  resumable runs. BMC-format figures and LaTeX tables come from `sim_figures.R`.

# SurvX 0.5.1

* Long titles ran off the page because `mtext()` in the outer margin ignores `par("cex")`: the title
  was measured at 0.85 x 1.1 but drawn at 1.1. Titles and captions are now drawn at exactly the size
  at which they are measured (checked against the Helvetica font metrics: the test title is 8.6
  inches wide at the default size, the page 10 inches).
* The validation script writes the package version on the first page of the PDF and in the first
  line of the output file.

# SurvX 0.5.0 (first public development version)

* README rewritten for GitHub: every function with its method and references, interpretation of
  the effects, plots, validation results, methods original to the package.
* `inst/CITATION`, GitHub Actions workflow for `R CMD check` on Linux, macOS and Windows.
* `cif_plot()`: values written by `mark.times` are stacked when several curves are drawn, so they no
  longer overlap; the note "sum of the risks > 1" of the stacked layout is drawn inside the panel
  so it does not collide with the panel title.

# SurvX 0.4.4

* Titles are measured in bold, the font they are drawn in, so a long title is wrapped instead of
  running off the page (0.4.3 measured them in the plain font, about 8% narrower); the `km()` plot
  wraps its title too.
* Maintainer: Habib Ezzatabadipour.

# SurvX 0.4.3

* All plots get a `cex` argument (size of all text); the default is now 0.85 instead of 1.
* Text widths are measured with the font metrics of the device instead of being estimated from
  the number of characters; titles, captions and the axis title of the comparison plot are wrapped
  so that they fit. ggplot2 titles are wrapped to the width of the open device.
* The check script stops with a clear message when an older SurvX is loaded.

# SurvX 0.4.2

* Plot titles and captions are wrapped to the width of the device (base graphics) or to about 100
  characters (ggplot2), so long covariate descriptions no longer run off the page; ggplot2 titles
  start at the left edge of the whole figure.
* `plot()` of a model comparison: labels written once at the left with their width measured, so
  neither the model names nor the axis title are cut; same model order in every panel, best model
  highlighted, dashed line at IPA = 0, new arguments `times`, `highlight`, `reference`, `col`, `main`.
* ggplot2 stacked layout: same stacking as the base engine (cause 1 at the bottom), stairs for step
  curves, running sum capped at 1 with a note when the risks add up to more than 1.
* `cif_data(type = "unadjusted")` for a cause-specific Cox fit now uses Breslow increments and the
  product-limit survival, so it is exactly the Aalen-Johansen estimate also with tied times
  (0.4.1 differed by up to 0.0014 on mgus2 because of Efron ties).

# SurvX 0.4.1

* `cif_plot()` and `cif_data()`: cumulative incidence curves of every cause for one model or a
  comparison of models (direct, Fine-Gray, cause-specific Cox, Aalen-Johansen), at given covariate
  values, at typical values (mean of numeric, most frequent level of categorical covariates), as the
  average over the sample, or from the models fitted without covariates. Pointwise confidence bands
  by the delta method (direct models, Aalen-Johansen) and by the bootstrap of the subjects (Fine-Gray,
  cause-specific Cox), computed in the C++ core. Layouts: one panel per cause, all causes in one
  panel, or stacked. Number at risk, legend, caption, bands, marked times, percent axis and grid can
  each be switched on or off. Base graphics, or ggplot2 with `engine = "ggplot2"` and `autoplot()`.
* `plot()` of `directcif()`, `fg()`, `csc()` and `aj()` fits is now `cif_plot()`; `plot()` of `km()`
  fits gains the number at risk, censoring marks and the same options.
* The printed long-run risk of direct models is given at typical covariate values (previously at the
  mean of the design matrix).

# SurvX 0.4.0

* New output layer for all models: readable reports from `print()` and `summary()` (effects with
  confidence intervals, factor levels against their reference, long-run risk of direct models,
  warnings in plain words), `estimates()`, `confint()`, confidence bands in `plot()`, ranked tables
  for `compare_models()`, and the S3 generics `formula`, `terms`, `model.matrix`, `getCall`,
  `update` and `as.data.frame`.
* `km()` and `aj()` gain confidence intervals (log scale for survival, complementary log-log for the
  cumulative incidence), the median survival time with its limits and numbers at risk in `summary()`.
* `summary()` of `km()` and `aj()` fits now returns an object; the tables are in `@tables`.

# SurvX 0.3.1

* `compare_models()` fits, cross-validates and scores weighted and constrained direct models in the C++ core.

# SurvX 0.3.0

* `directcif(estimation = "weighted")`: weighted likelihood for one cause with censoring weights.
* `directcif(estimation = "constrained")`: the model of Shi, Cheng and Jeong (2013).
* `cv_compare()` for any list of models.

# SurvX 0.2.1

* First version with the R interface: Kaplan-Meier, Aalen-Johansen, log-rank and Gray tests, Cox,
  cause-specific Cox, Fine-Gray, direct parametric models with eight baselines, prediction measures.
