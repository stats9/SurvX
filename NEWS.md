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
