# model_scan(): fit a set of candidate models to one data set and compare them on criteria that
# are valid for the comparison at hand. Orchestration only: the fits and the measures are the
# functions of the package (cox, parsurv, cure, csc, fg, directcif, brier_score, auc_t,
# calibration, ph_test).

# ----------------------------------------------------------------------------------------------
# candidate models
# ----------------------------------------------------------------------------------------------
# one candidate: name, family (for the grouping of the likelihood comparison), likelihood type
# ("full", "partial", "pseudo"), fit = function(formula, data), risk = function(fit, newdata, times)
.scan_cand <- function(name, family, lik, fit, risk, kind) list(name = name, family = family, lik = lik, fit = fit, risk = risk, kind = kind)

.risk_surv <- function(fit, nd, times) 1 - as.matrix(predict(fit, nd, times = times, type = "survival"))

.scan_candidates <- function(type, cause, depth, ev, plateau, K) {
  cand <- list()
  add <- function(x) cand[[length(cand) + 1L]] <<- x
  rich <- ev >= 50L                                   # enough events for the flexible models
  if (type == "single") {
    add(.scan_cand("Cox", "semi-parametric", "partial", function(f, d) cox(f, d), .risk_surv, "cox"))
    par <- list(c("exponential", "PH"), c("weibull", "PH"), c("gompertz", "PH"), c("lognormal", "AFT"), c("loglogistic", "AFT"))
    if (depth != "quick" && rich) par <- c(par, list(c("gengamma", "AFT"), c("rp", "PH")))
    if (depth == "full" && rich) par <- c(par, list(c("pwexp", "PH")))
    for (p in par) {
      local({
        b <- p[1L]; s <- p[2L]
        add(.scan_cand(paste0("parametric ", .par_label[[b]], " (", s, ")"), "parametric (full likelihood)", "full",
                       function(f, d) parsurv(f, d, baseline = b, scale = s), .risk_surv, "parsurv"))
      })
    }
    if (plateau) {
      cb <- if (depth == "quick") "weibull" else c("weibull", "lognormal")
      for (b in cb) for (ty in if (depth == "full") c("mixture", "nonmixture") else "mixture") {
        local({
          bb <- b; tt <- ty
          add(.scan_cand(paste0("cure ", tt, " ", .par_label[[bb]]), "parametric (full likelihood)", "full",
                         function(f, d) cure(f, d, baseline = bb, type = tt), .risk_surv, "cure"))
        })
      }
    }
    return(cand)
  }
  rc <- function(fit, nd, times) as.matrix(predict(fit, nd, times = times, cause = cause))
  add(.scan_cand("cause-specific Cox", "semi-parametric", "partial", function(f, d) csc(f, d), rc, "csc"))
  add(.scan_cand("Fine-Gray", "semi-parametric", "pseudo", function(f, d) fg(f, d, cause = cause),
                 function(fit, nd, times) as.matrix(predict(fit, nd, times = times)), "fg"))
  hz <- if (depth == "quick") "weibull" else c("weibull", "lognormal", if (rich) "rp")
  for (b in hz) local({
    bb <- b
    add(.scan_cand(paste0("cause-specific ", .par_label[[bb]]), "parametric (full likelihood)", "full",
                   function(f, d) csc(f, d, baseline = bb), rc, "csc_par"))
  })
  dir <- list(c("gompertz", "PH", "joint"), c("gompertz", "PO", "joint"), c("gompertz", "PH", "weighted"), c("gompertz", "PO", "weighted"))
  if (depth != "quick") {
    dir <- c(dir, list(c("gompertz", "PH", "constrained")))
    if (rich) dir <- c(dir, list(c("gompertz", "GOR", "weighted"), c("mix_weibull", "PH", "weighted")))
  }
  if (depth == "full" && rich) dir <- c(dir, list(c("spline", "PH", "weighted"), c("spline", "PO", "weighted"), c("logistic3", "PH", "constrained")))
  for (s in dir) local({
    b <- s[1L]; l <- s[2L]; e <- s[3L]
    nm <- paste0("direct ", b, " ", l, ", ", e)
    fam <- if (e == "weighted") "direct, weighted" else "parametric (full likelihood)"
    add(.scan_cand(nm, fam, if (e == "weighted") "pseudo" else "full",
                   function(f, d) directcif(f, d, baseline = b, link = l, estimation = e, cause = cause,
                                            cond = if (e == "constrained") "weibull" else "logistic"),
                   rc, paste0("direct_", e)))
  })
  cand
}

# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------
.scan_converged <- function(fit) {
  if (S7::S7_inherits(fit, CSCFit) || S7::S7_inherits(fit, ParCSCFit))
    return(all(vapply(fit@models, function(m) isTRUE(if (is.list(m)) m$converged %||% TRUE else m@converged), logical(1))))
  v <- tryCatch(fit@converged, error = function(e) TRUE)
  isTRUE(v)
}

.scan_loglik <- function(fit) {
  if (S7::S7_inherits(fit, DirectCIFFit)) return(c(fit@loglik, length(fit@coefficients)))
  ll <- tryCatch(logLik(fit), error = function(e) NULL)
  if (is.null(ll)) return(c(NA_real_, NA_real_))
  c(as.numeric(ll), attr(ll, "df"))
}

# trapezoidal integral of y(t) on the grid (with y(0) = 0), divided by the last grid point
.scan_integrate <- function(grid, y) {
  if (any(!is.finite(y))) return(NA_real_)
  x <- c(0, grid); v <- c(0, y)
  sum(diff(x) * (v[-1L] + v[-length(v)]) / 2) / max(grid)
}

.scan_brier_curve <- function(time, status, risk, grid, cause) {
  out <- matrix(NA_real_, length(grid), 2L)
  for (j in seq_along(grid)) {
    b <- tryCatch(brier_score(time, status, risk[, j], grid[j], cause), error = function(e) NULL)
    if (!is.null(b)) out[j, ] <- b[1:2]
  }
  out
}

# the non-parametric estimate of the event probability of the cause at the times (no covariates)
.scan_np <- function(time, status, cause, at) {
  d <- data.frame(time = time, status = status)
  if (max(status) <= 1L) {
    cv <- km(Event(time, status) ~ 1, d)@curves[[1L]]
    1 - .step_at(cv$time, cv$surv, at, 1)
  } else {
    cv <- aj(Event(time, status) ~ 1, d)@curves[[1L]]
    .step_at(cv$time, cv[[paste0("cif.", cause)]], at, 0)
  }
}

# a plateau in the Kaplan-Meier curve: a screening rule, not a test. The curve is flat over the
# last fifth of follow-up or more, with at least 10% of the subjects event-free at its end.
.scan_plateau <- function(time, status) {
  ev <- time[status > 0L]
  if (!length(ev)) return(FALSE)
  d <- data.frame(time = time, status = as.integer(status > 0L))
  cv <- km(Event(time, status) ~ 1, d)@curves[[1L]]
  s_end <- utils::tail(cv$surv, 1L)
  (max(time) - max(ev)) >= 0.2 * max(time) && s_end >= 0.1 && sum(time > max(ev)) >= 10L
}

.scan_folds <- function(status, folds, seed) {
  n <- length(status)
  .with_seed(seed, {
    fold <- integer(n)
    for (s in unique(status)) { i <- which(status == s); fold[i] <- sample(rep_len(seq_len(folds), length(i))) }
    fold
  })
}

.scan_checks <- function(fit, cand, n) {
  msg <- character(0)
  if (!.scan_converged(fit)) msg <- c(msg, "did not converge")
  if (cand$kind %in% c("cox", "csc")) {
    pt <- tryCatch(ph_test(fit)@table, error = function(e) NULL)
    if (!is.null(pt)) {
      g <- pt[pt$term == "GLOBAL", , drop = FALSE]
      if (nrow(g) && any(g$p < 0.05)) msg <- c(msg, sprintf("PH doubtful (global p = %s)", .pval(min(g$p))))
    }
  }
  if (cand$kind == "direct_joint") {
    d <- fit@diagnostics
    if (isTRUE(d$n_exceed_range > 0)) msg <- c(msg, sprintf("sum of CIFs > 1 for %s%% within follow-up", .num(100 * d$n_exceed_range / n, 0L)))
  }
  if (!length(msg)) "ok" else paste(msg, collapse = "; ")
}

# ----------------------------------------------------------------------------------------------
# model_scan
# ----------------------------------------------------------------------------------------------
#' Fit and compare all suitable models for a data set
#'
#' @description
#' Fits every model of the package that suits the data (one event type or competing risks) and
#' compares them on criteria that are valid for the comparison at hand: cross-validated prediction
#' for all models, the likelihood only within the models that share one, and checks of the
#' assumptions of every model. The result is a ranking with the models that predict equally well
#' grouped together, a reading in plain words, figures and, with \code{report()}, a complete HTML
#' report with the methods and references.
#'
#' \strong{Use it when} you want an overview of how the models of the package fit your data before
#' choosing one, or want to show that the conclusion does not depend on the model.
#'
#' \strong{It does not choose for you.} Models whose prediction error differs by less than its
#' standard error are equally good for prediction; choose among them on the question (rates or
#' risks), the assumptions and simplicity. Confidence intervals of a model chosen from the table
#' do not account for the choice.
#' @param formula \code{Event(time, status) ~ covariates}, with at least one covariate. Status
#'   0/1 for one event type; 0, 1, ..., K for competing risks. Delayed entry is not supported.
#' @param data a data frame. Rows with missing values in the variables of \code{formula} are
#'   removed.
#' @param cause competing risks: the cause whose cumulative incidence is predicted and compared.
#' @param goal \code{"prediction"} (default): the ranking is the answer. \code{"effect"}: in
#'   addition, \code{summary()} and the report show the effect of every covariate in every model,
#'   side by side, on the scale of each model (hazard ratio, subdistribution hazard ratio, odds
#'   ratio, time ratio).
#' @param times times for the AUC and the calibration error. Default: the quartiles of the event
#'   times of the cause.
#' @param horizon upper limit of the integrated Brier score. Default: the smaller of the 90th
#'   percentiles of the event times of the cause and of all observed times.
#' @param folds number of cross-validation folds (at least 2).
#' @param seed seed of the fold assignment; the state of the random number generator of the session
#'   is restored afterwards.
#' @param depth how many models: \code{"quick"} (the basic models), \code{"standard"} (default), or
#'   \code{"full"} (also the more flexible ones). Flexible models (generalized gamma, splines, GOR
#'   link) are tried only with at least 50 events of the cause.
#' @param models optional named list of additional models, each a function
#'   \code{function(formula, data)} returning a fit with a \code{predict(fit, newdata, times, ...)}
#'   method (the probability of the event or cumulative incidence of the cause, or for one event
#'   type \code{type = "survival"}).
#' @param verbose print the name of every model while it is fitted.
#' @param x a result of \code{model_scan()}.
#' @param file the HTML file to write.
#' @param ... \code{report()}: \code{width} and \code{height} of the figures in inches.
#' @details
#' \subsection{The models}{
#' \describe{
#'   \item{One event type}{\code{\link{cox}}; \code{\link{parsurv}} with the exponential, Weibull,
#'     Gompertz (PH), log-normal and log-logistic (AFT) distributions, with \code{depth = "standard"}
#'     also the generalized gamma and the Royston-Parmar spline, with \code{"full"} the piecewise
#'     exponential; \code{\link{cure}} models (mixture; with \code{"full"} also non-mixture) when the
#'     Kaplan-Meier curve levels off. The screening rule for a plateau: the curve is flat over at
#'     least the last fifth of follow-up, with at least 10\% of the subjects event-free and at least
#'     10 subjects censored after the last event. It is not a test; look at \code{plot(km(...))}.}
#'   \item{Competing risks}{cause-specific Cox and parametric (Weibull, log-normal, spline) models
#'     (\code{\link{csc}}); \code{\link{fg}}; direct models (\code{\link{directcif}}): joint Gompertz
#'     with PH and PO links, weighted Gompertz with PH, PO and GOR links and a Weibull mixture
#'     baseline, the constrained model, and with \code{"full"} spline and logistic baselines.}
#' }
#' }
#' \subsection{The criteria}{
#' Three kinds of evidence are kept apart:
#' \enumerate{
#'   \item \strong{Prediction, for every model.} Every model is fitted in K-fold cross-validation and
#'     predicts the probability of the event (or cumulative incidence of the cause) of the left-out
#'     subjects. On the pooled out-of-fold predictions: the integrated Brier score (IBS) over
#'     \code{(0, horizon]} with inverse probability of censoring weights, the index of prediction
#'     accuracy (IPA, 1 minus the ratio to the model without covariates), the time-dependent AUC and
#'     the calibration error at \code{times}. These measures compare any models, whatever their
#'     likelihood.
#'   \item \strong{Likelihood, within comparable models only.} AIC is reported for the models with a
#'     full likelihood of the same data: the parametric and cure models for one event type; the
#'     parametric cause-specific, joint and constrained direct models for competing risks (their
#'     likelihoods are all the joint density of the time and the cause). The Cox and cause-specific
#'     Cox models have a partial likelihood, Fine-Gray and weighted direct models a weighted
#'     pseudo-likelihood: no AIC.
#'   \item \strong{Checks.} Convergence; the test of proportional hazards (\code{\link{ph_test}}) for the
#'     Cox models; for joint direct models, subjects whose predicted incidences add up to more than
#'     one; and for every model the largest distance between its average predicted risk and the
#'     non-parametric estimate (Kaplan-Meier or Aalen-Johansen) on the evaluation grid.
#' }
#' }
#' \subsection{Ranking}{
#' The models are ranked on the cross-validated IBS. The standard error of the difference from the
#' best model is the standard deviation of the fold-wise differences in IBS divided by the square
#' root of the number of folds (a paired comparison: the same folds for every model). Every model
#' within one standard error of the best gets rank 1: the one-standard-error rule of Hastie,
#' Tibshirani and Friedman (2009, Section 7.10). The other models are numbered after them in the
#' order of their IBS. The reading also says when no model predicts clearly better than the model
#' without covariates (largest IPA below 0.01), and when AIC clearly disfavours (dAIC > 10) a model
#' of rank 1.
#' }
#' @section Interpreting the result:
#' \describe{
#'   \item{\code{IBS}, \code{diff (SE)}}{prediction error, lower is better, and its difference from the
#'     best model with the standard error of that difference. A difference smaller than its SE is
#'     not evidence that one model predicts better.}
#'   \item{\code{IPA}}{share of the prediction error of the model without covariates that the model
#'     removes; 0 = no better than no covariates, negative = worse.}
#'   \item{\code{AUC}}{discrimination at \code{times} (averaged): the probability that a subject who
#'     has the event by \eqn{t} has a higher predicted risk than one who does not.}
#'   \item{\code{calib}}{calibration error at \code{times} (averaged): mean absolute difference
#'     between predicted and observed risk in ten risk groups. Small is good.}
#'   \item{\code{dAIC}}{difference from the smallest AIC among the comparable models. Differences
#'     below 2 are small, above 10 large (Burnham and Anderson 2002, Section 2.6).}
#'   \item{\code{checks}}{"ok", or the problems found. A model with a problem can still predict well;
#'     it should not be used to interpret effects without addressing the problem.}
#'   \item{Rates or risks}{with competing risks, the cause-specific models describe the rate of the
#'     event among those still event-free; Fine-Gray and the direct models describe the cumulative
#'     incidence. They answer different questions, and a covariate can have different effects on
#'     the two (Putter, Fiocco and Geskus 2007; Latouche et al. 2013). \code{goal = "effect"} shows
#'     both side by side.}
#' }
#' A sentence for a report: \emph{"We compared twelve candidate models in 5-fold cross-validation.
#' Four models, the cause-specific Cox model, the Fine-Gray model and two direct parametric models,
#' had integrated Brier scores within one standard error of the best (0.112); we report the
#' Fine-Gray model, whose subdistribution hazard ratios refer to the cumulative incidence."}
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{Independent censoring}{as for every measure that uses censoring weights; the Kaplan-Meier
#'     estimate of the censoring distribution has no covariates.}
#'   \item{Enough events}{with fewer than about 10 events per fold and cause, the fold-wise measures
#'     are unstable and many flexible models fail; use \code{folds = 3} or \code{depth = "quick"}.}
#'   \item{The candidate set}{the ranking only compares the models that were fitted. Models that
#'     failed are listed under "Not fitted".}
#' }
#' @section Messages and common problems:
#' \describe{
#'   \item{\code{Not fitted: ...}}{the model failed on the full data or did not give finite
#'     predictions; with \code{verbose = TRUE} the models are named as they are fitted.}
#'   \item{\code{IBS} is \code{-}}{the model failed in at least one fold (\code{folds_failed} in
#'     \code{as.data.frame(x)}).}
#'   \item{Long computing times}{every model is fitted \code{folds + 1} times; the weighted GOR and
#'     spline models are the slowest. Use \code{depth = "quick"} first.}
#' }
#' @section Validation:
#' \code{model_scan()} calls the functions of the package, which are validated against reference
#' implementations (see their help pages and \code{vignette("validation", package = "SurvX")}),
#' and the measures of \code{\link{brier_score}}, \code{\link{auc_t}} and \code{\link{calibration}}.
#' The fold assignment, ranking and standard errors are tested in the test suite of the package.
#' @return
#' An S7 object of class \code{ModelScan}. \code{x@table} (also \code{as.data.frame(x)}) has one
#' row per model: \code{model}, \code{family}, \code{likelihood} (\code{"full"}, \code{"partial"},
#' \code{"pseudo"}), \code{fitted}, \code{converged}, \code{loglik}, \code{npar}, \code{aic},
#' \code{delta_aic}, \code{ibs_cv}, \code{ipa_cv}, \code{auc_cv}, \code{calib_cv},
#' \code{np_distance}, \code{folds_failed}, \code{checks}, \code{rank}, \code{diff_best},
#' \code{se_diff}, \code{equivalent}. Also \code{x@per_time} (measures at each of \code{times}),
#' \code{x@brier} (cross-validated Brier score on the grid), \code{x@effects} (with
#' \code{goal = "effect"}) and \code{x@fits} (the models fitted to the full data, for further use).
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(x)} \tab the ten best models and the reading \cr
#'   \code{summary(x)} \tab all models; with \code{goal = "effect"} the effects of the covariates \cr
#'   \code{plot(x, which)} \tab \code{"ranking"} (IBS with its SE), \code{"fit"} (average prediction
#'     against the non-parametric estimate), \code{"brier"} (prediction error over time),
#'     \code{"calibration"}; \code{top} models in the last three \cr
#'   \code{report(x, file)} \tab a self-contained HTML report: tables, figures, methods, references \cr
#'   \code{as.data.frame(x)} \tab the table \cr
#'   \code{\link{sx_references}(x)} \tab the references of the methods \cr
#' }
#' @references
#' \emph{Prediction error and its estimation under censoring}
#'
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999) Assessment and comparison of prognostic
#' classification schemes for survival data. Stat Med 18:2529-2545.
#'
#' Gerds TA, Schumacher M (2006) Consistent estimation of the expected Brier score in general
#' survival models with right-censored event times. Biom J 48:1029-1040.
#'
#' Kattan MW, Gerds TA (2018) The index of prediction accuracy: an intuitive measure useful for
#' evaluating risk prediction models. Diagn Progn Res 2:7.
#'
#' \emph{Discrimination and calibration}
#'
#' Blanche P, Dartigues JF, Jacqmin-Gadda H (2013) Estimating and comparing time-dependent areas
#' under receiver operating characteristic curves for censored event times with competing risks.
#' Stat Med 32:5381-5397.
#'
#' Gerds TA, Andersen PK, Kattan MW (2014) Calibration plots for risk prediction models in the
#' presence of competing risks. Stat Med 33:3191-3203.
#'
#' \emph{Cross-validation and model choice}
#'
#' Hastie T, Tibshirani R, Friedman J (2009) The Elements of Statistical Learning, 2nd ed.
#' Springer, Section 7.10.
#'
#' Akaike H (1974) A new look at the statistical model identification. IEEE Trans Automat Contr
#' 19:716-723.
#'
#' Burnham KP, Anderson DR (2002) Model Selection and Multimodel Inference: A Practical
#' Information-Theoretic Approach, 2nd ed. Springer.
#'
#' \emph{Rates and risks with competing risks}
#'
#' Putter H, Fiocco M, Geskus RB (2007) Tutorial in biostatistics: competing risks and multi-state
#' models. Stat Med 26:2389-2430.
#'
#' Latouche A, Allignol A, Beyersmann J, Labopin M, Fine JP (2013) A competing risks analysis should
#' report results on all cause-specific hazards and cumulative incidence functions. J Clin
#' Epidemiol 66:648-653.
#'
#' All references: \code{sx_references("model_scan")}.
#' @seealso
#' \code{\link{compare_models}} and \code{\link{cv_compare}} (a chosen set of models at chosen times),
#' \code{\link{brier_score}}, \code{\link{ph_test}}, \code{\link{km}}, \code{\link{aj}}.
#' @examples
#' ## one event type
#' set.seed(1)
#' n <- 300
#' d <- data.frame(x = rbinom(n, 1, 0.5), z = rnorm(n))
#' t <- rweibull(n, 1.4, exp(1 - 0.5 * d$x + 0.3 * d$z)); cens <- runif(n, 0, 5)
#' d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
#' sc <- model_scan(Event(time, status) ~ x + z, d, depth = "quick", folds = 3)
#' sc
#' plot(sc, which = "ranking")
#'
#' \donttest{
#' ## competing risks, with the effects of the covariates in every model
#' data(fertility)
#' sc2 <- model_scan(Event(time, status) ~ education + age, fertility, cause = 1,
#'                   goal = "effect", depth = "quick", folds = 3)
#' summary(sc2)
#' plot(sc2)
#' report(sc2, file = tempfile(fileext = ".html"))
#' }
#' @export
model_scan <- function(formula, data, cause = 1L, goal = c("prediction", "effect"), times = NULL, horizon = NULL,
                       folds = 5L, seed = 1L, depth = c("standard", "quick", "full"), models = NULL, verbose = FALSE) {
  cl <- match.call()
  goal <- match.arg(goal); depth <- match.arg(depth)
  cause <- as.integer(cause); folds <- as.integer(folds)
  if (missing(data) || is.null(data)) stop("data must be given", call. = FALSE)
  vars <- intersect(all.vars(formula), names(data))
  data <- data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  f <- .sx_frame(formula, data, what = "model_scan()")
  n <- length(f$time); p <- ncol(f$X); K <- max(f$status)
  if (nrow(data) != n) stop("could not align the response with the data", call. = FALSE)
  if (p < 1L) stop("model_scan() needs at least one covariate", call. = FALSE)
  if (K < 1L) stop("no events", call. = FALSE)
  if (K == 1L && cause != 1L) stop("with one event type, cause must be 1", call. = FALSE)
  if (cause < 1L || cause > K || !any(f$status == cause)) stop("no events of cause ", cause, call. = FALSE)
  if (folds < 2L) stop("folds must be at least 2", call. = FALSE)
  type <- if (K == 1L) "single" else "competing"
  ev <- sum(f$status == cause)
  tc <- f$time[f$status == cause]
  if (is.null(times)) times <- signif(as.numeric(stats::quantile(tc, c(0.25, 0.5, 0.75))), 3L)
  times <- sort(unique(as.double(times)))
  if (is.null(horizon)) horizon <- min(as.numeric(stats::quantile(tc, 0.9)), as.numeric(stats::quantile(f$time, 0.9)))
  grid <- seq(horizon / 10, horizon, length.out = 10L)
  at <- sort(unique(c(grid, times)))
  ig <- match(grid, at); it <- match(times, at)
  plateau <- type == "single" && .scan_plateau(f$time, f$status)

  cand <- .scan_candidates(type, cause, depth, ev, plateau, K)
  if (!is.null(models)) {
    if (!is.list(models) || is.null(names(models)) || any(names(models) == ""))
      stop("models must be a named list of functions function(formula, data)", call. = FALSE)
    for (nm in names(models)) local({
      fun <- models[[nm]]
      risk <- if (type == "single") function(fit, nd, times) {
        r <- tryCatch(as.matrix(predict(fit, nd, times = times, cause = 1L)), error = function(e) NULL)
        if (is.null(r)) .risk_surv(fit, nd, times) else r
      } else function(fit, nd, times) as.matrix(predict(fit, nd, times = times, cause = cause))
      cand[[length(cand) + 1L]] <<- .scan_cand(nm, "user", "unknown", fun, risk, "user")
    })
  }
  M <- length(cand)
  fold <- .scan_folds(f$status, folds, seed)
  fits <- vector("list", M); names(fits) <- vapply(cand, `[[`, "", "name")
  rows <- list(); per_time <- list(); bcurve <- list(); foldibs <- matrix(NA_real_, M, folds)
  np_grid <- .scan_np(f$time, f$status, cause, at)
  risk_ok <- function(r, nn) !is.null(r) && is.matrix(r) && nrow(r) == nn && ncol(r) == length(at) && all(is.finite(r))
  for (m in seq_len(M)) {
    cd <- cand[[m]]
    if (verbose) message("model_scan: ", cd$name)
    full <- tryCatch(suppressWarnings(cd$fit(formula, data)), error = function(e) NULL)
    rfull <- if (is.null(full)) NULL else tryCatch(cd$risk(full, data, at), error = function(e) NULL)
    if (!risk_ok(rfull, n)) {
      rows[[m]] <- data.frame(model = cd$name, family = cd$family, likelihood = cd$lik, fitted = FALSE, converged = NA,
                              loglik = NA_real_, npar = NA_real_, ibs_cv = NA_real_, ipa_cv = NA_real_, auc_cv = NA_real_,
                              calib_cv = NA_real_, np_distance = NA_real_, folds_failed = NA_integer_,
                              checks = "could not be fitted", stringsAsFactors = FALSE)
      next
    }
    fits[[m]] <- full
    pooled <- matrix(NA_real_, n, length(at)); failed <- 0L
    for (k in seq_len(folds)) {
      te <- fold == k
      fk <- tryCatch(suppressWarnings(cd$fit(formula, data[!te, , drop = FALSE])), error = function(e) NULL)
      rk <- if (is.null(fk)) NULL else tryCatch(cd$risk(fk, data[te, , drop = FALSE], at), error = function(e) NULL)
      if (!risk_ok(rk, sum(te))) { failed <- failed + 1L; next }
      pooled[te, ] <- rk
      bk <- .scan_brier_curve(f$time[te], f$status[te], rk[, ig, drop = FALSE], grid, cause)
      foldibs[m, k] <- .scan_integrate(grid, bk[, 1L])
    }
    ll <- .scan_loglik(full)
    ibs <- ipa <- aucm <- calm <- NA_real_
    if (failed == 0L) {
      bc <- .scan_brier_curve(f$time, f$status, pooled[, ig, drop = FALSE], grid, cause)
      ibs <- .scan_integrate(grid, bc[, 1L]); ibs0 <- .scan_integrate(grid, bc[, 2L])
      ipa <- 1 - ibs / ibs0
      bcurve[[m]] <- data.frame(model = cd$name, time = grid, brier_cv = bc[, 1L], brier_null = bc[, 2L])
      pt <- NULL
      for (j in seq_along(times)) {
        a <- tryCatch(auc_t(f$time, f$status, pooled[, it[j]], times[j], cause), error = function(e) NA_real_)
        ce <- tryCatch(attr(calibration(f$time, f$status, pooled[, it[j]], times[j], cause, 10L), "error"), error = function(e) NA_real_)
        b <- tryCatch(brier_score(f$time, f$status, pooled[, it[j]], times[j], cause), error = function(e) c(NA_real_, NA_real_, NA_real_))
        pt <- rbind(pt, data.frame(model = cd$name, time = times[j], brier_cv = b[[1L]], ipa_cv = b[[3L]], auc_cv = a, calib_cv = ce))
      }
      per_time[[m]] <- pt
      aucm <- mean(pt$auc_cv, na.rm = TRUE); calm <- mean(pt$calib_cv, na.rm = TRUE)
    }
    npd <- max(abs(colMeans(rfull[, ig, drop = FALSE]) - np_grid[ig]))
    rows[[m]] <- data.frame(model = cd$name, family = cd$family, likelihood = cd$lik, fitted = TRUE, converged = .scan_converged(full),
                            loglik = ll[1L], npar = ll[2L], ibs_cv = ibs, ipa_cv = ipa, auc_cv = aucm, calib_cv = calm,
                            np_distance = npd, folds_failed = failed, checks = .scan_checks(full, cd, n), stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  # AIC within the group of models with a full likelihood of the same data
  full_lik <- tab$likelihood == "full" & tab$fitted & is.finite(tab$loglik)
  tab$aic <- ifelse(full_lik, -2 * tab$loglik + 2 * tab$npar, NA_real_)
  tab$delta_aic <- if (any(full_lik)) tab$aic - min(tab$aic, na.rm = TRUE) else NA_real_
  # ranking on the cross-validated integrated Brier score, with the paired standard error of the
  # difference from the best model over the folds (one-standard-error rule, Hastie et al. 2009)
  ok <- is.finite(tab$ibs_cv) & tab$fitted
  tab$rank <- NA_integer_; tab$diff_best <- NA_real_; tab$se_diff <- NA_real_; tab$equivalent <- NA
  if (any(ok)) {
    best <- which(ok)[which.min(tab$ibs_cv[ok])]
    for (m in which(ok)) {
      dk <- foldibs[m, ] - foldibs[best, ]
      dk <- dk[is.finite(dk)]
      tab$diff_best[m] <- tab$ibs_cv[m] - tab$ibs_cv[best]
      tab$se_diff[m] <- if (m == best) 0 else if (length(dk) >= 2L) stats::sd(dk) / sqrt(length(dk)) else NA_real_
      tab$equivalent[m] <- m == best || (is.finite(tab$se_diff[m]) && tab$diff_best[m] <= tab$se_diff[m])
    }
    eqv <- ok & tab$equivalent %in% TRUE
    tab$rank[eqv] <- 1L
    rest <- which(ok & !eqv)
    tab$rank[rest[order(tab$ibs_cv[rest])]] <- sum(eqv) + seq_along(rest)
  }
  o <- order(is.na(tab$rank), tab$rank, tab$ibs_cv)
  tab <- tab[o, , drop = FALSE]; rownames(tab) <- NULL
  effects <- NULL
  if (goal == "effect") {
    for (m in seq_len(M)) {
      if (is.null(fits[[m]])) next
      e <- tryCatch(estimates(fits[[m]]), error = function(e) NULL)
      if (is.null(e) || !"role" %in% names(e)) next
      e <- e[e$role == "beta" & (is.na(e$cause) | e$cause == cause), , drop = FALSE]
      if (!nrow(e)) next
      effects <- rbind(effects, data.frame(model = cand[[m]]$name, term = e$term, scale = sub("^log ", "", e$scale), estimate = e$exp_estimate,
                                           lower = e$exp_lower, upper = e$exp_upper, p = e$p, stringsAsFactors = FALSE))
    }
  }
  counts <- table(factor(f$status, levels = 0:K))
  ModelScan(model_call = cl, type = type, K = as.integer(K), cause = cause, goal = goal, depth = depth,
            times = times, grid = grid, folds = folds, seed = as.integer(seed), n = as.integer(n),
            counts = as.integer(counts), plateau = plateau, table = tab,
            per_time = if (length(per_time)) do.call(rbind, per_time) else NULL,
            brier = if (length(bcurve)) do.call(rbind, bcurve) else NULL, effects = effects, fits = fits,
            nonparametric = data.frame(time = at, risk = np_grid), data = data, time = as.double(f$time), status = as.integer(f$status))
}

# ----------------------------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------------------------
.scan_lines <- function(x, top = 10L, detail = FALSE) {
  tab <- x@table
  what <- if (x@type == "single") "one event type" else paste0("competing risks, ", x@K, " causes, cause of interest ", x@cause)
  cnt <- x@counts
  ev_text <- if (x@type == "single") paste0("events: ", cnt[2L]) else paste0(paste0("cause ", seq_len(x@K), ": ", cnt[-1L]), collapse = "   ")
  out <- c(paste0("Model scan: ", what),
           paste0("  Subjects: ", x@n, "   censored: ", cnt[1L], "   ", ev_text),
           paste0("  Models fitted: ", sum(tab$fitted), " of ", nrow(tab), "   depth: ", x@depth, "   ", x@folds, "-fold cross-validation"),
           sprintf("  Integrated Brier score over (0, %s]; AUC and calibration error at t = %s", format(signif(max(x@grid), 3L)),
                   paste(format(x@times), collapse = ", ")))
  if (x@type == "single" && x@plateau) out <- c(out, "  The Kaplan-Meier curve levels off: cure models were included.")
  sh <- tab[tab$fitted, , drop = FALSE]
  if (!detail && nrow(sh) > top) sh <- sh[seq_len(top), , drop = FALSE]
  if (nrow(sh)) {
    df <- data.frame(rank = ifelse(is.na(sh$rank), "-", as.character(sh$rank)), model = sh$model,
                     IBS = .num(sh$ibs_cv, 4L), `diff (SE)` = ifelse(is.finite(sh$diff_best), paste0(.num(sh$diff_best, 4L), " (", .num(sh$se_diff, 4L), ")"), "-"),
                     IPA = .num(sh$ipa_cv, 3L), AUC = .num(sh$auc_cv, 3L), calib = .num(sh$calib_cv, 3L),
                     dAIC = .num(sh$delta_aic, 1L), checks = sh$checks, check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[1L] <- "rank"
    out <- c(out, "", "Cross-validated prediction (lower IBS and calibration error, higher IPA and AUC are better)", .table_lines(df))
    if (!detail && sum(tab$fitted) > top) out <- c(out, paste0("  ... and ", sum(tab$fitted) - top, " more (summary() shows all)"))
  }
  nf <- tab$model[!tab$fitted]
  if (length(nf)) out <- c(out, paste0("  Not fitted: ", paste(nf, collapse = ", ")))
  eq <- tab$model[tab$equivalent %in% TRUE]
  out <- c(out, "", "Reading")
  if (length(eq) > 1L)
    out <- c(out, paste0("  ", length(eq), " models predict equally well (rank 1: within one standard error of the best):"),
             paste0("    ", paste(eq, collapse = "; ")),
             "  Choose among them on interpretation, assumptions and simplicity, not on the small differences.")
  else if (length(eq) == 1L) out <- c(out, paste0("  Best prediction: ", eq, "; the others are more than one standard error worse."))
  best_ipa <- suppressWarnings(max(tab$ipa_cv, na.rm = TRUE))
  if (is.finite(best_ipa) && best_ipa < 0.01)
    out <- c(out, sprintf(paste0("  No model predicts clearly better than no covariates (largest IPA %s): the covariates carry little ",
                                 "prognostic information here. Compare the models on their effects (goal = \"effect\"), not on prediction."),
                          .num(best_ipa, 3L)))
  eqa <- tab[tab$equivalent %in% TRUE & is.finite(tab$delta_aic), , drop = FALSE]
  if (nrow(eqa) > 1L && any(eqa$delta_aic > 10))
    out <- c(out, paste0("  Within rank 1, AIC clearly disfavours ", paste(eqa$model[eqa$delta_aic > 10], collapse = "; "),
                         " (dAIC > 10): similar predictions, but a worse description of the data."))
  if (any(is.finite(tab$aic)))
    out <- c(out, paste0("  dAIC compares only the models with a full likelihood of the same data (",
                         if (x@type == "single") "parametric and cure models" else "parametric cause-specific, joint and constrained direct models",
                         "); Cox, Fine-Gray and weighted direct models have no comparable likelihood."))
  if (x@type == "competing")
    out <- c(out, "  Cause-specific models describe rates; Fine-Gray and direct models the cumulative incidence itself.",
             "  For the effect of a covariate on risk use the second, on the rate of the event the first (goal = \"effect\").")
  if (!is.null(x@effects) && detail) {
    e <- x@effects
    out <- c(out, "", "Effects of the covariates by model (exp(coefficient), 95% CI)",
             .table_lines(data.frame(model = e$model, term = e$term, scale = e$scale, estimate = .num(e$estimate, 3L),
                                     ci = .ci_text(e$lower, e$upper, 3L), p = .pval(e$p), stringsAsFactors = FALSE)))
  }
  out <- c(out, "", "NOTE: confidence intervals of a model chosen from this table do not account for the choice.")
  out
}

S7::method(print, ModelScan) <- function(x, ...) { cat(.scan_lines(x, detail = FALSE), sep = "\n"); invisible(x) }
S7::method(summary, ModelScan) <- function(object, ...) { cat(.scan_lines(object, detail = TRUE), sep = "\n"); invisible(object) }
S7::method(as.data.frame, ModelScan) <- function(x, ...) x@table

.scan_top <- function(x, top) {
  tab <- x@table[x@table$fitted & is.finite(x@table$ibs_cv), , drop = FALSE]
  utils::head(tab$model, top)
}

S7::method(plot, ModelScan) <- function(x, which = c("ranking", "fit", "brier", "calibration"), top = 4L, ...) {
  which <- match.arg(which, several.ok = TRUE)
  pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  if (length(which) > 1L) graphics::par(mfrow = c(ceiling(length(which) / 2), min(2L, length(which))))
  base_mar <- c(4.1, 4.1, 2.6, 1.1)
  short <- function(s, n = 30L) ifelse(nchar(s) > n, paste0(substr(s, 1L, n - 3L), "..."), s)
  leg <- function(pos, labels, ...) graphics::legend(pos, short(labels), bty = "n", cex = 0.7, ...)
  ylab_risk <- if (x@type == "single") "Probability of the event" else paste("Cumulative incidence of cause", x@cause)
  for (w in which) {
    if (w == "ranking") {
      tab <- x@table[x@table$fitted & is.finite(x@table$ibs_cv), , drop = FALSE]
      tab <- tab[order(tab$ibs_cv, decreasing = TRUE), , drop = FALSE]
      if (!nrow(tab)) next
      lo <- tab$ibs_cv - ifelse(is.finite(tab$se_diff), tab$se_diff, 0); hi <- tab$ibs_cv + ifelse(is.finite(tab$se_diff), tab$se_diff, 0)
      lab <- short(tab$model, 34L)
      graphics::par(mar = c(4.1, min(16, max(6, max(nchar(lab)) * 0.42)), 2.6, 1.1))
      graphics::plot(tab$ibs_cv, seq_len(nrow(tab)), xlim = range(c(lo, hi)), yaxt = "n", pch = 19,
                     col = ifelse(tab$equivalent %in% TRUE, pal[1L], "grey40"), xlab = "Integrated Brier score",
                     ylab = "", main = "Ranking")
      graphics::segments(lo, seq_len(nrow(tab)), hi, seq_len(nrow(tab)), col = "grey60")
      graphics::axis(2, at = seq_len(nrow(tab)), labels = lab, las = 1, cex.axis = 0.7)
      graphics::abline(v = min(tab$ibs_cv), lty = 3, col = "grey60")
      graphics::par(mar = base_mar)
    } else if (w == "fit") {
      graphics::par(mar = base_mar)
      mods <- .scan_top(x, top)
      tg <- seq(0, max(x@grid), length.out = 60L)[-1L]
      npv <- .scan_np(x@time, x@status, x@cause, tg)
      avg <- lapply(mods, function(m) { rk <- .scan_risk_any(x, x@fits[[m]], m, tg); if (is.null(rk)) NULL else colMeans(rk) })
      ymax <- max(c(npv, unlist(avg), 1e-6), na.rm = TRUE) * 1.1
      graphics::plot(tg, npv, type = "n", xlim = c(0, max(x@grid)), ylim = c(0, ymax),
                     xlab = "Time", ylab = ylab_risk, main = "Fit: average prediction")
      graphics::lines(c(0, tg), c(0, npv), type = "s", lwd = 2, col = "black")
      for (i in seq_along(mods)) if (!is.null(avg[[i]])) graphics::lines(tg, avg[[i]], col = pal[i], lwd = 2, lty = i + 1L)
      leg("topleft", c(if (x@type == "single") "Kaplan-Meier" else "Aalen-Johansen", mods),
          col = c("black", pal[seq_along(mods)]), lty = c(1, seq_along(mods) + 1L), lwd = 2)
    } else if (w == "brier") {
      graphics::par(mar = base_mar)
      b <- x@brier; mods <- .scan_top(x, top)
      if (is.null(b)) next
      b <- b[b$model %in% mods, , drop = FALSE]
      b0 <- stats::aggregate(brier_null ~ time, data = b, FUN = mean)
      graphics::plot(b0$time, b0$brier_null, type = "l", lty = 3, lwd = 2, col = "grey40", ylim = c(0, max(c(b$brier_cv, b0$brier_null), na.rm = TRUE)),
                     xlab = "Time", ylab = "Brier score (cross-validated)", main = "Prediction error")
      for (i in seq_along(mods)) { bi <- b[b$model == mods[i], ]; graphics::lines(bi$time, bi$brier_cv, col = pal[i], lwd = 2, lty = i + 1L) }
      leg("bottomright", c("no covariates", mods), col = c("grey40", pal[seq_along(mods)]), lty = c(3, seq_along(mods) + 1L), lwd = 2)
    } else if (w == "calibration") {
      graphics::par(mar = base_mar)
      mods <- .scan_top(x, top); t0 <- x@times[ceiling(length(x@times) / 2)]
      graphics::plot(c(0, 1), c(0, 1), type = "n", xlab = "Predicted risk", ylab = "Observed risk",
                     main = paste0("Calibration, t = ", format(t0)))
      graphics::abline(0, 1, lty = 2, col = "grey50")
      for (i in seq_along(mods)) {
        rk <- .scan_risk_any(x, x@fits[[mods[i]]], mods[i], t0)
        if (is.null(rk)) next
        cal <- tryCatch(calibration(x@time, x@status, rk[, 1L], t0, x@cause, 10L), error = function(e) NULL)
        if (!is.null(cal)) graphics::points(cal$predicted, cal$observed, col = pal[i], pch = 14 + i, type = "b")
      }
      leg("topleft", mods, col = pal[seq_along(mods)], pch = 14 + seq_along(mods))
    }
  }
  invisible(x)
}

# risk of the scanned cause for the subjects of the scan, at the times tg
.scan_risk_any <- function(x, fit, name, tg) {
  if (is.null(fit)) return(NULL)
  nd <- x@data
  cand <- .scan_candidates(x@type, x@cause, "full", 1e6L, TRUE, x@K)
  cd <- Filter(function(c) c$name == name, cand)
  risk <- if (length(cd)) cd[[1L]]$risk
          else if (x@type == "single") function(fit, nd, times) .risk_surv(fit, nd, times)
          else function(fit, nd, times) as.matrix(predict(fit, nd, times = times, cause = x@cause))
  tryCatch(risk(fit, nd, tg), error = function(e) NULL)
}

# ----------------------------------------------------------------------------------------------
# report(): one self-contained HTML file
# ----------------------------------------------------------------------------------------------
#' @rdname model_scan
#' @export
report <- S7::new_generic("report", "x", function(x, file = "model_scan_report.html", ...) S7::S7_dispatch())

.b64 <- function(raw) {
  tb <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "")[[1L]]
  v <- as.integer(raw); n <- length(v); pad <- (3L - n %% 3L) %% 3L
  v <- c(v, rep(0L, pad)); m <- matrix(v, nrow = 3L)
  b <- m[1L, ] * 65536L + m[2L, ] * 256L + m[3L, ]
  ch <- rbind(tb[b %/% 262144L + 1L], tb[(b %/% 4096L) %% 64L + 1L], tb[(b %/% 64L) %% 64L + 1L], tb[b %% 64L + 1L])
  s <- paste(ch, collapse = "")
  if (pad) substr(s, nchar(s) - pad + 1L, nchar(s)) <- strrep("=", pad)
  s
}

.scan_refs <- function(keys) {
  full <- vapply(.bib[intersect(keys, names(.bib))], `[[`, character(1), 2L)
  sort(unname(full))
}

.html_table <- function(df) {
  esc <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", as.character(s))))
  h <- paste0("<tr>", paste0("<th>", esc(names(df)), "</th>", collapse = ""), "</tr>")
  b <- vapply(seq_len(nrow(df)), function(i) paste0("<tr>", paste0("<td>", esc(unlist(df[i, ])), "</td>", collapse = ""), "</tr>"), "")
  paste0("<table>", h, paste(b, collapse = ""), "</table>")
}

S7::method(report, ModelScan) <- function(x, file = "model_scan_report.html", width = 7, height = 5, ...) {
  img <- function(which) {
    tf <- tempfile(fileext = ".png")
    grDevices::png(tf, width = width, height = height, units = "in", res = 110)
    ok <- tryCatch({ plot(x, which = which); TRUE }, error = function(e) FALSE)
    grDevices::dev.off()
    if (!ok || !file.exists(tf)) return("")
    on.exit(unlink(tf))
    paste0("<img alt='", which, "' src='data:image/png;base64,", .b64(readBin(tf, "raw", file.info(tf)$size)), "'/>")
  }
  tab <- x@table
  show <- data.frame(rank = ifelse(is.na(tab$rank), "-", tab$rank), model = tab$model, family = tab$family,
                     IBS = .num(tab$ibs_cv, 4L), diff = .num(tab$diff_best, 4L), SE = .num(tab$se_diff, 4L),
                     IPA = .num(tab$ipa_cv, 3L), AUC = .num(tab$auc_cv, 3L), calibration = .num(tab$calib_cv, 3L),
                     `distance to non-parametric` = .num(tab$np_distance, 3L), dAIC = .num(tab$delta_aic, 1L),
                     checks = tab$checks, check.names = FALSE)
  refs <- .scan_refs(c("graf1999", "gerds2006", "kattan2018", "blanche2013", "gerds2014", "hastie2009", "akaike1974", "burnham2002",
                         if (x@type == "competing") c("putter2007", "latouche2013")))
  methods <- paste0(
    "<p>Every candidate model was fitted to the full data and in ", x@folds, "-fold cross-validation (folds stratified by ",
    "event status; seed ", x@seed, "). Prediction was assessed on the pooled out-of-fold predictions of the probability of ",
    if (x@type == "single") "the event" else paste0("cause ", x@cause), ": the integrated Brier score over (0, ",
    format(signif(max(x@grid), 3L)), "] with inverse probability of censoring weights (Graf et al. 1999; Gerds and Schumacher 2006), ",
    "the index of prediction accuracy (IPA, 1 minus the ratio to the model without covariates; Kattan and Gerds 2018), the ",
    "time-dependent AUC (Blanche et al. 2013) and the calibration error (size-weighted mean absolute difference between ",
    "predicted and observed risk in ten risk groups; Gerds et al. 2014) at t = ", paste(format(x@times), collapse = ", "), ". ",
    "The standard error of the difference from the best model is the standard deviation of the fold-wise differences ",
    "divided by the square root of the number of folds; models within one standard error of the best were given rank 1 ",
    "(one-standard-error rule; Hastie et al. 2009, Section 7.10). AIC (Akaike 1974) is reported only for the models with a full ",
    "likelihood of the same data. The distance to the non-parametric estimate is the largest absolute difference between the ",
    "average predicted risk and the ", if (x@type == "single") "Kaplan-Meier" else "Aalen-Johansen", " estimate on the evaluation grid.</p>")
  html <- c("<!DOCTYPE html><html><head><meta charset='utf-8'><title>SurvX model scan</title>",
            "<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;color:#222;line-height:1.5}",
            "table{border-collapse:collapse;font-size:13px;margin:1em 0}th,td{border-bottom:1px solid #ddd;padding:4px 8px;text-align:right}",
            "th:nth-child(2),td:nth-child(2),td:last-child,th:last-child{text-align:left}pre{background:#f6f6f6;padding:1em;overflow-x:auto}",
            "img{max-width:100%}</style></head><body>",
            "<h1>Model scan</h1>", paste0("<p>SurvX ", as.character(utils::packageVersion("SurvX")), ", ", format(Sys.time(), "%Y-%m-%d %H:%M"), "</p>"),
            "<h2>Summary</h2>", paste0("<pre>", paste(gsub("<", "&lt;", .scan_lines(x, detail = TRUE)), collapse = "\n"), "</pre>"),
            "<h2>All models</h2>", .html_table(show),
            "<h2>Figures</h2>", img("ranking"), img("fit"), img("brier"), img("calibration"),
            if (!is.null(x@per_time)) c("<h2>Measures at each time</h2>", .html_table(within(x@per_time, {
              brier_cv <- .num(brier_cv, 4L); ipa_cv <- .num(ipa_cv, 3L); auc_cv <- .num(auc_cv, 3L); calib_cv <- .num(calib_cv, 3L) }))),
            if (!is.null(x@effects)) c("<h2>Effects of the covariates</h2>", .html_table(within(x@effects, {
              estimate <- .num(estimate, 3L); lower <- .num(lower, 3L); upper <- .num(upper, 3L); p <- .pval(p) }))),
            "<h2>Methods</h2>", methods,
            "<h2>References</h2>", paste0("<p>", refs, "</p>"),
            "</body></html>")
  writeLines(html, file, useBytes = TRUE)
  message("report written to ", normalizePath(file, mustWork = FALSE))
  invisible(file)
}
