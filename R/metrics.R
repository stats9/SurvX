#' Prediction performance at a time horizon
#'
#' @description
#' Measures how good predicted risks are at a time horizon \eqn{t}, for one event type or one cause
#' with competing risks, under right censoring:
#' \itemize{
#'   \item \code{brier_score()}: the Brier score (mean squared error of the predicted risk) and the
#'     index of prediction accuracy (IPA), which compares it with a model without covariates;
#'   \item \code{auc_t()}: the time-dependent area under the ROC curve (discrimination: are
#'     subjects who have the event by \eqn{t} given higher risks?);
#'   \item \code{calibration()}: predicted against observed risk in risk groups (do predicted
#'     risks of 20\% come true for 20\% of subjects?).
#' }
#' \strong{Use them} to validate a prediction model, ideally on new data or with cross-validation
#' (\code{\link{compare_models}}, \code{\link{cv_compare}}), which call these functions.
#' @param time,status follow-up time and status (0 = censored, 1..K = cause) of the subjects whose
#'   risks are assessed.
#' @param risk predicted probability of the event (or of the cause of interest) by time \code{t}, one
#'   per subject: for example \code{predict(fit, newdata, times = t)} of a Fine-Gray or
#'   cause-specific model, or \code{1 - predict(fit, newdata, times = t, type = "survival")} of a
#'   survival model.
#' @param t the time horizon (a single number).
#' @param cause cause of interest.
#' @param groups number of risk groups for calibration.
#' @details
#' \subsection{Brier score}{
#' \eqn{BS(t) = n^{-1} \sum_i \hat W_i(t) \{I(T_i \le t, \mathrm{cause}_i = 1) - \hat r_i\}^2},
#' with inverse-probability-of-censoring weights \eqn{\hat W_i(t)} from the Kaplan-Meier estimate of
#' the censoring distribution: \eqn{1/\hat G(T_i-)} for subjects with an event (of any cause) by
#' \eqn{t}, \eqn{1/\hat G(t)} for those still at risk, and 0 for those censored before \eqn{t}
#' (Graf et al. 1999; Gerds and Schumacher 2006). The null model predicts the Aalen-Johansen
#' (Kaplan-Meier) estimate for everybody, and \eqn{IPA = 1 - BS/BS_{null}} (Kattan and Gerds 2018).
#' }
#' \subsection{AUC}{
#' The probability that a subject with the event of interest by \eqn{t} has a higher predicted risk
#' than a subject without it, with the same weights; subjects with a competing event before \eqn{t}
#' count as controls ("definition 2" of Blanche et al. 2013). Ties in the predicted risk count one
#' half.
#' }
#' \subsection{Calibration}{
#' Subjects are grouped by predicted risk (\code{groups} quantile groups); in each group the mean
#' predicted risk is compared with the Aalen-Johansen estimate of the observed risk at \eqn{t}
#' (Gerds, Andersen and Kattan 2014). The calibration error is the size-weighted mean absolute
#' difference.
#' }
#' @section Interpreting the result:
#' \itemize{
#'   \item Brier score: 0 is perfect; a model that predicts 50\% for everybody scores 0.25. It
#'     depends on how common the event is, so compare it with the null model: IPA > 0 means the
#'     model is better than ignoring the covariates, IPA = 0.10 means 10\% of the null model's
#'     error is removed. Negative IPA: the model is worse than no model.
#'   \item AUC: 0.5 is no discrimination, 1 perfect; 0.7 to 0.8 is common for clinical models. AUC
#'     ignores calibration: a model can rank well and still predict wrong risks.
#'   \item Calibration: points near the diagonal (mean predicted = observed) and a small error.
#' }
#' Measured on the data used for fitting, all three are optimistic; use cross-validation or
#' external data for an honest assessment.
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}), the Brier score, that of the null model and the AUC agree with
#' \code{riskRegression::Score} to \eqn{2 \times 10^{-14}} on 82 data sets with and without competing
#' risks. They are about 10 times faster than \code{Score} on \eqn{n = 10000}.
#' @return
#' \code{brier_score}: named vector with \code{brier}, \code{brier_null} (Aalen-Johansen
#'   prediction for everyone) and \code{ipa = 1 - brier / brier_null}. \code{auc_t}: the AUC.
#'   \code{calibration}: a data frame with mean predicted and observed (Aalen-Johansen) risk per
#'   group; its attribute \code{"error"} is the size-weighted mean absolute difference.
#' @references
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999) Assessment and comparison of prognostic
#'   classification schemes for survival data. Statistics in Medicine 18:2529-2545.
#'
#' Gerds TA, Schumacher M (2006) Consistent estimation of the expected Brier score in general survival
#'   models with right-censored event times. Biometrical Journal 48:1029-1040.
#'
#' Kattan MW, Gerds TA (2018) The index of prediction accuracy: an intuitive measure useful for
#'   evaluating risk prediction models. Diagnostic and Prognostic Research 2:7.
#'
#' Blanche P, Dartigues J-F, Jacqmin-Gadda H (2013) Estimating and comparing time-dependent areas
#'   under receiver operating characteristic curves for censored event times with competing risks.
#'   Statistics in Medicine 32:5381-5397.
#'
#' Gerds TA, Andersen PK, Kattan MW (2014) Calibration plots for risk prediction models in the
#'   presence of competing risks. Statistics in Medicine 33:3191-3203.
#'
#' All references: \code{sx_references("brier_score")}, \code{sx_references("auc_t")},
#' \code{sx_references("calibration")}.
#' @seealso
#' \code{\link{compare_models}} and \code{\link{cv_compare}} (cross-validated comparison),
#' \code{\link{predict.SurvX}}.
#' @examples
#' set.seed(3)
#' n <- 400
#' d <- data.frame(x = rnorm(n))
#' t1 <- rexp(n, 0.3 * exp(0.8 * d$x)); t2 <- rexp(n, 0.2); cens <- runif(n, 0, 5)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#'
#' ## predicted 2-year risk of cause 1 from a Fine-Gray model
#' fit <- fg(Event(time, status) ~ x, d)
#' risk <- predict(fit, d, times = 2)[, 1]
#'
#' brier_score(d$time, d$status, risk, t = 2)     # brier, brier_null, ipa
#' auc_t(d$time, d$status, risk, t = 2)
#' cal <- calibration(d$time, d$status, risk, t = 2, groups = 5)
#' cal
#' plot(cal$predicted, cal$observed, xlim = c(0, 1), ylim = c(0, 1),
#'      xlab = "predicted risk", ylab = "observed risk"); abline(0, 1, lty = 2)
#' @export
brier_score <- function(time, status, risk, t, cause = 1L) {
  n <- length(time)
  if (length(risk) != n) stop("risk must have one value per subject", call. = FALSE)
  r <- .C("SX_brier", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.double(risk), as.double(t),
          brier = double(1), null = double(1), ipa = double(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "brier_score")
  c(brier = r$brier, brier_null = r$null, ipa = r$ipa)
}

#' @rdname brier_score
#' @export
auc_t <- function(time, status, risk, t, cause = 1L) {
  n <- length(time)
  if (length(risk) != n) stop("risk must have one value per subject", call. = FALSE)
  r <- .C("SX_auc", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.double(risk), as.double(t),
          auc = double(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "auc_t")
  r$auc
}

#' @rdname brier_score
#' @export
calibration <- function(time, status, risk, t, cause = 1L, groups = 10L) {
  n <- length(time)
  if (length(risk) != n) stop("risk must have one value per subject", call. = FALSE)
  r <- .C("SX_calibration", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.double(risk), as.double(t),
          as.integer(groups), pred = double(groups), obs = double(groups), size = double(groups), error = double(1),
          m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "calibration")
  k <- seq_len(r$m)
  structure(data.frame(group = k, n = r$size[k], predicted = r$pred[k], observed = r$obs[k]), error = r$error)
}

#' Compare direct, Fine-Gray and cause-specific Cox models
#'
#' @description
#' Fits a set of competing-risks models to the same data and compares how well they predict the
#' cumulative incidence of one cause: direct parametric models (\code{\link{directcif}}, any
#' estimator), the Fine-Gray model (\code{\link{fg}}), cause-specific Cox models
#' (\code{\link{csc}}) and, optionally, parametric cause-specific hazard models. For each model it
#' reports the apparent and K-fold cross-validated Brier score, index of prediction accuracy (IPA),
#' time-dependent AUC and calibration error at the chosen times, and the AIC where it is comparable.
#' Fitting, cross-validation and scoring are done in the C++ core, so the comparison is fast.
#'
#' \strong{Use it when} choosing a model for predicting the risk of a cause.
#' \strong{Use} \code{\link{cv_compare}} for models that \code{compare_models} does not include
#' (other packages, other options).
#' @param formula \code{Event(time, status) ~ covariates}.
#' @param data a data frame.
#' @param times time horizons.
#' @param direct a list of direct models, each a character vector \code{c(baseline, link)},
#'   \code{c(baseline, link, estimation)} or \code{c(baseline, link, "constrained", cond)}; link
#'   \code{"PH"}, \code{"PO"} or \code{"GOR"}; estimation \code{"joint"} (default), \code{"weighted"} or
#'   \code{"constrained"}; see \code{\link{directcif}}. \code{NULL}: Gompertz baseline with PH and PO
#'   links, joint and weighted, and with two causes the constrained PH model with a Weibull \eqn{G}.
#' @param cause cause of interest.
#' @param folds number of cross-validation folds (less than 2: apparent performance only).
#' @param seed seed of the fold assignment.
#' @param df degrees of freedom of spline baselines.
#' @param hazard optional list of parametric cause-specific hazard models, fitted by \code{\link{csc}}
#'   with a parametric \code{baseline}: each a distribution name or \code{c("rp", df)}, e.g.
#'   \code{list("weibull", "lognormal", c("rp", 4))}. They are cross-validated with the same folds as
#'   the other models, and their AIC (full likelihood of the data) is comparable with that of the joint
#'   and constrained direct models.
#' @section Interpreting the result:
#' Judge the models on the cross-validated columns (ending in \code{_cv}): lower Brier score and
#' calibration error, higher IPA and AUC are better. Apparent measures, computed on the data used for
#' fitting, favour complex models. Differences in the third decimal of the Brier score are usually
#' within noise; prefer the simpler model when the measures are close. AIC can be compared only among
#' models with a full likelihood of the same data (joint and constrained direct models, parametric
#' cause-specific models); see Value.
#' @return
#' an S7 object of class \code{ModelComparison}; \code{x@table} has one row per model and
#'   time, \code{x@fit} one row per model. AIC is reported for the joint and constrained direct models
#'   (full likelihoods of the same data, comparable with each other); the weighted direct model,
#'   Fine-Gray and Cox maximise pseudo- or partial likelihoods.
#' @references
#' \emph{Cross-validation:}
#' Stone M (1974). Cross-validatory choice and assessment of statistical predictions. Journal of the
#'   Royal Statistical Society, Series B 36:111-147.
#' Hastie T, Tibshirani R, Friedman J (2009). The Elements of Statistical Learning, 2nd ed.
#'   Springer, Section 7.10.
#' Gerds TA, Schumacher M (2007). Efron-type measures of prediction error for survival analysis.
#'   Biometrics 63:1283-1287.
#'
#' \emph{Measures:} see \code{\link{brier_score}}.
#'
#' All references: \code{sx_references("compare_models")}.
#' @seealso
#' \code{\link{cv_compare}}, \code{\link{brier_score}}, \code{\link{directcif}}, \code{\link{fg}},
#' \code{\link{csc}}.
#' @examples
#' set.seed(3)
#' n <- 300
#' d <- data.frame(x = rbinom(n, 1, 0.5))
#' t1 <- rexp(n, 0.4 * exp(0.6 * d$x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#'
#' cm <- compare_models(Event(time, status) ~ x, d, times = c(1, 2), folds = 3,
#'                      direct = list(c("gompertz", "PH"), c("gompertz", "PO", "weighted")),
#'                      hazard = list("weibull"))
#' cm
#' as.data.frame(cm)
#' @export
compare_models <- function(formula, data, times, direct = NULL, cause = 1L, folds = 5L, seed = 1L, df = 3L, hazard = NULL) {
  cl <- match.call()
  f <- .sx_frame(formula, data)
  n <- length(f$time); p <- ncol(f$X); K <- max(f$status); nt <- length(times)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  if (is.null(direct)) {
    direct <- list(c("gompertz", "PH"), c("gompertz", "PO"), c("gompertz", "PH", "weighted"), c("gompertz", "PO", "weighted"))
    if (K == 2L) direct <- c(direct, list(c("gompertz", "PH", "constrained", "weibull")))
  }
  nd <- length(direct)
  part <- function(i, default) vapply(direct, function(d) if (length(d) >= i) as.character(d[[i]]) else default, character(1))
  base <- part(1L, "gompertz"); lk <- part(2L, "PH"); est <- part(3L, "joint"); cd <- part(4L, "logistic")
  if (any(!base %in% names(.base_codes)) || any(!lk %in% c("PH", "PO", "GOR")) || any(!est %in% c("joint", "weighted", "constrained")) ||
      any(!cd %in% names(.cond_codes)))
    stop("direct: each entry must be c(baseline, link), c(baseline, link, estimation) or c(baseline, link, \"constrained\", cond)", call. = FALSE)
  if (any(est == "constrained") && K != 2L) stop("constrained models need exactly two causes", call. = FALSE)
  if (any(est == "constrained" & base %in% c("weibull", "spline"))) stop("constrained models need a baseline with a plateau (not weibull or spline)", call. = FALSE)
  M <- nd + 2L
  r <- .C("SX_compare", as.integer(n), as.integer(p), as.integer(K), f$time, f$status, .rm(f$X), as.integer(cause),
          as.integer(nt), as.double(times), as.integer(folds), as.integer(seed), as.integer(nd),
          as.integer(unname(.base_codes[base])), as.integer(match(lk, c("PH", "PO", "GOR")) - 1L), rep(as.integer(df), nd),
          as.integer(match(est, c("joint", "weighted", "constrained")) - 1L), as.integer(unname(.cond_codes[cd])),
          metrics = double(M * nt * 8L), fit = double(M * 5L), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "compare_models")
  suffix <- ifelse(est == "joint", "", ifelse(est == "weighted", ", weighted", paste0(", constrained (G ", cd, ")")))
  names_m <- c(paste0("direct ", base, " ", lk, suffix), "Fine-Gray", "cause-specific Cox")
  a <- array(r$metrics, dim = c(2L, 4L, nt, M))                      # [apparent/cv, metric, time, model]
  tab <- NULL
  for (m in seq_len(M)) for (j in seq_len(nt))
    tab <- rbind(tab, data.frame(model = names_m[m], time = times[j],
                                 brier = a[1L, 1L, j, m], brier_cv = a[2L, 1L, j, m], ipa = a[1L, 2L, j, m], ipa_cv = a[2L, 2L, j, m],
                                 auc = a[1L, 3L, j, m], auc_cv = a[2L, 3L, j, m], calib_error = a[1L, 4L, j, m], calib_error_cv = a[2L, 4L, j, m]))
  ft <- .matrix_rm(r$fit, M, 5L)
  fit <- data.frame(model = names_m, loglik = ft[, 1L], npar = ft[, 2L], aic = ft[, 3L], flags = ft[, 4L], folds_not_converged = ft[, 5L])
  fit$aic[!is.finite(fit$aic)] <- NA_real_
  if (!is.null(hazard)) {
    hz <- .hazard_specs(hazard)
    vars <- intersect(all.vars(formula), names(data))
    dd <- data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
    if (nrow(dd) != n) stop("compare_models: could not align the data for the parametric hazard models", call. = FALSE)
    fold <- rep(1L, n)
    if (folds >= 2L) {                      # the folds of the C++ comparison above
      rf <- .C("SX_cv_folds", as.integer(n), as.integer(folds), as.integer(seed), out = integer(n), rc = integer(1), PACKAGE = "SurvX")
      .sx_rc(rf$rc, "compare_models")
      fold <- rf$out + 1L
    }
    funs <- lapply(hz, function(h) { force(h); function(fm, d) csc(fm, d, baseline = h$baseline, df = h$df) })
    names(funs) <- vapply(hz, `[[`, "", "name")
    res <- .cv_run(formula, dd, f, funs, times, cause, folds, fold, 10L)
    tab <- rbind(tab, res$table)
    for (m in seq_along(funs)) {
      ft <- res$fits[[m]]
      ll <- if (is.null(ft)) NA_real_ else as.numeric(logLik(ft))
      np <- if (is.null(ft)) NA_real_ else attr(logLik(ft), "df")
      ok <- !is.null(ft) && all(vapply(ft@models, function(x) isTRUE(x@converged), logical(1)))
      fit <- rbind(fit, data.frame(model = names(funs)[m], loglik = ll, npar = np, aic = -2 * ll + 2 * np,
                                   flags = if (ok) 0 else 1, folds_not_converged = res$failed[m]))
    }
  }
  ModelComparison(call = cl, table = tab, fit = fit, folds = as.integer(folds), cause = as.integer(cause))
}

#' Cross-validated comparison of any set of models
#'
#' @description
#' Fits every model to the full data (apparent performance) and in K-fold cross-validation, and
#' reports the Brier score, the index of prediction accuracy, the time-dependent AUC and the
#' calibration error of the predicted risk of one cause at the chosen times.
#' \code{\link{compare_models}} does this in the C++ core for the models of this package;
#' \code{cv_compare} is for anything else that has a \code{predict(fit, newdata, times, cause)} method
#' (a model from another package wrapped in a function, a model with its own options).
#' @param formula \code{Event(time, status) ~ covariates}; passed to every model.
#' @param data a data frame; rows with missing values in the variables of \code{formula} are dropped.
#' @param models a named list of functions \code{function(formula, data)} that return a fitted model
#'   with a \code{predict(fit, newdata, times, cause)} method giving the cumulative incidence of
#'   the cause. Wrap any model, with its own options, in such
#'   a function.
#' @param times time horizons.
#' @param cause cause of interest.
#' @param folds number of folds (less than 2: apparent performance only).
#' @param seed seed of the fold assignment (the state of the random number generator is restored).
#' @param groups number of risk groups for the calibration error.
#' @section Interpreting the result:
#' As for \code{\link{compare_models}}: judge the models on the cross-validated measures (lower Brier
#' score and calibration error, higher IPA and AUC). The same folds are used for every model, so the
#' comparison is paired.
#' @return
#' an S7 object of class \code{ModelComparison}: \code{x@table} has one row per model and time
#'   (columns ending in \code{_cv} are computed from the pooled out-of-fold predictions), \code{x@fit}
#'   one row per model. A model that fails in a fold has \code{NA} cross-validated measures.
#' @references
#' Gerds TA, Schumacher M (2007) Efron-type measures of prediction error for survival analysis.
#'   Biometrics 63:1283-1287.
#'
#' Hastie T, Tibshirani R, Friedman J (2009) The Elements of Statistical Learning, 2nd ed. Springer,
#'   Section 7.10 (cross-validation).
#'
#' Akaike H (1974) A new look at the statistical model identification. IEEE Transactions on Automatic
#'   Control 19:716-723.
#'
#' All references: \code{sx_references("cv_compare")}.
#' @seealso
#' \code{\link{compare_models}}, \code{\link{brier_score}}.
#' @examples
#' set.seed(3)
#' n <- 300
#' x <- rbinom(n, 1, 0.5)
#' t1 <- rexp(n, 0.4 * exp(0.6 * x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
#' d <- data.frame(time = pmin(t1, t2, cens), x = x)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#' cv_compare(Event(time, status) ~ x, d, times = c(1, 2), folds = 3, models = list(
#'   "direct PH, weighted" = function(f, d) directcif(f, d, estimation = "weighted"),
#'   "Fine-Gray" = function(f, d) fg(f, d),
#'   "cause-specific Cox" = function(f, d) csc(f, d)))
#' @export
cv_compare <- function(formula, data, models, times, cause = 1L, folds = 5L, seed = 1L, groups = 10L) {
  cl <- match.call()
  if (!is.list(models) || is.null(names(models)) || any(names(models) == "")) stop("models must be a named list of functions", call. = FALSE)
  vars <- intersect(all.vars(formula), names(data))
  data <- data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  f <- .sx_frame(formula, data)
  n <- length(f$time); nt <- length(times); M <- length(models)
  if (n != nrow(data)) stop("could not align the response with the data", call. = FALSE)
  fold <- rep(1L, n)
  if (folds >= 2L) {
    had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed)
    fold <- sample(rep_len(seq_len(folds), n))
    if (had) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
  }
  res <- .cv_run(formula, data, f, models, times, cause, folds, fold, groups)
  tab <- res$table; fit <- NULL
  for (m in seq_len(M)) {
    ft <- res$fits[[m]]
    is_direct <- !is.null(ft) && S7::S7_inherits(ft, DirectCIFFit)
    is_par <- !is.null(ft) && (S7::S7_inherits(ft, ParCSCFit) || S7::S7_inherits(ft, ParSurvFit))
    ll <- if (is_direct) ft@loglik else if (is_par) as.numeric(logLik(ft)) else NA_real_
    np <- if (is_direct) length(ft@coefficients) else if (is_par) attr(logLik(ft), "df") else NA_integer_
    fit <- rbind(fit, data.frame(
      model = names(models)[m],
      loglik = ll,
      npar = np,
      aic = if ((is_direct && ft@estimation != "weighted") || is_par) -2 * ll + 2 * np else NA_real_,
      converged = if (is_direct) ft@converged else if (S7::S7_inherits(ft, ParCSCFit)) all(vapply(ft@models, function(x) isTRUE(x@converged), logical(1)))
                  else if (is_par) ft@converged else NA,
      folds_failed = res$failed[m]))
  }
  ModelComparison(call = cl, table = tab, fit = fit, folds = as.integer(folds), cause = as.integer(cause))
}

# Apparent and K-fold cross-validated measures of a named list of fitting functions; fold: the fold of
# every row of data (1..folds). Returns the table, the full-data fits and the number of failed folds.
.cv_run <- function(formula, data, f, models, times, cause, folds, fold, groups = 10L) {
  n <- length(f$time); nt <- length(times); M <- length(models)
  risk_of <- function(fun, train, test) {
    fit <- tryCatch(fun(formula, train), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    r <- tryCatch(stats::predict(fit, newdata = test, times = times, cause = cause), error = function(e) NULL)
    if (is.null(r) || !is.matrix(r) || nrow(r) != nrow(test) || any(!is.finite(r))) return(NULL)
    list(fit = fit, risk = r)
  }
  measures <- function(risk) {
    out <- matrix(NA_real_, nt, 4L)
    if (is.null(risk)) return(out)
    for (j in seq_len(nt)) {
      b <- tryCatch(brier_score(f$time, f$status, risk[, j], times[j], cause), error = function(e) c(NA_real_, NA_real_, NA_real_))
      a <- tryCatch(auc_t(f$time, f$status, risk[, j], times[j], cause), error = function(e) NA_real_)
      ce <- tryCatch(attr(calibration(f$time, f$status, risk[, j], times[j], cause, groups), "error"), error = function(e) NA_real_)
      out[j, ] <- c(b[[1L]], b[[3L]], a, ce)
    }
    out
  }
  tab <- NULL; fits <- vector("list", M); failed <- integer(M)
  for (m in seq_len(M)) {
    full <- risk_of(models[[m]], data, data)
    app <- measures(if (is.null(full)) NULL else full$risk)
    cvm <- matrix(NA_real_, nt, 4L)
    if (folds >= 2L) {
      pooled <- matrix(NA_real_, n, nt)
      for (k in seq_len(folds)) {
        te <- fold == k
        r <- risk_of(models[[m]], data[!te, , drop = FALSE], data[te, , drop = FALSE])
        if (is.null(r)) failed[m] <- failed[m] + 1L else pooled[te, ] <- r$risk
      }
      if (failed[m] == 0L) cvm <- measures(pooled)
    }
    tab <- rbind(tab, data.frame(model = names(models)[m], time = times,
                                 brier = app[, 1L], brier_cv = cvm[, 1L], ipa = app[, 2L], ipa_cv = cvm[, 2L],
                                 auc = app[, 3L], auc_cv = cvm[, 3L], calib_error = app[, 4L], calib_error_cv = cvm[, 4L]))
    if (!is.null(full)) fits[[m]] <- full$fit
  }
  list(table = tab, fits = fits, failed = failed)
}

# hazard = list("weibull", "lognormal", c("rp", 4), "pwexp"): one parametric cause-specific model each
.hazard_specs <- function(hazard) {
  if (is.character(hazard)) hazard <- as.list(hazard)
  lapply(hazard, function(h) {
    h <- as.character(h)
    b <- .par_dist(h[1L])
    dfv <- if (length(h) >= 2L) as.integer(h[2L]) else 3L
    list(baseline = b, df = dfv,
         name = paste0("parametric ", .par_label[[b]], if (b == "rp") paste0(" (df ", dfv, ")") else "", " cause-specific"))
  })
}
