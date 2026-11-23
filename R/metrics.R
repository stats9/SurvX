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
#'     risks of 20\% come true for 20\% of subjects?), the observed/expected ratio and the
#'     integrated calibration index of a smooth calibration curve.
#' }
#' With \code{se = TRUE} the Brier score and the AUC come with standard errors. Related:
#' \code{\link{cindex}} (concordance over the whole follow-up) and \code{\link{decision_curve}}
#' (net benefit).
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
#' @param se also return the standard error (see Details).
#' @param level confidence level of the interval of the AUC and of the observed/expected ratio.
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
#' \subsection{Standard errors}{
#' With \code{se = TRUE} the standard error of the Brier score and of the AUC is
#' \eqn{\mathrm{sd}(\psi_i)/\sqrt{n}}, where \eqn{\psi_i} is the influence function of subject
#' \eqn{i} with the censoring weights taken as known (Blanche et al. 2013): the "conservative"
#' standard errors of \code{riskRegression::Score(conservative = TRUE)}. \code{\link{compare_models}}
#' and \code{\link{cv_compare}} use the same influence functions for the paired difference of two
#' models.
#' }
#' \subsection{Calibration in the large and smooth calibration}{
#' \code{calibration()} also returns the observed/expected ratio: the Aalen-Johansen (Kaplan-Meier)
#' risk at \eqn{t} over the mean predicted risk, with the interval
#' \eqn{\exp\{\log(O/E) \pm z\, se(O)/O\}}; and the integrated calibration index of Austin, Harrell
#' and van Klaveren (2020): a Cox model (with competing risks, a Fine-Gray model) of the outcome on a
#' natural cubic spline of \eqn{\log\{-\log(1 - \hat r)\}} with knots at its 10th, 50th and 90th
#' percentiles gives a calibrated risk for every subject; ICI, E50 and E90 are the mean, median and
#' 90th percentile of its absolute difference from the predicted risk.
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
#'   prediction for everyone) and \code{ipa = 1 - brier / brier_null}, and with \code{se = TRUE}
#'   \code{se}. \code{auc_t}: the AUC; with \code{se = TRUE} a vector with \code{auc}, \code{se},
#'   \code{lower} and \code{upper}. \code{calibration}: a data frame with mean predicted and
#'   observed (Aalen-Johansen) risk per group; its attribute \code{"error"} is the size-weighted mean
#'   absolute difference, \code{"oe"} the observed/expected ratio with its interval, \code{"ici"},
#'   \code{"e50"} and \code{"e90"} the smooth calibration measures and \code{"smooth"} the predicted
#'   and calibrated risks (\code{NULL} when the calibration model cannot be fitted).
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
#' Austin PC, Harrell FE, van Klaveren D (2020) Graphical calibration curves and the integrated
#'   calibration index (ICI) for survival models. Statistics in Medicine 39:2714-2742.
#'
#' All references: \code{sx_references("brier_score")}, \code{sx_references("auc_t")},
#' \code{sx_references("calibration")}.
#' @seealso
#' \code{\link{compare_models}} and \code{\link{cv_compare}} (cross-validated comparison),
#' \code{\link{cindex}}, \code{\link{decision_curve}}, \code{\link{predict.SurvX}},
#' \code{vignette("prediction", "SurvX")}.
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
#' brier_score(d$time, d$status, risk, t = 2, se = TRUE)
#' auc_t(d$time, d$status, risk, t = 2, se = TRUE)   # auc, se, lower, upper
#' cal <- calibration(d$time, d$status, risk, t = 2, groups = 5)
#' cal
#' attr(cal, "oe")                                # observed / expected, with its interval
#' c(ICI = attr(cal, "ici"), E50 = attr(cal, "e50"), E90 = attr(cal, "e90"))
#' plot(cal$predicted, cal$observed, xlim = c(0, 1), ylim = c(0, 1),
#'      xlab = "predicted risk", ylab = "observed risk"); abline(0, 1, lty = 2)
#' @export
brier_score <- function(time, status, risk, t, cause = 1L, se = FALSE) {
  n <- length(time)
  if (length(risk) != n) stop("risk must have one value per subject", call. = FALSE)
  r <- .C("SX_brier", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.double(risk), as.double(t),
          brier = double(1), null = double(1), ipa = double(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "brier_score")
  out <- c(brier = r$brier, brier_null = r$null, ipa = r$ipa)
  if (isTRUE(se)) out <- c(out, se = unname(.sx_assess(time, status, list(risk), t, cause)$a["brier_se", 1L, 1L]))
  out
}

#' @rdname brier_score
#' @export
auc_t <- function(time, status, risk, t, cause = 1L, se = FALSE, level = 0.95) {
  n <- length(time)
  if (length(risk) != n) stop("risk must have one value per subject", call. = FALSE)
  r <- .C("SX_auc", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.double(risk), as.double(t),
          auc = double(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "auc_t")
  if (!isTRUE(se)) return(r$auc)
  a <- .sx_assess(time, status, list(risk), t, cause, level = level)$a[, 1L, 1L]
  c(auc = r$auc, se = unname(a["auc_se"]), lower = unname(a["auc_lower"]), upper = unname(a["auc_upper"]))
}

#' @rdname brier_score
#' @export
calibration <- function(time, status, risk, t, cause = 1L, groups = 10L, level = 0.95) {
  n <- length(time)
  if (length(risk) != n) stop("risk must have one value per subject", call. = FALSE)
  r <- .C("SX_calibration", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.double(risk), as.double(t),
          as.integer(groups), pred = double(groups), obs = double(groups), size = double(groups), error = double(1),
          m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "calibration")
  k <- seq_len(r$m)
  a <- .sx_assess(time, status, list(risk), t, cause, groups, level)$a[, 1L, 1L]
  sm <- .C("SX_smooth_calibration", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.double(risk),
           as.double(t), observed = double(n), out = double(3), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  if (sm$rc < 0L) .sx_rc(sm$rc, "calibration")
  o <- order(risk)
  smooth <- if (sm$rc == 0L) data.frame(predicted = as.double(risk)[o], observed = sm$observed[o]) else NULL
  structure(data.frame(group = k, n = r$size[k], predicted = r$pred[k], observed = r$obs[k]), error = r$error,
            oe = unname(a[c("oe", "oe_lower", "oe_upper")]), ici = unname(a["ici"]), e50 = unname(a["e50"]), e90 = unname(a["e90"]),
            smooth = smooth)
}

#' Compare direct, Fine-Gray and cause-specific Cox models
#'
#' @description
#' Fits a set of competing-risks models to the same data and compares how well they predict the
#' cumulative incidence of one cause: direct parametric models (\code{\link{directcif}}, any
#' estimator), the Fine-Gray model (\code{\link{fg}}), cause-specific Cox models
#' (\code{\link{csc}}) and, optionally, parametric cause-specific hazard models. For each model it
#' reports the apparent and K-fold cross-validated Brier score, index of prediction accuracy (IPA),
#' time-dependent AUC and calibration error at the chosen times, and the AIC where it is comparable,
#' with the standard errors and intervals of the Brier score and AUC, the paired difference of every
#' model from the best one with its p-value, the observed/expected ratio and the integrated
#' calibration index (see the section "Columns of the table").
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
#'   links, joint and weighted, and with two or more causes the constrained PH model with a Weibull
#'   \eqn{G}. Constrained models with three or more causes are fitted by \code{directcif()} with the
#'   same folds (as those of \code{constrained}).
#' @param cause cause of interest.
#' @param folds number of cross-validation folds (less than 2: apparent performance only).
#' @param seed seed of the fold assignment: a whole number, or \code{NULL} to take the seed from the
#'   random numbers of the session (so that \code{set.seed()} before the call governs the result).
#'   All random draws come from R's generator (see \code{RNGkind()}); with a number the state of the
#'   session's generator is left as it was.
#' @param df degrees of freedom of spline baselines.
#' @param hazard optional list of parametric cause-specific hazard models, fitted by \code{\link{csc}}
#'   with a parametric \code{baseline}: each a distribution name or \code{c("rp", df)}, e.g.
#'   \code{list("weibull", "lognormal", c("rp", 4))}. They are cross-validated with the same folds as
#'   the other models, and their AIC (full likelihood of the data) is comparable with that of the joint
#'   and constrained direct models.
#' @param constrained optional list of constrained direct models (\code{\link{directcif}} with
#'   \code{estimation = "constrained"}) with options that \code{direct} cannot express: each a list
#'   with any of \code{baseline} (default \code{"gompertz"}), \code{link} (\code{"PH"}), \code{cond}
#'   (\code{"weibull"}, one value or one per competing cause), \code{timing} and \code{shape} (the
#'   covariates of the time distributions of the competing causes), e.g.
#'   \code{list(list(), list(timing = ~ age))}. Any number of causes. They are fitted by
#'   \code{directcif()} and cross-validated with the same folds; their AIC is comparable with that of
#'   the joint and constrained direct models.
#' @section Interpreting the result:
#' Judge the models on the cross-validated columns (ending in \code{_cv}): lower Brier score and
#' calibration error, higher IPA and AUC are better. Apparent measures, computed on the data used for
#' fitting, favour complex models. Differences in the third decimal of the Brier score are usually
#' within noise; prefer the simpler model when the measures are close. AIC can be compared only among
#' models with a full likelihood of the same data (joint and constrained direct models, parametric
#' cause-specific models); see Value.
#' @section Columns of the table:
#' Every measure is computed twice, from the predictions of the full-data fit (apparent) and from the
#' pooled out-of-fold predictions (columns ending in \code{_cv}): \code{brier}, \code{ipa}, \code{auc},
#' \code{calib_error} (see \code{\link{brier_score}}); \code{brier_se} and \code{auc_se}, standard
#' errors from influence functions with the censoring weights taken as known (as
#' \code{riskRegression::Score(conservative = TRUE)}), with the 95\% intervals \code{brier_lower},
#' \code{brier_upper}, \code{auc_lower}, \code{auc_upper}; \code{oe}, \code{oe_lower}, \code{oe_upper},
#' the observed (Aalen-Johansen) over the mean predicted risk with its interval; \code{ici},
#' \code{e50}, \code{e90}, the smooth calibration measures of Austin, Harrell and van Klaveren (2020);
#' \code{brier_diff}, \code{brier_diff_se}, \code{brier_diff_p} and \code{auc_diff},
#' \code{auc_diff_se}, \code{auc_diff_p}, the paired difference from the model with the lowest Brier
#' score at that time (named in \code{reference} and \code{reference_cv}), its standard error and the
#' p-value of the test of no difference. The paired tests use the same subjects for every model, so
#' they are much more sensitive than comparing two intervals. Cross-validated differences are an
#' honest guide; their p-values are approximate because the pooled predictions are not independent.
#'
#' The concordance index (\code{\link{cindex}}) is deliberately not reported: it is not a proper
#' measure for risks at a time horizon (Blanche, Kattan and Gerds 2019), so it can prefer a wrong
#' model. The AUC at the horizon is its proper counterpart.
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
#' Austin PC, Harrell FE, van Klaveren D (2020). Graphical calibration curves and the integrated
#'   calibration index (ICI) for survival models. Statistics in Medicine 39:2714-2742.
#' Blanche P, Kattan MW, Gerds TA (2019). The c-index is not proper for the evaluation of t-year
#'   predicted risks. Biostatistics 20:347-357.
#'
#' All references: \code{sx_references("compare_models")}.
#' @seealso
#' \code{\link{cv_compare}}, \code{\link{model_scan}}, \code{\link{brier_score}}, \code{\link{cindex}},
#' \code{\link{decision_curve}}, \code{\link{directcif}}, \code{\link{fg}}, \code{\link{csc}},
#' \code{vignette("prediction", "SurvX")}.
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
#'                      hazard = list("weibull"), constrained = list(list(timing = ~ x)))
#' cm
#' as.data.frame(cm)
#' @export
compare_models <- function(formula, data, times, direct = NULL, cause = 1L, folds = 5L, seed = 1L, df = 3L, hazard = NULL,
                           constrained = NULL) {
  seed <- .resolve_seed(seed)
  cl <- match.call()
  f <- .sx_frame(formula, data, what = "compare_models() (for models with a shared frailty use model_scan() or cv_compare())")
  n <- length(f$time); p <- ncol(f$X); K <- max(f$status); nt <- length(times)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  if (is.null(direct)) {
    direct <- list(c("gompertz", "PH"), c("gompertz", "PO"), c("gompertz", "PH", "weighted"), c("gompertz", "PO", "weighted"))
    if (K >= 2L) direct <- c(direct, list(c("gompertz", "PH", "constrained", "weibull")))
  }
  nd <- length(direct)
  part <- function(i, default) vapply(direct, function(d) if (length(d) >= i) as.character(d[[i]]) else default, character(1))
  base <- part(1L, "gompertz"); lk <- part(2L, "PH"); est <- part(3L, "joint"); cd <- part(4L, "logistic")
  if (any(!base %in% names(.base_codes)) || any(!lk %in% c("PH", "PO", "GOR")) || any(!est %in% c("joint", "weighted", "constrained")) ||
      any(!cd %in% names(.cond_codes)))
    stop("direct: each entry must be c(baseline, link), c(baseline, link, estimation) or c(baseline, link, \"constrained\", cond)", call. = FALSE)
  if (any(est == "constrained") && K < 2L) stop("constrained models need at least two causes", call. = FALSE)
  if (any(est == "constrained" & base %in% c("weibull", "spline"))) stop("constrained models need a baseline with a plateau (not weibull or spline)", call. = FALSE)
  # constrained models with K >= 3 causes, and those of `constrained`, are fitted by directcif() in R
  cspec <- .constrained_specs(constrained)
  if (K != 2L && any(est == "constrained")) {
    mv <- which(est == "constrained")
    cspec <- c(.constrained_specs(lapply(mv, function(j) list(baseline = base[j], link = lk[j], cond = cd[j]))), cspec)
    keep <- setdiff(seq_len(nd), mv)
    direct <- direct[keep]; base <- base[keep]; lk <- lk[keep]; est <- est[keep]; cd <- cd[keep]; nd <- length(direct)
  }
  if (length(cspec) && K < 2L) stop("constrained models need at least two causes", call. = FALSE)
  M <- nd + 2L
  r <- .with_seed(seed, .C("SX_compare3", as.integer(n), as.integer(p), as.integer(K), f$time, f$status, .rm(f$X), as.integer(cause),
          as.integer(nt), as.double(times), as.integer(folds), as.integer(nd),
          as.integer(unname(.base_codes[base])), as.integer(match(lk, c("PH", "PO", "GOR")) - 1L), rep(as.integer(df), nd),
          as.integer(match(est, c("joint", "weighted", "constrained")) - 1L), as.integer(unname(.cond_codes[cd])),
          metrics = double(M * nt * 8L), fit = double(M * 5L), risk = double(M * 2L * nt * n), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX"))
  .sx_rc(r$rc, "compare_models")
  suffix <- ifelse(est == "joint", "", ifelse(est == "weighted", ", weighted", paste0(", constrained (G ", cd, ")")))
  names_m <- c(if (nd) paste0("direct ", base, " ", lk, suffix), "Fine-Gray", "cause-specific Cox")
  # the predicted risks that the core scored, [subject, time, apparent/cross-validated, model]; all
  # measures are computed from them by .assess_table (with those of the models fitted in R below)
  rk <- array(r$risk, dim = c(n, nt, 2L, M))
  keep_risk <- function(x) if (all(is.finite(x))) x else NULL
  app <- lapply(seq_len(M), function(m) keep_risk(matrix(rk[, , 1L, m], n, nt)))
  cvr <- lapply(seq_len(M), function(m) if (folds >= 2L) keep_risk(matrix(rk[, , 2L, m], n, nt)) else NULL)
  ft <- .matrix_rm(r$fit, M, 5L)
  fit <- data.frame(model = names_m, loglik = ft[, 1L], npar = ft[, 2L], aic = ft[, 3L], flags = ft[, 4L], folds_not_converged = ft[, 5L])
  fit$aic[!is.finite(fit$aic)] <- NA_real_
  if (!is.null(hazard) || length(cspec)) {
    hz <- if (is.null(hazard)) list() else .hazard_specs(hazard)
    vars <- intersect(all.vars(formula), names(data))
    dd <- data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
    if (nrow(dd) != n) stop("compare_models: could not align the data for the parametric hazard models", call. = FALSE)
    fold <- rep(1L, n)
    if (folds >= 2L) {                      # the folds of the C++ comparison above
      rf <- .with_seed(seed, .C("SX_cv_folds", as.integer(n), as.integer(folds), out = integer(n), rc = integer(1), PACKAGE = "SurvX"))
      .sx_rc(rf$rc, "compare_models")
      fold <- rf$out + 1L
    }
    funs <- lapply(hz, function(h) { force(h); function(fm, d) csc(fm, d, baseline = h$baseline, df = h$df) })
    names(funs) <- vapply(hz, `[[`, "", "name")
    cf <- lapply(cspec, function(s) {
      force(s)
      function(fm, d) directcif(fm, d, baseline = s$baseline, link = s$link, estimation = "constrained", cause = cause,
                                cond = s$cond, timing = s$timing, shape = s$shape)
    })
    names(cf) <- vapply(cspec, `[[`, "", "name")
    funs <- c(funs, cf)
    res <- .cv_run(formula, dd, f, funs, times, cause, folds, fold)
    names_m <- c(names_m, names(funs)); app <- c(app, res$app); cvr <- c(cvr, res$cv)
    for (m in seq_along(funs)) {
      ft <- res$fits[[m]]
      ll <- if (is.null(ft)) NA_real_ else as.numeric(logLik(ft))
      np <- if (is.null(ft)) NA_real_ else attr(logLik(ft), "df")
      ok <- !is.null(ft) && (if (S7::S7_inherits(ft, DirectCIFFit)) isTRUE(ft@converged)
                             else all(vapply(ft@models, function(x) isTRUE(x@converged), logical(1))))
      fit <- rbind(fit, data.frame(model = names(funs)[m], loglik = ll, npar = np, aic = .sx_aic(ll, np),
                                   flags = if (ok) 0 else 1, folds_not_converged = res$failed[m]))
    }
  }
  tab <- .assess_table(names_m, f$time, f$status, times, cause, app, if (folds >= 2L) cvr else NULL)
  ModelComparison(model_call = cl, table = tab, fit = fit, folds = as.integer(folds), cause = as.integer(cause))
}

#' Cross-validated comparison of any set of models
#'
#' @description
#' Fits every model to the full data (apparent performance) and in K-fold cross-validation, and
#' reports the Brier score, the index of prediction accuracy, the time-dependent AUC and the
#' calibration error of the predicted risk of one cause at the chosen times, with standard errors,
#' paired differences from the best model, the observed/expected ratio and the integrated calibration
#' index, as \code{\link{compare_models}}.
#' \code{\link{compare_models}} does this in the C++ core for the models of this package;
#' \code{cv_compare} is for anything else that has a \code{predict(fit, newdata, times, cause)} method
#' (a model from another package wrapped in a function, a model with its own options).
#' @param formula \code{Event(time, status) ~ covariates}; passed to every model.
#' @param data a data frame; rows with missing values in the variables of \code{formula} are dropped.
#' @param models a named list of functions \code{function(formula, data)} that return a fitted model
#'   with a \code{predict(fit, newdata, times, cause)} method giving the cumulative incidence of
#'   the cause. Wrap any model, with its own options, in such
#'   a function. Fits of \code{cox()}, \code{parsurv()} and \code{cure()} for one event type are scored
#'   by \code{1 - predict(fit, newdata, times, type = "survival")}. With \code{+ frailty(id)} in
#'   \code{formula}, the folds keep every cluster whole and frailty models predict for a new cluster.
#' @param times time horizons.
#' @param cause cause of interest.
#' @param folds number of folds (less than 2: apparent performance only).
#' @param seed seed of the fold assignment: a whole number, or \code{NULL} to take the seed from the
#'   random numbers of the session (so that \code{set.seed()} before the call governs the result).
#'   All random draws come from R's generator (see \code{RNGkind()}); with a number the state of the
#'   session's generator is left as it was.
#' @param groups number of risk groups for the calibration error.
#' @section Interpreting the result:
#' As for \code{\link{compare_models}}: judge the models on the cross-validated measures (lower Brier
#' score and calibration error, higher IPA and AUC). The same folds are used for every model, so the
#' comparison is paired. The table has the same columns as that of \code{compare_models} (standard
#' errors, paired differences from the best model, O/E and ICI; see the section "Columns of the
#' table" there).
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
#' \code{\link{compare_models}}, \code{\link{brier_score}}, \code{\link{cindex}},
#' \code{\link{decision_curve}}, \code{vignette("prediction", "SurvX")}.
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
  seed <- .resolve_seed(seed)
  cl <- match.call()
  if (!is.list(models) || is.null(names(models)) || any(names(models) == "")) stop("models must be a named list of functions", call. = FALSE)
  vars <- intersect(all.vars(formula), names(data))
  data <- data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  f <- .sx_frame(formula, data, frailty_ok = TRUE)
  n <- length(f$time); nt <- length(times); M <- length(models)
  if (n != nrow(data)) stop("could not align the response with the data", call. = FALSE)
  fold <- rep(1L, n)
  if (folds >= 2L) {
    if (!is.null(f$frailty)) fold <- .scan_folds_cluster(f$frailty, folds, seed)   # whole clusters
    else fold <- .sx_folds(rep(0L, n), folds, seed)
  }
  res <- .cv_run(formula, data, f, models, times, cause, folds, fold)
  tab <- .assess_table(names(models), f$time, f$status, times, cause, res$app, if (folds >= 2L) res$cv else NULL, groups)
  fit <- NULL
  for (m in seq_len(M)) {
    ft <- res$fits[[m]]
    is_direct <- !is.null(ft) && S7::S7_inherits(ft, DirectCIFFit)
    is_par <- !is.null(ft) && (S7::S7_inherits(ft, ParCSCFit) || S7::S7_inherits(ft, ParSurvFit) ||
                               (S7::S7_inherits(ft, FrailtyFit) && ft@kind == "parametric"))
    ll <- if (is_direct) ft@loglik else if (is_par) as.numeric(logLik(ft)) else NA_real_
    np <- if (is_direct) length(ft@coefficients) else if (is_par) attr(logLik(ft), "df") else NA_integer_
    fit <- rbind(fit, data.frame(
      model = names(models)[m],
      loglik = ll,
      npar = np,
      aic = if ((is_direct && ft@estimation != "weighted") || is_par) .sx_aic(ll, np) else NA_real_,
      converged = if (is_direct) ft@converged else if (S7::S7_inherits(ft, ParCSCFit)) all(vapply(ft@models, function(x) isTRUE(x@converged), logical(1)))
                  else if (is_par) ft@converged else NA,
      folds_failed = res$failed[m]))
  }
  ModelComparison(model_call = cl, table = tab, fit = fit, folds = as.integer(folds), cause = as.integer(cause))
}

# predicted risk of the event (one event type) or cumulative incidence of the cause: the models of
# this package for one event type predict survival (frailty models: for a new cluster)
.cv_predict <- function(fit, nd, times, cause) {
  if (S7::S7_inherits(fit, CoxFit) || S7::S7_inherits(fit, ParSurvFit) || S7::S7_inherits(fit, CureFit) ||
      S7::S7_inherits(fit, FrailtyFit))
    return(.sx_ew("complement", as.matrix(stats::predict(fit, nd, times = times, type = "survival"))))
  stats::predict(fit, newdata = nd, times = times, cause = cause)
}

# Apparent and K-fold cross-validated predictions of a named list of fitting functions; fold: the
# fold of every row of data (1..folds). Returns the full-data fits, the number of failed folds and
# the predicted risks (n x nt) of the full-data fits (app) and pooled out of fold (cv; NULL when a
# fold failed or folds < 2), to be scored by .assess_table.
.cv_run <- function(formula, data, f, models, times, cause, folds, fold) {
  n <- length(f$time); nt <- length(times); M <- length(models)
  risk_of <- function(fun, train, test) {
    fit <- tryCatch(fun(formula, train), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    r <- tryCatch(.cv_predict(fit, test, times, cause), error = function(e) NULL)
    if (is.null(r) || !is.matrix(r) || nrow(r) != nrow(test) || ncol(r) != nt || any(!is.finite(r))) return(NULL)
    list(fit = fit, risk = r)
  }
  fits <- vector("list", M); failed <- integer(M); app <- vector("list", M); cv <- vector("list", M)
  for (m in seq_len(M)) {
    full <- risk_of(models[[m]], data, data)
    if (!is.null(full)) { fits[[m]] <- full$fit; app[[m]] <- unname(full$risk) }
    if (folds >= 2L) {
      pooled <- matrix(NA_real_, n, nt)
      for (k in seq_len(folds)) {
        te <- fold == k
        r <- risk_of(models[[m]], data[!te, , drop = FALSE], data[te, , drop = FALSE])
        if (is.null(r)) failed[m] <- failed[m] + 1L else pooled[te, ] <- r$risk
      }
      if (failed[m] == 0L) cv[[m]] <- pooled
    }
  }
  list(fits = fits, failed = failed, app = app, cv = cv)
}

# constrained = list(list(), list(timing = ~ age), list(cond = "lognormal", shape = ~ x)): one
# constrained direct model each, with a readable name
.constrained_specs <- function(constrained) {
  if (is.null(constrained)) return(list())
  if (!is.list(constrained) || (length(constrained) && !all(vapply(constrained, is.list, logical(1)))))
    stop("constrained must be a list of lists, such as list(list(), list(timing = ~ age))", call. = FALSE)
  ok <- c("baseline", "link", "cond", "timing", "shape")
  lapply(constrained, function(s) {
    bad <- setdiff(names(s), ok)
    if (length(bad)) stop("constrained: unknown element ", paste(bad, collapse = ", "), "; use ", paste(ok, collapse = ", "), call. = FALSE)
    out <- list(baseline = s$baseline %||% "gompertz", link = s$link %||% "PH", cond = s$cond %||% "weibull",
                timing = s$timing, shape = s$shape)
    lab <- function(x) if (is.null(x)) NULL
                       else if (inherits(x, "formula")) paste(deparse(x), collapse = "")
                       else paste0(names(x), ": ", vapply(x, function(v) paste(deparse(v), collapse = ""), ""), collapse = "; ")
    out$name <- paste0("direct ", out$baseline, " ", out$link, ", constrained (G ", paste(out$cond, collapse = "/"), ")",
                       if (!is.null(out$timing)) paste0(", timing ", lab(out$timing)),
                       if (!is.null(out$shape)) paste0(", shape ", lab(out$shape)))
    out
  })
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
