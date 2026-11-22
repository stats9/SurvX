# Further measures of prediction performance (C++ core: survcore/assess.hpp): standard errors and
# paired differences of the Brier score and AUC, observed/expected ratio, smooth calibration (ICI,
# E50, E90), concordance index, decision curves, and likelihood-ratio tests of nested models.

.assess_cols <- c("brier", "brier_null", "ipa", "auc", "calib_error", "brier_se", "auc_se", "oe", "oe_lower", "oe_upper",
                  "ici", "e50", "e90", "brier_diff", "brier_diff_se", "brier_diff_p", "auc_diff", "auc_diff_se", "auc_diff_p",
                  "brier_lower", "brier_upper", "auc_lower", "auc_upper")

# all measures of M models at the times; risks: list of n x nt matrices (NULL for a model without
# predictions). Returns a [measure, time, model] array and the reference model (lowest Brier score)
# of every time (1-based, NA if none).
.sx_assess <- function(time, status, risks, times, cause, groups = 10L, level = 0.95) {
  n <- length(time); nt <- length(times); M <- length(risks)
  rk <- unlist(lapply(risks, function(r) if (is.null(r)) rep(NA_real_, n * nt) else as.double(matrix(r, n, nt))), use.names = FALSE)
  r <- .C("SX_assess", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.integer(nt), as.double(times),
          as.integer(M), rk, as.integer(groups), as.double(.sx_qnorm(.sx_upper_prob(level))),
          out = double(length(.assess_cols) * nt * M), ref = integer(nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "prediction measures")
  a <- array(r$out, c(length(.assess_cols), nt, M))
  a[is.nan(a)] <- NA_real_
  dimnames(a) <- list(.assess_cols, NULL, NULL)
  list(a = a, ref = ifelse(r$ref >= 0L, r$ref + 1L, NA_integer_))
}

# the measures added by .sx_assess to the four of the first versions of the comparison tables
.assess_new <- c("brier_se", "brier_lower", "brier_upper", "auc_se", "auc_lower", "auc_upper", "oe", "oe_lower", "oe_upper",
                 "ici", "e50", "e90",
                 "brier_diff", "brier_diff_se", "brier_diff_p", "auc_diff", "auc_diff_se", "auc_diff_p")

# comparison table (one row per model and time) from the apparent and cross-validated predictions
.assess_table <- function(names_m, time, status, times, cause, app, cv, groups = 10L) {
  M <- length(names_m); nt <- length(times)
  A <- .sx_assess(time, status, app, times, cause, groups)
  C <- if (!is.null(cv)) .sx_assess(time, status, cv, times, cause, groups) else NULL
  na <- rep(NA_real_, nt)
  rows <- vector("list", M)
  for (m in seq_len(M)) {
    g <- function(X, nm) if (is.null(X)) na else X$a[nm, , m]
    cols <- list(model = rep(names_m[m], nt), time = times,
                 brier = g(A, "brier"), brier_cv = g(C, "brier"), ipa = g(A, "ipa"), ipa_cv = g(C, "ipa"),
                 auc = g(A, "auc"), auc_cv = g(C, "auc"), calib_error = g(A, "calib_error"), calib_error_cv = g(C, "calib_error"))
    for (nm in .assess_new) { cols[[nm]] <- g(A, nm); cols[[paste0(nm, "_cv")]] <- g(C, nm) }
    cols$reference <- ifelse(is.na(A$ref), NA_character_, names_m[A$ref])
    cols$reference_cv <- if (is.null(C)) rep(NA_character_, nt) else ifelse(is.na(C$ref), NA_character_, names_m[C$ref])
    rows[[m]] <- .df(cols)
  }
  do.call(rbind, rows)
}

#' Concordance index
#'
#' @description
#' Harrell's concordance index C of a risk score, Uno's version with inverse probability of
#' censoring weights, and, with competing risks, the version of Wolbers et al. (2009) for one
#' cause, with the standard error of the infinitesimal jackknife.
#'
#' \strong{Use it} to describe how well a score ranks subjects by the time to the event over the
#' whole follow-up. \strong{Do not use it} to compare models that predict the risk by a fixed time
#' \eqn{t}: the C index is not a proper measure for \eqn{t}-year risks (Blanche, Kattan and Gerds
#' 2019), so a wrong model can have a higher C. For that use the AUC at \eqn{t} and the Brier score
#' (\code{\link{auc_t}}, \code{\link{brier_score}}, \code{\link{compare_models}}), which is why the
#' comparison tables of the package do not report C.
#' @param object a fitted model of the package, or a numeric risk score (higher = earlier event),
#'   one per subject of \code{time} and \code{status}.
#' @param time,status follow-up time and status (0 = censored, 1..K = cause) when \code{object} is a
#'   score; for a fitted model the data it was fitted to (or \code{newdata}).
#' @param newdata for a fitted model: data on which to compute C (default: the data of the fit;
#'   must contain the outcome variables).
#' @param cause cause of interest.
#' @param method \code{"harrell"} (all comparable pairs have weight 1) or \code{"uno"} (a pair whose
#'   first event is at \eqn{T_i} has weight \eqn{1/\hat G(T_i-)^2}, \eqn{\hat G} the Kaplan-Meier
#'   estimate of the censoring distribution).
#' @param tau only events at times \eqn{\le} \code{tau} define comparable pairs (Uno et al. 2011
#'   recommend a \code{tau} where follow-up is still adequate).
#' @param times for a model that does not have a linear predictor (Fine-Gray, cause-specific and
#'   direct models, cure models), the risk score is the predicted risk at this one time; default:
#'   the median follow-up time.
#' @param level confidence level.
#' @details
#' A pair \eqn{(i, j)} is comparable when \eqn{i} has the event of interest at \eqn{T_i \le \tau}
#' and \eqn{j} is known to be event-free at \eqn{T_i}: \eqn{T_j > T_i}, or \eqn{T_j = T_i} with
#' \eqn{j} censored, or (competing risks) \eqn{j} failed from another cause at any time. Events of
#' interest at the same time are not compared. The pair is concordant when \eqn{i} has the higher
#' score; tied scores count 1/2. \eqn{C} is the weighted share of concordant pairs.
#'
#' The standard error is \eqn{(\sum_u d_u^2)^{1/2}} with \eqn{d_u = (N_u - C D_u)/D}, where
#' \eqn{N_u} and \eqn{D_u} are the weighted concordance and the weight of the comparable pairs that
#' contain subject \eqn{u} and \eqn{D} the total weight: the infinitesimal jackknife of
#' \code{survival::concordance}. With Uno's weights the estimation of \eqn{\hat G} is ignored, as
#' there. The computation sorts the subjects once and uses binary indexed trees, so it takes
#' \eqn{O(n \log n)} time.
#' @section Interpreting the result:
#' 0.5 is no better than chance, 1 is perfect ranking; 0.6 to 0.75 is common for clinical scores.
#' Harrell's C depends on the censoring distribution when the score is not perfectly
#' proportional; Uno's version removes that dependence up to \code{tau}.
#' @section Validation:
#' Harrell's C and its standard error equal \code{survival::concordance(..., reverse = TRUE)} and,
#' with \code{method = "uno"}, \code{timewt = "n/G2"} with \code{ymax = tau} (validation suite,
#' \code{validation/R/07_assess.R}); with competing risks the estimate is reported next to \code{pec::cindex}. The C++ code
#' is tested against a direct computation over all pairs.
#' @return an object of class \code{SurvXConcordance}: a list with \code{concordance}, \code{se},
#'   \code{lower}, \code{upper}, \code{counts} (weighted numbers of concordant, discordant and tied
#'   pairs and the number of comparable pairs), \code{dfbeta}, \code{method}, \code{tau}, \code{cause}.
#' @references
#' Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA (1982) Evaluating the yield of medical tests.
#'   JAMA 247:2543-2546.
#'
#' Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011) On the C-statistics for evaluating overall
#'   adequacy of risk prediction procedures with censored survival data. Statistics in Medicine
#'   30:1105-1117.
#'
#' Wolbers M, Koller MT, Witteman JCM, Steyerberg EW (2009) Prognostic models with competing risks:
#'   methods and application to coronary risk prediction. Epidemiology 20:555-561.
#'
#' Blanche P, Kattan MW, Gerds TA (2019) The c-index is not proper for the evaluation of t-year
#'   predicted risks. Biostatistics 20:347-357.
#'
#' Therneau T, Atkinson E (2024) Concordance. Vignette of the survival package.
#' @seealso \code{\link{auc_t}}, \code{\link{compare_models}}.
#' @examples
#' set.seed(3)
#' n <- 300
#' d <- data.frame(x = rnorm(n))
#' t1 <- rexp(n, 0.3 * exp(0.8 * d$x)); t2 <- rexp(n, 0.2); cens <- runif(n, 0, 5)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#' cindex(d$x, d$time, d$status)                  # Wolbers' C for cause 1
#' cindex(d$x, d$time, d$status, method = "uno", tau = 3)
#' cindex(fg(Event(time, status) ~ x, d), times = 2)
#' @export
cindex <- function(object, time = NULL, status = NULL, newdata = NULL, cause = 1L, method = c("harrell", "uno"),
                   tau = Inf, times = NULL, level = 0.95) {
  method <- match.arg(method)
  if (is.numeric(object) && is.null(dim(object))) {
    score <- as.double(object)
    if (is.null(time) || is.null(status)) stop("cindex: give time and status with a numeric score", call. = FALSE)
  } else {
    sc <- .cindex_score(object, newdata, cause, times)
    score <- sc$score; time <- sc$time; status <- sc$status
  }
  time <- as.double(time); status <- as.integer(status); n <- length(time)
  if (length(score) != n || length(status) != n) stop("cindex: score, time and status must have the same length", call. = FALSE)
  ok <- is.finite(score) & is.finite(time) & !is.na(status)
  if (!all(ok)) { score <- score[ok]; time <- time[ok]; status <- status[ok]; n <- length(time) }
  r <- .C("SX_concordance", as.integer(n), time, status, as.integer(cause), score, as.double(tau), as.integer(method == "uno"),
          as.double(.sx_qnorm(.sx_upper_prob(level))), out = double(8), dfbeta = double(n), rc = integer(1), NAOK = TRUE,
          PACKAGE = "SurvX")
  .sx_rc(r$rc, "cindex")
  o <- r$out; o[is.nan(o)] <- NA_real_
  structure(list(concordance = o[1L], se = o[2L], lower = o[7L], upper = o[8L],
                 counts = c(concordant = o[3L], discordant = o[4L], tied = o[5L], pairs = o[6L]),
                 dfbeta = r$dfbeta, method = method, tau = tau, cause = as.integer(cause), level = level,
                 competing = any(status > 0L & status != cause), n = n),
            class = "SurvXConcordance")
}

# risk score of a fitted model: the linear predictor where the model has one for the cause (higher
# = earlier event), else the predicted risk of the cause at one time
.cindex_score <- function(object, newdata, cause, times) {
  nd <- newdata
  if (is.null(nd)) { time <- object@time; status <- object@status }
  else {
    vars <- intersect(all.vars(stats::formula(object)), names(nd))
    nd <- nd[stats::complete.cases(nd[, vars, drop = FALSE]), , drop = FALSE]
    f <- .sx_frame(stats::formula(object), nd, entry_ok = TRUE, specials_ok = TRUE, frailty_ok = TRUE)
    time <- f$time; status <- f$status
  }
  if (S7::S7_inherits(object, CoxFit) && is.null(object@tt))
    return(list(score = as.double(stats::predict(object, newdata = nd, type = "lp")), time = time, status = status))
  if (is.null(times)) times <- stats::median(time)
  if (length(times) != 1L) stop("cindex: times must be one time", call. = FALSE)
  d <- if (is.null(nd)) .cindex_data(object) else nd
  r <- .cv_predict(object, d, times, cause)
  list(score = as.double(r[, 1L]), time = time, status = status)
}

# the data of a fit (for predictions at the fitted subjects)
.cindex_data <- function(object) {
  cl <- object@model_call
  d <- tryCatch(eval(cl$data, envir = environment(stats::formula(object))), error = function(e) NULL)
  if (is.null(d)) stop("cindex: give newdata (the data of the fit could not be found)", call. = FALSE)
  vars <- intersect(all.vars(stats::formula(object)), names(d))
  d[stats::complete.cases(d[, vars, drop = FALSE]), , drop = FALSE]
}

#' @export
print.SurvXConcordance <- function(x, digits = 4L, ...) {
  lab <- if (x$method == "uno") "Uno's concordance index" else "Harrell's concordance index"
  if (x$competing) lab <- paste0(lab, " for cause ", x$cause, " with competing risks (Wolbers et al. 2009)")
  cat(lab, if (is.finite(x$tau)) paste0(", events up to time ", format(x$tau)), "\n", sep = "")
  cat("  C = ", format(round(x$concordance, digits), nsmall = digits), " (se ", format(round(x$se, digits), nsmall = digits), "), ",
      format(100 * x$level), "% interval ", format(round(x$lower, digits), nsmall = digits), " to ",
      format(round(x$upper, digits), nsmall = digits), "\n", sep = "")
  cat("  ", format(x$counts[["pairs"]], big.mark = ","), " comparable pairs; n = ", x$n, "\n", sep = "")
  cat("  C ranks over the whole follow-up; to compare predicted risks at a time use auc_t() and brier_score().\n")
  invisible(x)
}

#' Decision curve analysis
#'
#' @description
#' Net benefit of deciding by predicted risks at a time horizon, against treating everybody or
#' nobody, over a range of risk thresholds (Vickers and Elkin 2006), with right censoring and
#' competing risks (Vickers et al. 2008).
#'
#' \strong{Use it} to judge whether using a model to decide (treat, screen, refer) would do more
#' good than harm, at the thresholds that matter for the decision.
#' @param risks a named list of predicted risks of the cause by \code{time} (numeric vectors, one per
#'   subject), or of fitted models of the package (their predictions for \code{data} at \code{time}).
#'   A single numeric vector or model is allowed.
#' @param time,status follow-up time and status (0 = censored, 1..K = cause); with fitted models and
#'   \code{data} they may be omitted (taken from the formula of the first model).
#' @param data data for the predictions of fitted models.
#' @param t the time horizon.
#' @param cause cause of interest.
#' @param thresholds risk thresholds \eqn{p} in \eqn{[0, 1)}.
#' @details
#' Treating the subjects whose predicted risk is at least \eqn{p} has net benefit
#' \eqn{NB(p) = \pi_p F_p(t) - \pi_p \{1 - F_p(t)\} p/(1-p)}, with \eqn{\pi_p} the share of subjects
#' above the threshold and \eqn{F_p(t)} the Kaplan-Meier (Aalen-Johansen) risk of the cause by
#' \eqn{t} among them: true positives minus false positives per subject, the false ones weighted by
#' the odds of the threshold. Treating everybody has \eqn{F(t) - \{1 - F(t)\} p/(1-p)}, treating
#' nobody 0.
#' @section Interpreting the result:
#' A model is useful at threshold \eqn{p} when its net benefit is above both "treat all" and
#' "treat none". A net benefit of 0.05 means the same as finding 5 more true cases per 100 subjects
#' with no false positives. Look only at the thresholds a decision maker would use.
#' @section Validation:
#' Equals the definition computed with \code{survival::survfit} and \code{dcurves::dca} with a
#' \code{Surv} outcome and the same thresholds (validation suite, \code{validation/R/07_assess.R}).
#' @return an object of class \code{DecisionCurve}: a list with \code{table} (columns
#'   \code{model}, \code{threshold}, \code{net_benefit}, including the rows of "treat all" and
#'   "treat none"), \code{t}, \code{cause}. Methods \code{print}, \code{plot} and
#'   \code{as.data.frame}.
#' @references
#' Vickers AJ, Elkin EB (2006) Decision curve analysis: a novel method for evaluating prediction
#'   models. Medical Decision Making 26:565-574.
#'
#' Vickers AJ, Cronin AM, Elkin EB, Gonen M (2008) Extensions to decision curve analysis, a novel
#'   method for evaluating diagnostic tests, prediction models and molecular markers. BMC Medical
#'   Informatics and Decision Making 8:53.
#' @seealso \code{\link{brier_score}}, \code{\link{compare_models}}.
#' @examples
#' set.seed(3)
#' n <- 400
#' d <- data.frame(x = rnorm(n))
#' t1 <- rexp(n, 0.3 * exp(0.8 * d$x)); t2 <- rexp(n, 0.2); cens <- runif(n, 0, 5)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#' fit <- fg(Event(time, status) ~ x, d)
#' dc <- decision_curve(list("Fine-Gray" = fit), data = d, t = 2, thresholds = seq(0, 0.6, 0.02))
#' dc
#' plot(dc)
#' @export
decision_curve <- function(risks, time = NULL, status = NULL, data = NULL, t, cause = 1L,
                           thresholds = seq(0, 0.5, by = 0.01)) {
  if (!is.list(risks) || S7::S7_inherits(risks)) risks <- list(model = risks)
  if (is.null(names(risks)) || any(names(risks) == "")) names(risks) <- paste0("model ", seq_along(risks))
  if (length(t) != 1L || !is.finite(t)) stop("decision_curve: t must be one finite time", call. = FALSE)
  thresholds <- as.double(thresholds)
  if (any(!is.finite(thresholds) | thresholds < 0 | thresholds >= 1)) stop("decision_curve: thresholds must be in [0, 1)", call. = FALSE)
  if (is.null(time) || is.null(status)) {
    first <- Filter(function(r) !is.numeric(r), risks)
    if (!length(first) || is.null(data)) stop("decision_curve: give time and status, or fitted models with data", call. = FALSE)
    f <- .sx_frame(stats::formula(first[[1L]]), data, frailty_ok = TRUE)
    time <- f$time; status <- f$status
    vars <- intersect(all.vars(stats::formula(first[[1L]])), names(data))
    data <- data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  }
  time <- as.double(time); status <- as.integer(status); n <- length(time); m <- length(thresholds)
  tab <- NULL; nb_all <- NULL
  for (k in seq_along(risks)) {
    r <- risks[[k]]
    pr <- if (is.numeric(r)) as.double(r) else as.double(.cv_predict(r, data, t, cause)[, 1L])
    if (length(pr) != n || any(!is.finite(pr))) stop("decision_curve: the risks of ", names(risks)[k], " must be finite, one per subject", call. = FALSE)
    o <- .C("SX_net_benefit", as.integer(n), time, status, as.integer(cause), pr, as.double(t), as.integer(m), thresholds,
            nb = double(m), all = double(m), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(o$rc, "decision_curve")
    tab <- rbind(tab, .df(list(model = rep(names(risks)[k], m), threshold = thresholds, net_benefit = o$nb)))
    nb_all <- o$all
  }
  tab <- rbind(tab, .df(list(model = rep("treat all", m), threshold = thresholds, net_benefit = nb_all)),
               .df(list(model = rep("treat none", m), threshold = thresholds, net_benefit = rep(0, m))))
  structure(list(table = tab, t = t, cause = as.integer(cause), n = n), class = "DecisionCurve")
}

#' @export
print.DecisionCurve <- function(x, ...) {
  tb <- x$table
  mods <- unique(tb$model)
  th <- unique(tb$threshold)
  show <- th[unique(round(seq(1, length(th), length.out = min(length(th), 8L))))]
  w <- sapply(mods, function(mm) { s <- tb[tb$model == mm, , drop = FALSE]; s$net_benefit[match(show, s$threshold)] })
  w <- matrix(w, nrow = length(show), dimnames = list(NULL, mods))
  df <- data.frame(threshold = format(show), stringsAsFactors = FALSE)
  for (j in seq_along(mods)) df[[mods[j]]] <- .num(w[, j], 4L)
  cat("Decision curve: net benefit at time ", format(x$t), " (cause ", x$cause, "; n = ", x$n, ")\n", sep = "")
  cat(.table_lines(df), sep = "\n")
  cat("A model helps at a threshold where its net benefit is above both 'treat all' and 'treat none'.\n")
  invisible(x)
}

#' @export
as.data.frame.DecisionCurve <- function(x, ...) x$table

#' @export
plot.DecisionCurve <- function(x, ylim = NULL, col = NULL, lwd = 2, xlab = "risk threshold", ylab = "net benefit",
                               main = NULL, legend = TRUE, ...) {
  tb <- x$table
  mods <- unique(tb$model); mm <- setdiff(mods, c("treat all", "treat none"))
  if (is.null(col)) col <- rep_len(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"), length(mm))
  nb <- tb$net_benefit[tb$model %in% c(mm, "treat none")]
  if (is.null(ylim)) {
    top <- max(c(nb, tb$net_benefit[tb$model == "treat all"]), na.rm = TRUE)
    ylim <- c(-0.05 * max(top, 0.01), max(top, 0.01) * 1.05)
  }
  graphics::plot(NA, xlim = range(tb$threshold), ylim = ylim, xlab = xlab, ylab = ylab,
                 main = if (is.null(main)) paste0("Decision curve at time ", format(x$t)) else main, ...)
  s <- tb[tb$model == "treat all", ]; graphics::lines(s$threshold, s$net_benefit, col = "grey45", lwd = 1, lty = 2)
  graphics::abline(h = 0, col = "grey20", lwd = 1)
  for (k in seq_along(mm)) { s <- tb[tb$model == mm[k], ]; graphics::lines(s$threshold, s$net_benefit, col = col[k], lwd = lwd) }
  if (legend) graphics::legend("topright", legend = c(mm, "treat all", "treat none"), col = c(col, "grey45", "grey20"),
                               lty = c(rep(1, length(mm)), 2, 1), lwd = c(rep(lwd, length(mm)), 1, 1), bty = "n", cex = 0.85)
  invisible(x)
}

# ------------------------------------------------------------------------------------------------
# likelihood-ratio tests of nested models
# ------------------------------------------------------------------------------------------------
.anova_ok <- function(fit) {
  if (S7::S7_inherits(fit, FineGrayFit))
    stop("anova: the Fine-Gray model maximises a weighted pseudo-likelihood, so likelihood-ratio tests are not valid; ",
         "use the Wald tests of summary()", call. = FALSE)
  if (S7::S7_inherits(fit, DirectCIFFit) && identical(fit@estimation, "weighted"))
    stop("anova: weighted direct models maximise a pseudo-likelihood; use the Wald tests of summary()", call. = FALSE)
  if (S7::S7_inherits(fit, CoxFit) && isTRUE(fit@robust))
    warning("anova: the likelihood-ratio test ignores the clustering of the robust variance", call. = FALSE)
  invisible(TRUE)
}

.anova_terms <- function(fit) attr(stats::terms(stats::formula(fit)), "term.labels")

.sx_anova <- function(object, ...) {
  fits <- c(list(object), list(...))
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (length(fits) < 2L) stop("anova: give two or more nested models fitted to the same data, smallest first", call. = FALSE)
  for (f in fits) .anova_ok(f)
  cls <- vapply(fits, function(f) class(f)[1L], "")
  if (length(unique(cls)) > 1L) stop("anova: all models must be of the same kind (", paste(unique(cls), collapse = ", "), ")", call. = FALSE)
  ll <- lapply(fits, logLik)
  lv <- vapply(ll, as.numeric, numeric(1)); df <- vapply(ll, function(l) as.numeric(attr(l, "df")), numeric(1))
  nn <- vapply(ll, function(l) as.numeric(attr(l, "nobs")), numeric(1))
  if (length(unique(nn)) > 1L) stop("anova: the models were fitted to different numbers of subjects (", paste(nn, collapse = ", "), ")", call. = FALSE)
  for (k in seq_along(fits)[-1L]) {
    a <- .anova_terms(fits[[k - 1L]]); b <- .anova_terms(fits[[k]])
    if (!all(a %in% b)) stop("anova: model ", k - 1L, " is not nested in model ", k, " (terms ", paste(setdiff(a, b), collapse = ", "),
                             " are missing from the larger model); order the models from the smallest", call. = FALSE)
    if (!(df[k] > df[k - 1L])) stop("anova: model ", k, " must have more parameters than model ", k - 1L, call. = FALSE)
    if (S7::S7_inherits(fits[[k]], SurvXFit) && length(fits[[k]]@time) == length(fits[[k - 1L]]@time) &&
        (!identical(fits[[k]]@time, fits[[k - 1L]]@time) || !identical(fits[[k]]@status, fits[[k - 1L]]@status)))
      stop("anova: models ", k - 1L, " and ", k, " were not fitted to the same outcome", call. = FALSE)
  }
  K <- length(fits)
  dd <- c(NA_real_, .sx_ew("sub", df[-1L], df[-K]))
  lr <- .sx_lr(lv[-1L], lv[-K], dd[-1L])                        # 2 (ll_k - ll_{k-1}) and its chi-square p-value
  chi <- c(NA_real_, lr$stat); p <- c(NA_real_, lr$p)
  forms <- vapply(fits, function(f) paste(deparse(stats::formula(f)[[3L]], width.cutoff = 500L), collapse = ""), "")
  tab <- data.frame(loglik = lv, Df = dd, Chisq = chi, `Pr(>|Chi|)` = p, check.names = FALSE,
                    row.names = paste("Model", seq_along(fits)))
  structure(tab, heading = c("Likelihood-ratio tests of nested models\n", paste0(" Model ", seq_along(fits), ": ~ ", forms, collapse = "\n"), "\n"),
            class = c("anova", "data.frame"))
}

# likelihood-ratio tests of nested models: documented in ?SurvX-methods (man/SurvX-methods.Rd)
S7::method(anova, SurvXFit) <- function(object, ...) .sx_anova(object, ...)
S7::method(anova, CSCFit) <- function(object, ...) .sx_anova(object, ...)
S7::method(anova, ParCSCFit) <- function(object, ...) .sx_anova(object, ...)

# one line of the summary of a fitted model: Harrell's C of its risk score on the data of the fit
# (nothing when it cannot be computed)
.cindex_line <- function(score, time, status, what = "") {
  r <- tryCatch(cindex(as.double(score), time, status, cause = 1L), error = function(e) NULL)
  if (is.null(r) || !is.finite(r$concordance)) return(character(0))
  sprintf("     concordance (Harrell's C%s): %.3f (se %.3f)", what, r$concordance, r$se)
}
