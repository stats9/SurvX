# Parametric survival regression, parametric cause-specific hazards and the non-parametric hazard.
# The likelihoods, fits, predictions and numerical integration are in the C++ core
# (src/survcore/parametric.hpp, src/survcore/hazard.hpp); this file only prepares data and output.

.par_dists <- c(exponential = 0L, weibull = 1L, gompertz = 2L, lognormal = 3L, loglogistic = 4L)
.par_label <- c(exponential = "exponential", weibull = "Weibull", gompertz = "Gompertz",
                lognormal = "log-normal", loglogistic = "log-logistic")
.par_types <- c(survival = 0L, hazard = 1L, cumhaz = 2L, quantile = 3L, rmst = 4L)

# the covariate scale of a distribution: PH or AFT, checked against what the distribution supports
.par_scale <- function(dist, scale) {
  allowed <- switch(dist, exponential = c("PH", "AFT"), weibull = c("PH", "AFT"), gompertz = "PH",
                    lognormal = "AFT", loglogistic = "AFT")
  if (is.null(scale) || is.na(scale)) return(allowed[1L])
  sc <- toupper(scale)
  if (!sc %in% c("PH", "AFT")) stop("scale must be \"PH\" or \"AFT\"", call. = FALSE)
  if (!sc %in% allowed)
    stop("the ", .par_label[[dist]], " distribution is available on the ", paste(allowed, collapse = " and "),
         " scale only (", if (sc == "PH") "proportional hazards" else "accelerated failure time", " was requested)", call. = FALSE)
  sc
}

.par_dist <- function(baseline) {
  b <- tolower(baseline)
  if (!b %in% names(.par_dists))
    stop("baseline must be one of ", paste0("\"", names(.par_dists), "\"", collapse = ", "), call. = FALSE)
  b
}

.shape_name <- function(dist, scale) {
  if (dist == "exponential") character(0) else if (scale == "AFT") "log(scale)" else if (dist == "gompertz") "shape" else "log(shape)"
}

# one fit in the C++ core; status01 is 0/1
.par_fit <- function(time, status01, entry, X, dist, scale, init = NULL, max_iter = 200L, tol = 1e-10) {
  n <- length(time); p <- ncol(X)
  P <- 1L + p + as.integer(dist != "exponential")
  r <- .C("SX_par_fit", as.integer(n), as.integer(p), .par_dists[[dist]], as.integer(scale == "AFT"),
          as.integer(!is.null(entry)), if (is.null(entry)) double(1) else as.double(entry),
          as.double(time), as.integer(status01), .rm(X), as.integer(!is.null(init)),
          if (is.null(init)) double(P) else as.double(init), as.integer(max_iter), as.double(tol),
          theta = double(P), var = double(P * P), info_d = double(3), info_i = integer(4), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "parametric fit")
  nm <- c("(Intercept)", colnames(X), .shape_name(dist, scale))
  list(theta = stats::setNames(r$theta, nm), var = matrix(r$var, P, P, byrow = TRUE, dimnames = list(nm, nm)),
       loglik = r$info_d[2L], iterations = r$info_i[1L], converged = as.logical(r$info_i[2L]),
       flags = r$info_i[3L], nevent = r$info_i[4L])
}

.parsurv_core <- function(cl, f, status01, dist, scale, init, max_iter, tol, label) {
  if (any(f$time <= 0)) stop("times must be positive for parametric models", call. = FALSE)
  r <- .par_fit(f$time, status01, f$entry, f$X, dist, scale, init, max_iter, tol)
  ln <- if (ncol(f$X)) tryCatch(.par_fit(f$time, status01, f$entry, f$X[, 0L, drop = FALSE], dist, scale, NULL, max_iter, tol)$loglik,
                                error = function(e) NA_real_) else r$loglik
  ParSurvFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
             coefficients = r$theta, vcov = r$var, loglik = r$loglik, n = length(f$time),
             nevent = as.integer(r$nevent), iterations = as.integer(r$iterations), converged = isTRUE(r$converged),
             flags = as.integer(r$flags), time = f$time, status = as.integer(status01), X = f$X,
             dist = dist, scale = scale, dist_code = .par_dists[[dist]], aft = as.integer(scale == "AFT"),
             entry = f$entry, loglik_null = ln, event_label = label)
}

#' Parametric survival regression
#'
#' Maximum likelihood regression for a single event time with an exponential, Weibull, Gompertz,
#' log-normal or log-logistic distribution, under right censoring and, optionally, delayed entry
#' (left truncation). Covariates act on the hazard (proportional hazards, \code{scale = "PH"}) or
#' on the time scale (accelerated failure time, \code{scale = "AFT"}).
#'
#' @details
#' Parametrisations, with \eqn{\eta = \beta_0 + z'\beta}:
#' \describe{
#'   \item{PH}{\eqn{H(t \mid z) = e^{\eta} H_0(t)} with \eqn{H_0(t) = t} (exponential),
#'     \eqn{t^{\kappa}} (Weibull, \code{log(shape)} \eqn{= \log\kappa}) or
#'     \eqn{(e^{a t} - 1)/a} (Gompertz, \code{shape} \eqn{= a}, any sign; \eqn{a < 0} gives a
#'     fraction that never fails). \eqn{e^{\beta_j}} is a hazard ratio. These are the
#'     \code{"exp"}, \code{"weibullPH"} and \code{"gompertz"} models of \pkg{flexsurv}.}
#'   \item{AFT}{\eqn{\log T = \eta + \sigma W}, \code{log(scale)} \eqn{= \log\sigma}, with
#'     \eqn{W} standard minimum extreme value (Weibull; exponential with \eqn{\sigma = 1}), standard
#'     normal (log-normal) or standard logistic (log-logistic). \eqn{e^{\beta_j}} is a time ratio
#'     (acceleration factor). These are the models of \code{survival::survreg}.}
#' }
#' The log-likelihood is \eqn{\sum_i d_i \log h(t_i) + \log S(t_i) - \log S(e_i)} with entry time
#' \eqn{e_i} (zero without delayed entry); gradient and Hessian are analytic and the maximum is
#' found by Newton-Raphson with step halving, starting from the exponential fit. Standard errors are
#' from the observed information.
#'
#' @param formula \code{Event(time, status) ~ covariates}, status 0/1; \code{Event(time, status,
#'   entry = )} for delayed entry.
#' @param data a data frame.
#' @param baseline \code{"weibull"} (default), \code{"exponential"}, \code{"gompertz"},
#'   \code{"lognormal"} or \code{"loglogistic"}.
#' @param scale \code{"PH"} or \code{"AFT"}; default \code{"PH"} where it exists (exponential,
#'   Weibull, Gompertz), else \code{"AFT"} (log-normal, log-logistic).
#' @param init optional starting values (intercept, covariates, shape).
#' @param max_iter,tol iteration limit and convergence tolerance (relative change of the
#'   log-likelihood).
#' @return an S7 object of class \code{ParSurvFit} with the usual methods (\code{print},
#'   \code{summary}, \code{estimates}, \code{coef}, \code{vcov}, \code{confint}, \code{logLik},
#'   \code{AIC}, \code{BIC}, \code{update}) and \code{predict(fit, newdata, times, type)} with
#'   \code{type} \code{"survival"}, \code{"hazard"}, \code{"cumhaz"}, \code{"quantile"} (at
#'   probabilities \code{p}) or \code{"rmst"} (restricted mean survival time up to \code{times}).
#' @references
#' Kalbfleisch JD, Prentice RL (2002) The Statistical Analysis of Failure Time Data, 2nd ed. Wiley,
#'   Chapters 2-3.
#' Lawless JF (2003) Statistical Models and Methods for Lifetime Data, 2nd ed. Wiley, Chapters 4-6.
#' Klein JP, Moeschberger ML (2003) Survival Analysis: Techniques for Censored and Truncated Data,
#'   2nd ed. Springer, Sections 3.5 and 12.
#' @seealso \code{\link{csc}} with a parametric \code{baseline} for competing risks,
#'   \code{\link{hazard_plot}}, \code{\link{hazard_np}}.
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- rbinom(n, 1, 0.5)
#' t <- rweibull(n, shape = 1.5, scale = 4 * exp(-0.5 * x)); cens <- runif(n, 0, 8)
#' d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x)
#' fit <- parsurv(Event(time, status) ~ x, d)
#' fit
#' summary(fit)
#' predict(fit, data.frame(x = 0:1), times = c(1, 2, 4), se = TRUE)$est
#' predict(fit, data.frame(x = 0:1), type = "quantile", p = 0.5)
#' AIC(fit, update(fit, baseline = "lognormal"))
#' @export
parsurv <- function(formula, data, baseline = "weibull", scale = NULL, init = NULL, max_iter = 200L, tol = 1e-10) {
  cl <- match.call()
  dist <- .par_dist(baseline)
  scale <- .par_scale(dist, scale)
  f <- .sx_frame(formula, data, entry_ok = TRUE)
  if (any(f$status > 1L))
    stop("parsurv() is for one event (status 0/1); for competing risks use csc(..., baseline = \"", dist, "\")", call. = FALSE)
  .parsurv_core(cl, f, f$status, dist, scale, init, max_iter, tol, "the event")
}

.pcsc <- function(cl, f, baseline, scale, max_iter, tol) {
  K <- max(f$status)
  if (K < 1L) stop("no events", call. = FALSE)
  dist <- rep_len(vapply(baseline, .par_dist, character(1)), K)
  sc <- rep_len(if (is.null(scale)) NA_character_ else scale, K)
  sc <- vapply(seq_len(K), function(k) .par_scale(dist[k], sc[k]), character(1))
  models <- lapply(seq_len(K), function(k)
    tryCatch(.parsurv_core(cl, f, as.integer(f$status == k), dist[k], sc[k], NULL, max_iter, tol, paste("cause", k)),
             error = function(e) stop("cause ", k, ": ", conditionMessage(e), call. = FALSE)))
  names(models) <- paste0("cause", seq_len(K))
  ParCSCFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts, models = models,
            K = as.integer(K), dist = unname(dist), scale = unname(sc), n = length(f$time), time = f$time,
            status = f$status, entry = f$entry, X = f$X)
}

# ----------------------------------------------------------------------------------------------
# prediction
# ----------------------------------------------------------------------------------------------

# pointwise limits: survival-type quantities on the complementary log-log scale, positive ones on the log scale
.pred_limits <- function(est, se, type, level) {
  zq <- .zq(level)
  if (type %in% c("survival", "cif")) {
    lim <- .cif_interval(est, se, zq)
    return(list(lower = array(lim$lower, dim(est), dimnames(est)), upper = array(lim$upper, dim(est), dimnames(est))))
  }
  e <- suppressWarnings(zq * se / est)
  lo <- suppressWarnings(est * exp(-e)); up <- suppressWarnings(est * exp(e))
  bad <- !is.finite(est) | !is.finite(se) | !(est > 0)
  lo[bad] <- NA_real_; up[bad] <- NA_real_
  list(lower = lo, upper = up)
}

# predictions of a parametric fit at design rows Xn (no intercept column)
.par_predict_X <- function(object, Xn, at, type, se, level, average = FALSE) {
  nn <- nrow(Xn); nt <- length(at)
  usevar <- isTRUE(se) && all(is.finite(object@vcov))
  cells <- if (average) nt else nn * nt
  r <- .C("SX_par_predict", object@dist_code, object@aft, as.integer(ncol(Xn)), as.double(object@coefficients),
          as.integer(usevar), if (usevar) as.double(object@vcov) else double(1), as.integer(nn), .rm(Xn),
          as.integer(nt), as.double(at), .par_types[[type]], as.integer(average),
          out = double(cells), se = double(cells), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "predict")
  cn <- if (type == "quantile") paste0("p=", at) else paste0("t=", at)
  nr <- if (average) 1L else nn
  est <- .matrix_rm(r$out, nr, nt, list(NULL, cn))
  if (!isTRUE(se)) return(est)
  s <- if (usevar) .matrix_rm(r$se, nr, nt, list(NULL, cn)) else est * NA_real_
  lim <- .pred_limits(est, s, type, level)
  list(est = est, se = s, lower = lim$lower, upper = lim$upper)
}

#' @rdname predict.SurvX
#' @param type for \code{parsurv}: \code{"survival"}, \code{"hazard"}, \code{"cumhaz"},
#'   \code{"quantile"} or \code{"rmst"}; for a parametric \code{csc}: \code{"cif"} or \code{"survival"}
#'   (event-free probability).
#' @param p probabilities for \code{type = "quantile"}: the time by which a fraction \code{p} has
#'   failed (\code{p = 0.5}: the median).
#' @param se also return delta-method standard errors and pointwise limits.
#' @param average return the average over the rows of \code{newdata} (or of the data).
S7::method(predict, ParSurvFit) <- function(object, newdata = NULL, times = NULL, type = c("survival", "hazard", "cumhaz", "quantile", "rmst"),
                                           p = 0.5, se = FALSE, level = 0.95, average = FALSE, ...) {
  type <- match.arg(type)
  at <- if (type == "quantile") p else times
  if (is.null(at) || !length(at)) stop("times must be given (p for type = \"quantile\")", call. = FALSE)
  if (type == "quantile" && any(!(at > 0 & at < 1))) stop("p must be between 0 and 1", call. = FALSE)
  .par_predict_X(object, .sx_newX(object, newdata), as.double(at), type, se, level, average)
}

# CIF of cause k (0: event-free probability) of a parametric cause-specific fit at design rows Xn
.parcsc_predict_X <- function(object, Xn, times, k, se, level, average = FALSE) {
  if (any(!is.finite(times)) || any(times < 0)) stop("times must be finite and non-negative", call. = FALSE)
  o <- order(times); ts <- as.double(times[o])
  K <- object@K; nn <- nrow(Xn); nt <- length(ts)
  th <- unlist(lapply(object@models, function(m) unname(m@coefficients)))
  V <- unlist(lapply(object@models, function(m) as.double(t(m@vcov))))
  usevar <- isTRUE(se) && all(is.finite(V))
  cells <- if (average) nt else nn * nt
  r <- .C("SX_parcsc_predict", as.integer(K), vapply(object@models, function(m) m@dist_code, integer(1)),
          vapply(object@models, function(m) m@aft, integer(1)), as.integer(ncol(Xn)), as.double(th),
          as.integer(usevar), if (usevar) as.double(V) else double(1), as.integer(k), as.integer(nn), .rm(Xn),
          as.integer(nt), ts, as.double(max(object@time, ts)), as.integer(average),
          out = double(cells), se = double(cells), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "predict")
  nr <- if (average) 1L else nn
  back <- order(o); cn <- paste0("t=", times)
  est <- .matrix_rm(r$out, nr, nt)[, back, drop = FALSE]; dimnames(est) <- list(NULL, cn)
  if (!isTRUE(se)) return(est)
  s <- if (usevar) .matrix_rm(r$se, nr, nt)[, back, drop = FALSE] else est * NA_real_
  dimnames(s) <- list(NULL, cn)
  lim <- .pred_limits(est, s, "cif", level)
  list(est = est, se = s, lower = lim$lower, upper = lim$upper)
}

S7::method(predict, ParCSCFit) <- function(object, newdata = NULL, times, cause = 1L, type = c("cif", "survival"),
                                          se = FALSE, level = 0.95, average = FALSE, ...) {
  type <- match.arg(type)
  k <- if (type == "survival") 0L else as.integer(cause)
  if (k > object@K) stop("cause ", k, " does not occur in the data", call. = FALSE)
  .parcsc_predict_X(object, .sx_newX(object, newdata), times, k, se, level, average)
}

# ----------------------------------------------------------------------------------------------
# extractors of the parametric cause-specific fit
# ----------------------------------------------------------------------------------------------
S7::method(coef, ParCSCFit) <- function(object, ...) lapply(object@models, function(m) m@coefficients)
S7::method(vcov, ParCSCFit) <- function(object, ...) lapply(object@models, function(m) m@vcov)
S7::method(nobs, ParCSCFit) <- function(object, ...) object@n
S7::method(logLik, ParCSCFit) <- function(object, ...) {
  structure(sum(vapply(object@models, function(m) m@loglik, numeric(1))),
            df = sum(vapply(object@models, function(m) length(m@coefficients), integer(1))), nobs = object@n, class = "logLik")
}
S7::method(confint, ParCSCFit) <- function(object, parm, level = 0.95, ...) lapply(object@models, function(m) confint(m, level = level))
S7::method(formula, ParCSCFit) <- function(x, ...) stats::formula(x@terms)
S7::method(getCall, ParCSCFit) <- function(x, ...) x@call
S7::method(terms, ParCSCFit) <- function(x, ...) x@terms
S7::method(model.matrix, ParCSCFit) <- function(object, ...) object@X

.par_estimates <- function(m, cause, level) {
  est <- m@coefficients; nm <- names(est); p <- ncol(m@X)
  role <- c("baseline", rep("beta", p), rep("shape", length(est) - 1L - p))
  sc <- ifelse(role == "beta", if (m@scale == "AFT") "log time ratio" else if (is.na(cause)) "log hazard ratio" else "log cause-specific hazard ratio",
               ifelse(role == "baseline", if (m@scale == "AFT") "log time scale" else "log baseline rate", nm))
  .est_frame(cause, nm, role, sc, est, .se_of(m), level, role == "beta")
}

# baseline parameters on their natural scale, with Wald limits computed on the estimation scale
.natural_baseline <- function(m, level) {
  e <- .par_estimates(m, NA_integer_, level)
  e <- e[e$role != "beta", , drop = FALSE]
  rows <- lapply(seq_len(nrow(e)), function(i) {
    r <- e[i, ]
    if (r$role == "baseline") {
      if (m@scale == "AFT") list("time scale exp(intercept), all covariates 0", exp(r$estimate), exp(r$lower), exp(r$upper))
      else list("baseline rate exp(intercept), all covariates 0", exp(r$estimate), exp(r$lower), exp(r$upper))
    } else if (r$term == "shape") list("Gompertz shape a (a < 0: a fraction never fails)", r$estimate, r$lower, r$upper)
    else if (r$term == "log(shape)") list("Weibull shape (1: constant hazard, > 1 rising, < 1 falling)", exp(r$estimate), exp(r$lower), exp(r$upper))
    else list(paste0("scale sigma of log time", if (m@dist == "weibull") " (Weibull shape = 1/sigma)" else ""), exp(r$estimate), exp(r$lower), exp(r$upper))
  })
  data.frame(parameter = vapply(rows, `[[`, "", 1L), estimate = .num(vapply(rows, `[[`, 0, 2L)),
             ci = .ci_text(vapply(rows, `[[`, 0, 3L), vapply(rows, `[[`, 0, 4L), 3L), stringsAsFactors = FALSE)
}

# ----------------------------------------------------------------------------------------------
# readable output
# ----------------------------------------------------------------------------------------------

# one sentence per covariate effect with p < 0.05: what the ratio means for this model
.reading_lines <- function(e, xlevels, scale, what, level, ratio_label = "HR") {
  e <- e[e$role == "beta" & is.finite(e$p), , drop = FALSE]
  if (!nrow(e)) return(character(0))
  sig <- e[e$p < 0.05, , drop = FALSE]
  if (!nrow(sig)) return(c("Reading:", "  No covariate effect reached p < 0.05."))
  lab <- .pretty_terms(sig$term, xlevels)
  is_level <- grepl(": ", lab, fixed = TRUE) & grepl(" vs ", lab, fixed = TRUE)
  who <- ifelse(is_level, sub(": ", " (", lab, fixed = TRUE), paste0(lab, " (per one-unit increase"))
  who <- paste0(who, ")")
  out <- "Reading (effects with p < 0.05):"
  for (i in seq_len(nrow(sig))) {
    r <- sig$exp_estimate[i]
    rt <- .ratio_text(r, sig$exp_lower[i], sig$exp_upper[i])
    ci <- paste0(.pct(level), " CI ", rt$ci)
    sent <- if (scale == "AFT") {
      size <- if (r >= 2) paste0(rt$ratio, " times as long a") else paste0(round(100 * abs(r - 1)), "% ", if (r > 1) "longer" else "shorter")
      paste0(who[i], ": ", size, " time to ", what, " (time ratio ", rt$ratio, ", ", ci, ").")
    } else {
      size <- if (r >= 2) paste0(rt$ratio, " times the") else paste0(round(100 * abs(r - 1)), "% ", if (r > 1) "higher" else "lower")
      paste0(who[i], ": ", size, " rate of ", what, " (", ratio_label, " ", rt$ratio, ", ", ci, ").")
    }
    out <- c(out, paste0("  - ", sent))
  }
  out
}

.fit_line <- function(m) {
  P <- length(m@coefficients)
  aic <- -2 * m@loglik + 2 * P; bic <- -2 * m@loglik + log(m@n) * P
  out <- sprintf("Fit: log-likelihood %.2f | AIC %.2f | BIC %.2f | %d parameters", m@loglik, aic, bic, P)
  df <- ncol(m@X)
  if (df > 0L && is.finite(m@loglik_null)) {
    lr <- 2 * (m@loglik - m@loglik_null)
    pv <- .pval(stats::pchisq(lr, df, lower.tail = FALSE))
    out <- c(out, sprintf("     likelihood ratio test against no covariates: chi-square %.2f on %d df, p %s",
                          lr, df, if (startsWith(pv, "<")) sub("<", "< ", pv) else paste("=", pv)))
  }
  out
}

.par_model_label <- function(dist, scale) {
  paste0(.par_label[[dist]], ", ", if (scale == "AFT") "accelerated failure time (AFT)" else "proportional hazards (PH)")
}

.render_par <- function(x, level, detail) {
  e <- estimates(x, level = level)
  cens <- sum(x@status == 0L)
  out <- c(paste0("Parametric survival regression: ", .par_model_label(x@dist, x@scale)),
           paste0("  Subjects: ", x@n, "   events: ", x@nevent, "   censored: ", cens,
                  if (!is.null(x@entry)) "   delayed entry: yes" else ""))
  eb <- e[e$role == "beta", , drop = FALSE]
  if (nrow(eb)) {
    lab <- if (x@scale == "AFT") c("TR", "Effect of covariates: time ratio (TR)") else c("HR", "Effect of covariates: hazard ratio (HR)")
    out <- c(out, "", lab[2L], .effect_lines(eb, lab[1L], level, x@xlevels), "",
             .reading_lines(e, x@xlevels, x@scale, "the event", level))
  }
  if (detail) {
    out <- c(out, "", "Baseline distribution (natural scale)", .table_lines(stats::setNames(.natural_baseline(x, level), c("", "estimate", paste(.pct(level), "CI")))))
  }
  out <- c(out, "", .fit_line(x))
  if (detail)
    out <- c(out, "",
             if (x@scale == "AFT") "TR > 1: the event happens later (times are stretched); TR < 1: earlier."
             else "HR > 1: the event occurs at a higher rate at every time; HR < 1: at a lower rate.",
             "Numeric covariates: effect of a one-unit increase; categorical covariates: level against the reference level.",
             "Compare distributions with AIC(fit1, fit2) and check the fit with hazard_plot(fit).")
  c(out, .status_lines(x@converged, x@flags))
}

.render_parcsc <- function(x, level, detail) {
  e <- estimates(x, level = level)
  out <- c("Parametric cause-specific hazards",
           paste0("  Subjects: ", x@n, "   censored: ", sum(x@status == 0L),
                  if (!is.null(x@entry)) "   delayed entry: yes" else ""))
  for (k in seq_len(x@K)) {
    m <- x@models[[k]]
    ek <- e[e$cause == k, , drop = FALSE]
    out <- c(out, "", paste0("Cause ", k, " (", m@nevent, " events): ", .par_model_label(m@dist, m@scale)))
    eb <- ek[ek$role == "beta", , drop = FALSE]
    if (nrow(eb)) {
      lab <- if (m@scale == "AFT") c("TR", "  time ratio (TR)") else c("csHR", "  cause-specific hazard ratio (csHR)")
      out <- c(out, lab[2L], .effect_lines(eb, lab[1L], level, x@xlevels),
               .reading_lines(ek, x@xlevels, m@scale, paste0("cause ", k, " among subjects still event-free"), level, "csHR"))
    }
    if (detail) out <- c(out, .table_lines(stats::setNames(.natural_baseline(m, level), c("", "estimate", paste(.pct(level), "CI")))))
    out <- c(out, .fit_line(m), .status_lines(m@converged, m@flags))
  }
  ll <- logLik(x)
  out <- c(out, "", sprintf("All causes: log-likelihood %.2f | AIC %.2f", as.numeric(ll), -2 * as.numeric(ll) + 2 * attr(ll, "df")),
           "Check: the cumulative incidences of all causes and the event-free probability add up to one (coherent by construction).")
  if (detail)
    out <- c(out, "",
             "csHR: effect on the rate of one cause among subjects who are still event-free.",
             "A covariate can raise the rate of a cause and still lower its cumulative incidence,",
             "because the incidence also depends on the other causes; use predict() or cif_plot() for absolute risks.")
  out
}

S7::method(print, ParCSCFit) <- function(x, level = 0.95, ...) {
  cat(.render_parcsc(x, level, FALSE), sep = "\n")
  invisible(x)
}
S7::method(summary, ParCSCFit) <- function(object, level = 0.95, ...) {
  FitSummary(call = object@call, coefficients = estimates(object, level = level), tables = list(),
             text = .render_parcsc(object, level, TRUE))
}

# ----------------------------------------------------------------------------------------------
# non-parametric hazard
# ----------------------------------------------------------------------------------------------

# Nelson-Aalen increments of cause k (delayed entry allowed): distinct event times, events, at risk
.na_increments <- function(time, status, k, entry = NULL) {
  ev <- status == k
  et <- sort(unique(time[ev]))
  if (!length(et)) return(NULL)
  d <- tabulate(match(time[ev], et), length(et))
  Y <- vapply(et, function(s) sum(time >= s & (if (is.null(entry)) TRUE else entry < s)), numeric(1))
  list(time = et, d = as.numeric(d), Y = Y)
}

.na_frame <- function(inc, level) {
  H <- cumsum(inc$d / inc$Y); se <- sqrt(cumsum(inc$d / inc$Y^2))
  zq <- .zq(level)
  data.frame(time = inc$time, n.risk = inc$Y, n.event = inc$d, cumhaz = H, se = se,
             lower = H * exp(-zq * se / H), upper = H * exp(zq * se / H))
}

.kernel_frame <- function(inc, bw, boundary, grid, level) {
  x <- seq(0, max(inc$time), length.out = grid)
  r <- .C("SX_kernel_hazard", as.integer(length(inc$time)), as.double(inc$time), as.double(inc$d), as.double(inc$Y),
          as.double(if (is.null(bw)) -1 else bw), as.integer(boundary), as.integer(grid), as.double(x),
          h = double(grid), se = double(grid), bw = double(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "hazard_np")
  zq <- .zq(level); h <- r$h; se <- r$se
  ok <- h > 0
  lo <- ifelse(ok, h * exp(-zq * se / h), NA_real_); up <- ifelse(ok, h * exp(zq * se / h), NA_real_)
  list(frame = data.frame(time = x, hazard = h, se = se, lower = lo, upper = up), bw = r$bw)
}

#' Non-parametric hazard: Nelson-Aalen and kernel-smoothed hazard
#'
#' The Nelson-Aalen estimate of the cumulative hazard and the kernel-smoothed hazard rate, for
#' every group of the right-hand side and, with competing risks, for every cause (cause-specific
#' hazards, the other causes treated as censoring). Delayed entry is allowed.
#'
#' @details The cumulative hazard is \eqn{\hat H(t) = \sum_{t_j \le t} d_j / Y_j} with variance
#' \eqn{\sum d_j / Y_j^2}; limits on the log scale. The smoothed hazard is
#' \eqn{\hat h(t) = b^{-1}\sum_j K((t - t_j)/b)\, d_j / Y_j} with the Epanechnikov kernel and, near
#' the ends of follow-up, the boundary kernels of Mueller and Wang (1994), so the estimate is not
#' biased downwards at time zero. Its variance is \eqn{b^{-2}\sum_j K^2 d_j / Y_j^2}; limits on the
#' log scale. The bandwidth \eqn{b} minimises the least-squares cross-validation criterion of
#' Patil (1993) unless \code{bw} is given.
#'
#' @param formula \code{Event(time, status) ~ 1} or \code{~ group}.
#' @param data a data frame.
#' @param smooth compute the kernel-smoothed hazard.
#' @param bw bandwidth; \code{NULL} chooses it by cross-validation, separately for every curve.
#' @param boundary use boundary kernels near the ends of follow-up.
#' @param level confidence level.
#' @param grid number of time points of the smoothed curve.
#' @return an S7 object of class \code{HazardFit}: \code{x@cumhaz} and \code{x@hazard} are lists of
#'   data frames, one per group and cause; \code{x@bandwidth} the bandwidths. Methods \code{print}
#'   and \code{plot}.
#' @references
#' Nelson W (1972) Theory and applications of hazard plotting for censored failure data.
#'   Technometrics 14:945-966.
#' Aalen O (1978) Nonparametric inference for a family of counting processes. Ann Stat 6:701-726.
#' Ramlau-Hansen H (1983) Smoothing counting process intensities by means of kernel functions.
#'   Ann Stat 11:453-466.
#' Mueller HG, Wang JL (1994) Hazard rate estimation under random censoring with varying kernels
#'   and bandwidths. Biometrics 50:61-76.
#' Patil PN (1993) On the least squares cross-validation bandwidth in hazard rate estimation.
#'   Ann Stat 21:1792-1810.
#' @examples
#' set.seed(3)
#' d <- data.frame(time = rweibull(300, 1.5, 3), status = 1, arm = rep(c("A", "B"), 150))
#' hz <- hazard_np(Event(time, status) ~ arm, d)
#' hz
#' plot(hz)
#' plot(hz, type = "cumhaz")
#' @export
hazard_np <- function(formula, data, smooth = TRUE, bw = NULL, boundary = TRUE, level = 0.95, grid = 201L) {
  cl <- match.call()
  f <- .sx_groups(formula, data, what = "hazard_np()", entry_ok = TRUE)
  K <- max(f$status, 1L)
  cum <- list(); haz <- list(); bws <- NULL
  for (g in levels(f$group)) for (k in seq_len(K)) {
    sel <- f$group == g
    inc <- .na_increments(f$time[sel], f$status[sel], k, if (is.null(f$entry)) NULL else f$entry[sel])
    if (is.null(inc)) next
    key <- if (K > 1L) paste0(g, " | cause ", k) else g
    cum[[key]] <- cbind(group = g, cause = k, .na_frame(inc, level))
    if (smooth && length(inc$time) >= 2L) {
      kf <- .kernel_frame(inc, bw, boundary, grid, level)
      haz[[key]] <- cbind(group = g, cause = k, kf$frame)
      bws <- rbind(bws, data.frame(group = g, cause = k, bandwidth = kf$bw, events = sum(inc$d)))
    }
  }
  HazardFit(call = cl, cumhaz = cum, hazard = haz, bandwidth = bws, K = as.integer(K), level = level,
            time = f$time, status = f$status, group = f$group)
}

S7::method(print, HazardFit) <- function(x, ...) {
  cat("Non-parametric hazard", if (x@K > 1L) "(cause-specific)" else "", "\n")
  cat("  Nelson-Aalen cumulative hazard; kernel-smoothed hazard (Epanechnikov, boundary-corrected)\n\n")
  if (!is.null(x@bandwidth)) {
    b <- x@bandwidth
    df <- data.frame(curve = if (x@K > 1L) paste0(b$group, ", cause ", b$cause) else b$group,
                     events = b$events, bandwidth = .num(b$bandwidth, 3L), stringsAsFactors = FALSE)
    names(df) <- c("", "events", "bandwidth")
    cat(.table_lines(df), sep = "\n")
    cat("\n  The bandwidth is the half-width of the time window over which events are averaged;\n")
    cat("  it was chosen by cross-validation unless given. plot(x) draws the hazard over time.\n")
  }
  invisible(x)
}

# ----------------------------------------------------------------------------------------------
# plots of hazards
# ----------------------------------------------------------------------------------------------
.band <- function(x, lo, hi, col) {
  ok <- is.finite(lo) & is.finite(hi)
  if (sum(ok) < 2L) return(invisible())
  graphics::polygon(c(x[ok], rev(x[ok])), c(lo[ok], rev(hi[ok])), col = grDevices::adjustcolor(col, 0.15), border = NA)
}

S7::method(plot, HazardFit) <- function(x, type = c("hazard", "cumhaz"), conf.int = TRUE, cause = NULL, col = NULL, lwd = 2,
                                       xlab = "Time", ylab = NULL, main = NULL, legend = TRUE, cex = 0.85, ...) {
  type <- match.arg(type)
  curves <- if (type == "hazard") x@hazard else x@cumhaz
  if (!length(curves)) stop("nothing to plot", if (type == "hazard") ": no smoothed hazard was computed" else "", call. = FALSE)
  causes <- if (is.null(cause)) seq_len(x@K) else as.integer(cause)
  groups <- levels(x@group)
  col <- if (is.null(col)) rep_len(.okabe_ito, length(groups)) else rep_len(col, length(groups))
  if (is.null(ylab)) ylab <- if (type == "hazard") "Hazard (events per unit time)" else "Cumulative hazard"
  op <- graphics::par(mfrow = c(1L, length(causes)), cex = cex, mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  for (k in causes) {
    cs <- Filter(function(d) d$cause[1L] == k, curves)
    if (!length(cs)) next
    yv <- unlist(lapply(cs, function(d) c(d[[type]], if (conf.int) d$upper)))
    graphics::plot(range(unlist(lapply(cs, `[[`, "time"))), c(0, max(yv[is.finite(yv)], 0)), type = "n", xlab = xlab, ylab = ylab,
                   main = if (is.null(main)) (if (x@K > 1L) paste("Cause", k) else "") else main, ...)
    for (d in cs) {
      i <- match(d$group[1L], groups)
      if (conf.int) .band(d$time, d$lower, d$upper, col[i])
      graphics::lines(d$time, d[[type]], col = col[i], lwd = lwd, type = if (type == "cumhaz") "s" else "l")
    }
    if (legend && length(groups) > 1L) graphics::legend("topleft", legend = groups, col = col, lwd = lwd, bty = "n")
  }
  invisible(x)
}

#' Fitted hazard, cumulative hazard and survival of parametric models
#'
#' Draws the fitted hazard (or cumulative hazard, or survival) of one or several parametric models
#' at typical covariate values or at chosen covariate values, with delta-method confidence bands,
#' and overlays the non-parametric estimate of the same data (\code{\link{hazard_np}}) as a
#' reference. Agreement of the fitted and the non-parametric curve is the visual check of the
#' chosen distribution.
#'
#' @param x a model fitted by \code{\link{parsurv}} or \code{\link{csc}} with a parametric baseline,
#'   or a named list of such models fitted to the same data.
#' @param newdata data frame of covariate values, one curve per row; default typical values (mean
#'   of numeric covariates, most frequent level of categorical covariates).
#' @param type \code{"hazard"}, \code{"cumhaz"} or \code{"survival"}.
#' @param cause causes to show (cause-specific models); default all.
#' @param reference overlay the non-parametric estimate of all subjects (Nelson-Aalen, smoothed
#'   hazard or Kaplan-Meier). Meaningful as a check at typical covariate values or without covariates.
#' @param conf.int draw pointwise confidence bands.
#' @param level confidence level.
#' @param times time grid; default 201 points up to the largest observed time.
#' @param col,lwd colours and line width of the models.
#' @param xlab,ylab,main axis labels and title.
#' @param legend draw a legend.
#' @param cex character expansion.
#' @param ... passed to \code{plot.default}.
#' @return invisibly, a data frame of the plotted curves.
#' @references Klein JP, Moeschberger ML (2003) Survival Analysis, 2nd ed. Springer, Chapter 12
#'   (checking the fit of parametric models).
#' @examples
#' set.seed(1)
#' n <- 300
#' t <- rweibull(n, shape = 1.5, scale = 4); cens <- runif(n, 0, 8)
#' d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = rnorm(n))
#' fits <- list(Weibull = parsurv(Event(time, status) ~ x, d),
#'              exponential = parsurv(Event(time, status) ~ x, d, baseline = "exponential"))
#' hazard_plot(fits)
#' hazard_plot(fits, type = "cumhaz")
#' @export
hazard_plot <- function(x, newdata = NULL, type = c("hazard", "cumhaz", "survival"), cause = NULL, reference = TRUE,
                        conf.int = TRUE, level = 0.95, times = NULL, col = NULL, lwd = 2, xlab = "Time", ylab = NULL,
                        main = NULL, legend = TRUE, cex = 0.85, ...) {
  type <- match.arg(type)
  models <- if (is.list(x) && !S7::S7_inherits(x)) x else list(x)
  ok <- vapply(models, function(m) S7::S7_inherits(m, ParSurvFit) || S7::S7_inherits(m, ParCSCFit), logical(1))
  if (!all(ok)) stop("hazard_plot() takes models fitted by parsurv() or csc() with a parametric baseline", call. = FALSE)
  if (is.null(names(models))) names(models) <- vapply(models, .par_name, character(1))
  ref <- models[[1L]]
  K <- if (S7::S7_inherits(ref, ParCSCFit)) ref@K else 1L
  causes <- if (type == "survival") 0L else if (is.null(cause)) seq_len(K) else as.integer(cause)
  tg <- if (is.null(times)) seq(0, max(ref@time), length.out = 201L) else sort(times)
  tg_pos <- if (type == "hazard") tg[tg > 0] else tg
  rows <- if (is.null(newdata)) { pf <- .profile(ref); list(X = pf$X, labels = "typical covariate values") }
          else list(X = .sx_newX(ref, newdata), labels = .newdata_labels(ref, newdata))
  curves <- list()
  for (nm in names(models)) {
    m <- models[[nm]]
    for (k in causes) for (i in seq_len(nrow(rows$X))) {
      Xi <- rows$X[i, , drop = FALSE]
      pr <- if (S7::S7_inherits(m, ParCSCFit)) {
        if (type == "survival") .parcsc_predict_X(m, Xi, tg_pos, 0L, conf.int, level)
        else .par_predict_X(m@models[[k]], Xi, tg_pos, type, conf.int, level)
      } else .par_predict_X(m, Xi, tg_pos, type, conf.int, level)
      est <- if (conf.int) pr$est else pr
      curves[[length(curves) + 1L]] <- data.frame(model = nm, cause = k, profile = rows$labels[i], time = tg_pos,
        value = drop(est), lower = if (conf.int) drop(pr$lower) else NA_real_, upper = if (conf.int) drop(pr$upper) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, curves)
  # reference: non-parametric estimate of all subjects
  refd <- NULL
  if (reference) {
    st <- if (S7::S7_inherits(ref, ParCSCFit)) ref@status else ref@status
    refd <- do.call(rbind, lapply(causes, function(k) {
      kk <- if (k == 0L) NULL else k
      if (type == "survival") {
        inc <- .na_increments(ref@time, as.integer(st > 0L), 1L, ref@entry)
        if (is.null(inc)) return(NULL)
        s <- cumprod(1 - inc$d / inc$Y)
        return(data.frame(cause = k, time = c(0, inc$time), value = c(1, s), step = TRUE))
      }
      inc <- .na_increments(ref@time, st, kk, ref@entry)
      if (is.null(inc)) return(NULL)
      if (type == "cumhaz") { nf <- .na_frame(inc, level); return(data.frame(cause = k, time = c(0, nf$time), value = c(0, nf$cumhaz), step = TRUE)) }
      if (length(inc$time) < 2L) return(NULL)
      kf <- .kernel_frame(inc, NULL, TRUE, 201L, level)$frame
      data.frame(cause = k, time = kf$time, value = kf$hazard, step = FALSE)
    }))
  }
  keys <- unique(paste(d$model, d$profile, sep = " | "))
  lab <- if (length(unique(d$profile)) == 1L) unique(d$model) else if (length(unique(d$model)) == 1L) unique(d$profile) else keys
  col <- if (is.null(col)) rep_len(.okabe_ito, length(keys)) else rep_len(col, length(keys))
  if (is.null(ylab)) ylab <- switch(type, hazard = "Hazard (events per unit time)", cumhaz = "Cumulative hazard", survival = "Event-free probability")
  op <- graphics::par(mfrow = c(1L, length(causes)), cex = cex, mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  for (k in causes) {
    dk <- d[d$cause == k, , drop = FALSE]; rk <- if (is.null(refd)) NULL else refd[refd$cause == k, , drop = FALSE]
    yv <- c(dk$value, if (conf.int) dk$upper, if (!is.null(rk)) rk$value)
    yv <- yv[is.finite(yv)]
    ylim <- if (type == "survival") c(0, 1) else c(0, if (length(yv)) stats::quantile(yv, 0.995, names = FALSE) * 1.05 else 1)
    ttl <- if (!is.null(main)) main else if (k > 0L && K > 1L) paste0("Cause ", k) else ""
    graphics::plot(range(tg), ylim, type = "n", xlab = xlab, ylab = ylab, main = ttl, ...)
    if (!is.null(rk) && nrow(rk)) graphics::lines(rk$time, rk$value, type = if (rk$step[1L]) "s" else "l", col = "grey45", lwd = lwd, lty = 2)
    for (j in seq_along(keys)) {
      dj <- dk[paste(dk$model, dk$profile, sep = " | ") == keys[j], , drop = FALSE]
      if (!nrow(dj)) next
      if (conf.int) .band(dj$time, dj$lower, dj$upper, col[j])
      graphics::lines(dj$time, dj$value, col = col[j], lwd = lwd)
    }
    if (legend) {
      lg <- c(lab, if (!is.null(rk) && nrow(rk)) "non-parametric (all subjects)")
      graphics::legend("topleft", legend = lg, col = c(col, "grey45")[seq_along(lg)], lwd = lwd,
                       lty = c(rep(1L, length(lab)), 2L)[seq_along(lg)], bty = "n")
    }
  }
  invisible(d)
}

.par_name <- function(m) {
  if (S7::S7_inherits(m, ParCSCFit)) paste0(paste(unique(.par_label[m@dist]), collapse = "/"), " cause-specific")
  else paste0(.par_label[[m@dist]], " ", m@scale)
}
