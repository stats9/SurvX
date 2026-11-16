#' Cox proportional hazards model
#'
#' Partial-likelihood fit by Newton-Raphson with step halving, with strata, delayed entry or
#' time-varying covariates in (start, stop] form, case weights and a robust variance for clustered
#' data.
#'
#' @param formula \code{Event(time, status) ~ covariates}; status 0/1. For (start, stop] data or
#'   delayed entry use \code{Event(stop, status, entry = start)}. The right-hand side may contain
#'   \code{strata(...)} (a separate baseline hazard per stratum) and \code{cluster(id)} (robust
#'   variance for correlated rows).
#' @param data a data frame.
#' @param ties \code{"efron"} (default) or \code{"breslow"}.
#' @param weights optional case weights, one per row of \code{data}.
#' @param robust use the robust (sandwich) variance; \code{TRUE} whenever \code{cluster()} is given.
#' @param max_iter,tol iteration limit and convergence tolerance.
#' @return an S7 object of class \code{CoxFit}.
#' @export
cox <- function(formula, data, ties = c("efron", "breslow"), weights = NULL, robust = FALSE, max_iter = 50L, tol = 1e-10) {
  cl <- match.call(); ties <- match.arg(ties)
  f <- .sx_frame(formula, data, entry_ok = TRUE, what = "cox()", specials_ok = TRUE)
  if (any(f$status > 1L)) stop("cox() needs a 0/1 status; use csc() for cause-specific Cox models", call. = FALSE)
  n <- length(f$time); p <- ncol(f$X)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  if (!is.null(weights)) {
    omit <- attr(f$mf, "na.action")
    if (!is.null(omit)) weights <- weights[-omit]
    if (length(weights) != n) stop("weights must have one value per row of data", call. = FALSE)
  }
  robust <- isTRUE(robust) || !is.null(f$cluster)
  r <- .cox_fit(f$time, f$status, f$X, weights, ties, max_iter, tol, f$entry, f$strata)
  vm <- r$vcov
  if (robust) r$vcov <- .cox_robust(r, f$time, f$status, f$X, weights, ties, f$entry, f$strata, f$cluster)
  CoxFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
         coefficients = r$coefficients, vcov = r$vcov, loglik = r$loglik, n = n, nevent = r$nevent,
         iterations = r$iterations, converged = r$converged, flags = r$flags,
         time = f$time, status = f$status, X = f$X, ties = ties, weights = weights, loglik_null = r$loglik_null,
         entry = f$entry, strata = f$strata, strata_levels = f$strata_levels, strata_vars = f$strata_vars,
         cluster = f$cluster, robust = robust, vcov_model = vm)
}

# arguments shared by the C routines of the extended Cox model
.cox_args <- function(time, status, X, weights, entry, strata) {
  n <- length(time)
  list(as.integer(n), as.integer(ncol(X)), as.integer(!is.null(entry)), if (is.null(entry)) double(n) else as.double(entry),
       as.double(time), as.integer(status), as.integer(!is.null(strata)), if (is.null(strata)) integer(n) else as.integer(strata),
       .rm(X), as.integer(!is.null(weights)), if (is.null(weights)) double(n) else as.double(weights))
}

.cox_fit <- function(time, status, X, weights, ties, max_iter = 50L, tol = 1e-10, entry = NULL, strata = NULL) {
  p <- ncol(X)
  r <- do.call(.C, c(list("SX_cox2"), .cox_args(time, status, X, weights, entry, strata),
                     list(as.integer(ties == "efron"), as.integer(max_iter), as.double(tol),
                          beta = double(p), var = double(p * p), info_d = double(3), info_i = integer(4), rc = integer(1),
                          NAOK = TRUE, PACKAGE = "SurvX")))
  .sx_rc(r$rc, "cox")
  nm <- colnames(X)
  list(coefficients = stats::setNames(r$beta, nm), vcov = matrix(r$var, p, p, dimnames = list(nm, nm)),
       loglik = r$info_d[2L], loglik_null = r$info_d[1L], nevent = r$info_i[4L], iterations = r$info_i[1L],
       converged = as.logical(r$info_i[2L]), flags = r$info_i[3L])
}

# score residuals (n x p), Efron or Breslow form; their sum is the score
.cox_score <- function(beta, time, status, X, weights, ties, entry, strata) {
  n <- length(time); p <- ncol(X)
  r <- do.call(.C, c(list("SX_cox_score"), .cox_args(time, status, X, weights, entry, strata),
                     list(as.integer(ties == "efron"), as.double(beta), U = double(n * p), rc = integer(1),
                          NAOK = TRUE, PACKAGE = "SurvX")))
  .sx_rc(r$rc, "score residuals")
  matrix(r$U, n, p, byrow = TRUE, dimnames = list(NULL, colnames(X)))
}

# robust variance of Lin & Wei (1989): V (sum_c U_c U_c') V, with U_c the score residuals summed
# within cluster c (each row its own cluster without cluster())
.cox_robust <- function(fit, time, status, X, weights, ties, entry, strata, cluster) {
  U <- .cox_score(fit$coefficients, time, status, X, weights, ties, entry, strata)
  if (!is.null(cluster)) U <- rowsum(U, cluster, reorder = FALSE)
  V <- fit$vcov
  if (!all(is.finite(V))) return(V)
  out <- V %*% crossprod(U) %*% V
  dimnames(out) <- dimnames(V)
  out
}

# Breslow/Efron baseline cumulative hazard, per stratum (column stratum: 0-based code)
.cox_baseline_raw <- function(beta, time, status, X, weights, ties, entry, strata) {
  n <- length(time)
  r <- do.call(.C, c(list("SX_cox2_baseline"), .cox_args(time, status, X, weights, entry, strata),
                     list(as.integer(ties == "efron"), as.double(beta), time = double(n), cumhaz = double(n),
                          stratum = integer(n), m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")))
  .sx_rc(r$rc, "baseline hazard")
  k <- seq_len(r$m)
  data.frame(time = r$time[k], cumhaz = r$cumhaz[k], stratum = r$stratum[k])
}

.cox_baseline <- function(object) {
  .cox_baseline_raw(object@coefficients, object@time, object@status, object@X, object@weights, object@ties,
                    object@entry, object@strata)
}

# cumulative hazard of baseline bh at times tt for rows in strata s (0-based; NULL: one stratum)
.cumhaz_at <- function(bh, tt, s = NULL) {
  if (is.null(s)) return(.step_at(bh$time, bh$cumhaz, tt, 0))
  out <- numeric(length(tt))
  for (k in unique(s)) {
    i <- s == k; b <- bh[bh$stratum == k, , drop = FALSE]
    out[i] <- if (nrow(b)) .step_at(b$time, b$cumhaz, tt[i], 0) else 0
  }
  out
}

#' Predictions from SurvX models
#'
#' \code{predict()} methods of the model classes.
#'
#' @name predict.SurvX
#' @param object a fitted model.
#' @param newdata data frame with the covariates; the data used for fitting when \code{NULL}.
#' @param times time points. For \code{directcif} fits \code{Inf} gives the plateau.
#' @param type for \code{cox}: \code{"lp"} (linear predictor) or \code{"survival"}.
#' @param cause cause of interest (\code{csc} and \code{directcif} fits).
#' @param se for \code{directcif} fits: also return delta-method standard errors.
#' @param ... not used.
#' @return a matrix with one row per row of \code{newdata} and one column per time point
#'   (a vector for \code{type = "lp"}); with \code{se = TRUE} a list with \code{cif} and \code{se}.
NULL

S7::method(predict, CoxFit) <- function(object, newdata = NULL, times = NULL, type = c("lp", "survival"), ...) {
  type <- match.arg(type)
  lp <- drop(.sx_newX(object, newdata) %*% object@coefficients)
  if (type == "lp") return(lp)
  if (is.null(times)) stop("times are needed for type = \"survival\"", call. = FALSE)
  bh <- .cox_baseline(object)
  sn <- .sx_newstrata(object, newdata)
  out <- if (is.null(sn)) exp(-outer(exp(lp), .step_at(bh$time, bh$cumhaz, times, 0)))
         else t(vapply(seq_along(lp), function(i) exp(-exp(lp[i]) * .cumhaz_at(bh, times, rep(sn[i], length(times)))),
                       numeric(length(times))))
  if (length(times) == 1L && !is.null(sn)) out <- matrix(out, ncol = 1L)
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}
