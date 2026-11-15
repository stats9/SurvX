#' Cox proportional hazards model
#'
#' Partial-likelihood fit by Newton-Raphson with step halving.
#'
#' @param formula \code{Event(time, status) ~ covariates}; status 0/1.
#' @param data a data frame.
#' @param ties \code{"efron"} (default) or \code{"breslow"}.
#' @param weights optional case weights, one per row of \code{data}.
#' @param max_iter,tol iteration limit and convergence tolerance.
#' @return an S7 object of class \code{CoxFit}. Methods: \code{print}, \code{summary}, \code{coef},
#'   \code{vcov}, \code{logLik}, \code{nobs}, \code{predict}.
#' @references Cox DR (1972) J R Stat Soc B 34:187-220. Breslow N (1974) Biometrics 30:89-99.
#'   Efron B (1977) J Am Stat Assoc 72:557-565.
#' @export
cox <- function(formula, data, ties = c("efron", "breslow"), weights = NULL, max_iter = 50L, tol = 1e-10) {
  cl <- match.call(); ties <- match.arg(ties)
  f <- .sx_frame(formula, data)
  if (any(f$status > 1L)) stop("cox() needs a 0/1 status; use csc() for cause-specific Cox models", call. = FALSE)
  n <- length(f$time); p <- ncol(f$X)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  if (!is.null(weights)) {
    omit <- attr(f$mf, "na.action")
    if (!is.null(omit)) weights <- weights[-omit]
    if (length(weights) != n) stop("weights must have one value per row of data", call. = FALSE)
  }
  r <- .cox_fit(f$time, f$status, f$X, weights, ties, max_iter, tol)
  CoxFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
         coefficients = r$coefficients, vcov = r$vcov, loglik = r$loglik, n = n, nevent = r$nevent,
         iterations = r$iterations, converged = r$converged, flags = r$flags,
         time = f$time, status = f$status, X = f$X, ties = ties, weights = weights, loglik_null = r$loglik_null)
}

.cox_fit <- function(time, status, X, weights, ties, max_iter = 50L, tol = 1e-10) {
  n <- length(time); p <- ncol(X)
  r <- .C("SX_cox", as.integer(n), as.integer(p), as.double(time), as.integer(status), .rm(X),
          as.integer(!is.null(weights)), if (is.null(weights)) double(n) else as.double(weights),
          as.integer(ties == "efron"), as.integer(max_iter), as.double(tol),
          beta = double(p), var = double(p * p), info_d = double(3), info_i = integer(4), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "cox")
  nm <- colnames(X)
  list(coefficients = stats::setNames(r$beta, nm), vcov = matrix(r$var, p, p, dimnames = list(nm, nm)),
       loglik = r$info_d[2L], loglik_null = r$info_d[1L], nevent = r$info_i[4L], iterations = r$info_i[1L],
       converged = as.logical(r$info_i[2L]), flags = r$info_i[3L])
}

.cox_baseline <- function(object) {
  n <- object@n; p <- length(object@coefficients); w <- object@weights
  r <- .C("SX_cox_baseline", as.integer(n), as.integer(p), object@time, object@status, .rm(object@X),
          as.integer(!is.null(w)), if (is.null(w)) double(n) else as.double(w), as.integer(object@ties == "efron"),
          as.double(object@coefficients), time = double(n), cumhaz = double(n), m = integer(1), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "baseline hazard")
  k <- seq_len(r$m)
  data.frame(time = r$time[k], cumhaz = r$cumhaz[k])
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
  H <- .step_at(bh$time, bh$cumhaz, times, 0)
  out <- exp(-outer(exp(lp), H))
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}
