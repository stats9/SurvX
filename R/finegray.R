#' Fine-Gray model for the subdistribution hazard
#'
#' Proportional subdistribution hazards regression with inverse-probability-of-censoring
#' weights and the robust sandwich variance of Fine and Gray (1999).
#'
#' @param formula \code{Event(time, status) ~ covariates}; status 0 = censored, 1..K = cause.
#' @param data a data frame.
#' @param cause cause of interest.
#' @param variance \code{"robust"} (default, Fine and Gray 1999, Section 4) or \code{"model"}
#'   (inverse information).
#' @param max_iter,tol iteration limit and convergence tolerance.
#' @return an S7 object of class \code{FineGrayFit}. \code{predict(fit, newdata, times)} gives the
#'   cumulative incidence \eqn{1 - \exp\{-\Lambda_0(t) e^{z'\beta}\}}.
#' @references Fine JP, Gray RJ (1999) J Am Stat Assoc 94:496-509.
#' @export
fg <- function(formula, data, cause = 1L, variance = c("robust", "model"), max_iter = 50L, tol = 1e-10) {
  cl <- match.call(); variance <- match.arg(variance)
  f <- .sx_frame(formula, data)
  n <- length(f$time); p <- ncol(f$X)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  r <- .C("SX_fg", as.integer(n), as.integer(p), f$time, f$status, .rm(f$X), as.integer(cause),
          as.integer(max_iter), as.double(tol), beta = double(p), var_model = double(p * p), var_robust = double(p * p),
          info_d = double(3), info_i = integer(4), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "fg")
  nm <- colnames(f$X)
  vm <- matrix(r$var_model, p, p, dimnames = list(nm, nm)); vr <- matrix(r$var_robust, p, p, dimnames = list(nm, nm))
  FineGrayFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
              coefficients = stats::setNames(r$beta, nm), vcov = if (variance == "robust") vr else vm,
              loglik = r$info_d[2L], n = n, nevent = r$info_i[4L], iterations = r$info_i[1L],
              converged = as.logical(r$info_i[2L]), flags = r$info_i[3L],
              time = f$time, status = f$status, X = f$X, cause = as.integer(cause), variance = variance,
              vcov_model = vm, vcov_robust = vr, loglik_null = r$info_d[1L])
}

S7::method(predict, FineGrayFit) <- function(object, newdata = NULL, times, ...) {
  n <- object@n; p <- length(object@coefficients)
  r <- .C("SX_fg_baseline", as.integer(n), as.integer(p), object@time, object@status, .rm(object@X), object@cause,
          as.double(object@coefficients), time = double(n), cumhaz = double(n), m = integer(1), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "baseline subdistribution hazard")
  k <- seq_len(r$m)
  H <- .step_at(r$time[k], r$cumhaz[k], times, 0)
  lp <- drop(.sx_newX(object, newdata) %*% object@coefficients)
  out <- 1 - exp(-outer(exp(lp), H))
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}

#' Cause-specific Cox models
#'
#' One Cox model per cause, with the other causes treated as censoring, combined into
#' cumulative incidence functions.
#'
#' @param formula \code{Event(time, status) ~ covariates}; status 0 = censored, 1..K = cause.
#' @param data a data frame.
#' @param ties \code{"efron"} (default) or \code{"breslow"}.
#' @param product_limit \code{TRUE} (default): overall survival as a product-limit, as in
#'   \pkg{riskRegression}; \code{FALSE}: \eqn{\exp\{-\sum_k \Lambda_k(t)\}}.
#' @return an S7 object of class \code{CSCFit}; \code{x@models} holds the fit of each cause.
#'   \code{predict(fit, newdata, times, cause)} gives the cumulative incidence.
#' @references Prentice RL et al. (1978) Biometrics 34:541-554.
#'   Benichou J, Gail MH (1990) Biometrics 46:813-826.
#'   Ozenne B et al. (2017) R Journal 9(2):440-460.
#' @export
csc <- function(formula, data, ties = c("efron", "breslow"), product_limit = TRUE) {
  cl <- match.call(); ties <- match.arg(ties)
  f <- .sx_frame(formula, data)
  K <- max(f$status)
  if (ncol(f$X) < 1L) stop("at least one covariate is needed", call. = FALSE)
  models <- lapply(seq_len(K), function(k) .cox_fit(f$time, as.integer(f$status == k), f$X, NULL, ties))
  names(models) <- paste0("cause", seq_len(K))
  CSCFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts, models = models, K = as.integer(K),
         ties = ties, product_limit = isTRUE(product_limit), n = length(f$time), time = f$time, status = f$status, X = f$X)
}

S7::method(predict, CSCFit) <- function(object, newdata = NULL, times, cause = 1L, ...) {
  Xn <- .sx_newX(object, newdata)
  n <- object@n; p <- ncol(object@X); K <- object@K; nn <- nrow(Xn)
  B <- do.call(rbind, lapply(object@models, function(m) m$coefficients))
  r <- .C("SX_csc_cif", as.integer(n), as.integer(p), object@time, object@status, .rm(object@X), as.integer(K),
          .rm(B), as.integer(object@ties == "efron"), as.integer(nn), .rm(Xn), as.integer(object@product_limit),
          time = double(n), cif = double(nn * n * K), m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "csc prediction")
  m <- r$m
  arr <- array(r$cif[seq_len(nn * m * K)], dim = c(K, m, nn))      # [cause, time, subject]
  idx <- findInterval(times, r$time[seq_len(m)])
  out <- matrix(0, nn, length(times), dimnames = list(NULL, paste0("t=", times)))
  for (j in seq_along(times)) if (idx[j] > 0L) out[, j] <- arr[cause, idx[j], ]
  out
}
