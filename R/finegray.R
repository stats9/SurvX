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
  f <- .sx_frame(formula, data, what = "fg()")
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

#' Cause-specific hazard models: Cox or parametric
#'
#' One hazard model per cause, with the other causes treated as censoring, combined into
#' cumulative incidence functions. With \code{baseline = "cox"} (default) every cause gets a Cox
#' model with an unspecified baseline hazard; with a parametric \code{baseline} every cause gets a
#' parametric hazard model fitted by \code{\link{parsurv}}, and the cumulative incidence
#' \eqn{F_k(t) = \int_0^t h_k(u) S(u)\,du} is computed by numerical integration with delta-method
#' standard errors. Both are coherent by construction: the incidences of all causes and the
#' event-free probability add up to one.
#'
#' @param formula \code{Event(time, status) ~ covariates}; status 0 = censored, 1..K = cause.
#'   With a parametric baseline, \code{Event(time, status, entry = )} allows delayed entry.
#' @param data a data frame.
#' @param ties \code{"efron"} (default) or \code{"breslow"} (Cox baseline only).
#' @param product_limit \code{TRUE} (default): overall survival as a product-limit, as in
#'   \pkg{riskRegression}; \code{FALSE}: \eqn{\exp\{-\sum_k \Lambda_k(t)\}} (Cox baseline only).
#' @param baseline \code{"cox"}, or one of the distributions of \code{\link{parsurv}}
#'   (\code{"exponential"}, \code{"weibull"}, \code{"gompertz"}, \code{"lognormal"},
#'   \code{"loglogistic"}); a vector gives one distribution per cause.
#' @param scale \code{"PH"} or \code{"AFT"} for parametric baselines, one value or one per cause;
#'   default as in \code{\link{parsurv}}.
#' @param max_iter,tol iteration limit and convergence tolerance of the parametric fits.
#' @return \code{baseline = "cox"}: an S7 object of class \code{CSCFit}; \code{x@models} holds the
#'   fit of each cause. Parametric: an S7 object of class \code{ParCSCFit}; \code{x@models} holds
#'   the \code{\link{parsurv}} fit of each cause. \code{predict(fit, newdata, times, cause)} gives
#'   the cumulative incidence.
#' @references Prentice RL, Kalbfleisch JD, Peterson AV, Flournoy N, Farewell VT, Breslow NE (1978)
#'   The analysis of failure times in the presence of competing risks. Biometrics 34:541-554.
#'   Benichou J, Gail MH (1990) Biometrics 46:813-826.
#'   Ozenne B et al. (2017) R Journal 9(2):440-460.
#'   Hinchliffe SR, Lambert PC (2013) Flexible parametric modelling of cause specific hazards to
#'   estimate cumulative incidence functions. BMC Med Res Methodol 13:13.
#' @examples
#' set.seed(2)
#' n <- 300
#' x <- rbinom(n, 1, 0.5)
#' t1 <- rweibull(n, 1.3, 3 * exp(-0.4 * x)); t2 <- rweibull(n, 0.9, 5); cens <- runif(n, 0, 6)
#' d <- data.frame(time = pmin(t1, t2, cens), x = x)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#' csc(Event(time, status) ~ x, d)
#' fw <- csc(Event(time, status) ~ x, d, baseline = "weibull")
#' fw
#' predict(fw, data.frame(x = 0:1), times = c(1, 3), cause = 1)
#' @export
csc <- function(formula, data, ties = c("efron", "breslow"), product_limit = TRUE, baseline = "cox", scale = NULL,
                max_iter = 200L, tol = 1e-10) {
  cl <- match.call(); ties <- match.arg(ties)
  parametric <- any(baseline != "cox")
  if (parametric && any(baseline == "cox")) stop("baseline: 'cox' cannot be mixed with parametric distributions", call. = FALSE)
  f <- .sx_frame(formula, data, entry_ok = parametric, what = "csc() with a Cox baseline")
  if (parametric) return(.pcsc(cl, f, baseline, scale, max_iter, tol))
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
