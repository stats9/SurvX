#' Shared gamma frailty models for clustered survival data
#'
#' @description
#' When subjects come in clusters (patients in hospitals, children in families, teeth in mouths,
#' recurrent events of one person, twins), members of a cluster share unmeasured factors and their
#' event times are correlated. A shared frailty model gives every cluster an unobserved factor
#' \eqn{z_c} that multiplies the hazards of all its members,
#' \deqn{h_{j}(t \mid z_c) = z_c\, h_0(t) \exp(x_j'\beta),}
#' with \eqn{z_c} gamma-distributed with mean 1 and variance \eqn{\theta}. Write
#' \code{+ frailty(id)} in the formula of \code{\link{cox}} (semi-parametric, Breslow baseline) or of
#' \code{\link{parsurv}} (any parametric family); this page describes both.
#'
#' \strong{Use it when} you want covariate effects within clusters and a measure of how much clusters
#' differ (\eqn{\theta}), or predictions for a new or for an observed cluster. \strong{Use
#' something else when} you only need population-average effects with valid standard errors:
#' \code{cox(... + cluster(id))} (robust variance) answers that question without a model for the
#' dependence.
#'
#' @section Usage:
#' \preformatted{
#' cox(Event(time, status) ~ x + frailty(id), data)
#' parsurv(Event(time, status) ~ x + frailty(id), data, baseline = "weibull")
#' }
#'
#' @details
#' \strong{Likelihood.} The frailties are integrated out: cluster \eqn{c} with \eqn{D_c} events
#' and \eqn{\Lambda_c = \sum_j H_0(t_j) e^{x_j'\beta}} contributes
#' \deqn{\prod_j h_j(t_j)^{d_j} \prod_{m=0}^{D_c - 1}(1 + m\theta)\,(1 + \theta\Lambda_c)^{-(1/\theta + D_c)}}
#' (Clayton 1978; Klein 1992). With a parametric baseline the marginal likelihood is maximised by
#' Newton-Raphson with exact first and second derivatives. With the Breslow baseline, the baseline
#' hazard is profiled out by the EM fixed point of Nielsen et al. (1992) and Klein (1992), and the
#' profile likelihood in \eqn{(\beta, \log\theta)} is maximised by Newton-Raphson; its curvature
#' gives standard errors of \eqn{\beta} that account for the estimation of \eqn{\theta} (Murphy and
#' van der Vaart 2000). Ties are handled as in Breslow's method. \code{frailtyEM::emfrail}
#' maximises the same likelihood; the penalized likelihood of
#' \code{survival::coxph(... + frailty(id, distribution = "gamma"))} gives close but not identical
#' values.
#'
#' \strong{Reading the results.} Hazard ratios are conditional on the frailty: they compare
#' subjects of the same cluster (or of clusters with equal frailty). The population-average
#' (marginal) hazard ratio is closer to 1 and changes over time. \eqn{\theta} is printed with its
#' confidence interval, the likelihood ratio test of \eqn{\theta = 0} (with the half chi-square
#' distribution of a parameter on the boundary; Self and Liang 1987), Kendall's \eqn{\tau =
#' \theta/(\theta + 2)} between the event times of two members of a cluster, and the ratio of the
#' hazards of a cluster at the 75th and at the 25th percentile of the frailty distribution.
#'
#' \strong{Frailties of the clusters.} Given the data, \eqn{z_c} has a gamma distribution with
#' shape \eqn{1/\theta + D_c} and rate \eqn{1/\theta + \Lambda_c}; its mean is the empirical Bayes
#' estimate, shown with a 95\% posterior interval by \code{as.data.frame(fit, what = "frailty")} and
#' \code{plot(fit)}.
#'
#' \strong{Prediction.} \code{predict(fit, newdata, times, frailty = )}: \code{"marginal"} (a new
#' cluster, frailty integrated out: \eqn{(1 + \theta H)^{-1/\theta}}), \code{"conditional"}
#' (frailty 1, the average cluster: \eqn{e^{-H}}) or \code{"cluster"} (an observed cluster, given
#' its data: \eqn{\{r_c/(r_c + H)\}^{s_c}}, with the cluster in a column of \code{newdata} named as
#' in the formula). \code{type = "survival"} or \code{"cumhaz"}.
#'
#' Not available with frailty: delayed entry, weights, \code{strata()}, \code{cluster()},
#' \code{tt()} and interval-censored data.
#' @section Messages:
#' \describe{
#'   \item{\code{the frailty variance is at the boundary (0)}}{the data show no heterogeneity
#'     between clusters; the fit is then the model without frailty.}
#' }
#' @return
#' An S7 object of class \code{FrailtyFit} with \code{print()}, \code{summary()}, \code{coef()}
#' (\eqn{\beta}; for a parametric baseline also the baseline parameters), \code{vcov()},
#' \code{logLik()} (marginal; the degrees of freedom include \eqn{\theta}), \code{AIC()},
#' \code{predict()}, \code{plot()} and \code{as.data.frame(fit, what = c("coefficients",
#' "frailty"))}. \code{fit@theta}, \code{fit@theta_se} and \code{fit@lrt} hold the frailty variance,
#' its standard error and the test.
#' @references
#' Clayton DG (1978). A model for association in bivariate life tables and its application in
#'   epidemiological studies of familial tendency in chronic disease incidence. Biometrika 65:141-151.
#'
#' Klein JP (1992). Semiparametric estimation of random effects using the Cox model based on the EM
#'   algorithm. Biometrics 48:795-806.
#'
#' Nielsen GG, Gill RD, Andersen PK, Sorensen TIA (1992). A counting process approach to maximum
#'   likelihood estimation in frailty models. Scandinavian Journal of Statistics 19:25-43.
#'
#' Murphy SA, van der Vaart AW (2000). On profile likelihood. Journal of the American Statistical
#'   Association 95:449-465.
#'
#' Self SG, Liang KY (1987). Asymptotic properties of maximum likelihood estimators and likelihood
#'   ratio tests under nonstandard conditions. Journal of the American Statistical Association
#'   82:605-610.
#'
#' Duchateau L, Janssen P (2008). The Frailty Model. Springer, New York.
#'
#' Balan TA, Putter H (2019). frailtyEM: an R package for estimating semiparametric shared frailty
#'   models. Journal of Statistical Software 90(7):1-29.
#' @seealso
#' \code{\link{cox}}, \code{\link{parsurv}}; \code{cox(... + cluster(id))} for robust standard errors
#' without a frailty.
#' @examples
#' set.seed(3)
#' C <- 80; m <- 4
#' id <- rep(seq_len(C), each = m)
#' z <- rgamma(C, shape = 2, rate = 2)[id]          # frailty variance 0.5
#' x <- rnorm(C * m); trt <- rbinom(C * m, 1, 0.5)
#' t <- rweibull(C * m, 1.3, (0.1 * z * exp(0.5 * x - 0.7 * trt))^(-1 / 1.3))
#' cens <- runif(C * m, 0, 15)
#' d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x, trt, id)
#'
#' fc <- cox(Event(time, status) ~ x + trt + frailty(id), data = d)
#' fc
#' fw <- parsurv(Event(time, status) ~ x + trt + frailty(id), data = d, baseline = "weibull")
#' fw
#' head(as.data.frame(fc, what = "frailty"))
#' predict(fc, data.frame(x = 0, trt = 0:1), times = c(5, 10), frailty = "marginal")
#' plot(fc)
#' @name frailty
NULL

.frailty_var_name <- function(v) sub("^frailty\\((.*)\\)$", "\\1", v)

.frailty_ids <- function(f) {
  g <- f$frailty
  if (anyNA(g)) stop("frailty(): missing cluster identifiers", call. = FALSE)
  g <- as.factor(g)
  g <- droplevels(g)
  if (nlevels(g) < 2L) stop("frailty(): at least two clusters are needed", call. = FALSE)
  g
}

.frailty_checks <- function(f, weights = NULL) {
  if (!is.null(f$entry)) stop("frailty() cannot be combined with delayed entry", call. = FALSE)
  if (!is.null(weights)) stop("frailty() cannot be combined with weights", call. = FALSE)
  if (!is.null(f$time2)) stop("frailty() is not available for interval-censored data", call. = FALSE)
  if (any(f$status > 1L)) stop("frailty() needs a 0/1 status", call. = FALSE)
}

.frailty_lrt <- function(ll, ll0, boundary) {
  st <- if (boundary) 0 else max(0, 2 * (ll - ll0))
  list(statistic = st, p.value = 0.5 * stats::pchisq(st, 1, lower.tail = FALSE))
}

# semi-parametric: Breslow baseline profiled out
.cox_frailty <- function(cl, f, max_iter, tol) {
  .frailty_checks(f)
  g <- .frailty_ids(f)
  n <- length(f$time); p <- ncol(f$X); C <- nlevels(g)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  c0 <- .cox_fit(f$time, f$status, f$X, NULL, "breslow", 50L, 1e-10, NULL, NULL)
  # Cox model with the Breslow baseline: the marginal likelihood at theta = 0
  dk <- table(f$time[f$status == 1L])
  ll0 <- c0$loglik + sum(dk * log(dk)) - sum(f$status)
  r <- .C("SX_cox_frailty", as.integer(n), as.integer(p), as.double(f$time), as.integer(f$status), .rm(f$X),
          as.integer(g) - 1L, as.integer(C), as.double(c(c0$coefficients, log(0.5))), as.integer(max_iter), as.double(tol),
          est = double(p + 1L), var = double((p + 1L)^2), info_d = double(3), info_i = integer(4),
          K = integer(1), tk = double(n), hk = double(n), shape = double(C), rate = double(C), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "frailty model")
  nm <- colnames(f$X)
  theta <- exp(r$est[p + 1L])
  V <- matrix(r$var, p + 1L, p + 1L, byrow = TRUE, dimnames = list(c(nm, "log(theta)"), c(nm, "log(theta)")))
  boundary <- theta < 1e-5
  ll <- r$info_d[2L]
  if (boundary) {
    beta <- c0$coefficients; Vb <- c0$vcov; theta <- 0; tse <- NA_real_; ll <- ll0
  } else {
    beta <- stats::setNames(r$est[seq_len(p)], nm); Vb <- V[seq_len(p), seq_len(p), drop = FALSE]
    tse <- theta * sqrt(V[p + 1L, p + 1L])
  }
  K <- r$K
  FrailtyFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
             coefficients = beta, vcov = Vb, loglik = ll, n = n, nevent = as.integer(sum(f$status)),
             iterations = as.integer(r$info_i[1L]), converged = as.logical(r$info_i[2L]) && (boundary || all(is.finite(diag(V)))),
             flags = as.integer(r$info_i[3L]), time = f$time, status = f$status, X = f$X,
             kind = "cox", base = NULL, theta = theta, theta_se = tse, vcov_full = V, boundary = boundary,
             loglik0 = ll0, lrt = .frailty_lrt(ll, ll0, boundary), cluster = levels(g), cluster_index = as.integer(g),
             cluster_var = .frailty_var_name(f$frailty_var), shape = r$shape, rate = r$rate,
             basehaz = data.frame(time = r$tk[seq_len(K)], hazard = r$hk[seq_len(K)]),
             loglik_null = c0$loglik_null)
}

# parametric baseline: full marginal likelihood
.par_frailty <- function(cl, f, dist, scale, init, max_iter, tol, df, knots, cuts) {
  .frailty_checks(f)
  g <- .frailty_ids(f)
  n <- length(f$time); p <- ncol(f$X); C <- nlevels(g)
  kn <- .par_knots(dist, f$time, f$status, df, knots, cuts)
  b0 <- .par_fit(f$time, f$status, NULL, f$X, dist, scale, NULL, max_iter, tol, kn)
  P0 <- length(b0$theta)
  start <- if (!is.null(init)) as.double(init) else c(unname(b0$theta), log(0.5))
  if (length(start) != P0 + 1L) stop("init must have length ", P0 + 1L, " (the parameters of the model and log(theta))", call. = FALSE)
  r <- .C("SX_par_frailty", as.integer(n), as.integer(p), .par_dists[[dist]], as.integer(scale == "AFT"),
          as.integer(length(kn)), if (length(kn)) as.double(kn) else double(1), as.double(f$time), as.integer(f$status),
          .rm(f$X), as.integer(g) - 1L, as.integer(C), as.double(start), as.integer(max_iter), as.double(tol),
          theta = double(P0 + 1L), var = double((P0 + 1L)^2), info_d = double(3), info_i = integer(4),
          shape = double(C), rate = double(C), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "frailty model")
  nm <- names(b0$theta)
  V <- matrix(r$var, P0 + 1L, P0 + 1L, byrow = TRUE, dimnames = list(c(nm, "log(theta)"), c(nm, "log(theta)")))
  theta <- exp(r$theta[P0 + 1L])
  boundary <- theta < 1e-5
  ll <- r$info_d[2L]; ll0 <- b0$loglik
  if (boundary) {
    est <- b0$theta; Vb <- b0$var; theta <- 0; tse <- NA_real_; ll <- ll0
  } else {
    est <- stats::setNames(r$theta[seq_len(P0)], nm); Vb <- V[seq_len(P0), seq_len(P0), drop = FALSE]
    tse <- theta * sqrt(V[P0 + 1L, P0 + 1L])
  }
  base <- ParSurvFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
                     coefficients = est, vcov = Vb, loglik = ll, n = n, nevent = as.integer(sum(f$status)),
                     iterations = as.integer(r$info_i[1L]), converged = as.logical(r$info_i[2L]),
                     flags = as.integer(r$info_i[3L]), time = f$time, status = f$status, X = f$X,
                     dist = dist, scale = scale, dist_code = .par_dists[[dist]], aft = as.integer(scale == "AFT"),
                     entry = NULL, loglik_null = NA_real_, event_label = "the event", knots = kn)
  FrailtyFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
             coefficients = est, vcov = Vb, loglik = ll, n = n, nevent = as.integer(sum(f$status)),
             iterations = as.integer(r$info_i[1L]),
             converged = as.logical(r$info_i[2L]) && (boundary || all(is.finite(diag(V)))),
             flags = as.integer(r$info_i[3L]), time = f$time, status = f$status, X = f$X,
             kind = "parametric", base = base, theta = theta, theta_se = tse, vcov_full = V, boundary = boundary,
             loglik0 = ll0, lrt = .frailty_lrt(ll, ll0, boundary), cluster = levels(g), cluster_index = as.integer(g),
             cluster_var = .frailty_var_name(f$frailty_var), shape = r$shape, rate = r$rate, basehaz = NULL,
             loglik_null = NA_real_)
}

# ----------------------------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------------------------
.frailty_estimates <- function(x, level) {
  if (x@kind == "parametric") return(.par_estimates(x@base, NA_integer_, level))
  .est_frame(NA_integer_, names(x@coefficients), "beta", "log hazard ratio", x@coefficients, .se_of(x), level, TRUE)
}

.frailty_model_label <- function(x) {
  if (x@kind == "cox") "Cox (Breslow baseline)" else .par_model_label(x@base@dist, x@base@scale)
}

.frailty_lines <- function(x, level = 0.95, detail = FALSE) {
  e <- .frailty_estimates(x, level)
  sz <- tabulate(x@cluster_index)
  aft <- x@kind == "parametric" && x@base@scale == "AFT"
  lab <- if (aft) "TR" else "HR"
  out <- c(paste0("Shared gamma frailty model: ", .frailty_model_label(x)),
           paste0("  Subjects: ", x@n, "   events: ", x@nevent, "   clusters (", x@cluster_var, "): ", length(x@cluster),
                  "   cluster size ", min(sz), " to ", max(sz), ", median ", stats::median(sz)),
           "",
           paste0("Effect of covariates given the frailty: ", if (aft) "time ratio (TR)" else "hazard ratio (HR)"),
           .effect_lines(e[e$role == "beta", , drop = FALSE], lab, level, x@xlevels),
           "",
           .reading_lines(e, x@xlevels, if (aft) "AFT" else "PH", "the event", level),
           "  These ratios compare subjects of the same cluster (conditional effects); for a new cluster",
           "  the population-average effect is closer to 1: predict(fit, newdata, times, frailty = \"marginal\").")
  zq <- .zq(level)
  out <- c(out, "", "Frailty: heterogeneity between clusters")
  if (x@boundary) {
    out <- c(out, "  The frailty variance is at the boundary (0): no heterogeneity between clusters is seen;",
             "  the fit is the model without frailty.")
  } else {
    th <- x@theta; se_l <- x@theta_se / th
    lo <- th * exp(-zq * se_l); hi <- th * exp(zq * se_l)
    a <- 1 / th
    iqr <- stats::qgamma(0.75, a, a) / stats::qgamma(0.25, a, a)
    out <- c(out,
             paste0("  variance theta = ", .num(th, 3L), " (", .pct(level), " CI ", .ci_text(lo, hi, 3L), ")"),
             paste0("  likelihood ratio test of theta = 0: chi-square = ", .num(x@lrt$statistic, 2L), ", p ",
                    { pv <- .pval(x@lrt$p.value); if (startsWith(pv, "<")) sub("<", "< ", pv) else paste("=", pv) },
                    " (half chi-square, Self & Liang 1987)"),
             paste0("  Kendall's tau between two members of a cluster: ", .num(th / (th + 2), 3L)),
             paste0("  A cluster at the 75th percentile of the frailty has ", .num(iqr, 2L),
                    " times the hazard of a cluster at the 25th percentile."))
    out <- c(out, if (x@lrt$p.value < 0.05)
      "  Reading: clusters differ in their risk beyond the covariates; members of a cluster are correlated."
      else "  Reading: no evidence of heterogeneity between clusters.")
  }
  k <- length(x@coefficients) + 1L
  out <- c(out, "", sprintf("Fit: marginal log-likelihood %.2f | AIC %.2f | %d parameters (with theta)", x@loglik, -2 * x@loglik + 2 * k, k),
           sprintf("     without frailty: log-likelihood %.2f", x@loglik0))
  if (x@kind == "cox") out <- c(out, "     ties: Breslow; standard errors from the profile likelihood (theta estimated)")
  else out <- c(out, "     standard errors from the observed information (theta estimated)")
  if (detail && x@kind == "parametric") out <- c(out, "", "Baseline (given frailty 1):", .table_lines(.natural_baseline(x@base, level)))
  if (detail) {
    fr <- .frailty_frame(x, level)
    o <- order(fr$frailty)
    sel <- unique(c(utils::head(o, 3L), utils::tail(o, 3L)))
    out <- c(out, "", "Clusters with the lowest and the highest estimated frailty (posterior mean, 95% interval):",
             .table_lines(data.frame(cluster = fr$cluster[sel], events = fr$events[sel], frailty = .num(fr$frailty[sel], 2L),
                                     interval = .ci_text(fr$lower[sel], fr$upper[sel], 2L), stringsAsFactors = FALSE), 4L))
  }
  c(out, .status_lines(x@converged, x@flags), "", .method_lines(.method_keys(x)))
}

.frailty_frame <- function(x, level = 0.95) {
  a <- (1 - level) / 2
  ev <- as.vector(tapply(x@status, factor(x@cluster_index, seq_along(x@cluster)), sum))
  nn <- tabulate(x@cluster_index, length(x@cluster))
  if (x@boundary)
    return(data.frame(cluster = x@cluster, n = nn, events = ev, frailty = 1, lower = 1, upper = 1, stringsAsFactors = FALSE))
  data.frame(cluster = x@cluster, n = nn, events = ev, frailty = x@shape / x@rate,
             lower = stats::qgamma(a, x@shape, x@rate), upper = stats::qgamma(1 - a, x@shape, x@rate),
             stringsAsFactors = FALSE)
}

S7::method(print, FrailtyFit) <- function(x, level = 0.95, ...) { cat(.frailty_lines(x, level, FALSE), sep = "\n"); invisible(x) }

S7::method(summary, FrailtyFit) <- function(object, level = 0.95, ...) {
  FitSummary(model_call = object@model_call, coefficients = .frailty_estimates(object, level), tables = list(),
             text = .frailty_lines(object, level, TRUE))
}

S7::method(logLik, FrailtyFit) <- function(object, ...) {
  structure(object@loglik, df = length(object@coefficients) + 1L, nobs = object@n, class = "logLik")
}

S7::method(as.data.frame, FrailtyFit) <- function(x, row.names = NULL, optional = FALSE, what = c("coefficients", "frailty"),
                                                  level = 0.95, ...) {
  what <- match.arg(what)
  if (what == "frailty") .frailty_frame(x, level) else .frailty_estimates(x, level)
}

# conditional cumulative hazard (frailty 1) of the rows of X at the times
.frailty_cumhaz <- function(x, X, times) {
  if (x@kind == "parametric") return(.par_predict_X(x@base, X, times, "cumhaz", FALSE, 0.95))
  bh <- x@basehaz
  H0 <- vapply(times, function(t) sum(bh$hazard[bh$time <= t]), numeric(1))
  e <- exp(drop(X %*% x@coefficients))
  out <- outer(e, H0)
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}

S7::method(predict, FrailtyFit) <- function(object, newdata = NULL, times, type = c("survival", "cumhaz"),
                                            frailty = c("marginal", "conditional", "cluster"), ...) {
  type <- match.arg(type); frailty <- match.arg(frailty)
  if (missing(times) || !is.numeric(times)) stop("times (numeric) is required", call. = FALSE)
  X <- .sx_newX(object, newdata)
  H <- .frailty_cumhaz(object, X, times)
  th <- object@theta
  if (frailty == "conditional" || object@boundary) {
    S <- exp(-H)
  } else if (frailty == "marginal") {
    S <- (1 + th * H)^(-1 / th)
  } else {
    if (is.null(newdata) || !object@cluster_var %in% names(newdata))
      stop("frailty = \"cluster\" needs the cluster in a column '", object@cluster_var, "' of newdata", call. = FALSE)
    k <- match(as.character(newdata[[object@cluster_var]]), object@cluster)
    if (anyNA(k)) stop("unknown cluster(s) in newdata: ", paste(unique(newdata[[object@cluster_var]][is.na(k)]), collapse = ", "), call. = FALSE)
    S <- (object@rate[k] / (object@rate[k] + H))^object@shape[k]
  }
  out <- if (type == "survival") S else -log(S)
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}

S7::method(plot, FrailtyFit) <- function(x, which = c("frailty", "survival"), newdata = NULL, level = 0.95, col = NULL, ...) {
  which <- match.arg(which)
  if (which == "frailty") {
    fr <- .frailty_frame(x, level)
    o <- order(fr$frailty); fr <- fr[o, ]
    k <- nrow(fr)
    graphics::plot(fr$frailty, seq_len(k), xlim = range(c(fr$lower, fr$upper)), yaxt = "n", pch = 19, cex = 0.6,
                   col = col %||% .okabe_ito[1L], log = "x", xlab = "Frailty (posterior mean and interval, log scale)",
                   ylab = paste0("Cluster (", x@cluster_var, "), ordered"), main = "Estimated frailty of the clusters", ...)
    graphics::segments(fr$lower, seq_len(k), fr$upper, seq_len(k), col = grDevices::adjustcolor(col %||% .okabe_ito[1L], 0.5))
    graphics::abline(v = 1, lty = 2, col = "grey45")
    if (k <= 40L) graphics::axis(2, at = seq_len(k), labels = fr$cluster, las = 1, cex.axis = 0.6)
    return(invisible(x))
  }
  nd <- if (is.null(newdata)) NULL else newdata
  X <- if (is.null(nd)) .profile(x)$X else .sx_newX(x, nd)
  tg <- seq(0, max(x@time), length.out = 201L)[-1L]
  H <- .frailty_cumhaz(x, X[1L, , drop = FALSE], tg)
  Sc <- exp(-H[1L, ]); Sm <- if (x@boundary) Sc else (1 + x@theta * H[1L, ])^(-1 / x@theta)
  cc <- col %||% .okabe_ito[1:2]
  graphics::plot(c(0, tg), c(1, Sc), type = "l", lwd = 2, col = cc[1L], ylim = c(0, 1), xlab = "Time",
                 ylab = "Survival probability", main = if (is.null(nd)) "Typical covariate values" else "Covariate values given in newdata", ...)
  graphics::lines(c(0, tg), c(1, Sm), lwd = 2, col = cc[2L], lty = 2)
  graphics::legend("topright", legend = c("average cluster (frailty 1)", "new cluster (marginal)"), col = cc, lwd = 2, lty = 1:2, bty = "n")
  invisible(x)
}
