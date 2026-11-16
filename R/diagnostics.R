# Residuals, the test of proportional hazards and plots for choosing a distribution.

# risk-set sums of a Cox model (Breslow): for every event row, the weighted mean x-bar(t) and covariance
# V(t) of the covariates over the subjects still at risk
.cox_risksets <- function(time, status, X, beta, w = NULL) {
  n <- length(time); p <- ncol(X)
  if (is.null(w)) w <- rep(1, n)
  r <- w * exp(drop(X %*% beta))
  o <- order(time)
  ts <- time[o]
  rc <- function(v) rev(cumsum(rev(v)))                 # sum over j >= i in time order
  S0 <- rc(r[o])
  S1 <- apply(X[o, , drop = FALSE] * r[o], 2L, rc)
  if (!is.matrix(S1)) S1 <- matrix(S1, ncol = p)
  ev <- which(status == 1L)
  ev <- ev[order(time[ev])]
  first <- match(time[ev], ts)                           # first position of each event time: everybody with time >= t
  s0 <- S0[first]; xbar <- S1[first, , drop = FALSE] / s0
  # second moments
  S2 <- array(0, c(length(ev), p, p))
  for (a in seq_len(p)) for (b in seq_len(a)) {
    v <- rc(X[o, a] * X[o, b] * r[o])[first] / s0 - xbar[, a] * xbar[, b]
    S2[, a, b] <- v; S2[, b, a] <- v
  }
  list(rows = ev, time = time[ev], xbar = xbar, V = S2, w = w[ev])
}

.resid_cox_core <- function(time, status, X, beta, vcov, ties, type, w = NULL, baseline = NULL) {
  if (type %in% c("martingale", "deviance", "coxsnell")) {
    H <- .step_at(baseline$time, baseline$cumhaz, time, 0) * exp(drop(X %*% beta))
    m <- status - H
    out <- switch(type, martingale = m, coxsnell = H,
                  deviance = sign(m) * sqrt(-2 * (m + ifelse(status == 1L, log(pmax(status - m, 1e-300)), 0))))
    return(out)
  }
  rs <- .cox_risksets(time, status, X, beta, w)
  sch <- X[rs$rows, , drop = FALSE] - rs$xbar
  dimnames(sch) <- list(format(rs$time, trim = TRUE), colnames(X))
  if (type == "schoenfeld") return(structure(sch, time = rs$time))
  nd <- length(rs$rows)
  sc <- sch %*% vcov * nd + matrix(beta, nd, length(beta), byrow = TRUE)
  dimnames(sc) <- dimnames(sch)
  structure(sc, time = rs$time)
}

.par_resid <- function(m, type) {
  Xn <- m@X; n <- m@n
  Ht <- vapply(seq_len(n), function(i) .par_predict_X(m, Xn[i, , drop = FALSE], m@time[i], "cumhaz", FALSE, 0.95)[1L, 1L], numeric(1))
  if (!is.null(m@entry)) Ht <- Ht - vapply(seq_len(n), function(i)
    if (m@entry[i] > 0) .par_predict_X(m, Xn[i, , drop = FALSE], m@entry[i], "cumhaz", FALSE, 0.95)[1L, 1L] else 0, numeric(1))
  mres <- m@status - Ht
  switch(type, coxsnell = Ht, martingale = mres,
         deviance = sign(mres) * sqrt(-2 * (mres + ifelse(m@status == 1L, log(pmax(m@status - mres, 1e-300)), 0))),
         stop("parametric models have martingale, deviance and Cox-Snell residuals", call. = FALSE))
}

#' Residuals of hazard models
#'
#' Martingale, deviance, Cox-Snell and (scaled) Schoenfeld residuals of Cox models, cause-specific Cox
#' models and parametric models.
#'
#' @details With \eqn{\hat H_i} the fitted cumulative hazard of subject \eqn{i} at its time (minus that
#' at entry), the Cox-Snell residual is \eqn{\hat H_i}, the martingale residual \eqn{d_i - \hat H_i}
#' and the deviance residual \eqn{\mathrm{sign}(M_i)\sqrt{-2\{M_i + d_i \log(d_i - M_i)\}}}. Schoenfeld
#' residuals are defined at every event: the covariate vector minus its risk-weighted mean over the
#' subjects at risk (Breslow weights); the scaled residuals are \eqn{\hat\beta + D\,V r}, with
#' \eqn{D} the number of events and \eqn{V} the covariance of \eqn{\hat\beta}, as in \pkg{survival},
#' and their smoothed trend against time estimates \eqn{\beta(t)} (Grambsch and Therneau 1994).
#' For cause-specific models a list with one element per cause is returned.
#'
#' @param object a model fitted by \code{\link{cox}}, \code{\link{csc}} or \code{\link{parsurv}}.
#' @param type \code{"martingale"}, \code{"deviance"}, \code{"coxsnell"}, \code{"schoenfeld"} or
#'   \code{"scaledsch"} (Cox models only for the last two).
#' @param ... not used.
#' @return a vector with one value per subject, or for Schoenfeld residuals a matrix with one row per
#'   event (attribute \code{"time"}); for cause-specific models a list of these.
#' @references
#' Cox DR, Snell EJ (1968) A general definition of residuals. J R Stat Soc B 30:248-275.
#' Schoenfeld D (1982) Partial residuals for the proportional hazards regression model. Biometrika
#'   69:239-241.
#' Therneau TM, Grambsch PM, Fleming TR (1990) Martingale-based residuals for survival models.
#'   Biometrika 77:147-160.
#' Grambsch PM, Therneau TM (1994) Proportional hazards tests and diagnostics based on weighted
#'   residuals. Biometrika 81:515-526.
#' @name residuals.SurvX
NULL

.resid_types <- c("martingale", "deviance", "coxsnell", "schoenfeld", "scaledsch")

S7::method(residuals, CoxFit) <- function(object, type = c("martingale", "deviance", "coxsnell", "schoenfeld", "scaledsch"), ...) {
  type <- match.arg(type)
  .resid_cox_core(object@time, object@status, object@X, object@coefficients, object@vcov, object@ties, type,
                  object@weights, if (type %in% c("martingale", "deviance", "coxsnell")) .cox_baseline(object) else NULL)
}

.csc_baseline <- function(object, k) {
  m <- object@models[[k]]; st <- as.integer(object@status == k)
  n <- object@n; p <- ncol(object@X)
  r <- .C("SX_cox_baseline", as.integer(n), as.integer(p), object@time, st, .rm(object@X), 0L, double(n),
          as.integer(object@ties == "efron"), as.double(m$coefficients), time = double(n), cumhaz = double(n),
          m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "baseline hazard")
  k2 <- seq_len(r$m)
  data.frame(time = r$time[k2], cumhaz = r$cumhaz[k2])
}

S7::method(residuals, CSCFit) <- function(object, type = c("martingale", "deviance", "coxsnell", "schoenfeld", "scaledsch"), ...) {
  type <- match.arg(type)
  out <- lapply(seq_len(object@K), function(k) {
    m <- object@models[[k]]
    .resid_cox_core(object@time, as.integer(object@status == k), object@X, m$coefficients, m$vcov, object@ties, type, NULL,
                    if (type %in% c("martingale", "deviance", "coxsnell")) .csc_baseline(object, k) else NULL)
  })
  names(out) <- paste0("cause", seq_len(object@K))
  out
}

S7::method(residuals, ParSurvFit) <- function(object, type = c("martingale", "deviance", "coxsnell"), ...) {
  .par_resid(object, match.arg(type))
}

S7::method(residuals, ParCSCFit) <- function(object, type = c("martingale", "deviance", "coxsnell"), ...) {
  type <- match.arg(type)
  out <- lapply(object@models, function(m) .par_resid(m, type))
  names(out) <- paste0("cause", seq_len(object@K))
  out
}

# ----------------------------------------------------------------------------------------------
# test of proportional hazards
# ----------------------------------------------------------------------------------------------
.time_transform <- function(time, status, t_event, transform) {
  switch(transform,
         identity = t_event,
         log = log(t_event),
         rank = rank(t_event),
         km = {                                                # 1 - Kaplan-Meier just before each event time
           inc <- .na_increments(time, status, 1L)
           S <- cumprod(1 - inc$d / inc$Y)
           1 - c(1, S)[match(t_event, inc$time)]
         })
}

.zph_core <- function(time, status, X, beta, vcov, transform, w = NULL) {
  rs <- .cox_risksets(time, status, X, beta, w)
  sch <- X[rs$rows, , drop = FALSE] - rs$xbar
  g <- .time_transform(time, status, rs$time, transform)
  g <- g - mean(g)
  p <- ncol(X); nd <- length(rs$rows)
  U <- colSums(sch * g)
  Ibb <- apply(rs$V, c(2L, 3L), sum)
  Ibg <- apply(rs$V * g, c(2L, 3L), sum)
  Igg <- apply(rs$V * g^2, c(2L, 3L), sum)
  Ibb_inv <- solve(Ibb)
  A <- Igg - t(Ibg) %*% Ibb_inv %*% Ibg                 # efficient information for the time-varying terms
  chi <- vapply(seq_len(p), function(k) U[k]^2 / A[k, k], numeric(1))
  glob <- drop(t(U) %*% solve(A, U))
  tab <- data.frame(chisq = c(chi, glob), df = c(rep(1L, p), p),
                    p = stats::pchisq(c(chi, glob), c(rep(1L, p), p), lower.tail = FALSE),
                    row.names = c(colnames(X), "GLOBAL"))
  scaled <- sch %*% vcov * nd + matrix(beta, nd, p, byrow = TRUE)
  colnames(scaled) <- colnames(X)
  list(table = tab, x = g + mean(.time_transform(time, status, rs$time, transform)), y = scaled, time = rs$time,
       transform = transform, beta = beta, vcov = vcov)
}

#' Test of proportional hazards
#'
#' Score test of the proportional hazards assumption of a Cox model (or of every cause of a
#' cause-specific Cox model), for each covariate and globally, with the plot of the scaled Schoenfeld
#' residuals against time.
#'
#' @details For each covariate the test adds a term \eqn{\gamma_j x_j g(t)} to the model and tests
#' \eqn{\gamma = 0} by the score test at the fitted \eqn{\hat\beta}: the score is
#' \eqn{\sum_i g(t_i) r_{ij}} with the Schoenfeld residuals \eqn{r}, and its variance is the
#' efficient information for \eqn{\gamma} given \eqn{\beta}. This is the test of Grambsch and Therneau
#' (1994) in its exact score form, as computed by \code{survival::cox.zph} since version 3.0; risk-set
#' moments use Breslow weights, so results agree with \pkg{survival} exactly when there are no tied
#' event times. The plot shows the scaled Schoenfeld residuals, whose smoothed trend estimates the
#' time-varying effect \eqn{\beta_j(t)}: a flat line at \eqn{\hat\beta_j} supports proportional
#' hazards.
#'
#' @param fit a model fitted by \code{\link{cox}} or by \code{\link{csc}} with a Cox baseline.
#' @param transform the function \eqn{g(t)} of time: \code{"km"} (default, one minus the Kaplan-Meier
#'   estimate just before \eqn{t}), \code{"rank"}, \code{"identity"} or \code{"log"}.
#' @return an object of class \code{"ph_test"} (a list per cause for cause-specific models) with the
#'   table of tests; \code{print} and \code{plot(x, var, df = 4)} methods.
#' @references Grambsch PM, Therneau TM (1994) Proportional hazards tests and diagnostics based on
#'   weighted residuals. Biometrika 81:515-526.
#'   Therneau TM, Grambsch PM (2000) Modeling Survival Data: Extending the Cox Model. Springer,
#'   Chapter 6.
#' @examples
#' set.seed(4)
#' n <- 300
#' x <- rbinom(n, 1, 0.5)
#' t <- rweibull(n, ifelse(x == 1, 0.7, 1.5), 3); cens <- runif(n, 0, 6)  # non-proportional hazards
#' d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = rnorm(n))
#' pt <- ph_test(cox(Event(time, status) ~ x + z, d))
#' pt
#' plot(pt, "x")
#' @export
ph_test <- function(fit, transform = c("km", "rank", "identity", "log")) {
  transform <- match.arg(transform)
  if (S7::S7_inherits(fit, CoxFit)) {
    r <- .zph_core(fit@time, fit@status, fit@X, fit@coefficients, fit@vcov, transform, fit@weights)
    r$label <- "the event"
    return(structure(list(r), class = "ph_test"))
  }
  if (S7::S7_inherits(fit, CSCFit)) {
    out <- lapply(seq_len(fit@K), function(k) {
      m <- fit@models[[k]]
      r <- .zph_core(fit@time, as.integer(fit@status == k), fit@X, m$coefficients, m$vcov, transform)
      r$label <- paste("cause", k); r
    })
    return(structure(out, class = "ph_test"))
  }
  stop("ph_test() takes a model fitted by cox() or csc() with a Cox baseline", call. = FALSE)
}

#' @export
print.ph_test <- function(x, ...) {
  cat("Test of proportional hazards (score test of Grambsch and Therneau; time scale: ", x[[1L]]$transform, ")\n", sep = "")
  for (r in x) {
    tab <- r$table
    df <- data.frame(term = rownames(tab), chisq = .num(tab$chisq, 2L), df = tab$df, p = .pval(tab$p), stringsAsFactors = FALSE)
    names(df) <- c("", "chi-square", "df", "p")
    cat("\n", if (length(x) > 1L) paste0("Cause ", sub("cause ", "", r$label), "\n") else "", sep = "")
    cat(.table_lines(df), sep = "\n")
    bad <- rownames(tab)[-nrow(tab)][tab$p[-nrow(tab)] < 0.05]
    cat(if (length(bad)) paste0("Reading: the effect of ", paste(bad, collapse = ", "),
                                " changes over time (p < 0.05); proportional hazards are doubtful for ",
                                if (length(bad) > 1L) "these covariates." else "this covariate.")
        else "Reading: no evidence against proportional hazards (all p >= 0.05).", "\n")
  }
  cat("\nplot(x, var) shows the scaled Schoenfeld residuals: a flat trend supports proportional hazards.\n")
  cat("\n", .method_lines(.method_keys(x)), sep = "\n")
  invisible(x)
}

#' @export
plot.ph_test <- function(x, var = NULL, cause = 1L, df = 4L, level = 0.95, col = "#0072B2", xlab = NULL, ylab = NULL,
                         main = NULL, cex = 0.85, ...) {
  r <- x[[min(cause, length(x))]]
  vars <- if (is.null(var)) colnames(r$y) else var
  op <- graphics::par(mfrow = c(1L, length(vars)), cex = cex, mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  for (v in vars) {
    y <- r$y[, v]; xx <- r$x
    ok <- is.finite(xx) & is.finite(y)
    graphics::plot(xx[ok], y[ok], pch = 16, cex = 0.5, col = grDevices::adjustcolor("grey30", 0.5),
                   xlab = if (is.null(xlab)) paste0("Time (", r$transform, " scale)") else xlab,
                   ylab = if (is.null(ylab)) paste0("Beta(t) for ", v) else ylab,
                   main = if (is.null(main)) (if (length(x) > 1L) r$label else "") else main, ...)
    # natural spline regression of the scaled residuals on g(t), with pointwise limits
    k <- min(as.integer(df), length(unique(xx[ok])) - 1L)
    if (k >= 1L) {
      fit <- stats::lm(y ~ splines::ns(x, df = k), data = data.frame(x = xx[ok], y = y[ok]))
      grid <- seq(min(xx[ok]), max(xx[ok]), length.out = 101L)
      pr <- stats::predict(fit, newdata = data.frame(x = grid), se.fit = TRUE)
      zq <- .zq(level)
      .band(grid, pr$fit - zq * pr$se.fit, pr$fit + zq * pr$se.fit, col)
      graphics::lines(grid, pr$fit, col = col, lwd = 2)
    }
    graphics::abline(h = r$beta[[v]], lty = 2, col = "grey40")
  }
  invisible(x)
}

# ----------------------------------------------------------------------------------------------
# choosing a distribution
# ----------------------------------------------------------------------------------------------

#' Plots for choosing a distribution and checking a fitted model
#'
#' With a formula: the classic linearising plots of the non-parametric estimate (Nelson-Aalen
#' cumulative hazard, by group): straight lines indicate the distribution of the panel. With a fitted
#' parametric or Cox model: the Cox-Snell residual plot.
#'
#' @details Panels for a formula, with \eqn{\hat H} the Nelson-Aalen estimate and
#' \eqn{\hat S = e^{-\hat H}}:
#' \describe{
#'   \item{exponential}{\eqn{\hat H(t)} against \eqn{t}: a straight line through the origin.}
#'   \item{Weibull}{\eqn{\log \hat H(t)} against \eqn{\log t}: a straight line (slope = shape;
#'     parallel lines across groups indicate proportional hazards).}
#'   \item{log-normal}{\eqn{\Phi^{-1}(1 - \hat S(t))} against \eqn{\log t}: a straight line.}
#'   \item{log-logistic}{\eqn{\log\{(1 - \hat S)/\hat S\}} against \eqn{\log t}: a straight line
#'     (parallel lines indicate proportional odds).}
#' }
#' For a fitted model, the Nelson-Aalen cumulative hazard of the Cox-Snell residuals \eqn{\hat H_i}
#' (with the event indicators) is plotted against the residuals; if the model fits, the points lie
#' on the line of identity, since the residuals then behave like a censored sample from the unit
#' exponential distribution. With competing risks, \code{cause} selects the cause-specific hazard.
#'
#' @param x a formula \code{Event(time, status) ~ 1} or \code{~ group}, or a model fitted by
#'   \code{\link{parsurv}}, \code{\link{cox}} or \code{\link{csc}} (parametric or Cox).
#' @param data a data frame (with a formula).
#' @param cause the cause for competing risks data.
#' @param col colours of the groups.
#' @param cex character expansion.
#' @param ... passed to \code{plot.default}.
#' @return invisibly, the plotted data.
#' @references Klein JP, Moeschberger ML (2003) Survival Analysis, 2nd ed. Springer, Sections 12.2-12.5.
#'   Collett D (2015) Modelling Survival Data in Medical Research, 3rd ed. CRC Press, Chapter 5.
#'   Cox DR, Snell EJ (1968) J R Stat Soc B 30:248-275.
#' @examples
#' set.seed(5)
#' d <- data.frame(time = rweibull(300, 1.5, 3), status = 1, arm = rep(c("A", "B"), 150))
#' dist_plot(Event(time, status) ~ arm, d)
#' dist_plot(parsurv(Event(time, status) ~ arm, d))
#' @export
dist_plot <- function(x, data = NULL, cause = 1L, col = NULL, cex = 0.85, ...) {
  if (inherits(x, "formula")) return(.dist_plot_formula(x, data, cause, col, cex, ...))
  .coxsnell_plot(x, cause, col, cex, ...)
}

.dist_plot_formula <- function(formula, data, cause, col, cex, ...) {
  f <- .sx_groups(formula, data, what = "dist_plot()", entry_ok = TRUE)
  groups <- levels(f$group)
  col <- if (is.null(col)) rep_len(.okabe_ito, length(groups)) else rep_len(col, length(groups))
  cur <- lapply(groups, function(g) {
    sel <- f$group == g
    inc <- .na_increments(f$time[sel], f$status[sel], cause, if (is.null(f$entry)) NULL else f$entry[sel])
    if (is.null(inc)) return(NULL)
    H <- cumsum(inc$d / inc$Y)
    data.frame(group = g, time = inc$time, H = H, S = exp(-H))
  })
  names(cur) <- groups
  panels <- list(
    exponential = list(x = function(d) d$time, y = function(d) d$H, xl = "Time", yl = "Cumulative hazard H(t)"),
    Weibull = list(x = function(d) log(d$time), y = function(d) log(d$H), xl = "log time", yl = "log H(t)"),
    "log-normal" = list(x = function(d) log(d$time), y = function(d) stats::qnorm(1 - d$S), xl = "log time", yl = "probit(1 - S(t))"),
    "log-logistic" = list(x = function(d) log(d$time), y = function(d) log((1 - d$S) / d$S), xl = "log time", yl = "logit(1 - S(t))"))
  op <- graphics::par(mfrow = c(2L, 2L), cex = cex, mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  for (nm in names(panels)) {
    pn <- panels[[nm]]
    xs <- unlist(lapply(cur, function(d) if (is.null(d)) NULL else pn$x(d))); ys <- unlist(lapply(cur, function(d) if (is.null(d)) NULL else pn$y(d)))
    ok <- is.finite(xs) & is.finite(ys)
    graphics::plot(range(xs[ok]), range(ys[ok]), type = "n", xlab = pn$xl, ylab = pn$yl,
                   main = paste0(nm, ": straight line?"), ...)
    for (i in seq_along(cur)) if (!is.null(cur[[i]])) {
      xv <- pn$x(cur[[i]]); yv <- pn$y(cur[[i]]); k <- is.finite(xv) & is.finite(yv)
      graphics::lines(xv[k], yv[k], type = "s", col = col[i], lwd = 1.5)
    }
    if (length(groups) > 1L) graphics::legend("topleft", legend = groups, col = col, lwd = 1.5, bty = "n")
  }
  invisible(do.call(rbind, cur))
}

.coxsnell_plot <- function(x, cause, col, cex, ...) {
  if (S7::S7_inherits(x, ParCSCFit) || S7::S7_inherits(x, CSCFit)) {
    r <- residuals(x, type = "coxsnell")[[cause]]; st <- as.integer(x@status == cause)
    entry <- if (S7::S7_inherits(x, ParCSCFit)) x@entry else NULL
  } else if (S7::S7_inherits(x, ParSurvFit) || S7::S7_inherits(x, CoxFit)) {
    r <- residuals(x, type = "coxsnell"); st <- x@status
  } else stop("dist_plot() takes a formula or a model fitted by parsurv(), cox() or csc()", call. = FALSE)
  inc <- .na_increments(r, st, 1L)
  H <- cumsum(inc$d / inc$Y)
  col <- if (is.null(col)) .okabe_ito[1L] else col[1L]
  op <- graphics::par(cex = cex, mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  lim <- c(0, max(c(inc$time, H)))
  graphics::plot(inc$time, H, type = "s", col = col, lwd = 2, xlim = lim, ylim = lim,
                 xlab = "Cox-Snell residual", ylab = "Cumulative hazard of the residuals",
                 main = "Cox-Snell residuals: close to the line if the model fits", ...)
  graphics::abline(0, 1, lty = 2, col = "grey40")
  invisible(data.frame(residual = inc$time, cumhaz = H))
}
