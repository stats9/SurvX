# Residuals, the test of proportional hazards and plots for choosing a distribution.

# risk-set sums of a Cox model (Breslow): for every event row, the weighted mean x-bar(t) and covariance
# V(t) of the covariates over the subjects still at risk
.cox_risksets <- function(time, status, X, beta, w = NULL, entry = NULL, strata = NULL) {
  if (!is.null(entry) || !is.null(strata)) return(.cox_risksets_general(time, status, X, beta, w, entry, strata))
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

# the same with (entry, time] intervals and strata: risk set of an event at t in stratum s is
# {j in s: entry_j < t <= time_j}; O(n x events), used only when entry or strata are present
.cox_risksets_general <- function(time, status, X, beta, w, entry, strata) {
  n <- length(time); p <- ncol(X)
  if (is.null(w)) w <- rep(1, n)
  if (is.null(entry)) entry <- rep(-Inf, n)
  if (is.null(strata)) strata <- integer(n)
  r <- w * exp(drop(X %*% beta) - max(drop(X %*% beta)))
  ev <- which(status == 1L)
  ev <- ev[order(strata[ev], time[ev])]
  xbar <- matrix(0, length(ev), p); V <- array(0, c(length(ev), p, p))
  for (q in seq_along(ev)) {
    k <- ev[q]; t <- time[k]
    at <- which(strata == strata[k] & entry < t & time >= t)
    ra <- r[at]; Xa <- X[at, , drop = FALSE]; s0 <- sum(ra)
    m <- colSums(Xa * ra) / s0
    xbar[q, ] <- m
    V[q, , ] <- crossprod(Xa * sqrt(ra)) / s0 - tcrossprod(m)
  }
  list(rows = ev, time = time[ev], xbar = xbar, V = V, w = w[ev])
}

# martingale residuals and Schoenfeld means from the C++ core (Efron or Breslow form, as survival)
.cox_resid_c <- function(time, status, X, w, ties, beta, entry, strata) {
  n <- length(time); p <- ncol(X)
  r <- do.call(.C, c(list("SX_cox_resid"), .cox_args(time, status, X, w, entry, strata),
                     list(as.integer(ties == "efron"), as.double(beta), mart = double(n), xbar = double(n * p), rc = integer(1),
                          NAOK = TRUE, PACKAGE = "SurvX")))
  .sx_rc(r$rc, "residuals")
  list(mart = r$mart, xbar = matrix(r$xbar, n, p, byrow = TRUE))
}

.resid_cox_core <- function(time, status, X, beta, vcov, ties, type, w = NULL, baseline = NULL, entry = NULL, strata = NULL) {
  rc <- .cox_resid_c(time, status, X, w, ties, beta, entry, strata)
  if (type %in% c("martingale", "deviance", "coxsnell")) {
    m <- rc$mart
    return(switch(type, martingale = m, coxsnell = status - m,
                  deviance = sign(m) * sqrt(-2 * (m + ifelse(status == 1L, log(pmax(status - m, 1e-300)), 0)))))
  }
  ev <- which(status == 1L)
  ev <- ev[if (is.null(strata)) order(time[ev]) else order(strata[ev], time[ev])]
  sch <- X[ev, , drop = FALSE] - rc$xbar[ev, , drop = FALSE]
  dimnames(sch) <- list(format(time[ev], trim = TRUE), colnames(X))
  if (type == "schoenfeld") return(structure(sch, time = time[ev]))
  nd <- length(ev)
  sc <- sch %*% vcov * nd + matrix(beta, nd, length(beta), byrow = TRUE)
  dimnames(sc) <- dimnames(sch)
  structure(sc, time = time[ev])
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
#' @description
#' Martingale, deviance, Cox-Snell and (scaled) Schoenfeld residuals of Cox models, cause-specific
#' models, parametric models and cure models: \code{residuals(fit, type)}. They are used to check
#' the functional form of covariates, to find poorly fitted subjects, to check the overall fit and
#' the proportional hazards assumption.
#' @details
#' With \eqn{\hat H_i} the fitted cumulative hazard of subject \eqn{i} at its time (minus that at
#' entry), the Cox-Snell residual is \eqn{\hat H_i}, the martingale residual
#' \eqn{M_i = d_i - \hat H_i} and the deviance residual
#' \eqn{\mathrm{sign}(M_i)\sqrt{-2\{M_i + d_i \log(d_i - M_i)\}}}. Schoenfeld residuals are defined at
#' every event: the covariate vector minus its risk-weighted mean over the subjects at risk. For Cox
#' models the cumulative hazard and the means take the form of the ties method of the fit: with Efron
#' ties a subject who fails at a tied time contributes the Efron-weighted hazard increments, and the
#' mean is averaged over the Efron terms, exactly as in \pkg{survival} (Therneau and Grambsch 2000,
#' Section 3.3). The martingale residuals of a Cox model then sum to zero. The scaled residuals are \eqn{\hat\beta + D\,V r}, with \eqn{D} the number of events and
#' \eqn{V} the covariance matrix of \eqn{\hat\beta}, as in \pkg{survival}; their smoothed trend
#' against time estimates \eqn{\beta(t)} (Grambsch and Therneau 1994), see \code{\link{ph_test}}.
#' @section Usage:
#' \code{residuals(object, type = c("martingale", "deviance", "coxsnell", "schoenfeld", "scaledsch"))}
#' @section Arguments:
#' \describe{
#' \item{\code{object}}{a model fitted by \code{\link{cox}}, \code{\link{csc}} (Cox or parametric
#'   baseline), \code{\link{parsurv}} or \code{\link{cure}}.}
#' \item{\code{type}}{\code{"martingale"} (default), \code{"deviance"}, \code{"coxsnell"},
#'   \code{"schoenfeld"} or \code{"scaledsch"}; the last two for Cox models only.}
#' }
#' @section What each residual is for:
#' \describe{
#'   \item{Martingale}{observed minus expected number of events, between \eqn{-\infty} and 1.
#'     Plotted against a numeric covariate that is \emph{not} in the model, a smooth curve shows the
#'     functional form it should have (a straight line: enter it linearly; a threshold or a curve:
#'     categorise, transform or use a spline) (Therneau, Grambsch and Fleming 1990).}
#'   \item{Deviance}{a symmetrised martingale residual; values beyond about \eqn{\pm 2.5} point to
#'     subjects poorly predicted by the model (an early death of a low-risk patient, for example).}
#'   \item{Cox-Snell}{if the model is correct they behave like a censored sample from the unit
#'     exponential distribution: \code{\link{dist_plot}(fit)} plots their cumulative hazard, which
#'     should follow the identity line (Cox and Snell 1968).}
#'   \item{Schoenfeld and scaled Schoenfeld}{one row per event; their trend over time estimates the
#'     time-varying effect \eqn{\beta(t)} and is the basis of \code{\link{ph_test}}.}
#' }
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}), martingale, deviance and scaled Schoenfeld residuals of Cox models agree with
#' \code{survival::residuals.coxph} to \eqn{10^{-8}}, with Efron and Breslow ties, strata and
#' weights; Efron martingale residuals sum to zero to \eqn{2 \times 10^{-13}}.
#' @return
#' a vector with one value per subject, or for Schoenfeld residuals a matrix with one row per event in
#'   time order (attribute \code{"time"}); for cause-specific models a list with one element per cause.
#' @references
#' Cox DR, Snell EJ (1968) A general definition of residuals. Journal of the Royal Statistical Society,
#'   Series B 30:248-275.
#'
#' Schoenfeld D (1982) Partial residuals for the proportional hazards regression model. Biometrika
#'   69:239-241.
#'
#' Therneau TM, Grambsch PM, Fleming TR (1990) Martingale-based residuals for survival models.
#'   Biometrika 77:147-160.
#'
#' Grambsch PM, Therneau TM (1994) Proportional hazards tests and diagnostics based on weighted
#'   residuals. Biometrika 81:515-526.
#' @seealso
#' \code{\link{ph_test}}, \code{\link{dist_plot}}, \code{\link{cox}}, \code{\link{parsurv}}.
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(age = runif(n, 40, 80), x = rbinom(n, 1, 0.5))
#' t <- rexp(n, 0.05 * exp(0.5 * d$x + 0.002 * (d$age - 60)^2)); cens <- runif(n, 0, 20)
#' d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
#'
#' ## functional form of age: martingale residuals of a model without age
#' f0 <- cox(Event(time, status) ~ x, d)
#' plot(d$age, residuals(f0), xlab = "age", ylab = "martingale residual")
#' lines(stats::lowess(d$age, residuals(f0)), lwd = 2)      # U-shaped: age is not linear
#'
#' ## poorly fitted subjects
#' f1 <- cox(Event(time, status) ~ x + splines::ns(age, 3), d)
#' which(abs(residuals(f1, "deviance")) > 2.5)
#'
#' ## overall fit: Cox-Snell residuals
#' dist_plot(f1)
#' @name residuals.SurvX
NULL

.resid_types <- c("martingale", "deviance", "coxsnell", "schoenfeld", "scaledsch")

S7::method(residuals, CoxFit) <- function(object, type = c("martingale", "deviance", "coxsnell", "schoenfeld", "scaledsch"), ...) {
  .no_tt(object, "residuals()")
  type <- match.arg(type)
  .resid_cox_core(object@time, object@status, object@X, object@coefficients, object@vcov_model %||% object@vcov, object@ties, type,
                  object@weights, if (type %in% c("martingale", "deviance", "coxsnell")) .cox_baseline(object) else NULL,
                  object@entry, object@strata)
}

.csc_baseline <- function(object, k) {
  .cox_baseline_raw(object@models[[k]]$coefficients, object@time, as.integer(object@status == k), object@X, NULL,
                    object@ties, object@entry, object@strata)
}

S7::method(residuals, CSCFit) <- function(object, type = c("martingale", "deviance", "coxsnell", "schoenfeld", "scaledsch"), ...) {
  type <- match.arg(type)
  out <- lapply(seq_len(object@K), function(k) {
    m <- object@models[[k]]
    .resid_cox_core(object@time, as.integer(object@status == k), object@X, m$coefficients, m$vcov_model %||% m$vcov,
                    object@ties, type, NULL, if (type %in% c("martingale", "deviance", "coxsnell")) .csc_baseline(object, k) else NULL,
                    object@entry, object@strata)
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
.time_transform <- function(time, status, t_event, transform, entry = NULL, rows = NULL) {
  switch(transform,
         identity = t_event,
         log = log(t_event),
         # rank of the event time among all observed (stop) times, censored included, as in
         # survival::cox.zph (version 3.0 and later)
         rank = if (is.null(rows)) rank(t_event) else rank(time)[rows],
         km = {                                                # 1 - Kaplan-Meier just before each event time
           inc <- .na_increments(time, status, 1L, entry)
           S <- cumprod(1 - inc$d / inc$Y)
           1 - c(1, S)[match(t_event, inc$time)]
         })
}

.zph_core <- function(time, status, X, beta, vcov, transform, w = NULL, entry = NULL, strata = NULL) {
  rs <- .cox_risksets(time, status, X, beta, w, entry, strata)
  sch <- X[rs$rows, , drop = FALSE] - rs$xbar
  g <- .time_transform(time, status, rs$time, transform, entry, rs$rows)
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
  list(table = tab, x = g + mean(.time_transform(time, status, rs$time, transform, entry, rs$rows)), y = scaled, time = rs$time,
       transform = transform, beta = beta, vcov = vcov)
}

#' Test and plot of the proportional hazards assumption
#'
#' @description
#' Tests whether the effect of each covariate of a Cox model is constant over time (proportional
#' hazards), for every covariate and globally, and plots the estimated time-varying effect. For a
#' cause-specific Cox model (\code{\link{csc}}) every cause is tested.
#'
#' \strong{Use it} after every \code{\link{cox}} or Cox-baseline \code{\link{csc}} fit: the hazard
#' ratios are averages over time when the assumption fails.
#' @param fit a model fitted by \code{\link{cox}} or by \code{\link{csc}} with a Cox baseline.
#' @param transform the function \eqn{g(t)} of time against which a change of the effect is sought:
#'   \code{"km"} (default, one minus the Kaplan-Meier estimate just before \eqn{t}; robust to a few
#'   long times), \code{"rank"} (rank of \eqn{t} among all observed times, censored included, as in
#'   \code{survival::cox.zph}), \code{"identity"} (time itself) or \code{"log"} (log time; sensitive
#'   to early changes).
#' @details
#' For each covariate the test adds a term \eqn{\gamma_j x_j g(t)} to the model and tests
#' \eqn{\gamma = 0} by the score test at the fitted \eqn{\hat\beta}: the score is
#' \eqn{\sum_i g(t_i) r_{ij}} with the Schoenfeld residuals \eqn{r}, and its variance is the efficient
#' information for \eqn{\gamma} given \eqn{\beta}. The global test uses all covariates together. This
#' is the test of Grambsch and Therneau (1994) in its exact score form, as computed by
#' \code{survival::cox.zph} since version 3.0. Risk-set moments use Breslow weights, so the results
#' agree with \pkg{survival} exactly when no event times are tied.
#'
#' The plot shows the scaled Schoenfeld residuals against \eqn{g(t)} with a natural spline smooth and
#' its pointwise confidence band. The smooth estimates the time-varying effect \eqn{\beta_j(t)}; a flat
#' line at \eqn{\hat\beta_j} (dashed) supports proportional hazards.
#' @section Plot:
#' \code{plot(x, var = NULL, cause = 1, df = 4, level = 0.95, col = "#0072B2", xlab = NULL,
#' ylab = NULL, main = NULL, cex = 0.85, ...)}, with:
#' \describe{
#' \item{\code{x}}{the result of \code{ph_test()}.}
#' \item{\code{var}}{covariates to plot; default all.}
#' \item{\code{cause}}{the cause to plot for cause-specific models.}
#' \item{\code{df}}{degrees of freedom of the natural spline that smooths the scaled residuals.}
#' \item{\code{level}}{confidence level of the pointwise band of the smooth.}
#' \item{\code{col}}{colour of the smooth.}
#' \item{\code{xlab, ylab, main}}{axis labels and title.}
#' \item{\code{cex}}{character expansion.}
#' \item{\code{...}}{passed to \code{plot.default}.}
#' }
#' @section Interpreting the result:
#' \itemize{
#'   \item A small p-value for a covariate means that its hazard ratio changes over time. The plot
#'     shows how: a rising smooth means the effect grows (the HR moves up), a falling one that it
#'     fades. The global test asks whether any covariate violates the assumption.
#'   \item With large data sets, tiny and clinically irrelevant departures become significant: judge
#'     the size of the change on the plot, not only the p-value.
#'   \item What to do when the assumption fails for a covariate: put it in \code{strata()} if its
#'     effect is not of interest; let its effect change with time with \code{tt()} in
#'     \code{\link{cox}}, or let \code{\link{time_effects}} compare forms of the change and refit the
#'     model with them; or report hazard ratios for separate periods of follow-up. Restricted mean
#'     survival time (\code{\link{rmst}}) needs no proportional hazards.
#' }
#' A sentence for a report: \emph{"The proportional hazards assumption was checked with the
#' Grambsch-Therneau test (global p = 0.41) and plots of scaled Schoenfeld residuals."}
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}), the score statistics agree with \code{survival::cox.zph} for all four
#' transforms, with and without strata, clusters and ties, to \eqn{10^{-8}}. In a Monte Carlo study
#' (1000 replications, \eqn{n} = 200 to 1000) the size at the 5\% level was 4.5 to 5.9\% and identical to
#' that of \code{cox.zph} on the same data; with crossing hazards the power was 100\%.
#' @return
#' An S7 object of class \code{PHTest} with \code{print}, \code{plot} and \code{as.data.frame}
#' methods. \code{x@table} is a data frame of all tests (cause, term, chi-square, df, p; term
#' \code{GLOBAL} for the global test); \code{x@tests} has one element per cause (one for a Cox model)
#' with the table and the scaled Schoenfeld residuals used by the plot.
#' @references
#' \emph{Test and residual plot:}
#' Grambsch PM, Therneau TM (1994). Proportional hazards tests and diagnostics based on weighted
#'   residuals. Biometrika 81:515-526.
#' Therneau TM, Grambsch PM (2000). Modeling Survival Data: Extending the Cox Model. Springer,
#'   Chapter 6.
#'
#' \emph{Schoenfeld residuals:}
#' Schoenfeld D (1982). Partial residuals for the proportional hazards regression model.
#'   Biometrika 69:239-241.
#'
#' All references: \code{sx_references("ph_test")}.
#' @seealso
#' \code{\link{cox}} (\code{strata()}, \code{tt()}), \code{\link{time_effects}}, \code{\link{residuals.SurvX}},
#' \code{\link{rmst}}.
#' @examples
#' set.seed(4)
#' n <- 300
#' x <- rbinom(n, 1, 0.5)
#' t <- rweibull(n, ifelse(x == 1, 0.7, 1.5), 3); cens <- runif(n, 0, 6)   # non-proportional hazards
#' d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = rnorm(n))
#'
#' pt <- ph_test(cox(Event(time, status) ~ x + z, d))
#' pt                       # x violates proportional hazards, z does not
#' plot(pt, "x")            # the effect of x changes over time
#' as.data.frame(pt)
#' @export
ph_test <- function(fit, transform = c("km", "rank", "identity", "log")) {
  transform <- match.arg(transform)
  if (S7::S7_inherits(fit, CoxFit)) {
    .no_tt(fit, "ph_test()")
    r <- .zph_core(fit@time, fit@status, fit@X, fit@coefficients, fit@vcov_model %||% fit@vcov, transform, fit@weights,
                   fit@entry, fit@strata)
    r$label <- "the event"
    return(.ph_object(list(r), transform))
  }
  if (S7::S7_inherits(fit, CSCFit)) {
    out <- lapply(seq_len(fit@K), function(k) {
      m <- fit@models[[k]]
      r <- .zph_core(fit@time, as.integer(fit@status == k), fit@X, m$coefficients, m$vcov_model %||% m$vcov, transform,
                     NULL, fit@entry, fit@strata)
      r$label <- paste("cause", k); r
    })
    return(.ph_object(out, transform))
  }
  stop("ph_test() takes a model fitted by cox() or csc() with a Cox baseline", call. = FALSE)
}

.ph_object <- function(tests, transform) {
  tab <- do.call(rbind, lapply(seq_along(tests), function(k) {
    t <- tests[[k]]$table
    data.frame(cause = if (length(tests) > 1L) k else NA_integer_, term = rownames(t), chisq = t$chisq, df = t$df, p = t$p,
               stringsAsFactors = FALSE)
  }))
  rownames(tab) <- NULL
  PHTest(tests = tests, transform = transform, table = tab)
}

S7::method(print, PHTest) <- function(x, ...) {
  cat("Test of proportional hazards (score test of Grambsch and Therneau; time scale: ", x@transform, ")\n", sep = "")
  for (r in x@tests) {
    tab <- r$table
    df <- data.frame(term = rownames(tab), chisq = .num(tab$chisq, 2L), df = tab$df, p = .pval(tab$p), stringsAsFactors = FALSE)
    names(df) <- c("", "chi-square", "df", "p")
    cat("\n", if (length(x@tests) > 1L) paste0("Cause ", sub("cause ", "", r$label), "\n") else "", sep = "")
    cat(.table_lines(df), sep = "\n")
    bad <- rownames(tab)[-nrow(tab)][tab$p[-nrow(tab)] < 0.05]
    cat(if (length(bad)) paste0("Reading: the effect of ", paste(bad, collapse = ", "),
                                " changes over time (p < 0.05); proportional hazards are doubtful for ",
                                if (length(bad) > 1L) "these covariates." else "this covariate.")
        else "Reading: no evidence against proportional hazards (all p >= 0.05).", "\n")
  }
  cat("\nplot(x, var) shows the scaled Schoenfeld residuals: a flat trend supports proportional hazards.\n")
  if (any(x@table$p[x@table$term != "GLOBAL"] < 0.05, na.rm = TRUE))
    cat("time_effects(fit) compares forms of the change over time and refits the model with them.\n")
  cat("\n", .method_lines(.method_keys(x)), sep = "\n")
  invisible(x)
}

S7::method(plot, PHTest) <- function(x, var = NULL, cause = 1L, df = 4L, level = 0.95, col = "#0072B2", xlab = NULL, ylab = NULL,
                         main = NULL, cex = 0.85, ...) {
  r <- x@tests[[min(cause, length(x@tests))]]
  vars <- if (is.null(var)) colnames(r$y) else var
  op <- graphics::par(mfrow = c(1L, length(vars)), cex = cex, mar = c(4.5, 4.5, 3, 1)); on.exit(graphics::par(op))
  for (v in vars) {
    y <- r$y[, v]; xx <- r$x
    ok <- is.finite(xx) & is.finite(y)
    graphics::plot(xx[ok], y[ok], pch = 16, cex = 0.5, col = grDevices::adjustcolor("grey30", 0.5),
                   xlab = if (is.null(xlab)) paste0("Time (", r$transform, " scale)") else xlab,
                   ylab = if (is.null(ylab)) paste0("Beta(t) for ", v) else ylab,
                   main = if (is.null(main)) (if (length(x@tests) > 1L) r$label else "") else main, ...)
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
#' @description
#' With a formula: the classic linearising plots of the non-parametric (Nelson-Aalen) cumulative
#' hazard, by group; a straight line in a panel indicates the distribution of that panel. With a
#' fitted model: the Cox-Snell residual plot.
#' @param x a formula \code{Event(time, status) ~ 1} or \code{~ group}, or a model fitted by
#'   \code{\link{parsurv}}, \code{\link{cox}} or \code{\link{csc}} (Cox or parametric baseline).
#' @param data a data frame (with a formula).
#' @param cause the cause for competing risks data (cause-specific hazard).
#' @param col colours of the groups.
#' @param cex character expansion.
#' @param ... passed to \code{plot.default}.
#' @details
#' Panels for a formula, with \eqn{\hat H} the Nelson-Aalen estimate and \eqn{\hat S = e^{-\hat H}}:
#' \describe{
#'   \item{exponential}{\eqn{\hat H(t)} against \eqn{t}: a straight line through the origin.}
#'   \item{Weibull}{\eqn{\log \hat H(t)} against \eqn{\log t}: a straight line (slope = shape;
#'     parallel lines across groups indicate proportional hazards).}
#'   \item{log-normal}{\eqn{\Phi^{-1}(1 - \hat S(t))} against \eqn{\log t}: a straight line.}
#'   \item{log-logistic}{\eqn{\log\{(1 - \hat S)/\hat S\}} against \eqn{\log t}: a straight line
#'     (parallel lines indicate proportional odds).}
#' }
#' For a fitted model, the Nelson-Aalen cumulative hazard of the Cox-Snell residuals (with the event
#' indicators) is plotted against the residuals. If the model fits, the residuals behave like a
#' censored sample from the unit exponential distribution and the curve follows the line of identity.
#' @section Interpreting the plot:
#' With a formula: look for the panel in which the points lie closest to a straight line; that
#' distribution is a good candidate for \code{\link{parsurv}}. Parallel lines across groups in the
#' Weibull panel indicate proportional hazards; in the log-logistic panel, proportional odds. With a
#' fitted model: Cox-Snell residuals that follow the identity line support the model; a curve that
#' bends away shows lack of fit, usually in the tail first.
#' @return
#' invisibly, the plotted data.
#' @references
#' Klein JP, Moeschberger ML (2003) Survival Analysis: Techniques for Censored and Truncated Data,
#'   2nd ed. Springer, Sections 12.2-12.5.
#'
#' Collett D (2015) Modelling Survival Data in Medical Research, 3rd ed. CRC Press, Chapter 5.
#'
#' Cox DR, Snell EJ (1968) A general definition of residuals. Journal of the Royal Statistical Society,
#'   Series B 30:248-275.
#' @seealso
#' \code{\link{hazard_plot}}, \code{\link{hazard_np}}, \code{\link{residuals.SurvX}},
#' \code{\link{parsurv}}.
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
