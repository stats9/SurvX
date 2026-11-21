#' Adjusted absolute risks and causal contrasts between treatments
#'
#' @description
#' Estimates, for each value of a treatment (or exposure), the absolute risk of the event, or the
#' cumulative incidence of a cause, that would have been seen had everyone in the data received that
#' value, and the differences and ratios between treatments, at chosen times and as curves over
#' follow-up. Two methods are available and can be compared: standardisation of a fitted outcome
#' model (the g-formula) and inverse probability of treatment weighting. Intervals come from a
#' bootstrap in which the models are refitted.
#'
#' \strong{Use it when} the question is "how much would the risk change if everyone were treated
#' rather than not", in an observational study or to turn the hazard ratios of a trial into absolute
#' risks adjusted for covariates. Hazard ratios do not answer that question; absolute risks do.
#'
#' \strong{Use something else when} you want the risk of a particular patient (\code{predict()} of
#' the model) or a summary that needs no model (\code{\link{km}}, \code{\link{aj}} by group).
#' @param fit the outcome model: a fit of \code{\link{cox}}, \code{\link{fg}}, \code{\link{csc}},
#'   \code{\link{directcif}}, \code{\link{parsurv}}, \code{\link{cure}} or the result of
#'   \code{\link{time_effects}}. It must contain the treatment as a covariate (or, for \code{cox} and
#'   \code{csc}, as a stratum) and should contain the confounders. It must have been fitted with
#'   \code{data = }, so that it can be refitted in bootstrap samples.
#' @param data the data frame the model was fitted to.
#' @param treatment the name of the treatment variable (a column of \code{data}): a factor, character,
#'   logical or 0/1 variable with at most 10 values.
#' @param times times at which risks and contrasts are reported. Default: the quartiles of the event
#'   times of the cause (two significant digits). Curves are computed on a grid over follow-up in
#'   addition.
#' @param cause competing risks: the cause whose cumulative incidence is compared (for \code{fg} the
#'   cause of the model).
#' @param method \code{"gformula"} (default): standardisation of the outcome model.
#'   \code{"ipw"}: inverse probability of treatment weighting with a model for the treatment.
#'   \code{"both"}: both, side by side.
#' @param ps_formula the model of the treatment for \code{"ipw"}: \code{treatment ~ covariates}.
#'   Default: the treatment on all other terms of the outcome model (main effects).
#' @param trim truncation of the weights at the \code{trim} and \code{1 - trim} quantiles (for example
#'   0.01) when some weights are extreme; 0 (default) for none.
#' @param B number of bootstrap samples; 0 for point estimates only.
#' @param seed seed of the bootstrap; the state of the random number generator of the session is
#'   restored afterwards.
#' @param level confidence level.
#' @param reference the treatment value against which the others are compared; default the first
#'   level.
#' @details
#' \subsection{Standardisation (g-formula)}{
#' For each treatment value \eqn{a}, the outcome model predicts the risk \eqn{F(t \mid a, z_i)} of
#' every subject with the treatment set to \eqn{a} and the other covariates \eqn{z_i} as observed, and
#' the predictions are averaged: \eqn{\hat F_a(t) = n^{-1} \sum_i \hat F(t \mid a, z_i)}
#' (Robins 1986; Hernan and Robins 2020, Chapter 13). For one event type the risk is one minus the
#' predicted survival; for competing risks the predicted cumulative incidence of the cause. The model
#' can contain interactions of the treatment with covariates and, through \code{\link{time_effects}},
#' time-varying effects.
#' }
#' \subsection{Inverse probability of treatment weighting}{
#' A logistic (two values) or multinomial logistic (more values) model gives the probability
#' \eqn{e_a(z)} of each treatment; subjects get the stabilised weight \eqn{P(A = a)/e_A(z)}, and the
#' weighted Kaplan-Meier (one event type) or Aalen-Johansen estimate (competing risks) is computed in
#' each treatment group (Cole and Hernan 2004; Neumann and Billionnet 2016). Delayed entry of the
#' outcome model is used. Balance is shown by standardised mean differences of the covariates before
#' and after weighting, and the weights by the effective sample size and the largest weight (Austin
#' and Stuart 2015).
#' }
#' \subsection{Contrasts and intervals}{
#' Risk differences \eqn{\hat F_a(t) - \hat F_{ref}(t)} and ratios \eqn{\hat F_a(t)/\hat F_{ref}(t)}.
#' With \code{B > 0} the subjects (clusters, if the model has them) are resampled, the outcome model is
#' refitted with its own call (and the treatment model refitted), and percentile intervals are
#' reported; the p-value of a difference uses the bootstrap standard error. Samples in which a
#' treatment value is missing or a model fails are skipped and counted.
#' }
#' @section Interpreting the result:
#' \describe{
#'   \item{Risk by treatment}{the risk that would have been observed by time \eqn{t} had everyone
#'     received that treatment, in a population with the covariate distribution of the data.}
#'   \item{Risk difference}{in percentage points: the change in absolute risk from treating everyone
#'     with one value rather than the reference. Its inverse is a number needed to treat.}
#'   \item{Competing risks}{the contrast is a total effect on the cumulative incidence, which includes
#'     the effect of treatment on the competing causes: a treatment that prevents the competing event
#'     can increase the incidence of the cause of interest.}
#'   \item{\code{method = "both"}}{the two methods rest on different models, the outcome model and the
#'     treatment model. Close agreement is reassuring; a large difference means that at least one model
#'     is wrong.}
#' }
#' A sentence for a report: \emph{"Standardised to the covariate distribution of the cohort, the
#' 5-year risk of death was 18.2\% had all patients been treated and 24.9\% had none been treated, a
#' difference of -6.7 percentage points (95\% bootstrap CI -10.9 to -2.5)."}
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{No unmeasured confounding}{every common cause of treatment and outcome is among the
#'     covariates. This cannot be checked from the data.}
#'   \item{Positivity}{every subject could have received every treatment: look at the range of the
#'     propensity scores and the largest weights (\code{method = "ipw"} or \code{"both"}).}
#'   \item{Correct models}{the outcome model for the g-formula (check it with \code{\link{ph_test}},
#'     \code{\link{time_effects}} and \code{\link{model_scan}}), the treatment model for IPW
#'     (balance after weighting).}
#'   \item{Baseline covariates}{one row per subject; the g-formula with covariates that change over
#'     follow-up needs methods not provided here.}
#' }
#' @section Messages and common problems:
#' \describe{
#'   \item{\code{data does not match the model}}{give the same data frame that the model was fitted
#'     to; rows with missing values in the variables of the model are dropped as in the fit.}
#'   \item{\code{the model call has no 'data' argument}}{fit the model with \code{data = } so that it
#'     can be refitted in bootstrap samples, or use \code{B = 0}.}
#'   \item{Long computing times}{the model is refitted \code{B} times; check with \code{B = 0} first.}
#' }
#' @section Validation:
#' With a binary treatment and no other covariates both methods reduce to the Kaplan-Meier or
#' Aalen-Johansen estimate in each group, which is checked in the test suite, as are the weighted
#' estimator against the weighted estimates of \pkg{survival} and the multinomial model against
#' \code{stats::glm} for two values.
#' @return
#' An S7 object of class \code{AdjustedRisk} with \code{x@risk} (method, level, time, risk, se, lower,
#' upper), \code{x@contrast} (method, level, reference, time, difference, se, lower, upper, ratio,
#' ratio_lower, ratio_upper, p), \code{x@curve} (the risks on a grid over follow-up) and \code{x@ps}
#' (propensity score ranges, weights and balance for IPW).
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(x)} \tab risks, contrasts and a reading \cr
#'   \code{summary(x)} \tab also the balance of every covariate \cr
#'   \code{plot(x, method = NULL)} \tab adjusted risk curves with their bands \cr
#'   \code{as.data.frame(x, what)} \tab \code{"risk"}, \code{"contrast"}, \code{"curve"} or \code{"balance"} \cr
#'   \code{report(x, file)} \tab a self-contained HTML report \cr
#'   \code{\link{sx_references}(x)} \tab the references of the methods \cr
#' }
#' @references
#' \emph{Standardisation and causal inference}
#'
#' Robins J (1986) A new approach to causal inference in mortality studies with a sustained exposure
#' period: application to control of the healthy worker survivor effect. Math Model 7:1393-1512.
#'
#' Hernan MA, Robins JM (2020) Causal Inference: What If. Chapman and Hall/CRC.
#'
#' Ozenne BMH, Scheike TH, Staerk L, Gerds TA (2020) On the estimation of average treatment effects
#' with right-censored time to event outcome and competing risks. Biom J 62:751-763.
#'
#' \emph{Inverse probability of treatment weighting}
#'
#' Cole SR, Hernan MA (2004) Adjusted survival curves with inverse probability weights. Comput Methods
#' Programs Biomed 75:45-49.
#'
#' Neumann A, Billionnet C (2016) Covariate adjustment of cumulative incidence functions for
#' competing risks data using inverse probability of treatment weighting. Comput Methods Programs
#' Biomed 129:63-70.
#'
#' Austin PC, Stuart EA (2015) Moving towards best practice when using inverse probability of
#' treatment weighting (IPTW) using the propensity score to estimate causal treatment effects in
#' observational studies. Stat Med 34:3661-3679.
#'
#' \emph{Bootstrap}
#'
#' Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap. Chapman and Hall.
#'
#' All references: \code{sx_references("adjusted_risk")}.
#' @seealso
#' \code{\link{predict.SurvX}} (risks of given subjects), \code{\link{time_effects}},
#' \code{\link{km}}, \code{\link{aj}}, \code{\link{rmst}}.
#' @examples
#' set.seed(2)
#' n <- 400
#' age <- rnorm(n)
#' trt <- rbinom(n, 1, plogis(0.8 * age))                 # older patients treated more often
#' t <- rexp(n, 0.2 * exp(-0.5 * trt + 0.7 * age)); cens <- runif(n, 0, 8)
#' d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens),
#'                 trt = factor(trt, 0:1, c("control", "treated")), age = age)
#'
#' fit <- cox(Event(time, status) ~ trt + age, data = d)
#' ar <- adjusted_risk(fit, d, "trt", times = c(2, 5), method = "both", B = 0)
#' ar
#' plot(ar)
#'
#' \donttest{
#' ## with bootstrap intervals
#' ar_b <- adjusted_risk(fit, d, "trt", times = c(2, 5), B = 200)
#' ar_b
#' report(ar_b, file = tempfile(fileext = ".html"))
#'
#' ## competing risks: cumulative incidence of live birth by education
#' data(fertility)
#' adjusted_risk(fg(Event(time, status) ~ education, data = fertility, cause = 1),
#'               fertility, "education", times = c(1, 2, 5), B = 50)
#' }
#' @export
adjusted_risk <- function(fit, data, treatment, times = NULL, cause = 1L, method = c("gformula", "ipw", "both"),
                          ps_formula = NULL, trim = 0, B = 200L, seed = 1L, level = 0.95, reference = NULL) {
  cl <- match.call()
  method <- match.arg(method)
  if (missing(data) || !is.data.frame(data)) stop("data: give the data frame the model was fitted to", call. = FALSE)
  if (missing(treatment) || !is.character(treatment) || length(treatment) != 1L || !treatment %in% names(data))
    stop("treatment: the name of a column of data, for example treatment = \"arm\"", call. = FALSE)
  if (!is.numeric(trim) || length(trim) != 1L || trim < 0 || trim >= 0.5) stop("trim must be a number in [0, 0.5)", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || !(level > 0 && level < 1)) stop("level must be between 0 and 1", call. = FALSE)
  B <- as.integer(B)
  if (length(B) != 1L || is.na(B) || B < 0L) stop("B must be 0 or a number of bootstrap samples", call. = FALSE)
  if (B > 0L && B < 20L) stop("B: use at least 20 bootstrap samples (200 or more for intervals)", call. = FALSE)
  info <- .ar_info(fit, as.integer(cause))
  vars <- unique(c(info$vars, treatment))
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("data does not have the variables of the model: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!treatment %in% info$vars) stop("the treatment '", treatment, "' is not a covariate of the model", call. = FALSE)
  keep <- stats::complete.cases(data[, vars, drop = FALSE])
  d <- data[keep, , drop = FALSE]
  if (nrow(d) != length(info$time)) stop("data does not match the model: ", nrow(d), " complete rows, ", length(info$time),
                                         " in the model. Give the data frame the model was fitted to.", call. = FALSE)
  a <- d[[treatment]]
  lev <- if (is.factor(a)) levels(droplevels(a)) else sort(unique(as.character(a)))
  if (length(lev) < 2L) stop("the treatment has fewer than two values", call. = FALSE)
  if (length(lev) > 10L) stop("the treatment has more than 10 values; adjusted_risk() is for a categorical treatment", call. = FALSE)
  if (!is.null(reference)) {
    if (!as.character(reference) %in% lev) stop("reference must be one of: ", paste(lev, collapse = ", "), call. = FALSE)
    lev <- c(as.character(reference), setdiff(lev, as.character(reference)))
  }
  ev <- info$time[info$status == info$cause]
  if (length(ev) < 5L) stop("too few events of the cause of interest", call. = FALSE)
  times <- sort(unique(times %||% signif(stats::quantile(ev, c(0.25, 0.5, 0.75), type = 1, names = FALSE), 2L)))
  if (any(!is.finite(times)) || any(times <= 0)) stop("times must be positive numbers", call. = FALSE)
  grid <- seq(0, max(ev), length.out = 61L)[-1L]
  tt <- sort(unique(c(times, grid)))
  psf <- if (method != "gformula") .ar_ps_formula(ps_formula, fit, info, treatment) else NULL
  env <- parent.frame()
  est <- .ar_estimate(fit, d, info, treatment, lev, tt, method, psf, trim)
  boot <- NULL
  if (B > 0L) boot <- .ar_bootstrap(fit, d, treatment, lev, tt, method, psf, trim, B, seed, env, info)
  out <- .ar_tables(est, boot, lev, times, tt, level)
  AdjustedRisk(model_call = cl, fit = fit, treatment = treatment, levels = lev, times = times, cause = info$cause,
               kind = info$kind, method = method, risk = out$risk, contrast = out$contrast, curve = out$curve,
               ps = est$ps, B = B, ok = if (is.null(boot)) 0L else boot$ok, level = level, n = nrow(d),
               ps_formula = psf)
}

# ----------------------------------------------------------------------------------------------
# what the model is and how it predicts the risk of the cause of interest
# ----------------------------------------------------------------------------------------------
.ar_info <- function(fit, cause) {
  base <- if (S7::S7_inherits(fit, TimeEffects)) fit@fit else fit
  ok <- S7::S7_inherits(base, CoxFit) || S7::S7_inherits(base, FineGrayFit) || S7::S7_inherits(base, CSCFit) ||
        S7::S7_inherits(base, DirectCIFFit) || S7::S7_inherits(base, ParSurvFit) || S7::S7_inherits(base, ParCSCFit) ||
        S7::S7_inherits(base, CureFit)
  if (!ok) stop("adjusted_risk() takes a model fitted by cox(), fg(), csc(), directcif(), parsurv(), cure() or time_effects()", call. = FALSE)
  single <- S7::S7_inherits(base, CoxFit) || S7::S7_inherits(base, ParSurvFit) || S7::S7_inherits(base, CureFit)
  if (S7::S7_inherits(base, CoxFit) && !is.null(base@tt))
    stop("models with tt() terms are not supported; use time_effects() and give its result instead", call. = FALSE)
  if (S7::S7_inherits(base, FineGrayFit)) cause <- base@cause
  if (single) cause <- 1L
  vars <- all.vars(base@terms)
  sv <- tryCatch(base@strata_vars, error = function(e) NULL)
  if (!is.null(sv)) vars <- c(vars, unlist(lapply(sv, function(s) all.vars(str2lang(s)))))
  ent <- tryCatch(base@entry, error = function(e) NULL)
  risk <- if (single) {
    if (S7::S7_inherits(fit, TimeEffects)) function(f, nd, t) 1 - as.matrix(predict(f, nd, times = t))
    else function(f, nd, t) 1 - as.matrix(predict(f, nd, times = t, type = "survival"))
  } else if (S7::S7_inherits(base, FineGrayFit) && !S7::S7_inherits(fit, TimeEffects)) {
    function(f, nd, t) as.matrix(predict(f, nd, times = t))
  } else function(f, nd, t) as.matrix(predict(f, nd, times = t, cause = cause))
  list(base = base, kind = if (single) "single" else "competing", cause = as.integer(cause), vars = unique(vars),
       time = base@time, status = base@status, entry = ent, risk = risk)
}

.ar_ps_formula <- function(ps_formula, fit, info, treatment) {
  if (!is.null(ps_formula)) {
    if (!inherits(ps_formula, "formula") || length(ps_formula) != 3L) stop("ps_formula must be a two-sided formula: treatment ~ covariates", call. = FALSE)
    if (!identical(all.vars(ps_formula[[2L]]), treatment)) stop("the left-hand side of ps_formula must be the treatment", call. = FALSE)
    return(ps_formula)
  }
  lab <- attr(stats::delete.response(info$base@terms), "term.labels")
  lab <- lab[!vapply(lab, function(l) treatment %in% all.vars(str2lang(l)), logical(1))]
  sv <- tryCatch(info$base@strata_vars, error = function(e) NULL)
  if (!is.null(sv)) lab <- c(lab, sv[!vapply(sv, function(s) treatment %in% all.vars(str2lang(s)), logical(1))])
  stats::reformulate(if (length(lab)) lab else "1", response = treatment)
}

# ----------------------------------------------------------------------------------------------
# estimation on one data set: standardisation (g-formula) and inverse probability weighting
# ----------------------------------------------------------------------------------------------
.ar_estimate <- function(fit, d, info, treatment, lev, tt, method, psf, trim) {
  out <- list()
  if (method %in% c("gformula", "both")) {
    out$gformula <- t(vapply(lev, function(l) {
      nd <- d
      nd[[treatment]] <- if (is.factor(d[[treatment]])) factor(l, levels = levels(d[[treatment]]))
                         else if (is.logical(d[[treatment]])) as.logical(l)
                         else if (is.numeric(d[[treatment]])) as.numeric(l) else l
      r <- info$risk(fit, nd, tt)
      colMeans(r, na.rm = TRUE)
    }, numeric(length(tt))))
  }
  if (method %in% c("ipw", "both")) {
    ps <- .ar_ps(d, treatment, lev, psf, trim)
    a <- as.character(d[[treatment]])
    out$ipw <- t(vapply(lev, function(l) {
      i <- a == l
      .ar_wcif(info$time[i], info$status[i], if (is.null(info$entry)) NULL else info$entry[i], ps$w[i], info$cause, tt)
    }, numeric(length(tt))))
    out$ps <- ps$diag
  }
  out
}

# multinomial (binary: logistic) model of the treatment and stabilised weights
.ar_ps <- function(d, treatment, lev, psf, trim) {
  mf <- stats::model.frame(psf, d, na.action = stats::na.pass)
  X <- stats::model.matrix(psf, mf)
  a <- factor(as.character(d[[treatment]]), levels = lev)
  P <- .ar_mlogit(X, a)
  pa <- P[cbind(seq_len(nrow(P)), as.integer(a))]
  marg <- as.numeric(table(a)) / length(a)
  w <- marg[as.integer(a)] / pa
  if (trim > 0) {
    q <- stats::quantile(w, c(trim, 1 - trim), names = FALSE)
    w <- pmin(pmax(w, q[1L]), q[2L])
  }
  Xc <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  smd <- if (ncol(Xc)) {
    ref <- a == lev[1L]
    do.call(rbind, lapply(lev[-1L], function(l) {
      o <- a == l
      s <- sqrt((apply(Xc[ref, , drop = FALSE], 2L, stats::var) + apply(Xc[o, , drop = FALSE], 2L, stats::var)) / 2)
      s[!is.finite(s) | s == 0] <- NA_real_
      wm <- function(i) colSums(Xc[i, , drop = FALSE] * w[i]) / sum(w[i])
      data.frame(covariate = colnames(Xc), level = l,
                 before = (colMeans(Xc[o, , drop = FALSE]) - colMeans(Xc[ref, , drop = FALSE])) / s,
                 after = (wm(o) - wm(ref)) / s, stringsAsFactors = FALSE, row.names = NULL)
    }))
  } else NULL
  ess <- vapply(lev, function(l) { i <- a == l; sum(w[i])^2 / sum(w[i]^2) }, numeric(1))
  diag <- list(ps_range = t(vapply(lev, function(l) range(P[, l]), numeric(2))),
               weights = data.frame(level = lev, n = as.integer(table(a)), ess = unname(ess),
                                    max_weight = vapply(lev, function(l) max(w[a == l]), numeric(1)), row.names = NULL),
               smd = smd, trim = trim)
  list(w = w, diag = diag)
}

# multinomial logistic regression by Newton-Raphson; returns the n x K matrix of probabilities
.ar_mlogit <- function(X, a) {
  K <- nlevels(a); n <- nrow(X); p <- ncol(X)
  Y <- matrix(0, n, K); Y[cbind(seq_len(n), as.integer(a))] <- 1
  b <- matrix(0, p, K - 1L)
  prob <- function(b) { e <- cbind(0, X %*% b); e <- exp(e - apply(e, 1L, max)); e / rowSums(e) }
  ll <- function(P) sum(Y * log(pmax(P, 1e-300)))
  P <- prob(b); l0 <- ll(P)
  for (it in seq_len(100L)) {
    g <- as.vector(crossprod(X, Y[, -1L, drop = FALSE] - P[, -1L, drop = FALSE]))
    H <- matrix(0, p * (K - 1L), p * (K - 1L))
    for (j in seq_len(K - 1L)) for (k in seq_len(K - 1L)) {
      wjk <- P[, j + 1L] * ((j == k) - P[, k + 1L])
      H[(j - 1L) * p + seq_len(p), (k - 1L) * p + seq_len(p)] <- crossprod(X * wjk, X)
    }
    step <- tryCatch(solve(H + diag(1e-8, nrow(H)), g), error = function(e) NULL)
    if (is.null(step)) break
    s <- 1
    repeat {
      bn <- b + matrix(step * s, p, K - 1L)
      Pn <- prob(bn); ln <- ll(Pn)
      if (is.finite(ln) && ln >= l0 - 1e-10) break
      s <- s / 2
      if (s < 1e-8) break
    }
    b <- bn; P <- Pn
    if (abs(ln - l0) < 1e-10 * (abs(l0) + 1)) break
    l0 <- ln
  }
  colnames(P) <- levels(a)
  P
}

# weighted Aalen-Johansen estimate of the cumulative incidence of cause k (one event type: 1 - KM)
.ar_wcif <- function(time, status, entry, w, k, tt) {
  et <- sort(unique(time[status > 0L]))
  if (!length(et)) return(rep(0, length(tt)))
  Y <- vapply(et, function(s) sum(w[time >= s & (if (is.null(entry)) TRUE else entry < s)]), numeric(1))
  dall <- vapply(et, function(s) sum(w[time == s & status > 0L]), numeric(1))
  dk <- vapply(et, function(s) sum(w[time == s & status == k]), numeric(1))
  ok <- Y > 0
  haz_all <- ifelse(ok, dall / Y, 0); haz_k <- ifelse(ok, dk / Y, 0)
  S <- cumprod(1 - haz_all)
  Sm <- c(1, S[-length(S)])
  F <- cumsum(Sm * haz_k)
  .step_at(et, F, tt, 0)
}

# ----------------------------------------------------------------------------------------------
# bootstrap: the outcome model (and the treatment model) refitted in every sample
# ----------------------------------------------------------------------------------------------
.ar_refit <- function(fit, d, env) {
  if (S7::S7_inherits(fit, TimeEffects)) {
    inner <- .ar_refit(fit@fit, d, env)
    cl <- fit@model_call
    cl$fit <- quote(.ar_inner); cl$bootstrap <- 0L
    e <- new.env(parent = env); e$.ar_inner <- inner
    return(eval(cl, e))
  }
  cl <- fit@model_call
  if (is.null(cl$data)) stop("the model call has no 'data' argument; refit the model with data = ", call. = FALSE)
  cl$data <- quote(.ar_data)
  e <- new.env(parent = env); e$.ar_data <- d
  eval(cl, e)
}

.ar_bootstrap <- function(fit, d, treatment, lev, tt, method, psf, trim, B, seed, env, info) {
  cluster <- tryCatch(info$base@cluster, error = function(e) NULL)
  cl <- cluster %||% seq_len(nrow(d))
  ids <- unique(cl)
  rows_of <- split(seq_along(cl), match(cl, ids))
  res <- .with_seed(seed, lapply(seq_len(B), function(b) {
    rows <- unlist(rows_of[sample.int(length(ids), length(ids), replace = TRUE)], use.names = FALSE)
    db <- d[rows, , drop = FALSE]
    if (any(table(factor(as.character(db[[treatment]]), levels = lev)) == 0L)) return(NULL)
    tryCatch(suppressWarnings({
      fb <- .ar_refit(fit, db, env)
      ib <- .ar_info(fb, info$cause)
      e <- .ar_estimate(fb, db, ib, treatment, lev, tt, method, psf, trim)
      e[c("gformula", "ipw")]
    }), error = function(e) NULL)
  }))
  res <- Filter(Negate(is.null), res)
  list(res = res, ok = length(res), B = B)
}

# ----------------------------------------------------------------------------------------------
# tables of risks, contrasts and curves
# ----------------------------------------------------------------------------------------------
.ar_tables <- function(est, boot, lev, times, tt, level) {
  z <- .zq(level); a <- (1 - level) / 2
  ms <- intersect(c("gformula", "ipw"), names(est))
  have_b <- !is.null(boot) && boot$ok >= max(10L, boot$B / 2)
  pick <- function(m) if (have_b) lapply(boot$res, `[[`, m) else NULL
  risk <- list(); con <- list(); curve <- list()
  for (m in ms) {
    R <- est[[m]]                                   # levels x tt
    bs <- pick(m)
    arr <- if (have_b) simplify2array(bs) else NULL # levels x tt x B
    for (i in seq_along(lev)) {
      sd <- if (have_b) apply(arr[i, , , drop = FALSE], 2L, stats::sd) else rep(NA_real_, length(tt))
      lo <- if (have_b) apply(arr[i, , , drop = FALSE], 2L, stats::quantile, probs = a, names = FALSE) else rep(NA_real_, length(tt))
      hi <- if (have_b) apply(arr[i, , , drop = FALSE], 2L, stats::quantile, probs = 1 - a, names = FALSE) else rep(NA_real_, length(tt))
      df <- data.frame(method = m, level = lev[i], time = tt, risk = R[i, ], se = sd, lower = lo, upper = hi, stringsAsFactors = FALSE)
      curve[[length(curve) + 1L]] <- df
      risk[[length(risk) + 1L]] <- df[df$time %in% times, , drop = FALSE]
    }
    for (i in seq_along(lev)[-1L]) {
      j <- tt %in% times
      dd <- R[i, j] - R[1L, j]; rr <- R[i, j] / R[1L, j]
      if (have_b) {
        D <- matrix(arr[i, j, ] - arr[1L, j, ], nrow = sum(j))
        Q <- matrix(log(arr[i, j, ] / arr[1L, j, ]), nrow = sum(j))
        sdd <- apply(D, 1L, stats::sd); sdl <- apply(Q, 1L, function(v) stats::sd(v[is.finite(v)]))
        dlo <- apply(D, 1L, stats::quantile, probs = a, names = FALSE); dhi <- apply(D, 1L, stats::quantile, probs = 1 - a, names = FALSE)
        rlo <- exp(apply(Q, 1L, function(v) stats::quantile(v[is.finite(v)], a, names = FALSE)))
        rhi <- exp(apply(Q, 1L, function(v) stats::quantile(v[is.finite(v)], 1 - a, names = FALSE)))
        p <- 2 * stats::pnorm(-abs(dd / sdd))
      } else { sdd <- dlo <- dhi <- rlo <- rhi <- p <- rep(NA_real_, sum(j)) }
      con[[length(con) + 1L]] <- data.frame(method = m, level = lev[i], reference = lev[1L], time = tt[j], difference = dd, se = sdd,
                                            lower = dlo, upper = dhi, ratio = rr, ratio_lower = rlo, ratio_upper = rhi, p = p,
                                            stringsAsFactors = FALSE)
    }
  }
  list(risk = do.call(rbind, risk), contrast = do.call(rbind, con), curve = do.call(rbind, curve))
}

# ----------------------------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------------------------
.ar_mname <- function(m) c(gformula = "standardisation (g-formula)", ipw = "inverse probability of treatment weighting")[[m]]

.ar_lines <- function(x, detail = FALSE) {
  what <- if (x@kind == "single") "risk of the event (1 - survival)" else paste0("cumulative incidence of cause ", x@cause)
  lv <- .pct(x@level)
  out <- c(paste0("Adjusted absolute risks: ", what),
           paste0("  Treatment: ", x@treatment, " (", paste(x@levels, collapse = ", "), "; reference ", x@levels[1L], ")   subjects: ", x@n),
           paste0("  Method: ", paste(vapply(intersect(c("gformula", "ipw"), unique(x@risk$method)), .ar_mname, character(1)), collapse = " and ")),
           if (x@B > 0L) paste0("  Intervals: percentile bootstrap, ", x@ok, " of ", x@B, " samples (models refitted in every sample)")
           else "  Intervals: none (B = 0); use B = 200 or more", "")
  for (m in unique(x@risk$method)) {
    r <- x@risk[x@risk$method == m, , drop = FALSE]
    tab <- data.frame(time = format(r$time), level = r$level, risk = .num(r$risk, 3L),
                      ci = if (x@B > 0L) .ci_text(r$lower, r$upper, 3L) else "-", check.names = FALSE)
    names(tab)[4L] <- paste0(lv, " CI")
    out <- c(out, paste0("Risk by treatment, ", .ar_mname(m)), .table_lines(tab), "")
    cc <- x@contrast[x@contrast$method == m, , drop = FALSE]
    ct <- data.frame(time = format(cc$time), contrast = paste0(cc$level, " vs ", cc$reference),
                     difference = .num(cc$difference, 3L), ci = if (x@B > 0L) .ci_text(cc$lower, cc$upper, 3L) else "-",
                     ratio = .num(cc$ratio, 2L), rci = if (x@B > 0L) .ci_text(cc$ratio_lower, cc$ratio_upper) else "-",
                     p = .pval(cc$p), check.names = FALSE)
    names(ct)[c(4L, 6L)] <- c(paste0(lv, " CI"), paste0(lv, " CI "))
    out <- c(out, "Risk difference and risk ratio", .table_lines(ct), "")
  }
  # reading at the last time
  tl <- max(x@times)
  m1 <- unique(x@risk$method)[1L]
  r <- x@risk[x@risk$method == m1 & x@risk$time == tl, , drop = FALSE]
  cc <- x@contrast[x@contrast$method == m1 & x@contrast$time == tl, , drop = FALSE]
  out <- c(out, "Reading:",
           paste0("- At t = ", format(tl), ", had everyone received ", r$level, ", the ", if (x@kind == "single") "risk" else "cumulative incidence",
                  " would have been ", .num(100 * r$risk, 1L), "%", collapse = NULL))
  for (i in seq_len(nrow(cc))) {
    sig <- if (x@B > 0L && is.finite(cc$lower[i])) {
      if (cc$lower[i] > 0 || cc$upper[i] < 0) " (the interval excludes no difference)" else " (the interval includes no difference)"
    } else ""
    out <- c(out, paste0("- ", cc$level[i], " vs ", cc$reference[i], ": ", if (cc$difference[i] >= 0) "+" else "",
                         .num(100 * cc$difference[i], 1L), " percentage points", if (x@B > 0L) paste0(" (", lv, " CI ",
                         .num(100 * cc$lower[i], 1L), " to ", .num(100 * cc$upper[i], 1L), ")") else "", sig, "."))
  }
  if (identical(x@method, "both")) {
    g <- x@risk[x@risk$method == "gformula", ]; w <- x@risk[x@risk$method == "ipw", ]
    dmax <- max(abs(g$risk - w$risk[match(paste(g$level, g$time), paste(w$level, w$time))]), na.rm = TRUE)
    out <- c(out, paste0("- The two methods differ by at most ", .num(100 * dmax, 1L), " percentage points",
                         if (dmax < 0.02) ": the conclusion does not hinge on the outcome or the treatment model."
                         else ": they rely on different models (outcome vs treatment); a large difference suggests that one is misspecified."))
  }
  if (!is.null(x@ps)) {
    wt <- x@ps$weights
    out <- c(out, "", "Treatment model (IPW): balance and weights",
             .table_lines(data.frame(level = wt$level, n = wt$n, `effective n` = .num(wt$ess, 0L), `largest weight` = .num(wt$max_weight, 1L),
                                     check.names = FALSE)))
    if (!is.null(x@ps$smd)) {
      s <- x@ps$smd
      out <- c(out, paste0("  largest standardised mean difference: ", .num(max(abs(s$before), na.rm = TRUE), 2L), " before, ",
                           .num(max(abs(s$after), na.rm = TRUE), 2L), " after weighting (below 0.1 is usually taken as balanced)"))
      if (detail) out <- c(out, .table_lines(data.frame(covariate = s$covariate, level = s$level, before = .num(s$before, 3L), after = .num(s$after, 3L))))
    }
    pr <- x@ps$ps_range
    if (any(pr[, 1L] < 0.025 | pr[, 2L] > 0.975))
      out <- c(out, "  Note: some propensity scores are near 0 or 1 (positivity is doubtful); consider trim = 0.01.")
  }
  out <- c(out, "", "Assumptions:",
           "- No unmeasured confounding: the covariates of the model include every common cause of treatment and outcome.",
           "- Positivity: every subject could have received every treatment.",
           "- A correct outcome model (g-formula) or treatment model (IPW); method = \"both\" compares the two.",
           if (x@kind == "competing") "- The contrast is on the cumulative incidence (total effect through all causes), not on a rate.")
  c(out, "", .method_lines(.method_keys(x)))
}

S7::method(print, AdjustedRisk) <- function(x, ...) { cat(.ar_lines(x, FALSE), sep = "\n"); invisible(x) }
S7::method(summary, AdjustedRisk) <- function(object, ...) { cat(.ar_lines(object, TRUE), sep = "\n"); invisible(object) }

S7::method(as.data.frame, AdjustedRisk) <- function(x, row.names = NULL, optional = FALSE, what = c("risk", "contrast", "curve", "balance"), ...) {
  what <- match.arg(what)
  switch(what, risk = x@risk, contrast = x@contrast, curve = x@curve,
         balance = { if (is.null(x@ps)) stop("no treatment model: use method = \"ipw\" or \"both\"", call. = FALSE); x@ps$smd })
}

S7::method(plot, AdjustedRisk) <- function(x, method = NULL, col = NULL, ...) {
  ms <- method %||% unique(x@curve$method)
  col <- col %||% c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442", "#999999", "#882255")
  op <- graphics::par(mfrow = c(1L, length(ms)), mar = c(4, 4.2, 2.6, 1)); on.exit(graphics::par(op))
  for (m in ms) {
    cv <- x@curve[x@curve$method == m, , drop = FALSE]
    yl <- c(0, max(c(cv$risk, cv$upper), na.rm = TRUE))
    graphics::plot(NA, xlim = c(0, max(cv$time)), ylim = yl, xlab = "Time",
                   ylab = if (x@kind == "single") "Adjusted risk" else paste0("Adjusted cumulative incidence, cause ", x@cause), ...)
    for (i in seq_along(x@levels)) {
      d <- cv[cv$level == x@levels[i], , drop = FALSE]
      tt <- c(0, d$time); rr <- c(0, d$risk)
      if (all(is.finite(d$lower))) graphics::polygon(c(tt, rev(tt)), c(0, d$lower, rev(c(0, d$upper))),
                                                     col = grDevices::adjustcolor(col[i], 0.15), border = NA)
      graphics::lines(tt, rr, lwd = 2, col = col[i])
    }
    graphics::abline(v = x@times, lty = 3, col = "grey70")
    graphics::title(.ar_mname(m), cex.main = 0.9)
    graphics::legend("topleft", paste0(x@treatment, " = ", x@levels), col = col[seq_along(x@levels)], lwd = 2, bty = "n", cex = 0.8)
  }
  invisible(x)
}

# report(): see zzz-methods.R (the generic is defined in scan.R, which is collated later)
