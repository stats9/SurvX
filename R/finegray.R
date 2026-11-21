#' Fine-Gray regression for the cumulative incidence of a competing risk
#'
#' @description
#' Fits the proportional subdistribution hazards model of Fine and Gray (1999): the effect of
#' covariates on the \emph{cumulative incidence} of one cause (the probability of having that event
#' by time \eqn{t} when other events compete with it). Inverse-probability-of-censoring weights, the
#' robust variance of Fine and Gray, and delayed entry (Geskus 2011) are supported.
#'
#' \strong{Use it when} there are several event types (status 0 = censored, 1, 2, ... = causes) and
#' the question is how covariates change the \emph{risk} (absolute probability) of one cause, for
#' prognosis or prediction.
#' \strong{Use something else when}
#' \itemize{
#'   \item the question is about the mechanism, the rate at which the cause occurs among those still
#'     event-free: \code{\link{csc}} (cause-specific hazards). Reporting both is recommended
#'     (Latouche et al. 2013);
#'   \item a parametric cumulative incidence is wanted, a proportional odds link, or estimates for
#'     all causes jointly: \code{\link{directcif}};
#'   \item there is only one event type: \code{\link{cox}};
#'   \item no covariates, only groups: \code{\link{aj}} and \code{\link{gray_test}}.
#' }
#' @param formula a formula with \code{\link{Event}()} on the left and the covariates on the right:
#'   \tabular{ll}{
#'   \code{Event(time, status)} \tab one row per subject, right-censored \cr
#'   \code{Event(time, status, entry = entry)} \tab delayed entry (left truncation) \cr
#'   }
#'   \code{status} is 0 for censored and 1, 2, ... for the event types (a factor with the first
#'   level meaning censored also works, see \code{\link{Event}}). Covariates as in
#'   \code{\link[stats]{lm}}: numeric, factors, interactions, transformations, splines.
#'   \code{strata()} and \code{cluster()} are not available in \code{fg()}.
#' @param data a data frame with the variables of the formula. Rows with missing values in any of
#'   them are removed.
#' @param cause the event type of interest (default 1). All other non-zero codes are competing
#'   events.
#' @param variance \code{"robust"} (default): the sandwich variance of Fine and Gray (1999,
#'   Section 4), which accounts for the estimated censoring weights; use it for inference without
#'   delayed entry. \code{"bootstrap"}: the variance of the estimates over \code{B} bootstrap samples
#'   of the subjects; recommended with delayed entry, where the robust variance treats the weights as
#'   known and is too small (see Validation). \code{"model"}: the inverse information of the
#'   weighted partial likelihood, kept for comparison; not valid for inference in general. The
#'   robust and model-based variances are always stored in the fit.
#' @param max_iter maximum number of Newton-Raphson iterations (default 50).
#' @param tol convergence tolerance on the relative change of the log pseudo-likelihood.
#' @param B \code{variance = "bootstrap"}: number of bootstrap samples (default 200; 500 or more for
#'   final results).
#' @param seed \code{variance = "bootstrap"}: seed of the resampling, so that the same call gives
#'   the same standard errors. The random number generator of the session is restored afterwards.
#' @details
#' \subsection{The model}{
#' The subdistribution hazard of cause 1 is \eqn{\lambda_1(t) = -d \log\{1 - F_1(t)\} / dt}, with
#' \eqn{F_1(t) = P(T \le t, \mathrm{cause} = 1)} the cumulative incidence. The model is
#' \deqn{\lambda_1(t \mid x) = \lambda_{10}(t) \exp(\beta^\top x), \qquad
#'       F_1(t \mid x) = 1 - \exp\{-\Lambda_{10}(t) e^{\beta^\top x}\}.}
#' Subjects who have a competing event stay in the risk set after their event, with a weight that
#' decreases as the probability of still being under observation decreases.
#' }
#' \subsection{Estimation}{
#' Without delayed entry: the weighted partial likelihood of Fine and Gray (1999) with weights
#' \eqn{\hat G(t-)/\hat G(T_i-)}, where \eqn{\hat G} is the Kaplan-Meier estimate of the censoring
#' distribution; Breslow handling of tied event times. The computation uses a forward-backward scan
#' in linear time (as proposed by Kawaguchi et al. 2021), which is why \code{fg()} is fast on large
#' data. The robust variance is the sandwich of Fine and Gray (1999), with the term for the estimated
#' censoring distribution; ties between censoring and event times follow \code{cmprsk::crr}.
#'
#' With delayed entry: the weights of Geskus (2011),
#' \eqn{\hat G(t-)\hat H(t-)/\{\hat G(T_i-)\hat H(T_i-)\}}, with \eqn{\hat G} the product-limit
#' estimate of the censoring distribution among subjects under observation and \eqn{\hat H} the
#' reverse-time product-limit estimate of the entry distribution. The model is then a weighted Cox
#' model on (start, stop] rows cut where the weights change, and the robust variance is the sandwich
#' over subjects (Lin and Wei 1989) with the weights taken as known, as in
#' \code{survival::finegray} followed by \code{coxph(..., cluster = id)}.
#' }
#' @section Standard errors with delayed entry:
#' With delayed entry the robust variance treats the Geskus weights as known. In the Monte Carlo
#' study of the package it was 5 to 11\% too small, with 95\% intervals covering 91 to 95\%. The
#' bootstrap of the subjects (\code{variance = "bootstrap"}; Efron and Tibshirani 1993) re-estimates
#' the weights in every sample and so includes their uncertainty, at the cost of \code{B} refits.
#' @section Interpreting the coefficients:
#' \itemize{
#'   \item The printed effect is the \strong{subdistribution hazard ratio} (sHR, \eqn{e^\beta}). Its
#'     direction is directly interpretable: \strong{sHR > 1 means a higher cumulative incidence of
#'     the cause at every time point}, sHR < 1 a lower one, for a one-unit increase of a numeric
#'     covariate or a factor level against the reference, other covariates fixed.
#'   \item Its \emph{size} is not a ratio of risks, and it is not the effect on the rate of the
#'     event among those still event-free (that is the cause-specific HR of \code{\link{csc}}). An
#'     sHR of 1.5 does not mean a 50\% higher risk. Report absolute risks next to it:
#'     \code{predict(fit, newdata, times)} gives \eqn{F_1(t \mid x)} (Austin and Fine 2017).
#'   \item A covariate can raise the cumulative incidence of cause 1 without any effect on its
#'     cause-specific hazard, only by lowering the hazard of the competing cause (subjects then live
#'     longer and have more time to experience cause 1). Fitting \code{\link{csc}} as well shows
#'     which mechanism is at work.
#'   \item Units, interactions and factors are read as in \code{\link{cox}}: the sHR for \eqn{k}
#'     units is \eqn{sHR^k}; with \code{a * b} the effect of \code{a} applies at \code{b = 0}.
#' }
#' A sentence for a report: \emph{"In a Fine-Gray model, treatment was associated with a lower
#' cumulative incidence of relapse (sHR 0.72, 95\% CI 0.58 to 0.90); the predicted 5-year
#' incidence was 24\% with and 31\% without treatment for a 60-year-old patient."}
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{Proportional subdistribution hazards}{the sHR is constant over time. Compare
#'     \code{plot(fit, newdata)} with the non-parametric \code{\link{aj}} curves of the same groups,
#'     or fit the model on early and late follow-up separately (\code{\link{cif_plot}} with
#'     \code{type = "unadjusted"} overlays them).}
#'   \item{Censoring independent of the event and of the covariates}{the Kaplan-Meier weights assume
#'     censoring does not depend on the covariates. When it does strongly, results can be biased.}
#'   \item{Coherence across causes}{Fine-Gray models fitted separately for each cause need not give
#'     cumulative incidences that sum to at most 1. When predictions for all causes are needed
#'     together, \code{\link{directcif}} (\code{estimation = "joint"}) models them jointly.}
#' }
#' @section Messages and common problems:
#' \describe{
#'   \item{\code{WARNING: the fit did NOT converge}}{usually a covariate with no events of the cause
#'     in one of its levels (separation) or collinear covariates.}
#'   \item{\code{no events of the requested type}}{\code{cause} does not occur in the data; check the
#'     codes with \code{table(status)}.}
#'   \item{\code{strata() and cluster() are available in cox() and csc()}}{\code{fg()} has no
#'     stratified or clustered form; include the variable as a covariate.}
#' }
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}; 82 simulated and real data sets, up to \eqn{n = 20000}): against
#' \code{cmprsk::crr} 2.2, coefficients to \eqn{10^{-7}}, the robust variance to \eqn{10^{-9}}
#' (including heavily tied data) and predicted cumulative incidences (\code{predict.crr}) to
#' \eqn{2 \times 10^{-8}}; with delayed entry, against \code{survival::finegray} + \code{coxph}, robust
#' variances to \eqn{2 \times 10^{-7}} and coefficients to \eqn{2 \times 10^{-5}}, where the remaining
#' difference comes from \code{coxph} stopping before the maximum (its log-likelihood is higher at the
#' \code{fg()} estimate). In a Monte Carlo study (1000 replications, \eqn{n} = 200, 500, 1000) the
#' robust standard errors without delayed entry matched the empirical standard deviation (ratio 0.96
#' to 1.02, coverage 94.0 to 95.6\%). With delayed entry they were 5 to 11\% too small (coverage 91
#' to 95\%), because the weights are treated as known, as in \pkg{survival}; with small samples and
#' heavy truncation, interpret the intervals with this in mind. On \eqn{n = 10000}, \code{fg()} takes
#' about 0.01 s and \code{crr} about 150 s.
#' @return
#' An S7 object of class \code{FineGrayFit}. Main slots: \code{coefficients}, \code{vcov} (robust, bootstrap
#' or model, as requested), \code{vcov_robust}, \code{vcov_model}, \code{vcov_boot}, \code{loglik} (log
#' pseudo-likelihood), \code{loglik_null}, \code{cause}, \code{n}, \code{nevent}, \code{converged}.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(fit)}, \code{summary(fit)} \tab sHRs with intervals and p-values and the
#'     \code{Reading:} lines \cr
#'   \code{\link{estimates}(fit)} \tab the table as a data frame \cr
#'   \code{coef}, \code{vcov}, \code{confint}, \code{logLik}, \code{nobs} \tab extractors \cr
#'   \code{predict(fit, newdata, times)} \tab cumulative incidence of the cause, one row per row of
#'     \code{newdata}, one column per time \cr
#'   \code{plot(fit, ...)}, \code{\link{cif_plot}}, \code{\link{cif_data}} \tab cumulative
#'     incidence curves, with bootstrap intervals and adjusted (marginal) curves \cr
#'   \code{\link{brier_score}}, \code{\link{auc_t}}, \code{\link{compare_models}} \tab predictive
#'     performance \cr
#' }
#' @references
#' \emph{Model, weights and robust variance:}
#' Fine JP, Gray RJ (1999). A proportional hazards model for the subdistribution of a competing
#'   risk. Journal of the American Statistical Association 94:496-509.
#'
#' \emph{Sandwich variance in general:}
#' White H (1982). Maximum likelihood estimation of misspecified models. Econometrica 50:1-25.
#'
#' \emph{Delayed entry:}
#' Geskus RB (2011). Cause-specific cumulative incidence estimation and the Fine and Gray model under
#'   both left truncation and right censoring. Biometrics 67:39-49.
#' Lin DY, Wei LJ (1989). The robust inference for the Cox proportional hazards model. Journal of
#'   the American Statistical Association 84:1074-1078.
#'
#' \emph{Bootstrap variance:}
#' Efron B, Tibshirani RJ (1993). An Introduction to the Bootstrap. Chapman and Hall, New York.
#'
#' \emph{Computation:}
#' Kawaguchi ES, Shen JI, Suchard MA, Li G (2021). Scalable algorithms for large competing risks
#'   data. Journal of Computational and Graphical Statistics 30:685-693.
#'
#' \emph{Interpretation and reporting:}
#' Latouche A, Allignol A, Beyersmann J, Labopin M, Fine JP (2013). A competing risks analysis should
#'   report results on all cause-specific hazards and cumulative incidence functions. Journal of
#'   Clinical Epidemiology 66:648-653.
#' Austin PC, Fine JP (2017). Practical recommendations for reporting Fine-Gray model analyses for
#'   competing risk data. Statistics in Medicine 36:4391-4400.
#'
#' The references of a particular fit: \code{sx_references(fit)}; all references behind
#' \code{fg()}: \code{sx_references("fg")}.
#' @seealso
#' \code{\link{csc}} (cause-specific hazards), \code{\link{directcif}} (parametric cumulative
#' incidence), \code{\link{aj}} and \code{\link{gray_test}} (non-parametric), \code{\link{cif_plot}},
#' \code{\link{predict.SurvX}}, \code{\link{Event}}.
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- data.frame(treat = factor(rbinom(n, 1, 0.5), labels = c("control", "new")),
#'                 age = round(rnorm(n, 60, 10)))
#' t1 <- rexp(n, 0.10 * exp(-0.4 * (d$treat == "new") + 0.02 * (d$age - 60)))   # relapse
#' t2 <- rexp(n, 0.05 * exp(0.04 * (d$age - 60)))                               # death without relapse
#' cens <- runif(n, 0, 12)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0L, ifelse(t1 < t2, 1L, 2L))   # 0 censored, 1 relapse, 2 death
#' table(d$status)
#'
#' ## 1. subdistribution hazard ratios for relapse
#' fit <- fg(Event(time, status) ~ treat + I(age / 10), d, cause = 1)
#' fit
#'
#' ## 2. absolute risks: cumulative incidence of relapse at 2 and 5 years
#' nd <- data.frame(treat = c("control", "new"), age = 60)
#' predict(fit, nd, times = c(2, 5))
#' plot(fit, newdata = nd, conf.int = FALSE)
#'
#' ## 3. the mechanism: cause-specific hazards of both causes
#' csc(Event(time, status) ~ treat + I(age / 10), d)
#'
#' ## 4. delayed entry, with bootstrap standard errors
#' d$entry <- runif(n, 0, 0.3) * d$time
#' fg(Event(time, status, entry = entry) ~ treat + I(age / 10), d, variance = "bootstrap", B = 50)
#' @export
fg <- function(formula, data, cause = 1L, variance = c("robust", "model", "bootstrap"), max_iter = 50L, tol = 1e-10,
               B = 200L, seed = 1L) {
  cl <- match.call(); variance <- match.arg(variance)
  f <- .sx_frame(formula, data, entry_ok = TRUE, what = "fg()")
  n <- length(f$time); p <- ncol(f$X)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  fit <- if (!is.null(f$entry)) .fg_entry(cl, f, as.integer(cause), variance, max_iter, tol) else {
    r <- .C("SX_fg", as.integer(n), as.integer(p), f$time, f$status, .rm(f$X), as.integer(cause),
            as.integer(max_iter), as.double(tol), beta = double(p), var_model = double(p * p), var_robust = double(p * p),
            info_d = double(3), info_i = integer(4), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "fg")
    nm <- colnames(f$X)
    vm <- matrix(r$var_model, p, p, dimnames = list(nm, nm)); vr <- matrix(r$var_robust, p, p, dimnames = list(nm, nm))
    FineGrayFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
                coefficients = stats::setNames(r$beta, nm), vcov = if (variance == "model") vm else vr,
                loglik = r$info_d[2L], n = n, nevent = r$info_i[4L], iterations = r$info_i[1L],
                converged = as.logical(r$info_i[2L]), flags = r$info_i[3L],
                time = f$time, status = f$status, X = f$X, cause = as.integer(cause), variance = variance,
                vcov_model = vm, vcov_robust = vr, loglik_null = r$info_d[1L])
  }
  if (variance == "bootstrap") {
    vb <- .fg_bootstrap(f, as.integer(cause), as.integer(B), seed, max_iter, tol)
    fit@vcov_boot <- vb
    fit@vcov <- vb
  }
  fit
}

# nonparametric bootstrap of the subjects: the variance then includes the estimation of the
# censoring (and, with delayed entry, truncation) weights
.fg_bootstrap <- function(f, cause, B, seed, max_iter, tol) {
  if (B < 2L) stop("B must be at least 2", call. = FALSE)
  n <- length(f$time); p <- ncol(f$X); nm <- colnames(f$X)
  est <- .with_seed(seed, {
    out <- matrix(NA_real_, B, p)
    for (b in seq_len(B)) {
      i <- sample.int(n, n, replace = TRUE)
      out[b, ] <- tryCatch(.fg_coef(f$time[i], f$status[i], f$X[i, , drop = FALSE], if (is.null(f$entry)) NULL else f$entry[i],
                                    cause, max_iter, tol), error = function(e) rep(NA_real_, p))
    }
    out
  })
  ok <- stats::complete.cases(est)
  if (sum(ok) < max(10L, B / 2)) stop("the bootstrap failed: only ", sum(ok), " of ", B, " samples could be fitted", call. = FALSE)
  if (sum(ok) < B) warning(B - sum(ok), " of ", B, " bootstrap samples could not be fitted and were left out", call. = FALSE)
  v <- stats::cov(est[ok, , drop = FALSE])
  structure(matrix(v, p, p, dimnames = list(nm, nm)), B = sum(ok))
}

# coefficients of one Fine-Gray fit (NA when it does not converge)
.fg_coef <- function(time, status, X, entry, cause, max_iter, tol) {
  p <- ncol(X)
  if (!is.null(entry)) {
    e <- .fg_expand(time, status, entry, cause)
    r <- .cox_fit(e$stop, e$status, X[e$id, , drop = FALSE], e$w, "breslow", max_iter, tol, e$start, NULL)
    return(if (isTRUE(r$converged)) unname(r$coefficients) else rep(NA_real_, p))
  }
  r <- .C("SX_fg", as.integer(length(time)), as.integer(p), as.double(time), as.integer(status), .rm(X), as.integer(cause),
          as.integer(max_iter), as.double(tol), beta = double(p), var_model = double(p * p), var_robust = double(p * p),
          info_d = double(3), info_i = integer(4), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  if (r$rc != 0L || !as.logical(r$info_i[2L])) rep(NA_real_, p) else r$beta
}

# Fine-Gray data with delayed entry (Geskus 2011): subjects with a competing event stay in the risk
# set after their event time T_i with weight w_i(t) = G(t-) H(t-) / {G(T_i-) H(T_i-)}, where G is
# the product-limit estimate of the censoring distribution (risk set entry < s <= time) and H the
# reverse-time product-limit estimate of the entry distribution, H(t) = prod_{l > t} {1 - e(l)/R(l)}
# with R(l) = #{entry_j <= l < time_j}. The weights are constant between the jump points of G and H,
# so the extended rows are cut there. Without entry the weights are those of Fine & Gray (1999).
.fg_expand <- function(time, status, entry, cause) {
  n <- length(time)
  if (is.null(entry)) entry <- rep(0, n)
  se <- sort(entry); st <- sort(time)
  ct <- sort(unique(time[status == 0L]))
  Y <- findInterval(ct, se, left.open = TRUE) - findInterval(ct, st, left.open = TRUE)     # #{entry < c <= time}
  dc <- tabulate(match(time[status == 0L], ct), length(ct))
  cumG <- c(1, cumprod(1 - dc / Y))
  lt <- sort(unique(entry[entry > 0]))
  R <- findInterval(lt, se) - findInterval(lt, st)                                          # #{entry <= l < time}
  el <- tabulate(match(entry[entry > 0], lt), length(lt))
  Hrev <- c(rev(cumprod(rev(1 - el / R))), 1)                                               # prod over l_k, l_k+1, ...
  G_right <- function(a) cumG[findInterval(a, ct) + 1L]                                     # G(a)
  G_left <- function(a) cumG[findInterval(a, ct, left.open = TRUE) + 1L]                    # G(a-)
  H_after <- function(a) Hrev[findInterval(a, lt) + 1L]                                     # prod over l > a
  H_left <- function(a) Hrev[findInterval(a, lt, left.open = TRUE) + 1L]                    # H(a-): prod over l >= a
  jumps <- sort(unique(c(ct, lt)))
  tmax <- max(c(0, time[status == cause]))
  base <- data.frame(id = seq_len(n), start = entry, stop = time, status = as.integer(status == cause), w = 1)
  comp <- which(status != 0L & status != cause & time < tmax)
  ext <- lapply(comp, function(i) {
    Ti <- time[i]
    cuts <- c(Ti, jumps[jumps > Ti & jumps < tmax])
    stops <- c(cuts[-1L], tmax)
    den <- G_left(Ti) * H_left(Ti)
    w <- if (den > 0) G_right(cuts) * H_after(cuts) / den else rep(0, length(cuts))
    data.frame(id = i, start = cuts, stop = stops, status = 0L, w = w)
  })
  d <- do.call(rbind, c(list(base), ext))
  d[d$w > 0 & d$stop > d$start, , drop = FALSE]
}

.fg_entry <- function(cl, f, cause, variance, max_iter, tol) {
  e <- .fg_expand(f$time, f$status, f$entry, cause)
  Xe <- f$X[e$id, , drop = FALSE]
  r <- .cox_fit(e$stop, e$status, Xe, e$w, "breslow", max_iter, tol, e$start, NULL)
  vm <- r$vcov
  vr <- .cox_robust(r, e$stop, e$status, Xe, e$w, "breslow", e$start, NULL, e$id)
  FineGrayFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
              coefficients = r$coefficients, vcov = if (variance == "model") vm else vr,
              loglik = r$loglik, n = length(f$time), nevent = as.integer(sum(f$status == cause)), iterations = r$iterations,
              converged = r$converged, flags = r$flags,
              time = f$time, status = f$status, X = f$X, cause = cause, variance = variance,
              vcov_model = vm, vcov_robust = vr, loglik_null = r$loglik_null, entry = f$entry, expanded = e)
}

S7::method(predict, FineGrayFit) <- function(object, newdata = NULL, times, ...) {
  if (!is.null(object@entry)) {
    e <- object@expanded
    bh <- .cox_baseline_raw(object@coefficients, e$stop, e$status, object@X[e$id, , drop = FALSE], e$w, "breslow", e$start, NULL)
    H <- .step_at(bh$time, bh$cumhaz, times, 0)
  } else {
    n <- object@n; p <- length(object@coefficients)
    r <- .C("SX_fg_baseline", as.integer(n), as.integer(p), object@time, object@status, .rm(object@X), object@cause,
            as.double(object@coefficients), time = double(n), cumhaz = double(n), m = integer(1), rc = integer(1),
            NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "baseline subdistribution hazard")
    k <- seq_len(r$m)
    H <- .step_at(r$time[k], r$cumhaz[k], times, 0)
  }
  lp <- drop(.sx_newX(object, newdata) %*% object@coefficients)
  out <- 1 - exp(-outer(exp(lp), H))
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}

#' Cause-specific hazards regression for competing risks
#'
#' @description
#' Fits one hazard model per cause, each describing the rate of that cause among subjects still free
#' of all events, and combines them into the cumulative incidence of every cause. With
#' \code{baseline = "cox"} (default) each cause gets a Cox model; with a parametric
#' \code{baseline} each cause gets a parametric hazard model as in \code{\link{parsurv}}.
#'
#' \strong{Use it when} there are competing events and the question is about the mechanism: how a
#' covariate changes the rate at which each cause occurs (aetiology), or when absolute risks of all
#' causes are needed that are coherent (they never add up to more than 1).
#' \strong{Use something else when}
#' \itemize{
#'   \item the effect of covariates on the cumulative incidence of one cause is wanted as a single
#'     number: \code{\link{fg}} (or \code{\link{directcif}}). Reporting both kinds of model is
#'     recommended (Latouche et al. 2013);
#'   \item there is a single event type: \code{\link{cox}} or \code{\link{parsurv}}.
#' }
#' @param formula \code{Event(time, status) ~ covariates}, status 0 = censored, 1, 2, ..., K =
#'   cause. The same covariates enter the model of every cause. \code{Event(time, status, entry =
#'   entry)} for delayed entry and, with a Cox baseline, for (start, stop] rows with time-varying
#'   covariates. With a Cox baseline the right-hand side may contain \code{strata(...)} (a baseline
#'   hazard per stratum for every cause) and \code{cluster(id)} (robust variance), as in
#'   \code{\link{cox}}.
#' @param data a data frame. Rows with missing values are removed.
#' @param ties Cox baseline: \code{"efron"} (default) or \code{"breslow"}, as in \code{\link{cox}}.
#' @param product_limit Cox baseline: how the probability of being free of all events is computed
#'   for the cumulative incidence. \code{TRUE} (default): as a product-limit
#'   \eqn{\prod (1 - d\Lambda)}, as \pkg{riskRegression} does; \code{FALSE}:
#'   \eqn{\exp\{-\sum_k \Lambda_k(t)\}}. The difference is small except with few subjects at risk.
#' @param baseline \code{"cox"}, or one of the distributions of \code{\link{parsurv}}
#'   (\code{"exponential"}, \code{"weibull"}, \code{"gompertz"}, \code{"lognormal"},
#'   \code{"loglogistic"}, \code{"gengamma"}, \code{"rp"}, \code{"pwexp"}). A vector gives one
#'   distribution per cause, for example \code{c("weibull", "gompertz")}. Cox and parametric
#'   baselines cannot be mixed.
#' @param scale parametric baselines: \code{"PH"} or \code{"AFT"}, one value or one per cause;
#'   default as in \code{\link{parsurv}}.
#' @param df, knots, cuts spline (\code{"rp"}) and piecewise exponential (\code{"pwexp"})
#'   baselines, as in \code{\link{parsurv}}; by default knots and cut points are chosen from the
#'   events of each cause.
#' @param max_iter, tol iteration limit and convergence tolerance.
#' @param robust Cox baseline: robust (sandwich) variance of Lin and Wei (1989) for every cause;
#'   always used with \code{cluster()}.
#' @details
#' \subsection{The model}{
#' The cause-specific hazard of cause \eqn{k},
#' \eqn{h_k(t) = \lim P(t \le T < t + dt, \mathrm{cause} = k \mid T \ge t)/dt}, is modelled as
#' \eqn{h_k(t \mid x) = h_{k0}(t) \exp(\beta_k^\top x)} (Prentice et al. 1978). Each cause is fitted
#' separately, with the subjects who had another cause censored at their event time: the likelihood
#' factorises over causes, so this is the maximum likelihood estimate of the joint model, and each
#' cause-specific Cox model is identical to \code{survival::coxph} on \code{Surv(time, status == k)}.
#' }
#' \subsection{Cumulative incidence}{
#' \deqn{F_k(t \mid x) = \int_0^t S(u- \mid x)\, h_k(u \mid x)\,du, \qquad
#'       S(t \mid x) = \exp\Big\{-\sum_j \Lambda_j(t \mid x)\Big\}}
#' (Kalbfleisch and Prentice 2002, Section 8.2). With a Cox baseline this is a sum over event times
#' with the Breslow baseline hazards (Benichou and Gail 1990), as in \code{riskRegression::CSC}; for
#' a subject in stratum \eqn{s}, the baseline hazards of that stratum are used. With a parametric
#' baseline the integral is computed numerically (20-point Gauss-Legendre on a fine partition of the
#' time axis), with delta-method standard errors (Hinchliffe and Lambert 2013). Either way the
#' cumulative incidences of all causes and the event-free probability add up to one.
#' }
#' @section Interpreting the coefficients:
#' \itemize{
#'   \item The printed effect for cause \eqn{k} is the \strong{cause-specific hazard ratio}: the
#'     factor by which a covariate multiplies the rate of cause \eqn{k} \emph{among subjects who are
#'     still free of every event}. It is read like the HR of \code{\link{cox}}.
#'   \item A cause-specific HR does not translate into the risk of the cause: a covariate that does
#'     not change the hazard of cause 1 but lowers that of cause 2 increases the cumulative incidence
#'     of cause 1 (subjects survive longer to experience it). Look at the effects on \emph{all}
#'     causes together, and at the predicted cumulative incidences:
#'     \code{predict(fit, newdata, times, cause)} or \code{plot(fit, newdata)}.
#'   \item In aetiological questions (does exposure X cause the disease?) the cause-specific HR is the
#'     quantity of interest; for prognosis (what is the risk of this patient?) the cumulative
#'     incidence is (Andersen et al. 2012; Latouche et al. 2013).
#' }
#' A sentence for a report: \emph{"Treatment reduced the rate of relapse (cause-specific HR 0.65,
#' 95\% CI 0.50 to 0.85) and did not change the rate of death without relapse (HR 1.05, 0.75 to
#' 1.47); the predicted 5-year cumulative incidence of relapse was 22\% versus 31\%."}
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{Proportional cause-specific hazards}{one Cox model per cause: check each with
#'     \code{ph_test(fit)}, which tests every cause.}
#'   \item{Parametric baselines}{\code{hazard_plot(fit)} compares the fitted and non-parametric
#'     cause-specific hazards; \code{AIC()} compares distributions (for all causes together).}
#'   \item{Independent censoring}{as in \code{\link{cox}}; the competing events themselves are
#'     not censoring, and need no independence assumption for the cause-specific hazards.}
#' }
#' @section Messages and common problems:
#' \describe{
#'   \item{\code{WARNING: the fit did NOT converge} (for a cause)}{usually a cause with very few
#'     events, or separation; report that cause with care or merge rare causes.}
#'   \item{\code{baseline: 'cox' cannot be mixed with parametric distributions}}{choose all Cox or all
#'     parametric baselines.}
#'   \item{\code{at least one covariate is needed}}{for no covariates use \code{\link{aj}}.}
#' }
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}; 82 data sets with and without ties, strata and delayed entry): coefficients
#' agree with \code{survival::coxph} for each cause to \eqn{6 \times 10^{-7}} and variances to
#' \eqn{4 \times 10^{-10}}; predicted cumulative incidences agree with \code{riskRegression::CSC} to
#' \eqn{4 \times 10^{-7}} (stratified: \eqn{5 \times 10^{-10}}). In a Monte Carlo study (1000
#' replications, \eqn{n} = 200 to 1000) the 95\% intervals covered 94.5 to 97.3\%. On
#' \eqn{n = 10000}, a fit takes about 0.01 s (\code{riskRegression::CSC}: 0.07 s).
#' @return
#' With \code{baseline = "cox"}: an S7 object of class \code{CSCFit}; with a parametric baseline:
#' class \code{ParCSCFit}. In both, \code{fit@models} holds the fit of each cause and \code{fit@K} the
#' number of causes; \code{coef(fit)} and \code{vcov(fit)} are lists with one element per cause;
#' \code{logLik()} and \code{AIC()} refer to all causes together.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(fit)}, \code{summary(fit)} \tab one table of cause-specific HRs per cause, with
#'     likelihood ratio tests and \code{Reading:} lines \cr
#'   \code{\link{estimates}(fit)} \tab all causes in one data frame, with a \code{cause} column \cr
#'   \code{predict(fit, newdata, times, cause)} \tab cumulative incidence of the cause; parametric
#'     fits also \code{type = "survival"} (free of all events) and \code{se = TRUE} \cr
#'   \code{plot(fit)}, \code{\link{cif_plot}}, \code{\link{cif_data}} \tab cumulative incidence
#'     curves of all causes, with bootstrap (Cox) or delta-method (parametric) intervals \cr
#'   \code{\link{ph_test}(fit)}, \code{residuals(fit)} \tab checks per cause \cr
#'   \code{\link{hazard_plot}(fit)} \tab parametric fits: fitted against non-parametric hazards \cr
#' }
#' @references
#' \emph{Cause-specific hazards:}
#' Prentice RL, Kalbfleisch JD, Peterson AV, Flournoy N, Farewell VT, Breslow NE (1978). The
#'   analysis of failure times in the presence of competing risks. Biometrics 34:541-554.
#' Kalbfleisch JD, Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd ed.
#'   Wiley, Chapter 8.
#'
#' \emph{Cox models per cause:}
#' Cox DR (1972). Regression models and life-tables. Journal of the Royal Statistical Society,
#'   Series B 34:187-220.
#' Breslow NE (1974). Covariance analysis of censored survival data. Biometrics 30:89-99.
#' Efron B (1977). The efficiency of Cox's likelihood function for censored data. Journal of the
#'   American Statistical Association 72:557-565.
#' Andersen PK, Gill RD (1982). Cox's regression model for counting processes: a large sample study.
#'   Annals of Statistics 10:1100-1120.
#' Lin DY, Wei LJ (1989). The robust inference for the Cox proportional hazards model. Journal of
#'   the American Statistical Association 84:1074-1078.
#'
#' \emph{Absolute risks:}
#' Benichou J, Gail MH (1990). Estimates of absolute cause-specific risk in cohort studies.
#'   Biometrics 46:813-826.
#' Hinchliffe SR, Lambert PC (2013). Flexible parametric modelling of cause specific hazards to
#'   estimate cumulative incidence functions. BMC Medical Research Methodology 13:13.
#'
#' \emph{Parametric models:}
#' Lawless JF (2003). Statistical Models and Methods for Lifetime Data, 2nd ed. Wiley.
#'
#' \emph{Interpretation:}
#' Latouche A, Allignol A, Beyersmann J, Labopin M, Fine JP (2013). A competing risks analysis should
#'   report results on all cause-specific hazards and cumulative incidence functions. Journal of
#'   Clinical Epidemiology 66:648-653.
#' Andersen PK, Geskus RB, de Witte T, Putter H (2012). Competing risks in epidemiology:
#'   possibilities and pitfalls. International Journal of Epidemiology 41:861-870.
#'
#' All references: \code{sx_references("csc")}.
#' @seealso
#' \code{\link{fg}} and \code{\link{directcif}} (models for the cumulative incidence),
#' \code{\link{aj}}, \code{\link{cif_plot}}, \code{\link{cox}}, \code{\link{parsurv}},
#' \code{\link{predict.SurvX}}.
#' @examples
#' set.seed(2)
#' n <- 300
#' d <- data.frame(treat = rbinom(n, 1, 0.5), age = round(rnorm(n, 60, 10)))
#' t1 <- rweibull(n, 1.3, 3 * exp(0.4 * d$treat))        # relapse, delayed by treatment
#' t2 <- rweibull(n, 0.9, 5 * exp(-0.02 * (d$age - 60)))  # death without relapse
#' cens <- runif(n, 0, 6)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#'
#' ## 1. cause-specific Cox models for both causes
#' fit <- csc(Event(time, status) ~ treat + age, d)
#' fit
#'
#' ## 2. absolute risk of relapse at 1 and 3 years
#' nd <- data.frame(treat = 0:1, age = 60)
#' predict(fit, nd, times = c(1, 3), cause = 1)
#'
#' ## 3. parametric cause-specific hazards, one distribution per cause
#' fw <- csc(Event(time, status) ~ treat + age, d, baseline = c("weibull", "weibull"))
#' predict(fw, nd, times = 3, cause = 1, se = TRUE)
#' AIC(fw, csc(Event(time, status) ~ treat + age, d, baseline = "lognormal"))
#' @export
csc <- function(formula, data, ties = c("efron", "breslow"), product_limit = TRUE, baseline = "cox", scale = NULL,
                df = 3L, knots = NULL, cuts = NULL, max_iter = 200L, tol = 1e-10, robust = FALSE) {
  cl <- match.call(); ties <- match.arg(ties)
  parametric <- any(baseline != "cox")
  if (parametric && any(baseline == "cox")) stop("baseline: 'cox' cannot be mixed with parametric distributions", call. = FALSE)
  f <- .sx_frame(formula, data, entry_ok = TRUE, what = "csc() with a parametric baseline", specials_ok = !parametric)
  if (parametric) return(.pcsc(cl, f, baseline, scale, max_iter, tol, df, knots, cuts))
  K <- max(f$status)
  if (ncol(f$X) < 1L) stop("at least one covariate is needed", call. = FALSE)
  robust <- isTRUE(robust) || !is.null(f$cluster)
  models <- lapply(seq_len(K), function(k) {
    st <- as.integer(f$status == k)
    m <- .cox_fit(f$time, st, f$X, NULL, ties, 50L, tol, f$entry, f$strata)
    m$vcov_model <- m$vcov
    if (robust) m$vcov <- .cox_robust(m, f$time, st, f$X, NULL, ties, f$entry, f$strata, f$cluster)
    m
  })
  names(models) <- paste0("cause", seq_len(K))
  CSCFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts, models = models, K = as.integer(K),
         ties = ties, product_limit = isTRUE(product_limit), n = length(f$time), time = f$time, status = f$status, X = f$X,
         entry = f$entry, strata = f$strata, strata_levels = f$strata_levels, strata_vars = f$strata_vars,
         cluster = f$cluster, robust = robust)
}

# cumulative incidence from cause-specific Cox models with strata or (entry, time] data:
# F_k(t) = sum_{u <= t} S(u-) dLambda_k(u | x), S from all causes (product-limit or exponential),
# with the baseline hazards of the subject's stratum (Benichou & Gail 1990; Kalbfleisch & Prentice 2002, Sec. 8.2)
.csc_cif_general <- function(object, Xn, sn, times, cause) {
  K <- object@K
  bh <- lapply(seq_len(K), function(k) .csc_baseline(object, k))
  B <- lapply(object@models, function(m) m$coefficients)
  out <- matrix(0, nrow(Xn), length(times), dimnames = list(NULL, paste0("t=", times)))
  for (i in seq_len(nrow(Xn))) {
    s <- if (is.null(sn)) 0L else sn[i]
    inc <- lapply(seq_len(K), function(k) {
      b <- bh[[k]][bh[[k]]$stratum == s, , drop = FALSE]
      data.frame(time = b$time, dH = diff(c(0, b$cumhaz)) * exp(sum(Xn[i, ] * B[[k]])))
    })
    grid <- sort(unique(unlist(lapply(inc, `[[`, "time"))))
    if (!length(grid)) next
    dH <- vapply(inc, function(d) { v <- numeric(length(grid)); v[match(d$time, grid)] <- d$dH; v }, numeric(length(grid)))
    dH <- matrix(dH, length(grid), K)
    tot <- rowSums(dH)
    Sstep <- if (object@product_limit) pmax(0, 1 - tot) else exp(-tot)
    Sminus <- c(1, cumprod(Sstep))[seq_along(grid)]
    F <- cumsum(Sminus * dH[, cause])
    out[i, ] <- .step_at(grid, F, times, 0)
  }
  out
}

S7::method(predict, CSCFit) <- function(object, newdata = NULL, times, cause = 1L, ...) {
  Xn <- .sx_newX(object, newdata)
  if (!is.null(object@strata) || !is.null(object@entry))
    return(.csc_cif_general(object, Xn, .sx_newstrata(object, newdata), times, as.integer(cause)))
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
