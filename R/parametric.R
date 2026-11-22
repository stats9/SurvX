# Parametric survival regression, parametric cause-specific hazards and the non-parametric hazard.
# The likelihoods, fits, predictions and numerical integration are in the C++ core
# (src/survcore/parametric.hpp, src/survcore/hazard.hpp); this file only prepares data and output.

.par_dists <- c(exponential = 0L, weibull = 1L, gompertz = 2L, lognormal = 3L, loglogistic = 4L, gengamma = 5L, rp = 6L, pwexp = 7L)
.par_label <- c(exponential = "exponential", weibull = "Weibull", gompertz = "Gompertz",
                lognormal = "log-normal", loglogistic = "log-logistic", gengamma = "generalized gamma", rp = "Royston-Parmar spline",
                pwexp = "piecewise exponential")
.par_types <- c(survival = 0L, hazard = 1L, cumhaz = 2L, quantile = 3L, rmst = 4L)

# the covariate scale of a distribution: PH or AFT, checked against what the distribution supports
.par_scale <- function(dist, scale) {
  allowed <- switch(dist, exponential = c("PH", "AFT"), weibull = c("PH", "AFT"), gompertz = "PH",
                    lognormal = "AFT", loglogistic = "AFT", gengamma = "AFT", rp = "PH", pwexp = "PH")
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

.shape_name <- function(dist, scale, knots = numeric(0)) {
  if (dist == "gengamma") return(c("log(scale)", "Q"))
  if (dist == "rp") return(paste0("gamma", seq_len(length(knots) - 1L)))
  if (dist == "pwexp") {
    if (!length(knots)) return(character(0))
    lo <- .sx_signif(knots, 4L); hi <- c(.sx_signif(knots[-1L], 4L), Inf)
    return(paste0("log rate ratio (", lo, ", ", hi, "]"))
  }
  if (dist == "exponential") character(0) else if (scale == "AFT") "log(scale)" else if (dist == "gompertz") "shape" else "log(shape)"
}

.nshape <- function(dist, knots) {
  switch(dist, exponential = 0L, gengamma = 2L, rp = length(knots) - 1L, pwexp = length(knots), 1L)
}

# knots of the spline model (log time, boundary knots at the smallest and largest log event time, internal
# knots at equally spaced quantiles of the log event times, as in flexsurv::flexsurvspline) and cut points of
# the piecewise exponential model (default: quartiles of the event times)
.par_knots <- function(dist, time, status01, df, knots, cuts) {
  if (dist == "rp") {
    lt <- .sx_ew("log", time[status01 == 1L])
    if (length(unique(lt)) < 2L) stop("the spline model needs at least two distinct event times", call. = FALSE)
    if (!is.null(knots)) inner <- sort(.sx_ew("log", as.double(knots)))
    else {
      df <- as.integer(df)
      if (is.na(df) || df < 1L) stop("df must be a whole number of at least 1", call. = FALSE)
      inner <- if (df > 1L) .sx_quantile(lt, .sx_seq(0, 1, df + 1L)[-c(1L, df + 1L)]) else numeric(0)
    }
    b <- range(lt)
    if (any(inner <= b[1L] | inner >= b[2L]))
      stop("internal knots must lie strictly between the smallest and the largest event time", call. = FALSE)
    return(c(b[1L], inner, b[2L]))
  }
  if (dist == "pwexp") {
    if (!is.null(cuts)) {
      cs <- sort(unique(as.double(cuts)))
      if (any(!is.finite(cs) | cs <= 0)) stop("cuts must be positive and finite", call. = FALSE)
      return(cs)
    }
    return(unique(.sx_quantile(time[status01 == 1L], c(0.25, 0.5, 0.75))))
  }
  numeric(0)
}

# one fit in the C++ core; status01 is 0/1
.par_fit <- function(time, status01, entry, X, dist, scale, init = NULL, max_iter = 200L, tol = 1e-10, knots = numeric(0)) {
  n <- length(time); p <- ncol(X)
  P <- 1L + p + .nshape(dist, knots)
  r <- .C("SX_par_fit", as.integer(n), as.integer(p), .par_dists[[dist]], as.integer(scale == "AFT"),
          as.integer(length(knots)), if (length(knots)) as.double(knots) else double(1),
          as.integer(!is.null(entry)), if (is.null(entry)) double(1) else as.double(entry),
          as.double(time), as.integer(status01), .rm(X), as.integer(!is.null(init)),
          if (is.null(init)) double(P) else as.double(init), as.integer(max_iter), as.double(tol),
          theta = double(P), var = double(P * P), info_d = double(3), info_i = integer(4), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "parametric fit")
  nm <- c("(Intercept)", colnames(X), .shape_name(dist, scale, knots))
  V <- matrix(r$var, P, P, byrow = TRUE, dimnames = list(nm, nm))
  # a fit whose estimates or standard errors are not finite is never reported as converged
  ok <- all(is.finite(r$theta)) && all(is.finite(diag(V))) && all(diag(V) >= 0) && is.finite(r$info_d[2L])
  list(theta = stats::setNames(r$theta, nm), var = V,
       loglik = r$info_d[2L], iterations = r$info_i[1L], converged = as.logical(r$info_i[2L]) && ok,
       flags = r$info_i[3L], nevent = r$info_i[4L])
}

.parsurv_core <- function(cl, f, status01, dist, scale, init, max_iter, tol, label, df = 3L, knots = NULL, cuts = NULL, kn = NULL) {
  if (any(f$time <= 0)) stop("times must be positive for parametric models", call. = FALSE)
  if (is.null(kn)) kn <- .par_knots(dist, f$time, status01, df, knots, cuts)
  r <- .par_fit(f$time, status01, f$entry, f$X, dist, scale, init, max_iter, tol, kn)
  ln <- if (ncol(f$X)) tryCatch(.par_fit(f$time, status01, f$entry, f$X[, 0L, drop = FALSE], dist, scale, NULL, max_iter, tol, kn)$loglik,
                                error = function(e) NA_real_) else r$loglik
  ParSurvFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
             coefficients = r$theta, vcov = r$var, loglik = r$loglik, n = length(f$time),
             nevent = as.integer(r$nevent), iterations = as.integer(r$iterations), converged = isTRUE(r$converged),
             flags = as.integer(r$flags), time = f$time, status = as.integer(status01), X = f$X,
             dist = dist, scale = scale, dist_code = .par_dists[[dist]], aft = as.integer(scale == "AFT"),
             entry = f$entry, loglik_null = ln, event_label = label, knots = kn)
}

#' Parametric survival regression
#'
#' @description
#' Fits a fully parametric model for the time to a single event: exponential, Weibull, Gompertz,
#' log-normal, log-logistic, generalized gamma, piecewise exponential, or a flexible spline model
#' (Royston-Parmar), under right censoring and, optionally, delayed entry. Covariates act on the
#' hazard (proportional hazards, \code{scale = "PH"}) or on the time scale (accelerated failure
#' time, \code{scale = "AFT"}).
#'
#' \strong{Use it when} the shape of the hazard over time matters, smooth predicted curves are
#' wanted, survival must be extrapolated beyond follow-up (health economics), or time ratios
#' (how much longer or shorter the time to the event is) are the effect of interest.
#' \strong{Use something else when}
#' \itemize{
#'   \item only the covariate effects matter and no assumption on the baseline is wanted:
#'     \code{\link{cox}};
#'   \item there are competing events: \code{\link{csc}} with the same \code{baseline} (one
#'     parametric model per cause), or \code{\link{directcif}};
#'   \item part of the population will never have the event: \code{\link{cure}}.
#' }
#' @param formula \code{Event(time, status) ~ covariates} with status 0/1;
#'   \code{Event(time, status, entry = entry)} for delayed entry (left truncation);
#'   \code{Interval(left, right) ~ covariates} for interval-censored times (see the section
#'   Interval-censored data); \code{+ frailty(id)} for a shared gamma frailty of the clusters
#'   \code{id} (see \code{\link{frailty}}). Covariates as in \code{\link[stats]{lm}}.
#' @param data a data frame. Rows with missing values in the variables of the formula are removed.
#' @param baseline the distribution:
#'   \tabular{lll}{
#'   \code{"weibull"} (default) \tab PH or AFT \tab monotone hazard (rising, falling or constant) \cr
#'   \code{"exponential"} \tab PH or AFT \tab constant hazard \cr
#'   \code{"gompertz"} \tab PH \tab hazard rising or falling exponentially (ageing, mortality) \cr
#'   \code{"lognormal"} \tab AFT \tab hazard rising then falling \cr
#'   \code{"loglogistic"} \tab AFT \tab hazard rising then falling, closed forms \cr
#'   \code{"gengamma"} \tab AFT \tab generalized gamma: contains Weibull, log-normal and gamma \cr
#'   \code{"rp"} \tab PH \tab Royston-Parmar spline: any smooth shape, set by \code{df} \cr
#'   \code{"pwexp"} \tab PH \tab piecewise constant hazard, cut at \code{cuts} \cr
#'   }
#'   Compare candidates with \code{AIC()} and with \code{\link{hazard_plot}} against the
#'   non-parametric hazard.
#' @param scale \code{"PH"} or \code{"AFT"}. Default \code{"PH"} where it exists (exponential,
#'   Weibull, Gompertz, spline, piecewise exponential), else \code{"AFT"}. For the Weibull and
#'   exponential the two are the same model with different coefficients (\eqn{\beta_{PH} =
#'   -\beta_{AFT}/\sigma}).
#' @param df spline baseline: degrees of freedom, the number of internal knots plus one (default 3;
#'   \code{df = 1} is the Weibull model). Internal knots are placed at equally spaced quantiles of the
#'   log event times. Values from 2 to 5 are usually enough; compare them with \code{AIC()}.
#' @param knots spline baseline: internal knots on the time scale, instead of \code{df}. They must
#'   lie strictly between the smallest and the largest event time.
#' @param cuts piecewise exponential baseline: cut points of the time axis (default the quartiles
#'   of the event times).
#' @param init optional starting values (intercept, covariates, shape parameters), for difficult
#'   fits.
#' @param max_iter maximum number of Newton-Raphson iterations (default 200).
#' @param tol convergence tolerance on the relative change of the log-likelihood.
#' @details
#' Parametrisations, with \eqn{\eta = \beta_0 + z'\beta}:
#' \describe{
#'   \item{PH}{\eqn{H(t \mid z) = e^{\eta} H_0(t)} with \eqn{H_0(t) = t} (exponential),
#'     \eqn{t^{\kappa}} (Weibull; coefficient \code{log(shape)} \eqn{= \log\kappa}) or
#'     \eqn{(e^{a t} - 1)/a} (Gompertz; coefficient \code{shape} \eqn{= a}, any sign; \eqn{a < 0}
#'     gives a fraction that never fails). \eqn{e^{\beta_j}} is a hazard ratio. These are the
#'     \code{"exp"}, \code{"weibullPH"} and \code{"gompertz"} models of \pkg{flexsurv}.}
#'   \item{AFT}{\eqn{\log T = \eta + \sigma W}, coefficient \code{log(scale)} \eqn{= \log\sigma},
#'     with \eqn{W} standard minimum extreme value (Weibull; exponential with \eqn{\sigma = 1}),
#'     standard normal (log-normal) or standard logistic (log-logistic). \eqn{e^{\beta_j}} is a time
#'     ratio (acceleration factor). These are the models of \code{survival::survreg}.}
#'   \item{generalized gamma (\code{"gengamma"})}{AFT, \eqn{\log T = \eta + \sigma W} where, for
#'     \eqn{Q \ne 0}, \eqn{Q^{-2} e^{Q W}} has a gamma distribution with shape \eqn{Q^{-2}} (the
#'     parametrisation of Prentice 1974 and of \pkg{flexsurv}); coefficients \code{log(scale)}
#'     \eqn{= \log\sigma} and \code{Q}. It contains the log-normal (\eqn{Q = 0}), Weibull
#'     (\eqn{Q = 1}) and gamma (\eqn{Q = \sigma}) distributions, so a Wald or likelihood ratio test of
#'     \eqn{Q} chooses between them (Cox et al. 2007). The log density is evaluated in a form that is
#'     free of cancellation near \eqn{Q = 0}; the survivor function uses the incomplete gamma function
#'     for \eqn{|Q| \ge 0.01} and numerical integration of the density for smaller \eqn{|Q|}, so the
#'     likelihood is smooth across the log-normal limit. Derivatives in \eqn{(\mu, \log\sigma, Q)}
#'     are central differences of the exact contribution of each subject.}
#'   \item{spline (\code{"rp"})}{\eqn{\log H(t \mid z) = s(\log t) + z'\beta} with a restricted
#'     cubic spline \eqn{s(x) = \gamma_0 + \gamma_1 x + \sum_j \gamma_{j+1} v_j(x)}, linear beyond
#'     the boundary knots at the smallest and largest log event time (Royston and Parmar 2002); the
#'     basis of \code{flexsurv::flexsurvspline}. \code{(Intercept)} is \eqn{\gamma_0}; the
#'     coefficients \code{gamma1}, \code{gamma2}, ... are the shape parameters; \eqn{e^{\beta_j}} is a
#'     hazard ratio. With \code{df = 1} it is the Weibull model.}
#'   \item{piecewise exponential (\code{"pwexp"})}{hazard \eqn{e^{\eta + a_j}} on the interval
#'     \eqn{(c_{j-1}, c_j]}, \eqn{a_1 = 0} (Friedman 1982); the coefficients \code{log rate ratio}
#'     compare each later interval with the first. The maximum likelihood estimates equal those of a
#'     Poisson regression of the events on the split follow-up time.}
#' }
#' The log-likelihood is \eqn{\sum_i d_i \log h(t_i) + \log S(t_i) - \log S(e_i)} with entry time
#' \eqn{e_i} (zero without delayed entry). Gradient and Hessian are analytic and the maximum is found
#' by Newton-Raphson with step halving (and a Levenberg-Marquardt ridge where the information is not
#' positive definite), starting from the exponential fit. Standard errors are from the observed
#' information.
#'
#' \code{print()} shows the effects on the ratio scale with a one-line reading of every effect with
#' \eqn{p < 0.05}; \code{summary()} adds the baseline distribution on its natural scale (rate or time
#' scale, shape) and notes on interpretation.
#'
#' Estimation is by maximum likelihood with analytic first and second derivatives (central
#' differences for the generalized gamma), by Newton-Raphson with step halving. Standard errors are
#' the inverse of the observed information. With delayed entry, each subject contributes
#' \eqn{\log f(t)} or \eqn{\log S(t)} minus \eqn{\log S(entry)} (Klein and Moeschberger 2003,
#' Section 3.4).
#' @section Interpreting the coefficients:
#' \describe{
#'   \item{PH models}{\eqn{e^\beta} is a \strong{hazard ratio}, read as in \code{\link{cox}}: HR 0.70
#'     means a 30\% lower event rate at every time.}
#'   \item{AFT models}{\eqn{e^\beta} is a \strong{time ratio} (acceleration factor): TR 1.40 means
#'     that the time to the event is stretched by 40\% (every quantile of survival time, the median
#'     included, is 1.4 times longer); TR < 1 means shorter times and worse prognosis. Note the
#'     opposite direction to the HR: a protective factor has HR < 1 but TR > 1.}
#'   \item{Shape parameters}{(\code{log(shape)}, \code{log(scale)}, \code{Q}, \code{gamma1}, ...)
#'     describe the baseline and are rarely interpreted on their own. For the Weibull, a shape
#'     \eqn{\kappa > 1} means a rising hazard; for the generalized gamma, \eqn{Q} near 1 means that
#'     the Weibull fits, near 0 the log-normal.}
#' }
#' Absolute quantities usually speak more clearly than coefficients:
#' \code{predict(fit, newdata, times)} (survival), \code{type = "quantile", p = 0.5} (median time),
#' \code{type = "rmst"} (restricted mean). A sentence for a report: \emph{"In a Weibull
#' accelerated failure time model, treatment prolonged the time to progression by 40\% (time ratio
#' 1.40, 95\% CI 1.15 to 1.70); predicted median 14.2 versus 10.1 months."}
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{The distribution fits}{compare the fitted hazard and cumulative hazard with the
#'     non-parametric estimates: \code{hazard_plot(fit)}; \code{\link{dist_plot}} draws the
#'     linearising plots of each distribution (for example \eqn{\log H} against \eqn{\log t}, a straight
#'     line for the Weibull); compare candidate distributions with \code{AIC()}.}
#'   \item{PH or AFT structure}{for PH models, \code{\link{ph_test}} on a \code{cox} fit with the
#'     same covariates checks proportional hazards.}
#'   \item{Extrapolation}{beyond follow-up, the curve depends on the distribution only; compare the
#'     extrapolations of several well-fitting distributions before relying on one.}
#' }
#' @section Messages and common problems:
#' \describe{
#'   \item{\code{WARNING: the fit did NOT converge}}{often a distribution that cannot describe the
#'     data (the generalized gamma on small samples), separation, or a covariate on a very large or
#'     small scale (rescale it). Try another baseline, rescale covariates, or give \code{init}.}
#'   \item{\code{the spline model needs at least two distinct event times}}{too few events for
#'     \code{"rp"}.}
#'   \item{\code{internal knots must lie strictly between the smallest and the largest event time}}{
#'     choose \code{knots} inside the range of the event times, or a smaller \code{df}.}
#'   \item{\code{parsurv() is for one event (status 0/1)}}{use \code{csc(..., baseline = )} for
#'     competing risks.}
#' }
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}; 82 data sets): AFT models agree with \code{survival::survreg} to \eqn{3 \times
#' 10^{-9}} for coefficients, \eqn{2 \times 10^{-10}} for standard errors and \eqn{3 \times 10^{-11}} for
#' the log-likelihood; all distributions agree with \code{flexsurv::flexsurvreg} and
#' \code{flexsurvspline} within the tolerance of their optimisers (log-likelihood to a relative
#' \eqn{10^{-7}}, covariate coefficients and standard errors to about \eqn{10^{-4}}), and where they
#' differ more, SurvX reaches the higher log-likelihood. In a Monte Carlo study (1000 replications,
#' \eqn{n} = 200 to 1000, with and without delayed entry) the 95\% intervals covered 92.9 to 96.3\%.
#' On \eqn{n = 10000}, a Weibull fit takes about 0.02 s (\code{flexsurvreg}: 0.17 s), a
#' generalized gamma fit 0.36 s (2.8 s) and a spline fit 0.03 s (\code{flexsurvspline}: 1.6 s).
#' @section Interval-censored data:
#' With \code{\link{Interval}(left, right)} on the left-hand side, every subject contributes
#' \eqn{S(L) - S(R)} (an event between two visits), \eqn{1 - S(R)} (left censoring), \eqn{S(L)}
#' (right censoring) or the density (an exact time) to the likelihood, for every family, with exact
#' first and second derivatives (Lindsey and Ryan 1998; Sun 2006). The fit starts from the fit to the
#' midpoints of the intervals, which is a starting value only. This is the parametric counterpart of
#' \code{\link[survival]{survreg}} with \code{Surv(left, right, type = "interval2")} (and covers
#' also the Gompertz, generalized gamma, Royston-Parmar and piecewise exponential models).
#' \code{print()}, \code{predict()}, \code{AIC()}, \code{estimates()} work as usual; \code{plot()}
#' and \code{\link{hazard_plot}} (survival or cumulative hazard) draw the Turnbull estimate
#' (\code{\link{turnbull}}) as the non-parametric reference. Residuals are not defined for these data.
#' Delayed entry cannot be combined with interval censoring.
#' @return
#' An S7 object of class \code{ParSurvFit}. Main slots: \code{coefficients}, \code{vcov},
#' \code{loglik}, \code{baseline}, \code{scale}, \code{n}, \code{nevent}, \code{converged}, and the
#' knots or cut points of the spline and piecewise baselines.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(fit)}, \code{summary(fit)} \tab HRs or time ratios with intervals and p-values,
#'     shape parameters, \code{Reading:} lines \cr
#'   \code{\link{estimates}(fit)} \tab the table as a data frame \cr
#'   \code{coef}, \code{vcov}, \code{confint}, \code{logLik}, \code{AIC}, \code{BIC}, \code{nobs}
#'     \tab extractors (\code{AIC(fit1, fit2)} compares distributions) \cr
#'   \code{predict(fit, newdata, times, type)} \tab \code{"survival"}, \code{"hazard"},
#'     \code{"cumhaz"}, \code{"quantile"} (at probabilities \code{p}), \code{"rmst"}; \code{se = TRUE}
#'     adds delta-method standard errors and limits (see \code{\link{predict.SurvX}}) \cr
#'   \code{plot(fit)}, \code{\link{hazard_plot}} \tab fitted survival, hazard or cumulative hazard,
#'     with the non-parametric estimate for comparison \cr
#'   \code{residuals(fit, type)} \tab Cox-Snell, martingale and deviance residuals \cr
#' }
#' @references
#' \emph{Parametric models (general):}
#' Kalbfleisch JD, Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd ed.
#'   Wiley, Chapters 2-3.
#' Lawless JF (2003). Statistical Models and Methods for Lifetime Data, 2nd ed. Wiley, Chapters 4-6.
#' Klein JP, Moeschberger ML (2003). Survival Analysis: Techniques for Censored and Truncated Data,
#'   2nd ed. Springer, Sections 3.4 (delayed entry), 3.5 and 12.
#'
#' \emph{Spline (Royston-Parmar) model:}
#' Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards and proportional-odds
#'   models for censored survival data, with application to prognostic modelling and estimation of
#'   treatment effects. Statistics in Medicine 21:2175-2197.
#'
#' \emph{Piecewise exponential model:}
#' Friedman M (1982). Piecewise exponential models for survival data with covariates. Annals of
#'   Statistics 10:101-113.
#'
#' \emph{Generalized gamma:}
#' Prentice RL (1974). A log gamma model and its maximum likelihood estimation. Biometrika
#'   61:539-544.
#' Cox C, Chu H, Schneider MF, Munoz A (2007). Parametric survival analysis and taxonomy of hazard
#'   functions for the generalized gamma distribution. Statistics in Medicine 26:4352-4374.
#'
#' \emph{Model comparison:}
#' Akaike H (1974). A new look at the statistical model identification. IEEE Transactions on
#'   Automatic Control 19:716-723.
#'
#' All references: \code{sx_references("parsurv")}.
#' @seealso
#' \code{\link{Interval}}, \code{\link{turnbull}} for interval-censored data;
#' \code{\link{cox}}; \code{\link{csc}} with a parametric \code{baseline} for competing risks;
#' \code{\link{cure}}; \code{\link{hazard_plot}}, \code{\link{hazard_np}} and \code{\link{dist_plot}}
#' to check the fit; \code{\link{predict.SurvX}}.
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(treat = rbinom(n, 1, 0.5), age = round(rnorm(n, 60, 10)))
#' t <- rweibull(n, shape = 1.5, scale = 4 * exp(0.4 * d$treat - 0.02 * (d$age - 60)))
#' cens <- runif(n, 0, 8)
#' d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
#'
#' ## 1. Weibull PH: hazard ratios
#' fit <- parsurv(Event(time, status) ~ treat + age, d)
#' fit
#'
#' ## 2. the same model as AFT: time ratios
#' parsurv(Event(time, status) ~ treat + age, d, scale = "AFT")
#'
#' ## 3. choosing the distribution
#' AIC(fit,
#'     parsurv(Event(time, status) ~ treat + age, d, baseline = "lognormal"),
#'     parsurv(Event(time, status) ~ treat + age, d, baseline = "rp", df = 3))
#' hazard_plot(fit)                       # fitted against non-parametric hazard
#'
#' ## 4. absolute predictions for two patients aged 60
#' nd <- data.frame(treat = 0:1, age = 60)
#' predict(fit, nd, times = c(2, 5), se = TRUE)$est        # survival
#' predict(fit, nd, type = "quantile", p = 0.5)            # median time
#' predict(fit, nd, times = 5, type = "rmst")              # restricted mean up to 5
#'
#' ## 5. generalized gamma: Q near 1 supports the Weibull
#' summary(parsurv(Event(time, status) ~ treat + age, d, baseline = "gengamma"))
#' @export
parsurv <- function(formula, data, baseline = "weibull", scale = NULL, df = 3L, knots = NULL, cuts = NULL,
                    init = NULL, max_iter = 200L, tol = 1e-10) {
  cl <- match.call()
  dist <- .par_dist(baseline)
  scale <- .par_scale(dist, scale)
  f <- .sx_frame(formula, data, entry_ok = TRUE, what = "parsurv()", interval_ok = TRUE, frailty_ok = TRUE)
  if (!is.null(f$frailty)) return(.par_frailty(cl, f, dist, scale, init, max_iter, tol, df, knots, cuts))
  if (!is.null(f$time2)) {
    if (any(f$status > 1L))
      stop("parsurv() is for one event; for interval-censored competing risks use directcif()", call. = FALSE)
    return(.parsurv_ic(cl, f, dist, scale, init, max_iter, tol, df, knots, cuts))
  }
  if (any(f$status > 1L))
    stop("parsurv() is for one event (status 0/1); for competing risks use csc(..., baseline = \"", dist, "\")", call. = FALSE)
  .parsurv_core(cl, f, f$status, dist, scale, init, max_iter, tol, "the event", df, knots, cuts)
}

.par_fit_ic <- function(L, R, type, X, dist, scale, init, max_iter, tol, knots) {
  n <- length(L); p <- ncol(X)
  P <- 1L + p + .nshape(dist, knots)
  r <- .C("SX_par_fit_ic", as.integer(n), as.integer(p), .par_dists[[dist]], as.integer(scale == "AFT"),
          as.integer(length(knots)), if (length(knots)) as.double(knots) else double(1),
          as.double(L), as.double(R), as.integer(type), .rm(X), as.double(init), as.integer(max_iter), as.double(tol),
          theta = double(P), var = double(P * P), info_d = double(3), info_i = integer(4), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "parametric fit (interval-censored)")
  nm <- c("(Intercept)", colnames(X), .shape_name(dist, scale, knots))
  V <- matrix(r$var, P, P, byrow = TRUE, dimnames = list(nm, nm))
  ok <- all(is.finite(r$theta)) && all(is.finite(diag(V))) && all(diag(V) >= 0) && is.finite(r$info_d[2L])
  list(theta = stats::setNames(r$theta, nm), var = V,
       loglik = r$info_d[2L], iterations = r$info_i[1L], converged = as.logical(r$info_i[2L]) && ok,
       flags = r$info_i[3L], nevent = r$info_i[4L])
}

# interval-censored parametric fit: the likelihood is maximised from the fit to the midpoints of the
# intervals (a start only; the estimates are those of the interval-censored likelihood)
.parsurv_ic <- function(cl, f, dist, scale, init, max_iter, tol, df, knots, cuts) {
  L <- f$time; R <- f$time2
  type <- .ic_type(L, R)
  if (!any(type > 0L)) stop("no events: every subject is right-censored", call. = FALSE)
  mid <- .ic_mid(L, R)
  ev <- as.integer(type > 0L)
  kn <- .par_knots(dist, mid, ev, df, knots, cuts)
  start_of <- function(X, given) {
    if (!is.null(given)) return(as.double(given))
    s0 <- tryCatch(.par_fit(mid, ev, NULL, X, dist, scale, NULL, max_iter, tol, kn), error = function(e) NULL)
    if (is.null(s0) || !all(is.finite(s0$theta)))
      stop("no starting values: the model could not be fitted to the midpoints of the intervals; give init", call. = FALSE)
    unname(s0$theta)
  }
  r <- .par_fit_ic(L, R, type, f$X, dist, scale, start_of(f$X, init), max_iter, tol, kn)
  X0 <- f$X[, 0L, drop = FALSE]
  ln <- if (ncol(f$X)) tryCatch(.par_fit_ic(L, R, type, X0, dist, scale, start_of(X0, NULL), max_iter, tol, kn)$loglik,
                                error = function(e) NA_real_) else r$loglik
  ParSurvFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
             coefficients = r$theta, vcov = r$var, loglik = r$loglik, n = length(L),
             nevent = as.integer(r$nevent), iterations = as.integer(r$iterations), converged = isTRUE(r$converged),
             flags = as.integer(r$flags), time = ifelse(is.finite(R), R, L), status = ev, X = f$X,
             dist = dist, scale = scale, dist_code = .par_dists[[dist]], aft = as.integer(scale == "AFT"),
             entry = NULL, loglik_null = ln, event_label = "the event", knots = kn,
             interval = list(left = L, right = R, type = type))
}

.pcsc <- function(cl, f, baseline, scale, max_iter, tol, df = 3L, knots = NULL, cuts = NULL, kn_list = NULL) {
  K <- max(f$status)
  if (K < 1L) stop("no events", call. = FALSE)
  dist <- rep_len(vapply(baseline, .par_dist, character(1)), K)
  sc <- rep_len(if (is.null(scale)) NA_character_ else scale, K)
  sc <- vapply(seq_len(K), function(k) .par_scale(dist[k], sc[k]), character(1))
  models <- lapply(seq_len(K), function(k)
    tryCatch(.parsurv_core(cl, f, as.integer(f$status == k), dist[k], sc[k], NULL, max_iter, tol, paste("cause", k),
                           df, knots, cuts, if (is.null(kn_list)) NULL else kn_list[[k]]),
             error = function(e) stop("cause ", k, ": ", conditionMessage(e), call. = FALSE)))
  names(models) <- paste0("cause", seq_len(K))
  ParCSCFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts, models = models,
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
  lim <- .sx_log_limits(est, se, zq)
  lo <- up <- est
  lo[] <- lim$lower; up[] <- lim$upper                 # keep the shape and names of est
  list(lower = lo, upper = up)
}

# predictions of a parametric fit at design rows Xn (no intercept column)
.par_predict_X <- function(object, Xn, at, type, se, level, average = FALSE) {
  nn <- nrow(Xn); nt <- length(at)
  usevar <- isTRUE(se) && all(is.finite(object@vcov))
  cells <- if (average) nt else nn * nt
  kn <- object@knots
  r <- .C("SX_par_predict", object@dist_code, object@aft, as.integer(length(kn)), if (length(kn)) as.double(kn) else double(1),
          as.integer(ncol(Xn)), as.double(object@coefficients),
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

# predictions of parametric fits: documented in man/predict.SurvX.Rd (roxygen block in R/cox.R)
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
  kn <- lapply(object@models, function(m) m@knots)
  r <- .C("SX_parcsc_predict", as.integer(K), vapply(object@models, function(m) m@dist_code, integer(1)),
          vapply(object@models, function(m) m@aft, integer(1)), vapply(kn, length, integer(1)),
          if (length(unlist(kn))) as.double(unlist(kn)) else double(1), as.integer(ncol(Xn)), as.double(th),
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
  structure(.sx_reduce("sum", vapply(object@models, function(m) m@loglik, numeric(1))),
            df = sum(vapply(object@models, function(m) length(m@coefficients), integer(1))), nobs = object@n, class = "logLik")
}
S7::method(confint, ParCSCFit) <- function(object, parm, level = 0.95, ...) lapply(object@models, function(m) confint(m, level = level))
S7::method(formula, ParCSCFit) <- function(x, ...) stats::formula(x@terms)
S7::method(getCall, ParCSCFit) <- function(x, ...) x@model_call
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
      if (m@dist == "rp") list("gamma0 (intercept of the spline), all covariates 0", r$estimate, r$lower, r$upper)
      else if (m@dist == "pwexp") list("rate in the first interval exp(intercept), all covariates 0", .sx_ew("exp", r$estimate), .sx_ew("exp", r$lower), .sx_ew("exp", r$upper))
      else if (m@scale == "AFT") list("time scale exp(intercept), all covariates 0", .sx_ew("exp", r$estimate), .sx_ew("exp", r$lower), .sx_ew("exp", r$upper))
      else list("baseline rate exp(intercept), all covariates 0", .sx_ew("exp", r$estimate), .sx_ew("exp", r$lower), .sx_ew("exp", r$upper))
    } else if (m@dist == "rp") list(paste0(r$term, " (spline coefficient)"), r$estimate, r$lower, r$upper)
    else if (m@dist == "pwexp") list(sub("log rate ratio", "rate ratio against the first interval,", r$term), .sx_ew("exp", r$estimate), .sx_ew("exp", r$lower), .sx_ew("exp", r$upper))
    else if (r$term == "Q") list("Q (0: log-normal, 1: Weibull, equal to sigma: gamma)", r$estimate, r$lower, r$upper)
    else if (r$term == "shape") list("Gompertz shape a (a < 0: a fraction never fails)", r$estimate, r$lower, r$upper)
    else if (r$term == "log(shape)") list("Weibull shape (1: constant hazard, > 1 rising, < 1 falling)", .sx_ew("exp", r$estimate), .sx_ew("exp", r$lower), .sx_ew("exp", r$upper))
    else list(paste0("scale sigma of log time", if (m@dist == "weibull") " (Weibull shape = 1/sigma)" else ""), .sx_ew("exp", r$estimate), .sx_ew("exp", r$lower), .sx_ew("exp", r$upper))
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
      size <- if (r >= 2) paste0(rt$ratio, " times as long a") else paste0(.sx_round(.sx_ew("mul", .sx_ew("absdiff", r, 1), 100)), "% ", if (r > 1) "longer" else "shorter")
      paste0(who[i], ": ", size, " time to ", what, " (time ratio ", rt$ratio, ", ", ci, ").")
    } else {
      size <- if (r >= 2) paste0(rt$ratio, " times the") else paste0(.sx_round(.sx_ew("mul", .sx_ew("absdiff", r, 1), 100)), "% ", if (r > 1) "higher" else "lower")
      paste0(who[i], ": ", size, " rate of ", what, " (", ratio_label, " ", rt$ratio, ", ", ci, ").")
    }
    out <- c(out, paste0("  - ", sent))
  }
  out
}

# readings of effects on the cumulative incidence scale (Fine-Gray and direct models)
.reading_cif <- function(e, xlevels, link, what, level) {
  e <- e[e$role == "beta" & is.finite(e$p), , drop = FALSE]
  if (!nrow(e)) return(character(0))
  sig <- e[e$p < 0.05, , drop = FALSE]
  if (!nrow(sig)) return(c("  Reading:", "    No covariate effect reached p < 0.05."))
  lab <- .pretty_terms(sig$term, xlevels)
  is_level <- grepl(": ", lab, fixed = TRUE) & grepl(" vs ", lab, fixed = TRUE)
  who <- paste0(ifelse(is_level, sub(": ", " (", lab, fixed = TRUE), paste0(lab, " (per one-unit increase")), ")")
  out <- "  Reading (effects with p < 0.05):"
  for (i in seq_len(nrow(sig))) {
    r <- sig$exp_estimate[i]
    rt <- .ratio_text(r, sig$exp_lower[i], sig$exp_upper[i])
    ci <- paste0(.pct(level), " CI ", rt$ci)
    dir <- if (r > 1) "higher" else "lower"
    sent <- if (identical(link, "PO")) {
      size <- if (r >= 2) paste0(rt$ratio, " times the") else paste0(.sx_round(.sx_ew("mul", .sx_ew("absdiff", r, 1), 100)), "% ", dir)
      paste0(who[i], ": ", size, " odds of having had ", what, " by any time (OR ", rt$ratio, ", ", ci, ").")
    } else if (identical(link, "PH")) {
      paste0(who[i], ": ", dir, " cumulative incidence of ", what, " at every time (sHR ", rt$ratio, ", ", ci, ").")
    } else {
      paste0(who[i], ": ", dir, " cumulative incidence of ", what, " at every time (exp(b) ", rt$ratio, ", ", ci,
             "; generalized odds rate scale).")
    }
    out <- c(out, paste0("    - ", sent))
  }
  if (identical(link, "PH")) out <- c(out, "    The size of an sHR is not a ratio of risks; predict() and cif_plot() give absolute risks.")
  out
}

.fit_line <- function(m) {
  P <- length(m@coefficients)
  ic <- .sx_ic(m@loglik, P, m@n); aic <- ic$aic; bic <- ic$bic
  out <- sprintf("Fit: log-likelihood %.2f | AIC %.2f | BIC %.2f | %d parameters", m@loglik, aic, bic, P)
  df <- ncol(m@X)
  if (df > 0L && is.finite(m@loglik_null)) {
    t <- .sx_lr(m@loglik, m@loglik_null, df); lr <- t$stat
    pv <- .pval(t$p)
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
           if (.is_ic(x)) .ic_count_line(x@interval$type)
           else paste0("  Subjects: ", x@n, "   events: ", x@nevent, "   censored: ", cens,
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
  out <- c(out, "", sprintf("All causes: log-likelihood %.2f | AIC %.2f", as.numeric(ll), .sx_aic(as.numeric(ll), attr(ll, "df"))),
           "Check: the cumulative incidences of all causes and the event-free probability add up to one (coherent by construction).")
  if (detail)
    out <- c(out, "",
             "csHR: effect on the rate of one cause among subjects who are still event-free.",
             "A covariate can raise the rate of a cause and still lower its cumulative incidence,",
             "because the incidence also depends on the other causes; use predict() or cif_plot() for absolute risks.")
  out
}

S7::method(print, ParCSCFit) <- function(x, level = 0.95, ...) {
  cat(.render_parcsc(x, level, FALSE), "", .method_lines(.method_keys(x)), sep = "\n")
  invisible(x)
}
S7::method(summary, ParCSCFit) <- function(object, level = 0.95, ...) {
  FitSummary(model_call = object@model_call, coefficients = estimates(object, level = level), tables = list(),
             text = c(.render_parcsc(object, level, TRUE), "", .method_lines(.method_keys(object))))
}

# ----------------------------------------------------------------------------------------------
# non-parametric hazard
# ----------------------------------------------------------------------------------------------

# Nelson-Aalen increments of cause k (delayed entry allowed): distinct event times, events, at risk
.na_increments <- function(time, status, k, entry = NULL) {
  n <- length(time)
  r <- .C("SX_event_counts", as.integer(n), as.double(time), as.integer(status), as.integer(k), as.integer(!is.null(entry)),
          if (is.null(entry)) double(1) else as.double(entry), et = double(n), d = double(n), Y = double(n), m = integer(1),
          rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "event counts")
  if (!r$m) return(NULL)
  j <- seq_len(r$m)
  list(time = r$et[j], d = r$d[j], Y = r$Y[j])
}

.na_frame <- function(inc, level) {
  m <- length(inc$time)
  r <- .C("SX_na_table", as.integer(m), as.double(inc$d), as.double(inc$Y), as.double(.zq(level)), H = double(m), se = double(m),
          lower = double(m), upper = double(m), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "Nelson-Aalen estimate")
  data.frame(time = inc$time, n.risk = inc$Y, n.event = inc$d, cumhaz = r$H, se = r$se, lower = r$lower, upper = r$upper)
}

# the smoothed hazard is shown up to the last event time with at least max(10, 10% of the initial risk
# set) subjects at risk (the end of the cross-validation region in the C++ core); beyond it the
# estimate rests on a handful of subjects
.kernel_end <- function(inc) {
  need <- max(10, .sx_ew("mul", 0.1, inc$Y[1L]))
  ok <- inc$Y >= need
  if (any(ok)) max(inc$time[ok]) else inc$time[1L]
}

.kernel_frame <- function(inc, bw, boundary, grid, level) {
  x <- .sx_seq(0, .kernel_end(inc), grid)
  r <- .C("SX_kernel_hazard", as.integer(length(inc$time)), as.double(inc$time), as.double(inc$d), as.double(inc$Y),
          as.double(if (is.null(bw)) -1 else bw), as.integer(boundary), as.integer(grid), as.double(x),
          h = double(grid), se = double(grid), bw = double(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "hazard_np")
  h <- r$h; se <- r$se
  lim <- .sx_log_limits(h, se, .zq(level))
  list(frame = data.frame(time = x, hazard = h, se = se, lower = lim$lower, upper = lim$upper), bw = r$bw)
}

#' Non-parametric hazard: Nelson-Aalen and kernel-smoothed hazard rate
#'
#' @description
#' Estimates, without any model, the cumulative hazard (Nelson-Aalen) and the hazard rate itself
#' (kernel smoothing), for every group and, with competing risks, for every cause (cause-specific
#' hazards). Delayed entry is allowed.
#'
#' \strong{Use it when} you want to see the \emph{shape} of the hazard over time (rising, falling,
#' peaking), to choose a parametric distribution for \code{\link{parsurv}} or \code{\link{csc}}, to
#' check a fitted parametric model (\code{\link{hazard_plot}} overlays both), or when there is
#' delayed entry and \code{\link{km}} cannot be used.
#' @param formula \code{Event(time, status) ~ 1} or \code{~ group};
#'   \code{Event(time, status, entry = entry)} for delayed entry. With status 1, 2, ... one set of
#'   curves per cause.
#' @param data a data frame.
#' @param smooth \code{TRUE} (default): also compute the kernel-smoothed hazard rate.
#' @param bw bandwidth of the smoother, in time units. \code{NULL} (default) chooses it by
#'   cross-validation, separately for every curve. A larger bandwidth gives a smoother but more biased
#'   curve; try a few values if the automatic choice looks too wiggly or too flat.
#' @param boundary \code{TRUE} (default): boundary kernels near the ends of follow-up, so that the
#'   estimate is not biased downwards at time zero.
#' @param level confidence level of the pointwise intervals.
#' @param grid number of time points at which the smoothed hazard is computed.
#' @details
#' The cumulative hazard is \eqn{\hat H(t) = \sum_{t_j \le t} d_j / Y_j} with variance
#' \eqn{\sum d_j / Y_j^2}; limits on the log scale. The smoothed hazard is
#' \eqn{\hat h(t) = b^{-1}\sum_j K((t - t_j)/b)\, d_j / Y_j} with the Epanechnikov kernel and, within
#' one bandwidth of either end of follow-up, the boundary kernels of Mueller and Wang (1994), so the
#' estimate is not biased downwards at time zero. Its variance is \eqn{b^{-2}\sum_j K^2 d_j / Y_j^2};
#' limits on the log scale. The smoothed hazard is shown up to the last event time with at least
#' max(10, 10\% of the initial risk set) subjects at risk; beyond it the estimate rests on a handful of
#' subjects. The bandwidth \eqn{b} is chosen on that range from 40 values (1/50 to 1/2 of it) as the
#' largest local minimiser of the least-squares cross-validation criterion of Patil (1993), which is
#' far more stable than its global minimiser (Hall and Marron 1991), unless \code{bw} is given.
#' @section Interpreting the result:
#' The hazard \eqn{h(t)} is the instantaneous event rate among subjects still at risk, in events per
#' unit of time (0.05 per year: about 5 events per 100 person-years at that time). Its shape guides
#' the choice of distribution: constant (exponential), monotone (Weibull, Gompertz), rising then
#' falling (log-normal, log-logistic), or more complex (\code{"rp"} spline). The cumulative hazard is
#' easier to estimate and its slope is the hazard. The intervals of the smoothed hazard are wide in
#' the tails; do not over-interpret wiggles there.
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}), the Nelson-Aalen cumulative hazard agrees with
#' \code{survival::survfit(..., ctype = 1)} to rounding error (\eqn{10^{-15}}) on 81 data sets.
#' @return
#' An S7 object of class \code{HazardFit}: \code{fit@cumhaz} and \code{fit@hazard} are lists of data
#' frames (one per group and cause) with the estimates and limits; \code{fit@bandwidth} holds the
#' bandwidths. Methods: \code{print}, and \code{plot(fit, type = c("hazard", "cumhaz"))}.
#' @references
#' Nelson W (1972) Theory and applications of hazard plotting for censored failure data.
#'   Technometrics 14:945-966.
#'
#' Aalen O (1978) Nonparametric inference for a family of counting processes. Annals of Statistics
#'   6:701-726.
#'
#' Ramlau-Hansen H (1983) Smoothing counting process intensities by means of kernel functions.
#'   Annals of Statistics 11:453-466.
#'
#' Mueller HG, Wang JL (1994) Hazard rate estimation under random censoring with varying kernels and
#'   bandwidths. Biometrics 50:61-76.
#'
#' Patil PN (1993) On the least squares cross-validation bandwidth in hazard rate estimation. Annals
#'   of Statistics 21:1792-1810.
#'
#' Hall P, Marron JS (1991) Local minima in cross-validation functions. Journal of the Royal
#'   Statistical Society, Series B 53:245-252.
#'
#' All references: \code{sx_references("hazard_np")}.
#' @seealso
#' \code{\link{hazard_plot}} (fitted against non-parametric), \code{\link{dist_plot}}, \code{\link{km}},
#' \code{\link{parsurv}}.
#' @examples
#' set.seed(3)
#' d <- data.frame(time = rweibull(300, 1.5, 3), status = 1, arm = rep(c("A", "B"), 150))
#' hz <- hazard_np(Event(time, status) ~ arm, d)
#' hz                                # bandwidths chosen
#' plot(hz)                          # a rising hazard: Weibull with shape > 1
#' plot(hz, type = "cumhaz")
#'
#' ## a fixed, wider bandwidth for a smoother curve
#' plot(hazard_np(Event(time, status) ~ 1, d, bw = 1.5))
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
  HazardFit(model_call = cl, cumhaz = cum, hazard = haz, bandwidth = bws, K = as.integer(K), level = level,
            time = f$time, status = f$status, group = f$group)
}

S7::method(print, HazardFit) <- function(x, ...) {
  cat(paste0("Non-parametric hazard", if (x@K > 1L) " (cause-specific)" else ""), "\n", sep = "")
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
  cat("\n", .method_lines(.method_keys(x)), sep = "\n")
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
#' @description
#' Draws the fitted hazard (or cumulative hazard, or survival) of one or several parametric models at
#' typical or chosen covariate values, with delta-method confidence bands, and overlays the
#' non-parametric estimate of the same data (\code{\link{hazard_np}}, Nelson-Aalen or
#' Kaplan-Meier) as a reference. Agreement of the fitted and the non-parametric curve is the visual
#' check of the chosen distribution.
#' @param x a model fitted by \code{\link{parsurv}} or by \code{\link{csc}} with a parametric
#'   baseline, or a named list of such models fitted to the same data.
#' @param newdata data frame of covariate values, one curve per row; default typical values (mean of
#'   numeric covariates, most frequent level of categorical covariates).
#' @param type \code{"hazard"}, \code{"cumhaz"} or \code{"survival"}.
#' @param cause causes to show (cause-specific models); default all.
#' @param reference overlay the non-parametric estimate of all subjects. It is a fair check of the
#'   fitted curve at typical covariate values or without covariates.
#' @param conf.int draw pointwise confidence bands.
#' @param level confidence level.
#' @param times time grid; default 201 points up to the largest observed time.
#' @param col,lwd colours and line width of the models.
#' @param xlab,ylab,main axis labels and title.
#' @param legend draw a legend.
#' @param cex character expansion.
#' @param ... passed to \code{plot.default}.
#' @section Interpreting the plot:
#' The fitted curve should run inside the confidence band of the non-parametric estimate over the
#' range where many subjects are at risk. Systematic departures (a fitted monotone hazard where the
#' non-parametric one rises and falls, for example) suggest another distribution: try a log-normal
#' or log-logistic, or the flexible \code{"rp"} spline, and compare them with \code{AIC()}. Departures
#' in the far tail, where few subjects remain, matter less.
#' @return
#' invisibly, a data frame of the plotted curves.
#' @references
#' Klein JP, Moeschberger ML (2003) Survival Analysis: Techniques for Censored and Truncated Data,
#'   2nd ed. Springer, Chapter 12.
#'
#' Nelson W (1972) Theory and applications of hazard plotting for censored failure data.
#'   Technometrics 14:945-966.
#'
#' Ramlau-Hansen H (1983) Smoothing counting process intensities by means of kernel functions. Annals
#'   of Statistics 11:453-466.
#'
#' Mueller HG, Wang JL (1994) Hazard rate estimation under random censoring with varying kernels and
#'   bandwidths. Biometrics 50:61-76.
#' @seealso
#' \code{\link{hazard_np}}, \code{\link{dist_plot}}, \code{\link{parsurv}}, \code{\link{csc}},
#' \code{\link{cure}}.
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
  ic <- .is_ic(ref)
  K <- if (S7::S7_inherits(ref, ParCSCFit)) ref@K else 1L
  causes <- if (type == "survival") 0L else if (is.null(cause)) seq_len(K) else as.integer(cause)
  tg <- if (is.null(times)) .sx_seq(0, max(ref@time), 201L) else sort(times)
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
  if (reference && ic && type != "hazard") {
    # interval-censored data: the Turnbull estimate (survival, or minus its log)
    tb <- .tb_fit(ref@interval$left, ref@interval$right)
    sc <- .tb_steps(tb)
    v <- if (type == "survival") sc$surv else .sx_ew("neglog_floor", sc$surv)
    ok <- is.finite(v) & (type == "survival" | sc$surv > 0)
    refd <- data.frame(cause = 0L + (type != "survival"), time = sc$time[ok], value = v[ok], step = FALSE)
  } else if (reference && !ic) {
    st <- if (S7::S7_inherits(ref, ParCSCFit)) ref@status else ref@status
    refd <- do.call(rbind, lapply(causes, function(k) {
      kk <- if (k == 0L) NULL else k
      if (type == "survival") {
        inc <- .na_increments(ref@time, as.integer(st > 0L), 1L, ref@entry)
        if (is.null(inc)) return(NULL)
        s <- .sx_np_curve_inc(inc$d, inc$Y, "km")
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
    ylim <- if (type == "survival") c(0, 1) else c(0, if (length(yv)) .sx_quantile(yv, 0.995) * 1.05 else 1)
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
      lg <- c(lab, if (!is.null(rk) && nrow(rk)) if (ic) "Turnbull NPMLE (all subjects)" else "non-parametric (all subjects)")
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
