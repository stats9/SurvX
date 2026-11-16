\name{csc}
\alias{csc}
\title{Cause-specific hazard models: Cox or parametric}
\description{
One hazard model per cause, with the other causes treated as censoring, combined into cumulative
incidence functions. With \code{baseline = "cox"} (default) every cause gets a Cox model with an
unspecified baseline hazard; with a parametric \code{baseline} every cause gets a parametric hazard
model fitted as in \code{\link{parsurv}}, and the cumulative incidence
\eqn{F_k(t) = \int_0^t h_k(u) S(u)\,du} is computed by numerical integration (20-point
Gauss-Legendre on a fine partition of the time axis) with delta-method standard errors. Both are
coherent by construction: the incidences of all causes and the event-free probability add up to one.
}
\usage{
csc(formula, data, ties = c("efron", "breslow"), product_limit = TRUE, baseline = "cox",
    scale = NULL, df = 3L, knots = NULL, cuts = NULL, max_iter = 200L, tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}; status 0 = censored, 1..K = cause. With a
  parametric baseline, \code{Event(time, status, entry = )} allows delayed entry.}
\item{data}{a data frame.}
\item{ties}{\code{"efron"} (default) or \code{"breslow"} (Cox baseline only).}
\item{product_limit}{\code{TRUE} (default): overall survival as a product-limit, as in
  \pkg{riskRegression}; \code{FALSE}: \eqn{\exp\{-\sum_k \Lambda_k(t)\}} (Cox baseline only).}
\item{baseline}{\code{"cox"}, or one of the distributions of \code{\link{parsurv}}
  (\code{"exponential"}, \code{"weibull"}, \code{"gompertz"}, \code{"lognormal"},
  \code{"loglogistic"}, \code{"gengamma"}, \code{"rp"}, \code{"pwexp"}); a vector gives one
  distribution per cause.}
\item{scale}{\code{"PH"} or \code{"AFT"} for parametric baselines, one value or one per cause;
  default as in \code{\link{parsurv}}.}
\item{df,knots,cuts}{spline (\code{"rp"}) and piecewise exponential (\code{"pwexp"}) baselines, as
  in \code{\link{parsurv}}; default knots and cut points are chosen from the events of each cause.}
\item{max_iter,tol}{iteration limit and convergence tolerance of the parametric fits.}
}
\value{
With \code{baseline = "cox"}: an S7 object of class \code{CSCFit}; \code{x@models} holds the fit
  of each cause. With a parametric baseline: an S7 object of class \code{ParCSCFit};
  \code{x@models} holds the \code{\link{parsurv}} fit of each cause, and \code{logLik()} and
  \code{AIC()} refer to all causes together. \code{predict(fit, newdata, times, cause)} gives the
  cumulative incidence; for a parametric fit \code{type = "survival"} gives the event-free
  probability and \code{se = TRUE} delta-method standard errors and limits. Parametric fits work
  with \code{\link{cif_plot}}, \code{\link{cif_data}}, \code{\link{hazard_plot}} and
  \code{\link{cv_compare}}.
}
\references{
Prentice RL, Kalbfleisch JD, Peterson AV, Flournoy N, Farewell VT, Breslow NE (1978) The analysis
  of failure times in the presence of competing risks. Biometrics 34:541-554.

Benichou J, Gail MH (1990) Estimates of absolute cause-specific risk in cohort studies.
  Biometrics 46:813-826.

Ozenne B, Sorensen AL, Scheike T, Torp-Pedersen C, Gerds TA (2017) riskRegression: predicting the
  risk of an event using Cox regression models. R Journal 9(2):440-460.

Hinchliffe SR, Lambert PC (2013) Flexible parametric modelling of cause specific hazards to
  estimate cumulative incidence functions. BMC Medical Research Methodology 13:13.
}
\examples{
set.seed(2)
n <- 300
x <- rbinom(n, 1, 0.5)
t1 <- rweibull(n, 1.3, 3 * exp(-0.4 * x)); t2 <- rweibull(n, 0.9, 5); cens <- runif(n, 0, 6)
d <- data.frame(time = pmin(t1, t2, cens), x = x)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
csc(Event(time, status) ~ x, d)
fw <- csc(Event(time, status) ~ x, d, baseline = "weibull")
fw
predict(fw, data.frame(x = 0:1), times = c(1, 3), cause = 1)
AIC(fw, csc(Event(time, status) ~ x, d, baseline = "lognormal"))
}
