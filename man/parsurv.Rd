\name{parsurv}
\alias{parsurv}
\title{Parametric survival regression}
\description{
Maximum likelihood regression for a single event time with an exponential, Weibull, Gompertz,
log-normal, log-logistic, piecewise exponential or flexible spline (Royston-Parmar) distribution,
under right censoring and, optionally, delayed entry (left truncation). Covariates act on the hazard (proportional hazards, \code{scale = "PH"}) or on
the time scale (accelerated failure time, \code{scale = "AFT"}).
}
\usage{
parsurv(formula, data, baseline = "weibull", scale = NULL, df = 3L, knots = NULL, cuts = NULL,
        init = NULL, max_iter = 200L, tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}, status 0/1; \code{Event(time, status,
  entry = )} for delayed entry.}
\item{data}{a data frame.}
\item{baseline}{\code{"weibull"} (default), \code{"exponential"}, \code{"gompertz"},
  \code{"lognormal"}, \code{"loglogistic"}, \code{"rp"} (Royston-Parmar spline) or
  \code{"pwexp"} (piecewise exponential).}
\item{scale}{\code{"PH"} or \code{"AFT"}; default \code{"PH"} where it exists (exponential,
  Weibull, Gompertz, spline, piecewise exponential), else \code{"AFT"} (log-normal, log-logistic).}
\item{df}{spline baseline: degrees of freedom of the spline, the number of internal knots plus one;
  internal knots at equally spaced quantiles of the log event times.}
\item{knots}{spline baseline: internal knots on the time scale, instead of \code{df}.}
\item{cuts}{piecewise exponential baseline: the cut points of the time axis; default the quartiles of
  the event times.}
\item{init}{optional starting values (intercept, covariates, shape parameter).}
\item{max_iter,tol}{iteration limit and convergence tolerance (relative change of the
  log-likelihood).}
}
\details{
Parametrisations, with \eqn{\eta = \beta_0 + z'\beta}:
\describe{
  \item{PH}{\eqn{H(t \mid z) = e^{\eta} H_0(t)} with \eqn{H_0(t) = t} (exponential),
    \eqn{t^{\kappa}} (Weibull; coefficient \code{log(shape)} \eqn{= \log\kappa}) or
    \eqn{(e^{a t} - 1)/a} (Gompertz; coefficient \code{shape} \eqn{= a}, any sign; \eqn{a < 0}
    gives a fraction that never fails). \eqn{e^{\beta_j}} is a hazard ratio. These are the
    \code{"exp"}, \code{"weibullPH"} and \code{"gompertz"} models of \pkg{flexsurv}.}
  \item{AFT}{\eqn{\log T = \eta + \sigma W}, coefficient \code{log(scale)} \eqn{= \log\sigma},
    with \eqn{W} standard minimum extreme value (Weibull; exponential with \eqn{\sigma = 1}),
    standard normal (log-normal) or standard logistic (log-logistic). \eqn{e^{\beta_j}} is a time
    ratio (acceleration factor). These are the models of \code{survival::survreg}.}
  \item{spline (\code{"rp"})}{\eqn{\log H(t \mid z) = s(\log t) + z'\beta} with a restricted
    cubic spline \eqn{s(x) = \gamma_0 + \gamma_1 x + \sum_j \gamma_{j+1} v_j(x)}, linear beyond
    the boundary knots at the smallest and largest log event time (Royston and Parmar 2002); the
    basis of \code{flexsurv::flexsurvspline}. \code{(Intercept)} is \eqn{\gamma_0}; the
    coefficients \code{gamma1}, \code{gamma2}, ... are the shape parameters; \eqn{e^{\beta_j}} is a
    hazard ratio. With \code{df = 1} it is the Weibull model.}
  \item{piecewise exponential (\code{"pwexp"})}{hazard \eqn{e^{\eta + a_j}} on the interval
    \eqn{(c_{j-1}, c_j]}, \eqn{a_1 = 0} (Friedman 1982); the coefficients \code{log rate ratio}
    compare each later interval with the first. The maximum likelihood estimates equal those of a
    Poisson regression of the events on the split follow-up time.}
}
The log-likelihood is \eqn{\sum_i d_i \log h(t_i) + \log S(t_i) - \log S(e_i)} with entry time
\eqn{e_i} (zero without delayed entry). Gradient and Hessian are analytic and the maximum is found
by Newton-Raphson with step halving (and a Levenberg-Marquardt ridge where the information is not
positive definite), starting from the exponential fit. Standard errors are from the observed
information.

\code{print()} shows the effects on the ratio scale with a one-line reading of every effect with
\eqn{p < 0.05}; \code{summary()} adds the baseline distribution on its natural scale (rate or time
scale, shape) and notes on interpretation.
}
\value{
an S7 object of class \code{ParSurvFit} with the methods \code{print}, \code{summary},
  \code{estimates}, \code{coef}, \code{vcov}, \code{confint}, \code{logLik}, \code{AIC},
  \code{BIC}, \code{update} and \code{predict(fit, newdata, times, type)} with \code{type}
  \code{"survival"}, \code{"hazard"}, \code{"cumhaz"}, \code{"quantile"} (at probabilities
  \code{p}) or \code{"rmst"} (restricted mean survival time up to \code{times}); see
  \code{\link{predict.SurvX}}.
}
\references{
Kalbfleisch JD, Prentice RL (2002) The Statistical Analysis of Failure Time Data, 2nd ed. Wiley,
  Chapters 2-3.

Lawless JF (2003) Statistical Models and Methods for Lifetime Data, 2nd ed. Wiley, Chapters 4-6.

Klein JP, Moeschberger ML (2003) Survival Analysis: Techniques for Censored and Truncated Data,
  2nd ed. Springer, Sections 3.5 and 12.

Royston P, Parmar MKB (2002) Flexible parametric proportional-hazards and proportional-odds models
  for censored survival data, with application to prognostic modelling and estimation of treatment
  effects. Statistics in Medicine 21:2175-2197.

Friedman M (1982) Piecewise exponential models for survival data with covariates. Annals of
  Statistics 10:101-113.
}
\seealso{
\code{\link{csc}} with a parametric \code{baseline} for competing risks; \code{\link{hazard_plot}}
  and \code{\link{hazard_np}} to check the fit.
}
\examples{
set.seed(1)
n <- 300
x <- rbinom(n, 1, 0.5)
t <- rweibull(n, shape = 1.5, scale = 4 * exp(-0.5 * x)); cens <- runif(n, 0, 8)
d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x)
fit <- parsurv(Event(time, status) ~ x, d)
fit
summary(fit)
predict(fit, data.frame(x = 0:1), times = c(1, 2, 4), se = TRUE)$est
predict(fit, data.frame(x = 0:1), type = "quantile", p = 0.5)
AIC(fit, parsurv(Event(time, status) ~ x, d, baseline = "lognormal"))
fs <- parsurv(Event(time, status) ~ x, d, baseline = "rp", df = 3)
fs
predict(fs, data.frame(x = 0:1), type = "quantile", p = 0.5)
parsurv(Event(time, status) ~ x, d, baseline = "pwexp", cuts = c(2, 4))
}
