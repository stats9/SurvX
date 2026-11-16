\name{parsurv}
\alias{parsurv}
\title{Parametric survival regression}
\description{
Maximum likelihood regression for a single event time with an exponential, Weibull, Gompertz,
log-normal or log-logistic distribution, under right censoring and, optionally, delayed entry
(left truncation). Covariates act on the hazard (proportional hazards, \code{scale = "PH"}) or on
the time scale (accelerated failure time, \code{scale = "AFT"}).
}
\usage{
parsurv(formula, data, baseline = "weibull", scale = NULL, init = NULL, max_iter = 200L,
        tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}, status 0/1; \code{Event(time, status,
  entry = )} for delayed entry.}
\item{data}{a data frame.}
\item{baseline}{\code{"weibull"} (default), \code{"exponential"}, \code{"gompertz"},
  \code{"lognormal"} or \code{"loglogistic"}.}
\item{scale}{\code{"PH"} or \code{"AFT"}; default \code{"PH"} where it exists (exponential,
  Weibull, Gompertz), else \code{"AFT"} (log-normal, log-logistic).}
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
}
