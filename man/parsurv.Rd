\name{parsurv}
\alias{parsurv}
\title{Parametric survival regression}
\description{
Fits a fully parametric model for the time to a single event: exponential, Weibull, Gompertz,
log-normal, log-logistic, generalized gamma, piecewise exponential, or a flexible spline model
(Royston-Parmar), under right censoring and, optionally, delayed entry. Covariates act on the
hazard (proportional hazards, \code{scale = "PH"}) or on the time scale (accelerated failure
time, \code{scale = "AFT"}).

\strong{Use it when} the shape of the hazard over time matters, smooth predicted curves are
wanted, survival must be extrapolated beyond follow-up (health economics), or time ratios
(how much longer or shorter the time to the event is) are the effect of interest.
\strong{Use something else when}
\itemize{
  \item only the covariate effects matter and no assumption on the baseline is wanted:
    \code{\link{cox}};
  \item there are competing events: \code{\link{csc}} with the same \code{baseline} (one
    parametric model per cause), or \code{\link{directcif}};
  \item part of the population will never have the event: \code{\link{cure}}.
}
}
\usage{
parsurv(formula, data, baseline = "weibull", scale = NULL, df = 3L, knots = NULL,
        cuts = NULL, init = NULL, max_iter = 200L, tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates} with status 0/1;
  \code{Event(time, status, entry = entry)} for delayed entry (left truncation). Covariates as
  in \code{\link[stats]{lm}}.}
\item{data}{a data frame. Rows with missing values in the variables of the formula are removed.}
\item{baseline}{the distribution:
  \tabular{lll}{
    \code{"weibull"} (default) \tab PH or AFT \tab monotone hazard (rising, falling or constant) \cr
    \code{"exponential"} \tab PH or AFT \tab constant hazard \cr
    \code{"gompertz"} \tab PH \tab hazard rising or falling exponentially (ageing, mortality) \cr
    \code{"lognormal"} \tab AFT \tab hazard rising then falling \cr
    \code{"loglogistic"} \tab AFT \tab hazard rising then falling, closed forms \cr
    \code{"gengamma"} \tab AFT \tab generalized gamma: contains Weibull, log-normal and gamma \cr
    \code{"rp"} \tab PH \tab Royston-Parmar spline: any smooth shape, set by \code{df} \cr
    \code{"pwexp"} \tab PH \tab piecewise constant hazard, cut at \code{cuts} \cr
  }
  Compare candidates with \code{AIC()} and with \code{\link{hazard_plot}} against the
  non-parametric hazard.}
\item{scale}{\code{"PH"} or \code{"AFT"}. Default \code{"PH"} where it exists (exponential,
  Weibull, Gompertz, spline, piecewise exponential), else \code{"AFT"}. For the Weibull and
  exponential the two are the same model with different coefficients (\eqn{\beta_{PH} =
  -\beta_{AFT}/\sigma}).}
\item{df}{spline baseline: degrees of freedom, the number of internal knots plus one (default 3;
  \code{df = 1} is the Weibull model). Internal knots are placed at equally spaced quantiles of the
  log event times. Values from 2 to 5 are usually enough; compare them with \code{AIC()}.}
\item{knots}{spline baseline: internal knots on the time scale, instead of \code{df}. They must
  lie strictly between the smallest and the largest event time.}
\item{cuts}{piecewise exponential baseline: cut points of the time axis (default the quartiles
  of the event times).}
\item{init}{optional starting values (intercept, covariates, shape parameters), for difficult
  fits.}
\item{max_iter}{maximum number of Newton-Raphson iterations (default 200).}
\item{tol}{convergence tolerance on the relative change of the log-likelihood.}
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
  \item{generalized gamma (\code{"gengamma"})}{AFT, \eqn{\log T = \eta + \sigma W} where, for
    \eqn{Q \ne 0}, \eqn{Q^{-2} e^{Q W}} has a gamma distribution with shape \eqn{Q^{-2}} (the
    parametrisation of Prentice 1974 and of \pkg{flexsurv}); coefficients \code{log(scale)}
    \eqn{= \log\sigma} and \code{Q}. It contains the log-normal (\eqn{Q = 0}), Weibull
    (\eqn{Q = 1}) and gamma (\eqn{Q = \sigma}) distributions, so a Wald or likelihood ratio test of
    \eqn{Q} chooses between them (Cox et al. 2007). The log density is evaluated in a form that is
    free of cancellation near \eqn{Q = 0}; the survivor function uses the incomplete gamma function
    for \eqn{|Q| \ge 0.01} and numerical integration of the density for smaller \eqn{|Q|}, so the
    likelihood is smooth across the log-normal limit. Derivatives in \eqn{(\mu, \log\sigma, Q)}
    are central differences of the exact contribution of each subject.}
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

Estimation is by maximum likelihood with analytic first and second derivatives (central
differences for the generalized gamma), by Newton-Raphson with step halving. Standard errors are
the inverse of the observed information. With delayed entry, each subject contributes
\eqn{\log f(t)} or \eqn{\log S(t)} minus \eqn{\log S(entry)} (Klein and Moeschberger 2003,
Section 3.4).
}
\section{Interpreting the coefficients}{
\describe{
  \item{PH models}{\eqn{e^\beta} is a \strong{hazard ratio}, read as in \code{\link{cox}}: HR 0.70
    means a 30\% lower event rate at every time.}
  \item{AFT models}{\eqn{e^\beta} is a \strong{time ratio} (acceleration factor): TR 1.40 means
    that the time to the event is stretched by 40\% (every quantile of survival time, the median
    included, is 1.4 times longer); TR < 1 means shorter times and worse prognosis. Note the
    opposite direction to the HR: a protective factor has HR < 1 but TR > 1.}
  \item{Shape parameters}{(\code{log(shape)}, \code{log(scale)}, \code{Q}, \code{gamma1}, ...)
    describe the baseline and are rarely interpreted on their own. For the Weibull, a shape
    \eqn{\kappa > 1} means a rising hazard; for the generalized gamma, \eqn{Q} near 1 means that
    the Weibull fits, near 0 the log-normal.}
}
Absolute quantities usually speak more clearly than coefficients:
\code{predict(fit, newdata, times)} (survival), \code{type = "quantile", p = 0.5} (median time),
\code{type = "rmst"} (restricted mean). A sentence for a report: \emph{"In a Weibull
accelerated failure time model, treatment prolonged the time to progression by 40\% (time ratio
1.40, 95\% CI 1.15 to 1.70); predicted median 14.2 versus 10.1 months."}
}
\section{Assumptions and how to check them}{
\describe{
  \item{The distribution fits}{compare the fitted hazard and cumulative hazard with the
    non-parametric estimates: \code{hazard_plot(fit)}; \code{\link{dist_plot}} draws the
    linearising plots of each distribution (for example \eqn{\log H} against \eqn{\log t}, a straight
    line for the Weibull); compare candidate distributions with \code{AIC()}.}
  \item{PH or AFT structure}{for PH models, \code{\link{ph_test}} on a \code{cox} fit with the
    same covariates checks proportional hazards.}
  \item{Extrapolation}{beyond follow-up, the curve depends on the distribution only; compare the
    extrapolations of several well-fitting distributions before relying on one.}
}
}
\section{Messages and common problems}{
\describe{
  \item{\code{WARNING: the fit did NOT converge}}{often a distribution that cannot describe the
    data (the generalized gamma on small samples), separation, or a covariate on a very large or
    small scale (rescale it). Try another baseline, rescale covariates, or give \code{init}.}
  \item{\code{the spline model needs at least two distinct event times}}{too few events for
    \code{"rp"}.}
  \item{\code{internal knots must lie strictly between the smallest and the largest event time}}{
    choose \code{knots} inside the range of the event times, or a smaller \code{df}.}
  \item{\code{parsurv() is for one event (status 0/1)}}{use \code{csc(..., baseline = )} for
    competing risks.}
}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}; 82 data sets): AFT models agree with \code{survival::survreg} to \eqn{3 \times
10^{-9}} for coefficients, \eqn{2 \times 10^{-10}} for standard errors and \eqn{3 \times 10^{-11}} for
the log-likelihood; all distributions agree with \code{flexsurv::flexsurvreg} and
\code{flexsurvspline} within the tolerance of their optimisers (log-likelihood to a relative
\eqn{10^{-7}}, covariate coefficients and standard errors to about \eqn{10^{-4}}), and where they
differ more, SurvX reaches the higher log-likelihood. In a Monte Carlo study (1000 replications,
\eqn{n} = 200 to 1000, with and without delayed entry) the 95\% intervals covered 92.9 to 96.3\%.
On \eqn{n = 10000}, a Weibull fit takes about 0.02 s (\code{flexsurvreg}: 0.17 s), a
generalized gamma fit 0.36 s (2.8 s) and a spline fit 0.03 s (\code{flexsurvspline}: 1.6 s).
}
\value{
An S7 object of class \code{ParSurvFit}. Main slots: \code{coefficients}, \code{vcov},
\code{loglik}, \code{baseline}, \code{scale}, \code{n}, \code{nevent}, \code{converged}, and the
knots or cut points of the spline and piecewise baselines.

Methods:
\tabular{ll}{
  \code{print(fit)}, \code{summary(fit)} \tab HRs or time ratios with intervals and p-values,
    shape parameters, \code{Reading:} lines \cr
  \code{\link{estimates}(fit)} \tab the table as a data frame \cr
  \code{coef}, \code{vcov}, \code{confint}, \code{logLik}, \code{AIC}, \code{BIC}, \code{nobs}
    \tab extractors (\code{AIC(fit1, fit2)} compares distributions) \cr
  \code{predict(fit, newdata, times, type)} \tab \code{"survival"}, \code{"hazard"},
    \code{"cumhaz"}, \code{"quantile"} (at probabilities \code{p}), \code{"rmst"}; \code{se = TRUE}
    adds delta-method standard errors and limits (see \code{\link{predict.SurvX}}) \cr
  \code{plot(fit)}, \code{\link{hazard_plot}} \tab fitted survival, hazard or cumulative hazard,
    with the non-parametric estimate for comparison \cr
  \code{residuals(fit, type)} \tab Cox-Snell, martingale and deviance residuals \cr
}
}
\references{
\emph{Parametric models (general):}
Kalbfleisch JD, Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd ed.
  Wiley, Chapters 2-3.
Lawless JF (2003). Statistical Models and Methods for Lifetime Data, 2nd ed. Wiley, Chapters 4-6.
Klein JP, Moeschberger ML (2003). Survival Analysis: Techniques for Censored and Truncated Data,
  2nd ed. Springer, Sections 3.4 (delayed entry), 3.5 and 12.

\emph{Spline (Royston-Parmar) model:}
Royston P, Parmar MKB (2002). Flexible parametric proportional-hazards and proportional-odds
  models for censored survival data, with application to prognostic modelling and estimation of
  treatment effects. Statistics in Medicine 21:2175-2197.

\emph{Piecewise exponential model:}
Friedman M (1982). Piecewise exponential models for survival data with covariates. Annals of
  Statistics 10:101-113.

\emph{Generalized gamma:}
Prentice RL (1974). A log gamma model and its maximum likelihood estimation. Biometrika
  61:539-544.
Cox C, Chu H, Schneider MF, Munoz A (2007). Parametric survival analysis and taxonomy of hazard
  functions for the generalized gamma distribution. Statistics in Medicine 26:4352-4374.

\emph{Model comparison:}
Akaike H (1974). A new look at the statistical model identification. IEEE Transactions on
  Automatic Control 19:716-723.

All references: \code{sx_references("parsurv")}.
}
\seealso{
\code{\link{cox}}; \code{\link{csc}} with a parametric \code{baseline} for competing risks;
\code{\link{cure}}; \code{\link{hazard_plot}}, \code{\link{hazard_np}} and \code{\link{dist_plot}}
to check the fit; \code{\link{predict.SurvX}}.
}
\examples{
set.seed(1)
n <- 300
d <- data.frame(treat = rbinom(n, 1, 0.5), age = round(rnorm(n, 60, 10)))
t <- rweibull(n, shape = 1.5, scale = 4 * exp(0.4 * d$treat - 0.02 * (d$age - 60)))
cens <- runif(n, 0, 8)
d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)

## 1. Weibull PH: hazard ratios
fit <- parsurv(Event(time, status) ~ treat + age, d)
fit

## 2. the same model as AFT: time ratios
parsurv(Event(time, status) ~ treat + age, d, scale = "AFT")

## 3. choosing the distribution
AIC(fit,
    parsurv(Event(time, status) ~ treat + age, d, baseline = "lognormal"),
    parsurv(Event(time, status) ~ treat + age, d, baseline = "rp", df = 3))
hazard_plot(fit)                       # fitted against non-parametric hazard

## 4. absolute predictions for two patients aged 60
nd <- data.frame(treat = 0:1, age = 60)
predict(fit, nd, times = c(2, 5), se = TRUE)$est        # survival
predict(fit, nd, type = "quantile", p = 0.5)            # median time
predict(fit, nd, times = 5, type = "rmst")              # restricted mean up to 5

## 5. generalized gamma: Q near 1 supports the Weibull
summary(parsurv(Event(time, status) ~ treat + age, d, baseline = "gengamma"))
}
