\name{cure}
\alias{cure}
\title{Parametric cure models}
\description{
Mixture and non-mixture (promotion time) cure models: a fraction of the subjects never experiences
the event, and this cured fraction depends on covariates through a logistic model. The survival of
the subjects who are not cured (the latency) follows any distribution of \code{\link{parsurv}}, as
a proportional hazards or accelerated failure time model, with its own covariates.
}
\usage{
cure(formula, data, cure = NULL, baseline = "weibull", scale = NULL,
     type = c("mixture", "nonmixture"), df = 3L, knots = NULL, cuts = NULL, init = NULL,
     max_iter = 200L, tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ latency covariates}; status 0/1;
  \code{Event(time, status, entry = )} for delayed entry.}
\item{data}{a data frame.}
\item{cure}{one-sided formula of the covariates of the cured fraction, for example
  \code{~ age + stage}; \code{NULL} (default) uses the covariates of \code{formula}, and \code{~ 1}
  gives one cured fraction for all subjects.}
\item{baseline}{distribution of the latency: \code{"weibull"} (default), \code{"exponential"},
  \code{"gompertz"}, \code{"lognormal"}, \code{"loglogistic"}, \code{"gengamma"}, \code{"rp"} or
  \code{"pwexp"}, as in \code{\link{parsurv}}.}
\item{scale}{\code{"PH"} or \code{"AFT"} for the latency, as in \code{\link{parsurv}}.}
\item{type}{\code{"mixture"} (default) or \code{"nonmixture"}.}
\item{df,knots,cuts}{spline and piecewise exponential latencies, as in \code{\link{parsurv}}.}
\item{init}{optional starting values, in the order of \code{coef(fit)}.}
\item{max_iter,tol}{iteration limit and convergence tolerance.}
}
\details{
With \eqn{\pi(z) = \mathrm{expit}(\gamma^\top z)} the cured fraction and \eqn{S_u(t \mid x)} the
survival of the uncured,
\describe{
  \item{mixture}{\eqn{S(t \mid x, z) = \pi(z) + \{1 - \pi(z)\} S_u(t \mid x)}: the population is a
    mixture of cured subjects and subjects who will fail (Boag 1949; Berkson and Gage 1952;
    Farewell 1982). Latency effects (HR or TR) apply to the uncured only.}
  \item{nonmixture}{\eqn{S(t \mid x, z) = \pi(z)^{F_u(t \mid x)}}, with \eqn{F_u = 1 - S_u}: the
    cumulative hazard \eqn{-\log\pi(z) F_u(t \mid x)} is bounded, as in the promotion time model of
    Yakovlev and Tsodikov (1996) and Chen, Ibrahim and Sinha (1999); the hazard is proportional
    to \eqn{-\log \pi(z)}. This parametrisation, in which \eqn{\pi} is again the cured fraction,
    is that of Lambert et al. (2007) and of \pkg{flexsurvcure}.}
}
The likelihood is maximised by Newton-Raphson in the C++ core with exact derivatives (the latency
derivatives are those of \code{\link{parsurv}}; numerical for the generalized gamma), and
standard errors are from the inverse observed information. With delayed entry every contribution
is divided by \eqn{S(entry \mid x, z)}. Starting values: the latency fitted to the event times
alone and the logit of the Kaplan-Meier estimate at the last time for the cured fraction.

A cured fraction is identified only when follow-up is long enough for the survival curve to level
off; \code{plot(fit)} compares the fitted marginal survival with the Kaplan-Meier estimate and
shows the average cured fraction (Maller and Zhou 1996; Amico and Van Keilegom 2018).

Validation: the log-likelihood and its maximum agree with independent code to 1e-12 (Weibull and
log-normal latency, both types, with and without delayed entry), and the tests compare the fit with
\code{flexsurvcure::flexsurvcure}.
}
\value{
an S7 object of class \code{CureFit}. \code{coef(fit)}: the latency coefficients (as in
  \code{\link{parsurv}}) followed by those of the cured fraction (\code{"cure: ..."}, log odds
  ratios of being cured). \code{predict(fit, newdata, times, type)} with \code{type} one of
  \code{"survival"} (population), \code{"hazard"}, \code{"cumhaz"}, \code{"uncured"} (survival of
  the uncured) and \code{"cured"} (the cured fraction; no \code{times}), with delta-method
  standard errors when \code{se = TRUE}. Also \code{print}, \code{summary}, \code{estimates},
  \code{confint}, \code{vcov}, \code{logLik}, \code{AIC}, \code{residuals}, \code{plot}.
}
\references{
Boag JW (1949) Maximum likelihood estimates of the proportion of patients cured by cancer therapy.
  Journal of the Royal Statistical Society, Series B 11:15-53.

Berkson J, Gage RP (1952) Survival curve for cancer patients following treatment. Journal of the
  American Statistical Association 47:501-515.

Farewell VT (1982) The use of mixture models for the analysis of survival data with long-term
  survivors. Biometrics 38:1041-1046.

Yakovlev AY, Tsodikov AD (1996) Stochastic Models of Tumor Latency and Their Biostatistical
  Applications. World Scientific.

Chen MH, Ibrahim JG, Sinha D (1999) A new Bayesian model for survival data with a surviving
  fraction. Journal of the American Statistical Association 94:909-919.

Lambert PC, Thompson JR, Weston CL, Dickman PW (2007) Estimating and modeling the cure fraction in
  population-based cancer survival analysis. Biostatistics 8:576-594.

Maller RA, Zhou X (1996) Survival Analysis with Long-Term Survivors. Wiley.

Amico M, Van Keilegom I (2018) Cure models in survival analysis. Annual Review of Statistics and
  Its Application 5:311-342.
}
\examples{
set.seed(6)
n <- 400
d <- data.frame(age = rnorm(n), stage = rbinom(n, 1, 0.5))
cured <- runif(n) < plogis(0.5 - 1.2 * d$stage)
t <- ifelse(cured, Inf, rweibull(n, 1.3, 2 * exp(-0.3 * d$age)))
cens <- runif(n, 0, 10)
d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
fit <- cure(Event(time, status) ~ age, d, cure = ~ stage)
fit
predict(fit, data.frame(age = 0, stage = 0:1), type = "cured", se = TRUE)
plot(fit)
}
