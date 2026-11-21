\name{predict.SurvX}
\alias{predict.SurvX}
\title{Predictions from SurvX models}
\description{
Predictions from every fitted model of SurvX: linear predictors, survival probabilities,
cumulative incidences, hazards, quantiles (median time), restricted mean survival time and cured
fractions, for new subjects or for the subjects of the data, optionally with standard errors and
confidence limits. Absolute predictions such as these are usually the clearest way to report a
model: "a 60-year-old patient on the new treatment has a 24\% risk of relapse within 5 years".
}
\section{Usage}{
\code{predict(object, newdata = NULL, times, type, cause = 1, se = FALSE, level = 0.95,
average = FALSE, p = 0.5)}; each model class accepts the arguments that apply to it.
}
\section{Arguments}{
\describe{
\item{\code{object}}{a fitted model.}
\item{\code{newdata}}{a data frame with one row per subject to predict for, containing every
  covariate of the model (and the \code{strata()} variables of a stratified Cox model). Factors may
  be given as character values of their levels. \code{NULL}: the subjects used for fitting.}
\item{\code{times}}{time points. For \code{directcif} fits \code{Inf} gives the plateau; for
  \code{parsurv} with \code{type = "rmst"} the upper limits of integration.}
\item{\code{type}}{for \code{cox}: \code{"lp"} (linear predictor) or \code{"survival"}. For
  \code{parsurv}: \code{"survival"}, \code{"hazard"}, \code{"cumhaz"}, \code{"quantile"} (at the
  probabilities \code{p}) or \code{"rmst"} (restricted mean survival time). For \code{csc} with a
  parametric baseline: \code{"cif"} or \code{"survival"} (event-free probability). For
  \code{cure}: \code{"survival"}, \code{"hazard"}, \code{"cumhaz"}, \code{"uncured"} (survival
  of the uncured) or \code{"cured"} (the cured fraction, without \code{times}).}
\item{\code{p}}{\code{parsurv}, \code{type = "quantile"}: probabilities; the time by which a
  fraction \code{p} has failed (\code{p = 0.5}: the median).}
\item{\code{cause}}{cause of interest (\code{csc} and \code{directcif} fits).}
\item{\code{se}}{for \code{directcif}, \code{parsurv} and parametric \code{csc} fits: also return
  delta-method standard errors (and, for \code{parsurv} and parametric \code{csc}, pointwise
  limits: complementary log-log scale for probabilities, log scale for positive quantities).}
\item{\code{level}}{confidence level of the limits.}
\item{\code{average}}{\code{parsurv} and parametric \code{csc}: the average over the rows of
  \code{newdata} (marginal curve), with its standard error.}
\item{\code{...}}{not used.}
}
}
\details{
Cox models: \eqn{S(t \mid x) = \exp\{-\hat H_0(t) e^{x'\hat\beta}\}} with the Breslow
estimate of the baseline cumulative hazard. Cause-specific Cox models: the cumulative incidence
\eqn{\int_0^t \hat S(u-)\, d\hat H_k(u \mid x)} from all cause-specific hazards (product-limit
or exponential form of \eqn{\hat S}, see \code{product_limit} of \code{\link{csc}}). Fine-Gray
models: \eqn{1 - \exp\{-\hat\Lambda_{10}(t) e^{x'\hat\beta}\}}. Direct and parametric models
evaluate the fitted distribution; parametric cause-specific models integrate the hazards
numerically (Gauss-Legendre). Standard errors of direct and parametric predictions use the delta
method with the covariance matrix of the parameters; limits are computed on a transformed scale
(log-log for probabilities, log for hazards and times) and transformed back.
}
\section{Choosing what to predict}{
\itemize{
  \item For one event type: survival probabilities at clinically meaningful times
    (\code{type = "survival"}), the median time (\code{parsurv}, \code{type = "quantile"}) or
    the restricted mean (\code{type = "rmst"}).
  \item With competing risks: the cumulative incidence of the cause of interest, from
    \code{fg}, \code{csc} or \code{directcif}; \code{csc} gives incidences of all causes that
    add up correctly with the event-free probability.
  \item \emph{Individual} predictions (one row of \code{newdata} per profile) describe typical
    patients; \emph{marginal} predictions (the average over the subjects of a population,
    \code{average = TRUE}, or \code{\link{cif_data}} with \code{type = "marginal"}) give
    standardised, population-level risks, for example to compare treatments as if everybody had
    received each.
  \item Cox and Fine-Gray baselines are step functions and are not defined beyond the last
    event time; parametric models give smooth curves and can extrapolate, with the caveat that the
    extrapolation rests on the distribution only.
}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}): Cox predictions agree with \code{survival::survfit.coxph}, Fine-Gray
predictions with \code{cmprsk::predict.crr} (\eqn{2 \times 10^{-8}}), cause-specific predictions
with \code{riskRegression::CSC} (\eqn{4 \times 10^{-7}}), and parametric predictions with
\code{flexsurv} within the tolerance of its optimiser. In a Monte Carlo study the 95\% intervals of
parametric survival predictions covered 93.5 to 94.9\%.
}
\references{
Breslow NE (1974) Covariance analysis of censored survival data. Biometrics 30:89-99.

Prentice RL, Kalbfleisch JD, Peterson AV, Flournoy N, Farewell VT, Breslow NE (1978) The analysis
  of failure times in the presence of competing risks. Biometrics 34:541-554.

Benichou J, Gail MH (1990) Estimates of absolute cause-specific risk in cohort studies. Biometrics
  46:813-826.

Fine JP, Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk.
  Journal of the American Statistical Association 94:496-509.

Kalbfleisch JD, Prentice RL (2002) The Statistical Analysis of Failure Time Data, 2nd ed. Wiley.

Hinchliffe SR, Lambert PC (2013) Flexible parametric modelling of cause-specific hazards to
  estimate cumulative incidence functions. BMC Medical Research Methodology 13:13.

Casella G, Berger RL (2002) Statistical Inference, 2nd ed. Duxbury, Section 5.5 (delta method).
}
\value{
a matrix with one row per row of \code{newdata} and one column per time point
  (a vector for \code{type = "lp"}); with \code{se = TRUE} a list with \code{cif} and \code{se}
  (\code{directcif}) or \code{est}, \code{se}, \code{lower} and \code{upper} (\code{parsurv},
  parametric \code{csc}).
}
\seealso{
\code{\link{cif_data}} and \code{\link{cif_plot}} (curves with bootstrap intervals),
\code{\link{hazard_plot}}, \code{\link{cox}}, \code{\link{fg}}, \code{\link{csc}},
\code{\link{parsurv}}, \code{\link{directcif}}, \code{\link{cure}}.
}
\examples{
set.seed(1)
n <- 300
d <- data.frame(treat = rbinom(n, 1, 0.5), age = round(rnorm(n, 60, 10)))
t <- rweibull(n, 1.4, 5 * exp(0.4 * d$treat - 0.02 * (d$age - 60))); cens <- runif(n, 0, 10)
d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
nd <- data.frame(treat = 0:1, age = 60)

## Cox: survival at 2 and 5 years for two profiles
fc <- cox(Event(time, status) ~ treat + age, d)
predict(fc, nd, times = c(2, 5), type = "survival")

## Weibull: survival with limits, median time, restricted mean up to 5 years
fw <- parsurv(Event(time, status) ~ treat + age, d)
predict(fw, nd, times = c(2, 5), se = TRUE)
predict(fw, nd, type = "quantile", p = 0.5)
predict(fw, nd, times = 5, type = "rmst")

## marginal (standardised) survival: everybody treated versus nobody treated
d1 <- transform(d, treat = 1); d0 <- transform(d, treat = 0)
c(treated = predict(fw, d1, times = 5, average = TRUE),
  untreated = predict(fw, d0, times = 5, average = TRUE))
}
