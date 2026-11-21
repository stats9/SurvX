\name{cure}
\alias{cure}
\title{Cure models: survival with a cured fraction}
\description{
Fits a parametric cure model: a fraction of the population will never have the event (they are
"cured"), and the others have it after a time described by a parametric distribution (the
\emph{latency}). Covariates can act on the cured fraction and on the latency separately. Mixture
(Boag 1949; Berkson and Gage 1952; Farewell 1982) and non-mixture or promotion-time models
(Yakovlev and Tsodikov 1996; Chen, Ibrahim and Sinha 1999) are available, with delayed entry.

\strong{Use it when} the Kaplan-Meier curve levels off at a clearly positive value with long
follow-up (many cancers after curative treatment, recurrence after surgery), and the questions are
"what proportion is cured?" and "when do the others fail?".
\strong{Use something else when} there is no plateau, or follow-up is too short to see one: the
cured fraction is then not identified and the model only extrapolates. Use \code{\link{parsurv}}
or \code{\link{cox}}.
}
\usage{
cure(formula, data, cure = NULL, baseline = "weibull", scale = NULL,
     type = c("mixture", "nonmixture"), df = 3L, knots = NULL, cuts = NULL, init = NULL,
     max_iter = 200L, tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ latency covariates}, status 0/1;
  \code{Event(time, status, entry = entry)} for delayed entry. The right-hand side gives the
  covariates of the latency (the timing of the event among those who will have it).}
\item{data}{a data frame. Rows with missing values are removed.}
\item{cure}{one-sided formula with the covariates of the cured fraction, for example
  \code{~ age + stage}. \code{NULL} (default) uses the covariates of \code{formula}; \code{~ 1}
  gives one cured fraction for everybody. The two sets may differ, overlap or coincide.}
\item{baseline}{distribution of the latency: \code{"weibull"} (default), \code{"exponential"},
  \code{"gompertz"}, \code{"lognormal"}, \code{"loglogistic"}, \code{"gengamma"}, \code{"rp"} or
  \code{"pwexp"}, as in \code{\link{parsurv}}.}
\item{scale}{\code{"PH"} or \code{"AFT"} for the latency covariates, as in \code{\link{parsurv}}.}
\item{type}{\code{"mixture"} (default): the population is a mixture of cured subjects and
  subjects who will fail. \code{"nonmixture"}: the promotion-time model, in which the cumulative
  hazard of the whole population is bounded. See Details.}
\item{df, knots, cuts}{spline and piecewise exponential latencies, as in \code{\link{parsurv}}.}
\item{init}{optional starting values, in the order of \code{coef(fit)}.}
\item{max_iter, tol}{iteration limit and convergence tolerance.}
}
\details{
\subsection{The models}{
With \eqn{\pi(z) = \mathrm{expit}(\gamma^\top z)} the cured fraction and \eqn{S_u(t \mid x)} the
survival of the uncured (the latency distribution),
\describe{
  \item{mixture}{\eqn{S(t \mid x, z) = \pi(z) + \{1 - \pi(z)\} S_u(t \mid x)}.}
  \item{non-mixture}{\eqn{S(t \mid x, z) = \pi(z)^{F_u(t \mid x)}} with \eqn{F_u = 1 - S_u}: the
    cumulative hazard \eqn{-\log\pi(z) F_u(t \mid x)} is bounded and the hazard is proportional to
    \eqn{-\log \pi(z)}, as in the promotion time model. In this parametrisation (Lambert et al.
    2007, and \pkg{flexsurvcure}) \eqn{\pi} is again the cured fraction.}
}
In both, \eqn{S(t) \to \pi(z)} as \eqn{t \to \infty}.
}
\subsection{Estimation}{
Maximum likelihood by Newton-Raphson with exact derivatives (those of \code{\link{parsurv}} for
the latency; numerical for the generalized gamma); standard errors from the inverse observed
information. With delayed entry every contribution is divided by \eqn{S(entry \mid x, z)}.
Starting values: the latency fitted to the event times alone, and the logit of the Kaplan-Meier
estimate at the last time for the cured fraction.
}
\subsection{Identifiability}{
A cured fraction can be estimated only when follow-up is long enough for the survival curve to
level off, and when the latency distribution is plausible beyond the last events: the estimate of
\eqn{\pi} rests on the plateau (Maller and Zhou 1996; Amico and Van Keilegom 2018).
\code{plot(fit)} compares the fitted population survival with the Kaplan-Meier estimate and marks
the average cured fraction.
}
}
\section{Interpreting the coefficients}{
\describe{
  \item{Cured fraction (\code{cure: ...})}{log odds ratios: \eqn{e^\gamma} is the \strong{odds
    ratio of being cured}. OR 2.0 for a treatment means twice the odds of cure; OR < 1 means less
    cure. \code{predict(fit, newdata, type = "cured")} gives the cured fraction itself.}
  \item{Latency (the covariates of \code{formula})}{hazard ratios (PH) or time ratios (AFT)
    \emph{among the uncured} (mixture model): they say whether those who will have the event have
    it earlier or later, not whether they have it.}
  \item{Non-mixture model}{the latency coefficients describe the timing distribution \eqn{F_u};
    the cure coefficients change both the cured fraction and, proportionally, the hazard of the
    whole population.}
}
A covariate can increase the chance of cure and at the same time shorten the time to the event
among the uncured; report both parts. A sentence for a report: \emph{"In a Weibull mixture cure
model, the estimated cured fraction was 45\% (95\% CI 38 to 52\%) for stage I and 21\% (15 to
28\%) for stage II (odds ratio of cure 0.32); among patients who relapsed, relapse occurred later
with treatment (HR 0.70)."}
}
\section{Assumptions and how to check them}{
\describe{
  \item{A plateau exists}{\code{plot(km(...))}: the Kaplan-Meier curve levels off with many
    subjects still at risk; \code{plot(fit)} shows whether the fitted curve follows it.}
  \item{Latency distribution}{compare baselines with \code{AIC()}; the cured fraction can depend
    strongly on the distribution when the plateau is short, so report the sensitivity.}
}
}
\section{Messages and common problems}{
\describe{
  \item{\code{WARNING: the fit did NOT converge}}{often no plateau (the cured fraction drifts to 0)
    or too few events. Check the Kaplan-Meier curve; try a simpler latency or \code{cure = ~ 1}.}
  \item{A cured fraction near 0 or 1 with a huge standard error}{the cured fraction is not
    identified by these data.}
}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}; 62 data sets): the maximised log-likelihood agrees with
\code{flexsurvcure::flexsurvcure} to a relative \eqn{3 \times 10^{-10}} (mixture) and \eqn{3 \times
10^{-8}} (non-mixture), and predicted survival within the tolerance of its optimiser (\eqn{2 \times
10^{-4}}); the log-likelihood also agrees with independent code to \eqn{10^{-12}}. In a Monte Carlo
study (1000 replications, \eqn{n} = 200 to 1000) the 95\% intervals of the cure and latency
parameters covered 93.1 to 96.0\%. On \eqn{n = 10000}, a Weibull mixture cure fit takes about
0.03 s (\code{flexsurvcure}: 2.9 s).
}
\value{
An S7 object of class \code{CureFit}. \code{coef(fit)}: the latency coefficients (as in
\code{\link{parsurv}}) followed by those of the cured fraction, named \code{"cure: ..."}.

Methods:
\tabular{ll}{
  \code{print(fit)}, \code{summary(fit)} \tab latency (HR or TR) and cure (odds ratio) tables,
    average cured fraction, \code{Reading:} lines \cr
  \code{\link{estimates}(fit)}, \code{coef}, \code{vcov}, \code{confint}, \code{logLik},
    \code{AIC} \tab extractors \cr
  \code{predict(fit, newdata, times, type)} \tab \code{"survival"} (population),
    \code{"hazard"}, \code{"cumhaz"}, \code{"uncured"} (survival of the uncured) and
    \code{"cured"} (the cured fraction; no \code{times}); \code{se = TRUE} adds delta-method
    standard errors \cr
  \code{plot(fit)} \tab fitted population survival against Kaplan-Meier, with the cured fraction \cr
  \code{residuals(fit, type)} \tab Cox-Snell, martingale, deviance \cr
}
}
\references{
\emph{Mixture cure model:}
Boag JW (1949). Maximum likelihood estimates of the proportion of patients cured by cancer
  therapy. Journal of the Royal Statistical Society, Series B 11:15-44.
Berkson J, Gage RP (1952). Survival curve for cancer patients following treatment. Journal of the
  American Statistical Association 47:501-515.
Farewell VT (1982). The use of mixture models for the analysis of survival data with long-term
  survivors. Biometrics 38:1041-1046.

\emph{Non-mixture (promotion time) model:}
Yakovlev AY, Tsodikov AD (1996). Stochastic Models of Tumor Latency and Their Biostatistical
  Applications. World Scientific.
Chen MH, Ibrahim JG, Sinha D (1999). A new Bayesian model for survival data with a surviving
  fraction. Journal of the American Statistical Association 94:909-919.
Lambert PC, Thompson JR, Weston CL, Dickman PW (2007). Estimating and modeling the cure fraction
  in population-based cancer survival analysis. Biostatistics 8:576-594.

\emph{Identifiability and review:}
Maller RA, Zhou X (1996). Survival Analysis with Long-Term Survivors. Wiley.
Amico M, Van Keilegom I (2018). Cure models in survival analysis. Annual Review of Statistics and
  Its Application 5:311-342.

All references: \code{sx_references("cure")}.
}
\seealso{
\code{\link{parsurv}} (the latency distributions), \code{\link{km}} (look for the plateau),
\code{\link{predict.SurvX}}.
}
\examples{
set.seed(6)
n <- 400
d <- data.frame(age = rnorm(n), stage = rbinom(n, 1, 0.5))
cured <- runif(n) < plogis(0.5 - 1.2 * d$stage)                 # stage 1: less cure
t <- ifelse(cured, Inf, rweibull(n, 1.3, 2 * exp(-0.3 * d$age)))
cens <- runif(n, 0, 10)
d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)

## 1. look for a plateau
plot(km(Event(time, status) ~ stage, d))

## 2. mixture cure model: latency depends on age, cure on stage
fit <- cure(Event(time, status) ~ age, d, cure = ~ stage)
fit

## 3. cured fraction by stage, and survival of the uncured
nd <- data.frame(age = 0, stage = 0:1)
predict(fit, nd, type = "cured", se = TRUE)
predict(fit, nd, times = c(1, 3), type = "uncured")
plot(fit)

## 4. non-mixture model for comparison
AIC(fit, cure(Event(time, status) ~ age, d, cure = ~ stage, type = "nonmixture"))
}
