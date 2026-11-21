\name{csc}
\alias{csc}
\title{Cause-specific hazards regression for competing risks}
\description{
Fits one hazard model per cause, each describing the rate of that cause among subjects still free
of all events, and combines them into the cumulative incidence of every cause. With
\code{baseline = "cox"} (default) each cause gets a Cox model; with a parametric
\code{baseline} each cause gets a parametric hazard model as in \code{\link{parsurv}}.

\strong{Use it when} there are competing events and the question is about the mechanism: how a
covariate changes the rate at which each cause occurs (aetiology), or when absolute risks of all
causes are needed that are coherent (they never add up to more than 1).
\strong{Use something else when}
\itemize{
  \item the effect of covariates on the cumulative incidence of one cause is wanted as a single
    number: \code{\link{fg}} (or \code{\link{directcif}}). Reporting both kinds of model is
    recommended (Latouche et al. 2013);
  \item there is a single event type: \code{\link{cox}} or \code{\link{parsurv}}.
}
}
\usage{
csc(formula, data, ties = c("efron", "breslow"), product_limit = TRUE, baseline = "cox",
    scale = NULL, df = 3L, knots = NULL, cuts = NULL, max_iter = 200L, tol = 1e-10,
    robust = FALSE)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}, status 0 = censored, 1, 2, ..., K =
  cause. The same covariates enter the model of every cause. \code{Event(time, status, entry =
  entry)} for delayed entry and, with a Cox baseline, for (start, stop] rows with time-varying
  covariates. With a Cox baseline the right-hand side may contain \code{strata(...)} (a baseline
  hazard per stratum for every cause) and \code{cluster(id)} (robust variance), as in
  \code{\link{cox}}.}
\item{data}{a data frame. Rows with missing values are removed.}
\item{ties}{Cox baseline: \code{"efron"} (default) or \code{"breslow"}, as in \code{\link{cox}}.}
\item{product_limit}{Cox baseline: how the probability of being free of all events is computed
  for the cumulative incidence. \code{TRUE} (default): as a product-limit
  \eqn{\prod (1 - d\Lambda)}, as \pkg{riskRegression} does; \code{FALSE}:
  \eqn{\exp\{-\sum_k \Lambda_k(t)\}}. The difference is small except with few subjects at risk.}
\item{baseline}{\code{"cox"}, or one of the distributions of \code{\link{parsurv}}
  (\code{"exponential"}, \code{"weibull"}, \code{"gompertz"}, \code{"lognormal"},
  \code{"loglogistic"}, \code{"gengamma"}, \code{"rp"}, \code{"pwexp"}). A vector gives one
  distribution per cause, for example \code{c("weibull", "gompertz")}. Cox and parametric
  baselines cannot be mixed.}
\item{scale}{parametric baselines: \code{"PH"} or \code{"AFT"}, one value or one per cause;
  default as in \code{\link{parsurv}}.}
\item{df, knots, cuts}{spline (\code{"rp"}) and piecewise exponential (\code{"pwexp"})
  baselines, as in \code{\link{parsurv}}; by default knots and cut points are chosen from the
  events of each cause.}
\item{max_iter, tol}{iteration limit and convergence tolerance.}
\item{robust}{Cox baseline: robust (sandwich) variance of Lin and Wei (1989) for every cause;
  always used with \code{cluster()}.}
}
\details{
\subsection{The model}{
The cause-specific hazard of cause \eqn{k},
\eqn{h_k(t) = \lim P(t \le T < t + dt, \mathrm{cause} = k \mid T \ge t)/dt}, is modelled as
\eqn{h_k(t \mid x) = h_{k0}(t) \exp(\beta_k^\top x)} (Prentice et al. 1978). Each cause is fitted
separately, with the subjects who had another cause censored at their event time: the likelihood
factorises over causes, so this is the maximum likelihood estimate of the joint model, and each
cause-specific Cox model is identical to \code{survival::coxph} on \code{Surv(time, status == k)}.
}
\subsection{Cumulative incidence}{
\deqn{F_k(t \mid x) = \int_0^t S(u- \mid x)\, h_k(u \mid x)\,du, \qquad
      S(t \mid x) = \exp\Big\{-\sum_j \Lambda_j(t \mid x)\Big\}}
(Kalbfleisch and Prentice 2002, Section 8.2). With a Cox baseline this is a sum over event times
with the Breslow baseline hazards (Benichou and Gail 1990), as in \code{riskRegression::CSC}; for
a subject in stratum \eqn{s}, the baseline hazards of that stratum are used. With a parametric
baseline the integral is computed numerically (20-point Gauss-Legendre on a fine partition of the
time axis), with delta-method standard errors (Hinchliffe and Lambert 2013). Either way the
cumulative incidences of all causes and the event-free probability add up to one.
}
}
\section{Interpreting the coefficients}{
\itemize{
  \item The printed effect for cause \eqn{k} is the \strong{cause-specific hazard ratio}: the
    factor by which a covariate multiplies the rate of cause \eqn{k} \emph{among subjects who are
    still free of every event}. It is read like the HR of \code{\link{cox}}.
  \item A cause-specific HR does not translate into the risk of the cause: a covariate that does
    not change the hazard of cause 1 but lowers that of cause 2 increases the cumulative incidence
    of cause 1 (subjects survive longer to experience it). Look at the effects on \emph{all}
    causes together, and at the predicted cumulative incidences:
    \code{predict(fit, newdata, times, cause)} or \code{plot(fit, newdata)}.
  \item In aetiological questions (does exposure X cause the disease?) the cause-specific HR is the
    quantity of interest; for prognosis (what is the risk of this patient?) the cumulative
    incidence is (Andersen et al. 2012; Latouche et al. 2013).
}
A sentence for a report: \emph{"Treatment reduced the rate of relapse (cause-specific HR 0.65,
95\% CI 0.50 to 0.85) and did not change the rate of death without relapse (HR 1.05, 0.75 to
1.47); the predicted 5-year cumulative incidence of relapse was 22\% versus 31\%."}
}
\section{Assumptions and how to check them}{
\describe{
  \item{Proportional cause-specific hazards}{one Cox model per cause: check each with
    \code{ph_test(fit)}, which tests every cause.}
  \item{Parametric baselines}{\code{hazard_plot(fit)} compares the fitted and non-parametric
    cause-specific hazards; \code{AIC()} compares distributions (for all causes together).}
  \item{Independent censoring}{as in \code{\link{cox}}; the competing events themselves are
    not censoring, and need no independence assumption for the cause-specific hazards.}
}
}
\section{Messages and common problems}{
\describe{
  \item{\code{WARNING: the fit did NOT converge} (for a cause)}{usually a cause with very few
    events, or separation; report that cause with care or merge rare causes.}
  \item{\code{baseline: 'cox' cannot be mixed with parametric distributions}}{choose all Cox or all
    parametric baselines.}
  \item{\code{at least one covariate is needed}}{for no covariates use \code{\link{aj}}.}
}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}; 82 data sets with and without ties, strata and delayed entry): coefficients
agree with \code{survival::coxph} for each cause to \eqn{6 \times 10^{-7}} and variances to
\eqn{4 \times 10^{-10}}; predicted cumulative incidences agree with \code{riskRegression::CSC} to
\eqn{4 \times 10^{-7}} (stratified: \eqn{5 \times 10^{-10}}). In a Monte Carlo study (1000
replications, \eqn{n} = 200 to 1000) the 95\% intervals covered 94.5 to 97.3\%. On
\eqn{n = 10000}, a fit takes about 0.01 s (\code{riskRegression::CSC}: 0.07 s).
}
\value{
With \code{baseline = "cox"}: an S7 object of class \code{CSCFit}; with a parametric baseline:
class \code{ParCSCFit}. In both, \code{fit@models} holds the fit of each cause and \code{fit@K} the
number of causes; \code{coef(fit)} and \code{vcov(fit)} are lists with one element per cause;
\code{logLik()} and \code{AIC()} refer to all causes together.

Methods:
\tabular{ll}{
  \code{print(fit)}, \code{summary(fit)} \tab one table of cause-specific HRs per cause, with
    likelihood ratio tests and \code{Reading:} lines \cr
  \code{\link{estimates}(fit)} \tab all causes in one data frame, with a \code{cause} column \cr
  \code{predict(fit, newdata, times, cause)} \tab cumulative incidence of the cause; parametric
    fits also \code{type = "survival"} (free of all events) and \code{se = TRUE} \cr
  \code{plot(fit)}, \code{\link{cif_plot}}, \code{\link{cif_data}} \tab cumulative incidence
    curves of all causes, with bootstrap (Cox) or delta-method (parametric) intervals \cr
  \code{\link{ph_test}(fit)}, \code{residuals(fit)} \tab checks per cause \cr
  \code{\link{hazard_plot}(fit)} \tab parametric fits: fitted against non-parametric hazards \cr
}
}
\references{
\emph{Cause-specific hazards:}
Prentice RL, Kalbfleisch JD, Peterson AV, Flournoy N, Farewell VT, Breslow NE (1978). The
  analysis of failure times in the presence of competing risks. Biometrics 34:541-554.
Kalbfleisch JD, Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd ed.
  Wiley, Chapter 8.

\emph{Cox models per cause:}
Cox DR (1972). Regression models and life-tables. Journal of the Royal Statistical Society,
  Series B 34:187-220.
Breslow NE (1974). Covariance analysis of censored survival data. Biometrics 30:89-99.
Efron B (1977). The efficiency of Cox's likelihood function for censored data. Journal of the
  American Statistical Association 72:557-565.
Andersen PK, Gill RD (1982). Cox's regression model for counting processes: a large sample study.
  Annals of Statistics 10:1100-1120.
Lin DY, Wei LJ (1989). The robust inference for the Cox proportional hazards model. Journal of
  the American Statistical Association 84:1074-1078.

\emph{Absolute risks:}
Benichou J, Gail MH (1990). Estimates of absolute cause-specific risk in cohort studies.
  Biometrics 46:813-826.
Hinchliffe SR, Lambert PC (2013). Flexible parametric modelling of cause specific hazards to
  estimate cumulative incidence functions. BMC Medical Research Methodology 13:13.

\emph{Parametric models:}
Lawless JF (2003). Statistical Models and Methods for Lifetime Data, 2nd ed. Wiley.

\emph{Interpretation:}
Latouche A, Allignol A, Beyersmann J, Labopin M, Fine JP (2013). A competing risks analysis should
  report results on all cause-specific hazards and cumulative incidence functions. Journal of
  Clinical Epidemiology 66:648-653.
Andersen PK, Geskus RB, de Witte T, Putter H (2012). Competing risks in epidemiology:
  possibilities and pitfalls. International Journal of Epidemiology 41:861-870.

All references: \code{sx_references("csc")}.
}
\seealso{
\code{\link{fg}} and \code{\link{directcif}} (models for the cumulative incidence),
\code{\link{aj}}, \code{\link{cif_plot}}, \code{\link{cox}}, \code{\link{parsurv}},
\code{\link{predict.SurvX}}.
}
\examples{
set.seed(2)
n <- 300
d <- data.frame(treat = rbinom(n, 1, 0.5), age = round(rnorm(n, 60, 10)))
t1 <- rweibull(n, 1.3, 3 * exp(0.4 * d$treat))        # relapse, delayed by treatment
t2 <- rweibull(n, 0.9, 5 * exp(-0.02 * (d$age - 60)))  # death without relapse
cens <- runif(n, 0, 6)
d$time <- pmin(t1, t2, cens)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))

## 1. cause-specific Cox models for both causes
fit <- csc(Event(time, status) ~ treat + age, d)
fit

## 2. absolute risk of relapse at 1 and 3 years
nd <- data.frame(treat = 0:1, age = 60)
predict(fit, nd, times = c(1, 3), cause = 1)

## 3. parametric cause-specific hazards, one distribution per cause
fw <- csc(Event(time, status) ~ treat + age, d, baseline = c("weibull", "weibull"))
predict(fw, nd, times = 3, cause = 1, se = TRUE)
AIC(fw, csc(Event(time, status) ~ treat + age, d, baseline = "lognormal"))
}
