\name{directcif}
\alias{directcif}
\title{Direct parametric regression for cumulative incidence functions}
\description{
Maximum likelihood fit of the direct parametric models of Jeong and Fine (2006, 2007): the
cumulative incidence function of every cause is modelled directly,
\deqn{F_k(t; z) = 1 - \{1 + \alpha_k e^{z'\beta_k} u_k(t)\}^{-1/\alpha_k},}
with \eqn{\alpha_k \to 0} the proportional (subdistribution) hazards model and \eqn{\alpha_k = 1}
the proportional odds model. Three estimators are available (\code{estimation}).
}
\usage{
directcif(formula, data, baseline = "gompertz", link = "PH",
    estimation = c("joint", "weighted", "constrained"), cause = 1L,
    cond = c("logistic", "weibull", "loglogistic", "lognormal"), df = 3L,
    knots = NULL, constraint = FALSE, censoring_grid = 0L, init = NULL,
    max_iter = 500L, tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates} (or \code{~ 1}); status 0 = censored, 1..K = cause.}
\item{data}{a data frame.}
\item{baseline}{baseline family, one for all causes or one per cause.}
\item{link}{\code{"PH"} (alpha = 0), \code{"PO"} (alpha = 1), \code{"GOR"} (alpha estimated), or a
  number giving a fixed alpha; one for all causes or one per cause.}
\item{estimation}{\code{"joint"}, \code{"weighted"} or \code{"constrained"}; see Details.}
\item{cause}{cause of interest for \code{estimation = "weighted"} and \code{"constrained"}.}
\item{cond}{distribution \eqn{G} of the constrained model: \code{"logistic"}, \code{"weibull"},
  \code{"loglogistic"} or \code{"lognormal"}.}
\item{df}{degrees of freedom of the spline baseline.}
\item{knots}{optional list with the knots (log time, boundary knots first and last) of every cause.}
\item{constraint}{joint estimation only: keep the sum of the fitted CIFs below 1 within follow-up.}
\item{censoring_grid}{weighted estimation only: coarsen the estimated censoring distribution to at
  most this many support points (0 = exact). With the PH link the exact fit is linear in the
  sample size; with other links its cost grows with the number of competing events times the
  number of support points (a few seconds at n = 10000), and a grid can shorten it.}
\item{init}{optional start values on the estimation scale.}
\item{max_iter,tol}{iteration limit and convergence tolerance.}
}
\details{
Estimators:
\describe{
  \item{\code{"joint"}}{all causes jointly by the full likelihood of Jeong and Fine (2006, 2007).
    Appropriate when a share of the subjects never fails (the CIFs level off below a total of one).}
  \item{\code{"weighted"}}{the cause of interest alone. Subjects who fail from a competing cause
    enter through the Kaplan-Meier estimate of the censoring distribution, as in Fine and Gray
    (1999), so that no model for the other causes is needed (Geskus 2011; Lambert, Wilkes and
    Crowther 2017). Standard errors are from the sandwich estimator with one score per subject;
    the estimation of the censoring distribution is ignored. \code{logLik} is a weighted
    pseudo-likelihood: use it to compare links or baselines of weighted fits only.}
  \item{\code{"constrained"}}{two causes with asymptotes that add up to one (Shi, Cheng and Jeong
    2013): the cause of interest as above and \eqn{F_2(t; z) = \{1 - F_1(\infty; z)\} G(t)} with
    \eqn{G} the distribution of the time to the competing event given that it occurs
    (\code{cond}: modified logistic as in the paper, Weibull, log-logistic or log-normal). Needs a
    baseline with a plateau for the cause of interest (not \code{"weibull"}, not \code{"spline"}).
    Appropriate when every subject eventually fails from one of the two causes.}
}

Baseline families (\code{baseline}) and their parameters on the scale used for estimation:
\describe{
  \item{\code{"gompertz"}}{\eqn{u(t) = \tau (e^{\rho t} - 1)/\rho}; \code{rho, log_tau}. Jeong and Fine (2006).}
  \item{\code{"gompertz3"}}{three-parameter Gompertz; \code{rho, log_tau, zeta}. Haile et al. (2016).}
  \item{\code{"weibull"}}{\eqn{u(t) = (\lambda t)^\kappa}; \code{log_lambda, log_kappa}.}
  \item{\code{"mix_weibull"}, \code{"mix_loglogistic"}, \code{"mix_lognormal"}}{baseline CIF
    \eqn{F_0(t) = p\,G(t)} with \eqn{G} a Weibull, log-logistic or log-normal distribution function;
    \code{logit_p} and two parameters of \eqn{G}. Cheng (2009), Larson and Dinse (1985).}
  \item{\code{"logistic3"}}{modified three-parameter logistic baseline CIF; \code{log_b, c, logit_p}.
    Cheng (2009), Shi, Cheng and Jeong (2013).}
  \item{\code{"spline"}}{\eqn{\log u(t)} a restricted cubic spline in \eqn{\log t} with \code{df}
    degrees of freedom; \code{gamma0, gamma1, ...}. Royston and Parmar (2002), Mozumder et al. (2018).}
}
For baselines given as a CIF the baseline is the CIF at \eqn{z = 0}; centre continuous covariates.

The published models do not force \eqn{\sum_k F_k \le 1}. By default the model is fitted as
published and \code{fit@diagnostics} reports how many subjects have a fitted sum above 1 within
follow-up and at infinity. \code{constraint = TRUE} keeps the sum below 1 within follow-up.

A fit is reported as converged only at a stationary point with a positive definite information
matrix. Families with a plateau are not identified when the data show no plateau; such fits
are flagged and their standard errors should not be used.
}
\value{
an S7 object of class \code{DirectCIFFit}. Coefficients are named \code{"<cause>:<term>"}.
  Methods: \code{print}, \code{summary}, \code{coef}, \code{vcov}, \code{logLik} (so \code{AIC} and
  \code{BIC} work), \code{nobs}, \code{predict}, \code{plot}.
}
\references{
Jeong J-H, Fine JP (2006) Direct parametric inference for the cumulative incidence function.
J R Stat Soc C 55:187-200.

Jeong J-H, Fine JP (2007) Parametric regression on cumulative incidence function.
Biostatistics 8:184-196.

Cheng Y (2009) Modeling cumulative incidences of dementia and dementia-free death using a novel
three-parameter logistic function. Int J Biostat 5(1), Article 29.

Shi H, Cheng Y, Jeong J-H (2013) Constrained parametric model for simultaneous inference of two
cumulative incidence functions. Biom J 55:82-96.

Haile SR, Jeong J-H, Chen X, Cheng Y (2016) A 3-parameter Gompertz distribution for survival data
with competing risks, with an application to breast cancer data. J Appl Stat 43:2239-2253.

Royston P, Parmar MKB (2002) Flexible parametric proportional-hazards and proportional-odds
models for censored survival data. Stat Med 21:2175-2197.

Fine JP, Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk.
J Am Stat Assoc 94:496-509.

Geskus RB (2011) Cause-specific cumulative incidence estimation and the Fine and Gray model under
both left truncation and right censoring. Biometrics 67:39-49.

Lambert PC, Wilkes SR, Crowther MJ (2017) Flexible parametric modelling of the cause-specific
cumulative incidence function. Stat Med 36:1429-1446.

Mozumder SI, Rutherford MJ, Lambert PC (2018) Direct likelihood inference on the cause-specific
cumulative incidence function: a flexible parametric regression modelling approach.
Stat Med 37:82-97.
}
\examples{
set.seed(2)
n <- 400
x <- rbinom(n, 1, 0.5)
t1 <- rexp(n, 0.4 * exp(0.6 * x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
d <- data.frame(time = pmin(t1, t2, cens), x = x)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
fit <- directcif(Event(time, status) ~ x, d, baseline = "gompertz", link = "PH")
fit
predict(fit, data.frame(x = 0:1), times = c(1, 2, Inf))
fw <- directcif(Event(time, status) ~ x, d, estimation = "weighted", cause = 1)
fc <- directcif(Event(time, status) ~ x, d, estimation = "constrained", cause = 1, cond = "weibull")
nd <- data.frame(x = 0:1)
predict(fc, nd, times = Inf, cause = 1) + predict(fc, nd, times = Inf, cause = 2)
}
