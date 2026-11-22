\name{adjusted_risk}
\alias{adjusted_risk}
\alias{AdjustedRisk}
\title{Adjusted absolute risks and causal contrasts between treatments}
\description{
Estimates, for each value of a treatment (or exposure), the absolute risk of the event, or the
cumulative incidence of a cause, that would have been seen had everyone in the data received that
value, and the differences and ratios between treatments, at chosen times and as curves over
follow-up. Two methods are available and can be compared: standardisation of a fitted outcome
model (the g-formula) and inverse probability of treatment weighting. Intervals come from a
bootstrap in which the models are refitted.

\strong{Use it when} the question is "how much would the risk change if everyone were treated
rather than not", in an observational study or to turn the hazard ratios of a trial into absolute
risks adjusted for covariates. Hazard ratios do not answer that question; absolute risks do.

\strong{Use something else when} you want the risk of a particular patient (\code{predict()} of
the model) or a summary that needs no model (\code{\link{km}}, \code{\link{aj}} by group).
}
\usage{
adjusted_risk(fit, data, treatment, times = NULL, cause = 1L,
    method = c("gformula", "ipw", "both"), ps_formula = NULL, trim = 0,
    B = 200L, seed = 1L, level = 0.95, reference = NULL)
}
\arguments{
\item{fit}{the outcome model: a fit of \code{\link{cox}}, \code{\link{fg}}, \code{\link{csc}},
  \code{\link{directcif}}, \code{\link{parsurv}}, \code{\link{cure}} or the result of
  \code{\link{time_effects}}. It must contain the treatment as a covariate (or, for \code{cox} and
  \code{csc}, as a stratum) and should contain the confounders. It must have been fitted with
  \code{data = }, so that it can be refitted in bootstrap samples. A \code{cox()} or \code{parsurv()}
  model with \code{+ frailty(id)} is standardised with the frailty integrated out (the marginal
  risk of a subject from a new cluster, \code{predict(fit, frailty = "marginal")}); in the bootstrap
  the clusters are resampled whole, and a cluster drawn twice enters the refit as two clusters.}
\item{data}{the data frame the model was fitted to.}
\item{treatment}{the name of the treatment variable (a column of \code{data}): a factor, character,
  logical or 0/1 variable with at most 10 values.}
\item{times}{times at which risks and contrasts are reported. Default: the quartiles of the event
  times of the cause (two significant digits). Curves are computed on a grid over follow-up in
  addition.}
\item{cause}{competing risks: the cause whose cumulative incidence is compared (for \code{fg} the
  cause of the model).}
\item{method}{\code{"gformula"} (default): standardisation of the outcome model.
  \code{"ipw"}: inverse probability of treatment weighting with a model for the treatment.
  \code{"both"}: both, side by side.}
\item{ps_formula}{the model of the treatment for \code{"ipw"}: \code{treatment ~ covariates}.
  Default: the treatment on all other terms of the outcome model (main effects).}
\item{trim}{truncation of the weights at the \code{trim} and \code{1 - trim} quantiles (for example
  0.01) when some weights are extreme; 0 (default) for none.}
\item{B}{number of bootstrap samples; 0 for point estimates only.}
\item{seed}{seed of the bootstrap: a whole number, or \code{NULL} to take the seed from the random numbers of the session (so that \code{set.seed()} before the call governs the result). All random draws come from R's generator (see \code{RNGkind()}); with a number the state of the session's generator is left as it was.}
\item{level}{confidence level.}
\item{reference}{the treatment value against which the others are compared; default the first
  level.}
}
\details{
\subsection{Standardisation (g-formula)}{
For each treatment value \eqn{a}, the outcome model predicts the risk \eqn{F(t \mid a, z_i)} of
every subject with the treatment set to \eqn{a} and the other covariates \eqn{z_i} as observed, and
the predictions are averaged: \eqn{\hat F_a(t) = n^{-1} \sum_i \hat F(t \mid a, z_i)}
(Robins 1986; Hernan and Robins 2020, Chapter 13). For one event type the risk is one minus the
predicted survival; for competing risks the predicted cumulative incidence of the cause. The model
can contain interactions of the treatment with covariates and, through \code{\link{time_effects}},
time-varying effects.
}
\subsection{Inverse probability of treatment weighting}{
A logistic (two values) or multinomial logistic (more values) model gives the probability
\eqn{e_a(z)} of each treatment; subjects get the stabilised weight \eqn{P(A = a)/e_A(z)}, and the
weighted Kaplan-Meier (one event type) or Aalen-Johansen estimate (competing risks) is computed in
each treatment group (Cole and Hernan 2004; Neumann and Billionnet 2016). Delayed entry of the
outcome model is used. Balance is shown by standardised mean differences of the covariates before
and after weighting, and the weights by the effective sample size and the largest weight (Austin
and Stuart 2015).
}
\subsection{Contrasts and intervals}{
Risk differences \eqn{\hat F_a(t) - \hat F_{ref}(t)} and ratios \eqn{\hat F_a(t)/\hat F_{ref}(t)}.
With \code{B > 0} the subjects (clusters, if the model has them) are resampled, the outcome model is
refitted with its own call (and the treatment model refitted), and percentile intervals are
reported; the p-value of a difference uses the bootstrap standard error. Samples in which a
treatment value is missing or a model fails are skipped and counted.
}
}
\section{Interpreting the result}{
\describe{
  \item{Risk by treatment}{the risk that would have been observed by time \eqn{t} had everyone
    received that treatment, in a population with the covariate distribution of the data.}
  \item{Risk difference}{in percentage points: the change in absolute risk from treating everyone
    with one value rather than the reference. Its inverse is a number needed to treat.}
  \item{Competing risks}{the contrast is a total effect on the cumulative incidence, which includes
    the effect of treatment on the competing causes: a treatment that prevents the competing event
    can increase the incidence of the cause of interest.}
  \item{\code{method = "both"}}{the two methods rest on different models, the outcome model and the
    treatment model. Close agreement is reassuring; a large difference means that at least one model
    is wrong.}
}
A sentence for a report: \emph{"Standardised to the covariate distribution of the cohort, the
5-year risk of death was 18.2\% had all patients been treated and 24.9\% had none been treated, a
difference of -6.7 percentage points (95\% bootstrap CI -10.9 to -2.5)."}
}
\section{Assumptions and how to check them}{
\describe{
  \item{No unmeasured confounding}{every common cause of treatment and outcome is among the
    covariates. This cannot be checked from the data.}
  \item{Positivity}{every subject could have received every treatment: look at the range of the
    propensity scores and the largest weights (\code{method = "ipw"} or \code{"both"}).}
  \item{Correct models}{the outcome model for the g-formula (check it with \code{\link{ph_test}},
    \code{\link{time_effects}} and \code{\link{model_scan}}), the treatment model for IPW
    (balance after weighting).}
  \item{Baseline covariates}{one row per subject; the g-formula with covariates that change over
    follow-up needs methods not provided here.}
}
}
\section{Messages and common problems}{
\describe{
  \item{\code{data does not match the model}}{give the same data frame that the model was fitted
    to; rows with missing values in the variables of the model are dropped as in the fit.}
  \item{\code{the model call has no 'data' argument}}{fit the model with \code{data = } so that it
    can be refitted in bootstrap samples, or use \code{B = 0}.}
  \item{Long computing times}{the model is refitted \code{B} times; check with \code{B = 0} first.}
}
}
\section{Validation}{
With a binary treatment and no other covariates both methods reduce to the Kaplan-Meier or
Aalen-Johansen estimate in each group, which is checked in the test suite, as are the weighted
estimator against the weighted estimates of \pkg{survival} and the multinomial model against
\code{stats::glm} for two values.
}
\value{
An S7 object of class \code{AdjustedRisk} with \code{x@risk} (method, level, time, risk, se, lower,
upper), \code{x@contrast} (method, level, reference, time, difference, se, lower, upper, ratio,
ratio_lower, ratio_upper, p), \code{x@curve} (the risks on a grid over follow-up) and \code{x@ps}
(propensity score ranges, weights and balance for IPW).

Methods:
\tabular{ll}{
  \code{print(x)} \tab risks, contrasts and a reading \cr
  \code{summary(x)} \tab also the balance of every covariate \cr
  \code{plot(x, method = NULL)} \tab adjusted risk curves with their bands \cr
  \code{as.data.frame(x, what)} \tab \code{"risk"}, \code{"contrast"}, \code{"curve"} or \code{"balance"} \cr
  \code{report(x, file)} \tab a self-contained HTML report \cr
  \code{\link{sx_references}(x)} \tab the references of the methods \cr
}
}
\references{
\emph{Standardisation and causal inference}

Robins J (1986) A new approach to causal inference in mortality studies with a sustained exposure
period: application to control of the healthy worker survivor effect. Math Model 7:1393-1512.

Hernan MA, Robins JM (2020) Causal Inference: What If. Chapman and Hall/CRC.

Ozenne BMH, Scheike TH, Staerk L, Gerds TA (2020) On the estimation of average treatment effects
with right-censored time to event outcome and competing risks. Biom J 62:751-763.

\emph{Inverse probability of treatment weighting}

Cole SR, Hernan MA (2004) Adjusted survival curves with inverse probability weights. Comput Methods
Programs Biomed 75:45-49.

Neumann A, Billionnet C (2016) Covariate adjustment of cumulative incidence functions for
competing risks data using inverse probability of treatment weighting. Comput Methods Programs
Biomed 129:63-70.

Austin PC, Stuart EA (2015) Moving towards best practice when using inverse probability of
treatment weighting (IPTW) using the propensity score to estimate causal treatment effects in
observational studies. Stat Med 34:3661-3679.

\emph{Bootstrap}

Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap. Chapman and Hall.

All references: \code{sx_references("adjusted_risk")}.
}
\seealso{
\code{\link{predict.SurvX}} (risks of given subjects), \code{\link{time_effects}},
\code{\link{km}}, \code{\link{aj}}, \code{\link{rmst}}.
}
\examples{
set.seed(2)
n <- 400
age <- rnorm(n)
trt <- rbinom(n, 1, plogis(0.8 * age))                 # older patients treated more often
t <- rexp(n, 0.2 * exp(-0.5 * trt + 0.7 * age)); cens <- runif(n, 0, 8)
d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens),
                trt = factor(trt, 0:1, c("control", "treated")), age = age)

fit <- cox(Event(time, status) ~ trt + age, data = d)
ar <- adjusted_risk(fit, d, "trt", times = c(2, 5), method = "both", B = 0)
ar
plot(ar)

\donttest{
## with bootstrap intervals
ar_b <- adjusted_risk(fit, d, "trt", times = c(2, 5), B = 200)
ar_b
report(ar_b, file = tempfile(fileext = ".html"))

## competing risks: cumulative incidence of live birth by education
data(fertility)
adjusted_risk(fg(Event(time, status) ~ education, data = fertility, cause = 1),
              fertility, "education", times = c(1, 2, 5), B = 50)
}
}
