\name{fg}
\alias{fg}
\title{Fine-Gray regression for the cumulative incidence of a competing risk}
\description{
Fits the proportional subdistribution hazards model of Fine and Gray (1999): the effect of
covariates on the \emph{cumulative incidence} of one cause (the probability of having that event
by time \eqn{t} when other events compete with it). Inverse-probability-of-censoring weights, the
robust variance of Fine and Gray, and delayed entry (Geskus 2011) are supported.

\strong{Use it when} there are several event types (status 0 = censored, 1, 2, ... = causes) and
the question is how covariates change the \emph{risk} (absolute probability) of one cause, for
prognosis or prediction.
\strong{Use something else when}
\itemize{
  \item the question is about the mechanism, the rate at which the cause occurs among those still
    event-free: \code{\link{csc}} (cause-specific hazards). Reporting both is recommended
    (Latouche et al. 2013);
  \item a parametric cumulative incidence is wanted, a proportional odds link, or estimates for
    all causes jointly: \code{\link{directcif}};
  \item there is only one event type: \code{\link{cox}};
  \item no covariates, only groups: \code{\link{aj}} and \code{\link{gray_test}}.
}
}
\usage{
fg(formula, data, cause = 1L, variance = c("robust", "model", "bootstrap"), max_iter = 50L,
    tol = 1e-10, B = 200L, seed = 1L)
}
\arguments{
\item{formula}{a formula with \code{\link{Event}()} on the left and the covariates on the right:
  \tabular{ll}{
    \code{Event(time, status)} \tab one row per subject, right-censored \cr
    \code{Event(time, status, entry = entry)} \tab delayed entry (left truncation) \cr
  }
  \code{status} is 0 for censored and 1, 2, ... for the event types (a factor with the first
  level meaning censored also works, see \code{\link{Event}}). Covariates as in
  \code{\link[stats]{lm}}: numeric, factors, interactions, transformations, splines.
  \code{strata()} and \code{cluster()} are not available in \code{fg()}.}
\item{data}{a data frame with the variables of the formula. Rows with missing values in any of
  them are removed.}
\item{cause}{the event type of interest (default 1). All other non-zero codes are competing
  events.}
\item{variance}{\code{"robust"} (default): the sandwich variance of Fine and Gray (1999,
  Section 4), which accounts for the estimated censoring weights; use it for inference without
  delayed entry. \code{"bootstrap"}: the variance of the estimates over \code{B} bootstrap samples
  of the subjects; recommended with delayed entry, where the robust variance treats the weights as
  known and is too small (see Validation). \code{"model"}: the inverse information of the
  weighted partial likelihood, kept for comparison; not valid for inference in general. The
  robust and model-based variances are always stored in the fit.}
\item{max_iter}{maximum number of Newton-Raphson iterations (default 50).}
\item{tol}{convergence tolerance on the relative change of the log pseudo-likelihood.}
\item{B}{\code{variance = "bootstrap"}: number of bootstrap samples (default 200; 500 or more for
  final results).}
\item{seed}{\code{variance = "bootstrap"}: seed of the resampling, so that the same call gives
  the same standard errors. The random number generator of the session is restored afterwards.}
}
\details{
\subsection{The model}{
The subdistribution hazard of cause 1 is \eqn{\lambda_1(t) = -d \log\{1 - F_1(t)\} / dt}, with
\eqn{F_1(t) = P(T \le t, \mathrm{cause} = 1)} the cumulative incidence. The model is
\deqn{\lambda_1(t \mid x) = \lambda_{10}(t) \exp(\beta^\top x), \qquad
      F_1(t \mid x) = 1 - \exp\{-\Lambda_{10}(t) e^{\beta^\top x}\}.}
Subjects who have a competing event stay in the risk set after their event, with a weight that
decreases as the probability of still being under observation decreases.
}
\subsection{Estimation}{
Without delayed entry: the weighted partial likelihood of Fine and Gray (1999) with weights
\eqn{\hat G(t-)/\hat G(T_i-)}, where \eqn{\hat G} is the Kaplan-Meier estimate of the censoring
distribution; Breslow handling of tied event times. The computation uses a forward-backward scan
in linear time (as proposed by Kawaguchi et al. 2021), which is why \code{fg()} is fast on large
data. The robust variance is the sandwich of Fine and Gray (1999), with the term for the estimated
censoring distribution; ties between censoring and event times follow \code{cmprsk::crr}.

With delayed entry: the weights of Geskus (2011),
\eqn{\hat G(t-)\hat H(t-)/\{\hat G(T_i-)\hat H(T_i-)\}}, with \eqn{\hat G} the product-limit
estimate of the censoring distribution among subjects under observation and \eqn{\hat H} the
reverse-time product-limit estimate of the entry distribution. The model is then a weighted Cox
model on (start, stop] rows cut where the weights change, and the robust variance is the sandwich
over subjects (Lin and Wei 1989) with the weights taken as known, as in
\code{survival::finegray} followed by \code{coxph(..., cluster = id)}.
}
}
\section{Standard errors with delayed entry}{
With delayed entry the robust variance treats the Geskus weights as known. In the Monte Carlo
study of the package it was 5 to 11\% too small, with 95\% intervals covering 91 to 95\%. The
bootstrap of the subjects (\code{variance = "bootstrap"}; Efron and Tibshirani 1993) re-estimates
the weights in every sample and so includes their uncertainty, at the cost of \code{B} refits.
}
\section{Interpreting the coefficients}{
\itemize{
  \item The printed effect is the \strong{subdistribution hazard ratio} (sHR, \eqn{e^\beta}). Its
    direction is directly interpretable: \strong{sHR > 1 means a higher cumulative incidence of
    the cause at every time point}, sHR < 1 a lower one, for a one-unit increase of a numeric
    covariate or a factor level against the reference, other covariates fixed.
  \item Its \emph{size} is not a ratio of risks, and it is not the effect on the rate of the
    event among those still event-free (that is the cause-specific HR of \code{\link{csc}}). An
    sHR of 1.5 does not mean a 50\% higher risk. Report absolute risks next to it:
    \code{predict(fit, newdata, times)} gives \eqn{F_1(t \mid x)} (Austin and Fine 2017).
  \item A covariate can raise the cumulative incidence of cause 1 without any effect on its
    cause-specific hazard, only by lowering the hazard of the competing cause (subjects then live
    longer and have more time to experience cause 1). Fitting \code{\link{csc}} as well shows
    which mechanism is at work.
  \item Units, interactions and factors are read as in \code{\link{cox}}: the sHR for \eqn{k}
    units is \eqn{sHR^k}; with \code{a * b} the effect of \code{a} applies at \code{b = 0}.
}
A sentence for a report: \emph{"In a Fine-Gray model, treatment was associated with a lower
cumulative incidence of relapse (sHR 0.72, 95\% CI 0.58 to 0.90); the predicted 5-year
incidence was 24\% with and 31\% without treatment for a 60-year-old patient."}
}
\section{Assumptions and how to check them}{
\describe{
  \item{Proportional subdistribution hazards}{the sHR is constant over time. Compare
    \code{plot(fit, newdata)} with the non-parametric \code{\link{aj}} curves of the same groups,
    or fit the model on early and late follow-up separately (\code{\link{cif_plot}} with
    \code{type = "unadjusted"} overlays them).}
  \item{Censoring independent of the event and of the covariates}{the Kaplan-Meier weights assume
    censoring does not depend on the covariates. When it does strongly, results can be biased.}
  \item{Coherence across causes}{Fine-Gray models fitted separately for each cause need not give
    cumulative incidences that sum to at most 1. When predictions for all causes are needed
    together, \code{\link{directcif}} (\code{estimation = "joint"}) models them jointly.}
}
}
\section{Messages and common problems}{
\describe{
  \item{\code{WARNING: the fit did NOT converge}}{usually a covariate with no events of the cause
    in one of its levels (separation) or collinear covariates.}
  \item{\code{no events of the requested type}}{\code{cause} does not occur in the data; check the
    codes with \code{table(status)}.}
  \item{\code{strata() and cluster() are available in cox() and csc()}}{\code{fg()} has no
    stratified or clustered form; include the variable as a covariate.}
}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}; 82 simulated and real data sets, up to \eqn{n = 20000}): against
\code{cmprsk::crr} 2.2, coefficients to \eqn{10^{-7}}, the robust variance to \eqn{10^{-9}}
(including heavily tied data) and predicted cumulative incidences (\code{predict.crr}) to
\eqn{2 \times 10^{-8}}; with delayed entry, against \code{survival::finegray} + \code{coxph}, robust
variances to \eqn{2 \times 10^{-7}} and coefficients to \eqn{2 \times 10^{-5}}, where the remaining
difference comes from \code{coxph} stopping before the maximum (its log-likelihood is higher at the
\code{fg()} estimate). In a Monte Carlo study (1000 replications, \eqn{n} = 200, 500, 1000) the
robust standard errors without delayed entry matched the empirical standard deviation (ratio 0.96
to 1.02, coverage 94.0 to 95.6\%). With delayed entry they were 5 to 11\% too small (coverage 91
to 95\%), because the weights are treated as known, as in \pkg{survival}; with small samples and
heavy truncation, interpret the intervals with this in mind. On \eqn{n = 10000}, \code{fg()} takes
about 0.01 s and \code{crr} about 150 s.
}
\value{
An S7 object of class \code{FineGrayFit}. Main slots: \code{coefficients}, \code{vcov} (robust, bootstrap
or model, as requested), \code{vcov_robust}, \code{vcov_model}, \code{vcov_boot}, \code{loglik} (log
pseudo-likelihood), \code{loglik_null}, \code{cause}, \code{n}, \code{nevent}, \code{converged}.

Methods:
\tabular{ll}{
  \code{print(fit)}, \code{summary(fit)} \tab sHRs with intervals and p-values and the
    \code{Reading:} lines \cr
  \code{\link{estimates}(fit)} \tab the table as a data frame \cr
  \code{coef}, \code{vcov}, \code{confint}, \code{logLik}, \code{nobs} \tab extractors \cr
  \code{predict(fit, newdata, times)} \tab cumulative incidence of the cause, one row per row of
    \code{newdata}, one column per time \cr
  \code{plot(fit, ...)}, \code{\link{cif_plot}}, \code{\link{cif_data}} \tab cumulative
    incidence curves, with bootstrap intervals and adjusted (marginal) curves \cr
  \code{\link{brier_score}}, \code{\link{auc_t}}, \code{\link{compare_models}} \tab predictive
    performance \cr
}
}
\references{
\emph{Model, weights and robust variance:}
Fine JP, Gray RJ (1999). A proportional hazards model for the subdistribution of a competing
  risk. Journal of the American Statistical Association 94:496-509.

\emph{Sandwich variance in general:}
White H (1982). Maximum likelihood estimation of misspecified models. Econometrica 50:1-25.

\emph{Delayed entry:}
Geskus RB (2011). Cause-specific cumulative incidence estimation and the Fine and Gray model under
  both left truncation and right censoring. Biometrics 67:39-49.
Lin DY, Wei LJ (1989). The robust inference for the Cox proportional hazards model. Journal of
  the American Statistical Association 84:1074-1078.

\emph{Bootstrap variance:}
Efron B, Tibshirani RJ (1993). An Introduction to the Bootstrap. Chapman and Hall, New York.

\emph{Computation:}
Kawaguchi ES, Shen JI, Suchard MA, Li G (2021). Scalable algorithms for large competing risks
  data. Journal of Computational and Graphical Statistics 30:685-693.

\emph{Interpretation and reporting:}
Latouche A, Allignol A, Beyersmann J, Labopin M, Fine JP (2013). A competing risks analysis should
  report results on all cause-specific hazards and cumulative incidence functions. Journal of
  Clinical Epidemiology 66:648-653.
Austin PC, Fine JP (2017). Practical recommendations for reporting Fine-Gray model analyses for
  competing risk data. Statistics in Medicine 36:4391-4400.

The references of a particular fit: \code{sx_references(fit)}; all references behind
\code{fg()}: \code{sx_references("fg")}.
}
\seealso{
\code{\link{csc}} (cause-specific hazards), \code{\link{directcif}} (parametric cumulative
incidence), \code{\link{aj}} and \code{\link{gray_test}} (non-parametric), \code{\link{cif_plot}},
\code{\link{predict.SurvX}}, \code{\link{Event}}.
}
\examples{
set.seed(1)
n <- 400
d <- data.frame(treat = factor(rbinom(n, 1, 0.5), labels = c("control", "new")),
                age = round(rnorm(n, 60, 10)))
t1 <- rexp(n, 0.10 * exp(-0.4 * (d$treat == "new") + 0.02 * (d$age - 60)))   # relapse
t2 <- rexp(n, 0.05 * exp(0.04 * (d$age - 60)))                               # death without relapse
cens <- runif(n, 0, 12)
d$time <- pmin(t1, t2, cens)
d$status <- ifelse(d$time == cens, 0L, ifelse(t1 < t2, 1L, 2L))   # 0 censored, 1 relapse, 2 death
table(d$status)

## 1. subdistribution hazard ratios for relapse
fit <- fg(Event(time, status) ~ treat + I(age / 10), d, cause = 1)
fit

## 2. absolute risks: cumulative incidence of relapse at 2 and 5 years
nd <- data.frame(treat = c("control", "new"), age = 60)
predict(fit, nd, times = c(2, 5))
plot(fit, newdata = nd, conf.int = FALSE)

## 3. the mechanism: cause-specific hazards of both causes
csc(Event(time, status) ~ treat + I(age / 10), d)

## 4. delayed entry, with bootstrap standard errors
d$entry <- runif(n, 0, 0.3) * d$time
fg(Event(time, status, entry = entry) ~ treat + I(age / 10), d, variance = "bootstrap", B = 50)
}
