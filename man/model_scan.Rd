\name{model_scan}
\alias{model_scan}
\alias{report}
\alias{ModelScan}
\title{Fit and compare all suitable models for a data set}
\description{
Fits every model of the package that suits the data (one event type or competing risks) and
compares them on criteria that are valid for the comparison at hand: cross-validated prediction
for all models, the likelihood only within the models that share one, and checks of the
assumptions of every model. The result is a ranking with the models that predict equally well
grouped together, a reading in plain words, figures and, with \code{report()}, a complete HTML
report with the methods and references.

\strong{Use it when} you want an overview of how the models of the package fit your data before
choosing one, or want to show that the conclusion does not depend on the model.

\strong{It does not choose for you.} Models whose prediction error differs by less than its
standard error are equally good for prediction; choose among them on the question (rates or
risks), the assumptions and simplicity. Confidence intervals of a model chosen from the table
do not account for the choice.
}
\usage{
model_scan(formula, data, cause = 1L, goal = c("prediction", "effect"),
    times = NULL, horizon = NULL, folds = 5L, seed = 1L,
    depth = c("standard", "quick", "full"), models = NULL, verbose = FALSE)

report(x, file = "model_scan_report.html", ...)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}, with at least one covariate. Status
  0/1 for one event type; 0, 1, ..., K for competing risks. Delayed entry is not supported.}
\item{data}{a data frame. Rows with missing values in the variables of \code{formula} are
  removed.}
\item{cause}{competing risks: the cause whose cumulative incidence is predicted and compared.}
\item{goal}{\code{"prediction"} (default): the ranking is the answer. \code{"effect"}: in
  addition, \code{summary()} and the report show the effect of every covariate in every model,
  side by side, on the scale of each model (hazard ratio, subdistribution hazard ratio, odds
  ratio, time ratio).}
\item{times}{times for the AUC and the calibration error. Default: the quartiles of the event
  times of the cause.}
\item{horizon}{upper limit of the integrated Brier score. Default: the smaller of the 90th
  percentiles of the event times of the cause and of all observed times.}
\item{folds}{number of cross-validation folds (at least 2).}
\item{seed}{seed of the fold assignment; the state of the random number generator of the session
  is restored afterwards.}
\item{depth}{how many models: \code{"quick"} (the basic models), \code{"standard"} (default), or
  \code{"full"} (also the more flexible ones). Flexible models (generalized gamma, splines, GOR
  link) are tried only with at least 50 events of the cause.}
\item{models}{optional named list of additional models, each a function
  \code{function(formula, data)} returning a fit with a \code{predict(fit, newdata, times, ...)}
  method (the probability of the event or cumulative incidence of the cause, or for one event
  type \code{type = "survival"}).}
\item{verbose}{print the name of every model while it is fitted.}
\item{x}{a result of \code{model_scan()}.}
\item{file}{the HTML file to write.}
\item{...}{\code{report()}: \code{width} and \code{height} of the figures in inches.}
}
\details{
\subsection{The models}{
\describe{
  \item{One event type}{\code{\link{cox}}; \code{\link{parsurv}} with the exponential, Weibull,
    Gompertz (PH), log-normal and log-logistic (AFT) distributions, with \code{depth = "standard"}
    also the generalized gamma and the Royston-Parmar spline, with \code{"full"} the piecewise
    exponential; \code{\link{cure}} models (mixture; with \code{"full"} also non-mixture) when the
    Kaplan-Meier curve levels off. The screening rule for a plateau: the curve is flat over at
    least the last fifth of follow-up, with at least 10\% of the subjects event-free and at least
    10 subjects censored after the last event. It is not a test; look at \code{plot(km(...))}.}
  \item{Competing risks}{cause-specific Cox and parametric (Weibull, log-normal, spline) models
    (\code{\link{csc}}); \code{\link{fg}}; direct models (\code{\link{directcif}}): joint Gompertz
    with PH and PO links, weighted Gompertz with PH, PO and GOR links and a Weibull mixture
    baseline, the constrained model, and with \code{"full"} spline and logistic baselines.}
}
}
\subsection{The criteria}{
Three kinds of evidence are kept apart:
\enumerate{
  \item \strong{Prediction, for every model.} Every model is fitted in K-fold cross-validation and
    predicts the probability of the event (or cumulative incidence of the cause) of the left-out
    subjects. On the pooled out-of-fold predictions: the integrated Brier score (IBS) over
    \code{(0, horizon]} with inverse probability of censoring weights, the index of prediction
    accuracy (IPA, 1 minus the ratio to the model without covariates), the time-dependent AUC and
    the calibration error at \code{times}. These measures compare any models, whatever their
    likelihood.
  \item \strong{Likelihood, within comparable models only.} AIC is reported for the models with a
    full likelihood of the same data: the parametric and cure models for one event type; the
    parametric cause-specific, joint and constrained direct models for competing risks (their
    likelihoods are all the joint density of the time and the cause). The Cox and cause-specific
    Cox models have a partial likelihood, Fine-Gray and weighted direct models a weighted
    pseudo-likelihood: no AIC.
  \item \strong{Checks.} Convergence; the test of proportional hazards (\code{\link{ph_test}}) for the
    Cox models; for joint direct models, subjects whose predicted incidences add up to more than
    one; and for every model the largest distance between its average predicted risk and the
    non-parametric estimate (Kaplan-Meier or Aalen-Johansen) on the evaluation grid.
}
}
\subsection{Ranking}{
The models are ranked on the cross-validated IBS. The standard error of the difference from the
best model is the standard deviation of the fold-wise differences in IBS divided by the square
root of the number of folds (a paired comparison: the same folds for every model). Every model
within one standard error of the best gets rank 1: the one-standard-error rule of Hastie,
Tibshirani and Friedman (2009, Section 7.10). The other models are numbered after them in the
order of their IBS. The reading also says when no model predicts clearly better than the model
without covariates (largest IPA below 0.01), and when AIC clearly disfavours (dAIC > 10) a model
of rank 1.
}
}
\section{Interpreting the result}{
\describe{
  \item{\code{IBS}, \code{diff (SE)}}{prediction error, lower is better, and its difference from the
    best model with the standard error of that difference. A difference smaller than its SE is
    not evidence that one model predicts better.}
  \item{\code{IPA}}{share of the prediction error of the model without covariates that the model
    removes; 0 = no better than no covariates, negative = worse.}
  \item{\code{AUC}}{discrimination at \code{times} (averaged): the probability that a subject who
    has the event by \eqn{t} has a higher predicted risk than one who does not.}
  \item{\code{calib}}{calibration error at \code{times} (averaged): mean absolute difference
    between predicted and observed risk in ten risk groups. Small is good.}
  \item{\code{dAIC}}{difference from the smallest AIC among the comparable models. Differences
    below 2 are small, above 10 large (Burnham and Anderson 2002, Section 2.6).}
  \item{\code{checks}}{"ok", or the problems found. A model with a problem can still predict well;
    it should not be used to interpret effects without addressing the problem.}
  \item{Rates or risks}{with competing risks, the cause-specific models describe the rate of the
    event among those still event-free; Fine-Gray and the direct models describe the cumulative
    incidence. They answer different questions, and a covariate can have different effects on
    the two (Putter, Fiocco and Geskus 2007; Latouche et al. 2013). \code{goal = "effect"} shows
    both side by side.}
}
A sentence for a report: \emph{"We compared twelve candidate models in 5-fold cross-validation.
Four models, the cause-specific Cox model, the Fine-Gray model and two direct parametric models,
had integrated Brier scores within one standard error of the best (0.112); we report the
Fine-Gray model, whose subdistribution hazard ratios refer to the cumulative incidence."}
}
\section{Assumptions and how to check them}{
\describe{
  \item{Independent censoring}{as for every measure that uses censoring weights; the Kaplan-Meier
    estimate of the censoring distribution has no covariates.}
  \item{Enough events}{with fewer than about 10 events per fold and cause, the fold-wise measures
    are unstable and many flexible models fail; use \code{folds = 3} or \code{depth = "quick"}.}
  \item{The candidate set}{the ranking only compares the models that were fitted. Models that
    failed are listed under "Not fitted".}
}
}
\section{Messages and common problems}{
\describe{
  \item{\code{Not fitted: ...}}{the model failed on the full data or did not give finite
    predictions; with \code{verbose = TRUE} the models are named as they are fitted.}
  \item{\code{IBS} is \code{-}}{the model failed in at least one fold (\code{folds_failed} in
    \code{as.data.frame(x)}).}
  \item{Long computing times}{every model is fitted \code{folds + 1} times; the weighted GOR and
    spline models are the slowest. Use \code{depth = "quick"} first.}
}
}
\section{Validation}{
\code{model_scan()} calls the functions of the package, which are validated against reference
implementations (see their help pages and \code{vignette("validation", package = "SurvX")}),
and the measures of \code{\link{brier_score}}, \code{\link{auc_t}} and \code{\link{calibration}}.
The fold assignment, ranking and standard errors are tested in the test suite of the package.
}
\value{
An S7 object of class \code{ModelScan}. \code{x@table} (also \code{as.data.frame(x)}) has one
row per model: \code{model}, \code{family}, \code{likelihood} (\code{"full"}, \code{"partial"},
\code{"pseudo"}), \code{fitted}, \code{converged}, \code{loglik}, \code{npar}, \code{aic},
\code{delta_aic}, \code{ibs_cv}, \code{ipa_cv}, \code{auc_cv}, \code{calib_cv},
\code{np_distance}, \code{folds_failed}, \code{checks}, \code{rank}, \code{diff_best},
\code{se_diff}, \code{equivalent}. Also \code{x@per_time} (measures at each of \code{times}),
\code{x@brier} (cross-validated Brier score on the grid), \code{x@effects} (with
\code{goal = "effect"}) and \code{x@fits} (the models fitted to the full data, for further use).

Methods:
\tabular{ll}{
  \code{print(x)} \tab the ten best models and the reading \cr
  \code{summary(x)} \tab all models; with \code{goal = "effect"} the effects of the covariates \cr
  \code{plot(x, which)} \tab \code{"ranking"} (IBS with its SE), \code{"fit"} (average prediction
    against the non-parametric estimate), \code{"brier"} (prediction error over time),
    \code{"calibration"}; \code{top} models in the last three \cr
  \code{report(x, file)} \tab a self-contained HTML report: tables, figures, methods, references \cr
  \code{as.data.frame(x)} \tab the table \cr
  \code{\link{sx_references}(x)} \tab the references of the methods \cr
}
}
\references{
\emph{Prediction error and its estimation under censoring}

Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999) Assessment and comparison of prognostic
classification schemes for survival data. Stat Med 18:2529-2545.

Gerds TA, Schumacher M (2006) Consistent estimation of the expected Brier score in general
survival models with right-censored event times. Biom J 48:1029-1040.

Kattan MW, Gerds TA (2018) The index of prediction accuracy: an intuitive measure useful for
evaluating risk prediction models. Diagn Progn Res 2:7.

\emph{Discrimination and calibration}

Blanche P, Dartigues JF, Jacqmin-Gadda H (2013) Estimating and comparing time-dependent areas
under receiver operating characteristic curves for censored event times with competing risks.
Stat Med 32:5381-5397.

Gerds TA, Andersen PK, Kattan MW (2014) Calibration plots for risk prediction models in the
presence of competing risks. Stat Med 33:3191-3203.

\emph{Cross-validation and model choice}

Hastie T, Tibshirani R, Friedman J (2009) The Elements of Statistical Learning, 2nd ed.
Springer, Section 7.10.

Akaike H (1974) A new look at the statistical model identification. IEEE Trans Automat Contr
19:716-723.

Burnham KP, Anderson DR (2002) Model Selection and Multimodel Inference: A Practical
Information-Theoretic Approach, 2nd ed. Springer.

\emph{Rates and risks with competing risks}

Putter H, Fiocco M, Geskus RB (2007) Tutorial in biostatistics: competing risks and multi-state
models. Stat Med 26:2389-2430.

Latouche A, Allignol A, Beyersmann J, Labopin M, Fine JP (2013) A competing risks analysis should
report results on all cause-specific hazards and cumulative incidence functions. J Clin
Epidemiol 66:648-653.

All references: \code{sx_references("model_scan")}.
}
\seealso{
\code{\link{compare_models}} and \code{\link{cv_compare}} (a chosen set of models at chosen times),
\code{\link{brier_score}}, \code{\link{ph_test}}, \code{\link{time_effects}}, \code{\link{km}},
\code{\link{aj}}.
}
\examples{
## one event type
set.seed(1)
n <- 300
d <- data.frame(x = rbinom(n, 1, 0.5), z = rnorm(n))
t <- rweibull(n, 1.4, exp(1 - 0.5 * d$x + 0.3 * d$z)); cens <- runif(n, 0, 5)
d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
sc <- model_scan(Event(time, status) ~ x + z, d, depth = "quick", folds = 3)
sc
plot(sc, which = "ranking")

\donttest{
## competing risks, with the effects of the covariates in every model
data(fertility)
sc2 <- model_scan(Event(time, status) ~ education, fertility, cause = 1,
                  goal = "effect", depth = "quick", folds = 3)
summary(sc2)
plot(sc2)
report(sc2, file = tempfile(fileext = ".html"))
}
}
