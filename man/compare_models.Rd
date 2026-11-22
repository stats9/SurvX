\name{compare_models}
\alias{compare_models}
\title{Compare direct, Fine-Gray and cause-specific Cox models}
\description{
Fits a set of competing-risks models to the same data and compares how well they predict the
cumulative incidence of one cause: direct parametric models (\code{\link{directcif}}, any
estimator), the Fine-Gray model (\code{\link{fg}}), cause-specific Cox models
(\code{\link{csc}}) and, optionally, parametric cause-specific hazard models. For each model it
reports the apparent and K-fold cross-validated Brier score, index of prediction accuracy (IPA),
time-dependent AUC and calibration error at the chosen times, and the AIC where it is comparable.
Fitting, cross-validation and scoring are done in the C++ core, so the comparison is fast.

\strong{Use it when} choosing a model for predicting the risk of a cause.
\strong{Use} \code{\link{cv_compare}} for models that \code{compare_models} does not include
(other packages, other options).
}
\usage{
compare_models(formula, data, times, direct = NULL, cause = 1L, folds = 5L,
    seed = 1L, df = 3L, hazard = NULL, constrained = NULL)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}.}
\item{data}{a data frame.}
\item{times}{time horizons.}
\item{direct}{a list of direct models, each a character vector \code{c(baseline, link)},
  \code{c(baseline, link, estimation)} or \code{c(baseline, link, "constrained", cond)}; link
  \code{"PH"}, \code{"PO"} or \code{"GOR"}; estimation \code{"joint"} (default), \code{"weighted"} or
  \code{"constrained"}; see \code{\link{directcif}}. \code{NULL}: Gompertz baseline with PH and PO
  links, joint and weighted, and with two or more causes the constrained PH model with a Weibull
  \eqn{G}. Constrained models with three or more causes are fitted by \code{directcif()} with the
  same folds (as those of \code{constrained}).}
\item{cause}{cause of interest.}
\item{folds}{number of cross-validation folds (less than 2: apparent performance only).}
\item{seed}{seed of the fold assignment: a whole number, or \code{NULL} to take the seed from the random numbers of the session (so that \code{set.seed()} before the call governs the result). All random draws come from R's generator (see \code{RNGkind()}); with a number the state of the session's generator is left as it was.}
\item{df}{degrees of freedom of spline baselines.}
\item{hazard}{optional list of parametric cause-specific hazard models, fitted by \code{\link{csc}}
  with a parametric \code{baseline}: each a distribution name or \code{c("rp", df)}, e.g.
  \code{list("weibull", "lognormal", c("rp", 4))}. They are cross-validated with the same folds as
  the other models, and their AIC (full likelihood of the data) is comparable with that of the joint
  and constrained direct models.}
\item{constrained}{optional list of constrained direct models (\code{\link{directcif}} with
  \code{estimation = "constrained"}) with options that \code{direct} cannot express: each a list
  with any of \code{baseline} (default \code{"gompertz"}), \code{link} (\code{"PH"}), \code{cond}
  (\code{"weibull"}, one value or one per competing cause), \code{timing} and \code{shape} (the
  covariates of the time distributions of the competing causes), e.g.
  \code{list(list(), list(timing = ~ age))}. Any number of causes. They are fitted by
  \code{directcif()} and cross-validated with the same folds; their AIC is comparable with that of
  the joint and constrained direct models.}
}
\section{Interpreting the result}{
Judge the models on the cross-validated columns (ending in \code{_cv}): lower Brier score and
calibration error, higher IPA and AUC are better. Apparent measures, computed on the data used for
fitting, favour complex models. Differences in the third decimal of the Brier score are usually
within noise; prefer the simpler model when the measures are close. AIC can be compared only among
models with a full likelihood of the same data (joint and constrained direct models, parametric
cause-specific models); see Value.
}
\section{Columns of the table}{
Every measure is computed twice, from the predictions of the full-data fit (apparent) and from the
pooled out-of-fold predictions (columns ending in \code{_cv}): \code{brier}, \code{ipa}, \code{auc},
\code{calib_error} (see \code{\link{brier_score}}); \code{brier_se} and \code{auc_se}, standard
errors from influence functions with the censoring weights taken as known (as
\code{riskRegression::Score(conservative = TRUE)}), with the 95\% intervals \code{brier_lower},
\code{brier_upper}, \code{auc_lower}, \code{auc_upper}; \code{oe}, \code{oe_lower}, \code{oe_upper},
the observed (Aalen-Johansen) over the mean predicted risk with its interval; \code{ici},
\code{e50}, \code{e90}, the smooth calibration measures of Austin, Harrell and van Klaveren (2020);
\code{brier_diff}, \code{brier_diff_se}, \code{brier_diff_p} and \code{auc_diff},
\code{auc_diff_se}, \code{auc_diff_p}, the paired difference from the model with the lowest Brier
score at that time (named in \code{reference} and \code{reference_cv}), its standard error and the
p-value of the test of no difference. The paired tests use the same subjects for every model, so
they are much more sensitive than comparing two intervals. Cross-validated differences are an
honest guide; their p-values are approximate because the pooled predictions are not independent.

The concordance index (\code{\link{cindex}}) is deliberately not reported: it is not a proper
measure for risks at a time horizon (Blanche, Kattan and Gerds 2019), so it can prefer a wrong
model. The AUC at the horizon is its proper counterpart.
}
\value{
an S7 object of class \code{ModelComparison}; \code{x@table} has one row per model and
  time, \code{x@fit} one row per model. AIC is reported for the joint and constrained direct models
  (full likelihoods of the same data, comparable with each other); the weighted direct model,
  Fine-Gray and Cox maximise pseudo- or partial likelihoods.
}
\references{
\emph{Cross-validation:}
Stone M (1974). Cross-validatory choice and assessment of statistical predictions. Journal of the
  Royal Statistical Society, Series B 36:111-147.
Hastie T, Tibshirani R, Friedman J (2009). The Elements of Statistical Learning, 2nd ed.
  Springer, Section 7.10.
Gerds TA, Schumacher M (2007). Efron-type measures of prediction error for survival analysis.
  Biometrics 63:1283-1287.

\emph{Measures:} see \code{\link{brier_score}}.
Austin PC, Harrell FE, van Klaveren D (2020). Graphical calibration curves and the integrated
  calibration index (ICI) for survival models. Statistics in Medicine 39:2714-2742.
Blanche P, Kattan MW, Gerds TA (2019). The c-index is not proper for the evaluation of t-year
  predicted risks. Biostatistics 20:347-357.

All references: \code{sx_references("compare_models")}.
}
\seealso{
\code{\link{cv_compare}}, \code{\link{brier_score}}, \code{\link{directcif}}, \code{\link{fg}},
\code{\link{csc}}.
}
\examples{
set.seed(3)
n <- 300
d <- data.frame(x = rbinom(n, 1, 0.5))
t1 <- rexp(n, 0.4 * exp(0.6 * d$x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
d$time <- pmin(t1, t2, cens)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))

cm <- compare_models(Event(time, status) ~ x, d, times = c(1, 2), folds = 3,
                     direct = list(c("gompertz", "PH"), c("gompertz", "PO", "weighted")),
                     hazard = list("weibull"), constrained = list(list(timing = ~ x)))
cm
as.data.frame(cm)
}
