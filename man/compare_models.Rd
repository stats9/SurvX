\name{compare_models}
\alias{compare_models}
\title{Compare direct, Fine-Gray and cause-specific Cox models}
\description{
Fits a list of direct parametric models (with any of the three estimators), the Fine-Gray model
and cause-specific Cox models to the same data and reports apparent and K-fold cross-validated
Brier score, index of prediction accuracy, time-dependent AUC and calibration error at the chosen
times. Fitting, cross-validation and scoring are done in the C++ core.
}
\usage{
compare_models(formula, data, times, direct = NULL, cause = 1L, folds = 5L,
    seed = 1L, df = 3L)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}.}
\item{data}{a data frame.}
\item{times}{time horizons.}
\item{direct}{a list of direct models, each a character vector \code{c(baseline, link)},
  \code{c(baseline, link, estimation)} or \code{c(baseline, link, "constrained", cond)}; link
  \code{"PH"}, \code{"PO"} or \code{"GOR"}; estimation \code{"joint"} (default), \code{"weighted"} or
  \code{"constrained"}; see \code{\link{directcif}}. \code{NULL}: Gompertz baseline with PH and PO
  links, joint and weighted, and with two causes the constrained PH model with a Weibull \eqn{G}.}
\item{cause}{cause of interest.}
\item{folds}{number of cross-validation folds (less than 2: apparent performance only).}
\item{seed}{seed of the fold assignment.}
\item{df}{degrees of freedom of spline baselines.}
}
\value{
an S7 object of class \code{ModelComparison}; \code{x@table} has one row per model and
  time, \code{x@fit} one row per model. AIC is reported for the joint and constrained direct models
  (full likelihoods of the same data, comparable with each other); the weighted direct model,
  Fine-Gray and Cox maximise pseudo- or partial likelihoods.
}
\references{
Stone M (1974) J R Stat Soc B 36:111-147.
  Gerds TA, Schumacher M (2007) Biometrics 63:1283-1287.
}
