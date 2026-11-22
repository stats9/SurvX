\name{cv_compare}
\alias{cv_compare}
\title{Cross-validated comparison of any set of models}
\description{
Fits every model to the full data (apparent performance) and in K-fold cross-validation, and
reports the Brier score, the index of prediction accuracy, the time-dependent AUC and the
calibration error of the predicted risk of one cause at the chosen times.
\code{\link{compare_models}} does this in the C++ core for the models of this package;
\code{cv_compare} is for anything else that has a \code{predict(fit, newdata, times, cause)} method
(a model from another package wrapped in a function, a model with its own options).
}
\usage{
cv_compare(formula, data, models, times, cause = 1L, folds = 5L, seed = 1L,
    groups = 10L)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}; passed to every model.}
\item{data}{a data frame; rows with missing values in the variables of \code{formula} are dropped.}
\item{models}{a named list of functions \code{function(formula, data)} that return a fitted model
  with a \code{predict(fit, newdata, times, cause)} method giving the cumulative incidence of
  the cause. Wrap any model, with its own options, in such
  a function. Fits of \code{cox()}, \code{parsurv()} and \code{cure()} for one event type are scored
  by \code{1 - predict(fit, newdata, times, type = "survival")}. With \code{+ frailty(id)} in
  \code{formula}, the folds keep every cluster whole and frailty models predict for a new cluster.}
\item{times}{time horizons.}
\item{cause}{cause of interest.}
\item{folds}{number of folds (less than 2: apparent performance only).}
\item{seed}{seed of the fold assignment: a whole number, or \code{NULL} to take the seed from the random numbers of the session (so that \code{set.seed()} before the call governs the result). All random draws come from R's generator (see \code{RNGkind()}); with a number the state of the session's generator is left as it was.}
\item{groups}{number of risk groups for the calibration error.}
}
\section{Interpreting the result}{
As for \code{\link{compare_models}}: judge the models on the cross-validated measures (lower Brier
score and calibration error, higher IPA and AUC). The same folds are used for every model, so the
comparison is paired. The table has the same columns as that of \code{compare_models} (standard
errors, paired differences from the best model, O/E and ICI; see the section "Columns of the
table" there).
}
\value{
an S7 object of class \code{ModelComparison}: \code{x@table} has one row per model and time
  (columns ending in \code{_cv} are computed from the pooled out-of-fold predictions), \code{x@fit}
  one row per model. A model that fails in a fold has \code{NA} cross-validated measures.
}
\references{
Gerds TA, Schumacher M (2007) Efron-type measures of prediction error for survival analysis.
  Biometrics 63:1283-1287.

Hastie T, Tibshirani R, Friedman J (2009) The Elements of Statistical Learning, 2nd ed. Springer,
  Section 7.10 (cross-validation).

Akaike H (1974) A new look at the statistical model identification. IEEE Transactions on Automatic
  Control 19:716-723.

All references: \code{sx_references("cv_compare")}.
}
\seealso{\code{\link{compare_models}}, \code{\link{brier_score}}.}
\examples{
set.seed(3)
n <- 300
x <- rbinom(n, 1, 0.5)
t1 <- rexp(n, 0.4 * exp(0.6 * x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
d <- data.frame(time = pmin(t1, t2, cens), x = x)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
cv_compare(Event(time, status) ~ x, d, times = c(1, 2), folds = 3, models = list(
  "direct PH, weighted" = function(f, d) directcif(f, d, estimation = "weighted"),
  "Fine-Gray" = function(f, d) fg(f, d),
  "cause-specific Cox" = function(f, d) csc(f, d)))
}
