\name{brier_score}
\alias{brier_score}
\alias{auc_t}
\alias{calibration}
\title{Prediction performance at a time horizon}
\description{
Measures how good predicted risks are at a time horizon \eqn{t}, for one event type or one cause
with competing risks, under right censoring:
\itemize{
  \item \code{brier_score()}: the Brier score (mean squared error of the predicted risk) and the
    index of prediction accuracy (IPA), which compares it with a model without covariates;
  \item \code{auc_t()}: the time-dependent area under the ROC curve (discrimination: are
    subjects who have the event by \eqn{t} given higher risks?);
  \item \code{calibration()}: predicted against observed risk in risk groups (do predicted
    risks of 20\% come true for 20\% of subjects?).
}
\strong{Use them} to validate a prediction model, ideally on new data or with cross-validation
(\code{\link{compare_models}}, \code{\link{cv_compare}}), which call these functions.
}
\usage{
brier_score(time, status, risk, t, cause = 1L)

auc_t(time, status, risk, t, cause = 1L)

calibration(time, status, risk, t, cause = 1L, groups = 10L)
}
\arguments{
\item{time,status}{follow-up time and status (0 = censored, 1..K = cause) of the subjects whose
  risks are assessed.}
\item{risk}{predicted probability of the event (or of the cause of interest) by time \code{t}, one
  per subject: for example \code{predict(fit, newdata, times = t)} of a Fine-Gray or
  cause-specific model, or \code{1 - predict(fit, newdata, times = t, type = "survival")} of a
  survival model.}
\item{t}{the time horizon (a single number).}
\item{cause}{cause of interest.}
\item{groups}{number of risk groups for calibration.}
}
\details{
\subsection{Brier score}{
\eqn{BS(t) = n^{-1} \sum_i \hat W_i(t) \{I(T_i \le t, \mathrm{cause}_i = 1) - \hat r_i\}^2},
with inverse-probability-of-censoring weights \eqn{\hat W_i(t)} from the Kaplan-Meier estimate of
the censoring distribution: \eqn{1/\hat G(T_i-)} for subjects with an event (of any cause) by
\eqn{t}, \eqn{1/\hat G(t)} for those still at risk, and 0 for those censored before \eqn{t}
(Graf et al. 1999; Gerds and Schumacher 2006). The null model predicts the Aalen-Johansen
(Kaplan-Meier) estimate for everybody, and \eqn{IPA = 1 - BS/BS_{null}} (Kattan and Gerds 2018).
}
\subsection{AUC}{
The probability that a subject with the event of interest by \eqn{t} has a higher predicted risk
than a subject without it, with the same weights; subjects with a competing event before \eqn{t}
count as controls ("definition 2" of Blanche et al. 2013). Ties in the predicted risk count one
half.
}
\subsection{Calibration}{
Subjects are grouped by predicted risk (\code{groups} quantile groups); in each group the mean
predicted risk is compared with the Aalen-Johansen estimate of the observed risk at \eqn{t}
(Gerds, Andersen and Kattan 2014). The calibration error is the size-weighted mean absolute
difference.
}
}
\section{Interpreting the result}{
\itemize{
  \item Brier score: 0 is perfect; a model that predicts 50\% for everybody scores 0.25. It
    depends on how common the event is, so compare it with the null model: IPA > 0 means the
    model is better than ignoring the covariates, IPA = 0.10 means 10\% of the null model's
    error is removed. Negative IPA: the model is worse than no model.
  \item AUC: 0.5 is no discrimination, 1 perfect; 0.7 to 0.8 is common for clinical models. AUC
    ignores calibration: a model can rank well and still predict wrong risks.
  \item Calibration: points near the diagonal (mean predicted = observed) and a small error.
}
Measured on the data used for fitting, all three are optimistic; use cross-validation or
external data for an honest assessment.
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}), the Brier score, that of the null model and the AUC agree with
\code{riskRegression::Score} to \eqn{2 \times 10^{-14}} on 82 data sets with and without competing
risks. They are about 10 times faster than \code{Score} on \eqn{n = 10000}.
}
\value{
\code{brier_score}: named vector with \code{brier}, \code{brier_null} (Aalen-Johansen
  prediction for everyone) and \code{ipa = 1 - brier / brier_null}. \code{auc_t}: the AUC.
  \code{calibration}: a data frame with mean predicted and observed (Aalen-Johansen) risk per
  group; its attribute \code{"error"} is the size-weighted mean absolute difference.
}
\references{
Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999) Assessment and comparison of prognostic
  classification schemes for survival data. Statistics in Medicine 18:2529-2545.

Gerds TA, Schumacher M (2006) Consistent estimation of the expected Brier score in general survival
  models with right-censored event times. Biometrical Journal 48:1029-1040.

Kattan MW, Gerds TA (2018) The index of prediction accuracy: an intuitive measure useful for
  evaluating risk prediction models. Diagnostic and Prognostic Research 2:7.

Blanche P, Dartigues J-F, Jacqmin-Gadda H (2013) Estimating and comparing time-dependent areas
  under receiver operating characteristic curves for censored event times with competing risks.
  Statistics in Medicine 32:5381-5397.

Gerds TA, Andersen PK, Kattan MW (2014) Calibration plots for risk prediction models in the
  presence of competing risks. Statistics in Medicine 33:3191-3203.

All references: \code{sx_references("brier_score")}, \code{sx_references("auc_t")},
\code{sx_references("calibration")}.
}
\seealso{
\code{\link{compare_models}} and \code{\link{cv_compare}} (cross-validated comparison),
\code{\link{predict.SurvX}}.
}
\examples{
set.seed(3)
n <- 400
d <- data.frame(x = rnorm(n))
t1 <- rexp(n, 0.3 * exp(0.8 * d$x)); t2 <- rexp(n, 0.2); cens <- runif(n, 0, 5)
d$time <- pmin(t1, t2, cens)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))

## predicted 2-year risk of cause 1 from a Fine-Gray model
fit <- fg(Event(time, status) ~ x, d)
risk <- predict(fit, d, times = 2)[, 1]

brier_score(d$time, d$status, risk, t = 2)     # brier, brier_null, ipa
auc_t(d$time, d$status, risk, t = 2)
cal <- calibration(d$time, d$status, risk, t = 2, groups = 5)
cal
plot(cal$predicted, cal$observed, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "predicted risk", ylab = "observed risk"); abline(0, 1, lty = 2)
}
