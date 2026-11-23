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
    risks of 20\% come true for 20\% of subjects?), the observed/expected ratio and the
    integrated calibration index of a smooth calibration curve.
}
With \code{se = TRUE} the Brier score and the AUC come with standard errors. Related:
\code{\link{cindex}} (concordance over the whole follow-up) and \code{\link{decision_curve}}
(net benefit).
\strong{Use them} to validate a prediction model, ideally on new data or with cross-validation
(\code{\link{compare_models}}, \code{\link{cv_compare}}), which call these functions.
}
\usage{
brier_score(time, status, risk, t, cause = 1L, se = FALSE)

auc_t(time, status, risk, t, cause = 1L, se = FALSE, level = 0.95)

calibration(time, status, risk, t, cause = 1L, groups = 10L, level = 0.95)
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
\item{se}{also return the standard error (see Details).}
\item{level}{confidence level of the interval of the AUC and of the observed/expected ratio.}
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
\subsection{Standard errors}{
With \code{se = TRUE} the standard error of the Brier score and of the AUC is
\eqn{\mathrm{sd}(\psi_i)/\sqrt{n}}, where \eqn{\psi_i} is the influence function of subject
\eqn{i} with the censoring weights taken as known (Blanche et al. 2013): the "conservative"
standard errors of \code{riskRegression::Score(conservative = TRUE)}. \code{\link{compare_models}}
and \code{\link{cv_compare}} use the same influence functions for the paired difference of two
models.
}
\subsection{Calibration in the large and smooth calibration}{
\code{calibration()} also returns the observed/expected ratio: the Aalen-Johansen (Kaplan-Meier)
risk at \eqn{t} over the mean predicted risk, with the interval
\eqn{\exp\{\log(O/E) \pm z\, se(O)/O\}}; and the integrated calibration index of Austin, Harrell
and van Klaveren (2020): a Cox model (with competing risks, a Fine-Gray model) of the outcome on a
natural cubic spline of \eqn{\log\{-\log(1 - \hat r)\}} with knots at its 10th, 50th and 90th
percentiles gives a calibrated risk for every subject; ICI, E50 and E90 are the mean, median and
90th percentile of its absolute difference from the predicted risk.
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
    O/E = 1 means the average risk is right (O/E = 1.2: the model predicts 20\% too few events);
    an ICI of 0.02 means the predicted risks are on average 2 percentage points from the calibrated
    risks.
}
Measured on the data used for fitting, all three are optimistic; use cross-validation or
external data for an honest assessment.
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}), the Brier score, that of the null model and the AUC agree with
\code{riskRegression::Score} to \eqn{2 \times 10^{-14}} on 82 data sets with and without competing
risks. They are about 10 times faster than \code{Score} on \eqn{n = 10000}. The standard errors
equal those of \code{Score(conservative = TRUE, se.fit = TRUE)}; the O/E equals the Aalen-Johansen
estimate of \code{survival::survfit}; ICI, E50 and E90 equal a calibration model fitted with
\code{survival::coxph} and \code{splines::ns} (validation suite, \code{validation/R/07_assess.R}).
}
\value{
\code{brier_score}: named vector with \code{brier}, \code{brier_null} (Aalen-Johansen
  prediction for everyone) and \code{ipa = 1 - brier / brier_null}, and with \code{se = TRUE}
  \code{se}. \code{auc_t}: the AUC; with \code{se = TRUE} a vector with \code{auc}, \code{se},
  \code{lower} and \code{upper}.
  \code{calibration}: a data frame with mean predicted and observed (Aalen-Johansen) risk per
  group; its attribute \code{"error"} is the size-weighted mean absolute difference, \code{"oe"}
  the observed/expected ratio with its interval (three numbers), \code{"ici"}, \code{"e50"} and
  \code{"e90"} the smooth calibration measures and \code{"smooth"} a data frame of the predicted
  and calibrated risks, ordered by predicted risk (\code{NULL} when the calibration model cannot be
  fitted, for example when all predicted risks are equal).
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

Austin PC, Harrell FE, van Klaveren D (2020) Graphical calibration curves and the integrated
  calibration index (ICI) for survival models. Statistics in Medicine 39:2714-2742.

All references: \code{sx_references("brier_score")}, \code{sx_references("auc_t")},
\code{sx_references("calibration")}.
}
\seealso{
\code{\link{compare_models}} and \code{\link{cv_compare}} (cross-validated comparison),
\code{\link{cindex}}, \code{\link{decision_curve}}, \code{\link{predict.SurvX}},
\code{vignette("prediction", "SurvX")}.
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
brier_score(d$time, d$status, risk, t = 2, se = TRUE)
auc_t(d$time, d$status, risk, t = 2, se = TRUE)   # auc, se, lower, upper
cal <- calibration(d$time, d$status, risk, t = 2, groups = 5)
cal
attr(cal, "oe")                                # observed / expected, with its interval
c(ICI = attr(cal, "ici"), E50 = attr(cal, "e50"), E90 = attr(cal, "e90"))
plot(cal$predicted, cal$observed, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "predicted risk", ylab = "observed risk"); abline(0, 1, lty = 2)
}
