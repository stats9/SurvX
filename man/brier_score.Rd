\name{brier_score}
\alias{brier_score}
\alias{auc_t}
\alias{calibration}
\title{Prediction performance at a time horizon}
\description{
Inverse-probability-of-censoring weighted Brier score with the index of prediction accuracy,
time-dependent AUC, and grouped calibration, for predicted risks of one cause by time \code{t}.
}
\usage{
brier_score(time, status, risk, t, cause = 1L)

auc_t(time, status, risk, t, cause = 1L)

calibration(time, status, risk, t, cause = 1L, groups = 10L)
}
\arguments{
\item{time,status}{follow-up time and status (0 = censored, 1..K = cause).}
\item{risk}{predicted probability of the cause of interest by time \code{t}, one per subject.}
\item{t}{time horizon.}
\item{cause}{cause of interest.}
\item{groups}{number of risk groups for calibration.}
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
}
