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
Graf E et al. (1999) Stat Med 18:2529-2545.
  Gerds TA, Schumacher M (2006) Biom J 48:1029-1040.
  Kattan MW, Gerds TA (2018) Diagn Progn Res 2:7.
  Blanche P, Dartigues J-F, Jacqmin-Gadda H (2013) Stat Med 32:5381-5397.
}
