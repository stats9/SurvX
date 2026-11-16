\name{fg}
\alias{fg}
\title{Fine-Gray model for the subdistribution hazard}
\description{
Proportional subdistribution hazards regression with inverse-probability-of-censoring
weights and the robust sandwich variance of Fine and Gray (1999).
}
\usage{
fg(formula, data, cause = 1L, variance = c("robust", "model"), max_iter = 50L,
    tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}; status 0 = censored, 1..K = cause.
  \code{Event(time, status, entry = )} for delayed entry (left truncation).}
\item{data}{a data frame.}
\item{cause}{cause of interest.}
\item{variance}{\code{"robust"} (default, Fine and Gray 1999, Section 4) or \code{"model"}
  (inverse information).}
\item{max_iter,tol}{iteration limit and convergence tolerance.}
}
\details{
Without delayed entry the model is fitted with the inverse-probability-of-censoring weights
\eqn{\hat G(t-)/\hat G(T_i-)} of Fine and Gray (1999), \eqn{\hat G} the Kaplan-Meier estimate of
the censoring distribution, and the default variance is their robust sandwich, which accounts for
the estimation of \eqn{\hat G}; results agree with \code{cmprsk::crr}.

With delayed entry the weights of Geskus (2011) are used: subjects with a competing event at
\eqn{T_i} stay in the risk set with weight \eqn{\hat G(t-)\hat H(t-)/\{\hat G(T_i-)\hat H(T_i-)\}},
where \eqn{\hat G} is now the product-limit estimate of the censoring distribution among subjects
under observation and \eqn{\hat H} the reverse-time product-limit estimate of the entry-time
distribution. The model is then a weighted Cox model on (start, stop] rows cut where the weights
change, fitted in the extended Cox core of \code{\link{cox}}; the robust variance is the sandwich
over subjects (Lin and Wei 1989) with the weights taken as known, as in
\code{survival::finegray} followed by \code{coxph(..., cluster = id)}, with which the results
agree.
}
\value{
an S7 object of class \code{FineGrayFit}. \code{predict(fit, newdata, times)} gives the
  cumulative incidence \eqn{1 - \exp\{-\Lambda_0(t) e^{z'\beta}\}}.
}
\references{
Fine JP, Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk.
  Journal of the American Statistical Association 94:496-509.

Geskus RB (2011) Cause-specific cumulative incidence estimation and the Fine and Gray model under
  both left truncation and right censoring. Biometrics 67:39-49.

Lin DY, Wei LJ (1989) The robust inference for the Cox proportional hazards model. Journal of the
  American Statistical Association 84:1074-1078.
}
