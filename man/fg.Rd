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
\item{formula}{\code{Event(time, status) ~ covariates}; status 0 = censored, 1..K = cause.}
\item{data}{a data frame.}
\item{cause}{cause of interest.}
\item{variance}{\code{"robust"} (default, Fine and Gray 1999, Section 4) or \code{"model"}
  (inverse information).}
\item{max_iter,tol}{iteration limit and convergence tolerance.}
}
\value{
an S7 object of class \code{FineGrayFit}. \code{predict(fit, newdata, times)} gives the
  cumulative incidence \eqn{1 - \exp\{-\Lambda_0(t) e^{z'\beta}\}}.
}
\references{
Fine JP, Gray RJ (1999) J Am Stat Assoc 94:496-509.
}
