\name{estimates}
\alias{estimates}
\title{Estimates as a data frame}
\description{
All coefficients of a fitted model in one data frame, with standard errors, Wald confidence
limits and tests, and for regression coefficients the exponentiated estimate (hazard ratio,
subdistribution hazard ratio or odds ratio, see the column \code{scale}) with its limits.
}
\usage{
estimates(object, level = 0.95, ...)
}
\arguments{
\item{object}{a model fitted by \code{\link{cox}}, \code{\link{fg}}, \code{\link{csc}} or \code{\link{directcif}}.}
\item{level}{confidence level.}
\item{...}{not used.}
}
\value{
a data frame with the columns \code{cause, term, role, scale, estimate, se, lower, upper,
  z, p, exp_estimate, exp_lower, exp_upper}. \code{role} is \code{"beta"} for regression
  coefficients; for direct models also \code{"baseline"}, \code{"alpha"} and \code{"conditional"}
  (on the scale used for estimation; \code{summary()} shows them on the natural scale).
}
\examples{
set.seed(4)
d <- data.frame(time = rexp(200), status = rbinom(200, 1, 0.7), x = rnorm(200))
estimates(cox(Event(time, status) ~ x, d))
}
