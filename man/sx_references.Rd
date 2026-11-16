\name{sx_references}
\alias{sx_references}
\title{References of the methods behind a result}
\description{
Every result of SurvX rests on published methods. \code{sx_references()} prints the full references
of the methods that produced a fitted object, or of a function of the package, so that they can be
cited and checked. The printed output of every model ends with the short form of the same list
(the \dQuote{Methods:} line).
}
\usage{
sx_references(x = NULL)
}
\arguments{
\item{x}{a fitted object (from \code{\link{cox}}, \code{\link{fg}}, \code{\link{csc}},
  \code{\link{directcif}}, \code{\link{parsurv}}, \code{\link{cure}}, \code{\link{km}}, \code{\link{aj}},
  \code{\link{hazard_np}}, \code{\link{compare_models}} or \code{\link{ph_test}}), or the name of a
  function of the package as a character string (for example \code{"brier_score"}); \code{NULL}
  lists all references of the package.}
}
\details{
For a fitted object the list depends on the options that were used: for example the Efron or
Breslow method for ties of a Cox model, the estimation method and baselines of a direct model, or
delayed entry and the distribution of a parametric model.
}
\value{
invisibly, a character vector of full references.
}
\examples{
sx_references("ph_test")
set.seed(1)
d <- data.frame(time = rexp(50), status = rbinom(50, 1, 0.8), x = rnorm(50))
sx_references(cox(Event(time, status) ~ x, d))
}
