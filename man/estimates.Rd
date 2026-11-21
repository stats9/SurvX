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
\section{Interpreting the table}{
\code{estimate} and \code{se} are on the scale of estimation (log hazard ratio, log odds, log
shape, ...); \code{exp_estimate} with \code{exp_lower} and \code{exp_upper} is the effect on the
scale named in \code{scale} (hazard ratio, subdistribution hazard ratio, time ratio, odds ratio of
cure, ...), the one to report. Confidence limits and p-values are Wald limits and tests; for small
samples, compare nested models with a likelihood ratio test (\code{anova}-style: twice the
difference of \code{logLik()}).
}
\value{
a data frame with the columns \code{cause, term, role, scale, estimate, se, lower, upper,
  z, p, exp_estimate, exp_lower, exp_upper}. \code{role} is \code{"beta"} for regression
  coefficients; for direct models also \code{"baseline"}, \code{"alpha"} and \code{"conditional"}
  (on the scale used for estimation; \code{summary()} shows them on the natural scale).
}
\details{
Standard errors are the square roots of the diagonal of the inverse observed information (for
\code{fg()} the robust sandwich variance by default, for weighted direct models the sandwich with
one score per subject). Limits are Wald intervals \eqn{\hat\beta \pm z_{1-(1-level)/2}\,SE} on
the estimation scale, and \code{exp_*} columns are their exponentials; \code{p} is the two-sided
Wald test of a zero coefficient.
}
\references{
Wald A (1943) Tests of statistical hypotheses concerning several parameters when the number of
  observations is large. Transactions of the American Mathematical Society 54:426-482.

Casella G, Berger RL (2002) Statistical Inference, 2nd ed. Duxbury, Section 10.3.

White H (1982) Maximum likelihood estimation of misspecified models. Econometrica 50:1-25.
}
\seealso{
\code{\link{SurvX-methods}} (print, summary, coef, confint), \code{\link{sx_references}}.
}
\examples{
set.seed(4)
d <- data.frame(time = rexp(200), status = rbinom(200, 1, 0.7), x = rnorm(200))
estimates(cox(Event(time, status) ~ x, d))
}
