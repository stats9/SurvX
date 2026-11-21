\name{Interval}
\alias{Interval}
\alias{SurvXInterval}
\title{Interval-censored event times}
\description{
Builds the outcome of a model formula when the event time is not seen exactly but is known to lie
in an interval: the subject was event-free at the visit at time \code{left} and had had the event
by the visit at time \code{right}. Typical sources: periodic examinations (tumour recurrence or
progression seen at scheduled scans, HIV or HPV infection found at a test, dental caries, onset of
a chronic condition between two surveys). Use it on the left-hand side of
\code{\link{turnbull}} (non-parametric), \code{\link{parsurv}} (parametric regression) and
\code{\link{directcif}} (competing risks), for example
\code{parsurv(Interval(last_neg, first_pos) ~ age, data)}.
}
\usage{
Interval(left, right, cause = NULL)
}
\arguments{
\item{left}{the last time the subject was known to be event-free. \code{0} or \code{NA}: the
  event was already there at the first visit (left censoring, the event lies in
  \code{(0, right]}).}
\item{right}{the first time the event was known to have happened. \code{NA} or \code{Inf}: no
  event was seen (right censoring at \code{left}). \code{right = left}: the event time is known
  exactly.}
\item{cause}{competing risks only: the type of event (1, 2, ...) of the subjects with an event
  (finite \code{right}); ignored for right-censored subjects. A factor is read with its levels as
  causes 1, 2, ... (no censoring level: censoring is given by \code{right}).}
}
\details{
Every subject contributes the probability of the interval in which the event is known to lie:
\tabular{lll}{
  \code{left < right < Inf} \tab interval-censored \tab \eqn{S(L) - S(R)} \cr
  \code{left = 0, right < Inf} \tab left-censored \tab \eqn{1 - S(R)} \cr
  \code{right = Inf} \tab right-censored \tab \eqn{S(L)} \cr
  \code{left = right} \tab exact \tab density \eqn{f(L)} \cr
}
(with competing risks \eqn{F_k(R) - F_k(L)}, \eqn{1 - \sum_k F_k(L)} and \eqn{f_k(L)}). Imputing
the midpoint or the right end point of the interval and then treating the times as exact biases the
estimates and understates their uncertainty (Lindsey and Ryan 1998; Law and Brookmeyer 1992). The
inspection process must be non-informative: the visit times may differ between subjects but must
not depend on the event itself (given the covariates). Rows with \code{left = 0} and
\code{right = Inf} carry no information and are dropped by the model functions.
}
\section{Messages}{
\describe{
  \item{\code{right must be at least left}}{an interval with its end before its start.}
  \item{\code{cause must be 1, 2, ...}}{a subject with an event and no (or a zero) cause.}
}
}
\value{
A numeric matrix of class \code{"SurvXInterval"} with the columns \code{left}, \code{right}
(\code{Inf} for right censoring) and \code{status} (0 right-censored, otherwise the cause). It
prints as \code{(2, 5]}, \code{3} (exact), \code{4+} (right-censored) and \code{(0, 5]}
(left-censored).
}
\references{
Turnbull BW (1976). The empirical distribution function with arbitrarily grouped, censored and
  truncated data. Journal of the Royal Statistical Society B 38:290-295.

Lindsey JC, Ryan LM (1998). Tutorial in biostatistics: methods for interval-censored data.
  Statistics in Medicine 17:219-238.

Law CG, Brookmeyer R (1992). Effects of mid-point imputation on the analysis of doubly censored
  data. Statistics in Medicine 11:1569-1578.

Sun J (2006). The Statistical Analysis of Interval-censored Failure Time Data. Springer.
}
\seealso{
\code{\link{turnbull}}, \code{\link{parsurv}}, \code{\link{directcif}}, \code{\link{Event}}.
}
\examples{
## event between visits, exact, right-censored, left-censored
Interval(left = c(2, 3, 4, 0), right = c(5, 3, NA, 6))

## competing risks: the type of event is known
Interval(c(2, 1, 4), c(5, 3, Inf), cause = c(1, 2, 0))
}
