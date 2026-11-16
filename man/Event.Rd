\name{Event}
\alias{Event}
\title{Event-time response}
\description{
The response of every SurvX model: follow-up time and status, optionally with an entry time for
delayed entry (left truncation).
}
\usage{
Event(time, status, entry = NULL)
}
\arguments{
\item{time}{follow-up time (non-negative): the time of the event or of censoring.}
\item{status}{0 = censored; 1 = event, or 1, 2, ..., K = cause of failure with competing risks.
  A logical is read as 0/1. A factor is read with its first level as censoring.}
\item{entry}{optional entry time for delayed entry (left truncation): a subject is under
  observation from \code{entry} to \code{time} and contributes to the risk sets only in that
  interval. Must be smaller than \code{time}. Supported by \code{\link{parsurv}},
  \code{\link{csc}} with a parametric baseline and \code{\link{hazard_np}}; the other functions
  stop with a message.}
}
\value{
a matrix of class \code{"Event"} with the columns \code{time}, \code{status} and, with delayed
  entry, \code{entry}, for the left-hand side of a model formula.
}
\references{
Klein JP, Moeschberger ML (2003) Survival Analysis: Techniques for Censored and Truncated Data,
  2nd ed. Springer, Section 3.4 (left truncation).
}
\examples{
Event(c(5, 8, 3), c(1, 0, 2))
Event(c(5, 8, 3), c(1, 0, 1), entry = c(0, 2, 1))
}
