\name{Event}
\alias{Event}
\title{Event-time response}
\description{
Event-time response
}
\usage{
Event(time, status)
}
\arguments{
\item{time}{follow-up time (non-negative).}
\item{status}{0 = censored; 1 = event, or 1, 2, ..., K = cause of failure with competing risks.
  A logical is read as 0/1. A factor is read with its first level as censoring.}
}
\value{
a two-column matrix of class \code{"Event"} for the left-hand side of a model formula.
}
