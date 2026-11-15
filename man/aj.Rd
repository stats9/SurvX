\name{aj}
\alias{aj}
\title{Aalen-Johansen estimate of the cumulative incidence functions}
\description{
Cumulative incidence of every cause with delta-method standard errors, one set of curves per
group, and Gray's test for every cause when there is more than one group.
}
\usage{
aj(formula, data, level = 0.95)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ 1} or \code{~ group}; status 0 = censored, 1..K = cause.}
\item{data}{a data frame.}
\item{level}{confidence level of the intervals.}
}
\value{
an S7 object of class \code{AJFit}; \code{x@curves} is a list of data frames with columns
  \code{time, surv, cif.k, se.k, lower.k, upper.k} for every cause k (intervals on the
  complementary log-log scale). Methods: \code{print}, \code{summary}, \code{plot}, \code{as.data.frame}.
}
\references{
Aalen OO, Johansen S (1978) Scand J Stat 5:141-150.
  Gray RJ (1988) Ann Stat 16:1141-1154.
}
