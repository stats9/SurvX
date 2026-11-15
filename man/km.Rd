\name{km}
\alias{km}
\title{Kaplan-Meier estimate}
\description{
Product-limit estimate of the survival function with Greenwood standard errors, one curve
per group, and the log-rank test when there is more than one group.
}
\usage{
km(formula, data, level = 0.95)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ 1} or \code{~ group}; several variables on the
  right-hand side are combined into one grouping factor.}
\item{data}{a data frame.}
\item{level}{confidence level of the intervals.}
}
\value{
an S7 object of class \code{KMFit}; \code{x@curves} is a list of data frames with columns
  \code{time, n.risk, n.event, surv, se, lower, upper} (intervals on the log scale), \code{x@counts}
  has the median survival time of every group with its confidence limits (Brookmeyer and Crowley 1982).
  Methods: \code{print}, \code{summary}, \code{plot}, \code{as.data.frame}.
}
\references{
Kaplan EL, Meier P (1958) J Am Stat Assoc 53:457-481.
  Greenwood M (1926) Reports on Public Health and Medical Subjects 33:1-26.
  Brookmeyer R, Crowley J (1982) Biometrics 38:29-41.
}
\examples{
set.seed(1)
d <- data.frame(time = rexp(100), status = rbinom(100, 1, 0.7), arm = gl(2, 50))
fit <- km(Event(time, status) ~ arm, d)
fit
}
