\name{logrank_test}
\alias{logrank_test}
\alias{gray_test}
\title{K-sample tests}
\description{
\code{logrank_test} is the (stratified) log-rank test and, with \code{rho > 0}, the
G-rho family of Harrington and Fleming (1982). \code{gray_test} is Gray's (1988) test for
equality of the cumulative incidence of one cause across groups.
}
\usage{
logrank_test(time, status, group, strata = NULL, rho = 0)

gray_test(time, status, group, cause = 1L, strata = NULL, rho = 0)
}
\arguments{
\item{time,status}{follow-up time and status (0 = censored; for \code{gray_test} 1..K = cause).}
\item{group}{grouping variable.}
\item{strata}{optional stratification variable.}
\item{rho}{exponent of the weight function.}
\item{cause}{cause of interest.}
}
\value{
a list with \code{statistic}, \code{df}, \code{p.value} and, for the log-rank test,
  observed and expected numbers of events per group.
}
\references{
Mantel N (1966) Cancer Chemother Rep 50:163-170.
  Harrington DP, Fleming TR (1982) Biometrika 69:553-566.
  Gray RJ (1988) Ann Stat 16:1141-1154.
}
