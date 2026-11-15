\name{rmst}
\alias{rmst}
\title{Restricted mean survival time}
\description{
Area under the Kaplan-Meier curve on [0, tau] with its standard error.
}
\usage{
rmst(time, status, tau)
}
\arguments{
\item{time,status}{follow-up time and 0/1 status.}
\item{tau}{upper limit of integration.}
}
\value{
a named vector with \code{rmst} and \code{se}.
}
\references{
Royston P, Parmar MKB (2013) BMC Med Res Methodol 13:152.
}
