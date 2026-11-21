\name{rmst}
\alias{rmst}
\title{Restricted mean survival time}
\description{
The restricted mean survival time (RMST) up to a horizon \eqn{\tau}: the expected time lived
without the event during the first \eqn{\tau} units of follow-up, estimated by the area under the
Kaplan-Meier curve on \eqn{[0, \tau]}, with its standard error.

\strong{Use it when} an effect measure is wanted that does not need proportional hazards and is
easy to communicate ("on average, patients on the new treatment lived 2.1 months longer during
the first 5 years"), in particular when the survival curves cross.
\strong{Use something else when} there are competing events (the restricted mean time lost to a
cause is the area under its cumulative incidence, from \code{\link{aj}}), or covariates must be
adjusted for (fit \code{\link{cox}} or \code{\link{parsurv}} and average the predicted curves).
}
\usage{
rmst(time, status, tau)
}
\arguments{
\item{time}{follow-up time.}
\item{status}{0 = censored, 1 = event.}
\item{tau}{the horizon, in the unit of \code{time}. Choose it before looking at the results (a
  clinically meaningful time such as 5 years), and not beyond the follow-up of the data: beyond
  the last observed time the Kaplan-Meier curve is not estimated, and when the last time is
  censored the area is extended with the last value of the curve.}
}
\details{
\deqn{\widehat{RMST}(\tau) = \int_0^\tau \hat S(t)\,dt = \sum_j \hat S(t_{j-1})(t_j - t_{j-1}),}
summed over the event times up to \eqn{\tau} (Irwin 1949). The variance is the Greenwood-type
estimator \eqn{\sum_j A_j^2 d_j / \{Y_j (Y_j - d_j)\}} with \eqn{A_j = \int_{t_j}^\tau \hat S(t)dt}
(Klein and Moeschberger 2003, Section 4.5), the one of \code{survival::survfit(..., rmean = tau)}.
For two groups, the difference of the RMSTs has standard error
\eqn{\sqrt{se_1^2 + se_2^2}} (independent groups); see the example.
}
\section{Interpreting the result}{
An RMST of 3.8 years with \eqn{\tau} = 5 years means that, on average, a patient lived 3.8 of the
first 5 years without the event (and lost 1.2 years to it). The difference between two groups is
in time units, the same at every point of the curve, and needs no proportional hazards (Royston
and Parmar 2013). Always state \eqn{\tau} together with the result.

A sentence for a report: \emph{"Over 5 years, the restricted mean survival time was 3.8 years
(95\% CI 3.5 to 4.1) with the new treatment and 3.2 years (2.9 to 3.5) with the standard; difference
0.6 years (0.2 to 1.0)."}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}): the RMST agrees with \code{survival::survfit(..., rmean = tau)} to
\eqn{5 \times 10^{-12}} and with \code{survRM2::rmst2} to \eqn{2 \times 10^{-15}}, its standard error
with \code{survfit} to \eqn{2 \times 10^{-14}}, on 81 data sets. In a Monte Carlo study (1000
replications, \eqn{n} = 200 to 1000) the standard error matched the empirical standard deviation
(ratio 0.99 to 1.03) and the 95\% intervals covered 94.9 to 95.3\%.
}
\value{
A named numeric vector: \code{rmst} (the estimate) and \code{se} (its standard error).
}
\references{
\emph{Restricted mean:}
Irwin JO (1949). The standard error of an estimate of expectation of life, with special reference
  to expectation of tumourless life in experiments with mice. Journal of Hygiene 47:188-189.
Royston P, Parmar MKB (2013). Restricted mean survival time: an alternative to the hazard ratio
  for the design and analysis of randomized trials with a time-to-event outcome. BMC Medical
  Research Methodology 13:152.

\emph{Variance:}
Klein JP, Moeschberger ML (2003). Survival Analysis: Techniques for Censored and Truncated Data,
  2nd ed. Springer, Section 4.5.

All references: \code{sx_references("rmst")}.
}
\seealso{
\code{\link{km}} (the curve whose area this is), \code{\link{logrank_test}}.
}
\examples{
set.seed(1)
n <- 300
arm <- rep(c("standard", "new"), each = n / 2)
t <- rexp(n, ifelse(arm == "new", 0.15, 0.22)); cens <- runif(n, 2, 8)
time <- pmin(t, cens); status <- as.integer(t <= cens)

## RMST over the first 5 years in each arm
a <- rmst(time[arm == "new"], status[arm == "new"], tau = 5)
b <- rmst(time[arm == "standard"], status[arm == "standard"], tau = 5)
rbind(new = a, standard = b)

## difference with a 95\% confidence interval
diff <- a[["rmst"]] - b[["rmst"]]
se <- sqrt(a[["se"]]^2 + b[["se"]]^2)
c(difference = diff, lower = diff - 1.96 * se, upper = diff + 1.96 * se)
}
