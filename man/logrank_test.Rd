\name{logrank_test}
\alias{logrank_test}
\alias{gray_test}
\title{Comparing groups: log-rank test and Gray's test}
\description{
Tests whether two or more groups differ, without any model for the covariates.
\itemize{
  \item \code{logrank_test()}: equality of the \emph{survival curves} (one event type), the
    log-rank test (Mantel 1966; Peto and Peto 1972) and the \eqn{G^\rho} family of Harrington and
    Fleming (1982), optionally stratified.
  \item \code{gray_test()}: equality of the \emph{cumulative incidence} of one cause when there are
    competing events (Gray 1988).
}
\code{\link{km}} and \code{\link{aj}} run these tests automatically when the formula has groups;
call these functions directly for the \eqn{G^\rho} weights, for a stratified test, or when the
data are in vectors.

\strong{Use something else when} the effect has to be adjusted for other covariates
(\code{\link{cox}}, \code{\link{csc}}, \code{\link{fg}}), or when you need an effect size
(a hazard ratio from \code{cox}, or a difference in restricted mean survival time from
\code{\link{rmst}}).
}
\usage{
logrank_test(time, status, group, strata = NULL, rho = 0, timefix = TRUE)
gray_test(time, status, group, cause = 1L, strata = NULL, rho = 0, timefix = TRUE)
}
\arguments{
\item{time}{follow-up time (numeric vector).}
\item{status}{\code{logrank_test}: 0 = censored, any positive value = event.
  \code{gray_test}: 0 = censored, 1, 2, ... = cause.}
\item{group}{the groups to compare (any vector; it is converted to a factor). At least two
  groups.}
\item{strata}{optional stratification variable: the comparison is made within each stratum and
  the results are combined, so that a difference between strata (centres, for example) is not
  mistaken for a difference between groups.}
\item{rho}{\code{logrank_test}: the \eqn{G^\rho} weight \eqn{\hat S(t-)^\rho}. \code{0} (default)
  is the log-rank test, which has the most power when the hazards are proportional; \code{1} is
  the Peto-Peto version of the Wilcoxon test, which gives more weight to early differences.
  \code{gray_test}: the exponent of Gray's weight function (0 by default, as in
  \code{cmprsk::cuminc}).}
\item{cause}{\code{gray_test}: the cause whose cumulative incidence is compared.}
\item{timefix}{treat times that differ only by rounding as tied (see \code{\link{Event}}).}
}
\details{
\subsection{Log-rank and \eqn{G^\rho} tests}{
At every event time \eqn{t}, the observed number of events in each group is compared with the
number expected if all groups had the same hazard, given the numbers at risk. With weights
\eqn{w(t) = \hat S(t-)^\rho} (\eqn{\hat S} the pooled Kaplan-Meier estimate), the statistic is
\eqn{U^\top V^{-} U}, where \eqn{U} sums \eqn{w(t)(O - E)} over the event times and \eqn{V} is
its hypergeometric variance; under the null hypothesis it has a chi-square distribution with
(number of groups - 1) degrees of freedom. Stratified: \eqn{U} and \eqn{V} are summed over strata.
}
\subsection{Gray's test}{
Compares weighted averages of the subdistribution hazards of the cause across groups (Gray 1988,
with the variance of that paper), as in \code{cmprsk::cuminc}. It tests the cumulative
incidence, not the cause-specific hazard: a treatment that reduces only competing deaths can
increase the cumulative incidence of the cause, and Gray's test then detects that increase.
}
}
\section{Interpreting the result}{
A small p-value means that the curves (survival, or cumulative incidence of the cause) are not
all equal. The test says nothing about the size of the difference: report the curves
(\code{plot(km(...))}, \code{plot(aj(...))}), the median survival per group, or an effect
estimate. \code{observed} and \code{expected} of the log-rank test show which group has more
events than expected. The log-rank test can miss curves that cross (the early and late
differences cancel); then look at the curves and consider \code{rho = 1} or \code{\link{rmst}}.

A sentence for a report: \emph{"Overall survival differed between the arms (log-rank test,
chi-square 7.9 on 1 df, p = 0.005)."}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}): the log-rank, Peto-Peto and stratified statistics agree with
\code{survival::survdiff} to a relative \eqn{3 \times 10^{-6}}, and Gray's statistics with
\code{cmprsk::cuminc} to \eqn{10^{-11}}, on 81 data sets. In a Monte Carlo study (1000
replications) the size of the log-rank test was 4.1 to 5.7\% and that of Gray's test 3.6 to 4.9\%
at the 5\% level.
}
\value{
A list with \code{method}, \code{statistic} (chi-square), \code{df}, \code{p.value}; for
\code{logrank_test} also \code{observed} and \code{expected} (number of events per group), for
\code{gray_test} also \code{cause}.
}
\references{
\emph{Log-rank test:}
Mantel N (1966). Evaluation of survival data and two new rank order statistics arising in its
  consideration. Cancer Chemotherapy Reports 50:163-170.
Peto R, Peto J (1972). Asymptotically efficient rank invariant test procedures. Journal of the
  Royal Statistical Society, Series A 135:185-207.

\emph{\eqn{G^\rho} family:}
Harrington DP, Fleming TR (1982). A class of rank test procedures for censored survival data.
  Biometrika 69:553-566.

\emph{Gray's test:}
Gray RJ (1988). A class of K-sample tests for comparing the cumulative incidence of a competing
  risk. Annals of Statistics 16:1141-1154.

All references: \code{sx_references("logrank_test")}, \code{sx_references("gray_test")}.
}
\seealso{
\code{\link{km}} and \code{\link{aj}} (curves with the tests), \code{\link{rmst}},
\code{\link{cox}}, \code{\link{fg}}.
}
\examples{
set.seed(1)
n <- 200
arm <- rep(c("A", "B"), each = n / 2)
centre <- sample(c("c1", "c2"), n, TRUE)
t <- rexp(n, ifelse(arm == "B", 0.7, 1)); cens <- runif(n, 0, 3)
time <- pmin(t, cens); status <- as.integer(t <= cens)

logrank_test(time, status, arm)                    # log-rank
logrank_test(time, status, arm, rho = 1)           # Peto-Peto (early differences)
logrank_test(time, status, arm, strata = centre)   # stratified by centre

## competing risks: does the cumulative incidence of cause 1 differ?
cause <- ifelse(status == 1, sample(1:2, n, TRUE, prob = c(0.6, 0.4)), 0)
gray_test(time, cause, arm, cause = 1)
}
