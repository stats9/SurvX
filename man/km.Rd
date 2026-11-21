\name{km}
\alias{km}
\title{Kaplan-Meier estimate of the survival function}
\description{
Estimates the probability of being free of the event over time, \eqn{S(t) = P(T > t)}, without
any model: one curve per group, with confidence intervals, the median survival time of every
group, and the log-rank test when there are groups.

\strong{Use it when} there is one type of event and you want to describe survival or compare
groups. It is usually the first step of an analysis.
\strong{Use something else when}
\itemize{
  \item there are competing events: \code{\link{aj}}. One minus the Kaplan-Meier estimate that
    treats competing events as censored is \emph{not} the probability of the event and
    overestimates it;
  \item there is delayed entry (left truncation): \code{km()} does not handle it;
    \code{\link{hazard_np}} estimates the cumulative hazard with delayed entry;
  \item covariates must be adjusted for: \code{\link{cox}} or \code{\link{parsurv}}.
}
}
\usage{
km(formula, data, level = 0.95)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ 1} for one curve, or \code{Event(time, status) ~ group}
  for one curve per group; several variables on the right (\code{~ sex + stage}) give one curve
  per combination. \code{status} must be 0/1.}
\item{data}{a data frame. Rows with a missing time, status or group are removed.}
\item{level}{confidence level of the intervals (default 0.95).}
}
\details{
\subsection{Estimator}{
\deqn{\hat S(t) = \prod_{t_j \le t} \left(1 - \frac{d_j}{Y_j}\right),}
with \eqn{d_j} the number of events and \eqn{Y_j} the number at risk at event time \eqn{t_j}
(Kaplan and Meier 1958). Its variance is Greenwood's (1926) formula; the confidence interval is
computed on the log scale, \eqn{\hat S \exp(\pm z\,se\{\log \hat S\})}, truncated at 1, the
default of \code{survival::survfit}.
}
\subsection{Median survival}{
The median is the first time at which \eqn{\hat S(t) \le 0.5}; its confidence limits are the first
times at which the lower and upper limits of the curve reach 0.5 (Brookmeyer and Crowley 1982).
When the curve does not fall to 0.5 during follow-up, the median is \code{NA}: report the
survival at fixed times instead.
}
}
\section{Interpreting the result}{
\itemize{
  \item \eqn{\hat S(t)} is the estimated proportion still free of the event at time \eqn{t}; it is
    valid when censoring is unrelated to the risk of the event.
  \item The curve is reliable only where enough subjects are at risk: the intervals widen at the
    end of follow-up. The risk table under \code{plot(fit)} shows the numbers at risk.
  \item The log-rank test in the print compares the whole curves; see \code{\link{logrank_test}}
    for when it can miss differences.
}
A sentence for a report: \emph{"Five-year survival was 62\% (95\% CI 55 to 70\%) in arm A and 48\%
(40 to 57\%) in arm B; median survival 6.8 and 4.9 years (log-rank p = 0.01)."}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}), on 81 data sets: survival at the event times and Greenwood standard errors
agree with \code{survival::survfit} to \eqn{4 \times 10^{-16}}, medians and their confidence limits
exactly. In a Monte Carlo study (1000 replications, \eqn{n} = 200 to 1000) the 95\% intervals
covered 93.7 to 96.1\%.
}
\value{
An S7 object of class \code{KMFit}. Slots: \code{curves} (a list with one data frame per group,
columns \code{time}, \code{n.risk}, \code{n.event}, \code{surv}, \code{se} (of \eqn{\log S}),
\code{lower}, \code{upper}), \code{counts} (per group: \code{n}, \code{events}, \code{median} with
\code{median.lower} and \code{median.upper}), \code{test} (the log-rank test) and \code{level}.

Methods:
\tabular{ll}{
  \code{print(fit)} \tab numbers, events, medians with intervals, log-rank test \cr
  \code{summary(fit, times = c(1, 3, 5))} \tab survival with intervals at chosen times \cr
  \code{plot(fit, ...)} \tab survival curves with a risk table; arguments \code{conf.int},
    \code{risk.table}, \code{median}, \code{mark.times}, \code{percent}, \code{labels}, \code{col},
    \code{xlim} and the usual graphical ones (see \code{\link{SurvX-methods}}) \cr
  \code{as.data.frame(fit)} \tab all curves in one data frame with a \code{group} column \cr
}
}
\references{
\emph{Estimator and variance:}
Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete observations. Journal of the
  American Statistical Association 53:457-481.
Greenwood M (1926). The natural duration of cancer. Reports on Public Health and Medical Subjects
  33:1-26.

\emph{Median and its confidence interval:}
Brookmeyer R, Crowley J (1982). A confidence interval for the median survival time. Biometrics
  38:29-41.

\emph{Log-rank test in the print:}
Mantel N (1966). Evaluation of survival data and two new rank order statistics arising in its
  consideration. Cancer Chemotherapy Reports 50:163-170.

All references: \code{sx_references("km")}.
}
\seealso{
\code{\link{aj}} (competing risks), \code{\link{logrank_test}}, \code{\link{rmst}},
\code{\link{hazard_np}}, \code{\link{cox}}, \code{\link{Event}}.
}
\examples{
set.seed(1)
n <- 200
d <- data.frame(arm = rep(c("A", "B"), each = n / 2))
t <- rexp(n, ifelse(d$arm == "B", 0.15, 0.1)); cens <- runif(n, 2, 12)
d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)

fit <- km(Event(time, status) ~ arm, d)
fit                                   # medians and the log-rank test
summary(fit, times = c(1, 3, 5))      # survival at 1, 3 and 5 years
plot(fit, median = TRUE)              # curves with a risk table
head(as.data.frame(fit))
}
