\name{aj}
\alias{aj}
\title{Aalen-Johansen estimate of the cumulative incidence with competing risks}
\description{
Estimates, without any model, the \emph{cumulative incidence} of every cause,
\eqn{F_k(t) = P(T \le t, \mathrm{cause} = k)}: the probability of having had event \eqn{k} by time
\eqn{t} when other events compete with it. One set of curves per group, with confidence intervals
and Gray's test for every cause when there are groups.

\strong{Use it when} there are competing events (status 0 = censored, 1, 2, ... = causes) and you
want to describe the risk of each cause or compare groups. It replaces "1 minus Kaplan-Meier",
which treats competing events as censoring and overestimates the risk.
\strong{Use something else when} covariates must be adjusted for: \code{\link{fg}} or
\code{\link{directcif}} (cumulative incidence) and \code{\link{csc}} (cause-specific hazards).
}
\usage{
aj(formula, data, level = 0.95)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ 1} or \code{Event(time, status) ~ group}; several
  variables on the right give one set of curves per combination. \code{status}: 0 = censored,
  1, 2, ..., K = cause. Delayed entry is not supported.}
\item{data}{a data frame. Rows with a missing time, status or group are removed.}
\item{level}{confidence level of the intervals (default 0.95).}
}
\details{
\subsection{Estimator}{
\deqn{\hat F_k(t) = \sum_{t_j \le t} \hat S(t_{j-1}) \frac{d_{kj}}{Y_j},}
with \eqn{\hat S} the Kaplan-Meier estimate of being free of all events, \eqn{d_{kj}} the number of
events of cause \eqn{k} and \eqn{Y_j} the number at risk at time \eqn{t_j} (Aalen and Johansen
1978; Kalbfleisch and Prentice 2002, Section 8.2). At every time the cumulative incidences of all
causes and \eqn{\hat S} add up to one.
}
\subsection{Variance and intervals}{
Delta-method variance (Marubini and Valsecchi 1995, Section 10.3). The confidence interval is
computed on the complementary log-log scale, \eqn{\hat F^{\exp(\pm z\,se / (\hat F \log \hat F))}},
which keeps it inside (0, 1) and behaves well for small incidences.
}
}
\section{Interpreting the result}{
\itemize{
  \item \eqn{\hat F_k(t)} is the estimated proportion of subjects who have had cause \eqn{k} by time
    \eqn{t}, in a population where the other causes also occur. It depends on the competing
    causes: if fewer patients die of other causes, more can have cause \eqn{k}.
  \item Gray's test in the print compares the cumulative incidence of each cause across groups
    (see \code{\link{gray_test}}).
}
A sentence for a report: \emph{"The 5-year cumulative incidence of relapse was 24\% (95\% CI 19 to
30\%), and of death without relapse 11\% (8 to 15\%) (Aalen-Johansen estimator)."}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}), on 82 data sets: cumulative incidences agree with \code{cmprsk::cuminc} and
with the multi-state \code{survival::survfit} to \eqn{10^{-14}}, and Gray's statistics with
\code{cuminc} to \eqn{10^{-11}}. The standard errors differ slightly from those of \code{survfit},
which uses the infinitesimal jackknife. In a Monte Carlo study (1000 replications, \eqn{n} = 200
to 1000) the 95\% intervals covered 93.6 to 95.6\%.
}
\value{
An S7 object of class \code{AJFit}. Slots: \code{curves} (one data frame per group with columns
\code{time}, \code{surv}, and for every cause \eqn{k}: \code{cif.k}, \code{se.k}, \code{lower.k},
\code{upper.k}), \code{counts} (numbers censored and per cause), \code{tests} (Gray's test per
cause), \code{K} (number of causes) and \code{level}.

Methods:
\tabular{ll}{
  \code{print(fit)} \tab counts and Gray's tests \cr
  \code{summary(fit, times = c(1, 3, 5))} \tab cumulative incidences with intervals at chosen
    times \cr
  \code{plot(fit, ...)} \tab cumulative incidence curves (\code{\link{cif_plot}}): panels per cause,
    overlaid or stacked, with risk tables \cr
  \code{as.data.frame(fit)} \tab all curves in one data frame \cr
}
}
\references{
\emph{Estimator:}
Aalen OO, Johansen S (1978). An empirical transition matrix for non-homogeneous Markov chains
  based on censored observations. Scandinavian Journal of Statistics 5:141-150.
Kalbfleisch JD, Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd ed.
  Wiley, Section 8.2.

\emph{Variance:}
Marubini E, Valsecchi MG (1995). Analysing Survival Data from Clinical Trials and Observational
  Studies. Wiley, Section 10.3.
Andersen PK, Borgan O, Gill RD, Keiding N (1993). Statistical Models Based on Counting Processes.
  Springer, Section IV.4.

\emph{Gray's test in the print:}
Gray RJ (1988). A class of K-sample tests for comparing the cumulative incidence of a competing
  risk. Annals of Statistics 16:1141-1154.

\emph{Why not 1 minus Kaplan-Meier:}
Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: competing risks and multi-state
  models. Statistics in Medicine 26:2389-2430.

All references: \code{sx_references("aj")}.
}
\seealso{
\code{\link{cif_plot}}, \code{\link{gray_test}}, \code{\link{fg}}, \code{\link{csc}},
\code{\link{km}}, \code{\link{Event}}.
}
\examples{
set.seed(1)
n <- 300
d <- data.frame(arm = rep(c("A", "B"), each = n / 2))
t1 <- rexp(n, ifelse(d$arm == "B", 0.06, 0.1))   # relapse
t2 <- rexp(n, 0.04)                               # death without relapse
cens <- runif(n, 2, 10)
d$time <- pmin(t1, t2, cens)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))

fit <- aj(Event(time, status) ~ arm, d)
fit                                  # Gray's test for each cause
summary(fit, times = c(1, 3, 5))
plot(fit)
}
