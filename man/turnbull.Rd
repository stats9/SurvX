\name{turnbull}
\alias{turnbull}
\alias{TurnbullFit}
\title{Non-parametric survival curves for interval-censored data (Turnbull NPMLE)}
\description{
The non-parametric maximum likelihood estimate (NPMLE) of the survival function when event times
are only known to lie in intervals (\code{\link{Interval}}): the interval-censored counterpart of
the Kaplan-Meier estimate, to which it reduces when all times are exact or right-censored. Curves
can be computed by group, with bootstrap confidence bands, medians, and a logrank-type test of equal
survival for interval-censored data.
}
\usage{
turnbull(formula, data, times = NULL, B = 200L, level = 0.95, test = TRUE,
    seed = 1L, max_iter = 100000L, tol = 1e-10)
}
\arguments{
\item{formula}{\code{Interval(left, right) ~ 1} for one curve, or \code{Interval(left, right) ~ g}
  for one curve per level of \code{g} (several variables: one curve per combination).}
\item{data}{a data frame.}
\item{times}{times at which the survival table is printed (default: about five round times over
  the follow-up).}
\item{B}{number of bootstrap samples for the confidence bands and the confidence interval of the
  median (resampling subjects within groups); \code{0} for none.}
\item{level}{confidence level.}
\item{test}{with groups, compute the test of equal survival.}
\item{seed}{seed of the random numbers of the bootstrap (the global random state is left as it was).}
\item{max_iter,tol}{maximum number of iterations and convergence tolerance of the algorithm.}
}
\details{
\strong{Support.} The NPMLE puts all its probability on the innermost (Turnbull) intervals:
the intervals between a left end point and the next right end point of the data with no other end
point in between (Peto 1973; Turnbull 1976). How the probability is spread within an innermost
interval is not identified by the data; the survival curve is therefore only determined outside
these intervals. \code{predict()} gives the lower and upper limit in addition to the estimate, and
the curve is drawn linearly within the innermost intervals (the usual convention); exact event
times are point masses, as in the Kaplan-Meier estimate. After the last left end point of the
right-censored subjects the curve is undefined (as the Kaplan-Meier curve after the last time).

\strong{Computation.} The self-consistency (EM) algorithm of Turnbull (1976) in compiled code,
accelerated with SQUAREM (Varadhan and Roland 2008). Convergence is checked on the Kuhn-Tucker
conditions, which are necessary and sufficient for the maximum of this concave problem (Gentleman
and Geyer 1994): the printed value is zero at the NPMLE.

\strong{Uncertainty.} The NPMLE converges at a rate slower than \eqn{\sqrt{n}} and has no simple
standard error; the confidence bands are pointwise bootstrap percentile intervals.

\strong{Test.} With groups, the weighted logrank test for interval-censored data with the scores
of Finkelstein (1986), \eqn{c_i = \{S(L_i)\log S(L_i) - S(R_i)\log S(R_i)\}/\{S(L_i) - S(R_i)\}}
from the NPMLE of all subjects together, and the variance from the permutation distribution (the
permutational central limit theorem, \code{method = "pclt"} in Fay and Shaw 2010). It reduces to
the ordinary logrank test (permutation form) for right-censored data. A positive group score means
that the events happen earlier in that group than in the groups together.
}
\section{Messages}{
\describe{
  \item{\code{did not converge}}{the iteration limit was reached before the Kuhn-Tucker condition
    was met; increase \code{max_iter}.}
}
}
\value{
An S7 object of class \code{TurnbullFit}. Use \code{print()}, \code{summary()} (adds the innermost
intervals and their probabilities), \code{plot()}, \code{predict(fit, times)} and
\code{as.data.frame(fit, what = c("survival", "intervals", "test"))}.
}
\references{
Turnbull BW (1976). The empirical distribution function with arbitrarily grouped, censored and
  truncated data. Journal of the Royal Statistical Society B 38:290-295.

Peto R (1973). Experimental survival curves for interval-censored data. Applied Statistics
  22:86-91.

Gentleman R, Geyer CJ (1994). Maximum likelihood for interval censored data: consistency and
  computation. Biometrika 81:618-623.

Varadhan R, Roland C (2008). Simple and globally convergent methods for accelerating the
  convergence of any EM algorithm. Scandinavian Journal of Statistics 35:335-353.

Finkelstein DM (1986). A proportional hazards model for interval-censored failure time data.
  Biometrics 42:845-854.

Fay MP, Shaw PA (2010). Exact and asymptotic weighted logrank tests for interval censored data:
  the interval R package. Journal of Statistical Software 36(2):1-34.
}
\seealso{
\code{\link{Interval}}, \code{\link{parsurv}} (parametric models for the same data; its
\code{plot()} draws this estimate for comparison), \code{\link{km}}.
}
\examples{
set.seed(2)
n <- 150
g <- factor(rep(c("A", "B"), each = n / 2))
t <- rweibull(n, 1.5, ifelse(g == "A", 10, 7))
visit <- 3 * runif(n)                     # first visit; then every 3 time units
L <- R <- numeric(n)
for (i in seq_len(n)) {
  v <- visit[i] + 3 * (0:5)
  L[i] <- max(c(0, v[v < t[i]]))
  R[i] <- if (any(v >= t[i])) min(v[v >= t[i]]) else Inf
}
d <- data.frame(L, R, g)
fit <- turnbull(Interval(L, R) ~ g, data = d, B = 50)
fit
plot(fit)
predict(fit, times = c(5, 10))
}
