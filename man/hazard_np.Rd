\name{hazard_np}
\alias{hazard_np}
\title{Non-parametric hazard: Nelson-Aalen and kernel-smoothed hazard}
\description{
The Nelson-Aalen estimate of the cumulative hazard and the kernel-smoothed hazard rate, for every
group of the right-hand side and, with competing risks, for every cause (cause-specific hazards,
the other causes treated as censoring). Delayed entry is allowed.
}
\usage{
hazard_np(formula, data, smooth = TRUE, bw = NULL, boundary = TRUE, level = 0.95,
          grid = 201L)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ 1} or \code{~ group}.}
\item{data}{a data frame.}
\item{smooth}{compute the kernel-smoothed hazard.}
\item{bw}{bandwidth; \code{NULL} chooses it by cross-validation, separately for every curve.}
\item{boundary}{use boundary kernels near the ends of follow-up.}
\item{level}{confidence level.}
\item{grid}{number of time points of the smoothed curve.}
}
\details{
The cumulative hazard is \eqn{\hat H(t) = \sum_{t_j \le t} d_j / Y_j} with variance
\eqn{\sum d_j / Y_j^2}; limits on the log scale. The smoothed hazard is
\eqn{\hat h(t) = b^{-1}\sum_j K((t - t_j)/b)\, d_j / Y_j} with the Epanechnikov kernel and, within
one bandwidth of either end of follow-up, the boundary kernels of Mueller and Wang (1994), so the
estimate is not biased downwards at time zero. Its variance is \eqn{b^{-2}\sum_j K^2 d_j / Y_j^2};
limits on the log scale. The smoothed hazard is shown up to the last event time with at least
max(10, 10\% of the initial risk set) subjects at risk; beyond it the estimate rests on a handful of
subjects. The bandwidth \eqn{b} is chosen on that range from 40 values (1/50 to 1/2 of it) as the
largest local minimiser of the least-squares cross-validation criterion of Patil (1993), which is
far more stable than its global minimiser (Hall and Marron 1991), unless \code{bw} is given.
}
\value{
an S7 object of class \code{HazardFit}: \code{x@cumhaz} and \code{x@hazard} are lists of data
  frames, one per group and cause; \code{x@bandwidth} the bandwidths. Methods \code{print} and
  \code{plot(x, type = c("hazard", "cumhaz"))}.
}
\references{
Nelson W (1972) Theory and applications of hazard plotting for censored failure data.
  Technometrics 14:945-966.

Aalen O (1978) Nonparametric inference for a family of counting processes. Annals of Statistics
  6:701-726.

Ramlau-Hansen H (1983) Smoothing counting process intensities by means of kernel functions.
  Annals of Statistics 11:453-466.

Mueller HG, Wang JL (1994) Hazard rate estimation under random censoring with varying kernels and
  bandwidths. Biometrics 50:61-76.

Patil PN (1993) On the least squares cross-validation bandwidth in hazard rate estimation. Annals
  of Statistics 21:1792-1810.

Hall P, Marron JS (1991) Local minima in cross-validation functions. Journal of the Royal
  Statistical Society, Series B 53:245-252.
}
\examples{
set.seed(3)
d <- data.frame(time = rweibull(300, 1.5, 3), status = 1, arm = rep(c("A", "B"), 150))
hz <- hazard_np(Event(time, status) ~ arm, d)
hz
plot(hz)
plot(hz, type = "cumhaz")
}
