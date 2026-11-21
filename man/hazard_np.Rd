\name{hazard_np}
\alias{hazard_np}
\title{Non-parametric hazard: Nelson-Aalen and kernel-smoothed hazard rate}
\description{
Estimates, without any model, the cumulative hazard (Nelson-Aalen) and the hazard rate itself
(kernel smoothing), for every group and, with competing risks, for every cause (cause-specific
hazards). Delayed entry is allowed.

\strong{Use it when} you want to see the \emph{shape} of the hazard over time (rising, falling,
peaking), to choose a parametric distribution for \code{\link{parsurv}} or \code{\link{csc}}, to
check a fitted parametric model (\code{\link{hazard_plot}} overlays both), or when there is
delayed entry and \code{\link{km}} cannot be used.
}
\usage{
hazard_np(formula, data, smooth = TRUE, bw = NULL, boundary = TRUE, level = 0.95,
          grid = 201L)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ 1} or \code{~ group};
  \code{Event(time, status, entry = entry)} for delayed entry. With status 1, 2, ... one set of
  curves per cause.}
\item{data}{a data frame.}
\item{smooth}{\code{TRUE} (default): also compute the kernel-smoothed hazard rate.}
\item{bw}{bandwidth of the smoother, in time units. \code{NULL} (default) chooses it by
  cross-validation, separately for every curve. A larger bandwidth gives a smoother but more biased
  curve; try a few values if the automatic choice looks too wiggly or too flat.}
\item{boundary}{\code{TRUE} (default): boundary kernels near the ends of follow-up, so that the
  estimate is not biased downwards at time zero.}
\item{level}{confidence level of the pointwise intervals.}
\item{grid}{number of time points at which the smoothed hazard is computed.}
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
\section{Interpreting the result}{
The hazard \eqn{h(t)} is the instantaneous event rate among subjects still at risk, in events per
unit of time (0.05 per year: about 5 events per 100 person-years at that time). Its shape guides
the choice of distribution: constant (exponential), monotone (Weibull, Gompertz), rising then
falling (log-normal, log-logistic), or more complex (\code{"rp"} spline). The cumulative hazard is
easier to estimate and its slope is the hazard. The intervals of the smoothed hazard are wide in
the tails; do not over-interpret wiggles there.
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}), the Nelson-Aalen cumulative hazard agrees with
\code{survival::survfit(..., ctype = 1)} to rounding error (\eqn{10^{-15}}) on 81 data sets.
}
\value{
An S7 object of class \code{HazardFit}: \code{fit@cumhaz} and \code{fit@hazard} are lists of data
frames (one per group and cause) with the estimates and limits; \code{fit@bandwidth} holds the
bandwidths. Methods: \code{print}, and \code{plot(fit, type = c("hazard", "cumhaz"))}.
}
\references{Nelson W (1972) Theory and applications of hazard plotting for censored failure data.
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

All references: \code{sx_references("hazard_np")}.
}
\seealso{
\code{\link{hazard_plot}} (fitted against non-parametric), \code{\link{dist_plot}}, \code{\link{km}},
\code{\link{parsurv}}.
}
\examples{
set.seed(3)
d <- data.frame(time = rweibull(300, 1.5, 3), status = 1, arm = rep(c("A", "B"), 150))
hz <- hazard_np(Event(time, status) ~ arm, d)
hz                                # bandwidths chosen
plot(hz)                          # a rising hazard: Weibull with shape > 1
plot(hz, type = "cumhaz")

## a fixed, wider bandwidth for a smoother curve
plot(hazard_np(Event(time, status) ~ 1, d, bw = 1.5))
}
