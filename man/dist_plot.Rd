\name{dist_plot}
\alias{dist_plot}
\title{Plots for choosing a distribution and checking a fitted model}
\description{
With a formula: the classic linearising plots of the non-parametric (Nelson-Aalen) cumulative
hazard, by group; a straight line in a panel indicates the distribution of that panel. With a
fitted model: the Cox-Snell residual plot.
}
\usage{
dist_plot(x, data = NULL, cause = 1L, col = NULL, cex = 0.85, ...)
}
\arguments{
\item{x}{a formula \code{Event(time, status) ~ 1} or \code{~ group}, or a model fitted by
  \code{\link{parsurv}}, \code{\link{cox}} or \code{\link{csc}} (Cox or parametric baseline).}
\item{data}{a data frame (with a formula).}
\item{cause}{the cause for competing risks data (cause-specific hazard).}
\item{col}{colours of the groups.}
\item{cex}{character expansion.}
\item{...}{passed to \code{plot.default}.}
}
\details{
Panels for a formula, with \eqn{\hat H} the Nelson-Aalen estimate and \eqn{\hat S = e^{-\hat H}}:
\describe{
  \item{exponential}{\eqn{\hat H(t)} against \eqn{t}: a straight line through the origin.}
  \item{Weibull}{\eqn{\log \hat H(t)} against \eqn{\log t}: a straight line (slope = shape;
    parallel lines across groups indicate proportional hazards).}
  \item{log-normal}{\eqn{\Phi^{-1}(1 - \hat S(t))} against \eqn{\log t}: a straight line.}
  \item{log-logistic}{\eqn{\log\{(1 - \hat S)/\hat S\}} against \eqn{\log t}: a straight line
    (parallel lines indicate proportional odds).}
}
For a fitted model, the Nelson-Aalen cumulative hazard of the Cox-Snell residuals (with the event
indicators) is plotted against the residuals. If the model fits, the residuals behave like a
censored sample from the unit exponential distribution and the curve follows the line of identity.
}
\section{Interpreting the plot}{
With a formula: look for the panel in which the points lie closest to a straight line; that
distribution is a good candidate for \code{\link{parsurv}}. Parallel lines across groups in the
Weibull panel indicate proportional hazards; in the log-logistic panel, proportional odds. With a
fitted model: Cox-Snell residuals that follow the identity line support the model; a curve that
bends away shows lack of fit, usually in the tail first.
}
\value{
invisibly, the plotted data.
}
\references{
Klein JP, Moeschberger ML (2003) Survival Analysis: Techniques for Censored and Truncated Data,
  2nd ed. Springer, Sections 12.2-12.5.

Collett D (2015) Modelling Survival Data in Medical Research, 3rd ed. CRC Press, Chapter 5.

Cox DR, Snell EJ (1968) A general definition of residuals. Journal of the Royal Statistical Society,
  Series B 30:248-275.
}
\seealso{
\code{\link{hazard_plot}}, \code{\link{hazard_np}}, \code{\link{residuals.SurvX}},
\code{\link{parsurv}}.
}
\examples{
set.seed(5)
d <- data.frame(time = rweibull(300, 1.5, 3), status = 1, arm = rep(c("A", "B"), 150))
dist_plot(Event(time, status) ~ arm, d)
dist_plot(parsurv(Event(time, status) ~ arm, d))
}
