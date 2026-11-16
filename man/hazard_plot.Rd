\name{hazard_plot}
\alias{hazard_plot}
\title{Fitted hazard, cumulative hazard and survival of parametric models}
\description{
Draws the fitted hazard (or cumulative hazard, or survival) of one or several parametric models at
typical or chosen covariate values, with delta-method confidence bands, and overlays the
non-parametric estimate of the same data (\code{\link{hazard_np}}, Nelson-Aalen or
Kaplan-Meier) as a reference. Agreement of the fitted and the non-parametric curve is the visual
check of the chosen distribution.
}
\usage{
hazard_plot(x, newdata = NULL, type = c("hazard", "cumhaz", "survival"),
            cause = NULL, reference = TRUE, conf.int = TRUE, level = 0.95, times = NULL,
            col = NULL, lwd = 2, xlab = "Time", ylab = NULL, main = NULL, legend = TRUE,
            cex = 0.85, ...)
}
\arguments{
\item{x}{a model fitted by \code{\link{parsurv}} or by \code{\link{csc}} with a parametric
  baseline, or a named list of such models fitted to the same data.}
\item{newdata}{data frame of covariate values, one curve per row; default typical values (mean of
  numeric covariates, most frequent level of categorical covariates).}
\item{type}{\code{"hazard"}, \code{"cumhaz"} or \code{"survival"}.}
\item{cause}{causes to show (cause-specific models); default all.}
\item{reference}{overlay the non-parametric estimate of all subjects. It is a fair check of the
  fitted curve at typical covariate values or without covariates.}
\item{conf.int}{draw pointwise confidence bands.}
\item{level}{confidence level.}
\item{times}{time grid; default 201 points up to the largest observed time.}
\item{col,lwd}{colours and line width of the models.}
\item{xlab,ylab,main}{axis labels and title.}
\item{legend}{draw a legend.}
\item{cex}{character expansion.}
\item{...}{passed to \code{plot.default}.}
}
\value{
invisibly, a data frame of the plotted curves.
}
\references{
Klein JP, Moeschberger ML (2003) Survival Analysis: Techniques for Censored and Truncated Data,
  2nd ed. Springer, Chapter 12.
}
\examples{
set.seed(1)
n <- 300
t <- rweibull(n, shape = 1.5, scale = 4); cens <- runif(n, 0, 8)
d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = rnorm(n))
fits <- list(Weibull = parsurv(Event(time, status) ~ x, d),
             exponential = parsurv(Event(time, status) ~ x, d, baseline = "exponential"))
hazard_plot(fits)
hazard_plot(fits, type = "cumhaz")
}
