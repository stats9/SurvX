\name{cif_data}
\alias{cif_data}
\alias{cif_plot}
\title{Cumulative incidence curves of fitted models}
\description{
Computes the cumulative incidence function of every cause for one or several models, at chosen
covariate values, at typical covariate values, as the average over the sample, or from the same
models fitted without covariates. \code{cif_plot()} draws the result; \code{cif_data()} returns it
as a data frame for tables or other graphics.
}
\usage{
cif_data(x, newdata = NULL, cause = NULL,
    type = c("profile", "marginal", "unadjusted"), times = NULL, conf.int = TRUE,
    level = 0.95, B = 200L, seed = 1L, interval = c("cloglog", "percentile"))

cif_plot(x, newdata = NULL, cause = NULL,
    type = c("profile", "marginal", "unadjusted"),
    layout = c("panels", "overlay", "stacked"), times = NULL, conf.int = TRUE,
    level = 0.95, B = 200L, seed = 1L, interval = c("cloglog", "percentile"),
    risk.table = TRUE, legend = TRUE, caption = TRUE, reference = FALSE,
    mark.times = NULL, percent = FALSE, grid = FALSE, labels = NULL,
    cause.names = NULL, col = NULL, lty = NULL, lwd = 2, xlab = "Time", ylab = NULL,
    main = NULL, xlim = NULL, ylim = NULL, cex = 0.85,
    engine = c("base", "ggplot2"), ...)
}
\arguments{
\item{x}{a model fitted by \code{\link{directcif}}, \code{\link{fg}} or \code{\link{csc}} (Cox or parametric baseline), an
  \code{\link{aj}} estimate, or a named list of these for a comparison.}
\item{newdata}{data frame of covariate values, one curve per row.}
\item{cause}{causes to show; default all.}
\item{type}{\code{"profile"}, \code{"marginal"} or \code{"unadjusted"}; see Details.}
\item{times}{time points; default a fine grid (direct models) or the event times (step functions).}
\item{conf.int}{compute pointwise confidence intervals.}
\item{level}{confidence level.}
\item{B}{number of bootstrap samples for Fine-Gray and cause-specific Cox models.}
\item{seed}{seed of the bootstrap (the same seed gives the same intervals).}
\item{interval}{\code{"cloglog"} or \code{"percentile"} for bootstrap intervals.}
\item{layout}{\code{"panels"} (one panel per cause), \code{"overlay"} (all causes in one panel) or
  \code{"stacked"} (the cumulative incidences of the causes stacked on each other, with the
  event-free probability on top; one panel per model and covariate profile).}
\item{risk.table}{show the number of subjects at risk under the time axis.}
\item{legend}{show the legend.}
\item{caption}{show which covariate values were used.}
\item{reference}{add the Aalen-Johansen estimate of the whole sample for comparison.}
\item{mark.times}{time points at which the values of the curves are written on the plot.}
\item{percent}{label the probability axis in percent.}
\item{grid}{draw light horizontal grid lines.}
\item{labels}{labels of the curves (model and profile combinations), in the order of the legend.}
\item{cause.names}{names of the causes, e.g. \code{c("relapse", "death")}.}
\item{col,lty,lwd}{colours, line types and width of the curves.}
\item{xlab,ylab,main,xlim,ylim}{axis labels, title and limits.}
\item{cex}{size of all text, relative to the device font (default 0.85). With
  \code{engine = "ggplot2"} the theme's base size is \code{12 * cex} points.}
\item{engine}{\code{"base"} (default, no dependencies) or \code{"ggplot2"} (returns a ggplot object;
  the table of numbers at risk is not drawn).}
\item{...}{passed to \code{plot.default} (base engine).}
}
\details{
\describe{
  \item{\code{type = "profile"} (default)}{one curve at typical covariate values: the mean of every
    numeric covariate and the most frequent level of every categorical covariate.}
  \item{\code{newdata}}{one curve for every row; \code{type} is then ignored.}
  \item{\code{type = "marginal"}}{the average of the predicted curves of all subjects.}
  \item{\code{type = "unadjusted"}}{the same model fitted without covariates (for a direct model the
    parametric curve of Jeong and Fine 2006; for Fine-Gray and cause-specific Cox the models with all
    coefficients zero).}
}
Curves of all causes are computed: a Fine-Gray model is refitted for every other cause and a
weighted direct model is refitted for every other cause, with the same covariates.

Pointwise confidence intervals: delta method for direct models, parametric cause-specific models and the Aalen-Johansen estimator,
nonparametric bootstrap of the subjects for Fine-Gray and cause-specific Cox models (Efron and
Tibshirani 1993); both on the complementary log-log scale, \eqn{F^{\exp(\pm z\,se/(F \log F))}},
so that the limits stay in [0, 1]. \code{interval = "percentile"} gives bootstrap percentile limits instead.
}
\section{Interpreting the curves and intervals}{
\itemize{
  \item A \emph{profile} curve is the risk of one kind of subject (the typical one, or the rows of
    \code{newdata}); a \emph{marginal} curve is the average risk of the population of the data
    (standardisation), the natural choice for comparing treatments "as if everybody had received
    each"; an \emph{unadjusted} curve ignores the covariates and shows how the model compares with
    the data.
  \item The intervals are pointwise: each covers the true incidence at its own time with the stated
    probability; they are not a band for the whole curve. Parametric models (\code{directcif},
    parametric \code{csc}) use the delta method; Fine-Gray and Cox-baseline \code{csc} fits use the
    bootstrap (\code{B} resamples of subjects, or of clusters with \code{cluster()}), with
    \code{seed} for reproducibility.
  \item \code{reference = TRUE} adds the Aalen-Johansen estimate: a model whose unadjusted curve
    departs from it fits the data poorly.
}
}
\value{
a data frame with the columns \code{model, cause, profile, time, cif, se, lower, upper, step,
  method}; attribute \code{"description"} says which covariate values were used.
}
\references{
Jeong J-H, Fine JP (2006) J R Stat Soc C 55:187-200.
Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap. Chapman and Hall, New York.
Kalbfleisch JD, Prentice RL (2002) The Statistical Analysis of Failure Time Data, 2nd ed. Wiley,
  Section 8.2 (cumulative incidence functions).
}
\seealso{
\code{\link{aj}}, \code{\link{fg}}, \code{\link{csc}}, \code{\link{directcif}},
\code{\link{predict.SurvX}}, \code{\link{hazard_plot}}.
}
\examples{
set.seed(6)
n <- 300
x <- rbinom(n, 1, 0.5)
t1 <- rexp(n, 0.4 * exp(0.6 * x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
d <- data.frame(time = pmin(t1, t2, cens), x = x)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
fit <- directcif(Event(time, status) ~ x, d)
head(cif_data(fit))
cif_plot(fit)
cif_plot(fit, newdata = data.frame(x = 0:1), cause = 1, risk.table = FALSE)
cif_plot(list("direct PH" = fit, "Fine-Gray" = fg(Event(time, status) ~ x, d)), B = 50)
cif_plot(fit, layout = "stacked", conf.int = FALSE)
}
