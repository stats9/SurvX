\name{ph_test}
\alias{ph_test}
\title{Test and plot of the proportional hazards assumption}
\description{
Tests whether the effect of each covariate of a Cox model is constant over time (proportional
hazards), for every covariate and globally, and plots the estimated time-varying effect. For a
cause-specific Cox model (\code{\link{csc}}) every cause is tested.

\strong{Use it} after every \code{\link{cox}} or Cox-baseline \code{\link{csc}} fit: the hazard
ratios are averages over time when the assumption fails.
}
\usage{
ph_test(fit, transform = c("km", "rank", "identity", "log"))
}
\arguments{
\item{fit}{a model fitted by \code{\link{cox}} or by \code{\link{csc}} with a Cox baseline.}
\item{transform}{the function \eqn{g(t)} of time against which a change of the effect is sought:
  \code{"km"} (default, one minus the Kaplan-Meier estimate just before \eqn{t}; robust to a few
  long times), \code{"rank"} (rank of \eqn{t} among all observed times, censored included, as in
  \code{survival::cox.zph}), \code{"identity"} (time itself) or \code{"log"} (log time; sensitive
  to early changes).}
}
\details{
For each covariate the test adds a term \eqn{\gamma_j x_j g(t)} to the model and tests
\eqn{\gamma = 0} by the score test at the fitted \eqn{\hat\beta}: the score is
\eqn{\sum_i g(t_i) r_{ij}} with the Schoenfeld residuals \eqn{r}, and its variance is the efficient
information for \eqn{\gamma} given \eqn{\beta}. The global test uses all covariates together. This
is the test of Grambsch and Therneau (1994) in its exact score form, as computed by
\code{survival::cox.zph} since version 3.0. Risk-set moments use Breslow weights, so the results
agree with \pkg{survival} exactly when no event times are tied.

The plot shows the scaled Schoenfeld residuals against \eqn{g(t)} with a natural spline smooth and
its pointwise confidence band. The smooth estimates the time-varying effect \eqn{\beta_j(t)}; a flat
line at \eqn{\hat\beta_j} (dashed) supports proportional hazards.
}
\section{Plot}{
\code{plot(x, var = NULL, cause = 1, df = 4, level = 0.95, col = "#0072B2", xlab = NULL,
ylab = NULL, main = NULL, cex = 0.85, ...)}, with:
\describe{
\item{\code{x}}{the result of \code{ph_test()}.}
\item{\code{var}}{covariates to plot; default all.}
\item{\code{cause}}{the cause to plot for cause-specific models.}
\item{\code{df}}{degrees of freedom of the natural spline that smooths the scaled residuals.}
\item{\code{level}}{confidence level of the pointwise band of the smooth.}
\item{\code{col}}{colour of the smooth.}
\item{\code{xlab, ylab, main}}{axis labels and title.}
\item{\code{cex}}{character expansion.}
\item{\code{...}}{passed to \code{plot.default}.}
}
}
\section{Interpreting the result}{
\itemize{
  \item A small p-value for a covariate means that its hazard ratio changes over time. The plot
    shows how: a rising smooth means the effect grows (the HR moves up), a falling one that it
    fades. The global test asks whether any covariate violates the assumption.
  \item With large data sets, tiny and clinically irrelevant departures become significant: judge
    the size of the change on the plot, not only the p-value.
  \item What to do when the assumption fails for a covariate: put it in \code{strata()} if its
    effect is not of interest; let its effect change with time with \code{tt()} in
    \code{\link{cox}}, or let \code{\link{time_effects}} compare forms of the change and refit the
    model with them; or report hazard ratios for separate periods of follow-up. Restricted mean
    survival time (\code{\link{rmst}}) needs no proportional hazards.
}
A sentence for a report: \emph{"The proportional hazards assumption was checked with the
Grambsch-Therneau test (global p = 0.41) and plots of scaled Schoenfeld residuals."}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}), the score statistics agree with \code{survival::cox.zph} for all four
transforms, with and without strata, clusters and ties, to \eqn{10^{-8}}. In a Monte Carlo study
(1000 replications, \eqn{n} = 200 to 1000) the size at the 5\% level was 4.5 to 5.9\% and identical to
that of \code{cox.zph} on the same data; with crossing hazards the power was 100\%.
}
\value{
An S7 object of class \code{PHTest} with \code{print}, \code{plot} and \code{as.data.frame}
methods. \code{x@table} is a data frame of all tests (cause, term, chi-square, df, p; term
\code{GLOBAL} for the global test); \code{x@tests} has one element per cause (one for a Cox model)
with the table and the scaled Schoenfeld residuals used by the plot.
}
\references{
\emph{Test and residual plot:}
Grambsch PM, Therneau TM (1994). Proportional hazards tests and diagnostics based on weighted
  residuals. Biometrika 81:515-526.
Therneau TM, Grambsch PM (2000). Modeling Survival Data: Extending the Cox Model. Springer,
  Chapter 6.

\emph{Schoenfeld residuals:}
Schoenfeld D (1982). Partial residuals for the proportional hazards regression model.
  Biometrika 69:239-241.

All references: \code{sx_references("ph_test")}.
}
\seealso{
\code{\link{cox}} (\code{strata()}, \code{tt()}), \code{\link{time_effects}}, \code{\link{residuals.SurvX}},
\code{\link{rmst}}.
}
\examples{
set.seed(4)
n <- 300
x <- rbinom(n, 1, 0.5)
t <- rweibull(n, ifelse(x == 1, 0.7, 1.5), 3); cens <- runif(n, 0, 6)   # non-proportional hazards
d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = rnorm(n))

pt <- ph_test(cox(Event(time, status) ~ x + z, d))
pt                       # x violates proportional hazards, z does not
plot(pt, "x")            # the effect of x changes over time
as.data.frame(pt)
}
