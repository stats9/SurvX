\name{ph_test}
\alias{ph_test}
\alias{print.ph_test}
\alias{plot.ph_test}
\title{Test of proportional hazards}
\description{
Score test of the proportional hazards assumption of a Cox model, or of every cause of a
cause-specific Cox model, for each covariate and globally, with the plot of the scaled Schoenfeld
residuals against time.
}
\usage{
ph_test(fit, transform = c("km", "rank", "identity", "log"))

\method{print}{ph_test}(x, ...)

\method{plot}{ph_test}(x, var = NULL, cause = 1L, df = 4L, level = 0.95, col = "#0072B2",
     xlab = NULL, ylab = NULL, main = NULL, cex = 0.85, ...)
}
\arguments{
\item{fit}{a model fitted by \code{\link{cox}} or by \code{\link{csc}} with a Cox baseline.}
\item{transform}{the function \eqn{g(t)} of time: \code{"km"} (default, one minus the Kaplan-Meier
  estimate just before \eqn{t}), \code{"rank"}, \code{"identity"} or \code{"log"}.}
\item{x}{the result of \code{ph_test()}.}
\item{var}{covariates to plot; default all.}
\item{cause}{the cause to plot for cause-specific models.}
\item{df}{degrees of freedom of the natural spline that smooths the scaled residuals.}
\item{level}{confidence level of the pointwise band of the smooth.}
\item{col}{colour of the smooth.}
\item{xlab,ylab,main}{axis labels and title.}
\item{cex}{character expansion.}
\item{...}{passed to \code{plot.default}; not used by \code{print}.}
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
\value{
an object of class \code{"ph_test"}: a list with one element per cause (one for a Cox model), each
  with \code{table} (chi-square, df and p of every covariate and \code{GLOBAL}) and the residuals
  used by the plot.
}
\references{
Grambsch PM, Therneau TM (1994) Proportional hazards tests and diagnostics based on weighted
  residuals. Biometrika 81:515-526.

Therneau TM, Grambsch PM (2000) Modeling Survival Data: Extending the Cox Model. Springer, Chapter 6.
}
\examples{
set.seed(4)
n <- 300
x <- rbinom(n, 1, 0.5)
t <- rweibull(n, ifelse(x == 1, 0.7, 1.5), 3); cens <- runif(n, 0, 6)   # non-proportional hazards
d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = rnorm(n))
pt <- ph_test(cox(Event(time, status) ~ x + z, d))
pt
plot(pt, "x")
}
