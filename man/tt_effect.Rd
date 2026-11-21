\name{tt_effect}
\alias{tt_effect}
\title{Hazard ratio at chosen times for time-varying coefficients}
\description{
For a Cox model with \code{tt()} terms (see \code{\link{cox}}), the hazard ratio of each such
covariate at the requested times, \eqn{\exp\{\beta + \gamma\, g(t)\}}, with a Wald confidence
interval. It turns the two coefficients of a time-varying effect into the quantity to report:
"the hazard ratio of treatment was 0.55 at 6 months and 0.85 at 3 years".
}
\usage{
tt_effect(fit, times = NULL, level = 0.95)
}
\arguments{
\item{fit}{a model fitted by \code{\link{cox}} with \code{tt()} terms.}
\item{times}{the times at which the effect is wanted; default the quartiles of the event times.}
\item{level}{confidence level.}
}
\details{
For a covariate \eqn{x} entered as \code{x + tt(x)} with \code{tt = g}, the log hazard ratio at
time \eqn{t} is \eqn{\beta + \gamma g(t)} with variance
\eqn{V_{\beta\beta} + g(t)^2 V_{\gamma\gamma} + 2 g(t) V_{\beta\gamma}}. When \code{x} itself is not
in the model, \eqn{\beta = 0}.
}
\value{
A data frame with \code{term}, \code{time}, \code{g} (the value of the function of time),
\code{log_hr}, \code{se}, \code{hr}, \code{lower} and \code{upper}.
}
\references{
Therneau TM, Grambsch PM (2000). Modeling Survival Data: Extending the Cox Model. Springer,
  Chapter 6.
}
\seealso{\code{\link{cox}}, \code{\link{ph_test}}.}
\examples{
set.seed(4)
n <- 300
x <- rbinom(n, 1, 0.5)
t <- rweibull(n, ifelse(x == 1, 0.7, 1.5), 3); cens <- runif(n, 0, 6)
d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x)
fit <- cox(Event(time, status) ~ x + tt(x), d, tt = function(t) log(t))
tt_effect(fit, times = c(0.5, 1, 2, 4))
}
