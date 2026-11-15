\name{csc}
\alias{csc}
\title{Cause-specific Cox models}
\description{
One Cox model per cause, with the other causes treated as censoring, combined into
cumulative incidence functions.
}
\usage{
csc(formula, data, ties = c("efron", "breslow"), product_limit = TRUE)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}; status 0 = censored, 1..K = cause.}
\item{data}{a data frame.}
\item{ties}{\code{"efron"} (default) or \code{"breslow"}.}
\item{product_limit}{\code{TRUE} (default): overall survival as a product-limit, as in
  \pkg{riskRegression}; \code{FALSE}: \eqn{\exp\{-\sum_k \Lambda_k(t)\}}.}
}
\value{
an S7 object of class \code{CSCFit}; \code{x@models} holds the fit of each cause.
  \code{predict(fit, newdata, times, cause)} gives the cumulative incidence.
}
\references{
Prentice RL et al. (1978) Biometrics 34:541-554.
  Benichou J, Gail MH (1990) Biometrics 46:813-826.
  Ozenne B et al. (2017) R Journal 9(2):440-460.
}
