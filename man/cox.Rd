\name{cox}
\alias{cox}
\title{Cox proportional hazards model}
\description{
Partial-likelihood fit by Newton-Raphson with step halving.
}
\usage{
cox(formula, data, ties = c("efron", "breslow"), weights = NULL, max_iter = 50L,
    tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}; status 0/1.}
\item{data}{a data frame.}
\item{ties}{\code{"efron"} (default) or \code{"breslow"}.}
\item{weights}{optional case weights, one per row of \code{data}.}
\item{max_iter,tol}{iteration limit and convergence tolerance.}
}
\value{
an S7 object of class \code{CoxFit}. Methods: \code{print}, \code{summary}, \code{coef},
  \code{vcov}, \code{logLik}, \code{nobs}, \code{predict}.
}
\references{
Cox DR (1972) J R Stat Soc B 34:187-220. Breslow N (1974) Biometrics 30:89-99.
  Efron B (1977) J Am Stat Assoc 72:557-565.
}
