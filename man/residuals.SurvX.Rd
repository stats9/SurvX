\name{residuals.SurvX}
\alias{residuals.SurvX}
\title{Residuals of hazard models}
\description{
Martingale, deviance, Cox-Snell and (scaled) Schoenfeld residuals of Cox models, cause-specific
models and parametric models: \code{residuals(fit, type)}.
}
\section{Arguments}{
\describe{
\item{\code{object}}{a model fitted by \code{\link{cox}}, \code{\link{csc}} (Cox or parametric
  baseline) or \code{\link{parsurv}}.}
\item{\code{type}}{\code{"martingale"} (default), \code{"deviance"}, \code{"coxsnell"},
  \code{"schoenfeld"} or \code{"scaledsch"}; the last two for Cox models only.}
}
}
\details{
With \eqn{\hat H_i} the fitted cumulative hazard of subject \eqn{i} at its time (minus that at
entry), the Cox-Snell residual is \eqn{\hat H_i}, the martingale residual
\eqn{M_i = d_i - \hat H_i} and the deviance residual
\eqn{\mathrm{sign}(M_i)\sqrt{-2\{M_i + d_i \log(d_i - M_i)\}}}. Schoenfeld residuals are defined at
every event: the covariate vector minus its risk-weighted mean over the subjects at risk (Breslow
weights). The scaled residuals are \eqn{\hat\beta + D\,V r}, with \eqn{D} the number of events and
\eqn{V} the covariance matrix of \eqn{\hat\beta}, as in \pkg{survival}; their smoothed trend
against time estimates \eqn{\beta(t)} (Grambsch and Therneau 1994), see \code{\link{ph_test}}.
}
\value{
a vector with one value per subject, or for Schoenfeld residuals a matrix with one row per event in
  time order (attribute \code{"time"}); for cause-specific models a list with one element per cause.
}
\references{
Cox DR, Snell EJ (1968) A general definition of residuals. Journal of the Royal Statistical Society,
  Series B 30:248-275.

Schoenfeld D (1982) Partial residuals for the proportional hazards regression model. Biometrika
  69:239-241.

Therneau TM, Grambsch PM, Fleming TR (1990) Martingale-based residuals for survival models.
  Biometrika 77:147-160.

Grambsch PM, Therneau TM (1994) Proportional hazards tests and diagnostics based on weighted
  residuals. Biometrika 81:515-526.
}
\seealso{\code{\link{ph_test}}, \code{\link{dist_plot}}.}
