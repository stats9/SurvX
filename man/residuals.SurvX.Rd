\name{residuals.SurvX}
\alias{residuals.SurvX}
\title{Residuals of hazard models}
\description{
Martingale, deviance, Cox-Snell and (scaled) Schoenfeld residuals of Cox models, cause-specific
models, parametric models and cure models: \code{residuals(fit, type)}. They are used to check
the functional form of covariates, to find poorly fitted subjects, to check the overall fit and
the proportional hazards assumption.
}
\section{Usage}{
\code{residuals(object, type = c("martingale", "deviance", "coxsnell", "schoenfeld", "scaledsch"))}
}
\section{Arguments}{
\describe{
\item{\code{object}}{a model fitted by \code{\link{cox}}, \code{\link{csc}} (Cox or parametric
  baseline), \code{\link{parsurv}} or \code{\link{cure}}.}
\item{\code{type}}{\code{"martingale"} (default), \code{"deviance"}, \code{"coxsnell"},
  \code{"schoenfeld"} or \code{"scaledsch"}; the last two for Cox models only.}
}
}
\details{
With \eqn{\hat H_i} the fitted cumulative hazard of subject \eqn{i} at its time (minus that at
entry), the Cox-Snell residual is \eqn{\hat H_i}, the martingale residual
\eqn{M_i = d_i - \hat H_i} and the deviance residual
\eqn{\mathrm{sign}(M_i)\sqrt{-2\{M_i + d_i \log(d_i - M_i)\}}}. Schoenfeld residuals are defined at
every event: the covariate vector minus its risk-weighted mean over the subjects at risk. For Cox
models the cumulative hazard and the means take the form of the ties method of the fit: with Efron
ties a subject who fails at a tied time contributes the Efron-weighted hazard increments, and the
mean is averaged over the Efron terms, exactly as in \pkg{survival} (Therneau and Grambsch 2000,
Section 3.3). The martingale residuals of a Cox model then sum to zero. The scaled residuals are \eqn{\hat\beta + D\,V r}, with \eqn{D} the number of events and
\eqn{V} the covariance matrix of \eqn{\hat\beta}, as in \pkg{survival}; their smoothed trend
against time estimates \eqn{\beta(t)} (Grambsch and Therneau 1994), see \code{\link{ph_test}}.
}
\section{What each residual is for}{
\describe{
  \item{Martingale}{observed minus expected number of events, between \eqn{-\infty} and 1.
    Plotted against a numeric covariate that is \emph{not} in the model, a smooth curve shows the
    functional form it should have (a straight line: enter it linearly; a threshold or a curve:
    categorise, transform or use a spline) (Therneau, Grambsch and Fleming 1990).}
  \item{Deviance}{a symmetrised martingale residual; values beyond about \eqn{\pm 2.5} point to
    subjects poorly predicted by the model (an early death of a low-risk patient, for example).}
  \item{Cox-Snell}{if the model is correct they behave like a censored sample from the unit
    exponential distribution: \code{\link{dist_plot}(fit)} plots their cumulative hazard, which
    should follow the identity line (Cox and Snell 1968).}
  \item{Schoenfeld and scaled Schoenfeld}{one row per event; their trend over time estimates the
    time-varying effect \eqn{\beta(t)} and is the basis of \code{\link{ph_test}}.}
}
}
\section{Validation}{
In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
\code{validation}), martingale, deviance and scaled Schoenfeld residuals of Cox models agree with
\code{survival::residuals.coxph} to \eqn{10^{-8}}, with Efron and Breslow ties, strata and
weights; Efron martingale residuals sum to zero to \eqn{2 \times 10^{-13}}.
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
\seealso{\code{\link{ph_test}}, \code{\link{dist_plot}}, \code{\link{cox}}, \code{\link{parsurv}}.}
\examples{
set.seed(1)
n <- 300
d <- data.frame(age = runif(n, 40, 80), x = rbinom(n, 1, 0.5))
t <- rexp(n, 0.05 * exp(0.5 * d$x + 0.002 * (d$age - 60)^2)); cens <- runif(n, 0, 20)
d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)

## functional form of age: martingale residuals of a model without age
f0 <- cox(Event(time, status) ~ x, d)
plot(d$age, residuals(f0), xlab = "age", ylab = "martingale residual")
lines(stats::lowess(d$age, residuals(f0)), lwd = 2)      # U-shaped: age is not linear

## poorly fitted subjects
f1 <- cox(Event(time, status) ~ x + splines::ns(age, 3), d)
which(abs(residuals(f1, "deviance")) > 2.5)

## overall fit: Cox-Snell residuals
dist_plot(f1)
}
