\name{predict.SurvX}
\alias{predict.SurvX}
\title{Predictions from SurvX models}
\description{
\code{predict()} methods of the model classes.
}
\section{Arguments}{
\describe{
\item{\code{object}}{a fitted model.}
\item{\code{newdata}}{data frame with the covariates; the data used for fitting when \code{NULL}.}
\item{\code{times}}{time points. For \code{directcif} fits \code{Inf} gives the plateau; for
  \code{parsurv} with \code{type = "rmst"} the upper limits of integration.}
\item{\code{type}}{for \code{cox}: \code{"lp"} (linear predictor) or \code{"survival"}. For
  \code{parsurv}: \code{"survival"}, \code{"hazard"}, \code{"cumhaz"}, \code{"quantile"} (at the
  probabilities \code{p}) or \code{"rmst"} (restricted mean survival time). For \code{csc} with a
  parametric baseline: \code{"cif"} or \code{"survival"} (event-free probability).}
\item{\code{p}}{\code{parsurv}, \code{type = "quantile"}: probabilities; the time by which a
  fraction \code{p} has failed (\code{p = 0.5}: the median).}
\item{\code{cause}}{cause of interest (\code{csc} and \code{directcif} fits).}
\item{\code{se}}{for \code{directcif}, \code{parsurv} and parametric \code{csc} fits: also return
  delta-method standard errors (and, for \code{parsurv} and parametric \code{csc}, pointwise
  limits: complementary log-log scale for probabilities, log scale for positive quantities).}
\item{\code{level}}{confidence level of the limits.}
\item{\code{average}}{\code{parsurv} and parametric \code{csc}: the average over the rows of
  \code{newdata} (marginal curve), with its standard error.}
\item{\code{...}}{not used.}
}
}
\details{
Cox models: \eqn{S(t \mid x) = \exp\{-\hat H_0(t) e^{x'\hat\beta}\}} with the Breslow
estimate of the baseline cumulative hazard. Cause-specific Cox models: the cumulative incidence
\eqn{\int_0^t \hat S(u-)\, d\hat H_k(u \mid x)} from all cause-specific hazards (product-limit
or exponential form of \eqn{\hat S}, see \code{product_limit} of \code{\link{csc}}). Fine-Gray
models: \eqn{1 - \exp\{-\hat\Lambda_{10}(t) e^{x'\hat\beta}\}}. Direct and parametric models
evaluate the fitted distribution; parametric cause-specific models integrate the hazards
numerically (Gauss-Legendre). Standard errors of direct and parametric predictions use the delta
method with the covariance matrix of the parameters; limits are computed on a transformed scale
(log-log for probabilities, log for hazards and times) and transformed back.
}
\references{
Breslow NE (1974) Covariance analysis of censored survival data. Biometrics 30:89-99.

Prentice RL, Kalbfleisch JD, Peterson AV, Flournoy N, Farewell VT, Breslow NE (1978) The analysis
  of failure times in the presence of competing risks. Biometrics 34:541-554.

Benichou J, Gail MH (1990) Estimates of absolute cause-specific risk in cohort studies. Biometrics
  46:813-826.

Fine JP, Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk.
  Journal of the American Statistical Association 94:496-509.

Kalbfleisch JD, Prentice RL (2002) The Statistical Analysis of Failure Time Data, 2nd ed. Wiley.

Hinchliffe SR, Lambert PC (2013) Flexible parametric modelling of cause-specific hazards to
  estimate cumulative incidence functions. BMC Medical Research Methodology 13:13.

Casella G, Berger RL (2002) Statistical Inference, 2nd ed. Duxbury, Section 5.5 (delta method).
}
\value{
a matrix with one row per row of \code{newdata} and one column per time point
  (a vector for \code{type = "lp"}); with \code{se = TRUE} a list with \code{cif} and \code{se}
  (\code{directcif}) or \code{est}, \code{se}, \code{lower} and \code{upper} (\code{parsurv},
  parametric \code{csc}).
}
