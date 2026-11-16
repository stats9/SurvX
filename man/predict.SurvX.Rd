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
\value{
a matrix with one row per row of \code{newdata} and one column per time point
  (a vector for \code{type = "lp"}); with \code{se = TRUE} a list with \code{cif} and \code{se}
  (\code{directcif}) or \code{est}, \code{se}, \code{lower} and \code{upper} (\code{parsurv},
  parametric \code{csc}).
}
