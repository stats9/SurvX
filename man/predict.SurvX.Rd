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
\item{\code{times}}{time points. For \code{directcif} fits \code{Inf} gives the plateau.}
\item{\code{type}}{for \code{cox}: \code{"lp"} (linear predictor) or \code{"survival"}.}
\item{\code{cause}}{cause of interest (\code{csc} and \code{directcif} fits).}
\item{\code{se}}{for \code{directcif} fits: also return delta-method standard errors.}
\item{\code{...}}{not used.}
}
}
\value{
a matrix with one row per row of \code{newdata} and one column per time point
  (a vector for \code{type = "lp"}); with \code{se = TRUE} a list with \code{cif} and \code{se}.
}
