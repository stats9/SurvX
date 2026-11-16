\name{SurvX-methods}
\alias{SurvX-methods}
\title{Methods for fitted models}
\description{
Every model of SurvX is an S7 object: its parts are read with \code{@} (\code{fit@coefficients},
\code{fit@vcov}, \code{fit@converged}; \code{S7::props(fit)} lists them all). The usual S3 generics
of R work as for \code{lm()} or \code{survival::coxph()}.

\describe{
  \item{\code{print(fit)}}{covariate effects as hazard ratios, subdistribution hazard ratios or odds
    ratios with confidence intervals and p-values; for direct models also the long-run risk
    (plateau) at typical covariate values and a warning when the fit did not converge.}
  \item{\code{summary(fit, level = 0.95)}}{everything in \code{print} and in addition the predicted
    cumulative incidence at the quartiles of the event times, the baseline parameters on their
    natural scale and a note on how to read the effects. The result is an object whose
    \code{@coefficients} is the table of \code{\link{estimates}} and \code{@tables} the other tables;
    \code{as.data.frame()} returns the coefficient table.}
  \item{\code{coef}, \code{vcov}, \code{confint}, \code{logLik}, \code{AIC}, \code{BIC}, \code{nobs}}{as
    usual; coefficients are on the scale used for estimation (log ratios, log or logit baseline
    parameters). For \code{csc()} fits \code{coef}, \code{vcov} and \code{confint} return one element per cause.}
  \item{\code{estimates(fit)}}{all of this in one data frame, including the exponentiated effects.}
  \item{\code{predict(fit, newdata, times, cause)}}{absolute risks; see \code{\link{predict.SurvX}}.}
  \item{\code{plot(fit, ...)}}{for \code{directcif()}, \code{fg()}, \code{csc()} and \code{aj()} the same as
    \code{\link{cif_plot}(fit, ...)}: cumulative incidence curves of every cause with confidence bands and
    the number at risk; every element can be switched off. For \code{cox()} predicted survival curves.}
  \item{\code{formula}, \code{terms}, \code{model.matrix}, \code{getCall}, \code{update}}{as for \code{lm()};
    \code{update(fit, link = "PO")} refits with a changed argument.}
  \item{\code{print}, \code{summary(x, times)}, \code{plot}, \code{as.data.frame} for \code{km()} and \code{aj()}}{
    median survival with its confidence limits, survival or cumulative incidence with confidence
    intervals and numbers at risk at chosen times, tests between groups, curves with confidence bands.}
  \item{\code{print}, \code{summary}, \code{plot}, \code{as.data.frame} for \code{compare_models()} and \code{cv_compare()}}{
    models ranked by the Brier score at every time horizon, with a short explanation of each measure.}
  \item{\code{plot(x, measure = "ipa", times = NULL, highlight = TRUE, reference = TRUE, col, main, cex = 0.85)} for a comparison}{
    one panel per time horizon with the models in the same order in every panel (best on top, by the
    mean rank over the horizons); the best model at each horizon is drawn in orange
    (\code{highlight = FALSE} switches this off) and, for the IPA, a dashed line marks no gain over the
    null model (\code{reference = FALSE} switches it off). \code{measure} is one of \code{"ipa"},
    \code{"brier"}, \code{"auc"}, \code{"calib_error"}.}
  \item{\code{cex}}{every plot (\code{plot()} of \code{km()}, of a comparison, and \code{cif_plot()}) has a
    \code{cex} argument for the size of all text, default 0.85; titles, captions and axis titles are
    wrapped to the width of the device so that they never run off the page.}
}
}
\references{
Akaike H (1974) A new look at the statistical model identification. IEEE Transactions on Automatic
  Control 19:716-723 (\code{AIC}).

Schwarz G (1978) Estimating the dimension of a model. Annals of Statistics 6:461-464 (\code{BIC}).

Wilks SS (1938) The large-sample distribution of the likelihood ratio for testing composite
  hypotheses. Annals of Mathematical Statistics 9:60-62 (likelihood ratio tests).

Casella G, Berger RL (2002) Statistical Inference, 2nd ed. Duxbury, Section 5.5 (delta method).

Efron B, Tibshirani RJ (1993) An Introduction to the Bootstrap. Chapman and Hall (bootstrap bands
  of \code{plot()} for Fine-Gray and cause-specific Cox models).

Levenberg K (1944) Quarterly of Applied Mathematics 2:164-168; Marquardt DW (1963) Journal of the
  Society for Industrial and Applied Mathematics 11:431-441 (optimiser of the parametric models).
}
\seealso{\code{\link{sx_references}} for the references of the method behind any fitted object.}
\examples{
set.seed(5)
n <- 300
x <- rbinom(n, 1, 0.5)
t1 <- rexp(n, 0.4 * exp(0.6 * x)); t2 <- rexp(n, 0.3); cens <- runif(n, 0, 4)
d <- data.frame(time = pmin(t1, t2, cens), x = x)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
fit <- directcif(Event(time, status) ~ x, d)
fit
summary(fit)
confint(fit)
estimates(fit)
fit_po <- update(fit, link = "PO")
AIC(fit, fit_po)
plot(fit, newdata = data.frame(x = 0:1), labels = c("x = 0", "x = 1"))
}
