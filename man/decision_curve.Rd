\name{decision_curve}
\alias{decision_curve}
\alias{print.DecisionCurve}
\alias{plot.DecisionCurve}
\alias{as.data.frame.DecisionCurve}
\title{Decision curve analysis}
\description{
Net benefit of deciding by predicted risks at a time horizon, against treating everybody or
nobody, over a range of risk thresholds (Vickers and Elkin 2006), with right censoring and
competing risks (Vickers et al. 2008).

\strong{Use it} to judge whether using a model to decide (treat, screen, refer) would do more
good than harm, at the thresholds that matter for the decision.
}
\usage{
decision_curve(risks, time = NULL, status = NULL, data = NULL, t, cause = 1L,
               thresholds = seq(0, 0.5, by = 0.01))

\method{print}{DecisionCurve}(x, ...)

\method{plot}{DecisionCurve}(x, ylim = NULL, col = NULL, lwd = 2, xlab = "risk threshold",
     ylab = "net benefit", main = NULL, legend = TRUE, ...)

\method{as.data.frame}{DecisionCurve}(x, ...)
}
\arguments{
\item{risks}{a named list of predicted risks of the cause by \code{t} (numeric vectors, one per
  subject), or of fitted models of the package (their predictions for \code{data} at \code{t}).
  A single numeric vector or model is allowed.}
\item{time,status}{follow-up time and status (0 = censored, 1..K = cause); with fitted models and
  \code{data} they may be omitted (taken from the formula of the first model).}
\item{data}{data for the predictions of fitted models.}
\item{t}{the time horizon.}
\item{cause}{cause of interest.}
\item{thresholds}{risk thresholds \eqn{p} in \eqn{[0, 1)}.}
\item{x}{the result of \code{decision_curve()}.}
\item{ylim,col,lwd,xlab,ylab,main}{graphical parameters; \code{ylim} defaults to a range from a
  little below 0 to the largest net benefit.}
\item{legend}{draw a legend.}
\item{...}{further graphical parameters (\code{plot}); not used otherwise.}
}
\details{
Treating the subjects whose predicted risk is at least \eqn{p} has net benefit
\eqn{NB(p) = \pi_p F_p(t) - \pi_p \{1 - F_p(t)\} p/(1-p)}, with \eqn{\pi_p} the share of subjects
above the threshold and \eqn{F_p(t)} the Kaplan-Meier (Aalen-Johansen) risk of the cause by
\eqn{t} among them: true positives minus false positives per subject, the false ones weighted by
the odds of the threshold. Treating everybody has \eqn{F(t) - \{1 - F(t)\} p/(1-p)}, treating
nobody 0.
}
\section{Interpreting the result}{
A model is useful at threshold \eqn{p} when its net benefit is above both "treat all" and
"treat none". A net benefit of 0.05 means the same as finding 5 more true cases per 100 subjects
with no false positives. Look only at the thresholds a decision maker would use.
}
\section{Validation}{
Equals the definition computed with \code{survival::survfit} and \code{dcurves::dca} with a
\code{Surv} outcome and the same thresholds (validation suite, \code{validation/R/07_assess.R}).
}
\value{
an object of class \code{DecisionCurve}: a list with \code{table} (columns \code{model},
  \code{threshold}, \code{net_benefit}, including the rows of "treat all" and "treat none"),
  \code{t}, \code{cause} and \code{n}.
}
\references{
Vickers AJ, Elkin EB (2006) Decision curve analysis: a novel method for evaluating prediction
  models. Medical Decision Making 26:565-574.

Vickers AJ, Cronin AM, Elkin EB, Gonen M (2008) Extensions to decision curve analysis, a novel
  method for evaluating diagnostic tests, prediction models and molecular markers. BMC Medical
  Informatics and Decision Making 8:53.

All references: \code{sx_references("decision_curve")}.
}
\seealso{
\code{\link{brier_score}}, \code{\link{compare_models}}, \code{vignette("prediction", "SurvX")}.
}
\examples{
set.seed(3)
n <- 400
d <- data.frame(x = rnorm(n))
t1 <- rexp(n, 0.3 * exp(0.8 * d$x)); t2 <- rexp(n, 0.2); cens <- runif(n, 0, 5)
d$time <- pmin(t1, t2, cens)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
fit <- fg(Event(time, status) ~ x, d)
dc <- decision_curve(list("Fine-Gray" = fit), data = d, t = 2, thresholds = seq(0, 0.6, 0.02))
dc
plot(dc)
}
