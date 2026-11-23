\name{cindex}
\alias{cindex}
\alias{print.SurvXConcordance}
\title{Concordance index}
\description{
Harrell's concordance index C of a risk score and Uno's version with inverse probability of
censoring weights; with competing risks, the C index of one cause of Wolbers et al. (2009) and its
weighted version (Wolbers et al. 2014). One function adapts to both kinds of data: without
competing events the competing-risks versions are exactly Harrell's and Uno's. Standard errors
from the infinitesimal jackknife.

\strong{Use it} to describe how well a score ranks subjects by the time to the event over the
whole follow-up. \strong{Do not use it} to compare models that predict the risk by a fixed time
\eqn{t}: the C index is not a proper measure for \eqn{t}-year risks (Blanche, Kattan and Gerds
2019), so a wrong model can have a higher C. For that use the AUC at \eqn{t} and the Brier score
(\code{\link{auc_t}}, \code{\link{brier_score}}, \code{\link{compare_models}}), which is why the
comparison tables of the package do not report C.
}
\usage{
cindex(object, time = NULL, status = NULL, newdata = NULL, cause = 1L,
       method = c("harrell", "uno"), tau = Inf, times = NULL, level = 0.95)

\method{print}{SurvXConcordance}(x, digits = 4L, ...)
}
\arguments{
\item{object}{a fitted model of the package, or a numeric risk score (higher = earlier event),
  one per subject of \code{time} and \code{status}.}
\item{time,status}{follow-up time and status (0 = censored, 1..K = cause) when \code{object} is a
  score.}
\item{newdata}{for a fitted model: data on which to compute C (default: the data of the fit; must
  contain the outcome variables).}
\item{cause}{cause of interest.}
\item{method}{\code{"harrell"} (all comparable pairs have weight 1) or \code{"uno"} (inverse
  probability of censoring weights, \eqn{\hat G} the Kaplan-Meier estimate of the censoring
  distribution: \eqn{1/\hat G(T_i-)^2} for a pair whose control is still event-free at the case's
  time \eqn{T_i}, and \eqn{1/\{\hat G(T_i-)\hat G(T_j-)\}} for a control that failed from a
  competing cause at \eqn{T_j \le T_i}).}
\item{tau}{only events at times \eqn{\le} \code{tau} define comparable pairs (Uno et al. 2011
  recommend a \code{tau} where follow-up is still adequate).}
\item{times}{for a model without a linear predictor for the cause (Fine-Gray, cause-specific,
  direct, parametric and cure models), the risk score is the predicted risk at this one time;
  default: the median follow-up time. A Cox model uses its linear predictor.}
\item{level}{confidence level.}
\item{x}{the result of \code{cindex()}.}
\item{digits}{number of decimals printed.}
\item{...}{not used.}
}
\details{
A pair \eqn{(i, j)} is comparable when \eqn{i} has the event of interest at \eqn{T_i \le \tau}
and \eqn{j} is known to be event-free at \eqn{T_i}: \eqn{T_j > T_i}, or \eqn{T_j = T_i} with
\eqn{j} censored, or (competing risks) \eqn{j} failed from another cause at any time. Events of
interest at the same time are not compared. The pair is concordant when \eqn{i} has the higher
score; tied scores count 1/2. \eqn{C} is the weighted share of concordant pairs.

The standard error is \eqn{(\sum_u d_u^2)^{1/2}} with \eqn{d_u = (N_u - C D_u)/D}, where
\eqn{N_u} and \eqn{D_u} are the weighted concordance and the weight of the comparable pairs that
contain subject \eqn{u} and \eqn{D} the total weight: the infinitesimal jackknife of
\code{survival::concordance}. With Uno's weights the estimation of \eqn{\hat G} is ignored, as
there. The computation sorts the subjects once and uses binary indexed trees, so it takes
\eqn{O(n \log n)} time.

The summaries of Cox, cause-specific Cox and parametric survival models show Harrell's C of the
fitted model on its own data (for cause-specific models, of each cause with the other causes as
censored).
}
\section{Interpreting the result}{
0.5 is no better than chance, 1 is perfect ranking; 0.6 to 0.75 is common for clinical scores.
Harrell's C depends on the censoring distribution when the score is not perfectly
proportional; Uno's version removes that dependence up to \code{tau}.
}
\section{Validation}{
Harrell's C and its standard error equal \code{survival::concordance(..., reverse = TRUE)} and,
with \code{method = "uno"}, \code{timewt = "n/G2"} with \code{ymax = tau} (validation suite,
\code{validation/R/07_assess.R}); with competing risks Uno's version equals \code{pec::cindex}. The C++ code
is tested against a direct computation over all pairs.
}
\value{
an object of class \code{SurvXConcordance}: a list with \code{concordance}, \code{se},
  \code{lower}, \code{upper}, \code{counts} (weighted numbers of concordant, discordant and tied
  pairs and the number of comparable pairs), \code{dfbeta}, \code{method}, \code{tau},
  \code{cause}, \code{level}, \code{competing} and \code{n}.
}
\references{
Harrell FE, Califf RM, Pryor DB, Lee KL, Rosati RA (1982) Evaluating the yield of medical tests.
  JAMA 247:2543-2546.

Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011) On the C-statistics for evaluating overall
  adequacy of risk prediction procedures with censored survival data. Statistics in Medicine
  30:1105-1117.

Wolbers M, Koller MT, Witteman JCM, Steyerberg EW (2009) Prognostic models with competing risks:
  methods and application to coronary risk prediction. Epidemiology 20:555-561.

Wolbers M, Blanche P, Koller MT, Witteman JCM, Gerds TA (2014) Concordance for prognostic models
  with competing risks. Biostatistics 15:526-539.

Blanche P, Kattan MW, Gerds TA (2019) The c-index is not proper for the evaluation of t-year
  predicted risks. Biostatistics 20:347-357.

Therneau T, Atkinson E (2024) Concordance. Vignette of the survival package.

All references: \code{sx_references("cindex")}.
}
\seealso{
\code{\link{auc_t}}, \code{\link{compare_models}}, \code{vignette("prediction", "SurvX")}.
}
\examples{
set.seed(3)
n <- 300
d <- data.frame(x = rnorm(n))
t1 <- rexp(n, 0.3 * exp(0.8 * d$x)); t2 <- rexp(n, 0.2); cens <- runif(n, 0, 5)
d$time <- pmin(t1, t2, cens)
d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
cindex(d$x, d$time, d$status)                  # Wolbers' C for cause 1
cindex(d$x, d$time, d$status, method = "uno", tau = 3)
cindex(fg(Event(time, status) ~ x, d), times = 2)
}
