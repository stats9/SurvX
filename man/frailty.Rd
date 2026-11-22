\name{frailty}
\alias{frailty}
\alias{FrailtyFit}
\title{Shared gamma frailty models for clustered survival data}
\description{
When subjects come in clusters (patients in hospitals, children in families, teeth in mouths,
recurrent events of one person, twins), members of a cluster share unmeasured factors and their
event times are correlated. A shared frailty model gives every cluster an unobserved factor
\eqn{z_c} that multiplies the hazards of all its members,
\deqn{h_{j}(t \mid z_c) = z_c\, h_0(t) \exp(x_j'\beta),}
with \eqn{z_c} gamma-distributed with mean 1 and variance \eqn{\theta}. Write
\code{+ frailty(id)} in the formula of \code{\link{cox}} (semi-parametric, Breslow baseline) or of
\code{\link{parsurv}} (any parametric family); this page describes both.

\strong{Use it when} you want covariate effects within clusters and a measure of how much clusters
differ (\eqn{\theta}), or predictions for a new or for an observed cluster. \strong{Use
something else when} you only need population-average effects with valid standard errors:
\code{cox(... + cluster(id))} (robust variance) answers that question without a model for the
dependence.
}
\section{Usage}{
\preformatted{
cox(Event(time, status) ~ x + frailty(id), data)
parsurv(Event(time, status) ~ x + frailty(id), data, baseline = "weibull")
}
}
\details{
\strong{Likelihood.} The frailties are integrated out: cluster \eqn{c} with \eqn{D_c} events
and \eqn{\Lambda_c = \sum_j H_0(t_j) e^{x_j'\beta}} contributes
\deqn{\prod_j h_j(t_j)^{d_j} \prod_{m=0}^{D_c - 1}(1 + m\theta)\,(1 + \theta\Lambda_c)^{-(1/\theta + D_c)}}
(Clayton 1978; Klein 1992). With a parametric baseline the marginal likelihood is maximised by
Newton-Raphson with exact first and second derivatives. With the Breslow baseline, the baseline
hazard is profiled out by the EM fixed point of Nielsen et al. (1992) and Klein (1992), and the
profile likelihood in \eqn{(\beta, \log\theta)} is maximised by Newton-Raphson; its curvature
gives standard errors of \eqn{\beta} that account for the estimation of \eqn{\theta} (Murphy and
van der Vaart 2000). Ties are handled as in Breslow's method. The estimates are those of
\code{survival::coxph(... + frailty(id, distribution = "gamma"), ties = "breslow")} and of
\code{frailtyEM::emfrail}.

\strong{Reading the results.} Hazard ratios are conditional on the frailty: they compare
subjects of the same cluster (or of clusters with equal frailty). The population-average
(marginal) hazard ratio is closer to 1 and changes over time. \eqn{\theta} is printed with its
confidence interval, the likelihood ratio test of \eqn{\theta = 0} (with the half chi-square
distribution of a parameter on the boundary; Self and Liang 1987), Kendall's \eqn{\tau =
\theta/(\theta + 2)} between the event times of two members of a cluster, and the ratio of the
hazards of a cluster at the 75th and at the 25th percentile of the frailty distribution.

\strong{Frailties of the clusters.} Given the data, \eqn{z_c} has a gamma distribution with
shape \eqn{1/\theta + D_c} and rate \eqn{1/\theta + \Lambda_c}; its mean is the empirical Bayes
estimate, shown with a 95\% posterior interval by \code{as.data.frame(fit, what = "frailty")} and
\code{plot(fit)}.

\strong{Prediction.} \code{predict(fit, newdata, times, frailty = )}: \code{"marginal"} (a new
cluster, frailty integrated out: \eqn{(1 + \theta H)^{-1/\theta}}), \code{"conditional"}
(frailty 1, the average cluster: \eqn{e^{-H}}) or \code{"cluster"} (an observed cluster, given
its data: \eqn{\{r_c/(r_c + H)\}^{s_c}}, with the cluster in a column of \code{newdata} named as
in the formula). \code{type = "survival"} or \code{"cumhaz"}.

Not available with frailty: delayed entry, weights, \code{strata()}, \code{cluster()},
\code{tt()} and interval-censored data.
}
\section{Messages}{
\describe{
  \item{\code{the frailty variance is at the boundary (0)}}{the data show no heterogeneity
    between clusters; the fit is then the model without frailty.}
}
}
\value{
An S7 object of class \code{FrailtyFit} with \code{print()}, \code{summary()}, \code{coef()}
(\eqn{\beta}; for a parametric baseline also the baseline parameters), \code{vcov()},
\code{logLik()} (marginal; the degrees of freedom include \eqn{\theta}), \code{AIC()},
\code{predict()}, \code{plot()} and \code{as.data.frame(fit, what = c("coefficients",
"frailty"))}. \code{fit@theta}, \code{fit@theta_se} and \code{fit@lrt} hold the frailty variance,
its standard error and the test.
}
\references{
Clayton DG (1978). A model for association in bivariate life tables and its application in
  epidemiological studies of familial tendency in chronic disease incidence. Biometrika 65:141-151.

Klein JP (1992). Semiparametric estimation of random effects using the Cox model based on the EM
  algorithm. Biometrics 48:795-806.

Nielsen GG, Gill RD, Andersen PK, Sorensen TIA (1992). A counting process approach to maximum
  likelihood estimation in frailty models. Scandinavian Journal of Statistics 19:25-43.

Murphy SA, van der Vaart AW (2000). On profile likelihood. Journal of the American Statistical
  Association 95:449-465.

Self SG, Liang KY (1987). Asymptotic properties of maximum likelihood estimators and likelihood
  ratio tests under nonstandard conditions. Journal of the American Statistical Association
  82:605-610.

Duchateau L, Janssen P (2008). The Frailty Model. Springer, New York.

Balan TA, Putter H (2019). frailtyEM: an R package for estimating semiparametric shared frailty
  models. Journal of Statistical Software 90(7):1-29.
}
\seealso{
\code{\link{cox}}, \code{\link{parsurv}}; \code{cox(... + cluster(id))} for robust standard errors
without a frailty.
}
\examples{
set.seed(3)
C <- 80; m <- 4
id <- rep(seq_len(C), each = m)
z <- rgamma(C, shape = 2, rate = 2)[id]          # frailty variance 0.5
x <- rnorm(C * m); trt <- rbinom(C * m, 1, 0.5)
t <- rweibull(C * m, 1.3, (0.1 * z * exp(0.5 * x - 0.7 * trt))^(-1 / 1.3))
cens <- runif(C * m, 0, 15)
d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x, trt, id)

fc <- cox(Event(time, status) ~ x + trt + frailty(id), data = d)
fc
fw <- parsurv(Event(time, status) ~ x + trt + frailty(id), data = d, baseline = "weibull")
fw
head(as.data.frame(fc, what = "frailty"))
predict(fc, data.frame(x = 0, trt = 0:1), times = c(5, 10), frailty = "marginal")
plot(fc)
}
