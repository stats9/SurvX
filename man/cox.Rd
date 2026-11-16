\name{cox}
\alias{cox}
\title{Cox proportional hazards model}
\description{
Partial-likelihood fit by Newton-Raphson with step halving. Supports strata, delayed entry and
time-varying covariates in (start, stop] form, case weights, and a robust variance for clustered or
repeated rows.
}
\usage{
cox(formula, data, ties = c("efron", "breslow"), weights = NULL, robust = FALSE,
    max_iter = 50L, tol = 1e-10)
}
\arguments{
\item{formula}{\code{Event(time, status) ~ covariates}; status 0/1. For delayed entry, or for
  (start, stop] rows with time-varying covariates, use \code{Event(stop, status, entry = start)}.
  The right-hand side may contain \code{strata(...)}, which gives every stratum its own baseline
  hazard, and \code{cluster(id)}, which requests the robust variance with rows of the same
  \code{id} treated as correlated.}
\item{data}{a data frame.}
\item{ties}{\code{"efron"} (default) or \code{"breslow"}.}
\item{weights}{optional case weights, one per row of \code{data}.}
\item{robust}{use the robust (sandwich) variance; always \code{TRUE} when \code{cluster()} is in
  the formula. Without \code{cluster()} each row is its own cluster.}
\item{max_iter,tol}{iteration limit and convergence tolerance.}
}
\details{
The partial likelihood is that of Cox (1972) with the Efron (1977) or Breslow (1974) handling of
tied event times. With \code{strata()} the partial likelihood is the product over strata, each with
its own risk sets and baseline hazard, and common coefficients (Kalbfleisch and Prentice 2002,
Section 4.4). With \code{entry}, row \eqn{i} is at risk at time \eqn{t} when
\eqn{entry_i < t \le time_i}: this is the counting-process form of Andersen and Gill (1982), which
covers delayed entry (left truncation) and covariates that change over follow-up (one row per
interval of constant covariate values).

The robust variance of Lin and Wei (1989) is \eqn{V (\sum_c U_c U_c^\top) V}, with \eqn{V} the
inverse information and \eqn{U_c} the score residuals summed over the rows of cluster \eqn{c}; the
score residuals take the Efron or Breslow form of the fit, as in \pkg{survival}. The model-based
variance stays available as \code{x@vcov_model}. The baseline hazard (Breslow estimator, Efron form
with Efron ties) is estimated within each stratum; \code{predict(type = "survival")} uses the
stratum of each row of \code{newdata}.

Results agree with \code{survival::coxph}: coefficients, variances (model-based and robust),
log partial likelihood, predicted survival and residuals, with and without strata, (start, stop]
data and clusters (see the tests of the package).
}
\value{
an S7 object of class \code{CoxFit}. Methods: \code{print}, \code{summary}, \code{coef},
  \code{vcov}, \code{confint}, \code{logLik}, \code{AIC}, \code{nobs}, \code{predict},
  \code{residuals}, \code{plot}; \code{\link{ph_test}} tests proportional hazards.
}
\references{
Cox DR (1972) Regression models and life-tables. Journal of the Royal Statistical Society, Series B
  34:187-220.

Breslow NE (1974) Covariance analysis of censored survival data. Biometrics 30:89-99.

Efron B (1977) The efficiency of Cox's likelihood function for censored data. Journal of the
  American Statistical Association 72:557-565.

Andersen PK, Gill RD (1982) Cox's regression model for counting processes: a large sample study.
  Annals of Statistics 10:1100-1120.

Lin DY, Wei LJ (1989) The robust inference for the Cox proportional hazards model. Journal of the
  American Statistical Association 84:1074-1078.

Kalbfleisch JD, Prentice RL (2002) The Statistical Analysis of Failure Time Data, 2nd ed. Wiley,
  Section 4.4 (stratification).

Therneau TM, Grambsch PM (2000) Modeling Survival Data: Extending the Cox Model. Springer,
  Chapters 3 and 7.
}
\examples{
set.seed(3)
n <- 300
d <- data.frame(x = rbinom(n, 1, 0.5), z = rnorm(n), center = sample(c("A", "B", "C"), n, TRUE),
                id = rep(1:150, each = 2))
t <- rexp(n, 0.3 * exp(0.5 * d$x - 0.2 * d$z) * c(A = 1, B = 2, C = 0.5)[d$center])
cens <- runif(n, 0, 5)
d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
cox(Event(time, status) ~ x + z + strata(center), d)
cox(Event(time, status) ~ x + z + cluster(id), d)       # robust variance
# delayed entry: subjects enter the risk set at 'entry'
d$entry <- runif(n, 0, 0.2) * d$time
cox(Event(time, status, entry = entry) ~ x + z, d)
}
