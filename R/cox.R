#' Cox proportional hazards regression
#'
#' @description
#' Fits the Cox proportional hazards model: the effect of covariates on the \emph{hazard} (the
#' instantaneous rate of the event) is estimated without any assumption about the shape of the
#' baseline hazard. Ties, strata, delayed entry, time-varying covariates, case weights and correlated
#' rows (robust variance) are supported.
#'
#' \strong{Use it when} there is one type of event (status 0/1) and the question is how covariates
#' change the rate at which the event happens.
#' \strong{Use something else when}
#' \itemize{
#'   \item there are competing events (death from another cause, for example): \code{\link{csc}} for
#'     cause-specific hazards, \code{\link{fg}} or \code{\link{directcif}} for the cumulative
#'     incidence (the probability of the event by time t);
#'   \item a parametric baseline is wanted, to extrapolate beyond follow-up or to report
#'     time ratios: \code{\link{parsurv}};
#'   \item part of the population will never have the event: \code{\link{cure}};
#'   \item there are no covariates, only groups to compare: \code{\link{km}} and
#'     \code{\link{logrank_test}}.
#' }
#' @param formula a formula with \code{\link{Event}()} on the left and the covariates on the right.
#'   The forms of the left-hand side:
#'   \tabular{ll}{
#'   \code{Event(time, status)} \tab right-censored data, one row per subject \cr
#'   \code{Event(time, status, entry = entry)} \tab delayed entry (left truncation): the subject is
#'   observed from \code{entry} on \cr
#'   \code{Event(stop, status, entry = start)} \tab (start, stop] rows, for covariates that change
#'   over follow-up \cr
#'   }
#'   \code{status} must be 0 (censored) or 1 (event); \code{time} must be positive and larger than
#'   \code{entry}. On the right-hand side, numeric covariates, factors, interactions (\code{a * b}),
#'   transformations (\code{log(x)}, \code{I(age / 10)}) and splines (\code{splines::ns(x, 3)}) work
#'   as in \code{\link[stats]{lm}}. Two special terms are recognised:
#'   \itemize{
#'   \item \code{strata(v1, v2, ...)}: a separate baseline hazard for every combination of the
#'   variables, with common covariate effects. Use it for a variable that violates proportional
#'   hazards and whose effect you do not need to estimate (a centre, for example).
#'   \item \code{cluster(id)}: rows with the same \code{id} are correlated (several rows per
#'   subject in (start, stop] data, family members, recurrent events); the robust variance is
#'   then used.
#'   }
#' @param data a data frame containing the variables of the formula. Rows with a missing value in
#'   any of these variables are removed (the number of rows used is shown by \code{nobs()} and in the
#'   print).
#' @param ties how tied event times are handled. \code{"efron"} (default, as in
#'   \code{survival::coxph}) is more accurate when there are many ties; \code{"breslow"} is the
#'   classical approximation and gives identical results when there are no ties.
#' @param weights optional non-negative case weights, one per row of \code{data}. Integer weights act
#'   as frequency weights (a row with weight 2 counts as two identical subjects). For sampling or
#'   inverse-probability weights also set \code{robust = TRUE}: the model-based standard errors then
#'   are not valid.
#' @param robust \code{TRUE} gives the robust (sandwich) variance of Lin and Wei (1989); it changes
#'   the standard errors, intervals and tests, never the estimates. It is used automatically when
#'   \code{cluster()} is in the formula; without \code{cluster()}, every row is its own cluster.
#' @param max_iter maximum number of Newton-Raphson iterations (default 50; a well-posed model
#'   usually needs fewer than 10).
#' @param tol convergence tolerance: iteration stops when the relative change of the log partial
#'   likelihood is below \code{tol}.
#' @details
#' \subsection{The model}{
#' The hazard of subject \eqn{i} at time \eqn{t} is
#' \deqn{h_i(t) = h_0(t) \exp(\beta_1 x_{i1} + \dots + \beta_p x_{ip}),}
#' where the baseline hazard \eqn{h_0(t)} is left unspecified (Cox 1972). The coefficients are
#' estimated by maximising the partial likelihood, with the Efron (1977) or Breslow (1974) handling
#' of tied event times, by Newton-Raphson with step halving. With \code{strata()}, each stratum
#' \eqn{s} has its own \eqn{h_{0s}(t)} and the partial likelihood is the product over strata
#' (Kalbfleisch and Prentice 2002, Section 4.4). With \code{entry}, row \eqn{i} is in the risk set at
#' time \eqn{t} when \eqn{entry_i < t \le time_i}: the counting-process formulation of Andersen and
#' Gill (1982), which covers left truncation and time-varying covariates.
#' }
#' \subsection{Standard errors}{
#' Model-based: the inverse of the observed information matrix. Robust: \eqn{V (\sum_c U_c
#' U_c^\top) V}, with \eqn{V} the inverse information and \eqn{U_c} the score residuals summed over
#' the rows of cluster \eqn{c} (Lin and Wei 1989), in the Efron or Breslow form of the fit as in
#' \pkg{survival}. When the robust variance is used, the model-based one stays available as
#' \code{fit@vcov_model}.
#' }
#' \subsection{Baseline hazard and predictions}{
#' The cumulative baseline hazard is the Breslow (1972) estimator (in Efron form for Efron ties),
#' computed within each stratum. It refers to a subject with \emph{all covariates equal to 0}
#' (numeric covariates are not centred), so \code{predict(fit, type = "lp")} is \eqn{x^\top \hat\beta}
#' itself and predicted survival is \eqn{\exp\{-\hat H_0(t) e^{x^\top\hat\beta}\}}.
#' }
#' @section Interpreting the coefficients:
#' \itemize{
#'   \item The printed effect is the \strong{hazard ratio} \eqn{HR = e^{\beta}}: the factor by which
#'     the event rate is multiplied, at every time point, for a one-unit increase of a numeric
#'     covariate, or for a factor level compared with the reference level, \emph{with the other
#'     covariates held fixed}.
#'   \item HR = 1.50: a 50\% higher event rate; HR = 0.70: a 30\% lower event rate. HR < 1 means the
#'     event happens later on average, but the HR is \emph{not} a ratio of probabilities or of
#'     survival times. For absolute effects use \code{predict(fit, newdata, times, type =
#'     "survival")}, \code{plot(fit, newdata)} or \code{\link{rmst}}.
#'   \item The unit matters: the HR of age per year is often close to 1. Enter
#'     \code{I(age / 10)} for an HR per 10 years; the HR for \eqn{k} units is \eqn{HR^k}.
#'   \item With an interaction \code{a * b}, the coefficient of \code{a} is its effect when
#'     \code{b = 0} (or at the reference level of \code{b}); the HR of \code{a} at another value of
#'     \code{b} is \eqn{\exp(\beta_a + \beta_{ab} b)}.
#'   \item A variable in \code{strata()} gets no coefficient: its effect is absorbed into the
#'     baseline hazards.
#'   \item Proportional hazards means the HR is the same at every time. If \code{\link{ph_test}}
#'     rejects this, the printed HR is a weighted average of a time-varying effect; consider
#'     \code{strata()} for that variable, a split of follow-up, or a time-varying coefficient.
#' }
#' A sentence for a report: \emph{"In a Cox model adjusted for age and stage, treatment was
#' associated with a 30\% lower hazard of death (HR 0.70, 95\% CI 0.55 to 0.89, p = 0.004)."}
#' The print gives the same reading for every covariate on its \code{Reading:} lines.
#' @section Assumptions and how to check them:
#' \describe{
#'   \item{Proportional hazards}{\code{ph_test(fit)} (Grambsch and Therneau 1994) tests every
#'     covariate and the model globally; \code{plot(ph_test(fit))} shows the estimated time-varying
#'     effect \eqn{\beta(t)}, which should be flat.}
#'   \item{Functional form of numeric covariates}{plot \code{residuals(fit, "martingale")} of a model
#'     without the covariate against that covariate; a curve suggests a transformation or a spline
#'     (Therneau, Grambsch and Fleming 1990).}
#'   \item{Outliers and fit}{\code{residuals(fit, "deviance")} (roughly symmetric around 0);
#'     \code{residuals(fit, "coxsnell")} should behave like a unit exponential sample.}
#'   \item{Independent censoring}{cannot be tested from the data: censoring must be unrelated to the
#'     risk of the event given the covariates.}
#'   \item{Enough events}{as a rule of thumb, at least about 10 events per estimated coefficient.}
#' }
#' @section Messages and common problems:
#' \describe{
#'   \item{\code{WARNING: the fit did NOT converge}}{shown in the print with the reason. Usual
#'     causes: a covariate that perfectly separates events from non-events (its HR tends to 0 or
#'     infinity; check a cross-table of the covariate against status), collinear covariates
#'     (\code{information matrix not positive definite}), or covariates on a very large scale
#'     (rescale them). Do not use the standard errors of such a fit.}
#'   \item{\code{cox() needs a 0/1 status}}{the data have several event types; use
#'     \code{\link{csc}} or \code{\link{fg}}.}
#'   \item{\code{entry must be smaller than time for every subject}}{rows with zero or negative
#'     length; remove them or correct the times.}
#'   \item{\code{newdata has strata that are not in the data}}{\code{predict()} needs strata levels
#'     that were present in the fit.}
#' }
#' @section Validation:
#' Agreement with \code{survival::coxph} (\pkg{survival} 3.8) in the validation suite of the package
#' (\url{https://github.com/stats9/SurvX}, directory \code{validation}), over simulated data sets with
#' and without ties, strata, (start, stop] rows, clusters and case weights, and the real data sets
#' \code{lung}, \code{veteran} and \code{colon}: coefficients to \eqn{2 \times 10^{-9}}, model-based
#' and robust variances to \eqn{2 \times 10^{-11}}, log partial likelihood to
#' \eqn{2 \times 10^{-10}}, predicted survival (\code{survfit.coxph}) to \eqn{10^{-9}}, martingale,
#' deviance and scaled Schoenfeld residuals to \eqn{10^{-8}} (largest absolute differences).
#' @return
#' An S7 object of class \code{CoxFit}. Its main slots (read with \code{@}): \code{coefficients},
#' \code{vcov} (robust when requested), \code{vcov_model}, \code{loglik}, \code{loglik_null} (no
#' covariates), \code{n} (rows used), \code{nevent}, \code{iterations}, \code{converged},
#' \code{ties}, \code{robust}, \code{strata_levels}, \code{call}.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(fit)}, \code{summary(fit)} \tab table of HRs with confidence intervals and p-values,
#'     \code{Reading:} lines, likelihood ratio test, AIC; \code{summary} adds the null likelihood and
#'     guidance \cr
#'   \code{\link{estimates}(fit)} \tab the same table as a data frame (for reports) \cr
#'   \code{coef}, \code{vcov}, \code{confint}, \code{logLik}, \code{AIC}, \code{nobs}
#'     \tab the usual extractors (log-scale coefficients) \cr
#'   \code{predict(fit, newdata, times, type)} \tab \code{"lp"}: linear predictor;
#'     \code{"survival"}: survival probabilities at \code{times} \cr
#'   \code{plot(fit, newdata)} \tab predicted survival curves; without \code{newdata}, the average
#'     curve over the subjects of the data \cr
#'   \code{residuals(fit, type)} \tab martingale, deviance, Cox-Snell, Schoenfeld, scaled Schoenfeld
#'     (see \code{\link{residuals.SurvX}}) \cr
#'   \code{\link{ph_test}(fit)} \tab test of proportional hazards \cr
#'   \code{\link{compare_models}}, \code{\link{brier_score}}, \code{\link{cv_compare}} \tab
#'     predictive performance \cr
#' }
#' @references
#' \emph{Model and partial likelihood:}
#' Cox DR (1972). Regression models and life-tables. Journal of the Royal Statistical Society,
#'   Series B 34:187-220.
#'
#' \emph{Tied event times:}
#' Breslow NE (1974). Covariance analysis of censored survival data. Biometrics 30:89-99.
#' Efron B (1977). The efficiency of Cox's likelihood function for censored data. Journal of the
#'   American Statistical Association 72:557-565.
#'
#' \emph{Baseline hazard:}
#' Breslow NE (1972). Discussion of the paper by D. R. Cox. Journal of the Royal Statistical
#'   Society, Series B 34:216-217.
#'
#' \emph{Delayed entry and time-varying covariates (counting-process form):}
#' Andersen PK, Gill RD (1982). Cox's regression model for counting processes: a large sample study.
#'   Annals of Statistics 10:1100-1120.
#'
#' \emph{Strata:}
#' Kalbfleisch JD, Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd ed. Wiley,
#'   Section 4.4.
#'
#' \emph{Robust variance:}
#' Lin DY, Wei LJ (1989). The robust inference for the Cox proportional hazards model. Journal of
#'   the American Statistical Association 84:1074-1078.
#'
#' \emph{Likelihood ratio test and AIC in the print:}
#' Wilks SS (1938). The large-sample distribution of the likelihood ratio for testing composite
#'   hypotheses. Annals of Mathematical Statistics 9:60-62.
#' Akaike H (1974). A new look at the statistical model identification. IEEE Transactions on
#'   Automatic Control 19:716-723.
#'
#' \emph{Residuals and model checking:}
#' Therneau TM, Grambsch PM, Fleming TR (1990). Martingale-based residuals for survival models.
#'   Biometrika 77:147-160.
#' Grambsch PM, Therneau TM (1994). Proportional hazards tests and diagnostics based on weighted
#'   residuals. Biometrika 81:515-526.
#' Therneau TM, Grambsch PM (2000). Modeling Survival Data: Extending the Cox Model. Springer,
#'   Chapters 3, 6 and 7.
#'
#' The references of the methods used in a particular fit: \code{sx_references(fit)}; all
#' references behind \code{cox()}: \code{sx_references("cox")}.
#' @seealso
#' \code{\link{Event}} for the outcome; \code{\link{ph_test}}, \code{\link{residuals.SurvX}} for
#' model checking; \code{\link{predict.SurvX}}; \code{\link{csc}} and \code{\link{fg}} for competing
#' risks; \code{\link{parsurv}} for parametric models; \code{\link{sx_references}}.
#' @examples
#' set.seed(3)
#' n <- 300
#' d <- data.frame(treat = factor(rbinom(n, 1, 0.5), labels = c("control", "new")),
#'                 age = round(rnorm(n, 60, 10)),
#'                 centre = sample(c("A", "B", "C"), n, TRUE))
#' rate <- 0.1 * exp(-0.4 * (d$treat == "new") + 0.03 * (d$age - 60)) * c(A = 1, B = 2, C = 0.5)[d$centre]
#' t <- rexp(n, rate); cens <- runif(n, 0, 15)
#' d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
#'
#' ## 1. fit and read the hazard ratios (age per 10 years)
#' fit <- cox(Event(time, status) ~ treat + I(age / 10) + strata(centre), d)
#' fit
#' estimates(fit)                     # the table as a data frame
#'
#' ## 2. check proportional hazards
#' ph_test(fit)
#'
#' ## 3. absolute risks: survival at 5 and 10 years for two patients of centre A
#' nd <- data.frame(treat = c("control", "new"), age = 60, centre = "A")
#' predict(fit, nd, times = c(5, 10), type = "survival")
#' plot(fit, nd, labels = c("control", "new"))
#'
#' ## 4. delayed entry: subjects observed from 'entry' on
#' d$entry <- runif(n, 0, 0.3) * d$time
#' cox(Event(time, status, entry = entry) ~ treat + I(age / 10), d)
#'
#' ## 5. correlated rows: robust standard errors (two rows per family here)
#' d$family <- rep(seq_len(n / 2), each = 2)
#' cox(Event(time, status) ~ treat + I(age / 10) + cluster(family), d)
#' @export
cox <- function(formula, data, ties = c("efron", "breslow"), weights = NULL, robust = FALSE, max_iter = 50L, tol = 1e-10) {
  cl <- match.call(); ties <- match.arg(ties)
  f <- .sx_frame(formula, data, entry_ok = TRUE, what = "cox()", specials_ok = TRUE)
  if (any(f$status > 1L)) stop("cox() needs a 0/1 status; use csc() for cause-specific Cox models", call. = FALSE)
  n <- length(f$time); p <- ncol(f$X)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  if (!is.null(weights)) {
    omit <- attr(f$mf, "na.action")
    if (!is.null(omit)) weights <- weights[-omit]
    if (length(weights) != n) stop("weights must have one value per row of data", call. = FALSE)
  }
  robust <- isTRUE(robust) || !is.null(f$cluster)
  r <- .cox_fit(f$time, f$status, f$X, weights, ties, max_iter, tol, f$entry, f$strata)
  vm <- r$vcov
  if (robust) r$vcov <- .cox_robust(r, f$time, f$status, f$X, weights, ties, f$entry, f$strata, f$cluster)
  CoxFit(call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
         coefficients = r$coefficients, vcov = r$vcov, loglik = r$loglik, n = n, nevent = r$nevent,
         iterations = r$iterations, converged = r$converged, flags = r$flags,
         time = f$time, status = f$status, X = f$X, ties = ties, weights = weights, loglik_null = r$loglik_null,
         entry = f$entry, strata = f$strata, strata_levels = f$strata_levels, strata_vars = f$strata_vars,
         cluster = f$cluster, robust = robust, vcov_model = vm)
}

# arguments shared by the C routines of the extended Cox model
.cox_args <- function(time, status, X, weights, entry, strata) {
  n <- length(time)
  list(as.integer(n), as.integer(ncol(X)), as.integer(!is.null(entry)), if (is.null(entry)) double(n) else as.double(entry),
       as.double(time), as.integer(status), as.integer(!is.null(strata)), if (is.null(strata)) integer(n) else as.integer(strata),
       .rm(X), as.integer(!is.null(weights)), if (is.null(weights)) double(n) else as.double(weights))
}

.cox_fit <- function(time, status, X, weights, ties, max_iter = 50L, tol = 1e-10, entry = NULL, strata = NULL) {
  p <- ncol(X)
  r <- do.call(.C, c(list("SX_cox2"), .cox_args(time, status, X, weights, entry, strata),
                     list(as.integer(ties == "efron"), as.integer(max_iter), as.double(tol),
                          beta = double(p), var = double(p * p), info_d = double(3), info_i = integer(4), rc = integer(1),
                          NAOK = TRUE, PACKAGE = "SurvX")))
  .sx_rc(r$rc, "cox")
  nm <- colnames(X)
  list(coefficients = stats::setNames(r$beta, nm), vcov = matrix(r$var, p, p, dimnames = list(nm, nm)),
       loglik = r$info_d[2L], loglik_null = r$info_d[1L], nevent = r$info_i[4L], iterations = r$info_i[1L],
       converged = as.logical(r$info_i[2L]), flags = r$info_i[3L])
}

# score residuals (n x p), Efron or Breslow form; their sum is the score
.cox_score <- function(beta, time, status, X, weights, ties, entry, strata) {
  n <- length(time); p <- ncol(X)
  r <- do.call(.C, c(list("SX_cox_score"), .cox_args(time, status, X, weights, entry, strata),
                     list(as.integer(ties == "efron"), as.double(beta), U = double(n * p), rc = integer(1),
                          NAOK = TRUE, PACKAGE = "SurvX")))
  .sx_rc(r$rc, "score residuals")
  matrix(r$U, n, p, byrow = TRUE, dimnames = list(NULL, colnames(X)))
}

# robust variance of Lin & Wei (1989): V (sum_c U_c U_c') V, with U_c the score residuals summed
# within cluster c (each row its own cluster without cluster())
.cox_robust <- function(fit, time, status, X, weights, ties, entry, strata, cluster) {
  U <- .cox_score(fit$coefficients, time, status, X, weights, ties, entry, strata)
  if (!is.null(cluster)) U <- rowsum(U, cluster, reorder = FALSE)
  V <- fit$vcov
  if (!all(is.finite(V))) return(V)
  out <- V %*% crossprod(U) %*% V
  dimnames(out) <- dimnames(V)
  out
}

# Breslow/Efron baseline cumulative hazard, per stratum (column stratum: 0-based code)
.cox_baseline_raw <- function(beta, time, status, X, weights, ties, entry, strata) {
  n <- length(time)
  r <- do.call(.C, c(list("SX_cox2_baseline"), .cox_args(time, status, X, weights, entry, strata),
                     list(as.integer(ties == "efron"), as.double(beta), time = double(n), cumhaz = double(n),
                          stratum = integer(n), m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")))
  .sx_rc(r$rc, "baseline hazard")
  k <- seq_len(r$m)
  data.frame(time = r$time[k], cumhaz = r$cumhaz[k], stratum = r$stratum[k])
}

.cox_baseline <- function(object) {
  .cox_baseline_raw(object@coefficients, object@time, object@status, object@X, object@weights, object@ties,
                    object@entry, object@strata)
}

# cumulative hazard of baseline bh at times tt for rows in strata s (0-based; NULL: one stratum)
.cumhaz_at <- function(bh, tt, s = NULL) {
  if (is.null(s)) return(.step_at(bh$time, bh$cumhaz, tt, 0))
  out <- numeric(length(tt))
  for (k in unique(s)) {
    i <- s == k; b <- bh[bh$stratum == k, , drop = FALSE]
    out[i] <- if (nrow(b)) .step_at(b$time, b$cumhaz, tt[i], 0) else 0
  }
  out
}

#' Predictions from SurvX models
#'
#' \code{predict()} methods of the model classes.
#'
#' @name predict.SurvX
#' @param object a fitted model.
#' @param newdata data frame with the covariates; the data used for fitting when \code{NULL}.
#' @param times time points. For \code{directcif} fits \code{Inf} gives the plateau.
#' @param type for \code{cox}: \code{"lp"} (linear predictor) or \code{"survival"}.
#' @param cause cause of interest (\code{csc} and \code{directcif} fits).
#' @param se for \code{directcif} fits: also return delta-method standard errors.
#' @param ... not used.
#' @return a matrix with one row per row of \code{newdata} and one column per time point
#'   (a vector for \code{type = "lp"}); with \code{se = TRUE} a list with \code{cif} and \code{se}.
NULL

S7::method(predict, CoxFit) <- function(object, newdata = NULL, times = NULL, type = c("lp", "survival"), ...) {
  type <- match.arg(type)
  lp <- drop(.sx_newX(object, newdata) %*% object@coefficients)
  if (type == "lp") return(lp)
  if (is.null(times)) stop("times are needed for type = \"survival\"", call. = FALSE)
  bh <- .cox_baseline(object)
  sn <- .sx_newstrata(object, newdata)
  out <- if (is.null(sn)) exp(-outer(exp(lp), .step_at(bh$time, bh$cumhaz, times, 0)))
         else t(vapply(seq_along(lp), function(i) exp(-exp(lp[i]) * .cumhaz_at(bh, times, rep(sn[i], length(times)))),
                       numeric(length(times))))
  if (length(times) == 1L && !is.null(sn)) out <- matrix(out, ncol = 1L)
  dimnames(out) <- list(NULL, paste0("t=", times))
  out
}
