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
#'   \item \code{tt(x)}: a time-varying coefficient for the numeric covariate \code{x}, together
#'   with the argument \code{tt}: \code{x + tt(x)} gives the log hazard ratio
#'   \eqn{\beta + \gamma g(t)} at time \eqn{t}. Use it when \code{\link{ph_test}} shows that the
#'   effect of \code{x} changes over time and that effect is of interest.
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
#' @param tt the function of time \eqn{g(t)} of the \code{tt()} terms, for example
#'   \code{function(t) log(t)} (an effect that changes linearly in log time; the usual first choice),
#'   \code{function(t) t}, or \code{function(t) as.numeric(t > 2)} (a different effect after time
#'   2). A list gives one function per \code{tt()} term. Needed only with \code{tt()} terms.
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
#' \subsection{Time-varying coefficients}{
#' With \code{tt()} terms the model is \eqn{h(t \mid x) = h_0(t)\exp\{\beta x + \gamma x g(t)\}}. Every
#' row is split at the event times inside its follow-up, and the \code{tt()} column of each piece is
#' \eqn{x\,g(t)} at the end of the piece; the partial likelihood of the split data is exactly that of
#' the model, as in \code{survival::coxph(tt = )} (Therneau and Grambsch 2000, Chapter 6). With many
#' distinct event times the split data can be large (about \eqn{n} times the number of events
#' rows); rounding the times reduces it. \code{\link{tt_effect}} gives the hazard ratio at chosen
#' times. Predictions, residuals and \code{ph_test()} are not defined for such models.
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
#'     \code{strata()} for that variable, or a time-varying coefficient with \code{tt()}.
#'   \item With \code{x + tt(x)}, the row of \code{x} is the HR where \eqn{g(t) = 0} (at \eqn{t = 1}
#'     for \eqn{g = \log}) and the row of \code{tt(x)} the factor by which that HR changes per unit
#'     of \eqn{g(t)}; the print and \code{\link{tt_effect}} give the HR at chosen times, which is what
#'     to report.
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
#' (\url{https://github.com/stats9/SurvX}, directory \code{validation}), over 82 simulated and real data
#' sets (\code{lung}, \code{veteran}, \code{colon}; up to \eqn{n = 20000}) with and without ties, strata,
#' (start, stop] rows, clusters and case weights, 438 fits in all: coefficients to
#' \eqn{7 \times 10^{-7}}, model-based and robust variances to \eqn{4 \times 10^{-9}}, log partial
#' likelihood to a relative \eqn{4 \times 10^{-8}}; predicted survival (\code{survfit.coxph}),
#' martingale, deviance and scaled Schoenfeld residuals to \eqn{10^{-8}}. In a Monte Carlo study (1000
#' replications, \eqn{n} = 200, 500, 1000), with strata, delayed entry and time-varying covariates,
#' the standard errors matched the empirical standard deviation of the estimates (ratio 0.95 to 1.05)
#' and the 95\% intervals covered 93.4 to 96.1\%. A fit on \eqn{n = 10000} takes about 0.01 s
#' (\code{coxph}: 0.03 s).
#' @return
#' An S7 object of class \code{CoxFit}. Its main slots (read with \code{@}): \code{coefficients},
#' \code{vcov} (robust when requested), \code{vcov_model}, \code{loglik}, \code{loglik_null} (no
#' covariates), \code{n} (rows used), \code{nevent}, \code{iterations}, \code{converged},
#' \code{ties}, \code{robust}, \code{strata_levels}, \code{model_call} (also \code{getCall(fit)}).
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
#'   \code{\link{tt_effect}(fit, times)} \tab models with \code{tt()}: HR at chosen times \cr
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
#' rate <- 0.1 * exp(-0.4 * (d$treat == "new") + 0.03 * (d$age - 60)) *
#'   c(A = 1, B = 2, C = 0.5)[d$centre]
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
#'
#' ## 6. an effect that changes over time: tt()
#' d$new <- as.numeric(d$treat == "new")
#' ft <- cox(Event(time, status) ~ new + tt(new) + I(age / 10) + strata(centre), d,
#'           tt = function(t) log(t))
#' ft                                     # HR of the new treatment at the quartiles of time
#' tt_effect(ft, times = c(1, 5, 10))
#' @export
cox <- function(formula, data, ties = c("efron", "breslow"), weights = NULL, robust = FALSE, max_iter = 50L, tol = 1e-10,
                tt = NULL) {
  cl <- match.call(); ties <- match.arg(ties)
  f <- .sx_frame(formula, data, entry_ok = TRUE, what = "cox()", specials_ok = TRUE, tt_ok = TRUE)
  if (any(f$status > 1L)) stop("cox() needs a 0/1 status; use csc() for cause-specific Cox models", call. = FALSE)
  n <- length(f$time); p <- ncol(f$X)
  if (p < 1L) stop("at least one covariate is needed", call. = FALSE)
  if (!is.null(weights)) {
    omit <- attr(f$mf, "na.action")
    if (!is.null(omit)) weights <- weights[-omit]
    if (length(weights) != n) stop("weights must have one value per row of data", call. = FALSE)
  }
  robust <- isTRUE(robust) || !is.null(f$cluster)
  if (length(f$tt_cols)) {
    if (is.null(tt)) stop("the formula has tt() terms: give the function of time, for example tt = function(t) log(t)", call. = FALSE)
    s <- .tt_split(f, weights, tt)
    r <- .cox_fit(s$time, s$status, s$X, s$weights, ties, max_iter, tol, s$entry, s$strata)
    vm <- r$vcov
    if (robust) r$vcov <- .cox_robust(r, s$time, s$status, s$X, s$weights, ties, s$entry, s$strata, s$cluster)
  } else {
    r <- .cox_fit(f$time, f$status, f$X, weights, ties, max_iter, tol, f$entry, f$strata)
    vm <- r$vcov
    if (robust) r$vcov <- .cox_robust(r, f$time, f$status, f$X, weights, ties, f$entry, f$strata, f$cluster)
  }
  CoxFit(model_call = cl, terms = f$terms, xlevels = f$xlevels, contrasts = f$contrasts,
         coefficients = r$coefficients, vcov = r$vcov, loglik = r$loglik, n = n, nevent = r$nevent,
         iterations = r$iterations, converged = r$converged, flags = r$flags,
         time = f$time, status = f$status, X = f$X, ties = ties, weights = weights, loglik_null = r$loglik_null,
         entry = f$entry, strata = f$strata, strata_levels = f$strata_levels, strata_vars = f$strata_vars,
         cluster = f$cluster, robust = robust, vcov_model = vm,
         tt = if (length(f$tt_cols)) list(fun = tt, cols = f$tt_cols) else NULL)
}

# Time-varying coefficients: every row is split at the event times inside its (entry, time]
# interval, and the tt() columns become x * g(t) evaluated at the end of each piece. The pieces
# contain exactly one event time each (their end), so the partial likelihood of the split data is
# that of the model log h(t | x) = log h0(t) + beta'x + gamma x g(t), as in survival::coxph(tt = ).
.tt_split <- function(f, weights, tt) {
  time <- f$time; status <- f$status; n <- length(time)
  entry <- f$entry %||% rep(0, n)
  et <- sort(unique(time[status == 1L]))
  lo <- findInterval(entry, et)                      # event times <= entry
  hi <- findInterval(time, et, left.open = TRUE)     # event times < time
  k <- pmax(hi - lo, 0L)                             # cut points strictly inside (entry, time)
  rows <- k + 1L
  N <- sum(as.numeric(rows))
  if (N > 5e7) stop("tt(): splitting at every event time would create ", format(N, big.mark = ","),
                    " rows; use fewer distinct times (round them) or split follow-up into periods yourself", call. = FALSE)
  idx <- rep.int(seq_len(n), rows)
  r <- sequence(rows)
  ci <- lo[idx] + r                                  # index in et of the end of a piece that is a cut
  last <- r == rows[idx]
  stop_t <- ifelse(last, time[idx], et[pmin(ci, length(et))])
  start_t <- ifelse(r == 1L, entry[idx], et[pmax(ci - 1L, 1L)])
  X <- f$X[idx, , drop = FALSE]
  funs <- if (is.function(tt)) rep(list(tt), length(f$tt_cols)) else tt
  if (!is.list(funs) || length(funs) != length(f$tt_cols) || !all(vapply(funs, is.function, logical(1))))
    stop("tt must be a function of time, or a list with one function per tt() term", call. = FALSE)
  for (j in seq_along(f$tt_cols)) {
    g <- funs[[j]](stop_t)
    if (length(g) != length(stop_t) || any(!is.finite(g)))
      stop("tt function ", j, " must return a finite value for every time (log(t) needs t > 0)", call. = FALSE)
    X[, f$tt_cols[j]] <- X[, f$tt_cols[j]] * g
  }
  list(time = stop_t, status = ifelse(last, status[idx], 0L), X = X, entry = start_t,
       weights = if (is.null(weights)) NULL else weights[idx],
       strata = if (is.null(f$strata)) NULL else f$strata[idx],
       cluster = if (!is.null(f$cluster)) f$cluster[idx] else idx)
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
#' @description
#' Predictions from every fitted model of SurvX: linear predictors, survival probabilities,
#' cumulative incidences, hazards, quantiles (median time), restricted mean survival time and cured
#' fractions, for new subjects or for the subjects of the data, optionally with standard errors and
#' confidence limits. Absolute predictions such as these are usually the clearest way to report a
#' model: "a 60-year-old patient on the new treatment has a 24\% risk of relapse within 5 years".
#' @details
#' Cox models: \eqn{S(t \mid x) = \exp\{-\hat H_0(t) e^{x'\hat\beta}\}} with the Breslow
#' estimate of the baseline cumulative hazard. Cause-specific Cox models: the cumulative incidence
#' \eqn{\int_0^t \hat S(u-)\, d\hat H_k(u \mid x)} from all cause-specific hazards (product-limit
#' or exponential form of \eqn{\hat S}, see \code{product_limit} of \code{\link{csc}}). Fine-Gray
#' models: \eqn{1 - \exp\{-\hat\Lambda_{10}(t) e^{x'\hat\beta}\}}. Direct and parametric models
#' evaluate the fitted distribution; parametric cause-specific models integrate the hazards
#' numerically (Gauss-Legendre). Standard errors of direct and parametric predictions use the delta
#' method with the covariance matrix of the parameters; limits are computed on a transformed scale
#' (log-log for probabilities, log for hazards and times) and transformed back.
#' @section Usage:
#' \code{predict(object, newdata = NULL, times, type, cause = 1, se = FALSE, level = 0.95,
#' average = FALSE, p = 0.5)}; each model class accepts the arguments that apply to it.
#' @section Arguments:
#' \describe{
#' \item{\code{object}}{a fitted model.}
#' \item{\code{newdata}}{a data frame with one row per subject to predict for, containing every
#'   covariate of the model (and the \code{strata()} variables of a stratified Cox model). Factors may
#'   be given as character values of their levels. \code{NULL}: the subjects used for fitting.}
#' \item{\code{times}}{time points. For \code{directcif} fits \code{Inf} gives the plateau; for
#'   \code{parsurv} with \code{type = "rmst"} the upper limits of integration.}
#' \item{\code{type}}{for \code{cox}: \code{"lp"} (linear predictor) or \code{"survival"}. For
#'   \code{parsurv}: \code{"survival"}, \code{"hazard"}, \code{"cumhaz"}, \code{"quantile"} (at the
#'   probabilities \code{p}) or \code{"rmst"} (restricted mean survival time). For \code{csc} with a
#'   parametric baseline: \code{"cif"} or \code{"survival"} (event-free probability). For
#'   \code{cure}: \code{"survival"}, \code{"hazard"}, \code{"cumhaz"}, \code{"uncured"} (survival
#'   of the uncured) or \code{"cured"} (the cured fraction, without \code{times}).}
#' \item{\code{p}}{\code{parsurv}, \code{type = "quantile"}: probabilities; the time by which a
#'   fraction \code{p} has failed (\code{p = 0.5}: the median).}
#' \item{\code{cause}}{cause of interest (\code{csc} and \code{directcif} fits).}
#' \item{\code{se}}{for \code{directcif}, \code{parsurv} and parametric \code{csc} fits: also return
#'   delta-method standard errors (and, for \code{parsurv} and parametric \code{csc}, pointwise
#'   limits: complementary log-log scale for probabilities, log scale for positive quantities).}
#' \item{\code{level}}{confidence level of the limits.}
#' \item{\code{average}}{\code{parsurv} and parametric \code{csc}: the average over the rows of
#'   \code{newdata} (marginal curve), with its standard error.}
#' \item{\code{...}}{not used.}
#' }
#' @section Choosing what to predict:
#' \itemize{
#'   \item For one event type: survival probabilities at clinically meaningful times
#'     (\code{type = "survival"}), the median time (\code{parsurv}, \code{type = "quantile"}) or
#'     the restricted mean (\code{type = "rmst"}).
#'   \item With competing risks: the cumulative incidence of the cause of interest, from
#'     \code{fg}, \code{csc} or \code{directcif}; \code{csc} gives incidences of all causes that
#'     add up correctly with the event-free probability.
#'   \item \emph{Individual} predictions (one row of \code{newdata} per profile) describe typical
#'     patients; \emph{marginal} predictions (the average over the subjects of a population,
#'     \code{average = TRUE}, or \code{\link{cif_data}} with \code{type = "marginal"}) give
#'     standardised, population-level risks, for example to compare treatments as if everybody had
#'     received each.
#'   \item Cox and Fine-Gray baselines are step functions and are not defined beyond the last
#'     event time; parametric models give smooth curves and can extrapolate, with the caveat that the
#'     extrapolation rests on the distribution only.
#' }
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}): Cox predictions agree with \code{survival::survfit.coxph}, Fine-Gray
#' predictions with \code{cmprsk::predict.crr} (\eqn{2 \times 10^{-8}}), cause-specific predictions
#' with \code{riskRegression::CSC} (\eqn{4 \times 10^{-7}}), and parametric predictions with
#' \code{flexsurv} within the tolerance of its optimiser. In a Monte Carlo study the 95\% intervals of
#' parametric survival predictions covered 93.5 to 94.9\%.
#' @return
#' a matrix with one row per row of \code{newdata} and one column per time point
#'   (a vector for \code{type = "lp"}); with \code{se = TRUE} a list with \code{cif} and \code{se}
#'   (\code{directcif}) or \code{est}, \code{se}, \code{lower} and \code{upper} (\code{parsurv},
#'   parametric \code{csc}).
#' @references
#' Breslow NE (1974) Covariance analysis of censored survival data. Biometrics 30:89-99.
#'
#' Prentice RL, Kalbfleisch JD, Peterson AV, Flournoy N, Farewell VT, Breslow NE (1978) The analysis
#'   of failure times in the presence of competing risks. Biometrics 34:541-554.
#'
#' Benichou J, Gail MH (1990) Estimates of absolute cause-specific risk in cohort studies. Biometrics
#'   46:813-826.
#'
#' Fine JP, Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk.
#'   Journal of the American Statistical Association 94:496-509.
#'
#' Kalbfleisch JD, Prentice RL (2002) The Statistical Analysis of Failure Time Data, 2nd ed. Wiley.
#'
#' Hinchliffe SR, Lambert PC (2013) Flexible parametric modelling of cause-specific hazards to
#'   estimate cumulative incidence functions. BMC Medical Research Methodology 13:13.
#'
#' Casella G, Berger RL (2002) Statistical Inference, 2nd ed. Duxbury, Section 5.5 (delta method).
#' @seealso
#' \code{\link{cif_data}} and \code{\link{cif_plot}} (curves with bootstrap intervals),
#' \code{\link{hazard_plot}}, \code{\link{cox}}, \code{\link{fg}}, \code{\link{csc}},
#' \code{\link{parsurv}}, \code{\link{directcif}}, \code{\link{cure}}.
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(treat = rbinom(n, 1, 0.5), age = round(rnorm(n, 60, 10)))
#' t <- rweibull(n, 1.4, 5 * exp(0.4 * d$treat - 0.02 * (d$age - 60))); cens <- runif(n, 0, 10)
#' d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
#' nd <- data.frame(treat = 0:1, age = 60)
#'
#' ## Cox: survival at 2 and 5 years for two profiles
#' fc <- cox(Event(time, status) ~ treat + age, d)
#' predict(fc, nd, times = c(2, 5), type = "survival")
#'
#' ## Weibull: survival with limits, median time, restricted mean up to 5 years
#' fw <- parsurv(Event(time, status) ~ treat + age, d)
#' predict(fw, nd, times = c(2, 5), se = TRUE)
#' predict(fw, nd, type = "quantile", p = 0.5)
#' predict(fw, nd, times = 5, type = "rmst")
#'
#' ## marginal (standardised) survival: everybody treated versus nobody treated
#' d1 <- transform(d, treat = 1); d0 <- transform(d, treat = 0)
#' c(treated = predict(fw, d1, times = 5, average = TRUE),
#'   untreated = predict(fw, d0, times = 5, average = TRUE))
#' @name predict.SurvX
NULL

S7::method(predict, CoxFit) <- function(object, newdata = NULL, times = NULL, type = c("lp", "survival"), ...) {
  .no_tt(object, "predict()")
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

#' Hazard ratio at chosen times for time-varying coefficients
#'
#' @description
#' For a Cox model with \code{tt()} terms (see \code{\link{cox}}), the hazard ratio of each such
#' covariate at the requested times, \eqn{\exp\{\beta + \gamma\, g(t)\}}, with a Wald confidence
#' interval. It turns the two coefficients of a time-varying effect into the quantity to report:
#' "the hazard ratio of treatment was 0.55 at 6 months and 0.85 at 3 years".
#' @param fit a model fitted by \code{\link{cox}} with \code{tt()} terms.
#' @param times the times at which the effect is wanted; default the quartiles of the event times.
#' @param level confidence level.
#' @details
#' For a covariate \eqn{x} entered as \code{x + tt(x)} with \code{tt = g}, the log hazard ratio at
#' time \eqn{t} is \eqn{\beta + \gamma g(t)} with variance
#' \eqn{V_{\beta\beta} + g(t)^2 V_{\gamma\gamma} + 2 g(t) V_{\beta\gamma}}. When \code{x} itself is not
#' in the model, \eqn{\beta = 0}.
#' @return
#' A data frame with \code{term}, \code{time}, \code{g} (the value of the function of time),
#' \code{log_hr}, \code{se}, \code{hr}, \code{lower} and \code{upper}.
#' @references
#' Therneau TM, Grambsch PM (2000). Modeling Survival Data: Extending the Cox Model. Springer,
#'   Chapter 6.
#' @seealso
#' \code{\link{cox}}, \code{\link{ph_test}}.
#' @examples
#' set.seed(4)
#' n <- 300
#' x <- rbinom(n, 1, 0.5)
#' t <- rweibull(n, ifelse(x == 1, 0.7, 1.5), 3); cens <- runif(n, 0, 6)
#' d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x)
#' fit <- cox(Event(time, status) ~ x + tt(x), d, tt = function(t) log(t))
#' tt_effect(fit, times = c(0.5, 1, 2, 4))
#' @export
tt_effect <- function(fit, times = NULL, level = 0.95) {
  if (!S7::S7_inherits(fit, CoxFit) || is.null(fit@tt)) stop("fit must be a cox() model with tt() terms", call. = FALSE)
  if (is.null(times)) times <- unname(stats::quantile(fit@time[fit@status == 1L], c(0.25, 0.5, 0.75), type = 1))
  b <- fit@coefficients; V <- fit@vcov; nm <- names(b)
  funs <- if (is.function(fit@tt$fun)) rep(list(fit@tt$fun), length(fit@tt$cols)) else fit@tt$fun
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- lapply(seq_along(fit@tt$cols), function(j) {
    jt <- fit@tt$cols[j]
    base <- sub("^tt\\((.*)\\)$", "\\1", nm[jt])
    jb <- match(base, nm)
    g <- funs[[j]](times)
    est <- (if (is.na(jb)) 0 else b[jb]) + b[jt] * g
    v <- V[jt, jt] * g^2 + if (is.na(jb)) 0 else V[jb, jb] + 2 * g * V[jb, jt]
    se <- sqrt(pmax(v, 0))
    data.frame(term = base, time = times, g = g, log_hr = unname(est), se = unname(se), hr = exp(unname(est)),
               lower = exp(unname(est - z * se)), upper = exp(unname(est + z * se)), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
