#' Comparing groups: log-rank test and Gray's test
#'
#' @description
#' Tests whether two or more groups differ, without any model for the covariates.
#' \itemize{
#'   \item \code{logrank_test()}: equality of the \emph{survival curves} (one event type), the
#'     log-rank test (Mantel 1966; Peto and Peto 1972) and the \eqn{G^\rho} family of Harrington and
#'     Fleming (1982), optionally stratified.
#'   \item \code{gray_test()}: equality of the \emph{cumulative incidence} of one cause when there are
#'     competing events (Gray 1988).
#' }
#' \code{\link{km}} and \code{\link{aj}} run these tests automatically when the formula has groups;
#' call these functions directly for the \eqn{G^\rho} weights, for a stratified test, or when the
#' data are in vectors.
#'
#' \strong{Use something else when} the effect has to be adjusted for other covariates
#' (\code{\link{cox}}, \code{\link{csc}}, \code{\link{fg}}), or when you need an effect size
#' (a hazard ratio from \code{cox}, or a difference in restricted mean survival time from
#' \code{\link{rmst}}).
#' @param time follow-up time (numeric vector).
#' @param status \code{logrank_test}: 0 = censored, any positive value = event.
#'   \code{gray_test}: 0 = censored, 1, 2, ... = cause.
#' @param group the groups to compare (any vector; it is converted to a factor). At least two
#'   groups.
#' @param strata optional stratification variable: the comparison is made within each stratum and
#'   the results are combined, so that a difference between strata (centres, for example) is not
#'   mistaken for a difference between groups.
#' @param rho \code{logrank_test}: the \eqn{G^\rho} weight \eqn{\hat S(t-)^\rho}. \code{0} (default)
#'   is the log-rank test, which has the most power when the hazards are proportional; \code{1} is
#'   the Peto-Peto version of the Wilcoxon test, which gives more weight to early differences.
#'   \code{gray_test}: the exponent of Gray's weight function (0 by default, as in
#'   \code{cmprsk::cuminc}).
#' @param cause \code{gray_test}: the cause whose cumulative incidence is compared.
#' @param timefix treat times that differ only by rounding as tied (see \code{\link{Event}}).
#' @details
#' \subsection{Log-rank and \eqn{G^\rho} tests}{
#' At every event time \eqn{t}, the observed number of events in each group is compared with the
#' number expected if all groups had the same hazard, given the numbers at risk. With weights
#' \eqn{w(t) = \hat S(t-)^\rho} (\eqn{\hat S} the pooled Kaplan-Meier estimate), the statistic is
#' \eqn{U^\top V^{-} U}, where \eqn{U} sums \eqn{w(t)(O - E)} over the event times and \eqn{V} is
#' its hypergeometric variance; under the null hypothesis it has a chi-square distribution with
#' (number of groups - 1) degrees of freedom. Stratified: \eqn{U} and \eqn{V} are summed over strata.
#' }
#' \subsection{Gray's test}{
#' Compares weighted averages of the subdistribution hazards of the cause across groups (Gray 1988,
#' with the variance of that paper), as in \code{cmprsk::cuminc}. It tests the cumulative
#' incidence, not the cause-specific hazard: a treatment that reduces only competing deaths can
#' increase the cumulative incidence of the cause, and Gray's test then detects that increase.
#' }
#' @section Interpreting the result:
#' A small p-value means that the curves (survival, or cumulative incidence of the cause) are not
#' all equal. The test says nothing about the size of the difference: report the curves
#' (\code{plot(km(...))}, \code{plot(aj(...))}), the median survival per group, or an effect
#' estimate. \code{observed} and \code{expected} of the log-rank test show which group has more
#' events than expected. The log-rank test can miss curves that cross (the early and late
#' differences cancel); then look at the curves and consider \code{rho = 1} or \code{\link{rmst}}.
#'
#' A sentence for a report: \emph{"Overall survival differed between the arms (log-rank test,
#' chi-square 7.9 on 1 df, p = 0.005)."}
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}): the log-rank, Peto-Peto and stratified statistics agree with
#' \code{survival::survdiff} to a relative \eqn{3 \times 10^{-6}}, and Gray's statistics with
#' \code{cmprsk::cuminc} to \eqn{10^{-11}}, on 81 data sets. In a Monte Carlo study (1000
#' replications) the size of the log-rank test was 4.1 to 5.7\% and that of Gray's test 3.6 to 4.9\%
#' at the 5\% level.
#' @return
#' A list with \code{method}, \code{statistic} (chi-square), \code{df}, \code{p.value}; for
#' \code{logrank_test} also \code{observed} and \code{expected} (number of events per group), for
#' \code{gray_test} also \code{cause}.
#' @references
#' \emph{Log-rank test:}
#' Mantel N (1966). Evaluation of survival data and two new rank order statistics arising in its
#'   consideration. Cancer Chemotherapy Reports 50:163-170.
#' Peto R, Peto J (1972). Asymptotically efficient rank invariant test procedures. Journal of the
#'   Royal Statistical Society, Series A 135:185-207.
#'
#' \emph{\eqn{G^\rho} family:}
#' Harrington DP, Fleming TR (1982). A class of rank test procedures for censored survival data.
#'   Biometrika 69:553-566.
#'
#' \emph{Gray's test:}
#' Gray RJ (1988). A class of K-sample tests for comparing the cumulative incidence of a competing
#'   risk. Annals of Statistics 16:1141-1154.
#'
#' All references: \code{sx_references("logrank_test")}, \code{sx_references("gray_test")}.
#' @seealso
#' \code{\link{km}} and \code{\link{aj}} (curves with the tests), \code{\link{rmst}},
#' \code{\link{cox}}, \code{\link{fg}}.
#' @examples
#' set.seed(1)
#' n <- 200
#' arm <- rep(c("A", "B"), each = n / 2)
#' centre <- sample(c("c1", "c2"), n, TRUE)
#' t <- rexp(n, ifelse(arm == "B", 0.7, 1)); cens <- runif(n, 0, 3)
#' time <- pmin(t, cens); status <- as.integer(t <= cens)
#'
#' logrank_test(time, status, arm)                    # log-rank
#' logrank_test(time, status, arm, rho = 1)           # Peto-Peto (early differences)
#' logrank_test(time, status, arm, strata = centre)   # stratified by centre
#'
#' ## competing risks: does the cumulative incidence of cause 1 differ?
#' cause <- ifelse(status == 1, sample(1:2, n, TRUE, prob = c(0.6, 0.4)), 0)
#' gray_test(time, cause, arm, cause = 1)
#' @export
logrank_test <- function(time, status, group, strata = NULL, rho = 0, timefix = TRUE) {
  if (isTRUE(timefix)) time <- .timefix_vec(time)
  g <- droplevels(as.factor(group)); K <- nlevels(g); n <- length(time)
  if (K < 2L) stop("at least two groups are needed", call. = FALSE)
  st <- if (is.null(strata)) integer(n) else as.integer(as.factor(strata)) - 1L
  r <- .C("SX_logrank", as.integer(n), as.double(time), as.integer(status > 0), as.integer(g) - 1L, as.integer(K),
          as.integer(!is.null(strata)), as.integer(st), as.double(rho),
          stat = double(1), df = integer(1), p = double(1), obs = double(K), exp = double(K), var = double(K * K),
          rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "logrank_test")
  list(method = if (rho == 0) "log-rank test" else paste0("G-rho test, rho = ", rho),
       statistic = r$stat, df = r$df, p.value = r$p,
       observed = stats::setNames(r$obs, levels(g)), expected = stats::setNames(r$exp, levels(g)))
}

#' @rdname logrank_test
#' @export
gray_test <- function(time, status, group, cause = 1L, strata = NULL, rho = 0, timefix = TRUE) {
  if (isTRUE(timefix)) time <- .timefix_vec(time)
  g <- droplevels(as.factor(group)); K <- nlevels(g); n <- length(time)
  if (K < 2L) stop("at least two groups are needed", call. = FALSE)
  st <- if (is.null(strata)) integer(n) else as.integer(as.factor(strata)) - 1L
  r <- .C("SX_gray", as.integer(n), as.double(time), as.integer(status), as.integer(g) - 1L, as.integer(K),
          as.integer(!is.null(strata)), as.integer(st), as.integer(cause), as.double(rho),
          stat = double(1), df = integer(1), p = double(1), score = double(K), var = double(K * K),
          rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "gray_test")
  list(method = "Gray's test", cause = as.integer(cause), statistic = r$stat, df = r$df, p.value = r$p)
}

#' Kaplan-Meier estimate of the survival function
#'
#' @description
#' Estimates the probability of being free of the event over time, \eqn{S(t) = P(T > t)}, without
#' any model: one curve per group, with confidence intervals, the median survival time of every
#' group, and the log-rank test when there are groups.
#'
#' \strong{Use it when} there is one type of event and you want to describe survival or compare
#' groups. It is usually the first step of an analysis.
#' \strong{Use something else when}
#' \itemize{
#'   \item there are competing events: \code{\link{aj}}. One minus the Kaplan-Meier estimate that
#'     treats competing events as censored is \emph{not} the probability of the event and
#'     overestimates it;
#'   \item there is delayed entry (left truncation): \code{km()} does not handle it;
#'     \code{\link{hazard_np}} estimates the cumulative hazard with delayed entry;
#'   \item covariates must be adjusted for: \code{\link{cox}} or \code{\link{parsurv}}.
#' }
#' @param formula \code{Event(time, status) ~ 1} for one curve, or \code{Event(time, status) ~ group}
#'   for one curve per group; several variables on the right (\code{~ sex + stage}) give one curve
#'   per combination. \code{status} must be 0/1.
#' @param data a data frame. Rows with a missing time, status or group are removed.
#' @param level confidence level of the intervals (default 0.95).
#' @details
#' \subsection{Estimator}{
#' \deqn{\hat S(t) = \prod_{t_j \le t} \left(1 - \frac{d_j}{Y_j}\right),}
#' with \eqn{d_j} the number of events and \eqn{Y_j} the number at risk at event time \eqn{t_j}
#' (Kaplan and Meier 1958). Its variance is Greenwood's (1926) formula; the confidence interval is
#' computed on the log scale, \eqn{\hat S \exp(\pm z\,se\{\log \hat S\})}, truncated at 1, the
#' default of \code{survival::survfit}.
#' }
#' \subsection{Median survival}{
#' The median is the first time at which \eqn{\hat S(t) \le 0.5}; its confidence limits are the first
#' times at which the lower and upper limits of the curve reach 0.5 (Brookmeyer and Crowley 1982).
#' When the curve does not fall to 0.5 during follow-up, the median is \code{NA}: report the
#' survival at fixed times instead.
#' }
#' @section Interpreting the result:
#' \itemize{
#'   \item \eqn{\hat S(t)} is the estimated proportion still free of the event at time \eqn{t}; it is
#'     valid when censoring is unrelated to the risk of the event.
#'   \item The curve is reliable only where enough subjects are at risk: the intervals widen at the
#'     end of follow-up. The risk table under \code{plot(fit)} shows the numbers at risk.
#'   \item The log-rank test in the print compares the whole curves; see \code{\link{logrank_test}}
#'     for when it can miss differences.
#' }
#' A sentence for a report: \emph{"Five-year survival was 62\% (95\% CI 55 to 70\%) in arm A and 48\%
#' (40 to 57\%) in arm B; median survival 6.8 and 4.9 years (log-rank p = 0.01)."}
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}), on 81 data sets: survival at the event times and Greenwood standard errors
#' agree with \code{survival::survfit} to \eqn{4 \times 10^{-16}}, medians and their confidence limits
#' exactly. In a Monte Carlo study (1000 replications, \eqn{n} = 200 to 1000) the 95\% intervals
#' covered 93.7 to 96.1\%.
#' @return
#' An S7 object of class \code{KMFit}. Slots: \code{curves} (a list with one data frame per group,
#' columns \code{time}, \code{n.risk}, \code{n.event}, \code{surv}, \code{se} (of \eqn{\log S}),
#' \code{lower}, \code{upper}), \code{counts} (per group: \code{n}, \code{events}, \code{median} with
#' \code{median.lower} and \code{median.upper}), \code{test} (the log-rank test) and \code{level}.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(fit)} \tab numbers, events, medians with intervals, log-rank test \cr
#'   \code{summary(fit, times = c(1, 3, 5))} \tab survival with intervals at chosen times \cr
#'   \code{plot(fit, ...)} \tab survival curves with a risk table; arguments \code{conf.int},
#'     \code{risk.table}, \code{median}, \code{mark.times}, \code{percent}, \code{labels}, \code{col},
#'     \code{xlim} and the usual graphical ones (see \code{\link{SurvX-methods}}) \cr
#'   \code{as.data.frame(fit)} \tab all curves in one data frame with a \code{group} column \cr
#' }
#' @references
#' \emph{Estimator and variance:}
#' Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete observations. Journal of the
#'   American Statistical Association 53:457-481.
#' Greenwood M (1926). The natural duration of cancer. Reports on Public Health and Medical Subjects
#'   33:1-26.
#'
#' \emph{Median and its confidence interval:}
#' Brookmeyer R, Crowley J (1982). A confidence interval for the median survival time. Biometrics
#'   38:29-41.
#'
#' \emph{Log-rank test in the print:}
#' Mantel N (1966). Evaluation of survival data and two new rank order statistics arising in its
#'   consideration. Cancer Chemotherapy Reports 50:163-170.
#'
#' All references: \code{sx_references("km")}.
#' @seealso
#' \code{\link{aj}} (competing risks), \code{\link{logrank_test}}, \code{\link{rmst}},
#' \code{\link{hazard_np}}, \code{\link{cox}}, \code{\link{Event}}.
#' @examples
#' set.seed(1)
#' n <- 200
#' d <- data.frame(arm = rep(c("A", "B"), each = n / 2))
#' t <- rexp(n, ifelse(d$arm == "B", 0.15, 0.1)); cens <- runif(n, 2, 12)
#' d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
#'
#' fit <- km(Event(time, status) ~ arm, d)
#' fit                                   # medians and the log-rank test
#' summary(fit, times = c(1, 3, 5))      # survival at 1, 3 and 5 years
#' plot(fit, median = TRUE)              # curves with a risk table
#' head(as.data.frame(fit))
#' @export
km <- function(formula, data, level = 0.95) {
  cl <- match.call()
  f <- .sx_groups(formula, data, what = "this estimator")
  if (any(f$status > 1L)) stop("km() needs a 0/1 status; use aj() for competing risks", call. = FALSE)
  zq <- .zq(level)
  curves <- list(); counts <- NULL
  for (lev in levels(f$group)) {
    i <- f$group == lev; n <- sum(i)
    r <- .C("SX_km", as.integer(n), f$time[i], f$status[i], time = double(n), nrisk = double(n), nevent = double(n),
            surv = double(n), se = double(n), m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "km")
    k <- seq_len(r$m)
    S <- r$surv[k]; se <- r$se[k]
    lim <- .C("SX_km_limits", length(S), as.double(S), as.double(se), as.double(zq), lower = double(length(S)),
              upper = double(length(S)), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")   # interval on the log scale
    .sx_rc(lim$rc, "km")
    cv <- data.frame(time = r$time[k], n.risk = r$nrisk[k], n.event = r$nevent[k], surv = S, se = se,
                     lower = ifelse(is.nan(lim$lower), NA_real_, lim$lower), upper = ifelse(is.nan(lim$upper), NA_real_, lim$upper))
    curves[[lev]] <- cv
    first <- function(v) {
      q <- .C("SX_first_below", length(v), as.double(cv$time), as.double(v), 0.5, out = double(1), rc = integer(1),
              NAOK = TRUE, PACKAGE = "SurvX")$out
      if (is.nan(q)) NA_real_ else q
    }
    counts <- rbind(counts, data.frame(group = lev, n = n, events = sum(f$status[i]),
                                       median = first(cv$surv), median.lower = first(cv$lower), median.upper = first(cv$upper),
                                       stringsAsFactors = FALSE))
  }
  test <- if (nlevels(f$group) > 1L) logrank_test(f$time, f$status, f$group) else NULL
  KMFit(model_call = cl, curves = curves, counts = counts, test = test, level = level, time = f$time, status = f$status, group = f$group)
}

#' Aalen-Johansen estimate of the cumulative incidence with competing risks
#'
#' @description
#' Estimates, without any model, the \emph{cumulative incidence} of every cause,
#' \eqn{F_k(t) = P(T \le t, \mathrm{cause} = k)}: the probability of having had event \eqn{k} by time
#' \eqn{t} when other events compete with it. One set of curves per group, with confidence intervals
#' and Gray's test for every cause when there are groups.
#'
#' \strong{Use it when} there are competing events (status 0 = censored, 1, 2, ... = causes) and you
#' want to describe the risk of each cause or compare groups. It replaces "1 minus Kaplan-Meier",
#' which treats competing events as censoring and overestimates the risk.
#' \strong{Use something else when} covariates must be adjusted for: \code{\link{fg}} or
#' \code{\link{directcif}} (cumulative incidence) and \code{\link{csc}} (cause-specific hazards).
#' @param formula \code{Event(time, status) ~ 1} or \code{Event(time, status) ~ group}; several
#'   variables on the right give one set of curves per combination. \code{status}: 0 = censored,
#'   1, 2, ..., K = cause. Delayed entry is not supported.
#' @param data a data frame. Rows with a missing time, status or group are removed.
#' @param level confidence level of the intervals (default 0.95).
#' @details
#' \subsection{Estimator}{
#' \deqn{\hat F_k(t) = \sum_{t_j \le t} \hat S(t_{j-1}) \frac{d_{kj}}{Y_j},}
#' with \eqn{\hat S} the Kaplan-Meier estimate of being free of all events, \eqn{d_{kj}} the number of
#' events of cause \eqn{k} and \eqn{Y_j} the number at risk at time \eqn{t_j} (Aalen and Johansen
#' 1978; Kalbfleisch and Prentice 2002, Section 8.2). At every time the cumulative incidences of all
#' causes and \eqn{\hat S} add up to one.
#' }
#' \subsection{Variance and intervals}{
#' Delta-method variance (Marubini and Valsecchi 1995, Section 10.3). The confidence interval is
#' computed on the complementary log-log scale, \eqn{\hat F^{\exp(\pm z\,se / (\hat F \log \hat F))}},
#' which keeps it inside (0, 1) and behaves well for small incidences.
#' }
#' @section Interpreting the result:
#' \itemize{
#'   \item \eqn{\hat F_k(t)} is the estimated proportion of subjects who have had cause \eqn{k} by time
#'     \eqn{t}, in a population where the other causes also occur. It depends on the competing
#'     causes: if fewer patients die of other causes, more can have cause \eqn{k}.
#'   \item Gray's test in the print compares the cumulative incidence of each cause across groups
#'     (see \code{\link{gray_test}}).
#' }
#' A sentence for a report: \emph{"The 5-year cumulative incidence of relapse was 24\% (95\% CI 19 to
#' 30\%), and of death without relapse 11\% (8 to 15\%) (Aalen-Johansen estimator)."}
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}), on 82 data sets: cumulative incidences agree with \code{cmprsk::cuminc} and
#' with the multi-state \code{survival::survfit} to \eqn{10^{-14}}, and Gray's statistics with
#' \code{cuminc} to \eqn{10^{-11}}. The standard errors differ slightly from those of \code{survfit},
#' which uses the infinitesimal jackknife. In a Monte Carlo study (1000 replications, \eqn{n} = 200
#' to 1000) the 95\% intervals covered 93.6 to 95.6\%.
#' @return
#' An S7 object of class \code{AJFit}. Slots: \code{curves} (one data frame per group with columns
#' \code{time}, \code{surv}, and for every cause \eqn{k}: \code{cif.k}, \code{se.k}, \code{lower.k},
#' \code{upper.k}), \code{counts} (numbers censored and per cause), \code{tests} (Gray's test per
#' cause), \code{K} (number of causes) and \code{level}.
#'
#' Methods:
#' \tabular{ll}{
#'   \code{print(fit)} \tab counts and Gray's tests \cr
#'   \code{summary(fit, times = c(1, 3, 5))} \tab cumulative incidences with intervals at chosen
#'     times \cr
#'   \code{plot(fit, ...)} \tab cumulative incidence curves (\code{\link{cif_plot}}): panels per cause,
#'     overlaid or stacked, with risk tables \cr
#'   \code{as.data.frame(fit)} \tab all curves in one data frame \cr
#' }
#' @references
#' \emph{Estimator:}
#' Aalen OO, Johansen S (1978). An empirical transition matrix for non-homogeneous Markov chains
#'   based on censored observations. Scandinavian Journal of Statistics 5:141-150.
#' Kalbfleisch JD, Prentice RL (2002). The Statistical Analysis of Failure Time Data, 2nd ed.
#'   Wiley, Section 8.2.
#'
#' \emph{Variance:}
#' Marubini E, Valsecchi MG (1995). Analysing Survival Data from Clinical Trials and Observational
#'   Studies. Wiley, Section 10.3.
#' Andersen PK, Borgan O, Gill RD, Keiding N (1993). Statistical Models Based on Counting Processes.
#'   Springer, Section IV.4.
#'
#' \emph{Gray's test in the print:}
#' Gray RJ (1988). A class of K-sample tests for comparing the cumulative incidence of a competing
#'   risk. Annals of Statistics 16:1141-1154.
#'
#' \emph{Why not 1 minus Kaplan-Meier:}
#' Putter H, Fiocco M, Geskus RB (2007). Tutorial in biostatistics: competing risks and multi-state
#'   models. Statistics in Medicine 26:2389-2430.
#'
#' All references: \code{sx_references("aj")}.
#' @seealso
#' \code{\link{cif_plot}}, \code{\link{gray_test}}, \code{\link{fg}}, \code{\link{csc}},
#' \code{\link{km}}, \code{\link{Event}}.
#' @examples
#' set.seed(1)
#' n <- 300
#' d <- data.frame(arm = rep(c("A", "B"), each = n / 2))
#' t1 <- rexp(n, ifelse(d$arm == "B", 0.06, 0.1))   # relapse
#' t2 <- rexp(n, 0.04)                               # death without relapse
#' cens <- runif(n, 2, 10)
#' d$time <- pmin(t1, t2, cens)
#' d$status <- ifelse(d$time == cens, 0, ifelse(t1 < t2, 1, 2))
#'
#' fit <- aj(Event(time, status) ~ arm, d)
#' fit                                  # Gray's test for each cause
#' summary(fit, times = c(1, 3, 5))
#' plot(fit)
#' @export
aj <- function(formula, data, level = 0.95) {
  cl <- match.call()
  f <- .sx_groups(formula, data, what = "this estimator")
  K <- max(f$status)
  if (K < 1L) stop("no events", call. = FALSE)
  zq <- .zq(level)
  curves <- list(); counts <- NULL
  for (lev in levels(f$group)) {
    i <- f$group == lev; n <- sum(i)
    r <- .C("SX_aj", as.integer(n), f$time[i], f$status[i], as.integer(K), time = double(n), surv = double(n),
            cif = double(n * K), se = double(n * K), m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "aj")
    m <- r$m; k <- seq_len(m)
    cif <- .matrix_rm(r$cif[seq_len(m * K)], m, K, list(NULL, paste0("cif.", seq_len(K))))
    se <- .matrix_rm(r$se[seq_len(m * K)], m, K, list(NULL, paste0("se.", seq_len(K))))
    ci <- .cif_interval(cif, se, zq)
    colnames(ci$lower) <- paste0("lower.", seq_len(K)); colnames(ci$upper) <- paste0("upper.", seq_len(K))
    curves[[lev]] <- data.frame(time = r$time[k], surv = r$surv[k], cif, se, ci$lower, ci$upper)
    cnt <- data.frame(group = lev, n = n, censored = sum(f$status[i] == 0L), stringsAsFactors = FALSE)
    for (j in seq_len(K)) cnt[[paste0("cause.", j)]] <- sum(f$status[i] == j)
    counts <- rbind(counts, cnt)
  }
  tests <- NULL
  if (nlevels(f$group) > 1L) {
    tests <- do.call(rbind, lapply(seq_len(K), function(j) {
      g <- tryCatch(gray_test(f$time, f$status, f$group, cause = j), error = function(e) NULL)
      if (is.null(g)) data.frame(cause = j, statistic = NA_real_, df = NA_integer_, p.value = NA_real_)
      else data.frame(cause = j, statistic = g$statistic, df = g$df, p.value = g$p.value)
    }))
  }
  AJFit(model_call = cl, curves = curves, K = as.integer(K), counts = counts, tests = tests, level = level,
        time = f$time, status = f$status, group = f$group)
}

#' Restricted mean survival time
#'
#' @description
#' The restricted mean survival time (RMST) up to a horizon \eqn{\tau}: the expected time lived
#' without the event during the first \eqn{\tau} units of follow-up, estimated by the area under the
#' Kaplan-Meier curve on \eqn{[0, \tau]}, with its standard error.
#'
#' \strong{Use it when} an effect measure is wanted that does not need proportional hazards and is
#' easy to communicate ("on average, patients on the new treatment lived 2.1 months longer during
#' the first 5 years"), in particular when the survival curves cross.
#' \strong{Use something else when} there are competing events (the restricted mean time lost to a
#' cause is the area under its cumulative incidence, from \code{\link{aj}}), or covariates must be
#' adjusted for (fit \code{\link{cox}} or \code{\link{parsurv}} and average the predicted curves).
#' @param time follow-up time.
#' @param status 0 = censored, 1 = event.
#' @param tau the horizon, in the unit of \code{time}. Choose it before looking at the results (a
#'   clinically meaningful time such as 5 years), and not beyond the follow-up of the data: beyond
#'   the last observed time the Kaplan-Meier curve is not estimated, and when the last time is
#'   censored the area is extended with the last value of the curve.
#' @details
#' \deqn{\widehat{RMST}(\tau) = \int_0^\tau \hat S(t)\,dt = \sum_j \hat S(t_{j-1})(t_j - t_{j-1}),}
#' summed over the event times up to \eqn{\tau} (Irwin 1949). The variance is the Greenwood-type
#' estimator \eqn{\sum_j A_j^2 d_j / \{Y_j (Y_j - d_j)\}} with \eqn{A_j = \int_{t_j}^\tau \hat S(t)dt}
#' (Klein and Moeschberger 2003, Section 4.5), the one of \code{survival::survfit(..., rmean = tau)}.
#' For two groups, the difference of the RMSTs has standard error
#' \eqn{\sqrt{se_1^2 + se_2^2}} (independent groups); see the example.
#' @section Interpreting the result:
#' An RMST of 3.8 years with \eqn{\tau} = 5 years means that, on average, a patient lived 3.8 of the
#' first 5 years without the event (and lost 1.2 years to it). The difference between two groups is
#' in time units, the same at every point of the curve, and needs no proportional hazards (Royston
#' and Parmar 2013). Always state \eqn{\tau} together with the result.
#'
#' A sentence for a report: \emph{"Over 5 years, the restricted mean survival time was 3.8 years
#' (95\% CI 3.5 to 4.1) with the new treatment and 3.2 years (2.9 to 3.5) with the standard; difference
#' 0.6 years (0.2 to 1.0)."}
#' @section Validation:
#' In the validation suite of the package (\url{https://github.com/stats9/SurvX}, directory
#' \code{validation}): the RMST agrees with \code{survival::survfit(..., rmean = tau)} to
#' \eqn{5 \times 10^{-12}} and with \code{survRM2::rmst2} to \eqn{2 \times 10^{-15}}, its standard error
#' with \code{survfit} to \eqn{2 \times 10^{-14}}, on 81 data sets. In a Monte Carlo study (1000
#' replications, \eqn{n} = 200 to 1000) the standard error matched the empirical standard deviation
#' (ratio 0.99 to 1.03) and the 95\% intervals covered 94.9 to 95.3\%.
#' @return
#' A named numeric vector: \code{rmst} (the estimate) and \code{se} (its standard error).
#' @references
#' \emph{Restricted mean:}
#' Irwin JO (1949). The standard error of an estimate of expectation of life, with special reference
#'   to expectation of tumourless life in experiments with mice. Journal of Hygiene 47:188-189.
#' Royston P, Parmar MKB (2013). Restricted mean survival time: an alternative to the hazard ratio
#'   for the design and analysis of randomized trials with a time-to-event outcome. BMC Medical
#'   Research Methodology 13:152.
#'
#' \emph{Variance:}
#' Klein JP, Moeschberger ML (2003). Survival Analysis: Techniques for Censored and Truncated Data,
#'   2nd ed. Springer, Section 4.5.
#'
#' All references: \code{sx_references("rmst")}.
#' @seealso
#' \code{\link{km}} (the curve whose area this is), \code{\link{logrank_test}}.
#' @examples
#' set.seed(1)
#' n <- 300
#' arm <- rep(c("standard", "new"), each = n / 2)
#' t <- rexp(n, ifelse(arm == "new", 0.15, 0.22)); cens <- runif(n, 2, 8)
#' time <- pmin(t, cens); status <- as.integer(t <= cens)
#'
#' ## RMST over the first 5 years in each arm
#' a <- rmst(time[arm == "new"], status[arm == "new"], tau = 5)
#' b <- rmst(time[arm == "standard"], status[arm == "standard"], tau = 5)
#' rbind(new = a, standard = b)
#'
#' ## difference with a 95\% confidence interval
#' diff <- a[["rmst"]] - b[["rmst"]]
#' se <- sqrt(a[["se"]]^2 + b[["se"]]^2)
#' c(difference = diff, lower = diff - 1.96 * se, upper = diff + 1.96 * se)
#' @export
rmst <- function(time, status, tau) {
  n <- length(time)
  r <- .C("SX_rmst", as.integer(n), as.double(time), as.integer(status > 0), as.double(tau), mean = double(1), se = double(1),
          rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "rmst")
  c(rmst = r$mean, se = r$se)
}

