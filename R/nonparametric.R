#' K-sample tests
#'
#' \code{logrank_test} is the (stratified) log-rank test and, with \code{rho > 0}, the
#' G-rho family of Harrington and Fleming (1982). \code{gray_test} is Gray's (1988) test for
#' equality of the cumulative incidence of one cause across groups.
#'
#' @param time,status follow-up time and status (0 = censored; for \code{gray_test} 1..K = cause).
#' @param group grouping variable.
#' @param strata optional stratification variable.
#' @param rho exponent of the weight function.
#' @param cause cause of interest.
#' @return a list with \code{statistic}, \code{df}, \code{p.value} and, for the log-rank test,
#'   observed and expected numbers of events per group.
#' @references Mantel N (1966) Cancer Chemother Rep 50:163-170.
#'   Harrington DP, Fleming TR (1982) Biometrika 69:553-566.
#'   Gray RJ (1988) Ann Stat 16:1141-1154.
#' @export
logrank_test <- function(time, status, group, strata = NULL, rho = 0) {
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
gray_test <- function(time, status, group, cause = 1L, strata = NULL, rho = 0) {
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

#' Kaplan-Meier estimate
#'
#' Product-limit estimate of the survival function with Greenwood standard errors, one curve
#' per group, and the log-rank test when there is more than one group.
#'
#' @param formula \code{Event(time, status) ~ 1} or \code{~ group}; several variables on the
#'   right-hand side are combined into one grouping factor.
#' @param data a data frame.
#' @param level confidence level of the intervals.
#' @return an S7 object of class \code{KMFit}; \code{x@curves} is a list of data frames with columns
#'   \code{time, n.risk, n.event, surv, se, lower, upper} (intervals on the log scale), \code{x@counts}
#'   has the median survival time of every group with its confidence limits (Brookmeyer and Crowley 1982).
#'   Methods: \code{print}, \code{summary}, \code{plot}, \code{as.data.frame}.
#' @references Kaplan EL, Meier P (1958) J Am Stat Assoc 53:457-481.
#'   Greenwood M (1926) Reports on Public Health and Medical Subjects 33:1-26.
#'   Brookmeyer R, Crowley J (1982) Biometrics 38:29-41.
#' @examples
#' set.seed(1)
#' d <- data.frame(time = rexp(100), status = rbinom(100, 1, 0.7), arm = gl(2, 50))
#' fit <- km(Event(time, status) ~ arm, d)
#' fit
#' @export
km <- function(formula, data, level = 0.95) {
  cl <- match.call()
  f <- .sx_groups(formula, data)
  if (any(f$status > 1L)) stop("km() needs a 0/1 status; use aj() for competing risks", call. = FALSE)
  zq <- stats::qnorm(1 - (1 - level) / 2)
  curves <- list(); counts <- NULL
  for (lev in levels(f$group)) {
    i <- f$group == lev; n <- sum(i)
    r <- .C("SX_km", as.integer(n), f$time[i], f$status[i], time = double(n), nrisk = double(n), nevent = double(n),
            surv = double(n), se = double(n), m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
    .sx_rc(r$rc, "km")
    k <- seq_len(r$m)
    S <- r$surv[k]; se <- r$se[k]
    w <- suppressWarnings(exp(zq * se / S))                 # interval on the log scale
    cv <- data.frame(time = r$time[k], n.risk = r$nrisk[k], n.event = r$nevent[k], surv = S, se = se,
                     lower = ifelse(S > 0, S / w, NA_real_), upper = ifelse(S > 0, pmin(1, S * w), NA_real_))
    curves[[lev]] <- cv
    first <- function(v) { j <- which(v <= 0.5); if (length(j)) cv$time[j[1L]] else NA_real_ }
    counts <- rbind(counts, data.frame(group = lev, n = n, events = sum(f$status[i]),
                                       median = first(cv$surv), median.lower = first(cv$lower), median.upper = first(cv$upper),
                                       stringsAsFactors = FALSE))
  }
  test <- if (nlevels(f$group) > 1L) logrank_test(f$time, f$status, f$group) else NULL
  KMFit(call = cl, curves = curves, counts = counts, test = test, level = level, time = f$time, status = f$status, group = f$group)
}

#' Aalen-Johansen estimate of the cumulative incidence functions
#'
#' Cumulative incidence of every cause with delta-method standard errors, one set of curves per
#' group, and Gray's test for every cause when there is more than one group.
#'
#' @param formula \code{Event(time, status) ~ 1} or \code{~ group}; status 0 = censored, 1..K = cause.
#' @param data a data frame.
#' @param level confidence level of the intervals.
#' @return an S7 object of class \code{AJFit}; \code{x@curves} is a list of data frames with columns
#'   \code{time, surv, cif.k, se.k, lower.k, upper.k} for every cause k (intervals on the
#'   complementary log-log scale). Methods: \code{print}, \code{summary}, \code{plot}, \code{as.data.frame}.
#' @references Aalen OO, Johansen S (1978) Scand J Stat 5:141-150.
#'   Gray RJ (1988) Ann Stat 16:1141-1154.
#' @export
aj <- function(formula, data, level = 0.95) {
  cl <- match.call()
  f <- .sx_groups(formula, data)
  K <- max(f$status)
  if (K < 1L) stop("no events", call. = FALSE)
  zq <- stats::qnorm(1 - (1 - level) / 2)
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
  AJFit(call = cl, curves = curves, K = as.integer(K), counts = counts, tests = tests, level = level,
        time = f$time, status = f$status, group = f$group)
}

#' Restricted mean survival time
#'
#' Area under the Kaplan-Meier curve on [0, tau] with its standard error.
#' @param time,status follow-up time and 0/1 status.
#' @param tau upper limit of integration.
#' @return a named vector with \code{rmst} and \code{se}.
#' @references Royston P, Parmar MKB (2013) BMC Med Res Methodol 13:152.
#' @export
rmst <- function(time, status, tau) {
  n <- length(time)
  r <- .C("SX_rmst", as.integer(n), as.double(time), as.integer(status > 0), as.double(tau), mean = double(1), se = double(1),
          rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "rmst")
  c(rmst = r$mean, se = r$se)
}

