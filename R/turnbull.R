#' Non-parametric survival curves for interval-censored data (Turnbull NPMLE)
#'
#' @description
#' The non-parametric maximum likelihood estimate (NPMLE) of the survival function when event times
#' are only known to lie in intervals (\code{\link{Interval}}): the interval-censored counterpart of
#' the Kaplan-Meier estimate, to which it reduces when all times are exact or right-censored. Curves
#' can be computed by group, with bootstrap confidence bands, medians, and a logrank-type test of equal
#' survival for interval-censored data.
#'
#' @param formula \code{Interval(left, right) ~ 1} for one curve, or \code{Interval(left, right) ~ g}
#'   for one curve per level of \code{g} (several variables: one curve per combination).
#' @param data a data frame.
#' @param times times at which the survival table is printed (default: about five round times over
#'   the follow-up).
#' @param B number of bootstrap samples for the confidence bands and the confidence interval of the
#'   median (resampling subjects within groups); \code{0} for none.
#' @param level confidence level.
#' @param test with groups, compute the test of equal survival.
#' @param seed seed of the bootstrap (group k uses seed + k - 1): a whole number, or \code{NULL} to
#'   take the seed from the random numbers of the session (so that \code{set.seed()} before the call
#'   governs the result). All random draws come from R's generator (see \code{RNGkind()}); with a
#'   number the state of the session's generator is left as it was.
#' @param max_iter,tol maximum number of iterations and convergence tolerance of the algorithm.
#'
#' @details
#' \strong{Support.} The NPMLE puts all its probability on the innermost (Turnbull) intervals:
#' the intervals between a left end point and the next right end point of the data with no other end
#' point in between (Peto 1973; Turnbull 1976). How the probability is spread within an innermost
#' interval is not identified by the data; the survival curve is therefore only determined outside
#' these intervals. \code{predict()} gives the lower and upper limit in addition to the estimate, and
#' the curve is drawn linearly within the innermost intervals (the usual convention); exact event
#' times are point masses, as in the Kaplan-Meier estimate. After the last left end point of the
#' right-censored subjects the curve is undefined (as the Kaplan-Meier curve after the last time).
#'
#' \strong{Computation.} In compiled code, steps of the self-consistency (EM) algorithm of Turnbull
#' (1976) alternate with steps of the iterative convex minorant algorithm (Groeneboom and Wellner
#' 1992) with a line search (Jongbloed 1998), the EM-ICM hybrid of Wellner and Zhan (1997), which
#' converges in a fraction of a second also with thousands of subjects. Convergence is checked on the Kuhn-Tucker
#' conditions, which are necessary and sufficient for the maximum of this concave problem (Gentleman
#' and Geyer 1994): the printed value is zero at the NPMLE.
#'
#' \strong{Uncertainty.} The NPMLE converges at a rate slower than \eqn{\sqrt{n}} and has no simple
#' standard error; the confidence bands are pointwise bootstrap percentile intervals.
#'
#' \strong{Test.} With groups, the weighted logrank test for interval-censored data with the scores
#' of Finkelstein (1986), \eqn{c_i = \{S(L_i)\log S(L_i) - S(R_i)\log S(R_i)\}/\{S(L_i) - S(R_i)\}}
#' from the NPMLE of all subjects together, and the variance from the permutation distribution (the
#' permutational central limit theorem, \code{method = "pclt"} in Fay and Shaw 2010). It reduces to
#' the ordinary logrank test (permutation form) for right-censored data. A positive group score means
#' that the events happen earlier in that group than in the groups together.
#'
#' @section Messages:
#' \describe{
#'   \item{\code{did not converge}}{the iteration limit was reached before the Kuhn-Tucker condition
#'     was met; increase \code{max_iter}.}
#' }
#' @return
#' An S7 object of class \code{TurnbullFit}. Use \code{print()}, \code{summary()} (adds the innermost
#' intervals and their probabilities), \code{plot()}, \code{predict(fit, times)} and
#' \code{as.data.frame(fit, what = c("survival", "intervals", "test"))}.
#' @references
#' Turnbull BW (1976). The empirical distribution function with arbitrarily grouped, censored and
#'   truncated data. Journal of the Royal Statistical Society B 38:290-295.
#'
#' Peto R (1973). Experimental survival curves for interval-censored data. Applied Statistics
#'   22:86-91.
#'
#' Gentleman R, Geyer CJ (1994). Maximum likelihood for interval censored data: consistency and
#'   computation. Biometrika 81:618-623.
#'
#' Groeneboom P, Wellner JA (1992). Information Bounds and Nonparametric Maximum Likelihood
#'   Estimation. Birkhauser, Basel.
#' 
#' Jongbloed G (1998). The iterative convex minorant algorithm for nonparametric estimation.
#'   Journal of Computational and Graphical Statistics 7:310-321.
#' 
#' Wellner JA, Zhan Y (1997). A hybrid algorithm for computation of the nonparametric maximum
#'   likelihood estimator from censored data. Journal of the American Statistical Association
#'   92:945-959.
#'
#' Finkelstein DM (1986). A proportional hazards model for interval-censored failure time data.
#'   Biometrics 42:845-854.
#'
#' Fay MP, Shaw PA (2010). Exact and asymptotic weighted logrank tests for interval censored data:
#'   the interval R package. Journal of Statistical Software 36(2):1-34.
#' @seealso
#' \code{\link{Interval}}, \code{\link{parsurv}} (parametric models for the same data; its
#' \code{plot()} draws this estimate for comparison), \code{\link{km}}.
#' @examples
#' set.seed(2)
#' n <- 150
#' g <- factor(rep(c("A", "B"), each = n / 2))
#' t <- rweibull(n, 1.5, ifelse(g == "A", 10, 7))
#' visit <- 3 * runif(n)                     # first visit; then every 3 time units
#' L <- R <- numeric(n)
#' for (i in seq_len(n)) {
#'   v <- visit[i] + 3 * (0:5)
#'   L[i] <- max(c(0, v[v < t[i]]))
#'   R[i] <- if (any(v >= t[i])) min(v[v >= t[i]]) else Inf
#' }
#' d <- data.frame(L, R, g)
#' fit <- turnbull(Interval(L, R) ~ g, data = d, B = 50)
#' fit
#' plot(fit)
#' predict(fit, times = c(5, 10))
#' @export
turnbull <- function(formula, data, times = NULL, B = 200L, level = 0.95, test = TRUE, seed = 1L,
                     max_iter = 100000L, tol = 1e-10) {
  seed <- .resolve_seed(seed)
  cl <- match.call()
  if (!inherits(formula, "formula")) stop("formula must be a formula, e.g. Interval(left, right) ~ group", call. = FALSE)
  if (missing(data)) stop("data is required", call. = FALSE)
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  y <- mf[[1L]]
  if (!inherits(y, "SurvXInterval"))
    stop("the left-hand side of the formula must be Interval(left, right); for exact or right-censored times use km()", call. = FALSE)
  y <- unclass(y)
  L <- as.double(y[, 1L]); R <- as.double(y[, 2L])
  if (any(y[, 3L] > 1)) stop("turnbull() is for one event; the cause is ignored only when it is 1", call. = FALSE)
  if (ncol(mf) > 1L) {
    vars <- names(mf)[-1L]
    g <- if (length(vars) == 1L) as.factor(mf[[2L]]) else interaction(mf[-1L], sep = ", ", drop = TRUE)
    g <- droplevels(g)
    variable <- paste(vars, collapse = ", ")
  } else {
    g <- factor(rep("all", length(L))); variable <- ""
  }
  keep <- !(L == 0 & !is.finite(R))
  L <- L[keep]; R <- R[keep]; g <- droplevels(g[keep])
  if (!length(L)) stop("no informative rows: every subject has left = 0 and right = Inf", call. = FALSE)
  B <- as.integer(B)
  if (is.na(B) || B < 0L) stop("B must be a whole number >= 0", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) stop("level must be in (0, 1)", call. = FALSE)
  levs <- levels(g)
  fits <- lapply(levs, function(lv) {
    i <- g == lv
    f <- .tb_fit(L[i], R[i], NULL, max_iter, tol)
    f$L <- L[i]; f$R <- R[i]
    f
  })
  names(fits) <- levs
  for (lv in levs) if (!fits[[lv]]$converged)
    warning("turnbull(): the algorithm did not converge for ", if (variable == "") "the data" else paste0("group ", lv),
            " (Kuhn-Tucker gap ", signif(fits[[lv]]$kkt, 2), "); increase max_iter", call. = FALSE)
  fin <- c(L, R[is.finite(R)])
  tmax <- max(fin)
  ptimes <- if (length(times)) as.double(times[is.finite(times)]) else .default_times(tmax)
  grid <- sort(unique(c(.sx_seq(0, tmax, 201L), ptimes)))
  boot <- NULL
  if (B > 0L) {
    # resamples of the subjects of every group from R's generator (group k: set.seed(seed + k - 1))
    boot <- lapply(seq_along(levs), function(k) {
      f <- fits[[levs[k]]]; n <- length(f$L); ng <- length(grid)
      r <- .with_seed(as.integer((as.double(seed) + k - 1) %% .Machine$integer.max), .C("SX_tb_bootstrap", as.integer(n), as.double(f$L), as.double(f$R), as.integer(ng), as.double(grid),
              as.integer(B), as.integer(max_iter), 1e-8, out = double(B * ng), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX"))
      .sx_rc(r$rc, "Turnbull bootstrap")
      out <- .matrix_rm(r$out, B, ng)
      out[is.nan(out)] <- NA_real_
      out
    })
    names(boot) <- levs
  }
  tst <- if (isTRUE(test) && length(levs) > 1L) .tb_test(L, R, g, max_iter, tol) else NULL
  TurnbullFit(model_call = cl, groups = levs, variable = variable, fits = fits, times = grid, boot = boot, B = B,
              level = level, test = tst, n = length(L), print_times = ptimes)
}

# NPMLE on the innermost intervals (compiled); w: case weights or NULL
.tb_fit <- function(L, R, w = NULL, max_iter = 100000L, tol = 1e-10) {
  n <- length(L)
  r <- .C("SX_turnbull", as.integer(n), as.double(L), as.double(R), as.integer(!is.null(w)),
          if (is.null(w)) double(1) else as.double(w), as.integer(max_iter), as.double(tol),
          n_int = integer(1), lo = double(n), hi = double(n), mass = double(n), a = integer(n), b = integer(n),
          info_d = double(2), info_i = integer(2), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "Turnbull estimate")
  J <- r$n_int
  list(lo = r$lo[seq_len(J)], hi = r$hi[seq_len(J)], p = r$mass[seq_len(J)], a = r$a + 1L, b = r$b + 1L,
       loglik = r$info_d[1L], kkt = r$info_d[2L], iterations = r$info_i[1L], converged = r$info_i[2L] == 1L)
}

# survival P(T > t) (left = FALSE) or P(T >= t) (left = TRUE): estimate with the mass of an innermost
# interval spread uniformly over it, and the lower and upper limits over all spreads (C++ core)
.tb_surv <- function(tb, t, left = FALSE) {
  J <- length(tb$lo); nt <- length(t)
  r <- .C("SX_tb_surv", as.integer(J), as.double(tb$lo), as.double(tb$hi), as.double(tb$p), as.integer(nt), as.double(t),
          as.integer(left), est = double(nt), lower = double(nt), upper = double(nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "Turnbull survival")
  list(estimate = ifelse(is.nan(r$est), NA_real_, r$est), lower = r$lower, upper = r$upper)
}

# the curve as a polyline: at every end point the left limit and the value
.tb_steps <- function(tb) {
  J <- length(tb$lo); mx <- 2L * (2L * J + 1L)
  r <- .C("SX_tb_steps", as.integer(J), as.double(tb$lo), as.double(tb$hi), as.double(tb$p), time = double(mx), surv = double(mx),
          m = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "Turnbull curve")
  k <- seq_len(r$m)
  data.frame(time = r$time[k], surv = r$surv[k])
}

# weighted logrank test, Finkelstein (1986) scores, permutation central limit variance (Fay & Shaw 2010)
.tb_test <- function(L, R, g, max_iter, tol) {
  G <- nlevels(g)
  r <- .C("SX_tb_test", length(L), as.double(L), as.double(R), as.integer(g) - 1L, as.integer(G), as.integer(max_iter),
          as.double(tol), stat = double(1), df = integer(1), p = double(1), size = double(G), U = double(G), sd = double(G),
          rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "Turnbull test")
  list(statistic = r$stat, df = r$df, p.value = r$p,
       table = data.frame(group = levels(g), n = as.integer(r$size), score = r$U, sd = r$sd))
}

.ic_count_line <- function(type) {
  ct <- tabulate(factor(type, levels = 0:3), 4L)
  paste0("  Subjects: ", length(type), "   interval-censored: ", ct[3L], "   left-censored: ", ct[4L],
         "   right-censored: ", ct[1L], "   exact: ", ct[2L])
}

.tb_table <- function(x, times) {
  out <- NULL
  for (lv in x@groups) {
    f <- x@fits[[lv]]
    s <- .tb_surv(f, times)
    lo <- hi <- rep(NA_real_, length(times))
    if (!is.null(x@boot)) {
      bm <- x@boot[[lv]]
      r <- .C("SX_grid_boot_limits", length(x@times), as.double(x@times), nrow(bm), .rm(bm), length(times), as.double(times),
              as.double(x@level), 10L, lower = double(length(times)), upper = double(length(times)), rc = integer(1),
              NAOK = TRUE, PACKAGE = "SurvX")
      .sx_rc(r$rc, "Turnbull bootstrap limits")
      lo <- ifelse(is.nan(r$lower), NA_real_, r$lower); hi <- ifelse(is.nan(r$upper), NA_real_, r$upper)
    }
    out <- rbind(out, data.frame(group = lv, time = times, survival = s$estimate, lower = lo, upper = hi,
                                 bound_low = s$lower, bound_high = s$upper, stringsAsFactors = FALSE))
  }
  out
}

# time at which a curve given on a grid first drops to p (linear interpolation)
.tb_quant <- function(time, surv, p = 0.5) {
  r <- .C("SX_tb_quantile", length(time), as.double(time), as.double(surv), as.double(p), out = double(1), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "Turnbull quantile")
  if (is.nan(r$out)) NA_real_ else r$out
}

.tb_medians <- function(x) {
  do.call(rbind, lapply(x@groups, function(lv) {
    st <- .tb_steps(x@fits[[lv]])
    md <- .tb_quant(st$time, st$surv)
    ci <- c(NA_real_, NA_real_)
    if (!is.null(x@boot)) {
      bm <- x@boot[[lv]]
      r <- .C("SX_tb_median_limits", length(x@times), as.double(x@times), nrow(bm), .rm(bm), as.double(x@level),
              lower = double(1), upper = double(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
      .sx_rc(r$rc, "Turnbull median")
      ci <- ifelse(is.nan(c(r$lower, r$upper)), NA_real_, c(r$lower, r$upper))
    }
    data.frame(group = lv, median = md, lower = ci[1L], upper = ci[2L], stringsAsFactors = FALSE)
  }))
}

.tb_lines <- function(x, detail = FALSE) {
  types <- unlist(lapply(x@fits, function(f) .ic_type(f$L, f$R)))
  grp <- x@variable != ""
  out <- c("Turnbull non-parametric estimate (NPMLE) of the survival function, interval-censored data",
           .ic_count_line(types),
           if (grp) paste0("  Groups (", x@variable, "): ", paste0(x@groups, " (n = ", vapply(x@fits, function(f) length(f$L), 1L), ")", collapse = ", ")),
           "")
  tmax <- max(vapply(x@fits, function(f) max(c(f$L, f$R[is.finite(f$R)])), 1))
  tt <- x@print_times
  tab <- .tb_table(x, tt)
  ci_lab <- if (!is.null(x@boot)) paste0(.pct(x@level), " CI (bootstrap, B = ", x@B, ")") else NULL
  shown <- data.frame(time = .num(tab$time, 3L), survival = .num(tab$survival, 3L), stringsAsFactors = FALSE)
  if (grp) shown <- cbind(group = tab$group, shown, stringsAsFactors = FALSE)
  if (!is.null(ci_lab)) shown[[ci_lab]] <- ifelse(is.finite(tab$lower), .ci_text(tab$lower, tab$upper, 3L), "-")
  shown[["identified range"]] <- ifelse(.sx_ew("absdiff", tab$bound_high, tab$bound_low) < 1e-10, "exact",
                                        .ci_text(tab$bound_low, tab$bound_high, 3L))
  out <- c(out, "Survival probability", .table_lines(shown), "")
  md <- .tb_medians(x)
  ml <- paste0(if (grp) paste0(md$group, ": ") else "", ifelse(is.finite(md$median), .num(md$median, 2L), "not reached"),
               ifelse(is.finite(md$lower), paste0(" (", .pct(x@level), " CI ", .ci_text(md$lower, md$upper, 2L), ")"), ""))
  out <- c(out, paste0("Median time to the event: ", paste(ml, collapse = "; ")))
  if (!is.null(x@test)) {
    ts <- x@test
    out <- c(out, "", "Test of equal survival (logrank scores for interval-censored data, Finkelstein 1986;",
             "permutational variance, Fay & Shaw 2010):",
             paste0("  chi-square = ", .num(ts$statistic, 2L), " on ", ts$df, " df, p = ", .pval(ts$p.value)))
    tb <- ts$table
    hi <- tb$group[which.max(tb$score)]; lo <- tb$group[which.min(tb$score)]
    out <- c(out, if (ts$p.value < 0.05)
      paste0("  Reading: the survival curves differ; events happen earliest in ", x@variable, " = ", hi,
             " and latest in ", x@variable, " = ", lo, ".")
      else "  Reading: no evidence that the survival curves differ.")
    if (detail) out <- c(out, .table_lines(data.frame(group = tb$group, n = tb$n, score = .num(tb$score, 2L),
                                                      sd = .num(tb$sd, 2L), stringsAsFactors = FALSE), 4L))
  }
  nint <- vapply(x@fits, function(f) sum(f$p > 0 & f$lo < f$hi), 1L)
  out <- c(out, "", strwrap(paste0("Note: the NPMLE puts its probability on ", sum(vapply(x@fits, function(f) sum(f$p > 0), 1L)),
                                   " innermost (Turnbull) intervals or exact times; within the ", sum(nint),
                                   " intervals of positive length the data do not say how the probability is spread ",
                                   "(the 'identified range' column), and the curve is drawn linearly there."), width = 100L, prefix = "", initial = ""))
  conv <- vapply(x@fits, function(f) f$converged, TRUE)
  out <- c(out, paste0("Algorithm: EM-ICM hybrid, ", if (all(conv)) "converged" else "NOT converged", " (",
                       paste0(vapply(x@fits, function(f) f$iterations, 1L), collapse = ", "), " iterations; Kuhn-Tucker gap ",
                       format(signif(max(vapply(x@fits, function(f) f$kkt, 1)), 2)), ")."))
  if (detail) {
    for (lv in x@groups) {
      f <- x@fits[[lv]]; k <- f$p > 0
      out <- c(out, "", paste0("Innermost intervals", if (grp) paste0(" (", lv, ")"), " and their probabilities:"),
               .table_lines(data.frame(interval = ifelse(f$lo[k] == f$hi[k], .num(f$lo[k], 3L),
                                                         paste0("(", .num(f$lo[k], 3L), ", ", .num(f$hi[k], 3L), "]")),
                                       probability = .num(f$p[k], 4L), stringsAsFactors = FALSE), 4L))
    }
  }
  c(out, "", .method_lines(.method_keys(x)))
}

S7::method(print, TurnbullFit) <- function(x, ...) { cat(.tb_lines(x, FALSE), sep = "\n"); invisible(x) }

S7::method(summary, TurnbullFit) <- function(object, ...) {
  FitSummary(model_call = object@model_call, coefficients = .tb_table(object, object@times), tables = list(),
             text = .tb_lines(object, TRUE))
}

S7::method(predict, TurnbullFit) <- function(object, times, ...) {
  if (missing(times) || !is.numeric(times)) stop("times (numeric) is required", call. = FALSE)
  x <- object
  if (!is.null(x@boot) && !all(times %in% x@times)) {
    # bootstrap curves at the new times, interpolated linearly on the stored grid
    g <- sort(unique(c(x@times, times)))
    nb <- lapply(x@boot, function(bm) {
      r <- .C("SX_approx_rows", nrow(bm), length(x@times), as.double(x@times), .rm(bm), length(g), as.double(g),
              out = double(nrow(bm) * length(g)), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
      .sx_rc(r$rc, "interpolation")
      out <- .matrix_rm(r$out, nrow(bm), length(g))
      out[is.nan(out)] <- NA_real_
      out
    })
    x@boot <- nb
    x@times <- g
  }
  out <- .tb_table(x, times)
  names(out)[names(out) == "bound_low"] <- "identified_low"
  names(out)[names(out) == "bound_high"] <- "identified_high"
  if (x@variable == "") out$group <- NULL
  out
}

S7::method(as.data.frame, TurnbullFit) <- function(x, row.names = NULL, optional = FALSE, what = c("survival", "intervals", "test"), ...) {
  what <- match.arg(what)
  switch(what,
         survival = .tb_table(x, x@times),
         intervals = do.call(rbind, lapply(x@groups, function(lv) {
           f <- x@fits[[lv]]
           data.frame(group = lv, left = f$lo, right = f$hi, probability = f$p, stringsAsFactors = FALSE)
         })),
         test = if (is.null(x@test)) NULL else cbind(x@test$table, statistic = x@test$statistic, df = x@test$df, p.value = x@test$p.value))
}

S7::method(plot, TurnbullFit) <- function(x, conf.int = !is.null(x@boot), col = NULL, lwd = 2, lty = 1, legend = TRUE,
                                          xlab = "Time", ylab = "Survival probability", main = NULL, xlim = NULL, ylim = c(0, 1), ...) {
  G <- length(x@groups)
  col <- rep_len(col %||% .okabe_ito, G)
  lty <- rep_len(lty, G)
  st <- lapply(x@fits, .tb_steps)
  xr <- xlim %||% c(0, max(unlist(lapply(st, `[[`, "time"))))
  graphics::plot(NA, xlim = xr, ylim = ylim, xlab = xlab, ylab = ylab, main = main %||% "Turnbull estimate (interval-censored data)", ...)
  if (isTRUE(conf.int) && !is.null(x@boot)) {
    a <- (1 - x@level) / 2
    for (k in seq_len(G)) {
      bm <- x@boot[[k]]
      r <- .C("SX_column_limits", nrow(bm), ncol(bm), .rm(bm), as.double(x@level), 10L, lower = double(ncol(bm)),
              upper = double(ncol(bm)), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
      .sx_rc(r$rc, "bootstrap band")
      .band(x@times, ifelse(is.nan(r$lower), NA_real_, r$lower), ifelse(is.nan(r$upper), NA_real_, r$upper), col[k])
    }
  }
  for (k in seq_len(G)) graphics::lines(st[[k]]$time, st[[k]]$surv, col = col[k], lwd = lwd, lty = lty[k])
  if (isTRUE(legend) && G > 1L)
    graphics::legend("topright", legend = paste0(x@variable, " = ", x@groups), col = col, lwd = lwd, lty = lty, bty = "n")
  invisible(x)
}
