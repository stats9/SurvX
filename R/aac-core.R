# Numerical summaries computed by the C++ core (src/capi_support.cpp), so that every front end of
# the core gives the same numbers. Vectorised over their arguments.

.sx_qnorm <- function(p) {
  p <- as.double(p)
  out <- .C("SX_qnorm", length(p), p, out = double(length(p)), NAOK = TRUE, PACKAGE = "SurvX")$out
  out[is.na(p)] <- NA_real_
  out
}

.sx_pnorm <- function(x) {
  x <- as.double(x)
  out <- .C("SX_pnorm", length(x), x, out = double(length(x)), NAOK = TRUE, PACKAGE = "SurvX")$out
  out[is.na(x)] <- NA_real_
  out
}

# B resamples (rows) of 1..n with replacement, from R's generator after set.seed(seed): row b is
# exactly the b-th sample.int(n, n, replace = TRUE)
.sx_resample <- function(n, B, seed) {
  r <- .with_seed(seed, .C("SX_resample", as.integer(n), as.integer(B), out = integer(n * B), rc = integer(1), PACKAGE = "SurvX"))
  .sx_rc(r$rc, "resampling")
  matrix(r$out + 1L, B, n, byrow = TRUE)
}

# upper tail of the chi-square distribution
.sx_pchisq_upper <- function(x, df) {
  n <- max(length(x), length(df))
  if (!n) return(numeric(0))
  x <- rep_len(as.double(x), n); df <- rep_len(as.double(df), n)
  out <- .C("SX_pchisq_upper", as.integer(n), x, df, out = double(n), NAOK = TRUE, PACKAGE = "SurvX")$out
  out[is.na(x) | is.na(df)] <- NA_real_
  out
}

.sx_qgamma <- function(p, shape, rate = 1) {
  n <- max(length(p), length(shape), length(rate))
  if (!n) return(numeric(0))
  .C("SX_qgamma", as.integer(n), rep_len(as.double(p), n), rep_len(as.double(shape), n), rep_len(as.double(rate), n),
     out = double(n), NAOK = TRUE, PACKAGE = "SurvX")$out
}

# sample quantiles as stats::quantile(x, probs, type, names = FALSE, na.rm = TRUE) for types 1 and 7
.sx_quantile <- function(x, probs, type = 7L) {
  x <- as.double(x); probs <- as.double(probs)
  r <- .C("SX_quantile", length(x), x, length(probs), probs, as.integer(type), out = double(length(probs)), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "quantile")
  r$out
}

# limits of cumulative incidences on the complementary log-log scale (NA where F is 0 or 1)
.cif_interval <- function(Fv, se, zq) {
  n <- length(Fv)
  r <- .C("SX_cif_interval", as.integer(n), as.double(Fv), rep_len(as.double(se), n), as.double(zq),
          lower = double(n), upper = double(n), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "confidence interval")
  lower <- upper <- Fv                                    # keeps dim and names
  storage.mode(lower) <- storage.mode(upper) <- "double"
  lower[] <- ifelse(is.nan(r$lower), NA_real_, r$lower); upper[] <- ifelse(is.nan(r$upper), NA_real_, r$upper)
  list(lower = lower, upper = upper)
}

# Brier score, Brier score without covariates, IPA, AUC and calibration error of the risks
# (n x length(times)) at every time; NA where a measure cannot be computed (C++ core)
.sx_measures <- function(time, status, risk, times, cause, groups = 10L) {
  n <- length(time); nt <- length(times)
  risk <- matrix(as.double(risk), n, nt)
  r <- .C("SX_measures", as.integer(n), as.double(time), as.integer(status), as.integer(cause), as.integer(nt), as.double(times),
          .rm(risk), as.integer(groups), brier = double(nt), brier_null = double(nt), ipa = double(nt), auc = double(nt),
          calib = double(nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "prediction measures")
  out <- cbind(brier = r$brier, brier_null = r$brier_null, ipa = r$ipa, auc = r$auc, calib = r$calib)
  out[is.nan(out)] <- NA_real_
  out
}

# folds 1..folds, balanced within strata (a vector of codes; one stratum for unstratified folds),
# from R's generator after set.seed(seed): the folds of each stratum, taken in order of first
# appearance, are sample(rep_len(1:folds, m))
.sx_folds <- function(strata, folds, seed) {
  st <- match(strata, unique(strata)) - 1L
  r <- .with_seed(seed, .C("SX_stratified_folds", length(st), as.integer(st), as.integer(folds), out = integer(length(st)), rc = integer(1),
                           PACKAGE = "SurvX"))
  .sx_rc(r$rc, "folds")
  r$out
}

# Holm-adjusted p-values (as stats::p.adjust(p, "holm"); NA stays NA)
.sx_p_holm <- function(p) {
  p <- as.double(p)
  r <- .C("SX_p_adjust_holm", length(p), p, out = double(length(p)), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "p-value adjustment")
  ifelse(is.nan(r$out), NA_real_, r$out)
}

# standard deviation and percentile limits (type 7) of every column of B x m draws, given in
# row-major order (draw after draw; for an array cells x B that is as.double(array)). finite_only:
# leave out non-finite draws (failed fits). NA where there are too few draws.
.sx_boot_summary <- function(draws, nb, m, level, finite_only = FALSE) {
  r <- .C("SX_boot_summary", as.integer(nb), as.integer(m), as.double(draws), as.integer(finite_only), as.double(level),
          sd = double(m), lower = double(m), upper = double(m), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "bootstrap summary")
  lapply(r[c("sd", "lower", "upper")], function(v) ifelse(is.nan(v), NA_real_, v))
}

# limits est * exp(-+ z se / est) of a positive quantity; NA where est or se is not finite or est <= 0
.sx_log_limits <- function(est, se, z) {
  n <- length(est)
  r <- .C("SX_log_limits", as.integer(n), as.double(est), as.double(rep_len(se, n)), as.double(z), lower = double(n), upper = double(n),
          rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "confidence limits")
  list(lower = ifelse(is.nan(r$lower), NA_real_, r$lower), upper = ifelse(is.nan(r$upper), NA_real_, r$upper))
}

# ---- toolkit of the core: every small computation on results ------------------------------------
.ew_codes <- c(exp = 1L, log = 2L, plogis = 3L, logit = 4L, qnorm = 5L, complement = 6L, neglog = 7L, sub = 8L, div = 9L,
               logratio = 10L, pvalue = 11L, exp_se = 12L, ic_mid = 13L, neglog_floor = 14L, mul = 15L, sqrt = 16L, add = 17L,
               qnorm_comp = 18L, logit_comp = 19L, exp_neg = 20L, absdiff = 21L, signif = 22L, round = 23L)

# f(x) or f(x, y) elementwise; the result keeps the dimensions and names of x; NA where an input is NA
.sx_ew <- function(op, x, y = NULL) {
  code <- .ew_codes[[op]]; n <- length(x)
  if (!is.null(y)) y <- rep_len(as.double(y), n)
  r <- .C("SX_elementwise", code, as.integer(n), as.double(x), as.integer(!is.null(y)), if (is.null(y)) double(1) else y,
          out = double(n), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "computation")
  out <- r$out
  na <- is.na(x) & !is.nan(x)
  if (!is.null(y)) na <- na | (is.na(y) & !is.nan(y))
  out[na] <- NA_real_
  attributes(out) <- attributes(x)[intersect(names(attributes(x)), c("dim", "dimnames", "names"))]
  out
}

.rd_codes <- c(sum = 0L, mean = 1L, maxabs = 2L, max = 3L, min = 4L)
.sx_reduce <- function(op, x, na.rm = FALSE) {
  r <- .C("SX_reduce", .rd_codes[[op]], length(x), as.double(x), as.integer(na.rm), out = double(1), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "computation")
  if (is.nan(r$out) && anyNA(x) && !na.rm) NA_real_ else r$out
}

# column (or row) means of a matrix, as colMeans / rowMeans
.sx_means <- function(M, by_row = FALSE, na.rm = FALSE) {
  M <- as.matrix(M); nr <- nrow(M); nc <- ncol(M)
  r <- .C("SX_means", as.integer(nr), as.integer(nc), .rm(M), as.integer(by_row), as.integer(na.rm),
          out = double(if (by_row) nr else nc), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "means")
  out <- r$out
  names(out) <- if (by_row) rownames(M) else colnames(M)
  out
}
.sx_colmeans <- function(M, na.rm = FALSE) .sx_means(M, FALSE, na.rm)
.sx_rowmeans <- function(M, na.rm = FALSE) .sx_means(M, TRUE, na.rm)

# means of x within the groups g (in order of the sorted unique values of g), NA left out
.sx_group_means <- function(x, g) {
  lv <- sort(unique(g)); gi <- match(g, lv) - 1L
  r <- .C("SX_group_means", length(x), as.integer(gi), as.double(x), length(lv), out = double(length(lv)), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "group means")
  list(group = lv, mean = r$out)
}

# running sums along the rows of M, capped at cap (NULL: not capped)
.sx_cumsum_rows <- function(M, cap = NULL) {
  M <- as.matrix(M)
  r <- .C("SX_cumsum_rows", nrow(M), ncol(M), .rm(M), as.double(if (is.null(cap)) -1 else cap), out = double(length(M)),
          rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "running sums")
  .matrix_rm(r$out, nrow(M), ncol(M), dimnames(M))
}

.sx_linpred <- function(X, beta) {
  X <- as.matrix(X)
  r <- .C("SX_linpred", nrow(X), ncol(X), .rm(X), as.double(beta), out = double(nrow(X)), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "linear predictor")
  stats::setNames(r$out, rownames(X))
}

# survival (type "survival") or cumulative incidence ("cif") from a cumulative baseline H: a vector
# of one value per time, or a rows x times matrix (one baseline per row)
.sx_ph_predict <- function(X, beta, H, type = c("survival", "cif")) {
  type <- match.arg(type); X <- as.matrix(X); nr <- nrow(X)
  byrow <- is.matrix(H); nt <- if (byrow) ncol(H) else length(H)
  r <- .C("SX_ph_predict", as.integer(type == "cif"), as.integer(nr), ncol(X), .rm(X), as.double(beta), as.integer(nt),
          if (byrow) .rm(H) else as.double(H), as.integer(byrow), out = double(nr * nt), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "prediction")
  .matrix_rm(r$out, nr, nt)
}

# curves from event counts d and numbers at risk Y: "na" (Nelson-Aalen), "na_surv" (exp(-H)), "km"
.sx_np_curve_inc <- function(d, Y, type = c("na", "na_surv", "km")) {
  mode <- match(match.arg(type), c("na", "na_surv", "km")) - 1L
  r <- .C("SX_np_curve_inc", length(d), as.double(d), as.double(Y), mode, out = double(length(d)), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "curve")
  r$out
}

# AIC and BIC (nobs NULL: BIC not computed)
.sx_ic <- function(ll, k, nobs = NULL) {
  n <- max(length(ll), length(k))
  r <- .C("SX_info_criteria", as.integer(n), as.double(rep_len(ll, n)), as.double(rep_len(k, n)), as.integer(!is.null(nobs)),
          if (is.null(nobs)) double(1) else as.double(rep_len(nobs, n)), aic = double(n), bic = double(n), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "information criteria")
  list(aic = ifelse(is.nan(r$aic), NA_real_, r$aic), bic = ifelse(is.nan(r$bic), NA_real_, r$bic))
}
.sx_aic <- function(ll, k) .sx_ic(ll, k)$aic
.sx_bic <- function(ll, k, nobs) .sx_ic(ll, k, nobs)$bic

# likelihood ratio test 2 (ll1 - ll0) on df; p NA when df < 1
.sx_lr <- function(ll1, ll0, df) {
  n <- max(length(ll1), length(ll0), length(df))
  r <- .C("SX_lr_test", as.integer(n), as.double(rep_len(ll1, n)), as.double(rep_len(ll0, n)), as.double(rep_len(df, n)),
          stat = double(n), p = double(n), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "likelihood ratio test")
  list(stat = ifelse(is.nan(r$stat), NA_real_, r$stat), p = ifelse(is.nan(r$p), NA_real_, r$p))
}

# standard errors from a covariance matrix (named as its rows)
.sx_se <- function(V) {
  V <- as.matrix(V); p <- nrow(V)
  r <- .C("SX_se_diag", as.integer(p), .rm(V), out = double(p), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "standard errors")
  out <- ifelse(is.nan(r$out), NA_real_, r$out); names(out) <- rownames(V)
  out
}

# Wald table: limits, z, p and, where expo, the exponentiated estimate and limits
.sx_wald <- function(est, se, level, expo) {
  n <- length(est)
  r <- .C("SX_wald_table", as.integer(n), as.double(est), as.double(rep_len(se, n)), as.double(.zq(level)),
          as.integer(rep_len(expo, n)), lo = double(n), hi = double(n), z = double(n), p = double(n), eest = double(n),
          elo = double(n), ehi = double(n), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "Wald table")
  lapply(r[c("lo", "hi", "z", "p", "eest", "elo", "ehi")], function(v) ifelse(is.nan(v), NA_real_, v))
}

# residuals from fitted cumulative hazards H: "coxsnell", "martingale" or "deviance"; from martingale
# residuals given as H: "deviance_from_m", "coxsnell_from_m"
.sx_hazard_resid <- function(H, status, type) {
  n <- length(H)
  r <- .C("SX_hazard_residuals", as.integer(n), as.double(H), as.integer(status),
          match(type, c("coxsnell", "martingale", "deviance", "deviance_from_m", "coxsnell_from_m")) - 1L,
          out = double(n), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "residuals")
  r$out
}

# grids as R: seq(from, to, length.out = n), seq(from, to, by) and pretty(c(from, to), n)
.sx_grid <- function(kind, from, to, by = 0, n = 0L) {
  call <- function(cap) .C("SX_grid", as.integer(kind), as.double(from), as.double(to), as.double(by), as.integer(n), as.integer(cap),
                           out = double(max(cap, 1L)), m = integer(1), rc = integer(1), PACKAGE = "SurvX")
  r <- call(0L); .sx_rc(r$rc, "grid")
  if (!r$m) return(numeric(0))
  r <- call(r$m); .sx_rc(r$rc, "grid")
  r$out[seq_len(r$m)]
}
.sx_seq <- function(from, to, length.out) .sx_grid(0L, from, to, n = length.out)
.sx_seq_by <- function(from, to, by) .sx_grid(1L, from, to, by = by)
.sx_pretty <- function(from, to, n = 5L) .sx_grid(2L, from, to, n = n)

# 1 - (1 - level) / 2, the upper probability of a two-sided interval
.sx_upper_prob <- function(level) .sx_ew("complement", .sx_ew("div", .sx_ew("complement", level), 2))

# hazard increments (rows x grid times tj) from a baseline data frame (time, cumhaz, stratum):
# dH0 exp(LP), LP a vector (one value per row) or a rows x times matrix
.sx_hazard_inc <- function(bh, rstr, tj, LP) {
  nr <- length(rstr); m <- length(tj); bymat <- is.matrix(LP)
  r <- .C("SX_hazard_increments", nrow(bh), as.double(bh$time), as.double(bh$cumhaz), as.integer(bh$stratum), as.integer(nr),
          as.integer(rstr), as.integer(m), as.double(tj), if (bymat) .rm(LP) else as.double(LP), as.integer(bymat),
          out = double(nr * m), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "hazard increments")
  .matrix_rm(r$out, nr, m)
}

# proportions as percentages (100 x)
.sx_pct <- function(x) .sx_ew("mul", x, 100)

# signif() and round() of R, in the core
.sx_signif <- function(x, digits) .sx_ew("signif", x, digits)
.sx_round <- function(x) .sx_ew("round", x)

# stop with a clear message when covariates are constant or linear combinations of the others (C++
# core: QR with column pivoting, tolerance 1e-7 as lm()); skip: columns not checked (tt() columns)
.sx_check_collinear <- function(X, what = "covariates", skip = integer(0)) {
  keep <- setdiff(seq_len(ncol(X)), skip)
  if (length(keep) < 1L || nrow(X) < 1L) return(invisible())
  Xk <- X[, keep, drop = FALSE]
  if (any(!is.finite(Xk))) return(invisible())
  r <- .C("SX_collinear", nrow(Xk), ncol(Xk), .rm(Xk), 1e-7, out = integer(ncol(Xk)), ndep = integer(1), rc = integer(1),
          NAOK = TRUE, PACKAGE = "SurvX")
  .sx_rc(r$rc, "collinearity check")
  if (!r$ndep) return(invisible())
  dep <- colnames(Xk)[r$out[seq_len(r$ndep)] + 1L]
  const <- dep[vapply(dep, function(v) length(unique(Xk[, v])) <= 1L, logical(1))]
  comb <- setdiff(dep, const)
  msg <- c(if (length(const)) paste0("constant in the data (no information on its effect): ", paste(const, collapse = ", ")),
           if (length(comb)) paste0("a linear combination of the other ", what, " (its effect cannot be separated from theirs): ",
                                    paste(comb, collapse = ", ")))
  stop(paste(msg, collapse = "; "), "; remove ", if (length(dep) == 1L) "it" else "them", " from the formula", call. = FALSE)
}
