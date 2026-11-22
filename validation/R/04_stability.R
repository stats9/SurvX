# Part 4. Stability: (a) difficult and invalid data, where every function must either give a finite,
# flagged result or stop with a clear message (never crash R, never return silent NaN);
# (b) invariances that every correct implementation satisfies exactly.

# outcome of one call: "ok" (finite result), "flagged" (finite but not converged / warned),
# "clear error" (stopped with a message), "NON-FINITE" (silent NaN / Inf: a defect)
classify <- function(x) {
  if (failed(x)) return(list(outcome = "clear error", message = fail_msg(x)))
  w <- attr(x, "warnings")
  conv <- tryCatch(x@converged, error = function(e) NA)
  fin <- tryCatch(all(is.finite(coef(x))) && all(is.finite(diag(as.matrix(vcov(x))))), error = function(e) NA)
  if (is.list(tryCatch(coef(x), error = function(e) NULL)))
    fin <- all(vapply(coef(x), function(v) all(is.finite(v)), logical(1)))
  if (isFALSE(fin) && !isFALSE(conv)) return(list(outcome = "NON-FINITE", message = "non-finite estimates without a flag"))
  if (isFALSE(conv) || length(w)) return(list(outcome = "flagged", message = paste(c(if (isFALSE(conv)) "not converged", w), collapse = "; ")))
  list(outcome = "ok", message = "")
}

stab_rec <- function(case, func, x, expected = "any") {
  cl <- classify(x)
  rec_add("stability", data.frame(kind = "difficult data", case = case, func = func, outcome = cl$outcome,
                                  expected = expected, message = substr(cl$message, 1, 160),
                                  status = if (cl$outcome == "NON-FINITE") "FAIL"
                                           else if (expected == "error" && cl$outcome != "clear error") "FAIL"
                                           else "PASS", stringsAsFactors = FALSE))
}
inv_rec <- function(case, func, diff, tol = 1e-7) {
  rec_add("stability", data.frame(kind = "invariance", case = case, func = func,
                                  outcome = if (is.na(diff)) "ERROR" else formatC(diff, format = "e", digits = 2),
                                  expected = paste("<=", tol), message = "",
                                  status = if (is.na(diff)) "FAIL" else if (diff <= tol) "PASS" else "FAIL", stringsAsFactors = FALSE))
}

# every regression function on one data set
fit_all <- function(d, cr = FALSE) {
  if (!cr) list(
    cox = attempt(cox(Event(time, status) ~ x1 + x2, d)),
    parsurv_weibull = attempt(parsurv(Event(time, status) ~ x1 + x2, d)),
    parsurv_lognormal = attempt(parsurv(Event(time, status) ~ x1 + x2, d, baseline = "lognormal")),
    parsurv_gengamma = attempt(parsurv(Event(time, status) ~ x1 + x2, d, baseline = "gengamma")),
    parsurv_rp = attempt(parsurv(Event(time, status) ~ x1 + x2, d, baseline = "rp")),
    cure = attempt(cure(Event(time, status) ~ x2, d, cure = ~ x1)),
    km = attempt(km(Event(time, status) ~ x1, d)))
  else list(
    fg = attempt(fg(Event(time, status) ~ x1 + x2, d)),
    csc = attempt(csc(Event(time, status) ~ x1 + x2, d)),
    csc_weibull = attempt(csc(Event(time, status) ~ x1 + x2, d, baseline = "weibull")),
    directcif_joint = attempt(directcif(Event(time, status) ~ x1 + x2, d)),
    directcif_weighted = attempt(directcif(Event(time, status) ~ x1 + x2, d, estimation = "weighted")),
    aj = attempt(aj(Event(time, status) ~ x1, d)))
}

run_stability <- function(cfg) {
  logmsg("Part 4: stability and invariances")
  set.seed(cfg$seed + 4L)
  base <- gen_surv(400, p = 2); basecr <- gen_cr(400, p = 2)
  cases <- list(
    "n = 15" = list(gen_surv(15, p = 2), gen_cr(15, p = 2)),
    "95% censored" = list(gen_surv(400, p = 2, cens_max = 0.25), gen_cr(400, p = 2)),
    "heavy ties (5 distinct times)" = list(gen_surv(400, p = 2, ties = 5), gen_cr(400, p = 2, ties = 5)),
    "covariate scale 1e6" = list(transform(base, x2 = x2 * 1e6), transform(basecr, x2 = x2 * 1e6)),
    "covariate scale 1e-6" = list(transform(base, x2 = x2 * 1e-6), transform(basecr, x2 = x2 * 1e-6)),
    "time scale 1e5" = list(transform(base, time = time * 1e5), transform(basecr, time = time * 1e5)),
    "near separation (all events with x1 = 1)" = list(transform(base, status = status * x1), transform(basecr, status = status * x1)),
    "one event" = list(transform(base, status = as.integer(seq_along(status) == which(status == 1)[1])),
                       transform(basecr, status = ifelse(seq_along(status) == which(status == 1)[1], 1L, ifelse(status == 2L, 2L, 0L)))),
    "collinear covariates (x2 = 2 x1)" = list(transform(base, x2 = 2 * x1), transform(basecr, x2 = 2 * x1)),
    "no plateau (cure model)" = list(gen_surv(400, p = 2, cens_max = 40), gen_cr(400, p = 2)))
  progress_start(length(cases) + 1L)
  for (nm in names(cases)) {
    a <- fit_all(cases[[nm]][[1]]); for (f in names(a)) stab_rec(nm, f, a[[f]])
    b <- fit_all(cases[[nm]][[2]], cr = TRUE); for (f in names(b)) stab_rec(nm, f, b[[f]])
    progress_tick()
  }
  # invalid input: must stop with a message
  bad <- list(
    "negative time" = transform(base, time = ifelse(seq_along(time) == 1, -1, time)),
    "missing time" = transform(base, time = ifelse(seq_along(time) == 1, NA, time)),
    "status coded 3 in cox()" = transform(base, status = ifelse(seq_along(status) == 1, 3L, status)),
    "no events" = transform(base, status = 0L),
    "entry >= time" = transform(base, entry = time + 1))
  for (nm in names(bad)) {
    d <- bad[[nm]]
    expected <- if (nm == "missing time") "any" else "error"     # missing rows may be dropped (na.omit)
    lhs <- if (nm == "entry >= time") "Event(time, status, entry = entry)" else "Event(time, status)"
    stab_rec(nm, "cox", attempt(cox(fml(lhs, c("x1", "x2")), d)), expected)
    stab_rec(nm, "parsurv", attempt(parsurv(fml(lhs, c("x1", "x2")), d)), expected)
  }
  # missing covariate values: rows dropped, n reported
  dm <- base; dm$x2[1:10] <- NA
  f <- attempt(cox(Event(time, status) ~ x1 + x2, dm))
  rec_add("stability", data.frame(kind = "difficult data", case = "10 missing covariate values", func = "cox",
                                  outcome = if (failed(f)) "clear error" else paste("n =", f@n), expected = "n = 390", message = "",
                                  status = if (!failed(f) && f@n == 390L) "PASS" else "FAIL", stringsAsFactors = FALSE))

  # ---- invariances ----
  md <- function(a, b) tryCatch({
    a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b))
    if (length(a) != length(b) || any(xor(is.na(a), is.na(b)))) NA_real_ else max(c(0, abs(a - b)), na.rm = TRUE)
  }, error = function(e) NA_real_)
  perm <- sample(nrow(base)); permcr <- sample(nrow(basecr))
  A <- fit_all(base); P <- fit_all(base[perm, ])
  for (f in setdiff(names(A), "km")) inv_rec("row order", f, md(coef(A[[f]]), coef(P[[f]])), 1e-6)
  A2 <- fit_all(basecr, TRUE); P2 <- fit_all(basecr[permcr, ], TRUE)
  for (f in setdiff(names(A2), "aj")) inv_rec("row order", f, md(coef(A2[[f]]), coef(P2[[f]])), 1e-6)
  # covariate rescaling: coefficient of x2 * 10 is coefficient / 10 (semi-parametric models; exact)
  for (f in c("cox")) inv_rec("x2 multiplied by 10", f,
                              md(coef(attempt(cox(Event(time, status) ~ x1 + x2, transform(base, x2 = 10 * x2))))[["x2"]] * 10, coef(A$cox)[["x2"]]), 1e-7)
  inv_rec("x2 multiplied by 10", "fg", md(coef(attempt(fg(Event(time, status) ~ x1 + x2, transform(basecr, x2 = 10 * x2))))[["x2"]] * 10,
                                          coef(A2$fg)[["x2"]]), 1e-7)
  # time rescaling: Cox, Fine-Gray and cause-specific coefficients unchanged; Weibull PH covariate
  # coefficients unchanged
  inv_rec("time multiplied by 1000", "cox", md(coef(attempt(cox(Event(time, status) ~ x1 + x2, transform(base, time = 1000 * time)))), coef(A$cox)), 1e-8)
  inv_rec("time multiplied by 1000", "fg", md(coef(attempt(fg(Event(time, status) ~ x1 + x2, transform(basecr, time = 1000 * time)))), coef(A2$fg)), 1e-8)
  inv_rec("time multiplied by 1000", "csc", md(coef(attempt(csc(Event(time, status) ~ x1 + x2, transform(basecr, time = 1000 * time)))), coef(A2$csc)), 1e-8)
  w <- attempt(parsurv(Event(time, status) ~ x1 + x2, transform(base, time = 1000 * time)))
  inv_rec("time multiplied by 1000", "parsurv weibull (covariate coefficients)", md(coef(w)[c("x1", "x2")], coef(A$parsurv_weibull)[c("x1", "x2")]), 1e-5)
  # duplicated data: same coefficients, half the variance
  # duplicated data (Breslow ties: duplicating creates tied times, and the Efron approximation is by
  # design not invariant to that)
  dd <- rbind(base, base)
  f1 <- attempt(cox(Event(time, status) ~ x1 + x2, base, ties = "breslow"))
  f2 <- attempt(cox(Event(time, status) ~ x1 + x2, dd, ties = "breslow"))
  inv_rec("data duplicated: coefficients (Breslow)", "cox", md(coef(f2), coef(f1)), 1e-7)
  inv_rec("data duplicated: variance halved (Breslow)", "cox", md(2 * vcov(f2), vcov(f1)), 1e-7)
  inv_rec("Efron martingale residuals sum to zero", "residuals(cox)", abs(sum(residuals(attempt(cox(Event(time, status) ~ x1 + x2,
                                                                                 gen_surv(400, ties = 8)))))), 1e-8)
  # case weights 2 = duplicated rows (Breslow ties, where the two are exactly equivalent)
  fw <- attempt(cox(Event(time, status) ~ x1 + x2, base, ties = "breslow", weights = rep(2, nrow(base))))
  fd <- attempt(cox(Event(time, status) ~ x1 + x2, dd, ties = "breslow"))
  inv_rec("weights 2 vs duplicated rows (Breslow)", "cox", md(coef(fw), coef(fd)), 1e-7)
  # one stratum = no strata; one row per subject in (start, stop] form = ordinary data
  inv_rec("strata with one level", "cox", md(coef(attempt(cox(Event(time, status) ~ x1 + x2 + strata(one), transform(base, one = 1)))), coef(A$cox)), 1e-8)
  inv_rec("entry = 0 for everyone", "cox", md(coef(attempt(cox(Event(time, status, entry = z0) ~ x1 + x2, transform(base, z0 = 0)))), coef(A$cox)), 1e-8)
  inv_rec("entry = 0 for everyone", "fg", md(coef(attempt(fg(Event(time, status, entry = z0) ~ x1 + x2, transform(basecr, z0 = 0)))), coef(A2$fg)), 1e-6)
  # reproducibility: the same call twice gives identical results (seeded bootstrap / folds)
  c1 <- attempt(cif_data(A2$csc, B = 20, seed = 5)); c2 <- attempt(cif_data(A2$csc, B = 20, seed = 5))
  inv_rec("same seed twice", "cif_data (bootstrap)", md(c1$lower, c2$lower), 0)
  m1 <- attempt(compare_models(Event(time, status) ~ x1 + x2, basecr, times = event_quantiles(basecr)[2], seed = 3))
  m2 <- attempt(compare_models(Event(time, status) ~ x1 + x2, basecr, times = event_quantiles(basecr)[2], seed = 3))
  inv_rec("same seed twice", "compare_models (cross-validation)", md(as.data.frame(m1)[sapply(as.data.frame(m1), is.numeric)],
                                                                      as.data.frame(m2)[sapply(as.data.frame(m2), is.numeric)]), 0)
  progress_tick()
  st <- rec_get("stability")
  logmsg("  ", nrow(st), " checks: ", sum(st$status == "PASS"), " PASS, ", sum(st$status == "FAIL"), " FAIL")
  invisible(st)
}
