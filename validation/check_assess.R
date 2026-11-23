# Quick check of the measures of prediction performance (cindex, standard errors and paired
# differences of the Brier score and AUC, O/E, ICI/E50/E90, decision curves, anova) before CI.
# Run from the repository root after installing dev:
#   R CMD INSTALL .
#   Rscript validation/check_assess.R
# Takes about a minute. Writes check_assess.txt: one line per comparison (PASS, FAIL, INFO or
# SKIP when a reference package is not installed) and, for every FAIL or INFO, the numbers of both
# packages, so that a disagreement can be traced at once. The C index with competing risks is also
# computed by a direct R implementation over all pairs (written from Wolbers et al. 2009, 2014),
# so it is checked even without pec; the pec section prints the variants needed to see where pec
# differs (weights, truncation time, direction of the score).
library(SurvX)
suppressPackageStartupMessages(library(survival))
sink("check_assess.txt", split = TRUE)
cat("SurvX", as.character(packageVersion("SurvX")), "|", R.version.string, "|", format(Sys.time()), "\n")
has <- function(p) requireNamespace(p, quietly = TRUE)
for (p in c("riskRegression", "prodlim", "pec", "dcurves", "splines"))
  cat(sprintf("%-15s %s\n", p, if (has(p)) as.character(packageVersion(p)) else "not installed"))
res <- data.frame(check = character(0), status = character(0), max_diff = numeric(0), stringsAsFactors = FALSE)
check <- function(label, ours, ref, tol = 1e-8, info = FALSE) {
  st <- tryCatch({
    a <- as.numeric(unlist(ours)); b <- as.numeric(unlist(ref))
    if (length(a) != length(b)) stop("lengths ", length(a), " and ", length(b))
    d <- max(abs(a - b) / pmax(1, abs(b)))
    s <- if (info) "INFO" else if (is.finite(d) && d <= tol) "PASS" else "FAIL"
    cat(sprintf("%-4s  %-78s  %.2e\n", s, label, d))
    if (s != "PASS") { cat("        SurvX:    ", format(signif(a, 10)), "\n"); cat("        reference:", format(signif(b, 10)), "\n") }
    c(s, d)
  }, error = function(e) { cat(sprintf("ERR   %-78s  %s\n", label, conditionMessage(e))); c("ERROR", NA) })
  res[nrow(res) + 1L, ] <<- list(label, st[1L], as.numeric(st[2L]))
  invisible(st[1L] == "PASS")
}
skip <- function(label, why) { cat(sprintf("SKIP  %-78s  %s\n", label, why)); res[nrow(res) + 1L, ] <<- list(label, "SKIP", NA_real_) }
section <- function(s) cat("\n=====", s, "=====\n")

# ---- data ---------------------------------------------------------------------------------------
set.seed(20261123)
sim_cr <- function(n, ties = FALSE) {
  x1 <- rbinom(n, 1, 0.5); x2 <- rnorm(n)
  t1 <- rweibull(n, 1.2, 1 / (0.35 * exp(0.6 * x1 + 0.4 * x2))^(1 / 1.2)); t2 <- rexp(n, 0.25 * exp(-0.3 * x1))
  cens <- runif(n, 0, 5); time <- pmin(t1, t2, cens)
  if (ties) time <- ceiling(time * 4) / 4
  data.frame(time = time, status = ifelse(pmin(t1, t2) > cens, 0L, ifelse(t1 < t2, 1L, 2L)), x1 = x1, x2 = x2)
}
dcr <- sim_cr(500); dcr_t <- sim_cr(300, ties = TRUE)
dsv <- dcr; dsv$status <- as.integer(dsv$status > 0)                # one event type
lung <- survival::lung; lung$status <- lung$status - 1
lung <- lung[complete.cases(lung[, c("time", "status", "age", "sex", "ph.ecog")]), ]

# ---- 1. concordance -----------------------------------------------------------------------------
section("1. cindex")
# C over all pairs from the definitions: kind A (control event-free after T_i; censored at T_i
# counts), kind B (competing failure at or before T_i); weights 1, or 1/G(T_i-)^2 (A) and
# 1/{G(T_i-) G(T_j-)} (B); G = Kaplan-Meier of the censoring times, events before censorings at ties
ref_c <- function(score, time, status, cause = 1, tau = Inf, uno = FALSE) {
  n <- length(time); G <- numeric(n); g <- 1
  for (u in sort(unique(time))) {                     # G(u-); at tied times the events leave first
    at <- time == u; G[at] <- g
    cens <- sum(at & status == 0); ev <- sum(at & status != 0); risk <- sum(time >= u)
    if (cens > 0) g <- g * (1 - cens / (risk - ev))
  }
  num <- den <- 0
  for (i in which(status == cause & time <= tau)) {
    A <- time > time[i] | (time == time[i] & status == 0)
    B <- !A & status != 0 & status != cause
    A[i] <- B[i] <- FALSE
    w <- ifelse(A, if (uno) 1 / G[i]^2 else 1, ifelse(B, if (uno) 1 / (G[i] * G) else 1, 0))
    s <- ifelse(score[i] > score, 1, ifelse(score[i] == score, 0.5, 0))
    num <- num + sum(w * s); den <- den + sum(w)
  }
  num / den
}
sc <- dcr$x1 + 0.5 * dcr$x2; sct <- round(dcr_t$x1 + 0.5 * dcr_t$x2, 1)
check("Harrell, one event: survival::concordance (C, se)",
      unlist(cindex(sc, dsv$time, dsv$status)[c("concordance", "se")]),
      { r <- survival::concordance(Surv(time, status) ~ sc, data = cbind(dsv, sc = sc), reverse = TRUE); c(r$concordance, sqrt(r$var)) })
tau <- unname(quantile(dsv$time[dsv$status == 1], 0.8))
check("Uno, one event, tau: survival::concordance(timewt = 'n/G2', ymax)",
      cindex(sc, dsv$time, dsv$status, method = "uno", tau = tau)$concordance,
      unname(survival::concordance(Surv(time, status) ~ sc, data = cbind(dsv, sc = sc), reverse = TRUE, timewt = "n/G2", ymax = tau)$concordance))
check("Harrell, one event, Cox fit: survival::coxph concordance",
      cindex(cox(Event(time, status) ~ age + sex + ph.ecog, lung))$concordance,
      unname(survival::coxph(Surv(time, status) ~ age + sex + ph.ecog, lung)$concordance["concordance"]))
for (tt in list(list(dcr, sc, "continuous"), list(dcr_t, sct, "tied times and scores"))) {
  dd <- tt[[1]]; s <- tt[[2]]
  check(paste0("Wolbers 2009 (unweighted), competing risks, ", tt[[3]], ": all pairs in R"),
        cindex(s, dd$time, dd$status)$concordance, ref_c(s, dd$time, dd$status))
  check(paste0("Wolbers 2014 (weighted), competing risks, ", tt[[3]], ": all pairs in R"),
        cindex(s, dd$time, dd$status, method = "uno")$concordance, ref_c(s, dd$time, dd$status, uno = TRUE))
  check(paste0("Wolbers 2014, tau = 2, ", tt[[3]], ": all pairs in R"),
        cindex(s, dd$time, dd$status, method = "uno", tau = 2)$concordance, ref_c(s, dd$time, dd$status, tau = 2, uno = TRUE))
}
check("without competing events, the competing-risks C is Harrell's (cause recoded)",
      cindex(sc, dsv$time, dsv$status)$concordance, cindex(sc, dsv$time, ifelse(dsv$status == 1, 1L, 0L))$concordance, tol = 0)
check("standard error: infinitesimal jackknife = sqrt(sum dfbeta^2)",
      { r <- cindex(sc, dcr$time, dcr$status, method = "uno"); r$se }, sqrt(sum(cindex(sc, dcr$time, dcr$status, method = "uno")$dfbeta^2)))
if (has("pec") && has("prodlim")) {
  Hist <- prodlim::Hist
  pc <- function(score, d, ...) {
    out <- utils::capture.output(v <- pec::cindex(list(m = matrix(score, ncol = 1L)), formula = Hist(time, status) ~ 1, data = d,
                                                  cause = 1, cens.model = "marginal", splitMethod = "none", ...))
    unname(unlist(v$AppCindex))[1L]
  }
  et <- max(dcr$time) + 1
  ours_w <- cindex(sc, dcr$time, dcr$status, method = "uno")$concordance
  ours_u <- cindex(sc, dcr$time, dcr$status)$concordance
  p1 <- tryCatch(pc(sc, dcr, eval.times = et), error = function(e) NA)
  check("Wolbers 2014: pec::cindex(cause = 1, eval.times beyond the data)", ours_w, p1, tol = 1e-8)
  cat("        diagnostics: SurvX weighted", format(ours_w, digits = 10), " unweighted", format(ours_u, digits = 10),
      " | pec", format(p1, digits = 10), " pec with -score", format(tryCatch(pc(-sc, dcr, eval.times = et), error = function(e) NA), digits = 10),
      " pec at eval.times = 2", format(tryCatch(pc(sc, dcr, eval.times = 2), error = function(e) NA), digits = 10),
      " SurvX tau = 2", format(cindex(sc, dcr$time, dcr$status, method = "uno", tau = 2)$concordance, digits = 10), "\n")
} else skip("Wolbers 2014: pec::cindex", "pec or prodlim not installed")

# ---- 2. standard errors and paired differences ---------------------------------------------------
section("2. standard errors and paired differences of the Brier score and AUC")
for (dd in list(list(dsv, "one event"), list(dcr, "competing risks"))) {
  d <- dd[[1]]; t <- median(d$time[d$status > 0])
  f1 <- if (dd[[2]] == "one event") cox(Event(time, status) ~ x1, d) else fg(Event(time, status) ~ x1, d)
  f2 <- if (dd[[2]] == "one event") cox(Event(time, status) ~ x1 + x2, d) else fg(Event(time, status) ~ x1 + x2, d)
  rk <- function(f) if (dd[[2]] == "one event") 1 - drop(predict(f, d, times = t, type = "survival")) else drop(predict(f, d, times = t))
  r1 <- rk(f1); r2 <- rk(f2)
  b2 <- brier_score(d$time, d$status, r2, t, se = TRUE); a2 <- auc_t(d$time, d$status, r2, t, se = TRUE)
  cc <- cv_compare(Event(time, status) ~ x1 + x2, d, times = t, folds = 0,
                   models = list(m1 = function(f, dd) if (all(dd$status <= 1)) cox(Event(time, status) ~ x1, dd) else fg(Event(time, status) ~ x1, dd),
                                 m2 = function(f, dd) if (all(dd$status <= 1)) cox(f, dd) else fg(f, dd)))
  tb <- cc@table
  if (has("riskRegression") && has("prodlim")) {
    Hist <- prodlim::Hist
    sc_r <- riskRegression::Score(list(m1 = matrix(r1), m2 = matrix(r2)), formula = Hist(time, status) ~ 1, data = d, times = t,
                                  cause = 1, metrics = c("auc", "brier"), null.model = FALSE, conservative = TRUE, se.fit = TRUE,
                                  cens.model = "km", contrasts = TRUE)
    A <- as.data.frame(sc_r$AUC$score); B <- as.data.frame(sc_r$Brier$score)
    ac <- as.data.frame(sc_r$AUC$contrasts); bc <- as.data.frame(sc_r$Brier$contrasts)
    check(paste0(dd[[2]], ": Brier and its se (brier_score) vs Score"), b2[c("brier", "se")], c(B$Brier[B$model == "m2"], B$se[B$model == "m2"]))
    check(paste0(dd[[2]], ": AUC and its se (auc_t) vs Score"), a2[c("auc", "se")], c(A$AUC[A$model == "m2"], A$se[A$model == "m2"]))
    sgn <- if (tb$reference[1L] == "m1") 1 else -1        # our differences are model minus the reference
    o <- tb[tb$model != tb$reference[1L], ]
    check(paste0(dd[[2]], ": paired difference m2 - m1 of the AUC, se, p vs Score"),
          c(sgn * o$auc_diff, o$auc_diff_se, o$auc_diff_p), c(ac$delta.AUC[1L], ac$se[1L], ac$p[1L]))
    check(paste0(dd[[2]], ": paired difference m2 - m1 of the Brier score, se, p vs Score"),
          c(sgn * o$brier_diff, o$brier_diff_se, o$brier_diff_p), c(bc$delta.Brier[1L], bc$se[1L], bc$p[1L]))
  } else skip(paste0(dd[[2]], ": Score"), "riskRegression or prodlim not installed")
  z <- qnorm(0.975)
  check(paste0(dd[[2]], ": intervals of the table = estimate -/+ z se"),
        with(tb[tb$model == "m2", ], c(brier_lower, brier_upper, auc_lower, auc_upper)),
        with(tb[tb$model == "m2", ], c(max(0, brier - z * brier_se), brier + z * brier_se, max(0, auc - z * auc_se), min(1, auc + z * auc_se))))
  check(paste0(dd[[2]], ": table measures = brier_score, auc_t"), with(tb[tb$model == "m2", ], c(brier, brier_se, auc, auc_se)),
        c(b2[["brier"]], b2[["se"]], a2[["auc"]], a2[["se"]]))
}

# ---- 3. calibration ---------------------------------------------------------------------------
section("3. calibration: O/E, ICI, E50, E90")
for (dd in list(list(dsv, FALSE, "one event"), list(dcr, TRUE, "competing risks"))) {
  d <- dd[[1]]; cr <- dd[[2]]; t <- median(d$time[d$status > 0])
  f <- if (cr) fg(Event(time, status) ~ x1 + x2, d) else cox(Event(time, status) ~ x1 + x2, d)
  r <- if (cr) drop(predict(f, d, times = t)) else 1 - drop(predict(f, d, times = t, type = "survival"))
  cal <- calibration(d$time, d$status, r, t)
  s <- survival::survfit(Surv(time, factor(status)) ~ 1, data = d)
  O <- summary(s, times = t, extend = TRUE)$pstate[, which(s$states == "1")]
  check(paste0(dd[[3]], ": O/E vs survival::survfit"), attr(cal, "oe")[1L], O / mean(r))
  if (has("splines")) {
    cl <- log(-log(1 - r)); kn <- quantile(cl, c(0.1, 0.5, 0.9))
    Bm <- splines::ns(cl, knots = kn[2L], Boundary.knots = kn[c(1L, 3L)]); colnames(Bm) <- c("b1", "b2")
    D <- data.frame(time = d$time, status = d$status, Bm)
    if (!cr) { cf <- coxph(Surv(time, status == 1) ~ b1 + b2, data = D, ties = "efron") }
    else {
      D$ev <- factor(D$status, 0:2, c("censor", "one", "two"))
      fgd <- survival::finegray(Surv(time, ev) ~ b1 + b2, data = D, etype = "one")
      cf <- coxph(Surv(fgstart, fgstop, fgstatus) ~ b1 + b2, data = fgd, weights = fgwt, ties = "breslow")
    }
    bh <- basehaz(cf, centered = FALSE); H <- max(c(0, bh$hazard[bh$time <= t]))
    a <- abs(1 - exp(-H * exp(drop(Bm %*% coef(cf)))) - r)
    check(paste0(dd[[3]], ": ICI, E50, E90 vs ", if (cr) "finegray + coxph" else "coxph", " + splines::ns"),
          c(attr(cal, "ici"), attr(cal, "e50"), attr(cal, "e90")), c(mean(a), median(a), unname(quantile(a, 0.9))), tol = if (cr) 1e-6 else 1e-7)
  }
}

# ---- 4. decision curve ------------------------------------------------------------------------
section("4. decision_curve")
th <- c(0.05, 0.1, 0.2, 0.3, 0.5)
for (dd in list(list(dsv, "one event"), list(dcr, "competing risks"))) {
  d <- dd[[1]]; t <- median(d$time[d$status > 0])
  r <- plogis(-1 + d$x1 + 0.5 * d$x2)
  dc <- decision_curve(list(m = r), d$time, d$status, t = t, thresholds = th)
  F_of <- function(i) { if (!any(i)) return(0); s <- survfit(Surv(d$time[i], factor(d$status[i], 0:2)) ~ 1)
                        summary(s, times = t, extend = TRUE)$pstate[, which(s$states == "1")] }
  nb <- vapply(th, function(p) { i <- r >= p; F <- F_of(i); mean(i) * F - mean(i) * (1 - F) * p / (1 - p) }, 0)
  Fa <- F_of(rep(TRUE, nrow(d)))
  check(paste0(dd[[2]], ": net benefit (model and treat all) vs survfit"),
        c(dc$table$net_benefit[dc$table$model == "m"], dc$table$net_benefit[dc$table$model == "treat all"]), c(nb, Fa - (1 - Fa) * th / (1 - th)))
  if (dd[[2]] == "one event") {
    if (has("dcurves")) {
      x <- dcurves::dca(Surv(time, status) ~ m, data = data.frame(time = d$time, status = d$status, m = r), time = t, thresholds = th)
      tb <- as.data.frame(x$dca)
      check("one event: net benefit vs dcurves::dca", c(dc$table$net_benefit[dc$table$model == "m"], dc$table$net_benefit[dc$table$model == "treat all"]),
            c(tb$net_benefit[tb$variable == "m"], tb$net_benefit[tb$variable == "all"]))
    } else skip("net benefit vs dcurves::dca", "dcurves not installed")
  }
}

# ---- 5. anova ---------------------------------------------------------------------------------
section("5. anova")
a <- anova(cox(Event(time, status) ~ age, lung), cox(Event(time, status) ~ age + sex + ph.ecog, lung))
r <- anova(coxph(Surv(time, status) ~ age, lung), coxph(Surv(time, status) ~ age + sex + ph.ecog, lung))
check("Cox: chi-square, df, p vs survival::anova.coxph", c(a$Chisq[2L], a$Df[2L], a[["Pr(>|Chi|)"]][2L]), c(r$Chisq[2L], r$Df[2L], r[["Pr(>|Chi|)"]][2L]))
a <- anova(parsurv(Event(time, status) ~ age, lung), parsurv(Event(time, status) ~ age + sex, lung))
st <- 2 * (survreg(Surv(time, status) ~ age + sex, lung)$loglik[2L] - survreg(Surv(time, status) ~ age, lung)$loglik[2L])
check("Weibull: chi-square vs survreg log-likelihoods", a$Chisq[2L], st, tol = 1e-6)
a <- anova(csc(Event(time, status) ~ x1, dcr), csc(Event(time, status) ~ x1 + x2, dcr))
ll <- function(v) sum(vapply(1:2, function(k) { dk <- dcr; dk$ev <- as.integer(dcr$status == k)
                               coxph(reformulate(v, "Surv(time, ev)"), data = dk)$loglik[2L] }, 0))
check("cause-specific Cox: chi-square vs coxph per cause", a$Chisq[2L], 2 * (ll(c("x1", "x2")) - ll("x1")))
cat("refused for Fine-Gray: ", inherits(tryCatch(anova(fg(Event(time, status) ~ x1, dcr), fg(Event(time, status) ~ x1 + x2, dcr)),
                                                 error = function(e) e), "error"), "\n")

# ---- 6. the tables and methods of the comparisons ----------------------------------------------
section("6. compare_models, model_scan: columns, print, plot")
cm <- compare_models(Event(time, status) ~ x1 + x2, dcr, times = c(1, 2), folds = 3, direct = list(c("gompertz", "PH")))
need <- c("brier_se_cv", "brier_lower_cv", "auc_diff_p_cv", "oe_cv", "ici_cv", "reference_cv")
check("compare_models: new columns present", all(need %in% names(cm@table)), TRUE, tol = 0)
for (t in c(1, 2)) {
  s <- cm@table[cm@table$time == t, ]
  check(paste0("compare_models t = ", t, ": the reference has the lowest cross-validated Brier score and difference 0"),
        c(s$brier_cv[s$model == s$reference_cv[1L]], s$brier_diff_cv[s$model == s$reference_cv[1L]]), c(min(s$brier_cv), 0))
}
invisible(utils::capture.output(print(cm)))
pdf(NULL); for (m in c("brier", "auc", "ipa", "calib_error", "ici", "oe")) plot(cm, measure = m); invisible(dev.off())
ms <- model_scan(Event(time, status) ~ x1 + x2, dsv, depth = "quick", folds = 3)
check("model_scan: iauc_cv in (0.5, 1) for every fitted model", all(ms@table$iauc_cv[ms@table$fitted] > 0.5 & ms@table$iauc_cv[ms@table$fitted] < 1, na.rm = TRUE), TRUE, tol = 0)
invisible(utils::capture.output(print(ms)))

# ---- summary ----------------------------------------------------------------------------------
section("summary")
print(table(res$status))
bad <- res[res$status %in% c("FAIL", "ERROR"), ]
if (nrow(bad)) { cat("\nNot passed:\n"); print(bad, row.names = FALSE) } else cat("Everything passed.\n")
sink()
