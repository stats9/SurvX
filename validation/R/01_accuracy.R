# Part 1. Accuracy: every SurvX estimate against the same quantity from an established
# implementation, on simulated designs and on real data sets. Each comparison is one row of
# accuracy.csv with the largest absolute and relative difference and PASS / FAIL.

# ---- non-parametric estimators -----------------------------------------------------------------------
acc_nonparametric <- function(d, label, group = "x1") {
  if (!has("survival")) return(invisible())
  d$grp <- d[[group]]
  ev <- as.integer(d$status > 0)
  # Kaplan-Meier: survival, standard error of log S, median
  o <- attempt(as.data.frame(km(Event(time, ev) ~ grp, cbind(d, ev = ev))))
  r <- attempt(survival::survfit(survival::Surv(time, ev) ~ grp, data = cbind(d, ev = ev)))
  if (!failed(o) && !failed(r)) {
    rs <- summary(r, censored = FALSE)
    key_r <- paste(sub("^grp=", "", as.character(rs$strata)), rs$time)
    key_o <- paste(o$group, o$time)
    m <- match(key_r, key_o)
    check("non-parametric", "km", "survival::survfit", label, "survival at event times", o$surv[m], rs$surv)
    check("non-parametric", "km", "survival::survfit", label, "std. error of log S (Greenwood)",
          ifelse(o$surv[m] > 0, o$se[m] / o$surv[m], NA), ifelse(rs$surv > 0, rs$std.err / rs$surv, NA))
    check("non-parametric", "km", "survival::survfit", label, "median and 95% CI",
          attempt(as.matrix(km(Event(time, ev) ~ grp, cbind(d, ev = ev))@counts[, c("median", "median.lower", "median.upper")])),
          attempt(unname(summary(r)$table[, c("median", "0.95LCL", "0.95UCL")])))
  }
  # log-rank and Fleming-Harrington G-rho tests
  for (rho in c(0, 1)) {
    check("non-parametric", "logrank_test", "survival::survdiff", label, paste0("chi-square, rho = ", rho),
          attempt(logrank_test(d$time, ev, d$grp, rho = rho)$statistic),
          attempt(survival::survdiff(survival::Surv(time, ev) ~ grp, data = cbind(d, ev = ev), rho = rho)$chisq))
  }
  check("non-parametric", "logrank_test", "survival::survdiff", label, "stratified chi-square",
        attempt(logrank_test(d$time, ev, d$grp, strata = d$g)$statistic),
        attempt(survival::survdiff(survival::Surv(time, ev) ~ grp + strata(g), data = cbind(d, ev = ev))$chisq))
  # restricted mean survival time
  tau <- unname(stats::quantile(d$time, 0.8))
  o <- attempt(rmst(d$time, ev, tau))
  r <- attempt(summary(survival::survfit(survival::Surv(time, ev) ~ 1, data = cbind(d, ev = ev)), rmean = tau)$table[c("rmean", "se(rmean)")])
  check("non-parametric", "rmst", "survival::survfit (rmean)", label, "RMST", attempt(o[["rmst"]]), attempt(r[[1]]))
  check("non-parametric", "rmst", "survival::survfit (rmean)", label, "standard error of RMST", attempt(o[["se"]]), attempt(r[[2]]),
        note = "survival uses the same Irwin/Kaplan-Meier variance")
  if (has("survRM2")) {
    r2 <- attempt(survRM2::rmst2(d$time, ev, as.integer(d$grp == d$grp[1]), tau = tau))
    check("non-parametric", "rmst", "survRM2::rmst2", label, "RMST of one arm",
          attempt(rmst(d$time[d$grp == d$grp[1]], ev[d$grp == d$grp[1]], tau)[["rmst"]]),
          attempt(r2$RMST.arm1$rmst[["Est."]]))
  }
  # Nelson-Aalen cumulative hazard
  o <- attempt(hazard_np(Event(time, ev) ~ 1, cbind(d, ev = ev), smooth = FALSE)@cumhaz[[1]])
  r <- attempt(survival::survfit(survival::Surv(time, ev) ~ 1, data = cbind(d, ev = ev), ctype = 1))
  if (!failed(o) && !failed(r)) {
    k <- r$n.event > 0
    check("non-parametric", "hazard_np", "survival::survfit (ctype = 1)", label, "Nelson-Aalen cumulative hazard",
          o$cumhaz[match(r$time[k], o$time)], r$cumhaz[k])
  }
  invisible()
}

acc_competing_np <- function(d, label, group = "x1") {
  d$grp <- d[[group]]
  K <- max(d$status)
  # cmprsk keeps times as they are; survival merges times that differ only by rounding (timefix)
  oa <- attempt(aj(Event(time, status) ~ grp, d))
  oa0 <- attempt(aj(Event(time, status, timefix = FALSE) ~ grp, d))
  if (has("cmprsk")) {
    ra <- attempt(cmprsk::cuminc(d$time, d$status, d$grp))
    if (!failed(oa0) && !failed(ra)) {
      tq <- event_quantiles(d)
      tp <- cmprsk::timepoints(ra, tq)
      tabs <- summary(oa0, times = tq)@tables
      for (k in seq_len(K)) {
        ours <- unlist(lapply(names(tabs), function(l) tabs[[l]][[paste0("cif.", k)]]))
        refs <- unlist(lapply(names(tabs), function(l) tp$est[paste(l, k), ]))
        check("non-parametric", "aj", "cmprsk::cuminc", label, paste("CIF of cause", k), ours, refs)
      }
      check("non-parametric", "gray_test", "cmprsk::cuminc", label, "Gray test statistics (all causes)",
            oa0@tests$statistic, ra$Tests[, "stat"])
    }
  }
  if (has("survival")) {
    rs <- attempt(survival::survfit(survival::Surv(time, factor(status, 0:K)) ~ grp, data = d))
    if (!failed(oa) && !failed(rs)) {
      o <- as.data.frame(oa)
      ss <- summary(rs, censored = FALSE)
      key_r <- paste(sub("^grp=", "", as.character(ss$strata)), ss$time)
      m <- match(key_r, paste(o$group, o$time))
      for (k in seq_len(K)) {
        check("non-parametric", "aj", "survival::survfit (multi-state)", label, paste("CIF of cause", k),
              o[[paste0("cif.", k)]][m], ss$pstate[, k + 1L])
        check("non-parametric", "aj", "survival::survfit (multi-state)", label, paste("std. error of CIF, cause", k),
              o[[paste0("se.", k)]][m], ss$std.err[, k + 1L], tol = "info",
              note = "SurvX: delta-method (Aalen-type) variance; survival: infinitesimal jackknife")
      }
    }
  }
  invisible()
}

# ---- Cox models --------------------------------------------------------------------------------------
cox_pair <- function(d, label, vars, variant) {
  if (!has("survival")) return(invisible())
  lhs_o <- if (variant$entry) "Event(time, status, entry = entry)" else "Event(time, status)"
  lhs_r <- if (variant$entry) "survival::Surv(entry, time, status)" else "survival::Surv(time, status)"
  extra_o <- c(if (variant$strata) "strata(g)", if (variant$cluster) "cluster(id)")
  extra_r <- c(if (variant$strata) "strata(g)")
  w <- if (variant$weights) d$w else NULL
  o <- attempt(cox(fml(lhs_o, c(vars, extra_o)), d, ties = variant$ties, weights = w))
  # survival::coxph uses the robust variance by default with non-integer case weights: ask for the
  # model-based one here (compared with cox()), and compare the robust one separately below
  r <- attempt(do.call(survival::coxph, c(list(fml(lhs_r, c(vars, extra_r)), data = d, ties = variant$ties),
                                          if (variant$weights) list(weights = w, robust = FALSE), if (variant$cluster) list(cluster = d$id))))
  if (variant$weights) {
    o_rob <- attempt(cox(fml(lhs_o, c(vars, extra_o)), d, ties = variant$ties, weights = w, robust = TRUE))
    r_rob <- attempt(survival::coxph(fml(lhs_r, c(vars, extra_r)), data = d, ties = variant$ties, weights = w))
    check("Cox", "cox", "survival::coxph", paste0(label, " [", variant$name, "]"), "robust variance (case weights)",
          attempt(vcov(o_rob)), attempt(vcov(r_rob)))
  }
  tag <- paste0(label, " [", variant$name, "]")
  check("Cox", "cox", "survival::coxph", tag, "coefficients", attempt(coef(o)), attempt(coef(r)))
  check("Cox", "cox", "survival::coxph", tag, if (variant$cluster) "robust variance" else "variance", attempt(vcov(o)), attempt(vcov(r)))
  check("Cox", "cox", "survival::coxph", tag, "log partial likelihood", attempt(as.numeric(logLik(o))), attempt(r$loglik[2]))
  if (!variant$weights) {
    check("Cox", "residuals", "survival::residuals.coxph", tag, "martingale residuals", attempt(residuals(o)), attempt(unname(residuals(r, "martingale"))))
    if (!variant$strata)
      check("Cox", "residuals", "survival::residuals.coxph", tag, "scaled Schoenfeld residuals",
            attempt(residuals(o, "scaledsch")), attempt(unname(residuals(r, "scaledsch"))))
    check("Cox", "residuals", "survival::residuals.coxph", tag, "deviance residuals", attempt(residuals(o, "deviance")),
          attempt(unname(residuals(r, "deviance"))))
  }
  if (!variant$weights && !variant$cluster) {
    nd <- d[seq_len(5), , drop = FALSE]; tq <- event_quantiles(d)
    ref_s <- attempt(t(vapply(seq_len(nrow(nd)), function(i) summary(survival::survfit(r, newdata = nd[i, , drop = FALSE]), times = tq, extend = TRUE)$surv,
                              numeric(length(tq)))))
    check("Cox", "predict(cox)", "survival::survfit.coxph", tag, "predicted survival", attempt(predict(o, nd, times = tq, type = "survival")), ref_s)
  }
  if (variant$ties == "breslow" && !variant$weights) {
    for (tr in c("identity", "log", "km", "rank")) {
      ours <- attempt(ph_test(o, transform = tr)@table$chisq)
      ref <- attempt(unname(survival::cox.zph(r, transform = tr, terms = FALSE)$table[, "chisq"]))
      check("Cox", "ph_test", "survival::cox.zph", tag, paste0("score test chi-square (", tr, ")"), ours, ref)
    }
  }
  invisible()
}

acc_cox <- function(d, label, vars) {
  d$w <- stats::runif(nrow(d), 0.5, 2)
  variants <- list(
    list(name = "Efron", ties = "efron", strata = FALSE, entry = FALSE, cluster = FALSE, weights = FALSE),
    list(name = "Breslow", ties = "breslow", strata = FALSE, entry = FALSE, cluster = FALSE, weights = FALSE),
    list(name = "case weights", ties = "efron", strata = FALSE, entry = FALSE, cluster = FALSE, weights = TRUE),
    list(name = "strata", ties = "efron", strata = TRUE, entry = FALSE, cluster = FALSE, weights = FALSE),
    list(name = "strata, Breslow", ties = "breslow", strata = TRUE, entry = FALSE, cluster = FALSE, weights = FALSE),
    list(name = "delayed entry", ties = "efron", strata = FALSE, entry = TRUE, cluster = FALSE, weights = FALSE),
    list(name = "entry + strata + cluster", ties = "efron", strata = TRUE, entry = TRUE, cluster = TRUE, weights = FALSE),
    list(name = "cluster", ties = "breslow", strata = FALSE, entry = FALSE, cluster = TRUE, weights = FALSE))
  for (v in variants) {
    if ((v$entry && is.null(d$entry)) || (v$strata && length(unique(d$g)) < 2L)) next
    cox_pair(d, label, vars, v)
  }
}

# (start, stop] data with a time-varying covariate
acc_cox_tv <- function(d, label) {
  if (!has("survival")) return(invisible())
  cuts <- stats::quantile(d$time, c(0.3, 0.6))
  rows <- lapply(seq_len(nrow(d)), function(i) {
    br <- c(0, cuts[cuts < d$time[i]], d$time[i]); k <- length(br) - 1L
    data.frame(start = br[-(k + 1L)], stop = br[-1L], event = c(rep(0L, k - 1L), d$status[i]), x1 = d$x1[i], x2 = d$x2[i], id = d$id[i])
  })
  s <- do.call(rbind, rows)
  s$xt <- s$x1 * (s$start >= cuts[1])
  o <- attempt(cox(Event(stop, event, entry = start) ~ xt + x2 + cluster(id), s))
  r <- attempt(survival::coxph(survival::Surv(start, stop, event) ~ xt + x2, data = s, cluster = id))
  check("Cox", "cox", "survival::coxph", paste0(label, " [time-varying covariate]"), "coefficients", attempt(coef(o)), attempt(coef(r)))
  check("Cox", "cox", "survival::coxph", paste0(label, " [time-varying covariate]"), "robust variance", attempt(vcov(o)), attempt(vcov(r)))
  invisible()
}

# time-varying coefficient x1 * log(t): cox(tt()) against coxph(tt = )
acc_cox_tt <- function(d, label) {
  if (!has("survival") || nrow(d) > 3000) return(invisible())     # the split data grow as n x events
  for (ties in c("efron", "breslow")) {
    o <- attempt(cox(Event(time, status) ~ x1 + x2 + tt(x1), d, ties = ties, tt = function(t) log(t)))
    r <- attempt(survival::coxph(survival::Surv(time, status) ~ x1 + x2 + tt(x1), data = d, ties = ties,
                                 tt = function(x, t, ...) x * log(t)))
    tag <- paste0(label, " [tt(x1) = x1 log t, ", ties, "]")
    check("Cox", "cox (tt)", "survival::coxph (tt)", tag, "coefficients", attempt(coef(o)), attempt(coef(r)))
    check("Cox", "cox (tt)", "survival::coxph (tt)", tag, "variance", attempt(vcov(o)), attempt(vcov(r)))
  }
  invisible()
}

# ---- competing risks regression ----------------------------------------------------------------------
acc_fg <- function(d, label, vars) {
  X <- stats::model.matrix(fml("~", vars), d)[, -1L, drop = FALSE]
  o <- attempt(fg(fml("Event(time, status)", vars), d))
  o0 <- attempt(fg(fml("Event(time, status, timefix = FALSE)", vars), d))   # as crr: times kept as they are
  tq <- event_quantiles(d)
  if (has("cmprsk")) {
    # crr stops at a gradient tolerance of 1e-6 by default; a tight tolerance makes the comparison exact
    r <- attempt(cmprsk::crr(d$time, d$status, X, failcode = 1, cencode = 0, gtol = 1e-12, maxiter = 100))
    check("Fine-Gray", "fg", "cmprsk::crr", label, "coefficients", attempt(coef(o0)), attempt(r$coef))
    check("Fine-Gray", "fg", "cmprsk::crr", label, "robust (Fine-Gray) variance", attempt(vcov(o0)), attempt(r$var))
    if (!failed(o0) && !failed(r)) {
      pr <- attempt(stats::predict(r, X[1:5, , drop = FALSE]))
      refp <- attempt(t(sapply(tq, function(s) pr[max(which(pr[, 1] <= s)), -1])))
      check("Fine-Gray", "predict(fg)", "cmprsk::predict.crr", label, "predicted CIF", attempt(t(predict(o0, d[1:5, ], times = tq))), refp)
    }
  }
  if (has("survival")) {
    dd <- d; dd$ev <- factor(dd$status, 0:max(dd$status))
    fd <- attempt(survival::finegray(fml("survival::Surv(time, ev)", vars), data = dd, etype = "1"))
    ctl <- survival::coxph.control(eps = 1e-12, toler.chol = 1e-13, iter.max = 100)
    tied <- anyDuplicated(d$time) > 0
    tl <- if (tied) "info" else "exact"
    nt <- if (tied) "tied times: finegray and crr order censoring and events at the same time differently" else ""
    r2 <- if (failed(fd)) fd else attempt(survival::coxph(fml("survival::Surv(fgstart, fgstop, fgstatus)", vars), data = fd, weights = fgwt,
                                                          ties = "breslow", control = ctl))
    # log partial likelihood of the reference model (the expanded, weighted Cox model) at SurvX's
    # estimate: if it is at least the reference's own maximum, both maximise the same function and
    # the reference stopped early
    ref_ll_at <- function(fd, beta, entry = FALSE) {
      fo <- fml("survival::Surv(fgstart, fgstop, fgstatus)", vars)
      tryCatch(survival::coxph(fo, data = fd, weights = fgwt, ties = "breslow", init = unname(beta),
                               control = survival::coxph.control(iter.max = 0))$loglik[2], error = function(e) NA_real_)
    }
    judge <- function(ours, ref, fd, tl0, nt0) {
      if (failed(ours) || failed(ref) || failed(fd)) return(list(tol = tl0, note = nt0))
      gain <- ref_ll_at(fd, coef(ours)) - ref$loglik[2]
      if (is.finite(gain) && gain >= -1e-9 * abs(ref$loglik[2]) && max(abs(coef(ours) - coef(ref))) > 1e-6)
        list(tol = "info", note = sprintf("same objective: at SurvX's estimate the reference's log-likelihood is higher by %.2g, the reference stopped early", gain))
      else list(tol = tl0, note = nt0)
    }
    jg <- judge(o, r2, fd, tl, nt)
    check("Fine-Gray", "fg", "survival::finegray + coxph", label, "coefficients", attempt(coef(o)), attempt(coef(r2)), tol = jg$tol, note = jg$note)
    if (!is.null(d$entry)) {
      oe <- attempt(fg(fml("Event(time, status, entry = entry)", vars), d))
      dd$rid <- seq_len(nrow(dd))
      fe <- attempt(survival::finegray(fml("survival::Surv(entry, time, ev)", c(vars, "rid")), data = dd, etype = "1", id = rid))
      re <- if (failed(fe)) fe else attempt(survival::coxph(fml("survival::Surv(fgstart, fgstop, fgstatus)", vars), data = fe, weights = fgwt,
                                                             ties = "breslow", cluster = rid, control = ctl))
      je <- judge(oe, re, fe, tl, nt)
      check("Fine-Gray", "fg (delayed entry)", "survival::finegray + coxph", label, "coefficients", attempt(coef(oe)), attempt(coef(re)), tol = je$tol, note = je$note)
      check("Fine-Gray", "fg (delayed entry)", "survival::finegray + coxph", label, "robust variance", attempt(vcov(oe)), attempt(vcov(re)), tol = tl, note = nt)
    }
  }
  invisible()
}

acc_csc <- function(d, label, vars) {
  o <- attempt(csc(fml("Event(time, status)", vars), d))
  K <- max(d$status)
  if (has("survival")) for (k in seq_len(K)) {
    r <- attempt(survival::coxph(fml(paste0("survival::Surv(time, status == ", k, ")"), vars), data = d))
    check("cause-specific", "csc", "survival::coxph", label, paste("coefficients, cause", k), attempt(coef(o)[[k]]), attempt(coef(r)))
    check("cause-specific", "csc", "survival::coxph", label, paste("variance, cause", k), attempt(vcov(o)[[k]]), attempt(vcov(r)))
  }
  if (has("riskRegression") && has("prodlim")) {
    Hist <- prodlim::Hist
    r <- attempt(riskRegression::CSC(fml("Hist(time, status)", vars), data = d))
    nd <- d[1:5, ]; tq <- event_quantiles(d)
    for (k in seq_len(K))
      check("cause-specific", "predict(csc)", "riskRegression::CSC", label, paste("predicted CIF, cause", k),
            attempt(predict(o, nd, times = tq, cause = k)), attempt(riskRegression::predictRisk(r, newdata = nd, times = tq, cause = k)))
    if (length(unique(d$g)) > 1L) {
      os <- attempt(csc(fml("Event(time, status)", c(vars, "strata(g)")), d))
      rs <- attempt(riskRegression::CSC(fml("Hist(time, status)", c(vars, "strata(g)")), data = d))
      check("cause-specific", "predict(csc, strata)", "riskRegression::CSC", label, "predicted CIF, cause 1, strata",
            attempt(predict(os, nd, times = tq, cause = 1)), attempt(riskRegression::predictRisk(rs, newdata = nd, times = tq, cause = 1)))
    }
  }
  invisible()
}

# ---- parametric models -------------------------------------------------------------------------------
acc_parametric <- function(d, label, vars) {
  if (!has("flexsurv")) return(invisible())
  entry <- !is.null(d$entry)
  lhs_o <- if (entry) "Event(time, status, entry = entry)" else "Event(time, status)"
  lhs_r <- if (entry) "survival::Surv(entry, time, status)" else "survival::Surv(time, status)"
  tq <- event_quantiles(d); nd <- d[1:5, , drop = FALSE]
  dists <- list(c("exponential", "PH", "exp"), c("weibull", "PH", "weibullPH"), c("weibull", "AFT", "weibull"),
                c("gompertz", "PH", "gompertz"), c("lognormal", "AFT", "lnorm"), c("loglogistic", "AFT", "llogis"),
                c("gengamma", "AFT", "gengamma"))
  for (ds in dists) {
    o <- attempt(parsurv(fml(lhs_o, vars), d, baseline = ds[1], scale = ds[2]))
    r <- attempt(flexsurv::flexsurvreg(fml(lhs_r, vars), data = d, dist = ds[3]))
    tag <- paste0(ds[1], " ", ds[2])
    check("parametric", "parsurv", "flexsurv::flexsurvreg", label, paste(tag, "log-likelihood"),
          attempt(as.numeric(logLik(o))), attempt(r$loglik), tol = "optimiser")
    # covariate effects and their standard errors (the covariate block of the inverse information
    # does not depend on how the baseline parameters are parametrised)
    cv <- if (failed(o) || failed(r)) character(0) else intersect(names(coef(o)), rownames(stats::vcov(r)))
    g0 <- tryCatch(as.numeric(logLik(o)) - r$loglik, error = function(e) NA_real_)
    tl_c <- if (is.finite(g0) && g0 > 1e-4) "info" else "optimiser"
    if (length(cv)) {
      check("parametric", "parsurv", "flexsurv::flexsurvreg", label, paste(tag, "covariate coefficients"),
            unname(coef(o)[cv]), unname(coef(r)[cv]), tol = tl_c)
      check("parametric", "parsurv", "flexsurv::flexsurvreg", label, paste(tag, "covariate standard errors"),
            unname(sqrt(diag(vcov(o)))[match(cv, names(coef(o)))]), unname(sqrt(diag(stats::vcov(r))[cv])), tol = tl_c)
    }
    rp <- attempt(t(vapply(seq_len(nrow(nd)), function(i) summary(r, newdata = nd[i, , drop = FALSE], t = tq, type = "survival", tidy = TRUE)$est,
                           numeric(length(tq)))))
    # when SurvX reached a clearly higher likelihood the reference stopped early: report, do not judge
    gain <- tryCatch(as.numeric(logLik(o)) - r$loglik, error = function(e) NA_real_)
    better <- is.finite(gain) && gain > 1e-4
    check("parametric", "predict(parsurv)", "flexsurv::flexsurvreg", label, paste(tag, "predicted survival"),
          attempt(predict(o, nd, times = tq)), rp, tol = if (better) "info" else "optimiser",
          note = if (better) sprintf("SurvX log-likelihood higher by %.2g: the reference did not reach the maximum", gain) else "")
  }
  # Royston-Parmar spline (df = 3: two internal knots)
  o <- attempt(parsurv(fml(lhs_o, vars), d, baseline = "rp", df = 3))
  r <- attempt(flexsurv::flexsurvspline(fml(lhs_r, vars), data = d, k = 2, scale = "hazard"))
  check("parametric", "parsurv", "flexsurv::flexsurvspline", label, "Royston-Parmar log-likelihood",
        attempt(as.numeric(logLik(o))), attempt(as.numeric(logLik(r))), tol = "optimiser")
  # AFT models against survival::survreg (no delayed entry there)
  if (!entry && has("survival")) for (ds in list(c("weibull", "weibull"), c("lognormal", "lognormal"), c("loglogistic", "loglogistic"))) {
    o <- attempt(parsurv(fml("Event(time, status)", vars), d, baseline = ds[1], scale = "AFT"))
    r <- attempt(survival::survreg(fml("survival::Surv(time, status)", vars), data = d, dist = ds[2]))
    check("parametric", "parsurv", "survival::survreg", label, paste(ds[1], "AFT coefficients"),
          attempt(unname(coef(o)[seq_len(1 + length(coef(r)) - 1)])), attempt(unname(coef(r))), tol = "optimiser")
    check("parametric", "parsurv", "survival::survreg", label, paste(ds[1], "AFT log-likelihood"),
          attempt(as.numeric(logLik(o))), attempt(as.numeric(logLik(r))), tol = "optimiser")
    k <- if (failed(r)) 0L else length(coef(r))
    if (k > 0 && !failed(o))
      check("parametric", "parsurv", "survival::survreg", label, paste(ds[1], "AFT standard errors"),
            unname(sqrt(diag(vcov(o)))[seq_len(k)]), unname(sqrt(diag(stats::vcov(r)))[seq_len(k)]), tol = "optimiser")
  }
  # cure models
  if (has("flexsurvcure") && !entry) {
    dc <- gen_surv(nrow(d), plateau = 0.35)
    for (mix in c(TRUE, FALSE)) {
      o <- attempt(cure(Event(time, status) ~ x2, dc, cure = ~ x1, type = if (mix) "mixture" else "nonmixture"))
      r <- attempt(flexsurvcure::flexsurvcure(survival::Surv(time, status) ~ x1, data = dc, dist = "weibullPH", mixture = mix,
                                              anc = list(scale = ~ x2)))
      tag <- if (mix) "mixture" else "non-mixture"
      check("parametric", "cure", "flexsurvcure", label, paste(tag, "cure log-likelihood"), attempt(as.numeric(logLik(o))), attempt(r$loglik), tol = "optimiser")
      ndc <- dc[1:5, ]; tc <- event_quantiles(dc)
      rp <- attempt(t(vapply(1:5, function(i) summary(r, newdata = ndc[i, , drop = FALSE], t = tc, type = "survival", tidy = TRUE)$est, numeric(length(tc)))))
      check("parametric", "cure", "flexsurvcure", label, paste(tag, "predicted survival"), attempt(predict(o, ndc, times = tc)), rp, tol = "optimiser")
    }
  }
  invisible()
}

# ---- prediction performance --------------------------------------------------------------------------
acc_metrics <- function(d, label, vars) {
  if (!(has("riskRegression") && has("prodlim"))) return(invisible())
  o <- attempt(fg(fml("Event(time, status)", vars), d))
  if (failed(o)) return(invisible())
  t0 <- event_quantiles(d)[2]
  risk <- drop(predict(o, times = t0))
  Hist <- prodlim::Hist
  sc <- attempt(riskRegression::Score(list(m = matrix(risk)), formula = Hist(time, status) ~ 1, data = d, times = t0, cause = 1,
                                      metrics = c("auc", "brier"), cens.model = "km", conf.int = FALSE, null.model = TRUE))
  if (failed(sc)) { check("prediction", "brier_score", "riskRegression::Score", label, "Brier score", 0, sc); return(invisible()) }
  bs <- as.data.frame(sc$Brier$score); au <- as.data.frame(sc$AUC$score)
  ob <- attempt(brier_score(d$time, d$status, risk, t0))
  check("prediction", "brier_score", "riskRegression::Score", label, "Brier score", attempt(ob[["brier"]]), bs$Brier[bs$model == "m"])
  check("prediction", "brier_score", "riskRegression::Score", label, "Brier score of the null model", attempt(ob[["brier_null"]]), bs$Brier[bs$model == "Null model"])
  check("prediction", "auc_t", "riskRegression::Score", label, "time-dependent AUC", attempt(auc_t(d$time, d$status, risk, t0)), au$AUC[au$model == "m"])
  if (has("timeROC"))
    check("prediction", "auc_t", "timeROC::timeROC (definition 2)", label, "time-dependent AUC",
          attempt(auc_t(d$time, d$status, risk, t0)),
          attempt(unname(utils::tail(timeROC::timeROC(d$time, d$status, risk, cause = 1, weighting = "marginal", times = t0, iid = FALSE)$AUC_2, 1))),
          tol = "info",
          note = "definition 2 (competing events before t are controls). Reported only: timeROC and riskRegression::Score differ on the same data; SurvX is judged against Score")
  invisible()
}

# ---- driver ------------------------------------------------------------------------------------------
run_accuracy <- function(cfg) {
  set.seed(cfg$seed + 1L)
  progress_start(cfg$acc_reps + 1L)
  designs <- list(
    standard = function() gen_surv(500, p = 3),
    ties = function() gen_surv(500, p = 3, ties = 15),
    entry_strata = function() gen_surv(600, p = 3, entry = TRUE, strata = TRUE),
    large = function() gen_surv(5000, p = 5))
  cr_designs <- list(
    cr_standard = function() gen_cr(500, p = 3),
    cr_ties = function() gen_cr(500, p = 3, ties = 15),
    cr_entry_strata = function() gen_cr(600, p = 3, entry = TRUE, strata = TRUE),
    cr_large = function() gen_cr(5000, p = 5))
  for (rep in seq_len(cfg$acc_reps)) {
    for (nm in names(designs)) {
      d <- designs[[nm]](); vars <- grep("^x", names(d), value = TRUE); lab <- paste0(nm, " #", rep)
      acc_nonparametric(d, lab)
      acc_cox(d, lab, vars)
      if (nm == "standard") { acc_cox_tv(d, lab); acc_cox_tt(d, lab) }
      acc_parametric(d, lab, vars)
    }
    for (nm in names(cr_designs)) {
      d <- cr_designs[[nm]](); vars <- grep("^x", names(d), value = TRUE); lab <- paste0(nm, " #", rep)
      acc_competing_np(d, lab)
      acc_fg(d, lab, vars)
      acc_csc(d, lab, vars)
      acc_metrics(d, lab, vars)
    }
    logmsg("  accuracy replication ", rep, " of ", cfg$acc_reps, " done")
    progress_tick()
  }
  # real data
  rd <- real_data()
  if (length(rd)) {
    acc_nonparametric(transform(rd$lung, x1 = sex, g = "A"), "lung", "x1")
    acc_cox(transform(rd$lung, g = cut(age, 3), id = seq_len(nrow(rd$lung))), "lung", c("age", "sex", "ecog"))
    acc_parametric(rd$lung, "lung", c("age", "sex", "ecog"))
    acc_cox(transform(rd$veteran, g = celltype, id = seq_len(nrow(rd$veteran))), "veteran", c("trt", "karno", "age"))
    acc_parametric(rd$veteran, "veteran", c("trt", "karno"))
    acc_cox(transform(rd$colon, g = rx, id = seq_len(nrow(rd$colon))), "colon", c("age", "nodes", "sex"))
    for (nm in c("pbc_cr", "mgus2_cr")) {
      d <- transform(rd[[nm]], g = "A")
      vars <- setdiff(names(d), c("time", "status", "g"))
      grp <- if (nm == "pbc_cr") "trt" else "sex"
      acc_competing_np(d, nm, grp)
      acc_fg(d, nm, vars)
      acc_csc(d, nm, vars)
      acc_metrics(d, nm, vars)
    }
  }
  progress_tick()
  acc <- rec_get("accuracy")
  logmsg("  ", nrow(acc), " comparisons: ", sum(acc$status == "PASS"), " PASS, ", sum(acc$status == "FAIL"), " FAIL, ",
         sum(acc$status == "ERROR"), " ERROR, ", sum(acc$status == "INFO"), " INFO (shown, no tolerance)")
  invisible(acc)
}

# the methods added after 1.0.0, as a stage of its own (seeded on its own, so that it can be run
# separately and gives the same data whether or not the first stage ran in the same session)
run_accuracy_new <- function(cfg) {
  set.seed(cfg$seed + 11L)
  acc_new_methods(cfg)
  acc <- rec_get("accuracy")
  logmsg("  ", nrow(acc), " comparisons: ", sum(acc$status == "PASS"), " PASS, ", sum(acc$status == "FAIL"), " FAIL, ",
         sum(acc$status == "ERROR"), " ERROR, ", sum(acc$status == "INFO"), " INFO (shown, no tolerance)")
  invisible(acc)
}
