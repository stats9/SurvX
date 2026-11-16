# Part 2. Speed: elapsed time of one fit against the sample size n and the number of covariates p,
# for SurvX and the established implementation of the same model, on the same data. A method that
# takes longer than speed_max_seconds is not run at larger n (recorded as "skipped").

speed_methods <- function() {
  S <- list()
  add <- function(model, impl, family, fun) S[[length(S) + 1L]] <<- list(model = model, impl = impl, family = family, fun = fun)
  X_of <- function(d) as.matrix(d[, grep("^x", names(d)), drop = FALSE])
  fx <- function(lhs, d) fml(lhs, grep("^x", names(d), value = TRUE))
  # one event
  add("Kaplan-Meier", "SurvX::km", "surv", function(d) km(Event(time, status) ~ x1, d))
  if (has("survival")) add("Kaplan-Meier", "survival::survfit", "surv", function(d) survival::survfit(survival::Surv(time, status) ~ x1, data = d))
  add("Cox", "SurvX::cox", "surv", function(d) cox(fx("Event(time, status)", d), d))
  if (has("survival")) add("Cox", "survival::coxph", "surv", function(d) survival::coxph(fx("survival::Surv(time, status)", d), data = d))
  add("Cox, strata + cluster", "SurvX::cox", "surv", function(d) cox(fml("Event(time, status)", c(grep("^x", names(d), value = TRUE), "strata(g)", "cluster(id)")), d))
  if (has("survival")) add("Cox, strata + cluster", "survival::coxph", "surv",
                           function(d) survival::coxph(fml("survival::Surv(time, status)", c(grep("^x", names(d), value = TRUE), "strata(g)")), data = d, cluster = id))
  add("Weibull PH", "SurvX::parsurv", "surv", function(d) parsurv(fx("Event(time, status)", d), d))
  if (has("flexsurv")) add("Weibull PH", "flexsurv::flexsurvreg", "surv", function(d) flexsurv::flexsurvreg(fx("survival::Surv(time, status)", d), data = d, dist = "weibullPH"))
  if (has("survival")) add("Weibull AFT", "survival::survreg", "surv", function(d) survival::survreg(fx("survival::Surv(time, status)", d), data = d, dist = "weibull"))
  add("Weibull AFT", "SurvX::parsurv", "surv", function(d) parsurv(fx("Event(time, status)", d), d, scale = "AFT"))
  add("generalized gamma", "SurvX::parsurv", "surv", function(d) parsurv(fx("Event(time, status)", d), d, baseline = "gengamma"))
  if (has("flexsurv")) add("generalized gamma", "flexsurv::flexsurvreg", "surv", function(d) flexsurv::flexsurvreg(fx("survival::Surv(time, status)", d), data = d, dist = "gengamma"))
  add("Royston-Parmar (df 3)", "SurvX::parsurv", "surv", function(d) parsurv(fx("Event(time, status)", d), d, baseline = "rp", df = 3))
  if (has("flexsurv")) add("Royston-Parmar (df 3)", "flexsurv::flexsurvspline", "surv", function(d) flexsurv::flexsurvspline(fx("survival::Surv(time, status)", d), data = d, k = 2))
  add("mixture cure (Weibull)", "SurvX::cure", "cure", function(d) cure(Event(time, status) ~ x2, d, cure = ~ x1))
  if (has("flexsurvcure")) add("mixture cure (Weibull)", "flexsurvcure::flexsurvcure", "cure",
                               function(d) flexsurvcure::flexsurvcure(survival::Surv(time, status) ~ x1, data = d, dist = "weibullPH", anc = list(scale = ~ x2)))
  # competing risks
  add("Aalen-Johansen", "SurvX::aj", "cr", function(d) aj(Event(time, status) ~ x1, d))
  if (has("cmprsk")) add("Aalen-Johansen", "cmprsk::cuminc", "cr", function(d) cmprsk::cuminc(d$time, d$status, d$x1))
  if (has("survival")) add("Aalen-Johansen", "survival::survfit", "cr", function(d) survival::survfit(survival::Surv(time, factor(status, 0:2)) ~ x1, data = d))
  add("Fine-Gray", "SurvX::fg", "cr", function(d) fg(fx("Event(time, status)", d), d))
  if (has("cmprsk")) add("Fine-Gray", "cmprsk::crr", "cr", function(d) cmprsk::crr(d$time, d$status, X_of(d), failcode = 1, cencode = 0))
  if (has("survival")) add("Fine-Gray", "survival::finegray + coxph", "cr", function(d) {
    dd <- d; dd$ev <- factor(dd$status, 0:2)
    fd <- survival::finegray(fx("survival::Surv(time, ev)", dd), data = dd, etype = "1")
    survival::coxph(fx("survival::Surv(fgstart, fgstop, fgstatus)", dd), data = fd, weights = fgwt)
  })
  add("cause-specific Cox", "SurvX::csc", "cr", function(d) csc(fx("Event(time, status)", d), d))
  if (has("riskRegression") && has("prodlim")) add("cause-specific Cox", "riskRegression::CSC", "cr", function(d) {
    Hist <- prodlim::Hist; riskRegression::CSC(fx("Hist(time, status)", d), data = d)
  })
  add("CIF prediction (csc, 100 x 10)", "SurvX::predict", "cr", function(d) {
    f <- csc(fx("Event(time, status)", d), d); predict(f, d[1:100, ], times = event_quantiles(d, seq(0.1, 0.9, length.out = 10)), cause = 1)
  })
  if (has("riskRegression") && has("prodlim")) add("CIF prediction (csc, 100 x 10)", "riskRegression::predictRisk", "cr", function(d) {
    Hist <- prodlim::Hist; f <- riskRegression::CSC(fx("Hist(time, status)", d), data = d)
    riskRegression::predictRisk(f, newdata = d[1:100, ], times = event_quantiles(d, seq(0.1, 0.9, length.out = 10)), cause = 1)
  })
  add("direct CIF (joint, Gompertz PH)", "SurvX::directcif", "cr", function(d) directcif(fx("Event(time, status)", d), d))
  add("direct CIF (weighted, Gompertz PH)", "SurvX::directcif", "cr", function(d) directcif(fx("Event(time, status)", d), d, estimation = "weighted"))
  add("Brier + AUC at one time", "SurvX::brier_score + auc_t", "cr", function(d) {
    r <- d$x1 * 0.2 + stats::plogis(d$x2); t0 <- event_quantiles(d)[2]
    c(brier_score(d$time, d$status, r, t0), auc_t(d$time, d$status, r, t0))
  })
  if (has("riskRegression") && has("prodlim")) add("Brier + AUC at one time", "riskRegression::Score", "cr", function(d) {
    Hist <- prodlim::Hist; r <- d$x1 * 0.2 + stats::plogis(d$x2)
    riskRegression::Score(list(m = matrix(r)), formula = Hist(time, status) ~ 1, data = d, times = event_quantiles(d)[2], cause = 1,
                          metrics = c("auc", "brier"), conf.int = FALSE, null.model = FALSE)
  })
  S
}

run_speed <- function(cfg) {
  logmsg("Part 2: speed")
  ms <- speed_methods()
  too_slow <- character(0)
  set.seed(cfg$seed + 2L)
  # (a) sample size, p = 5
  for (n in cfg$speed_n) {
    for (rep in seq_len(cfg$speed_reps)) {
      dat <- list(surv = gen_surv(n, p = 5, strata = TRUE), cure = gen_surv(n, p = 2, plateau = 0.3), cr = gen_cr(n, p = 5))
      for (m in ms) {
        key <- paste(m$model, m$impl)
        if (key %in% too_slow) { if (rep == 1L) rec_add("speed", data.frame(design = "n", model = m$model, impl = m$impl, n = n, p = 5L, rep = rep, seconds = NA_real_, status = "skipped (too slow)")); next }
        tm <- time_it(function() m$fun(dat[[m$family]]))
        st <- if (failed(tm$value)) paste("error:", substr(fail_msg(tm$value), 1, 80)) else "ok"
        rec_add("speed", data.frame(design = "n", model = m$model, impl = m$impl, n = n, p = 5L, rep = rep, seconds = tm$seconds, status = st))
        if (is.finite(tm$seconds) && tm$seconds > cfg$speed_max_seconds) too_slow <- union(too_slow, key)
      }
    }
    logmsg("  speed: n = ", n, " done")
  }
  # (b) number of covariates at n = 5000, for the regression models
  reg <- Filter(function(m) m$model %in% c("Cox", "Fine-Gray", "cause-specific Cox", "Weibull PH", "direct CIF (joint, Gompertz PH)"), ms)
  for (p in cfg$speed_p) for (rep in seq_len(cfg$speed_reps)) {
    dat <- list(surv = gen_surv(5000, p = p), cr = gen_cr(5000, p = p))
    for (m in reg) {
      tm <- time_it(function() m$fun(dat[[m$family]]))
      rec_add("speed", data.frame(design = "p", model = m$model, impl = m$impl, n = 5000L, p = as.integer(p), rep = rep, seconds = tm$seconds,
                                  status = if (failed(tm$value)) paste("error:", substr(fail_msg(tm$value), 1, 80)) else "ok"))
    }
  }
  logmsg("  speed: covariates done")
  invisible(rec_get("speed"))
}
