# Part 1 (continued). Accuracy of the further performance measures: concordance index (cindex),
# standard errors and paired differences of the Brier score and AUC, observed/expected ratio, smooth
# calibration (ICI, E50, E90), decision curves (decision_curve) and likelihood-ratio tests of nested
# models (anova). Same records as 01_accuracy.R; called at the end of acc_new_methods().

# ICI, E50, E90 written from Austin, Harrell & van Klaveren (2020): Cox (one event type) or Fine-Gray
# (survival::finegray + coxph) model of the outcome on a natural spline of cloglog(risk) with knots at
# the 10th, 50th and 90th percentiles; the calibrated risk at t of every subject
ref_ici <- function(d, risk, t, cr) {
  cl <- log(-log(1 - pmin(pmax(risk, 1e-10), 1 - 1e-10)))
  kn <- stats::quantile(cl, c(0.1, 0.5, 0.9))
  B <- splines::ns(cl, knots = kn[2L], Boundary.knots = kn[c(1L, 3L)])
  colnames(B) <- c("b1", "b2")
  dd <- data.frame(time = d$time, status = d$status, B)
  if (!cr) {
    cf <- survival::coxph(survival::Surv(time, status == 1) ~ b1 + b2, data = dd, ties = "efron")
    bh <- survival::basehaz(cf, centered = FALSE)
  } else {
    dd$ev <- factor(dd$status, 0:2, c("censor", "one", "two"))
    fgd <- survival::finegray(survival::Surv(time, ev) ~ b1 + b2, data = dd, etype = "one")
    cf <- survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ b1 + b2, data = fgd, weights = fgwt, ties = "breslow")
    bh <- survival::basehaz(cf, centered = FALSE)
  }
  H <- max(c(0, bh$hazard[bh$time <= t]))
  obs <- 1 - exp(-H * exp(drop(B %*% stats::coef(cf))))
  a <- abs(obs - risk)
  c(mean(a), stats::median(a), unname(stats::quantile(a, 0.9)))
}

acc_assess_one <- function(d, label, cr) {
  if (!has("survival")) return(invisible())
  fo <- Event(time, status) ~ x1 + x2
  t <- stats::median(d$time[d$status > 0])
  # two models to compare: with both covariates and with x1 only
  f2 <- attempt(if (cr) fg(fo, d) else cox(fo, d))
  f1 <- attempt(if (cr) fg(Event(time, status) ~ x1, d) else cox(Event(time, status) ~ x1, d))
  rk <- function(f) if (cr) drop(predict(f, d, times = t)) else 1 - drop(predict(f, d, times = t, type = "survival"))
  r2 <- if (failed(f2)) f2 else attempt(rk(f2)); r1 <- if (failed(f1)) f1 else attempt(rk(f1))
  if (failed(r1) || failed(r2)) {
    check("prediction measures", "compare measures", "-", label, "models", r2, r1); return(invisible())
  }
  # ---- concordance --------------------------------------------------------------------------------
  sc <- d$x1 - 0.5 * d$x2 + 0.3 * round(d$x2)            # a score with ties
  if (!cr) {
    o <- attempt(cindex(sc, d$time, d$status))
    r <- attempt(survival::concordance(survival::Surv(time, status) ~ sc, data = data.frame(time = d$time, status = d$status, sc = sc),
                                       reverse = TRUE))
    check("concordance", "cindex (Harrell)", "survival::concordance(reverse = TRUE)", label, "C and standard error",
          if (failed(o)) o else c(o$concordance, o$se), if (failed(r)) r else c(r$concordance, sqrt(r$var)))
    tau <- unname(stats::quantile(d$time[d$status == 1], 0.8, type = 1))
    o <- attempt(cindex(sc, d$time, d$status, method = "uno", tau = tau))
    r <- attempt(survival::concordance(survival::Surv(time, status) ~ sc, data = data.frame(time = d$time, status = d$status, sc = sc),
                                       reverse = TRUE, timewt = "n/G2", ymax = tau))
    check("concordance", "cindex (Uno, tau)", "survival::concordance(timewt = \"n/G2\", ymax)", label, "C",
          if (failed(o)) o else o$concordance, if (failed(r)) r else unname(r$concordance))
  } else if (has("pec")) {
    o <- attempt(cindex(sc, d$time, d$status, cause = 1L, method = "uno"))
    r <- attempt({
      dd <- data.frame(time = d$time, status = d$status, sc = sc)
      Hist <- prodlim::Hist                                # pec and riskRegression read the name Hist from the formula
      pc <- pec::cindex(list(m = matrix(sc, ncol = 1L)), formula = Hist(time, status) ~ 1, data = dd, cause = 1,
                        eval.times = max(dd$time) + 1, cens.model = "marginal", splitMethod = "none")
      unname(unlist(pc$AppCindex))[1L]
    })
    check("concordance", "cindex (Uno, competing risks: Wolbers et al. 2014)", "pec::cindex(cause = 1)", label, "C",
          if (failed(o)) o else o$concordance, r)
  }
  # ---- standard errors and paired differences (riskRegression::Score, conservative) ---------------
  ours <- attempt(.sx_assess_ext(d, list(m1 = r1, m2 = r2), t))
  if (has("riskRegression") && has("prodlim")) {
    ref <- attempt({
      suppressPackageStartupMessages(library(survival))
      Hist <- prodlim::Hist                                # Score reads the name Hist from the formula
      scr <- riskRegression::Score(list(m1 = matrix(r1), m2 = matrix(r2)), formula = Hist(time, status) ~ 1, data = d,
                                  times = t, cause = 1, metrics = c("auc", "brier"), null.model = FALSE, conservative = TRUE,
                                  se.fit = TRUE, cens.model = "km", contrasts = TRUE)
      a <- as.data.frame(scr$AUC$score); b <- as.data.frame(scr$Brier$score)
      ac <- as.data.frame(scr$AUC$contrasts); bc <- as.data.frame(scr$Brier$contrasts)
      list(auc_se = a$se[match(c("m1", "m2"), a$model)], brier_se = b$se[match(c("m1", "m2"), b$model)],
           dauc = ac$delta.AUC[1L], dauc_se = ac$se[1L], dauc_p = ac$p[1L], dbr = bc$delta.Brier[1L], dbr_se = bc$se[1L], dbr_p = bc$p[1L])
    })
    if (!failed(ours) && !failed(ref)) {
      check("prediction measures", "compare measures (se)", "riskRegression::Score(conservative = TRUE)", label, "se of the AUC", ours$auc_se, ref$auc_se)
      check("prediction measures", "compare measures (se)", "riskRegression::Score(conservative = TRUE)", label, "se of the Brier score", ours$brier_se, ref$brier_se)
      check("prediction measures", "compare measures (difference)", "riskRegression::Score(conservative = TRUE)", label,
            "m2 - m1: AUC, se, p", c(ours$dauc, ours$dauc_se, ours$dauc_p), c(ref$dauc, ref$dauc_se, ref$dauc_p))
      check("prediction measures", "compare measures (difference)", "riskRegression::Score(conservative = TRUE)", label,
            "m2 - m1: Brier, se, p", c(ours$dbr, ours$dbr_se, ours$dbr_p), c(ref$dbr, ref$dbr_se, ref$dbr_p))
    } else check("prediction measures", "compare measures (se)", "riskRegression::Score", label, "se", ours, ref)
  }
  # ---- O/E and smooth calibration -----------------------------------------------------------------
  cal <- attempt(calibration(d$time, d$status, r2, t))
  O <- attempt({
    s <- survival::survfit(survival::Surv(time, factor(status)) ~ 1, data = d)
    ps <- summary(s, times = t, extend = TRUE)$pstate
    ps[, which(s$states == "1")]
  })
  check("calibration", "calibration (O/E)", "survival::survfit (Aalen-Johansen)", label, "O/E",
        if (failed(cal)) cal else attr(cal, "oe")[1L], if (failed(O)) O else O / mean(r2))
  if (has("splines")) {
    check("calibration", "calibration (ICI, E50, E90)", if (cr) "survival::finegray + coxph + splines::ns" else "survival::coxph + splines::ns",
          label, "ICI, E50, E90", if (failed(cal)) cal else c(attr(cal, "ici"), attr(cal, "e50"), attr(cal, "e90")),
          attempt(ref_ici(d, r2, t, cr)), tol = if (cr) "loose" else "exact",
          note = if (cr) "finegray() weights the expanded rows by its own censoring estimate" else "")
  }
  # ---- decision curve -----------------------------------------------------------------------------
  th <- c(0.05, 0.1, 0.2, 0.3, 0.4)
  dc <- attempt(decision_curve(list(m = r2), d$time, d$status, t = t, thresholds = th))
  ours_nb <- if (failed(dc)) dc else c(dc$table$net_benefit[dc$table$model == "m"], dc$table$net_benefit[dc$table$model == "treat all"])
  km_nb <- attempt({
    risk_of <- function(i) {
      if (!any(i)) return(0)
      s <- survival::survfit(survival::Surv(d$time[i], factor(d$status[i], 0:2)) ~ 1)
      ps <- summary(s, times = t, extend = TRUE)$pstate
      ps[, which(s$states == "1")]
    }
    nb <- vapply(th, function(p) { i <- r2 >= p; pi <- mean(i); F <- risk_of(i); pi * F - pi * (1 - F) * p / (1 - p) }, 0)
    Fa <- risk_of(rep(TRUE, nrow(d)))
    c(nb, Fa - (1 - Fa) * th / (1 - th))
  })
  check("decision curve", "decision_curve", "survival::survfit (definition of Vickers et al. 2008)", label, "net benefit", ours_nb, km_nb)
  if (!cr && has("dcurves")) {
    ref <- attempt({
      dd <- data.frame(time = d$time, status = d$status, m = r2)
      x <- dcurves::dca(survival::Surv(time, status) ~ m, data = dd, time = t, thresholds = th)
      tb <- as.data.frame(x$dca)
      c(tb$net_benefit[tb$variable == "m"], tb$net_benefit[tb$variable == "all"])
    })
    check("decision curve", "decision_curve", "dcurves::dca", label, "net benefit", ours_nb, ref)
  }
  # ---- likelihood-ratio tests ---------------------------------------------------------------------
  if (!cr) {
    a <- attempt(anova(cox(Event(time, status) ~ x1, d), cox(fo, d)))
    r <- attempt(anova(survival::coxph(survival::Surv(time, status) ~ x1, d), survival::coxph(survival::Surv(time, status) ~ x1 + x2, d)))
    check("likelihood-ratio test", "anova (Cox)", "survival::anova.coxph", label, "chi-square and p",
          if (failed(a)) a else c(a$Chisq[2L], a[["Pr(>|Chi|)"]][2L]), if (failed(r)) r else c(r$Chisq[2L], r[["Pr(>|Chi|)"]][2L]))
    a <- attempt(anova(parsurv(Event(time, status) ~ x1, d), parsurv(fo, d)))
    r <- attempt({
      st <- 2 * (survival::survreg(survival::Surv(time, status) ~ x1 + x2, d)$loglik[2L] - survival::survreg(survival::Surv(time, status) ~ x1, d)$loglik[2L])
      c(st, stats::pchisq(st, 1, lower.tail = FALSE))
    })
    check("likelihood-ratio test", "anova (Weibull)", "survival::survreg log-likelihoods", label, "chi-square and p",
          if (failed(a)) a else c(a$Chisq[2L], a[["Pr(>|Chi|)"]][2L]), r, tol = "optimiser")
  } else {
    a <- attempt(anova(csc(Event(time, status) ~ x1, d), csc(fo, d)))
    r <- attempt({
      ll <- function(vars) sum(vapply(1:2, function(k) {
        dk <- d; dk$ev <- as.integer(d$status == k)          # cause k as the event, the others censored
        survival::coxph(stats::reformulate(vars, "survival::Surv(time, ev)"), data = dk)$loglik[2L]
      }, 0))
      st <- 2 * (ll(c("x1", "x2")) - ll("x1"))
      c(st, stats::pchisq(st, 2, lower.tail = FALSE))
    })
    check("likelihood-ratio test", "anova (cause-specific Cox)", "survival::coxph per cause", label, "chi-square and p",
          if (failed(a)) a else c(a$Chisq[2L], a[["Pr(>|Chi|)"]][2L]), r)
  }
  invisible()
}

# the measures of two models from the C++ core: standard errors and the paired difference m2 - m1
.sx_assess_ext <- function(d, risks, t) {
  A <- SurvX:::.sx_assess(d$time, d$status, risks, t, 1L)
  a <- A$a[, 1L, ]
  ref <- A$ref[1L]; other <- 3L - ref
  sgn <- if (other == 2L) 1 else -1                       # the core reports model - reference
  list(auc_se = unname(a["auc_se", ]), brier_se = unname(a["brier_se", ]),
       dauc = sgn * a["auc_diff", other], dauc_se = a["auc_diff_se", other], dauc_p = a["auc_diff_p", other],
       dbr = sgn * a["brier_diff", other], dbr_se = a["brier_diff_se", other], dbr_p = a["brier_diff_p", other])
}

acc_assess <- function(cfg) {
  logmsg("  further performance measures: cindex, standard errors and differences, O/E, ICI, decision curves, anova")
  for (rep in seq_len(cfg$acc_reps)) {
    lab <- paste0("#", rep)
    acc_assess_one(gen_surv(500, p = 2), paste("standard", lab), FALSE)
    acc_assess_one(gen_cr(500, p = 2), paste("cr_standard", lab), TRUE)
  }
  invisible()
}
