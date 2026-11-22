# Part 1 (continued). Accuracy of the methods added after 1.0.0: time-varying effects
# (time_effects), adjusted risks (adjusted_risk), interval censoring (Interval, turnbull, parsurv,
# directcif), shared gamma frailty (cox and parsurv with frailty()) and covariates of the time
# distributions of the constrained direct model (directcif(timing = )). Same records as
# 01_accuracy.R; called at the end of run_accuracy().

# ---- data generators ---------------------------------------------------------------------------------
# interval censoring: Weibull PH times seen at irregular visits (gaps 0.3 to 1.2), up to 10 visits
gen_ic <- function(n, p = 2L) {
  d <- gen_surv(n, p = p, cens_max = 1e8)              # no censoring of the latent times
  lat <- d$time
  L <- R <- numeric(n)
  for (i in seq_len(n)) {
    v <- cumsum(runif(10, 0.3, 1.2))
    L[i] <- max(c(0, v[v < lat[i]]))
    R[i] <- if (any(v >= lat[i])) min(v[v >= lat[i]]) else Inf
  }
  d$left <- L; d$right <- R
  d$left2 <- ifelse(L > 0, L, NA); d$right2 <- ifelse(is.finite(R), R, NA)
  d
}

# shared gamma frailty: C clusters of size m, frailty variance theta, Weibull PH baseline
gen_frailty <- function(C, m, theta = 0.5) {
  id <- rep(seq_len(C), each = m); n <- C * m
  z <- rgamma(C, 1 / theta, 1 / theta)[id]
  x1 <- rbinom(n, 1, 0.5); x2 <- rnorm(n)
  t <- 2 * (-log(runif(n)) / (z * exp(0.5 * x1 - 0.4 * x2)))^(1 / 1.3)
  cens <- runif(n, 0, 6)
  data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x1 = x1, x2 = x2, id = id)
}

# constrained model with time ratios of the competing causes: cause 1 Gompertz PH (rho -0.8, tau 0.6,
# effects 0.5 x1 - 0.3 x2); the others share the rest (logit of cause 3: -0.2 + 0.4 x1) with Weibull
# times and time ratios exp(0.5 x2) (cause 2, shape 1.6, rate 1) and exp(-0.4 x2) (cause 3, 1.3, 0.6)
gen_aft <- function(n) {
  x1 <- rbinom(n, 1, 0.5); x2 <- rnorm(n)
  lp <- 0.5 * x1 - 0.3 * x2
  p1 <- 1 - exp(-exp(lp) * 0.6 / 0.8)
  t <- numeric(n); s <- integer(n)
  one <- runif(n) < p1
  H <- -log(1 - runif(sum(one)) * p1[one]) / exp(lp[one])
  t[one] <- log1p(H * (-0.8) / 0.6) / (-0.8); s[one] <- 1L
  three <- !one & runif(n) < stats::plogis(-0.2 + 0.4 * x1); two <- !one & !three
  t[two] <- (-log(runif(sum(two))))^(1 / 1.6) * exp(0.5 * x2[two]); s[two] <- 2L
  t[three] <- (-log(runif(sum(three))))^(1 / 1.3) / 0.6 * exp(-0.4 * x2[three]); s[three] <- 3L
  cens <- runif(n, 0.5, 6)
  data.frame(time = pmin(t, cens), status = ifelse(t <= cens, s, 0L), x1 = x1, x2 = x2)
}
# true F_2(t | x1, x2) of gen_aft
aft_true_F2 <- function(t, x1, x2) {
  q <- exp(-exp(0.5 * x1 - 0.3 * x2) * 0.75)
  q * (1 - stats::plogis(-0.2 + 0.4 * x1)) * (1 - exp(-(t * exp(-0.5 * x2))^1.6))
}

# independent log-likelihood of the constrained model with three causes, Gompertz PH for cause 1,
# Weibull G_k with time ratios on the covariates tvars (written from the model definition)
aft_loglik <- function(b, d, vars, tvars) {
  X <- as.matrix(d[vars])
  lp <- drop(X %*% b[paste0("1:", vars)])
  rho <- b[["1:rho"]]; tau <- exp(b[["1:log_tau"]])
  u <- tau * expm1(rho * d$time) / rho
  F1 <- 1 - exp(-exp(lp) * u)
  lf1 <- lp + log(tau) + rho * d$time - exp(lp) * u
  q <- exp(-exp(lp) * tau / (-rho))
  l3 <- b[["3:share_(Intercept)"]] + drop(X %*% b[paste0("3:share_", vars)])
  pi <- cbind(1 - stats::plogis(l3), stats::plogis(l3))
  G <- lg <- matrix(0, nrow(d), 2)
  for (j in 1:2) {
    k <- j + 1L
    u_k <- drop(as.matrix(d[tvars]) %*% b[paste0(k, ":time_", tvars)])
    lam <- exp(b[[paste0(k, ":log_lambda")]] - u_k); kap <- exp(b[[paste0(k, ":log_kappa")]])
    y <- (lam * d$time)^kap
    G[, j] <- -expm1(-y); lg[, j] <- log(kap) + log(y) - y - log(d$time)
  }
  s <- d$status
  sum(lf1[s == 1L]) + sum(log(q[s == 2L] * pi[s == 2L, 1]) + lg[s == 2L, 1]) +
    sum(log(q[s == 3L] * pi[s == 3L, 2]) + lg[s == 3L, 2]) +
    sum(log(1 - F1[s == 0L] - q[s == 0L] * rowSums(pi[s == 0L, , drop = FALSE] * G[s == 0L, , drop = FALSE])))
}

col_of <- function(df, names) { k <- intersect(names, colnames(df))[1]; if (is.na(k)) stop("no column ", paste(names, collapse = "/")); df[[k]] }

# ---- time_effects ------------------------------------------------------------------------------------
acc_time_effects <- function() {
  if (!has("survival")) return(invisible())
  vet <- survival::veteran
  fv <- attempt(cox(Event(time, status) ~ trt + celltype + karno, vet))
  x <- if (failed(fv)) fv else attempt(time_effects(fv, forms = c("constant", "log"),
                                                    fix = list(trt = "constant", celltype = "constant", karno = "log")))
  et <- sort(unique(vet$time[vet$status == 1])); tm <- stats::median(et)
  r <- attempt(survival::coxph(survival::Surv(time, status) ~ trt + celltype + karno + tt(karno), data = vet,
                               tt = function(x, t, ...) x * log(t / tm)))
  check("time-varying effects", "time_effects (log form)", "survival::coxph(tt = )", "veteran", "coefficients",
        if (failed(x)) x else attempt(unname(x@units[[1]]$final$coefficients)), if (failed(r)) r else attempt(unname(coef(r))))
  invisible()
}

# ---- constrained model with covariates of the time distributions G_k (timing) ----------------------
acc_constrained_aft <- function(d, label) {
  fo <- Event(time, status) ~ x1 + x2
  f0 <- attempt(directcif(fo, d, estimation = "constrained", cause = 1, cond = "weibull"))
  f <- attempt(directcif(fo, d, estimation = "constrained", cause = 1, cond = "weibull", timing = ~ x2))
  check("constrained CIF, AFT timing", "directcif(timing = ~ x2)", "independent R likelihood", label, "log-likelihood",
        if (failed(f)) f else attempt(as.numeric(logLik(f))),
        if (failed(f)) f else attempt(aft_loglik(coef(f), d, c("x1", "x2"), "x2")))
  # the maximum: a general optimiser of the independent likelihood, started away from it
  opt <- if (failed(f)) f else attempt({
    b0 <- coef(f)
    nll <- function(v) { b <- stats::setNames(v, names(b0)); r <- -aft_loglik(b, d, c("x1", "x2"), "x2"); if (is.finite(r)) r else 1e10 }
    o <- stats::optim(b0 + 0.05, nll, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14))
    o$par
  })
  check("constrained CIF, AFT timing", "directcif(timing = ~ x2)", "optim (BFGS) of the independent likelihood", label,
        "coefficients", if (failed(f)) f else attempt(unname(coef(f))), if (failed(opt)) opt else attempt(unname(opt)), tol = "optimiser")
  check("constrained CIF, AFT timing", "directcif(timing = ~ x2)", "model without timing (nested)", label,
        "log-likelihood not below the nested model",
        if (failed(f)) f else if (failed(f0)) f0 else attempt(as.numeric(as.numeric(logLik(f)) >= as.numeric(logLik(f0)) - 1e-8)), 1)
  invisible()
}

# ---- adjusted_risk ---# ---- adjusted_risk -----------------------------------------------------------------------------------
acc_adjusted_risk <- function(d, label) {
  if (!has("survival")) return(invisible())
  d$trt <- factor(d$x1, 0:1, c("control", "treated"))
  tq <- event_quantiles(d)
  o <- attempt(cox(Event(time, status) ~ trt + x2, d, ties = "breslow"))
  ar <- if (failed(o)) o else attempt(adjusted_risk(o, d, "trt", times = tq, method = "gformula", B = 0))
  ours <- if (failed(ar)) ar else attempt(ar@risk$risk[order(ar@risk$level, ar@risk$time)])
  if (has("riskRegression")) {
    rc <- attempt(survival::coxph(survival::Surv(time, status) ~ trt + x2, data = d, ties = "breslow", x = TRUE))
    a <- if (failed(rc)) rc else attempt(riskRegression::ate(rc, treatment = "trt", times = tq, data = d, se = FALSE, verbose = FALSE))
    ref <- if (failed(a)) a else attempt({
      mr <- as.data.frame(a$meanRisk)
      v <- col_of(mr, c("estimate", "meanRisk")); lv <- as.character(col_of(mr, c("treatment", "level"))); tt <- col_of(mr, "time")
      v[order(lv, tt)]
    })
    check("adjusted risk", "adjusted_risk (g-formula, Cox)", "riskRegression::ate", label, "risk at the event-time quartiles", ours, ref)
  }
  # inverse probability weighting: weighted Kaplan-Meier by group with logistic propensity weights
  ai <- if (failed(o)) o else attempt(adjusted_risk(o, d, "trt", times = tq, method = "ipw", ps_formula = trt ~ x2, B = 0))
  oi <- if (failed(ai)) ai else attempt(ai@risk$risk[order(ai@risk$level, ai@risk$time)])
  ref <- attempt({
    ps <- stats::fitted(stats::glm(trt ~ x2, data = d, family = stats::binomial()))
    w <- ifelse(d$trt == "treated", mean(d$trt == "treated") / ps, mean(d$trt == "control") / (1 - ps))
    unlist(lapply(c("control", "treated"), function(l) {
      k <- d$trt == l
      1 - summary(survival::survfit(survival::Surv(time, status) ~ 1, data = d[k, ], weights = w[k]), times = tq)$surv
    }))
  })
  check("adjusted risk", "adjusted_risk (IPW)", "survival::survfit (weighted)", label, "risk at the event-time quartiles", oi, ref)
  invisible()
}

# ---- interval censoring ------------------------------------------------------------------------------
acc_interval <- function(d, label, vars = c("x1", "x2")) {
  fo <- function(lhs) fml(lhs, vars)
  # Turnbull NPMLE
  tb <- attempt(turnbull(Interval(left, right) ~ 1, data = d, B = 0))
  if (has("icenReg")) {
    r <- attempt(icenReg::ic_np(cbind(d$left, d$right)))
    check("interval censoring", "turnbull", "icenReg::ic_np", label, "NPMLE log-likelihood",
          if (failed(tb)) tb else attempt(tb@fits[[1]]$loglik), if (failed(r)) r else attempt(r$llk), tol = "optimiser")
  }
  if (has("interval")) {
    # interval's default stopping rule (epsilon 1e-6) stops the EM before the NPMLE; with its
    # defaults the curves differed by up to 2e-3 while the log-likelihood of SurvX matched icenReg
    # to 1e-11. The reference is therefore run to the NPMLE with a tight tolerance.
    ctl <- interval::icfitControl(epsilon = 1e-12, maxit = 1000000)
    ic <- attempt(interval::icfit(survival::Surv(left2, right2, type = "interval2") ~ 1, data = d, control = ctl))
    tt <- if (failed(ic)) NULL else ic$intmap[2, is.finite(ic$intmap[2, ])]
    check("interval censoring", "turnbull", "interval::icfit", label, "survival at the right ends of the innermost intervals",
          if (failed(tb) || failed(ic)) tb else attempt(.subset2(SurvX:::.tb_surv(tb@fits[[1]], tt), "estimate")),
          if (failed(ic)) ic else attempt((1 - cumsum(ic$pf))[is.finite(ic$intmap[2, ])]), tol = "optimiser")
    d$grp <- factor(d$x1)
    tg <- attempt(turnbull(Interval(left, right) ~ grp, data = d, B = 0))
    it <- attempt(interval::ictest(survival::Surv(left2, right2, type = "interval2") ~ grp, data = d, scores = "logrank2", method = "pclt",
                                   icontrol = ctl))
    check("interval censoring", "turnbull (test)", "interval::ictest (logrank2, pclt)", label, "chi-square statistic",
          if (failed(tg)) tg else attempt(tg@test$statistic), if (failed(it)) it else attempt(unname(it$statistic)^2), tol = "optimiser")
  }
  # parametric regression
  if (has("survival")) for (ds in c("weibull", "lognormal", "loglogistic", "exponential")) {
    o <- attempt(parsurv(fo("Interval(left, right)"), d, baseline = ds, scale = "AFT"))
    r <- attempt(survival::survreg(fo("survival::Surv(left2, right2, type = 'interval2')"), data = d, dist = ds))
    k <- length(vars) + 1L
    check("interval censoring", "parsurv(Interval)", "survival::survreg(interval2)", label, paste(ds, "AFT coefficients"),
          if (failed(o)) o else attempt(unname(coef(o)[seq_len(k)])), if (failed(r)) r else attempt(unname(coef(r))), tol = "optimiser")
    check("interval censoring", "parsurv(Interval)", "survival::survreg(interval2)", label, paste(ds, "log-likelihood"),
          if (failed(o)) o else attempt(as.numeric(logLik(o))), if (failed(r)) r else attempt(as.numeric(logLik(r))), tol = "optimiser")
    check("interval censoring", "parsurv(Interval)", "survival::survreg(interval2)", label, paste(ds, "standard errors"),
          if (failed(o)) o else attempt(unname(sqrt(diag(vcov(o)))[seq_len(k)])), if (failed(r)) r else attempt(unname(sqrt(diag(stats::vcov(r)))[seq_len(k)])),
          tol = "optimiser")
  }
  if (has("flexsurv")) for (ds in c("gompertz", "gengamma")) {
    o <- attempt(parsurv(fo("Interval(left, right)"), d, baseline = ds))
    r <- attempt(flexsurv::flexsurvreg(fo("survival::Surv(left2, right2, type = 'interval2')"), data = d, dist = ds))
    check("interval censoring", "parsurv(Interval)", "flexsurv::flexsurvreg(interval2)", label, paste(ds, "log-likelihood"),
          if (failed(o)) o else attempt(as.numeric(logLik(o))), if (failed(r)) r else attempt(r$loglik), tol = "optimiser")
  }
  invisible()
}

# exact and right-censored times written as Interval(): the fits must not change
acc_interval_exact <- function(d, label) {
  d$L <- d$time; d$R <- ifelse(d$status > 0, d$time, Inf)
  a <- attempt(parsurv(Event(time, status) ~ x1 + x2, transform(d, status = as.integer(status > 0))))
  b <- attempt(parsurv(Interval(L, R) ~ x1 + x2, d))
  check("interval censoring", "parsurv(Interval)", "parsurv(Event) on the same data", label, "coefficients (exact times as intervals)",
        if (failed(b)) b else attempt(unname(coef(b))), if (failed(a)) a else attempt(unname(coef(a))), tol = "optimiser")
  if (max(d$status) > 1) {
    a <- attempt(directcif(Event(time, status) ~ x1 + x2, d))
    b <- attempt(directcif(Interval(L, R, cause = status) ~ x1 + x2, d))
    check("interval censoring", "directcif(Interval)", "directcif(Event) on the same data", label, "coefficients (exact times as intervals)",
          if (failed(b)) b else attempt(unname(coef(b))), if (failed(a)) a else attempt(unname(coef(a))), tol = "optimiser")
  }
  invisible()
}

# ---- shared gamma frailty ----------------------------------------------------------------------------
acc_frailty <- function(d, label, ties = FALSE) {
  tem <- if (ties) "info" else "optimiser"
  nem <- if (ties) "tied times: the reference may break ties differently; reported only" else ""
  o <- attempt(cox(Event(time, status) ~ x1 + x2 + frailty(id), d))
  if (has("frailtyEM")) {
    em <- attempt(frailtyEM::emfrail(survival::Surv(time, status) ~ x1 + x2 + cluster(id), data = d))
    check("frailty", "cox + frailty()", "frailtyEM::emfrail", label, "coefficients",
          if (failed(o)) o else attempt(unname(coef(o))), if (failed(em)) em else attempt(unname(em$coefficients)), tol = tem, note = nem)
    check("frailty", "cox + frailty()", "frailtyEM::emfrail", label, "frailty variance",
          if (failed(o)) o else attempt(o@theta), if (failed(em)) em else attempt(1 / exp(em$logtheta)), tol = tem, note = nem)
    check("frailty", "cox + frailty()", "frailtyEM::emfrail", label, "likelihood ratio statistic of no frailty",
          if (failed(o)) o else attempt(o@lrt$statistic), if (failed(em)) em else attempt(2 * diff(em$loglik)), tol = tem, note = nem)
  }
  if (has("survival")) {
    r <- attempt(survival::coxph(survival::Surv(time, status) ~ x1 + x2 + survival::frailty(id, distribution = "gamma"),
                                 data = d, ties = "breslow"))
    check("frailty", "cox + frailty()", "survival::coxph(frailty(), penalized)", label, "coefficients",
          if (failed(o)) o else attempt(unname(coef(o))), if (failed(r)) r else attempt(unname(coef(r))[1:2]), tol = "info",
          note = "coxph maximises a penalized likelihood: close, not identical; reported only")
  }
  if (has("parfm")) for (ds in c("weibull", "gompertz", "exponential", "lognormal", "loglogistic")) {
    p <- attempt(parsurv(Event(time, status) ~ x1 + x2 + frailty(id), d, baseline = ds))
    r <- attempt(parfm::parfm(survival::Surv(time, status) ~ x1 + x2, cluster = "id", data = d, dist = ds, frailty = "gamma"))
    judged <- ds %in% c("weibull", "gompertz", "exponential")
    gain <- if (failed(p) || failed(r)) NA_real_ else tryCatch(p@loglik - attr(r, "loglik"), error = function(e) NA_real_)
    tl <- if (judged && !(is.finite(gain) && gain > 1e-3)) "optimiser" else "info"
    nt <- if (is.finite(gain) && gain > 1e-3) sprintf("SurvX log-likelihood higher by %.2g: the reference did not reach the maximum", gain) else ""
    check("frailty", "parsurv + frailty()", "parfm::parfm", label, paste(ds, "log-likelihood"),
          if (failed(p)) p else attempt(p@loglik), if (failed(r)) r else attempt(attr(r, "loglik")), tol = tl, note = nt)
    check("frailty", "parsurv + frailty()", "parfm::parfm", label, paste(ds, "frailty variance"),
          if (failed(p)) p else attempt(p@theta), if (failed(r)) r else attempt(unname(r["theta", 1])), tol = tl, note = nt)
  }
  invisible()
}

# ---- driver ------------------------------------------------------------------------------------------
acc_new_methods <- function(cfg) {
  logmsg("  methods added after 1.0.0: time_effects, adjusted_risk, interval censoring, frailty, AFT timing of the constrained model")
  progress_start(cfg$acc_reps + 2L)
  acc_time_effects()
  progress_tick()
  for (rep in seq_len(cfg$acc_reps)) {
    lab <- paste0("#", rep)
    acc_adjusted_risk(gen_surv(500, p = 2), paste("standard", lab))
    acc_interval(gen_ic(400), paste("interval", lab))
    acc_interval_exact(gen_surv(400, p = 2), paste("standard", lab))
    acc_interval_exact(gen_cr(400, p = 2), paste("cr_standard", lab))
    acc_frailty(gen_frailty(100, 4), paste("frailty 100 x 4", lab))
    acc_constrained_aft(gen_aft(800), paste("three causes", lab))
    progress_tick()
  }
  if (has("survival")) {
    acc_frailty(transform(survival::kidney, x1 = age, x2 = sex), "kidney", ties = TRUE)
    if (has("interval")) {
      e <- new.env(); utils::data("bcos", package = "interval", envir = e)
      b <- e$bcos
      b <- transform(b, left2 = ifelse(left > 0, left, NA), right2 = ifelse(is.finite(right), right, NA),
                     x1 = as.integer(treatment != levels(treatment)[1]))
      acc_interval(b, "bcos", vars = "x1")
    }
  }
  invisible()
}
