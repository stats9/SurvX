# Part 3. Monte Carlo study of the statistical properties: for every model, data are generated
# from the model itself (known true values), the model is fitted R times, and the study reports the
# bias, the empirical standard deviation of the estimates, the mean model-based standard error,
# their ratio, the coverage of the 95% Wald interval, and the share of fits that converged.
# For tests: rejection rate at the 5% level under the null (size) and under an alternative (power).
# Coverage is judged against the Monte Carlo band 95% +/- 1.96 sqrt(0.95 * 0.05 / R).

# ---- extract estimates -------------------------------------------------------------------------------
# rows of estimates(fit) for the given terms (and cause), as data.frame(term, est, se)
grab <- function(fit, terms, cause = NULL, role = NULL) {
  e <- estimates(fit)
  sel <- e$term %in% terms
  if (!is.null(cause)) sel <- sel & !is.na(e$cause) & e$cause == cause
  if (!is.null(role)) sel <- sel & e$role %in% role
  e <- e[sel, , drop = FALSE]
  e[match(terms, e$term), c("term", "estimate", "se")]
}

# one replication result: rows of (scenario, method, parameter, true, est, se, converged, judged)
res_rows <- function(scenario, method, truth, fit, cause = NULL, role = NULL, judged = TRUE) {
  if (failed(fit)) return(data.frame(scenario = scenario, method = method, parameter = names(truth), true = unname(truth),
                                     est = NA_real_, se = NA_real_, converged = FALSE, judged = judged, stringsAsFactors = FALSE))
  g <- grab(fit, names(truth), cause, role)
  conv <- tryCatch(isTRUE(fit@converged), error = function(e) TRUE)
  data.frame(scenario = scenario, method = method, parameter = names(truth), true = unname(truth), est = g$estimate, se = g$se,
             converged = conv, judged = judged, stringsAsFactors = FALSE)
}
val_rows <- function(scenario, method, parameter, true, est, se = NA_real_, judged = TRUE)
  data.frame(scenario = scenario, method = method, parameter = parameter, true = true, est = est, se = se,
             converged = is.finite(est), judged = judged, stringsAsFactors = FALSE)

# ---- generators with known truth ---------------------------------------------------------------------
# Weibull PH, x1 ~ Bernoulli(0.5), x2 ~ N(0, 1); S(t | x) = exp{-(t / 2)^1.3 exp(b'x)}
sim_weibull <- function(n, beta = c(x1 = 0.5, x2 = -0.4), cens_max = 6) {
  x1 <- rbinom(n, 1, 0.5); x2 <- rnorm(n)
  t <- 2 * (-log(runif(n)) / exp(beta[1] * x1 + beta[2] * x2))^(1 / 1.3)
  cens <- runif(n, 0, cens_max)
  data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x1 = x1, x2 = x2, latent = t)
}

# left truncation: subjects are seen only if T > L, L ~ U(0, 1.5); sampled until n are observed
sim_truncated <- function(n, gen) {
  out <- NULL
  while (is.null(out) || nrow(out) < n) {
    d <- gen(2L * n); L <- runif(nrow(d), 0, 1.5)
    keep <- d$latent > L
    d <- d[keep, , drop = FALSE]; d$entry <- L[keep]
    out <- rbind(out, d)
  }
  out <- out[seq_len(n), ]
  out[out$time > out$entry, , drop = FALSE]
}

# time-varying covariate: z(t) = 1 after a switching time s ~ U(0, 2) for half of the subjects;
# constant hazard 0.4 exp(b z(t) + b2 x2)
sim_tv <- function(n, b = 0.7, b2 = -0.3) {
  x2 <- rnorm(n); sw <- ifelse(runif(n) < 0.5, runif(n, 0, 2), Inf)
  h0 <- 0.4 * exp(b2 * x2); E <- -log(runif(n))
  t <- ifelse(E < h0 * sw, E / h0, sw + (E - h0 * sw) / (h0 * exp(b)))
  cens <- runif(n, 0, 5); time <- pmin(t, cens); status <- as.integer(t <= cens)
  rows <- lapply(seq_len(n), function(i) {
    if (time[i] <= sw[i]) data.frame(id = i, start = 0, stop = time[i], event = status[i], z = 0, x2 = x2[i])
    else data.frame(id = i, start = c(0, sw[i]), stop = c(sw[i], time[i]), event = c(0L, status[i]), z = c(0, 1), x2 = x2[i])
  })
  do.call(rbind, rows)
}

# Fine-Gray model as in Fine & Gray (1999, Sec. 5): F1(t | z) = 1 - [1 - p{1 - exp(-t)}]^exp(b1'z),
# cause 2 given cause 2 exponential with rate exp(b2'z); z1 ~ Bernoulli(0.5), z2 ~ N(0, 1)
sim_fg <- function(n, b1 = c(x1 = 0.5, x2 = -0.5), b2 = c(-0.5, 0.5), p = 0.4, cens_max = 4) {
  x1 <- rbinom(n, 1, 0.5); x2 <- rnorm(n)
  e1 <- exp(b1[1] * x1 + b1[2] * x2)
  pinf <- 1 - (1 - p)^e1
  U <- runif(n)
  c1 <- U < pinf
  t <- numeric(n)
  t[c1] <- -log(1 - (1 - (1 - U[c1])^(1 / e1[c1])) / p)
  t[!c1] <- rexp(sum(!c1), exp(b2[1] * x1[!c1] + b2[2] * x2[!c1]))
  cens <- runif(n, 0, cens_max)
  data.frame(time = pmin(t, cens), status = ifelse(t <= cens, ifelse(c1, 1L, 2L), 0L), x1 = x1, x2 = x2, latent = t)
}

# cause-specific exponential hazards
sim_csc <- function(n, b1 = c(x1 = 0.5, x2 = -0.3), b2 = c(x1 = -0.4, x2 = 0.3)) {
  x1 <- rbinom(n, 1, 0.5); x2 <- rnorm(n)
  t1 <- rexp(n, 0.4 * exp(b1[1] * x1 + b1[2] * x2)); t2 <- rexp(n, 0.25 * exp(b2[1] * x1 + b2[2] * x2))
  cens <- runif(n, 0, 4); time <- pmin(t1, t2, cens)
  data.frame(time = time, status = ifelse(time == cens, 0L, ifelse(t1 < t2, 1L, 2L)), x1 = x1, x2 = x2)
}

# Direct CIF models of Jeong & Fine with Gompertz baselines u(t) = tau (exp(rho t) - 1) / rho (rho < 0):
# F_k(t | z) = 1 - {1 + a e^{b'z} u_k(t)}^{-1/a} (a = 0: PH, a = 1: PO). x2 ~ U(-1, 1) keeps
# F_1(inf) + F_2(inf) < 1, the rest never fail. design "joint": both causes follow the model;
# design "complete": every subject fails, cause 2 given cause 2 is Weibull(1.5, 1) (the constrained
# model of Shi, Cheng & Jeong with cond = "weibull").
gomp_u <- function(t, rho, tau) tau * (exp(rho * t) - 1) / rho
gomp_inv <- function(u, rho, tau) log(1 + rho * u / tau) / rho
link_F <- function(u, eta, a) if (a == 0) 1 - exp(-exp(eta) * u) else 1 - (1 + a * exp(eta) * u)^(-1 / a)
link_Finv_u <- function(Fv, eta, a) if (a == 0) -log(1 - Fv) / exp(eta) else ((1 - Fv)^(-a) - 1) / (a * exp(eta))
sim_direct <- function(n, a, design = "joint", b1 = c(x1 = 0.5, x2 = -0.3), b2 = c(x1 = -0.3, x2 = 0.2),
                       g1 = c(rho = -0.6, tau = 0.5), g2 = c(rho = -0.8, tau = 0.25)) {
  x1 <- rbinom(n, 1, 0.5); x2 <- runif(n, -1, 1)
  eta1 <- b1[1] * x1 + b1[2] * x2; eta2 <- b2[1] * x1 + b2[2] * x2
  F1inf <- link_F(g1[["tau"]] / -g1[["rho"]], eta1, a)
  F2inf <- if (design == "joint") link_F(g2[["tau"]] / -g2[["rho"]], eta2, a) else 1 - F1inf
  U <- runif(n); t <- rep(Inf, n); cause <- integer(n)
  c1 <- U < F1inf; c2 <- !c1 & U < F1inf + F2inf
  t[c1] <- gomp_inv(link_Finv_u(U[c1], eta1[c1], a), g1[["rho"]], g1[["tau"]])
  if (design == "joint") {
    V <- (U[c2] - F1inf[c2])                                           # F2(t) = V
    t[c2] <- gomp_inv(link_Finv_u(V, eta2[c2], a), g2[["rho"]], g2[["tau"]])
  } else t[c2] <- rweibull(sum(c2), 1.5, 1)
  cause[c1] <- 1L; cause[c2] <- 2L
  cens <- runif(n, 0, 6)
  data.frame(time = pmin(t, cens), status = ifelse(t <= cens, cause, 0L), x1 = x1, x2 = x2)
}

# generalized gamma AFT (Prentice 1974): log T = mu + b x + sigma W, Q^-2 exp(Q W) ~ Gamma(Q^-2)
sim_gengamma <- function(n, mu = 0.5, b = c(x1 = -0.4, x2 = 0.3), sigma = 0.8, Q = 0.5) {
  x1 <- rbinom(n, 1, 0.5); x2 <- rnorm(n)
  G <- rgamma(n, 1 / Q^2); W <- log(Q^2 * G) / Q
  t <- exp(mu + b[1] * x1 + b[2] * x2 + sigma * W)
  cens <- runif(n, 0, quantile(t, 0.95))
  data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x1 = x1, x2 = x2)
}

# cure models: cured fraction expit(g0 + g1 x1), latency Weibull PH with b on x2 (shape 1.3, scale 1.5)
sim_cure <- function(n, type = "mixture", g = c(-0.3, 0.8), b = -0.4) {
  x1 <- rbinom(n, 1, 0.5); x2 <- rnorm(n)
  pi <- plogis(g[1] + g[2] * x1)
  U <- runif(n); t <- rep(Inf, n)
  if (type == "mixture") {
    unc <- U >= pi                                                      # uncured: Weibull latency
    t[unc] <- 1.5 * (-log(runif(sum(unc))) / exp(b * x2[unc]))^(1 / 1.3)
  } else {
    fail <- U < 1 - pi                                                  # S = pi^{F_u}: F_u(t) = log(1 - U) / log(pi)
    Fu <- log(1 - U[fail]) / log(pi[fail])
    t[fail] <- 1.5 * (-log(1 - Fu) / exp(b * x2[fail]))^(1 / 1.3)
  }
  cens <- runif(n, 0, 10)
  data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x1 = x1, x2 = x2)
}

# ---- one replication of every scenario ---------------------------------------------------------------
sim_one <- function(rep, n, seed) {
  set.seed(seed + 7919L * rep + n)
  out <- list()
  add <- function(x) out[[length(out) + 1L]] <<- x
  b <- c(x1 = 0.5, x2 = -0.4)

  # Cox
  d <- sim_weibull(n)
  add(res_rows("Cox PH", "cox (Efron)", b, attempt(cox(Event(time, status) ~ x1 + x2, d))))
  add(res_rows("Cox PH", "cox (Breslow)", b, attempt(cox(Event(time, status) ~ x1 + x2, d, ties = "breslow"))))
  ds <- d; ds$g <- sample(c("A", "B", "C"), n, TRUE)
  ds$time <- ds$time * c(A = 1, B = 0.5, C = 2)[ds$g]                # different baseline by stratum, same b
  add(res_rows("Cox PH, strata", "cox + strata()", b, attempt(cox(Event(time, status) ~ x1 + x2 + strata(g), ds))))
  dl <- sim_truncated(n, sim_weibull)
  add(res_rows("Cox PH, left truncation", "cox, entry", b, attempt(cox(Event(time, status, entry = entry) ~ x1 + x2, dl))))
  add(res_rows("Cox PH, left truncation", "cox ignoring entry (naive)", b, attempt(cox(Event(time, status) ~ x1 + x2, dl)), judged = FALSE))
  dt <- sim_tv(n)
  add(res_rows("Cox, time-varying covariate", "cox, (start, stop] + cluster(id)", c(z = 0.7, x2 = -0.3),
               attempt(cox(Event(stop, event, entry = start) ~ z + x2 + cluster(id), dt))))
  # robust variance with correlated pairs (gamma frailty, variance 0.5): only the ratio of the robust
  # SE to the empirical SD is judged (the marginal coefficient differs from the conditional one)
  dp <- gen_surv(n, p = 2, beta = c(0.5, -0.4), frailty = 0.5)
  f <- attempt(cox(Event(time, status) ~ x1 + x2 + cluster(id), dp))
  add(transform(res_rows("Cox, clustered pairs", "cox + cluster(id): robust SE", b, f), judged = FALSE))
  if (!failed(f)) add(val_rows("Cox, clustered pairs", "cox: model-based SE", c("x1", "x2"), b, coef(f), sqrt(diag(f@vcov_model)), judged = FALSE))

  # Fine-Gray
  d <- sim_fg(n); b1 <- c(x1 = 0.5, x2 = -0.5)
  add(res_rows("Fine-Gray", "fg (robust SE)", b1, attempt(fg(Event(time, status) ~ x1 + x2, d))))
  add(res_rows("Fine-Gray", "fg (model SE)", b1, attempt(fg(Event(time, status) ~ x1 + x2, d, variance = "model")), judged = FALSE))
  dl <- sim_truncated(n, sim_fg)
  add(res_rows("Fine-Gray, left truncation", "fg, entry (Geskus)", b1, attempt(fg(Event(time, status, entry = entry) ~ x1 + x2, dl))))
  # direct PH model (Gompertz) on Fine-Gray data: misspecified baseline, reported for comparison
  add(res_rows("Fine-Gray", "directcif weighted, Gompertz PH (misspecified baseline)", b1,
               attempt(directcif(Event(time, status) ~ x1 + x2, d, estimation = "weighted")), cause = 1, role = "beta", judged = FALSE))

  # cause-specific Cox
  d <- sim_csc(n)
  f <- attempt(csc(Event(time, status) ~ x1 + x2, d))
  add(res_rows("cause-specific Cox", "csc, cause 1", c(x1 = 0.5, x2 = -0.3), f, cause = 1))
  add(res_rows("cause-specific Cox", "csc, cause 2", c(x1 = -0.4, x2 = 0.3), f, cause = 2))
  f <- attempt(csc(Event(time, status) ~ x1 + x2, d, baseline = "exponential"))
  add(res_rows("cause-specific Cox", "csc exponential, cause 1", c(x1 = 0.5, x2 = -0.3), f, cause = 1))

  # direct CIF models
  for (a in c(0, 1)) {
    lk <- if (a == 0) "PH" else "PO"
    d <- sim_direct(n, a, "joint")
    t1 <- c(x1 = 0.5, x2 = -0.3)
    add(res_rows(paste("direct CIF", lk), "directcif joint", t1, attempt(directcif(Event(time, status) ~ x1 + x2, d, link = lk)), cause = 1, role = "beta"))
    add(res_rows(paste("direct CIF", lk), "directcif joint, cause 2", c(x1 = -0.3, x2 = 0.2),
                 attempt(directcif(Event(time, status) ~ x1 + x2, d, link = lk)), cause = 2, role = "beta"))
    add(res_rows(paste("direct CIF", lk), "directcif weighted", t1,
                 attempt(directcif(Event(time, status) ~ x1 + x2, d, link = lk, estimation = "weighted")), cause = 1, role = "beta"))
    d <- sim_direct(n, a, "complete")
    add(res_rows(paste("direct CIF", lk, "(all fail)"), "directcif constrained, cond = weibull", t1,
                 attempt(directcif(Event(time, status) ~ x1 + x2, d, link = lk, estimation = "constrained", cond = "weibull")), cause = 1, role = "beta"))
    add(res_rows(paste("direct CIF", lk, "(all fail)"), "directcif weighted", t1,
                 attempt(directcif(Event(time, status) ~ x1 + x2, d, link = lk, estimation = "weighted")), cause = 1, role = "beta"))
  }

  # parametric models
  d <- sim_weibull(n)
  add(res_rows("Weibull PH", "parsurv weibull", c(x1 = 0.5, x2 = -0.4, "log(shape)" = log(1.3)),
               attempt(parsurv(Event(time, status) ~ x1 + x2, d))))
  add(res_rows("Weibull PH", "parsurv rp (df 3), coefficients only", b,
               attempt(parsurv(Event(time, status) ~ x1 + x2, d, baseline = "rp", df = 3))))
  pf <- attempt(parsurv(Event(time, status) ~ x1 + x2, d))
  if (!failed(pf)) {
    s1 <- attempt(predict(pf, data.frame(x1 = 1, x2 = 0), times = 1.5, se = TRUE))
    if (!failed(s1)) add(val_rows("Weibull PH", "predict(parsurv): S(1.5 | x1 = 1)", "S(1.5)", exp(-(1.5 / 2)^1.3 * exp(0.5)), s1$est[1, 1], s1$se[1, 1]))
  }
  dl <- sim_truncated(n, sim_weibull)
  add(res_rows("Weibull PH, left truncation", "parsurv weibull, entry", c(x1 = 0.5, x2 = -0.4, "log(shape)" = log(1.3)),
               attempt(parsurv(Event(time, status, entry = entry) ~ x1 + x2, dl))))
  d <- sim_gengamma(n)
  add(res_rows("generalized gamma", "parsurv gengamma", c(x1 = -0.4, x2 = 0.3, "log(scale)" = log(0.8), Q = 0.5),
               attempt(parsurv(Event(time, status) ~ x1 + x2, d, baseline = "gengamma"))))

  # cure models
  for (ty in c("mixture", "nonmixture")) {
    d <- sim_cure(n, ty)
    add(res_rows(paste("cure,", ty), paste("cure", ty), c("cure: (Intercept)" = -0.3, "cure: x1" = 0.8, x2 = -0.4, "log(shape)" = log(1.3)),
                 attempt(cure(Event(time, status) ~ x2, d, cure = ~ x1, type = ty))))
  }

  # non-parametric estimators: coverage of the pointwise intervals (one group, Weibull(1.3, 2))
  t0 <- 2 * (-log(runif(n)))^(1 / 1.3); cc <- runif(n, 0, 6)
  dn <- data.frame(time = pmin(t0, cc), status = as.integer(t0 <= cc))
  k <- attempt(summary(km(Event(time, status) ~ 1, dn), times = c(1, 2, 3))@tables[[1]])
  if (!failed(k)) for (j in 1:3) {
    tt <- c(1, 2, 3)[j]; S <- exp(-(tt / 2)^1.3)
    add(data.frame(scenario = "Kaplan-Meier", method = "km, log-scale interval", parameter = paste0("S(", tt, ")"), true = S, est = k$surv[j],
                   se = k$se[j], converged = TRUE, judged = TRUE, lower = k$lower[j], upper = k$upper[j]))
  }
  rm4 <- attempt(rmst(dn$time, dn$status, 4))
  true_rmst <- stats::integrate(function(t) exp(-(t / 2)^1.3), 0, 4)$value
  if (!failed(rm4)) add(val_rows("RMST", "rmst (tau = 4)", "RMST(4)", true_rmst, rm4[["rmst"]], rm4[["se"]]))
  t1 <- rexp(n, 0.3); t2 <- rexp(n, 0.2); cc <- runif(n, 0, 5); tm <- pmin(t1, t2, cc)
  da <- data.frame(time = tm, status = ifelse(tm == cc, 0L, ifelse(t1 < t2, 1L, 2L)))
  a <- attempt(summary(aj(Event(time, status) ~ 1, da), times = c(1, 2, 3))@tables[[1]])
  if (!failed(a)) for (j in 1:3) {
    tt <- c(1, 2, 3)[j]; F1 <- 0.3 / 0.5 * (1 - exp(-0.5 * tt))
    add(data.frame(scenario = "Aalen-Johansen", method = "aj, cloglog interval", parameter = paste0("F1(", tt, ")"), true = F1, est = a$cif.1[j],
                   se = a$se.1[j], converged = TRUE, judged = TRUE, lower = a$lower.1[j], upper = a$upper.1[j]))
  }

  # tests: rejection at 5% (true = 0.05 under the null; power reported)
  g <- rbinom(n, 1, 0.5)
  h0 <- rexp(n, 0.3); h1 <- rexp(n, 0.3 * exp(0.4 * g)); cc <- runif(n, 0, 5)
  rej <- function(p) as.numeric(p < 0.05)
  lr0 <- attempt(logrank_test(pmin(h0, cc), as.integer(h0 <= cc), g)$p.value)
  lr1 <- attempt(logrank_test(pmin(h1, cc), as.integer(h1 <= cc), g)$p.value)
  if (!failed(lr0)) add(val_rows("tests", "logrank_test: size", "rejection rate", 0.05, rej(lr0)))
  if (!failed(lr1)) add(val_rows("tests", "logrank_test: power (HR 1.5)", "rejection rate", NA_real_, rej(lr1), judged = FALSE))
  dc <- sim_csc(n); dc$g <- rbinom(n, 1, 0.5)
  gt <- attempt(gray_test(dc$time, dc$status, dc$g)$p.value)
  if (!failed(gt)) add(val_rows("tests", "gray_test: size", "rejection rate", 0.05, rej(gt)))
  dw <- sim_weibull(n)
  pt <- attempt(ph_test(cox(Event(time, status) ~ x1 + x2, dw))@table)
  if (!failed(pt)) add(val_rows("tests", "ph_test (km): size, global", "rejection rate", 0.05, rej(pt$p[pt$term == "GLOBAL"])))
  x <- rbinom(n, 1, 0.5); tn <- rweibull(n, ifelse(x == 1, 0.7, 1.6), 2); cc <- runif(n, 0, 5)
  dnph <- data.frame(time = pmin(tn, cc), status = as.integer(tn <= cc), x = x)
  pt <- attempt(ph_test(cox(Event(time, status) ~ x, dnph))@table)
  if (!failed(pt)) add(val_rows("tests", "ph_test (km): power, crossing hazards", "rejection rate", NA_real_, rej(pt$p[pt$term == "x"]), judged = FALSE))

  r <- do.call(rbind, lapply(out, function(x) { for (c in c("lower", "upper")) if (is.null(x[[c]])) x[[c]] <- NA_real_; x }))
  r$rep <- rep; r$n <- n
  r
}

run_simulation <- function(cfg) {
  logmsg("Part 3: Monte Carlo study (", cfg$sim_R, " replications per scenario, n = ", paste(cfg$sim_n, collapse = ", "),
         ", ", cfg$cores, " cores)")
  all <- list()
  for (n in cfg$sim_n) {
    t0 <- proc.time()[["elapsed"]]
    res <- plapply(seq_len(cfg$sim_R), function(r) tryCatch(sim_one(r, n, cfg$seed + 3L),
                                                            error = function(e) data.frame(scenario = "ERROR", method = conditionMessage(e), rep = r, n = n)),
                   cfg$cores)
    cols <- unique(unlist(lapply(res, names)))
    res <- do.call(rbind, lapply(res, function(x) { for (c in setdiff(cols, names(x))) x[[c]] <- NA; x[cols] }))
    all[[length(all) + 1L]] <- res
    logmsg("  simulation n = ", n, " done in ", round(proc.time()[["elapsed"]] - t0), " s")
  }
  raw <- do.call(rbind, all)
  assign("simulation", list(raw), envir = .REC)
  invisible(raw)
}

# summary table of the Monte Carlo study
summarise_simulation <- function(raw) {
  raw <- raw[raw$scenario != "ERROR", , drop = FALSE]
  key <- interaction(raw$scenario, raw$method, raw$parameter, raw$n, drop = TRUE, lex.order = TRUE)
  out <- lapply(split(raw, key), function(s) {
    R <- nrow(s); ok <- is.finite(s$est)
    est <- s$est[ok]; se <- s$se[ok]; tr <- s$true[1]
    is_rate <- s$parameter[1] == "rejection rate"
    lo <- ifelse(is.finite(s$lower[ok]), s$lower[ok], est - 1.959964 * se)
    hi <- ifelse(is.finite(s$upper[ok]), s$upper[ok], est + 1.959964 * se)
    cover <- if (!is_rate && is.finite(tr) && any(is.finite(se))) mean(lo <= tr & tr <= hi, na.rm = TRUE) else NA_real_
    emp_sd <- if (length(est) > 1L) stats::sd(est) else NA_real_
    data.frame(scenario = s$scenario[1], method = s$method[1], parameter = s$parameter[1], n = s$n[1], R = R,
               converged_pct = 100 * mean(s$converged & ok), true = tr, mean_est = mean(est),
               bias = if (is_rate) NA_real_ else mean(est) - tr,
               rel_bias_pct = if (is_rate || !is.finite(tr) || tr == 0) NA_real_ else 100 * (mean(est) - tr) / abs(tr),
               emp_sd = if (is_rate) NA_real_ else emp_sd, mean_se = if (is_rate) NA_real_ else mean(se, na.rm = TRUE),
               se_ratio = if (is_rate) NA_real_ else mean(se, na.rm = TRUE) / emp_sd,
               coverage_pct = 100 * cover, rate_pct = if (is_rate) 100 * mean(est) else NA_real_,
               mc_band = 100 * 1.959964 * sqrt(0.95 * 0.05 / sum(ok)), judged = s$judged[1], stringsAsFactors = FALSE)
  })
  sm <- do.call(rbind, out)
  band <- sm$mc_band
  sm$verdict <- ifelse(!sm$judged, "reported",
                ifelse(!is.na(sm$rate_pct), ifelse(abs(sm$rate_pct - 5) <= 100 * 1.959964 * sqrt(0.05 * 0.95 / sm$R), "OK", "CHECK"),
                ifelse(is.na(sm$coverage_pct), "reported",
                ifelse(abs(sm$coverage_pct - 95) <= band & sm$converged_pct >= 99 & abs(sm$se_ratio - 1) <= 0.1, "OK", "CHECK"))))
  sm[order(sm$scenario, sm$method, sm$parameter, sm$n), ]
}
