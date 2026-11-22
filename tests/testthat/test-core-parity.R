# The C++ core gives the numbers of the R code it replaced (helper-reference.R).

test_that("distribution functions and sample quantiles agree with R", {
  set.seed(1)
  p <- c(runif(500), 1e-12, 0.025, 0.5, 0.975, 1 - 1e-12)
  expect_equal(SurvX:::.sx_qnorm(p), stats::qnorm(p), tolerance = 1e-13)
  x <- rexp(500, 0.2); df <- sample(1:20, 500, TRUE)
  expect_equal(SurvX:::.sx_pchisq_upper(x, df), stats::pchisq(x, df, lower.tail = FALSE), tolerance = 1e-12)
  sh <- exp(rnorm(500, 1, 1.5)); rt <- exp(rnorm(500))
  pp <- runif(500, 0.001, 0.999)
  expect_equal(SurvX:::.sx_qgamma(pp, sh, rt), stats::qgamma(pp, sh, rt), tolerance = 1e-11)
  for (i in 1:200) {
    v <- if (i %% 2) round(rnorm(sample(1:30, 1)), 1) else rnorm(sample(1:30, 1))
    pr <- c(runif(5), 0, 0.25, 0.5, 0.75, 1)
    expect_equal(SurvX:::.sx_quantile(v, pr, 7L), unname(stats::quantile(v, pr, type = 7)), tolerance = 1e-15)
    expect_equal(SurvX:::.sx_quantile(v, pr, 1L), unname(stats::quantile(v, pr, type = 1)), tolerance = 0)
  }
})

test_that("intervals, step functions and counts agree with R", {
  set.seed(2)
  Fv <- matrix(c(runif(40), 0, 1, NA, 0.5), 11, 4); se <- matrix(runif(44, 0, 0.1), 11, 4)
  expect_equal(SurvX:::.cif_interval(Fv, se, 1.96), ref_cif_interval(Fv, se, 1.96), tolerance = 1e-14)
  x <- sort(round(runif(30, 0, 10), 1)); y <- cumsum(runif(30)); tt <- c(-1, 0, x[5], runif(50, 0, 12), NA)
  expect_identical(SurvX:::.step_at(x, y, tt, 0.25), ref_step_at(x, y, tt, 0.25))
  time <- round(rexp(300, 0.3), 1) + 0.1; status <- rbinom(300, 2, 0.4); entry <- ifelse(runif(300) < 0.3, time * runif(300), 0)
  for (k in 1:2) for (e in list(NULL, entry)) {
    a <- SurvX:::.na_increments(time, status, k, e); b <- ref_na_increments(time, status, k, e)
    expect_equal(a, b, tolerance = 0)
    expect_equal(SurvX:::.na_frame(a, 0.9), ref_na_frame(b, 0.9), tolerance = 1e-14)
  }
  S <- c(1, 0.9, 0.5, 0.2, 0); se <- c(0, 0.03, 0.05, 0.04, 0)
  fit <- km(Event(time, status) ~ 1, data.frame(time = time, status = as.integer(status > 0)))
  cv <- fit@curves[[1L]]
  ref <- ref_km_limits(cv$surv, cv$se, stats::qnorm(0.975))
  expect_equal(cv$lower, ref$lower, tolerance = 1e-14)
  expect_equal(cv$upper, ref$upper, tolerance = 1e-14)
})

test_that("frailty summaries agree with R", {
  set.seed(3)
  C <- 40; m <- 4; id <- rep(seq_len(C), each = m)
  z <- rgamma(C, 2, 2)[id]; x <- rnorm(C * m)
  t <- rweibull(C * m, 1.3, (0.2 * z * exp(0.5 * x))^(-1 / 1.3)); cens <- runif(C * m, 0, 10)
  d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, id = id)
  f <- cox(Event(time, status) ~ x + frailty(id), d)
  expect_equal(SurvX:::.frailty_frame(f, 0.9), ref_frailty_frame(f, 0.9), tolerance = 1e-11)
  expect_equal(SurvX:::.frailty_lrt(f@loglik, f@loglik0, f@boundary), ref_frailty_lrt(f@loglik, f@loglik0, f@boundary),
               tolerance = 1e-13)
  X <- cbind(x = c(-1, 0, 2)); tt <- c(0.5, 2, 8)
  expect_equal(unname(SurvX:::.frailty_cumhaz(f, X, tt)), ref_breslow_cumhaz(f@basehaz, X, f@coefficients, tt), tolerance = 1e-14)
  H <- ref_breslow_cumhaz(f@basehaz, X, f@coefficients, tt); th <- f@theta
  expect_equal(unname(predict(f, data.frame(x = X[, 1]), times = tt)), (1 + th * H)^(-1 / th), tolerance = 1e-14)
  expect_equal(unname(predict(f, data.frame(x = X[, 1]), times = tt, frailty = "conditional")), exp(-H), tolerance = 1e-14)
})

test_that("adjusted_risk numerics agree with R: treatment model, weights, balance, weighted AJ, bootstrap summaries", {
  set.seed(4)
  n <- 300
  d <- data.frame(x = rnorm(n), z = rbinom(n, 1, 0.4), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$arm <- factor(ifelse(runif(n) < plogis(0.5 * d$x - 0.3 * d$z), "B", ifelse(runif(n) < 0.5, "A", "C")), levels = c("A", "B", "C"))
  psf <- arm ~ x + z + g
  for (trim in c(0, 0.05)) {
    a <- SurvX:::.ar_ps(d, "arm", levels(d$arm), psf, trim); b <- ref_ar_ps(d, "arm", levels(d$arm), psf, trim)
    expect_equal(a$w, b$w, tolerance = 1e-10)
    expect_equal(a$diag$smd, b$diag$smd, tolerance = 1e-10)
    expect_equal(a$diag$weights, b$diag$weights, tolerance = 1e-10)
    expect_equal(unname(a$diag$ps_range), unname(b$diag$ps_range), tolerance = 1e-10)
  }
  time <- round(rexp(n, 0.3), 1) + 0.1; status <- rbinom(n, 2, 0.4); w <- runif(n, 0.3, 3)
  entry <- ifelse(runif(n) < 0.3, time * runif(n), 0); tt <- sort(runif(20, 0, 8))
  for (k in 1:2) for (e in list(NULL, entry))
    expect_equal(SurvX:::.ar_wcif(time, status, e, w, k, tt), ref_ar_wcif(time, status, e, w, k, tt), tolerance = 1e-13)
  lev <- c("A", "B"); tt <- c(1, 2, 3); est <- list(gformula = matrix(c(0.2, 0.3, 0.25, 0.4, 0.3, 0.5), 2))
  res <- lapply(1:50, function(b) list(gformula = est$gformula * exp(rnorm(6, 0, 0.1))))
  boot <- list(res = res, ok = 50L, B = 50L)
  expect_equal(SurvX:::.ar_tables(est, boot, lev, c(1, 3), tt, 0.9), ref_ar_tables(est, boot, lev, c(1, 3), tt, 0.9), tolerance = 1e-12)
})

test_that("Turnbull summaries agree with R: survival and bounds, curve, test, tables, medians", {
  set.seed(5)
  n <- 150
  lat <- rweibull(n, 1.3, 2); g <- factor(sample(c("a", "b", "c"), n, TRUE))
  L <- R <- numeric(n)
  for (i in seq_len(n)) {
    v <- cumsum(runif(8, 0.3, 1.2))
    L[i] <- max(c(0, v[v < lat[i]])); R[i] <- if (any(v >= lat[i])) min(v[v >= lat[i]]) else Inf
  }
  d <- data.frame(L = L, R = R, g = g)
  tb <- SurvX:::.tb_fit(L, R)
  tt <- c(0, sort(runif(30, 0, 10)), tb$lo[1:3])
  for (left in c(FALSE, TRUE)) expect_equal(SurvX:::.tb_surv(tb, tt, left), ref_tb_surv(tb, tt, left), tolerance = 1e-14)
  expect_equal(SurvX:::.tb_steps(tb), ref_tb_steps(tb), tolerance = 1e-14)
  expect_equal(SurvX:::.tb_test(L, R, g, 100000L, 1e-10), ref_tb_test(L, R, g, 100000L, 1e-10), tolerance = 1e-12)
  fit <- turnbull(Interval(L, R) ~ g, d, B = 30, seed = 2)
  for (times in list(fit@times, c(0.5, 1.37, 2.5, 30))) expect_equal(SurvX:::.tb_table(fit, times), ref_tb_table(fit, times), tolerance = 1e-12)
  expect_equal(SurvX:::.tb_medians(fit), ref_tb_medians(fit), tolerance = 1e-12)
  st <- SurvX:::.tb_steps(tb)
  for (p in c(0.25, 0.5, 0.9)) expect_equal(SurvX:::.tb_quant(st$time, st$surv, p), ref_tb_quant(st$time, st$surv, p), tolerance = 1e-14)
})

test_that("model comparison helpers agree with R: measures, curves, integral, ranking", {
  set.seed(6)
  n <- 250; time <- round(rexp(n, 0.3), 1) + 0.05; status <- rbinom(n, 2, 0.45)
  risk <- matrix(runif(n * 4), n, 4); tt <- c(1, 2, 3, 4)
  for (cause in 1:2) expect_equal(unname(SurvX:::.sx_measures(time, status, risk, tt, cause)), unname(ref_measures(time, status, risk, tt, cause)),
                                  tolerance = 1e-14)
  at <- c(0, sort(runif(20, 0, 10)))
  expect_equal(SurvX:::.scan_np(time, status, 2L, at), ref_scan_np(time, status, 2L, at), tolerance = 1e-14)
  s01 <- as.integer(status > 0)
  expect_equal(SurvX:::.scan_np(time, s01, 1L, at), ref_scan_np(time, s01, 1L, at), tolerance = 1e-14)
  y <- runif(10); g <- seq(0.5, 5, length.out = 10)
  expect_equal(SurvX:::.scan_integrate(g, y), ref_scan_integrate(g, y), tolerance = 1e-14)
  expect_identical(SurvX:::.scan_integrate(g, c(y[-1], NA)), NA_real_)
  M <- 7; folds <- 5
  tab <- data.frame(fitted = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE), ibs_cv = c(0.15, 0.14, NA, 0.149, 0.2, 0.14, 0.145))
  foldibs <- matrix(rep(tab$ibs_cv, folds) + rnorm(M * folds, 0, 0.01), M, folds); foldibs[2, 3] <- NA
  new <- tab
  rr <- .C("SX_rank_one_se", as.integer(M), as.integer(folds), as.double(ifelse(tab$fitted, tab$ibs_cv, NA)), SurvX:::.rm(foldibs),
           rank = integer(M), diff = double(M), se = double(M), equivalent = integer(M), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  ref <- ref_rank(tab, foldibs)
  fin <- is.finite(tab$ibs_cv)
  expect_equal(ifelse(fin, rr$rank, NA_integer_), ref$rank)
  expect_equal(ifelse(fin, rr$diff, NA_real_), ref$diff_best, tolerance = 1e-14)
  expect_equal(ifelse(fin & !is.nan(rr$se), rr$se, NA_real_), ref$se_diff, tolerance = 1e-13)
  expect_equal(ifelse(fin, rr$equivalent == 1L, NA), ref$equivalent)
})

test_that("tests of proportional hazards agree with R (risk sets with and without entry, strata, weights)", {
  set.seed(7)
  n <- 200
  d <- data.frame(x = rnorm(n), z = rbinom(n, 1, 0.5), s = sample(1:2, n, TRUE))
  d$time <- round(rexp(n, 0.2 * exp(0.5 * d$x)), 1) + 0.1; d$status <- rbinom(n, 1, 0.7)
  d$entry <- ifelse(runif(n) < 0.3, d$time * runif(n), 0); w <- runif(n, 0.5, 2)
  X <- as.matrix(d[c("x", "z")]); b <- c(0.4, -0.2); V <- diag(c(0.01, 0.02))
  for (tr in c("km", "rank", "identity", "log"))
    for (cfg in list(list(), list(w = w), list(entry = d$entry), list(strata = d$s - 1L), list(entry = d$entry, strata = d$s - 1L, w = w))) {
      a <- SurvX:::.zph_core(d$time, d$status, X, b, V, tr, cfg$w, cfg$entry, cfg$strata)
      r <- ref_zph_core(d$time, d$status, X, b, V, tr, cfg$w, cfg$entry, cfg$strata)
      expect_equal(a$table, r$table, tolerance = 1e-10)
      expect_equal(a$x, r$x, tolerance = 1e-12)
      expect_equal(unname(a$y), unname(r$y), tolerance = 1e-10)
      expect_equal(a$time, r$time, tolerance = 0)
    }
  for (weighted in c(FALSE, TRUE)) {
    a <- SurvX:::.te_zph(d$time, d$status, X, b, V, SurvX:::.te_gfun("km", d$time, d$status, NULL, w, weighted), w, NULL, NULL, c("x", "z"))
    r <- ref_te_zph(d$time, d$status, X, b, V, ref_te_gfun("km", d$time, d$status, NULL, w, weighted), w, NULL, NULL, c("x", "z"))
    expect_equal(a$terms, r$terms, tolerance = 1e-10)
    expect_equal(a$columns, r$columns, tolerance = 1e-10)
    expect_equal(a$global, r$global, tolerance = 1e-10)
  }
  p <- c(0.01, NA, 0.04, 0.03, 0.2, 0.03)
  expect_equal(SurvX:::.sx_p_holm(p), stats::p.adjust(p, "holm"), tolerance = 1e-15)
})

test_that("time_effects numerics agree with R: time functions (splines::ns), splitting, effects", {
  set.seed(8)
  et <- sort(unique(round(rexp(60, 0.3), 2) + 0.05)); tt <- c(et[1] / 4, sort(runif(40, 0, max(et) * 1.3)), range(et))
  for (fs in list(list(name = "log", type = "log"), list(name = "linear", type = "linear"), list(name = "spline", type = "spline")))
    for (df in 1:4) {
      a <- SurvX:::.te_basis(fs, et, df, NULL); r <- ref_te_basis(fs, et, df, NULL)
      expect_equal(unname(a$B(tt)), unname(r$B(tt)), tolerance = 1e-12)
    }
  fs <- list(name = "piecewise", type = "piecewise")
  a <- SurvX:::.te_basis(fs, et, 3, c(1, 2.5)); r <- ref_te_basis(fs, et, 3, c(1, 2.5))
  expect_equal(unname(a$B(tt)), unname(r$B(tt)), tolerance = 0)
  n <- 80; stop <- round(rexp(n, 0.3), 1) + 0.1; start <- ifelse(runif(n) < 0.3, stop * runif(n), 0); status <- rbinom(n, 1, 0.6)
  ev <- sort(unique(stop[status == 1]))
  expect_equal(SurvX:::.te_split(start, stop, status, ev), ref_te_split(start, stop, status, ev))
  expect_equal(SurvX:::.te_split_size(start, stop, ev), ref_te_split_size(start, stop, ev))
  b <- c(x = 0.3, z = -0.2, "x:log" = 0.1); V <- crossprod(matrix(rnorm(9), 3)) / 50
  bl <- SurvX:::.te_basis(list(name = "log", type = "log"), et, 3, NULL)
  map <- list(list(col = 1L, form = "log", tv = 3L, basis = bl), list(col = 2L, form = "constant", tv = integer(0), basis = NULL))
  fin <- list(coefficients = b)
  expect_equal(SurvX:::.te_effects(fin, V, map, c("x", "z"), c(0.5, 1, 3), 0.9),
               ref_te_effects(fin, V, map, c("x", "z"), c(0.5, 1, 3), 0.9), tolerance = 1e-12)
})

test_that("Fine-Gray expansion, timefix and robust variance agree with R", {
  set.seed(9)
  for (rep in 1:20) {
    n <- 40; time <- round(rexp(n, 0.3) + 0.1, 1); status <- sample(0:3, n, TRUE)
    entry <- if (rep %% 2) pmin(ifelse(runif(n) < 0.4, round(time * runif(n), 1), 0), time - 0.05) else NULL
    a <- SurvX:::.fg_expand(time, status, entry, 1L); r <- ref_fg_expand(time, status, entry, 1L)
    expect_equal(as.list(a), as.list(r), tolerance = 1e-14)
  }
  y <- Event(c(1, 1 + 1e-12, 2, 3, 3 + 5e-9, 4), c(1, 0, 1, 1, 0, 1), timefix = FALSE)
  expect_equal(unclass(SurvX:::.timefix(y)), unclass(ref_timefix(y)))
  d <- data.frame(time = round(rexp(60), 2) + 0.01, status = rbinom(60, 1, 0.7), x = rnorm(60), z = rnorm(60), cl = rep(1:15, 4))
  X <- as.matrix(d[, c("x", "z")]); f <- SurvX:::.cox_fit(d$time, d$status, X, NULL, "efron", 50L, 1e-9, NULL, NULL)
  for (cl in list(NULL, d$cl))
    expect_equal(SurvX:::.cox_robust(f, d$time, d$status, X, NULL, "efron", NULL, NULL, cl),
                 ref_cox_robust(f, d$time, d$status, X, NULL, "efron", NULL, NULL, cl), tolerance = 1e-12)
})

test_that("cause-specific cumulative incidence with strata agrees with R", {
  d <- sim_cr(300, 4); d$g <- rep(c("a", "b", "c"), 100)
  for (pl in c(TRUE, FALSE)) {
    f <- csc(Event(time, status) ~ x + z + strata(g), d, product_limit = pl)
    Xn <- f@X[1:25, , drop = FALSE]; sn <- f@strata[1:25]; tt <- c(0, 0.5, 1, 2, 5, 50)
    for (k in 1:2)
      expect_equal(SurvX:::.csc_cif_general(f, Xn, sn, tt, k), ref_csc_cif_general(f, Xn, sn, tt, k), tolerance = 1e-12)
  }
})

test_that("tt() splitting agrees with R", {
  set.seed(10); n <- 50
  time <- round(rexp(n, 0.4), 1) + 0.1; status <- rbinom(n, 1, 0.7)
  for (ent in list(NULL, ifelse(runif(n) < 0.3, round(time * runif(n), 1), 0))) {
    f <- list(time = time, status = status, X = cbind(x = rnorm(n), z = rnorm(n)), entry = ent, tt_cols = 2L, strata = NULL, cluster = NULL)
    expect_equal(SurvX:::.tt_split(f, NULL, function(t) log(t)), ref_tt_split(f, NULL, function(t) log(t)))
  }
})

test_that("every random draw comes from R's generator", {
  # a resample is sample.int(n, n, TRUE), row after row, after set.seed(seed)
  set.seed(5); ref <- t(replicate(4, sample.int(9, 9, TRUE)))
  set.seed(99); before <- .Random.seed
  expect_identical(SurvX:::.sx_resample(9L, 4L, 5L), ref)
  expect_identical(.Random.seed, before)                       # a numeric seed leaves the session alone
  # folds: sample(rep_len(1:folds, m)) for each stratum in order of first appearance
  st <- c(1, 0, 1, 1, 0, 2, 1, 0, 2, 1, 1)
  set.seed(3); ref <- integer(length(st))
  for (s in unique(st)) { i <- which(st == s); ref[i] <- sample(rep_len(1:3, length(i))) }
  expect_identical(SurvX:::.sx_folds(st, 3L, 3L), ref)
  # seed = NULL: set.seed() before the call governs the draws
  set.seed(11); a <- SurvX:::.resolve_seed(NULL); set.seed(11); b <- SurvX:::.resolve_seed(NULL)
  expect_identical(a, b)
  expect_error(SurvX:::.resolve_seed(1.5), "whole number")
  d <- data.frame(L = c(0, 1, 2, 2, 3, 0, 4, 1), R = c(2, 3, Inf, 4, 5, 1, Inf, 2), g = rep(c("a", "b"), 4))
  set.seed(7); f1 <- turnbull(Interval(L, R) ~ g, d, B = 20, seed = NULL, test = FALSE)
  set.seed(7); f2 <- turnbull(Interval(L, R) ~ g, d, B = 20, seed = NULL, test = FALSE)
  expect_identical(f1@boot, f2@boot)
  # the C++ bootstrap of cif() draws the same subjects as sample.int
  dd <- sim_cr(120, 2); fit <- csc(Event(time, status) ~ x, dd)
  set.seed(1); a <- cif_data(fit, times = 1, B = 30, seed = 4); b <- cif_data(fit, times = 1, B = 30, seed = 4)
  expect_identical(a, b)
})

test_that("residuals of parametric and cure models agree with R", {
  set.seed(12); n <- 150
  x <- rnorm(n); z <- rbinom(n, 1, 0.5); t0 <- rweibull(n, 1.4, exp(0.3 * x)); cens <- runif(n, 0, 3)
  cured <- runif(n) < plogis(-0.5 + z); t0[cured] <- Inf
  d <- data.frame(time = pmin(t0, cens), status = as.integer(t0 <= cens), x = x, z = z)
  d$entry <- ifelse(runif(n) < 0.3, d$time * runif(n) * 0.5, 0)
  for (b in c("weibull", "lognormal", "rp")) for (en in c(FALSE, TRUE)) {
    f <- if (en) parsurv(Event(time, status, entry = entry) ~ x + z, d, baseline = b) else parsurv(Event(time, status) ~ x + z, d, baseline = b)
    for (ty in c("martingale", "deviance", "coxsnell")) expect_equal(residuals(f, ty), ref_par_resid(f, ty), tolerance = 1e-12)
  }
  for (ty in c("mixture", "nonmixture")) {
    f <- cure(Event(time, status) ~ x, d, cure = ~ z, type = ty)
    for (rt in c("martingale", "deviance", "coxsnell")) expect_equal(residuals(f, rt), ref_cure_resid(f, rt), tolerance = 1e-12)
  }
})

test_that("bootstrap summaries, covariance, smoother, Fine-Gray predictions and log limits agree with R", {
  set.seed(13)
  bs <- array(runif(3 * 4 * 40), c(3, 4, 40)); bs[, , c(5, 17)] <- NA
  s <- SurvX:::.sx_boot_summary(bs, 40L, 12L, 0.9, finite_only = TRUE)
  expect_equal(matrix(s$sd, 3, 4), apply(bs, c(1L, 2L), stats::sd, na.rm = TRUE), tolerance = 1e-12)
  expect_equal(matrix(s$lower, 3, 4), apply(bs, c(1L, 2L), stats::quantile, probs = 0.05, na.rm = TRUE, names = FALSE), tolerance = 1e-12)
  est <- matrix(rnorm(60), 20); est[c(3, 9), 2] <- NA
  r <- .C("SX_cov_complete", 20L, 3L, as.double(t(est)), out = double(9), nok = integer(1), rc = integer(1), NAOK = TRUE, PACKAGE = "SurvX")
  expect_equal(matrix(r$out, 3), stats::cov(est[stats::complete.cases(est), ]), tolerance = 1e-12)
  expect_identical(r$nok, 18L)
  x <- runif(80, 0, 3); y <- sin(2 * x) + rnorm(80, 0, 0.3)
  for (k in 1:4) {
    sm <- SurvX:::.ns_smooth(x, y, k, 0.95)
    fit <- stats::lm(y ~ splines::ns(x, df = k)); pr <- stats::predict(fit, data.frame(x = sm$x), se.fit = TRUE)
    expect_equal(sm$fit, unname(pr$fit), tolerance = 1e-10)
    expect_equal(sm$upper, unname(pr$fit + stats::qnorm(0.975) * pr$se.fit), tolerance = 1e-10)
  }
  d <- sim_cr(200, 5); f <- fg(Event(time, status) ~ x + z, d)
  expect_equal(predict(f, d[1:15, ], times = c(0.5, 1, 3)), ref_fg_predict(f, d[1:15, ], times = c(0.5, 1, 3)), tolerance = 1e-12)
  e <- matrix(c(0.2, 1, 3, -1, NA, 2), 2, dimnames = list(NULL, c("a", "b", "c"))); se <- matrix(c(0.1, 0.3, NA, 0.2, 0.1, 0.5), 2)
  expect_equal(SurvX:::.pred_limits(e, se, "hazard", 0.9), ref_pred_limits(e, se, "hazard", 0.9))
})

test_that("the toolkit of the core gives R's own numbers", {
  set.seed(14)
  x <- c(rnorm(20), NA, 0); y <- c(runif(20) + 0.1, 2, NA); p <- c(runif(20), 0, 1)
  ew <- SurvX:::.sx_ew
  expect_equal(ew("exp", x), exp(x)); expect_equal(ew("log", y), log(y)); expect_equal(ew("plogis", x), stats::plogis(x))
  expect_equal(ew("logit", p), log(p / (1 - p))); expect_equal(ew("qnorm", p), stats::qnorm(p), tolerance = 1e-14)
  expect_equal(ew("complement", x), 1 - x); expect_equal(ew("neglog", y), -log(y)); expect_equal(ew("sub", x, y), x - y)
  expect_equal(ew("div", x, y), x / y); expect_equal(ew("logratio", y, y + 1), log(y / (y + 1)))
  expect_equal(ew("pvalue", x, y), 2 * stats::pnorm(-abs(x / y)), tolerance = 1e-14)
  expect_equal(ew("exp_se", x, y), exp(x) * sqrt(y)); expect_equal(ew("neglog_floor", p), -log(pmax(p, 1e-300)))
  expect_equal(ew("mul", x, y), x * y); expect_equal(ew("add", x, y), x + y); expect_equal(ew("sqrt", y), sqrt(y))
  expect_equal(ew("qnorm_comp", p), stats::qnorm(1 - p), tolerance = 1e-14); expect_equal(ew("logit_comp", p), log((1 - p) / p))
  expect_equal(ew("exp_neg", x), exp(-x)); expect_equal(ew("absdiff", x, y), abs(x - y))
  L <- c(0, 1, 2, 3); R <- c(2, Inf, 4, 3)
  expect_equal(ew("ic_mid", L, R), ifelse(!is.finite(R), L, ifelse(L == 0, R / 2, (L + R) / 2)))
  M <- matrix(c(runif(29), NA), 10, dimnames = list(NULL, c("a", "b", "c")))
  expect_identical(SurvX:::.sx_reduce("sum", x[1:20]), sum(x[1:20])); expect_identical(SurvX:::.sx_reduce("mean", x, TRUE), mean(x, na.rm = TRUE))
  expect_identical(SurvX:::.sx_reduce("maxabs", x, TRUE), max(abs(x), na.rm = TRUE)); expect_identical(SurvX:::.sx_reduce("min", x, TRUE), min(x, na.rm = TRUE))
  expect_identical(SurvX:::.sx_colmeans(M), colMeans(M)); expect_identical(SurvX:::.sx_colmeans(M, TRUE), colMeans(M, na.rm = TRUE))
  expect_identical(SurvX:::.sx_rowmeans(M, TRUE), rowMeans(M, na.rm = TRUE))
  g <- sample(c(2, 5, 9), 30, TRUE); v <- c(runif(29), NA); a <- stats::aggregate(v ~ g, FUN = mean)
  gm <- SurvX:::.sx_group_means(v, g); expect_equal(gm$group, a$g); expect_identical(gm$mean, a$v)
  M2 <- matrix(runif(12), 4); expect_identical(SurvX:::.sx_cumsum_rows(M2), t(apply(M2, 1L, cumsum)))
  expect_identical(SurvX:::.sx_cumsum_rows(M2, 1), pmin(t(apply(M2, 1L, cumsum)), 1))
  X <- matrix(rnorm(20), 10); b <- c(0.3, -0.7); H <- c(0, 0.2, 1.5)
  expect_equal(SurvX:::.sx_linpred(X, b), drop(X %*% b), tolerance = 1e-15)
  expect_equal(SurvX:::.sx_ph_predict(X, b, H), exp(-outer(exp(drop(X %*% b)), H)), tolerance = 1e-14)
  expect_equal(SurvX:::.sx_ph_predict(X, b, H, "cif"), 1 - exp(-outer(exp(drop(X %*% b)), H)), tolerance = 1e-14)
  d <- c(1, 2, 1, 1); Y <- c(10, 8, 5, 2)
  expect_identical(SurvX:::.sx_np_curve_inc(d, Y, "na"), cumsum(d / Y)); expect_identical(SurvX:::.sx_np_curve_inc(d, Y, "km"), cumprod(1 - d / Y))
  ic <- SurvX:::.sx_ic(-120.5, 4, 80); expect_identical(ic$aic, -2 * -120.5 + 2 * 4); expect_identical(ic$bic, -2 * -120.5 + log(80) * 4)
  lr <- SurvX:::.sx_lr(-100, -104.2, 2); expect_equal(lr$stat, 8.4); expect_equal(lr$p, stats::pchisq(8.4, 2, lower.tail = FALSE), tolerance = 1e-13)
  V <- crossprod(matrix(rnorm(9), 3)); expect_equal(SurvX:::.sx_se(V), sqrt(diag(V)))
  w <- SurvX:::.sx_wald(c(0.2, -1), c(0.1, 0.4), 0.95, c(TRUE, FALSE)); zq <- stats::qnorm(0.975)
  expect_equal(w$lo, c(0.2, -1) - zq * c(0.1, 0.4), tolerance = 1e-14); expect_equal(w$eest, c(exp(0.2), NA))
  for (up in c(1, 3.65, 7.3, 10, 59.9, 137, 2500, 0.034)) expect_equal(SurvX:::.sx_pretty(0, up, 5L), pretty(c(0, up), n = 5L))
  expect_identical(SurvX:::.sx_seq(0, 7.3, 201L), seq(0, 7.3, length.out = 201L))
  expect_identical(SurvX:::.sx_seq_by(0.2, 0.8, 0.1), seq(0.2, 0.8, by = 0.1))
  expect_identical(SurvX:::.sx_upper_prob(0.95), 1 - (1 - 0.95) / 2)
  st <- rbinom(15, 1, 0.6); m <- ifelse(st == 1L, runif(15, -2, 0.9), -runif(15, 0, 2))   # censored: M <= 0
  expect_equal(SurvX:::.sx_hazard_resid(m, st, "deviance_from_m"), sign(m) * sqrt(-2 * (m + ifelse(st == 1L, log(pmax(st - m, 1e-300)), 0))))
  bh <- data.frame(time = c(0.5, 1, 2, 0.7, 1.5), cumhaz = c(0.1, 0.3, 0.4, 0.2, 0.6), stratum = c(0L, 0L, 0L, 1L, 1L))
  tj <- sort(bh$time); rs <- c(0L, 1L, 1L, 0L); lp <- rnorm(4)
  ref <- matrix(0, 4, 5)
  for (s in 0:1) { bs <- bh[bh$stratum == s, ]; i <- which(rs == s); ref[i, match(bs$time, tj)] <- outer(exp(lp[i]), diff(c(0, bs$cumhaz))) }
  expect_equal(SurvX:::.sx_hazard_inc(bh, rs, tj, lp), ref, tolerance = 1e-15)
})

test_that("predictions, residuals and summaries built on the toolkit agree with the earlier R code", {
  d <- sim_cr(200, 6); d$g <- rep(c("a", "b"), 100); d$dead <- as.integer(d$status > 0)
  f <- cox(Event(time, dead) ~ x + z + strata(g), d)
  bh <- SurvX:::.cox_baseline(f); X <- f@X[1:12, ]; sn <- f@strata[1:12]; tt <- c(0.3, 1, 2.5)
  ref <- t(vapply(1:12, function(i) exp(-exp(sum(X[i, ] * f@coefficients)) * SurvX:::.cumhaz_at(bh, tt, rep(sn[i], 3))), numeric(3)))
  expect_equal(unname(predict(f, d[1:12, ], times = tt, type = "survival")), ref, tolerance = 1e-13)
  f0 <- cox(Event(time, dead) ~ x + z, d)
  for (ty in c("martingale", "deviance", "coxsnell")) {
    m <- residuals(f0, "martingale")
    ref <- switch(ty, martingale = m, coxsnell = d$dead - m, deviance = sign(m) * sqrt(-2 * (m + ifelse(d$dead == 1L, log(pmax(d$dead - m, 1e-300)), 0))))
    expect_equal(residuals(f0, ty), ref, tolerance = 1e-14)
  }
  sch <- residuals(f0, "schoenfeld"); sc <- residuals(f0, "scaledsch"); nd <- nrow(sch)
  expect_equal(unclass(sc)[, ], (sch %*% f0@vcov * nd + matrix(f0@coefficients, nd, 2, byrow = TRUE))[, ], tolerance = 1e-13)
})

test_that("signif() and round() of the core are R's", {
  x <- c(1.23456, 0.0012345, 123456, 2.5, -0.015, 0.15, 9.995, 1e-310, 0, NA, -7.77e5, 3.14159e-7)
  for (d in 1:5) expect_identical(SurvX:::.sx_signif(x, d), signif(x, d))
  r <- c(0.5, 1.5, 2.5, -0.5, -1.5, 3.4999, 12.5, NA)
  expect_identical(SurvX:::.sx_round(r), round(r))
})

test_that(".df() gives the same object as data.frame(), and km/aj tables keep their form", {
  m <- matrix(1:6 + 0.5, 3, dimnames = list(NULL, c("a.1", "a.2")))
  expect_identical(SurvX:::.df(list(x = c(1, 2, 3), g = c("u", "v", "w"), m = m)),
                   data.frame(x = c(1, 2, 3), g = c("u", "v", "w"), m, stringsAsFactors = FALSE))
  d <- sim_cr(150, 3); d$grp <- rep(c("a", "b", "c"), 50); d$dead <- as.integer(d$status > 0)
  k <- km(Event(time, dead) ~ grp, d)
  expect_identical(names(k@counts), c("group", "n", "events", "median", "median.lower", "median.upper"))
  expect_identical(k@counts$n, c(50L, 50L, 50L))
  a <- aj(Event(time, status) ~ grp, d)
  expect_identical(names(a@counts), c("group", "n", "censored", paste0("cause.", 1:2)))
  expect_identical(a@counts$censored + a@counts$cause.1 + a@counts$cause.2, a@counts$n)
  expect_identical(names(a@curves[[1]])[1:4], c("time", "surv", "cif.1", "cif.2"))
})
