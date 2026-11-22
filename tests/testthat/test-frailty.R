# Shared gamma frailty: cox(... + frailty(id)) and parsurv(... + frailty(id))

fr_data <- function(C = 80, m = 4, theta = 0.5, seed = 1) {
  set.seed(seed)
  id <- rep(seq_len(C), each = m)
  z <- if (theta > 0) rgamma(C, shape = 1 / theta, rate = 1 / theta)[id] else rep(1, C * m)
  x <- rnorm(C * m); trt <- rbinom(C * m, 1, 0.5)
  t <- rweibull(C * m, 1.3, (0.1 * z * exp(0.5 * x - 0.7 * trt))^(-1 / 1.3))
  cens <- runif(C * m, 0, 15)
  data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, trt = trt, id = id)
}

# marginal log-likelihood of a parametric frailty fit, computed from its predictions
fr_loglik <- function(f, d) {
  X <- f@X; th <- f@theta
  H <- vapply(seq_len(nrow(d)), function(i) SurvX:::.par_predict_X(f@base, X[i, , drop = FALSE], d$time[i], "cumhaz", FALSE, 0.95)[1, 1], 1)
  h <- vapply(seq_len(nrow(d)), function(i) SurvX:::.par_predict_X(f@base, X[i, , drop = FALSE], d$time[i], "hazard", FALSE, 0.95)[1, 1], 1)
  ll <- sum(d$status * log(h))
  for (c in unique(d$id)) {
    k <- d$id == c; D <- sum(d$status[k]); L <- sum(H[k])
    ll <- ll + D * log(th) + lgamma(1 / th + D) - lgamma(1 / th) - (1 / th + D) * log1p(th * L)
  }
  ll
}

test_that("cox() with frailty() is the marginal maximum likelihood estimate (EM with coxph offsets)", {
  skip_if_not_installed("survival")
  d <- fr_data()
  f <- cox(Event(time, status) ~ x + trt + frailty(id), data = d)
  expect_true(S7::S7_inherits(f, SurvX:::FrailtyFit))
  expect_true(f@converged)
  th <- f@theta
  # EM of Klein (1992) at the estimated theta, with survival::coxph for the M-step
  a <- rep(1, max(d$id)); D <- as.vector(tapply(d$status, d$id, sum))
  for (it in 1:1000) {
    r <- survival::coxph(survival::Surv(time, status) ~ x + trt + offset(log(a[id])), data = d, ties = "breslow")
    e <- exp(drop(as.matrix(d[, c("x", "trt")]) %*% coef(r)) + log(a[d$id]))
    tk <- sort(unique(d$time[d$status == 1]))
    hk <- vapply(tk, function(u) sum(d$status[d$time == u]) / sum(e[d$time >= u]), 1)
    H <- vapply(d$time, function(u) sum(hk[tk <= u]), 1)
    L <- as.vector(tapply(exp(drop(as.matrix(d[, c("x", "trt")]) %*% coef(r))) * H, d$id, sum))
    an <- (1 + th * D) / (1 + th * L)
    if (max(abs(an - a)) < 1e-10) { a <- an; break }
    a <- an
  }
  expect_close(coef(f), coef(r), 1e-5)
  # penalized likelihood of survival::coxph: close to the marginal estimates
  rp <- survival::coxph(survival::Surv(time, status) ~ x + trt + survival::frailty(id, distribution = "gamma"),
                        data = d, ties = "breslow")
  expect_lt(max(abs(coef(f) - coef(rp)[1:2])), 0.1)
  expect_lt(abs(log(f@theta / rp$history[[1]]$theta)), 0.3)
  expect_gt(f@lrt$statistic, 0)
  expect_output(print(f), "Shared gamma frailty")
  expect_output(print(summary(f)), "lowest and the highest")
})

test_that("the Cox frailty fit without heterogeneity is close to the Cox model", {
  d <- fr_data(theta = 0, seed = 2)
  f <- cox(Event(time, status) ~ x + trt + frailty(id), data = d)
  c0 <- cox(Event(time, status) ~ x + trt, data = d, ties = "breslow")
  expect_lt(f@theta, 0.2)
  expect_close(coef(f), coef(c0), 0.05)
  expect_gte(f@loglik, f@loglik0 - 1e-6)
})

test_that("parsurv() with frailty() maximises the marginal likelihood", {
  d <- fr_data(60, 4, 0.5, 3)
  for (sp in list(c("weibull", "PH"), c("weibull", "AFT"), c("lognormal", "AFT"), c("gompertz", "PH"))) {
    f <- parsurv(Event(time, status) ~ x + trt + frailty(id), data = d, baseline = sp[1], scale = sp[2])
    expect_true(f@converged)
    expect_close(fr_loglik(f, d), as.numeric(logLik(f)), 1e-7)
    expect_identical(attr(logLik(f), "df"), length(coef(f)) + 1L)
    # no better point nearby: the likelihood decreases along every coordinate
    par <- c(coef(f), log(f@theta))
    for (j in seq_along(par)) for (s in c(-1, 1)) {
      q <- par; q[j] <- q[j] + s * 1e-3
      b <- f@base; b@coefficients <- q[-length(q)]
      g <- f; g@base <- b; g@theta <- unname(exp(q[length(q)]))
      expect_lte(fr_loglik(g, d), as.numeric(logLik(f)) + 1e-8)
    }
  }
  # Weibull PH and AFT are the same model
  a <- parsurv(Event(time, status) ~ x + trt + frailty(id), data = d, baseline = "weibull", scale = "PH")
  b <- parsurv(Event(time, status) ~ x + trt + frailty(id), data = d, baseline = "weibull", scale = "AFT")
  expect_close(as.numeric(logLik(a)), as.numeric(logLik(b)), 1e-8)
  expect_close(a@theta, b@theta, 1e-6)
})

test_that("frailty predictions, cluster estimates and messages", {
  d <- fr_data(50, 4, 0.8, 4)
  f <- cox(Event(time, status) ~ x + trt + frailty(id), data = d)
  nd <- data.frame(x = 0, trt = 0:1)
  m <- predict(f, nd, times = c(2, 5, 10), frailty = "marginal")
  cnd <- predict(f, nd, times = c(2, 5, 10), frailty = "conditional")
  expect_true(all(m >= cnd - 1e-12))                       # Jensen: the marginal curve lies above
  ndc <- data.frame(x = 0, trt = 0, id = c(1, 2))
  pc <- predict(f, ndc, times = 5, frailty = "cluster")
  expect_identical(dim(pc), c(2L, 1L))
  fr <- as.data.frame(f, what = "frailty")
  expect_identical(nrow(fr), 50L)
  expect_true(all(fr$lower <= fr$frailty & fr$frailty <= fr$upper))
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(f))
  expect_silent(plot(f, which = "survival"))
  expect_error(fg(Event(time, status) ~ x + frailty(id), data = d), "frailty")
  expect_error(cox(Event(time, status) ~ x + frailty(id), data = d, weights = rep(1, nrow(d))), "weights")
  expect_error(cox(Event(time, status) ~ x + frailty(id) + cluster(id), data = d), "cannot be combined")
})
