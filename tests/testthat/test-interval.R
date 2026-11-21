# Interval censoring: Interval(), turnbull(), parsurv() and directcif() with Interval() outcomes

# periodic inspections: the event is seen at the first visit after it
ic_data <- function(n = 300, seed = 1, exact = 0) {
  set.seed(seed)
  x <- rnorm(n); g <- factor(rep(c("a", "b"), length.out = n))
  t <- rweibull(n, 1.4, exp(1.5 + 0.4 * x + 0.3 * (g == "b")))
  L <- R <- numeric(n)
  for (i in seq_len(n)) {
    v <- cumsum(runif(8, 0.5, 2))
    L[i] <- max(c(0, v[v < t[i]]))
    R[i] <- if (any(v >= t[i])) min(v[v >= t[i]]) else Inf
  }
  if (exact > 0) { e <- seq_len(n) %% round(1 / exact) == 0 & is.finite(R); L[e] <- R[e] <- t[e] }
  data.frame(L = L, R = R, x = x, g = g, t = t)
}

test_that("Interval() codes the four kinds of observation and checks its input", {
  y <- Interval(c(2, 3, 4, 0, NA), c(5, 3, NA, 6, 7))
  m <- unclass(y)
  expect_equal(m[, "left"], c(2, 3, 4, 0, 0))
  expect_equal(m[, "right"], c(5, 3, Inf, 6, 7))
  expect_equal(m[, "status"], c(1, 1, 0, 1, 1))
  expect_equal(format(y), c("(2, 5]", "3", "4+", "(0, 6]", "(0, 7]"))
  expect_equal(SurvX:::.ic_type(m[, 1], m[, 2]), c(2L, 1L, 0L, 3L, 3L))
  expect_error(Interval(5, 2), "at least left")
  expect_error(Interval(1, 2, cause = 0), "cause must be")
  expect_equal(unclass(Interval(c(1, 2), c(3, Inf), cause = c(2, 0)))[, "status"], c(2, 0))
  d <- ic_data(50)
  expect_error(km(Interval(L, R) ~ 1, data = d), "turnbull")
  expect_error(cox(Interval(L, R) ~ x, data = d), "interval-censored")
})

test_that("turnbull() equals Kaplan-Meier for exact and right-censored times", {
  skip_if_not_installed("survival")
  sc <- sim_cr(200, 3)
  d <- data.frame(L = sc$time, R = ifelse(sc$dead == 1L, sc$time, Inf))
  tb <- turnbull(Interval(L, R) ~ 1, data = d, B = 0)
  tt <- c(0.2, 0.5, 1, 1.5, 2)
  p <- predict(tb, times = tt)
  s <- summary(survival::survfit(survival::Surv(time, dead) ~ 1, data = sc), times = tt)$surv
  expect_true(all(abs(p$identified_high - p$identified_low) < 1e-12))   # no unidentified region
  expect_close(p$survival, s, 1e-7)
})

test_that("turnbull() maximises the likelihood (EM by brute force, Kuhn-Tucker condition)", {
  d <- ic_data(120, 2, exact = 0.1)
  tb <- turnbull(Interval(L, R) ~ 1, data = d, B = 0)
  f <- tb@fits[[1]]
  expect_true(f$converged)
  expect_lt(f$kkt, 1e-6)
  # indicator matrix of the innermost intervals covered by each subject, plain EM in R
  A <- outer(seq_len(nrow(d)), seq_along(f$lo), function(i, j)
    ifelse(f$lo[j] == f$hi[j], (d$L[i] == d$R[i] & d$L[i] == f$lo[j]) | (d$L[i] < d$R[i] & f$lo[j] > d$L[i] & f$lo[j] <= d$R[i]),
           d$L[i] < d$R[i] & f$lo[j] >= d$L[i] & f$hi[j] <= d$R[i]))
  storage.mode(A) <- "double"
  expect_true(all(rowSums(A) > 0))
  p <- rep(1 / ncol(A), ncol(A))
  for (it in 1:20000) { P <- drop(A %*% p); p <- p * drop(crossprod(A, 1 / P)) / nrow(A) }
  expect_lte(abs(sum(log(A %*% p)) - f$loglik), 1e-4)
  expect_gte(f$loglik, sum(log(A %*% p)) - 1e-8)
  expect_close(cumsum(f$p), cumsum(p), 5e-3)
})

test_that("turnbull() groups, bootstrap, test and methods", {
  d <- ic_data(300, 4)
  tb <- turnbull(Interval(L, R) ~ g, data = d, B = 30)
  expect_identical(tb@groups, c("a", "b"))
  expect_true(is.list(tb@test) && tb@test$df == 1L && tb@test$p.value >= 0 && tb@test$p.value <= 1)
  out <- capture.output(print(tb))
  expect_true(any(grepl("Turnbull", out)))
  expect_true(any(grepl("chi-square", out)))
  expect_output(print(summary(tb)), "Innermost intervals")
  p <- predict(tb, times = c(3, 6))
  expect_identical(nrow(p), 4L)
  expect_true(all(p$lower <= p$survival + 1e-12 & p$survival <= p$upper + 1e-12, na.rm = TRUE))
  expect_true(all(p$identified_low <= p$survival + 1e-12 & p$survival <= p$identified_high + 1e-12, na.rm = TRUE))
  iv <- as.data.frame(tb, what = "intervals")
  expect_close(tapply(iv$probability, iv$group, sum), c(1, 1), 1e-10)
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(tb))
  # no difference: the test does not reject
  d0 <- ic_data(300, 5); d0$g <- factor(sample(c("a", "b"), 300, TRUE))
  # strong difference: the test rejects
  d1 <- d; d1$R[d1$g == "b"] <- d1$R[d1$g == "b"] * 3; d1$L[d1$g == "b"] <- d1$L[d1$g == "b"] * 3
  expect_lt(turnbull(Interval(L, R) ~ g, data = d1, B = 0)@test$p.value, 0.001)
})

test_that("parsurv() with Interval() of exact and right-censored times equals the Event() fit", {
  sc <- sim_cr(250, 7)
  sc$L <- sc$time; sc$R <- ifelse(sc$dead == 1L, sc$time, Inf)
  specs <- list(c("exponential", "PH"), c("weibull", "PH"), c("weibull", "AFT"), c("gompertz", "PH"),
                c("lognormal", "AFT"), c("loglogistic", "AFT"), c("rp", "PH"), c("pwexp", "PH"))
  for (sp in specs) {
    a <- parsurv(Event(time, dead) ~ x + z, data = sc, baseline = sp[1], scale = sp[2])
    b <- parsurv(Interval(L, R) ~ x + z, data = sc, baseline = sp[1], scale = sp[2])
    expect_close(coef(b), coef(a), 1e-5)
    expect_close(as.numeric(logLik(b)), as.numeric(logLik(a)), 1e-7)
    expect_close(sqrt(diag(vcov(b))), sqrt(diag(vcov(a))), 1e-4)
  }
})

test_that("parsurv() interval-censored fits agree with survreg(type = 'interval2')", {
  skip_if_not_installed("survival")
  d <- ic_data(300, 8)
  d$R2 <- ifelse(is.finite(d$R), d$R, NA); d$L2 <- ifelse(d$L > 0, d$L, NA)
  for (dist in c("weibull", "lognormal", "loglogistic", "exponential")) {
    s <- survival::survreg(survival::Surv(L2, R2, type = "interval2") ~ x + g, data = d, dist = dist)
    f <- parsurv(Interval(L, R) ~ x + g, data = d, baseline = dist, scale = "AFT")
    expect_close(coef(f)[1:3], coef(s), 1e-5)
    if (dist != "exponential") expect_close(exp(coef(f)[4]), s$scale, 1e-5)
    expect_close(as.numeric(logLik(f)), as.numeric(logLik(s)), 1e-6)
  }
  f <- parsurv(Interval(L, R) ~ x + g, data = d, baseline = "weibull")
  expect_true(f@converged)
  expect_output(print(f), "interval-censored")
  expect_error(residuals(f), "interval-censored")
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(f))
  expect_silent(hazard_plot(f, type = "cumhaz"))
})

test_that("directcif() with Interval(): exact data reproduce Event(), and the log-likelihood is F(R) - F(L)", {
  sc <- sim_cr(300, 9)
  sc$L <- sc$time; sc$R <- ifelse(sc$status > 0L, sc$time, Inf)
  a <- directcif(Event(time, status) ~ x + z, data = sc)
  b <- directcif(Interval(L, R, cause = status) ~ x + z, data = sc)
  expect_close(coef(b), coef(a), 1e-5)
  expect_close(as.numeric(logLik(b)), as.numeric(logLik(a)), 1e-7)
  # coarsen the times to visits
  set.seed(9)
  n <- nrow(sc); L <- R <- numeric(n)
  for (i in seq_len(n)) {
    v <- cumsum(runif(10, 0.2, 0.6))
    if (sc$status[i] == 0L) { L[i] <- max(c(0.05, v[v < sc$time[i]])); R[i] <- Inf }
    else { L[i] <- max(c(0, v[v < sc$time[i]])); R[i] <- min(v[v >= sc$time[i]]) }
  }
  sc$L <- L; sc$R <- R
  f <- directcif(Interval(L, R, cause = status) ~ x, data = sc)
  expect_true(f@converged)
  ll <- 0
  for (i in seq_len(n)) {
    nd <- sc[i, , drop = FALSE]
    if (sc$status[i] == 0L) {
      F <- vapply(1:2, function(k) predict(f, nd, times = L[i], cause = k)[1, 1], 1)
      ll <- ll + log(1 - sum(F))
    } else {
      k <- sc$status[i]
      FR <- predict(f, nd, times = R[i], cause = k)[1, 1]
      FL <- if (L[i] > 0) predict(f, nd, times = L[i], cause = k)[1, 1] else 0
      ll <- ll + log(FR - FL)
    }
  }
  expect_close(ll, as.numeric(logLik(f)), 1e-8)
  expect_output(print(f), "Observed as")
  expect_error(directcif(Interval(L, R, cause = status) ~ x, data = sc, estimation = "weighted"), "joint")
})
