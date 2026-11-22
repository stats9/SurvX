# covariates on very large or very small scales, and collinear covariates

sc_data <- function(n = 400, seed = 21) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  t1 <- rexp(n, 0.3 * exp(0.5 * x1 - 0.4 * x2)); t2 <- rexp(n, 0.2); cen <- runif(n, 0, 6)
  time <- pmin(t1, t2, cen); status <- ifelse(time == cen, 0L, ifelse(time == t1, 1L, 2L))
  data.frame(time = time, status = status, dead = as.integer(status > 0), x1 = x1, x2 = x2)
}

test_that("fits do not depend on the scale of a covariate", {
  d <- sc_data()
  for (k in c(1e6, 1e-6, 37)) {
    dk <- transform(d, x1 = x1 * k)
    a <- cox(Event(time, dead) ~ x1 + x2, d); b <- cox(Event(time, dead) ~ x1 + x2, dk)
    expect_true(b@converged)
    expect_equal(unname(coef(b)) * c(k, 1), unname(coef(a)), tolerance = 1e-8)
    expect_equal(unname(sqrt(diag(vcov(b)))) * c(k, 1), unname(sqrt(diag(vcov(a)))), tolerance = 1e-8)
    a <- fg(Event(time, status) ~ x1 + x2, d); b <- fg(Event(time, status) ~ x1 + x2, dk)
    expect_equal(unname(coef(b)) * c(k, 1), unname(coef(a)), tolerance = 1e-8)
    for (bl in c("weibull", "lognormal", "gengamma")) {
      a <- parsurv(Event(time, dead) ~ x1 + x2, d, baseline = bl); b <- parsurv(Event(time, dead) ~ x1 + x2, dk, baseline = bl)
      expect_true(b@converged)
      expect_equal(b@loglik, a@loglik, tolerance = 1e-8)
      expect_equal(predict(b, dk[1:8, ], times = c(1, 2)), predict(a, d[1:8, ], times = c(1, 2)), tolerance = 1e-6)
    }
    a <- cure(Event(time, dead) ~ x1, d, cure = ~ x1); b <- cure(Event(time, dead) ~ x1, dk, cure = ~ x1)
    expect_equal(b@loglik, a@loglik, tolerance = 1e-8)
    expect_equal(predict(b, dk[1:8, ], times = c(1, 2)), predict(a, d[1:8, ], times = c(1, 2)), tolerance = 1e-6)
  }
})

test_that("collinear or constant covariates stop with a clear message", {
  d <- sc_data(); d$x3 <- 2 * d$x1; d$one <- 1
  expect_error(cox(Event(time, dead) ~ x1 + x3, d), "linear combination.*x3")
  expect_error(parsurv(Event(time, dead) ~ x1 + x3 + x2, d), "x3")
  expect_error(cox(Event(time, dead) ~ x1 + one, d), "constant.*one")
  expect_error(cure(Event(time, dead) ~ x1, d, cure = ~ x1 + x3), "cured fraction.*x3")
  expect_no_error(cox(Event(time, dead) ~ x1 + x2, d))
})
