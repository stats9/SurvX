# Time-varying coefficients: cox(... + tt(x), tt = g) against survival::coxph(tt = function(x, t, ...) x * g(t))

tt_data <- function(seed = 4, n = 300) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.5); z <- rnorm(n)
  t <- rweibull(n, ifelse(x == 1, 0.7, 1.5), 3 * exp(-0.3 * z)); cens <- runif(n, 0, 6)
  data.frame(time = round(pmin(t, cens), 2) + 0.01, status = as.integer(t <= cens), x = x, z = z,
             g = sample(c("a", "b"), n, TRUE))
}

test_that("tt() agrees with coxph(tt = ) for Efron and Breslow ties", {
  skip_if_not_installed("survival")
  d <- tt_data()
  for (ties in c("efron", "breslow")) {
    f <- cox(Event(time, status) ~ x + z + tt(x), d, tt = function(t) log(t), ties = ties)
    r <- survival::coxph(survival::Surv(time, status) ~ x + z + tt(x), data = d, ties = ties,
                         tt = function(x, t, ...) x * log(t))
    expect_close(coef(f), coef(r), 1e-7)
    expect_close(vcov(f), vcov(r), 1e-7)
    expect_close(as.numeric(logLik(f)), r$loglik[2], 1e-7)
  }
})

test_that("tt() works with strata and with another function of time", {
  skip_if_not_installed("survival")
  d <- tt_data(5)
  f <- cox(Event(time, status) ~ x + tt(x) + strata(g), d, tt = function(t) t)
  r <- survival::coxph(survival::Surv(time, status) ~ x + tt(x) + survival::strata(g), data = d,
                       tt = function(x, t, ...) x * t)
  expect_close(coef(f), coef(r), 1e-7)
  expect_close(vcov(f), vcov(r), 1e-7)
})

test_that("tt_effect() gives exp(b + gamma g(t)) and methods refuse tt fits", {
  d <- tt_data(6)
  f <- cox(Event(time, status) ~ x + tt(x), d, tt = function(t) log(t))
  te <- tt_effect(f, times = c(1, 2))
  b <- coef(f)
  expect_close(te$hr, exp(b[["x"]] + b[["tt(x)"]] * log(c(1, 2))), 1e-12)
  expect_true(all(te$lower < te$hr & te$hr < te$upper))
  expect_error(predict(f, d[1:2, ], times = 1, type = "survival"), "time-varying")
  expect_error(ph_test(f), "time-varying")
  expect_error(cox(Event(time, status) ~ x + tt(x), d), "tt = function")
  expect_output(print(f), "Time-varying effects")
})
