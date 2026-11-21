test_that("Event() builds, subsets and prints", {
  y <- Event(c(1, 2, 3, 4, 5), c(0L, 1L, 2L, 1L, 0L))
  expect_s3_class(y, "SurvXEvent")
  expect_s3_class(y[2:4], "SurvXEvent")
  expect_identical(nrow(unclass(y[2:4])), 3L)
  expect_output(print(y))
  expect_false(any(is.na(y)))
})

test_that("Event() in a formula gives the same fit as the data columns", {
  d <- sim_cr()
  f1 <- fg(Event(time, status) ~ x + z, d)
  d2 <- d; names(d2)[1:4] <- c("t", "a", "b", "s")
  f2 <- fg(Event(t, s) ~ a + b, d2)
  expect_close(coef(f1), coef(f2), 1e-12)
})

test_that("times that differ only by rounding are tied, as in survival (timefix)", {
  set.seed(5)
  n <- 200
  t0 <- round(rexp(n), 3)
  t <- t0 * (1 + sample(c(0, 1e-12), n, replace = TRUE))      # rounding-level differences
  d <- data.frame(time = t, status = rbinom(n, 1, 0.7), x = rnorm(n))
  y <- unclass(Event(d$time, d$status))
  expect_identical(length(unique(y[, "time"])), length(unique(t0)))
  y0 <- unclass(Event(d$time, d$status, timefix = FALSE))
  expect_identical(y0[, "time"], d$time)
  # survival merges the same values
  k <- km(Event(time, status) ~ 1, d)
  s <- survival::survfit(survival::Surv(time, status) ~ 1, data = d)
  expect_identical(sum(as.data.frame(k)$n.event > 0), sum(s$n.event > 0))
  f <- cox(Event(time, status) ~ x, d)
  r <- survival::coxph(survival::Surv(time, status) ~ x, data = d)
  expect_close(coef(f), coef(r), 1e-8)
  expect_close(residuals(f), residuals(r), 1e-8)
})

test_that("timefix keeps entry < time", {
  expect_error(Event(c(1, 2), c(1, 1), entry = c(1 - 1e-15, 0)), "rounding error|smaller than time")
})
