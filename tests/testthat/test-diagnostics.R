# residuals, ph_test and dist_plot

sim_cont <- function(n = 300, seed = 21) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.5); z <- rnorm(n)
  t <- rweibull(n, ifelse(x == 1, 0.8, 1.4), 3 * exp(0.3 * z)); cens <- runif(n, 0, 7)
  data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = z)
}

test_that("Cox residuals agree with survival (Breslow, no ties)", {
  skip_if_not_installed("survival")
  d <- sim_cont()
  expect_identical(anyDuplicated(d$time), 0L)
  f <- cox(Event(time, status) ~ x + z, d, ties = "breslow")
  r <- survival::coxph(survival::Surv(time, status) ~ x + z, d, ties = "breslow")
  expect_close(residuals(f), residuals(r, "martingale"), 1e-6)
  expect_close(residuals(f, "deviance"), residuals(r, "deviance"), 1e-6)
  expect_close(residuals(f, "coxsnell"), d$status - residuals(r, "martingale"), 1e-6)
  expect_close(residuals(f, "schoenfeld"), residuals(r, "schoenfeld"), 1e-6)
  expect_close(residuals(f, "scaledsch"), residuals(r, "scaledsch"), 1e-6)
  # martingale residuals sum to zero
  expect_lt(abs(sum(residuals(f))), 1e-8)
})

test_that("martingale residuals with tied times agree with survival (Breslow)", {
  skip_if_not_installed("survival")
  d <- sim_cont(); d$time <- round(d$time, 1) + 0.1
  f <- cox(Event(time, status) ~ x + z, d, ties = "breslow")
  r <- survival::coxph(survival::Surv(time, status) ~ x + z, d, ties = "breslow")
  expect_close(residuals(f), residuals(r, "martingale"), 1e-6)
})

test_that("ph_test() agrees with survival::cox.zph", {
  skip_if_not_installed("survival")
  d <- sim_cont()
  f <- cox(Event(time, status) ~ x + z, d, ties = "breslow")
  r <- survival::coxph(survival::Surv(time, status) ~ x + z, d, ties = "breslow")
  for (tr in c("identity", "log")) {
    ours <- ph_test(f, transform = tr)[[1]]$table
    ref <- survival::cox.zph(r, transform = tr, terms = FALSE)$table
    expect_close(ours$chisq, ref[, "chisq"], 1e-5)
    expect_close(ours$p, ref[, "p"], 1e-5)
  }
  pt <- ph_test(f)
  expect_lt(pt[[1]]$table["x", "p"], 0.05)        # x has non-proportional hazards
  expect_output(print(pt), "Reading: the effect of x")
  with_null_device(expect_silent(plot(pt)))
})

test_that("cause-specific residuals and tests use the cause-specific events", {
  skip_if_not_installed("survival")
  d <- sim_cr(300, 5)
  f <- csc(Event(time, status) ~ x + z, d, ties = "breslow")
  r <- survival::coxph(survival::Surv(time, status == 2) ~ x + z, d, ties = "breslow")
  expect_close(residuals(f)$cause2, residuals(r, "martingale"), 1e-6)
  pt <- ph_test(f, transform = "identity")
  expect_length(pt, 2L)
  expect_close(pt[[2]]$table$chisq, survival::cox.zph(r, transform = "identity", terms = FALSE)$table[, "chisq"], 1e-5)
  expect_error(ph_test(csc(Event(time, status) ~ x, d, baseline = "weibull")), "Cox baseline")
})

test_that("parametric residuals: Cox-Snell is the fitted cumulative hazard", {
  d <- sim_cont()
  f <- parsurv(Event(time, status) ~ x + z, d, baseline = "weibull")
  cs <- residuals(f, "coxsnell")
  sh <- exp(coef(f)[4]); eta <- coef(f)[1] + coef(f)[2] * d$x + coef(f)[3] * d$z
  expect_close(cs, exp(eta) * d$time^sh, 1e-8)
  expect_close(residuals(f), d$status - cs, 1e-10)
  g <- csc(Event(time, status) ~ x, sim_cr(), baseline = "exponential")
  expect_length(residuals(g, "deviance"), 2L)
})

test_that("dist_plot() draws for formulas and fitted models", {
  d <- sim_cont()
  with_null_device({
    expect_s3_class(dist_plot(Event(time, status) ~ x, d), "data.frame")
    expect_s3_class(dist_plot(parsurv(Event(time, status) ~ x, d)), "data.frame")
    expect_s3_class(dist_plot(cox(Event(time, status) ~ x, d)), "data.frame")
    expect_s3_class(dist_plot(csc(Event(time, status) ~ x, sim_cr(), baseline = "weibull"), cause = 2), "data.frame")
  })
  expect_error(dist_plot(1), "formula or a model")
})

test_that("printed models carry Reading and Fit lines", {
  d <- sim_cr(400, 2)
  expect_output(print(cox(Event(time, dead) ~ x + z, d)), "Fit: log partial likelihood")
  expect_output(print(fg(Event(time, status) ~ x + z, d)), "Reading")
  expect_output(print(csc(Event(time, status) ~ x + z, d)), "Reading")
})

test_that("every result names its sources", {
  d <- sim_cr(200, 8)
  f <- cox(Event(time, dead) ~ x, d)
  expect_output(print(f), "Methods: Cox \\(1972\\); Efron \\(1977\\)")
  expect_output(print(cox(Event(time, dead) ~ x, d, ties = "breslow")), "Breslow \\(1974\\)")
  expect_output(print(fg(Event(time, status) ~ x, d)), "Fine & Gray \\(1999\\)")
  expect_output(print(csc(Event(time, status) ~ x, d)), "Prentice et al. \\(1978\\)")
  expect_output(print(km(Event(time, dead) ~ x, d)), "Kaplan & Meier")
  expect_output(print(aj(Event(time, status) ~ 1, d)), "Aalen & Johansen")
  g <- parsurv(Event(time, dead) ~ x, d, baseline = "gengamma")
  expect_output(print(g), "Prentice \\(1974\\)")
  r <- sx_references(f)
  expect_true(any(grepl("Regression models and life-tables", r)))
  expect_output(sx_references("brier_score"), "Graf E")
  expect_error(sx_references("nothing"), "no references")
  expect_gt(length(suppressWarnings(capture.output(all <- sx_references()))), 100)
})
