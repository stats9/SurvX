# Every generic function a user may call on a SurvX object must dispatch to a SurvX method on every
# platform and R version (including R-devel). One test per class; each generic is called and the
# type of its result is checked.

d_cr <- sim_cr(250, 31)

check_model_generics <- function(f, multi = FALSE) {
  expect_output(print(f))
  s <- summary(f)
  expect_output(print(s))
  expect_s3_class(estimates(f), "data.frame")
  if (multi) {
    expect_type(coef(f), "list"); expect_type(vcov(f), "list"); expect_type(confint(f), "list")
  } else {
    expect_true(is.numeric(coef(f)) && length(coef(f)) >= 1L)
    expect_true(is.matrix(vcov(f)))
    expect_true(is.matrix(confint(f)))
  }
  expect_s3_class(logLik(f), "logLik")
  expect_true(is.finite(AIC(f)) || is.na(AIC(f)))
  expect_true(is.finite(BIC(f)) || is.na(BIC(f)))
  expect_true(is.numeric(nobs(f)))
  expect_s3_class(formula(f), "formula")
  expect_true(is.call(getCall(f)))
  expect_s3_class(terms(f), "terms")
  expect_true(is.matrix(model.matrix(f)))
  g <- update(f)
  expect_identical(class(g), class(f))
  invisible(TRUE)
}

test_that("cox(): all generics", {
  f <- cox(Event(time, dead) ~ x + z, d_cr)
  check_model_generics(f)
  expect_true(is.numeric(predict(f, d_cr[1:3, ])))
  expect_length(residuals(f), nrow(d_cr))
  with_null_device(plot(f, newdata = data.frame(x = 0:1, z = 0)))
})

test_that("fg(): all generics", {
  f <- fg(Event(time, status) ~ x + z, d_cr)
  check_model_generics(f)
  expect_true(all(is.finite(predict(f, d_cr[1:3, ], times = 1))))
  with_null_device(plot(f, newdata = data.frame(x = 0:1, z = 0)))
})

test_that("csc() with a Cox baseline: all generics", {
  f <- csc(Event(time, status) ~ x + z, d_cr)
  check_model_generics(f, multi = TRUE)
  expect_true(all(is.finite(predict(f, d_cr[1:3, ], times = 1))))
  expect_type(residuals(f), "list")
  with_null_device(plot(f, newdata = data.frame(x = 0:1, z = 0)))
})

test_that("directcif(): all generics", {
  f <- directcif(Event(time, status) ~ x + z, d_cr)
  check_model_generics(f)
  expect_true(all(is.finite(predict(f, d_cr[1:3, ], times = 1))))
  with_null_device(plot(f, newdata = data.frame(x = 0:1, z = 0)))
})

test_that("parsurv(): all generics", {
  f <- parsurv(Event(time, dead) ~ x + z, d_cr)
  check_model_generics(f)
  expect_true(all(is.finite(predict(f, d_cr[1:3, ], times = 1))))
  expect_length(residuals(f), nrow(d_cr))
  with_null_device(plot(f))
})

test_that("csc() with a parametric baseline: all generics", {
  f <- csc(Event(time, status) ~ x + z, d_cr, baseline = "weibull")
  check_model_generics(f, multi = TRUE)
  expect_true(all(is.finite(predict(f, d_cr[1:3, ], times = 1))))
  expect_type(residuals(f), "list")
  with_null_device(plot(f, newdata = data.frame(x = 0:1, z = 0)))
})

test_that("km(), aj(), hazard_np(): print, summary, plot, as.data.frame", {
  k <- km(Event(time, dead) ~ x, d_cr)
  expect_output(print(k)); expect_no_error(summary(k)); expect_s3_class(as.data.frame(k), "data.frame")
  with_null_device(plot(k))
  a <- aj(Event(time, status) ~ x, d_cr)
  expect_output(print(a)); expect_no_error(summary(a)); expect_s3_class(as.data.frame(a), "data.frame")
  with_null_device(plot(a))
  h <- hazard_np(Event(time, status) ~ 1, d_cr)
  expect_output(print(h))
  with_null_device(plot(h))
})

test_that("compare_models() and ph_test(): print, summary/as.data.frame, plot", {
  cm <- compare_models(Event(time, status) ~ x + z, d_cr, times = 1, folds = 0)
  expect_output(print(cm)); expect_s3_class(summary(cm), "data.frame"); expect_s3_class(as.data.frame(cm), "data.frame")
  with_null_device(plot(cm))
  pt <- ph_test(cox(Event(time, dead) ~ x + z, d_cr))
  expect_output(print(pt)); expect_s3_class(as.data.frame(pt), "data.frame")
  with_null_device(plot(pt))
})

test_that("Event(): S3 methods", {
  y <- Event(c(2, 3, NA), c(1, 0, 1))
  expect_output(print(y))
  expect_s3_class(y[1:2], "Event")
  expect_identical(is.na(y), c(FALSE, FALSE, TRUE))
})

set.seed(62)
dc <- data.frame(x = rnorm(300), z = rbinom(300, 1, 0.5))
dc$tt <- ifelse(runif(300) < plogis(-0.3 + dc$z), Inf, rweibull(300, 1.2, 2)); dc$cc <- runif(300, 0, 8)
dc$time <- pmin(dc$tt, dc$cc); dc$status <- as.integer(dc$tt <= dc$cc)

test_that("cure(): all generics", {
  f <- cure(Event(time, status) ~ x, dc, cure = ~ z)
  check_model_generics(f)
  expect_true(all(is.finite(predict(f, dc[1:3, ], times = 1))))
  expect_length(residuals(f), nrow(dc))
  with_null_device(plot(f))
})
