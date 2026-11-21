scan_single_data <- function(n = 240, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rbinom(n, 1, 0.5), z = rnorm(n))
  t <- rweibull(n, 1.4, exp(1 - 0.5 * d$x + 0.3 * d$z)); cens <- runif(n, 0, 5)
  d$time <- pmin(t, cens); d$status <- as.integer(t <= cens)
  d
}

test_that("model_scan() for one event type: table, ranking and methods", {
  d <- scan_single_data()
  sc <- model_scan(Event(time, status) ~ x + z, d, depth = "quick", folds = 3)
  expect_true(S7::S7_inherits(sc, ModelScan))
  tab <- as.data.frame(sc)
  expect_true(all(c("model", "ibs_cv", "ipa_cv", "auc_cv", "calib_cv", "aic", "delta_aic", "rank", "se_diff", "equivalent", "checks") %in% names(tab)))
  expect_true("Cox" %in% tab$model)
  expect_gte(sum(tab$fitted), 5L)
  ok <- tab[is.finite(tab$ibs_cv), ]
  expect_equal(min(ok$rank), 1L)
  expect_true(all(ok$ibs_cv > 0 & ok$ibs_cv < 0.25))
  expect_true(all(is.na(tab$aic[tab$model == "Cox"])))                 # partial likelihood: no AIC
  expect_true(any(is.finite(tab$aic)))
  expect_true(all(ok$equivalent[ok$rank == 1L]))
  expect_output(print(sc), "Model scan: one event type")
  expect_output(summary(sc), "NOTE: confidence intervals")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(plot(sc), NA)
})

test_that("model_scan() restores the random number generator and is reproducible", {
  d <- scan_single_data()
  set.seed(42); u0 <- runif(1)
  set.seed(42); a <- model_scan(Event(time, status) ~ x + z, d, depth = "quick", folds = 3, seed = 7); u1 <- runif(1)
  expect_identical(u0, u1)
  b <- model_scan(Event(time, status) ~ x + z, d, depth = "quick", folds = 3, seed = 7)
  expect_equal(a@table$ibs_cv, b@table$ibs_cv)
})

test_that("model_scan() for competing risks and goal = 'effect'", {
  sc <- model_scan(Event(time, status) ~ education, fertility, cause = 1, goal = "effect", depth = "quick", folds = 3)
  tab <- as.data.frame(sc)
  expect_identical(sc@type, "competing")
  expect_true(all(c("cause-specific Cox", "Fine-Gray") %in% tab$model))
  expect_true(all(is.na(tab$aic[tab$model %in% c("cause-specific Cox", "Fine-Gray")])))
  expect_true(!is.null(sc@effects) && nrow(sc@effects) > 0)
  expect_output(summary(sc), "Effects of the covariates")
  f <- tempfile(fileext = ".html")
  report(sc, file = f)
  expect_true(file.exists(f))
  h <- readLines(f)
  expect_true(any(grepl("data:image/png;base64", h)))
  expect_true(any(grepl("Graf E", h)))
})

test_that("model_scan() checks its input", {
  d <- scan_single_data()
  expect_error(model_scan(Event(time, status) ~ 1, d), "at least one covariate")
  expect_error(model_scan(Event(time, status) ~ x, d, folds = 1), "at least 2")
  expect_error(model_scan(Event(time, status) ~ x, d, cause = 2), "cause must be 1")
})

test_that("the base64 encoder of report() is correct", {
  expect_identical(SurvX:::.b64(charToRaw("Man")), "TWFu")
  expect_identical(SurvX:::.b64(charToRaw("Ma")), "TWE=")
  expect_identical(SurvX:::.b64(charToRaw("M")), "TQ==")
  expect_identical(SurvX:::.b64(charToRaw("hello world")), "aGVsbG8gd29ybGQ=")
})
