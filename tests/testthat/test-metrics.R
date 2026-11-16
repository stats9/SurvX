test_that("brier_score() and auc_t() agree with riskRegression::Score", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed")
  skip_if_not_installed("riskRegression"); skip_if_not_installed("prodlim")
  # riskRegression::Score calls Surv() through prodlim without loading survival, so attach it
  suppressPackageStartupMessages(library(survival))
  Hist <- prodlim::Hist
  f <- directcif(Event(time, status) ~ sex + age + hgb + mspike, d)
  risk <- drop(predict(f, times = 100))
  b <- brier_score(d$time, d$status, risk, 100); a <- auc_t(d$time, d$status, risk, 100)
  sc <- riskRegression::Score(list(m = matrix(risk)), formula = Hist(time, status) ~ 1, data = d, times = 100,
                              cause = 1, metrics = c("auc", "brier"), summary = "ipa", cens.model = "km", conf.int = FALSE)
  bs <- as.data.frame(sc$Brier$score); au <- as.data.frame(sc$AUC$score)
  expect_close(b[["brier"]], bs$Brier[bs$model == "m"])
  expect_close(a, au$AUC[au$model == "m"])
})

test_that("calibration() returns one row per risk group", {
  d <- sim_cr(400)
  f <- fg(Event(time, status) ~ x + z, d)
  cal <- calibration(d$time, d$status, drop(predict(f, times = 1)), 1)
  expect_identical(nrow(cal), 10L)
})

test_that("compare_models() and cv_compare() score every model and agree with brier_score()", {
  skip_on_cran()
  d <- sim_cr(400); fm <- Event(time, status) ~ x + z
  cm <- compare_models(fm, d, times = c(1, 2), folds = 3)
  expect_output(print(cm))
  expect_identical(nrow(cm@table), 2L * 7L)
  f <- directcif(fm, d)
  expect_close(cm@table$brier[cm@table$model == "direct gompertz PH" & cm@table$time == 2],
               brier_score(d$time, d$status, drop(predict(f, times = 2)), 2)[["brier"]], 1e-9)
  ff <- fg(fm, d)
  cc <- cv_compare(fm, d, times = c(1, 2), folds = 3, models = list(
    "Fine-Gray" = function(f, d) fg(f, d), "cause-specific Cox" = function(f, d) csc(f, d)))
  expect_identical(nrow(cc@table), 4L)
  expect_close(cc@table$brier[cc@table$model == "Fine-Gray" & cc@table$time == 2],
               brier_score(d$time, d$status, drop(predict(ff, times = 2)), 2)[["brier"]], 1e-9)
  expect_true(all(is.finite(cc@table$brier_cv)))
  with_null_device(plot(cm))
})
