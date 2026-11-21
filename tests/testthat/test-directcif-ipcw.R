# Variance of the weighted estimator with the term for the estimated censoring distribution.
# Reference values: an independent Python implementation of the influence function (the exact
# derivative of the total score with respect to case weights, the Kaplan-Meier estimate of the
# censoring distribution recomputed), used during development.
ipcw_data <- function() {
  i <- 1:120
  data.frame(time = ((i * 37) %% 97) / 10 + 0.1,
             status = c(0L, 1L, 2L, 1L)[(i * 7) %% 4 + 1L],
             x = as.numeric(i %% 3 == 0))
}

test_that("variance = 'ipcw' matches the independent implementation (PH link)", {
  d <- ipcw_data()
  f <- directcif(Event(time, status) ~ x, d, estimation = "weighted", cause = 1)
  expect_equal(unname(coef(f)), c(0.161842213, -3.1746247365, 0.0965235834), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))), c(0.0488571177, 0.2884656942, 0.274513599), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(f@diagnostics$vcov_robust))), c(0.0495987245, 0.2900894336, 0.2746867101), tolerance = 1e-6)
  expect_identical(f@diagnostics$variance, "ipcw")
})

test_that("variance = 'ipcw' matches the independent implementation (PO link)", {
  d <- ipcw_data()
  f <- directcif(Event(time, status) ~ x, d, link = "PO", estimation = "weighted", cause = 1)
  expect_equal(unname(coef(f)), c(0.2541516371, -3.2874393941, 0.1943584652), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))), c(0.0538029708, 0.3307788014, 0.3832716016), tolerance = 1e-6)
})

test_that("the variance argument selects the matrix", {
  d <- ipcw_data()
  fr <- directcif(Event(time, status) ~ x, d, estimation = "weighted", variance = "robust")
  fm <- directcif(Event(time, status) ~ x, d, estimation = "weighted", variance = "model")
  expect_equal(vcov(fr), fr@diagnostics$vcov_robust)
  expect_equal(vcov(fm), fm@vcov_model)
  expect_equal(fr@diagnostics$vcov_ipcw, directcif(Event(time, status) ~ x, d, estimation = "weighted")@vcov)
})

test_that("the censoring term is computed on the exact distribution with a censoring grid", {
  d <- ipcw_data()
  fg5 <- directcif(Event(time, status) ~ x, d, estimation = "weighted", censoring_grid = 5L)
  expect_true(all(is.finite(vcov(fg5))))
  expect_identical(fg5@diagnostics$variance, "ipcw")
})
