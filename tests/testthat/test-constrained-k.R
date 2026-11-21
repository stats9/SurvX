# Constrained model for K >= 3 causes. Reference values: the C++ fit checked against an
# independent log-likelihood (paper2/proto/check_ck2.py of the development repository; the
# log-likelihoods agree to 1e-12) and, for K = 2, identical to the two-cause model.
ck_data <- function() {
  i <- 1:300
  data.frame(time = ((i * 37) %% 113) / 10 + 0.05,
             status = c(0L, 1L, 2L, 3L, 1L, 2L)[(i * 7) %% 6 + 1L],
             x = ((i * 5) %% 11) / 10 - 0.5)
}

test_that("constrained fit with three causes matches the reference", {
  d <- ck_data()
  f <- directcif(Event(time, status) ~ x, d, estimation = "constrained", cause = 1, cond = "weibull")
  expect_true(f@converged)
  expect_equal(names(coef(f)), c("1:rho", "1:log_tau", "1:x", "2:log_lambda", "2:log_kappa",
                                 "3:log_lambda", "3:log_kappa", "3:share_(Intercept)", "3:share_x"))
  expect_equal(unname(coef(f)), c(-0.1187750884, -2.8159744639, 0.0125199365, -1.8687368833, 0.4173793309,
                                  -2.0444170452, 0.7859263298, -0.6484394133, -0.185202035), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))), c(0.0143490838, 0.1319075531, 0.3078145888, 0.0736045287, 0.085121726,
                                              0.0721656529, 0.1198299087, 0.1762973872, 0.5443678209), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), -971.5675904875836, tolerance = 1e-8)
})

test_that("the long-term incidences of the K causes add up to one for every covariate value", {
  d <- ck_data()
  f <- directcif(Event(time, status) ~ x, d, estimation = "constrained", cause = 1, cond = "weibull")
  nd <- data.frame(x = seq(-3, 3, by = 1))
  tot <- predict(f, nd, times = Inf, cause = 1) + predict(f, nd, times = Inf, cause = 2) + predict(f, nd, times = Inf, cause = 3)
  expect_equal(as.numeric(tot), rep(1, nrow(nd)), tolerance = 1e-10)
  expect_equal(as.numeric(predict(f, data.frame(x = 0), times = Inf, cause = 2)), 0.39674803451441426, tolerance = 1e-6)
})

test_that("one distribution per competing cause, and K = 2 is the two-cause model", {
  d <- ck_data()
  f2 <- directcif(Event(time, status) ~ x, d, estimation = "constrained", cause = 1, cond = c("weibull", "lognormal"))
  expect_equal(names(coef(f2))[6:7], c("3:mu", "3:log_sigma"))
  expect_error(directcif(Event(time, status) ~ x, d, estimation = "constrained", cond = c("weibull", "weibull", "weibull")),
               "length 1 or 2")
  d2 <- d; d2$status[d2$status == 3L] <- 0L
  a <- directcif(Event(time, status) ~ x, d2, estimation = "constrained", cond = "weibull")
  expect_equal(length(coef(a)), 5L)
  expect_output(print(f2), "share")
})

test_that("cif_data works for the K-cause model", {
  d <- ck_data()
  f <- directcif(Event(time, status) ~ x, d, estimation = "constrained", cause = 1, cond = "weibull")
  cd <- cif_data(f, cause = 1:3, times = c(1, 5))
  expect_true(all(is.finite(cd$cif)))
})
