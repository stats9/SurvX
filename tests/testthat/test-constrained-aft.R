# Constrained model with covariates of the time distributions G_k (timing = accelerated failure
# time, shape). Reference values: the C++ core fitted to the same deterministic data during
# development, where the likelihood was checked against independent code written from the model
# definition (1e-15), the score and information against numerical differentiation (1e-8) and the
# maximum against a general optimiser. The R log-likelihood below is a second independent check.
aft_data <- function() {
  frac <- function(v) v - floor(v)
  i <- 1:600
  x <- 2 * frac(i * 0.7548776662) - 1; age <- 2 * frac(i * 0.5698402910) - 1
  u1 <- frac(i * 0.6180339887 + 0.1); u2 <- frac(i * 0.4142135624 + 0.3)
  u3 <- frac(i * 0.3247179572 + 0.7); uc <- frac(i * 0.2360679775 + 0.5)
  lp <- 0.5 * x + 0.3 * age
  p1 <- 1 - exp(-exp(lp) * 0.75)
  time <- numeric(600); status <- integer(600)
  one <- u1 < p1                                        # cause 1: Gompertz, rho = -0.8, tau = 0.6
  H <- -log(1 - u2[one] * p1[one]) / exp(lp[one])
  time[one] <- log1p(H * (-0.8) / 0.6) / (-0.8); status[one] <- 1L
  pi3 <- 1 / (1 + exp(-(-0.2 + 0.4 * x)))
  three <- !one & u3 < pi3; two <- !one & !three       # Weibull times, time ratios exp(0.5 age), exp(-0.4 age)
  time[two] <- (-log(1 - u2[two]))^(1 / 1.6) / 1.0 * exp(0.5 * age[two]); status[two] <- 2L
  time[three] <- (-log(1 - u2[three]))^(1 / 1.3) / 0.6 * exp(-0.4 * age[three]); status[three] <- 3L
  cz <- 0.5 + 5.5 * uc
  cen <- cz < time
  time[cen] <- cz[cen]; status[cen] <- 0L
  data.frame(time = time, status = status, x = x, age = age)
}

# independent log-likelihood: Gompertz PH for cause 1, Weibull G_k with time ratios exp(eta_k age)
aft_loglik <- function(b, d) {
  lp <- b[["1:x"]] * d$x + b[["1:age"]] * d$age
  rho <- b[["1:rho"]]; tau <- exp(b[["1:log_tau"]])
  u <- tau * expm1(rho * d$time) / rho
  F1 <- 1 - exp(-exp(lp) * u)
  lf1 <- lp + log(tau) + rho * d$time - exp(lp) * u
  q <- exp(-exp(lp) * tau / (-rho))
  l3 <- b[["3:share_(Intercept)"]] + b[["3:share_x"]] * d$x + b[["3:share_age"]] * d$age
  pi <- cbind(1 - stats::plogis(l3), stats::plogis(l3))
  G <- lg <- matrix(0, nrow(d), 2)
  for (j in 1:2) {
    k <- j + 1L
    lam <- exp(b[[paste0(k, ":log_lambda")]] - b[[paste0(k, ":time_age")]] * d$age)
    kap <- exp(b[[paste0(k, ":log_kappa")]])
    y <- (lam * d$time)^kap
    G[, j] <- -expm1(-y)
    lg[, j] <- log(kap) + log(y) - y - log(d$time)
  }
  s <- d$status
  sum(lf1[s == 1L]) + sum(log(q[s == 2L] * pi[s == 2L, 1] ) + lg[s == 2L, 1]) +
    sum(log(q[s == 3L] * pi[s == 3L, 2]) + lg[s == 3L, 2]) +
    sum(log(1 - F1[s == 0L] - q[s == 0L] * rowSums(pi[s == 0L, , drop = FALSE] * G[s == 0L, , drop = FALSE])))
}

test_that("time ratios of the competing causes: reference values and independent likelihood", {
  d <- aft_data()
  expect_equal(as.integer(table(d$status)), c(77L, 293L, 144L, 86L))
  f0 <- directcif(Event(time, status) ~ x + age, d, estimation = "constrained", cause = 1, cond = "weibull")
  f <- directcif(Event(time, status) ~ x + age, d, estimation = "constrained", cause = 1, cond = "weibull", timing = ~ age)
  expect_true(f@converged)
  expect_equal(names(coef(f)), c("1:rho", "1:log_tau", "1:x", "1:age", "2:log_lambda", "2:log_kappa",
                                 "3:log_lambda", "3:log_kappa", "3:share_(Intercept)", "3:share_x", "3:share_age",
                                 "2:time_age", "3:time_age"))
  expect_equal(unname(coef(f)), c(-0.7599255477, -0.5046967157, 0.5301847312, 0.3029407996, -0.0054082373, 0.5287947972,
                                  -0.5063144567, 0.1387710159, -0.3392000498, 0.265471914, -0.3752679015, 0.5344616795,
                                  -0.3940407642), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))), c(0.0608404674, 0.0802416438, 0.1012949889, 0.0996612039, 0.0545215067,
                                              0.0666640254, 0.115612115, 0.0910987957, 0.1421203217, 0.2415522622,
                                              0.2356877323, 0.0957003054, 0.1850185219), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), -1035.2816072623514, tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f)), aft_loglik(coef(f), d), tolerance = 1e-10)
  expect_equal(as.numeric(logLik(f0)), -1050.7041066543761, tolerance = 1e-8)
  expect_gt(as.numeric(logLik(f)), as.numeric(logLik(f0)))       # nested: never below the model without timing
  e <- estimates(f)
  expect_equal(e$role[e$term == "time_age"], c("timing", "timing"))
  expect_true(all(is.finite(e$exp_estimate[e$role == "timing"])))
  # predictions: values, standard errors and the coherence of the long-run risks
  nd <- data.frame(x = c(0, 0.5), age = c(0.5, -0.5))
  p2 <- predict(f, nd, times = c(1, 3), cause = 2, se = TRUE)
  expect_equal(as.numeric(t(p2$cif)), c(0.11656106259002476, 0.24518884102456237, 0.16368151114959378, 0.20725325719159235),
               tolerance = 1e-6)
  expect_equal(as.numeric(t(p2$se)), c(0.01533644634151418, 0.024747385218861873, 0.022618513005950757, 0.027332764741398514),
               tolerance = 1e-5)
  nd2 <- data.frame(x = seq(-2, 2, 1), age = seq(2, -2, -1))
  tot <- predict(f, nd2, times = Inf, cause = 1) + predict(f, nd2, times = Inf, cause = 2) + predict(f, nd2, times = Inf, cause = 3)
  expect_equal(as.numeric(tot), rep(1, nrow(nd2)), tolerance = 1e-10)
  # the timing does not change the long-run risks of a subject, only when the events occur
  expect_equal(as.numeric(predict(f, nd, times = Inf, cause = 2) + predict(f, nd, times = Inf, cause = 3)),
               as.numeric(1 - predict(f, nd, times = Inf, cause = 1)), tolerance = 1e-12)
})

test_that("timing = ~ 1 is the model without covariates of G, digit for digit", {
  d <- aft_data()
  f0 <- directcif(Event(time, status) ~ x + age, d, estimation = "constrained", cause = 1, cond = "weibull")
  f1 <- directcif(Event(time, status) ~ x + age, d, estimation = "constrained", cause = 1, cond = "weibull", timing = ~ 1)
  expect_identical(coef(f1), coef(f0))
  expect_identical(vcov(f1), vcov(f0))
  expect_null(f1@cond_cov)
})

test_that("covariates per cause, shape effects and two causes", {
  d <- aft_data()
  f <- directcif(Event(time, status) ~ x + age, d, estimation = "constrained", cause = 1, cond = c("weibull", "lognormal"),
                 timing = list("2" = ~ age), shape = list("3" = ~ x))
  expect_true(f@converged)
  expect_equal(names(coef(f))[12:13], c("2:time_age", "3:shape_x"))
  expect_equal(unname(coef(f))[12:13], c(0.5241607414, -0.1320958971), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f)), -1043.4329783538053, tolerance = 1e-8)
  d2 <- d; d2$status[d2$status == 3L] <- 2L
  g <- directcif(Event(time, status) ~ x + age, d2, estimation = "constrained", cause = 1, cond = "logistic",
                 timing = ~ age, shape = ~ x)
  expect_true(g@converged)
  expect_equal(names(coef(g)), c("1:rho", "1:log_tau", "1:x", "1:age", "2:log_b", "2:c", "2:time_age", "2:shape_x"))
  expect_equal(unname(coef(g)), c(-0.7402225335, -0.5196866878, 0.5361495763, 0.300176706, 0.2374909461, 0.224792298,
                                  0.0323719686, 0.0528936673), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(g)), -914.0463947663557, tolerance = 1e-8)
  nd <- data.frame(x = c(-1, 1), age = c(1, -1))
  expect_equal(as.numeric(predict(g, nd, times = Inf, cause = 1) + predict(g, nd, times = Inf, cause = 2)), c(1, 1), tolerance = 1e-10)
  expect_output(print(f), "time ratio")
  expect_output(print(f), "shape of G")
  expect_output(print(summary(g)), "Parameters of G \\(covariates of G at zero\\)")
  cd <- cif_data(f, cause = 1:3, times = c(1, 3))
  expect_true(all(is.finite(cd$cif)))
})

test_that("timing and shape are checked", {
  d <- aft_data()
  expect_error(directcif(Event(time, status) ~ x + age, d, timing = ~ age), "estimation = \"constrained\"")
  expect_error(directcif(Event(time, status) ~ x, d, estimation = "constrained", cond = "weibull", timing = ~ age),
               "must also be a term of the model formula")
  expect_error(directcif(Event(time, status) ~ x + age, d, estimation = "constrained", cond = "weibull",
                         timing = list("4" = ~ age)), "not a competing cause")
  expect_error(directcif(Event(time, status) ~ x + age, d, estimation = "constrained", cond = "weibull",
                         timing = list(~ age)), "named by the competing causes")
  expect_error(directcif(Event(time, status) ~ 1, d, estimation = "constrained", cond = "weibull", timing = ~ age),
               "needs covariates")
  small <- d[1:120, ]
  expect_warning(directcif(Event(time, status) ~ x + age, small, estimation = "constrained", cond = "weibull",
                           timing = ~ x + age, shape = ~ x + age), "may be poorly estimated")
})

test_that("compare_models and model_scan include the constrained model with K causes and timing", {
  d <- aft_data()
  cm <- compare_models(Event(time, status) ~ x + age, d, times = c(1, 2), folds = 2,
                       direct = list(c("gompertz", "PH"), c("gompertz", "PH", "constrained", "weibull")),
                       constrained = list(list(timing = ~ age)))
  nm <- cm@fit$model
  expect_true("direct gompertz PH, constrained (G weibull)" %in% nm)        # three causes: fitted in R
  expect_true("direct gompertz PH, constrained (G weibull), timing ~age" %in% nm)
  ft <- cm@fit[nm == "direct gompertz PH, constrained (G weibull), timing ~age", ]
  expect_equal(ft$loglik, -1035.2816072623514, tolerance = 1e-8)
  expect_true(is.finite(ft$aic))
  expect_true(all(is.finite(cm@table$brier_cv[cm@table$model == ft$model])))
  expect_error(compare_models(Event(time, status) ~ x + age, d, times = 1, constrained = list(list(timming = ~ age))),
               "unknown element")
  cand <- .scan_candidates("competing", 1L, "standard", 100L, FALSE, 3L)
  expect_true("direct gompertz PH, constrained, timing" %in% vapply(cand, `[[`, "", "name"))
  expect_false("direct gompertz PH, constrained, timing" %in% vapply(.scan_candidates("competing", 1L, "quick", 100L, FALSE, 3L), `[[`, "", "name"))
})
