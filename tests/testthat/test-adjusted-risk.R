# adjusted_risk(): standardisation, inverse probability weighting, bootstrap and output

ar_data <- function(seed = 2, n = 400) {
  set.seed(seed)
  age <- rnorm(n)
  trt <- rbinom(n, 1, plogis(0.8 * age))
  t1 <- rexp(n, 0.2 * exp(-0.5 * trt + 0.7 * age)); t2 <- rexp(n, 0.1); cens <- runif(n, 0, 8)
  data.frame(time = pmin(t1, t2, cens), status = ifelse(cens < pmin(t1, t2), 0L, ifelse(t1 < t2, 1L, 2L)),
             dead = as.integer(pmin(t1, cens) == t1), time1 = pmin(t1, cens),
             trt = factor(trt, 0:1, c("control", "treated")), age = age)
}

test_that("the g-formula is the average of the predictions with the treatment set", {
  d <- ar_data()
  fit <- cox(Event(time1, dead) ~ trt + age, data = d)
  ar <- adjusted_risk(fit, d, "trt", times = c(1, 3), B = 0)
  for (l in c("control", "treated")) {
    nd <- d; nd$trt <- factor(l, levels = levels(d$trt))
    ref <- colMeans(1 - predict(fit, nd, times = c(1, 3), type = "survival"))
    expect_close(ar@risk$risk[ar@risk$level == l], ref, 1e-10)
  }
  cc <- ar@contrast
  expect_close(cc$difference, ar@risk$risk[ar@risk$level == "treated"] - ar@risk$risk[ar@risk$level == "control"], 1e-12)
  expect_true(all(cc$difference < 0))                 # treatment lowers the risk
})

test_that("IPW without covariates is Kaplan-Meier / Aalen-Johansen by group; weighted estimates agree with survival", {
  skip_if_not_installed("survival")
  d <- ar_data(3)
  fit <- cox(Event(time1, dead) ~ trt + age, data = d)
  ar <- adjusted_risk(fit, d, "trt", times = c(1, 3), method = "ipw", ps_formula = trt ~ 1, B = 0)
  for (l in c("control", "treated")) {
    s <- survival::survfit(survival::Surv(time1, dead) ~ 1, data = d[d$trt == l, ])
    expect_close(ar@risk$risk[ar@risk$level == l], 1 - summary(s, times = c(1, 3))$surv, 1e-10)
  }
  set.seed(1); w <- runif(nrow(d), 0.5, 2)
  s <- survival::survfit(survival::Surv(time1, dead) ~ 1, data = d, weights = w)
  expect_close(SurvX:::.ar_wcif(d$time1, d$dead, NULL, w, 1L, c(0.5, 1, 2, 4)), 1 - summary(s, times = c(0.5, 1, 2, 4))$surv, 1e-10)
  s2 <- survival::survfit(survival::Surv(time, factor(status, 0:2)) ~ 1, data = d, weights = w)
  ss <- summary(s2, times = c(0.5, 1, 2, 4))
  for (k in 1:2) expect_close(SurvX:::.ar_wcif(d$time, d$status, NULL, w, k, c(0.5, 1, 2, 4)), ss$pstate[, k + 1L], 1e-10)
})

test_that("the treatment model is logistic regression for two values", {
  d <- ar_data(4)
  X <- model.matrix(~ age, d)
  P <- SurvX:::.ar_mlogit(X, d$trt)
  g <- glm(trt ~ age, binomial, data = d)
  expect_close(P[, "treated"], fitted(g), 1e-6)
  expect_close(rowSums(P), rep(1, nrow(d)), 1e-12)
  d$arm <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  P3 <- SurvX:::.ar_mlogit(model.matrix(~ age, d), d$arm)
  expect_close(colMeans(P3), as.numeric(table(d$arm)) / nrow(d), 1e-6)   # score equations: fitted = observed shares
})

test_that("both methods, competing risks, time_effects fits and balance", {
  d <- ar_data(5)
  fit <- cox(Event(time1, dead) ~ trt + age, data = d)
  ar <- adjusted_risk(fit, d, "trt", times = c(1, 3), method = "both", B = 0)
  expect_setequal(unique(ar@risk$method), c("gformula", "ipw"))
  expect_true(max(abs(ar@risk$risk[ar@risk$method == "gformula"] - ar@risk$risk[ar@risk$method == "ipw"])) < 0.1)
  b <- as.data.frame(ar, what = "balance")
  expect_true(abs(b$after[b$covariate == "age"]) < abs(b$before[b$covariate == "age"]))
  expect_output(print(ar), "Risk difference")
  expect_output(summary(ar), "No unmeasured confounding")
  for (f in list(fg(Event(time, status) ~ trt + age, data = d, cause = 1), csc(Event(time, status) ~ trt + age, data = d))) {
    a2 <- adjusted_risk(f, d, "trt", times = c(1, 3), B = 0)
    expect_identical(a2@kind, "competing")
    expect_true(all(a2@risk$risk >= 0 & a2@risk$risk <= 1))
  }
  te <- time_effects(fit, fix = list(trt = "log"))
  a3 <- adjusted_risk(te, d, "trt", times = c(1, 3), B = 0)
  nd <- d; nd$trt <- factor("treated", levels = levels(d$trt))
  expect_close(a3@risk$risk[a3@risk$level == "treated"], colMeans(1 - predict(te, nd, times = c(1, 3))), 1e-10)
  with_null_device(expect_no_error(plot(ar)))
  html <- tempfile(fileext = ".html")
  expect_message(report(ar, file = html), "report written")
  expect_true(any(grepl("g-formula", readLines(html))))
})

test_that("the bootstrap refits the model, is reproducible and leaves the RNG alone", {
  d <- ar_data(6, 250)
  fit <- cox(Event(time1, dead) ~ trt + age, data = d)
  set.seed(9); before <- runif(1); set.seed(9)
  a <- adjusted_risk(fit, d, "trt", times = c(1, 3), B = 30, seed = 4)
  expect_identical(runif(1), before)
  b <- adjusted_risk(fit, d, "trt", times = c(1, 3), B = 30, seed = 4)
  expect_identical(a@risk, b@risk)
  expect_identical(a@ok, 30L)
  expect_true(all(a@risk$lower <= a@risk$risk & a@risk$risk <= a@risk$upper))
  expect_true(all(is.finite(a@contrast$p)))
})

test_that("adjusted_risk() refuses what it cannot do", {
  d <- ar_data()
  fit <- cox(Event(time1, dead) ~ trt + age, data = d)
  expect_error(adjusted_risk(fit, d, "arm"), "column of data")
  expect_error(adjusted_risk(fit, d[1:100, ], "trt"), "does not match")
  expect_error(adjusted_risk(fit, d, "age"), "more than 10 values")
  expect_error(adjusted_risk(fit, d, "trt", B = 5), "at least 20")
  expect_error(adjusted_risk(km(Event(time1, dead) ~ 1, d), d, "trt"), "takes a model")
  d$other <- factor(sample(c("x", "y"), nrow(d), TRUE))
  expect_error(adjusted_risk(fit, d, "other"), "not a covariate")
})
