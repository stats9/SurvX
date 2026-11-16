test_that("cox() agrees with survival::coxph (Efron and Breslow)", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed")
  f <- cox(Event(time, dead) ~ sex + age + hgb + mspike, d)
  r <- survival::coxph(survival::Surv(time, dead) ~ sex + age + hgb + mspike, d)
  expect_close(coef(f), coef(r))
  expect_close(sqrt(diag(vcov(f))), sqrt(diag(vcov(r))))
  expect_close(as.numeric(logLik(f)), r$loglik[2])
  nd <- d[1:5, ]
  s <- summary(survival::survfit(r, newdata = nd), times = c(50, 100))$surv
  expect_close(predict(f, nd, times = c(50, 100), type = "survival"), t(s))
  fb <- cox(Event(time, dead) ~ age, d, ties = "breslow")
  rb <- survival::coxph(survival::Surv(time, dead) ~ age, d, ties = "breslow")
  expect_close(coef(fb), coef(rb))
})

test_that("fg() agrees with cmprsk::crr", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed"); skip_if_not_installed("cmprsk")
  f <- fg(Event(time, status) ~ sex + age + hgb + mspike, d)
  X <- model.matrix(~ sex + age + hgb + mspike, d)[, -1]
  r <- cmprsk::crr(d$time, d$status, X, failcode = 1, cencode = 0, gtol = 1e-10)
  expect_close(coef(f), r$coef)
  expect_close(sqrt(diag(vcov(f))), sqrt(diag(r$var)), 1e-5)
  pr <- predict(r, X[1:5, ]); tp <- c(50, 100)
  refp <- sapply(tp, function(s) pr[max(which(pr[, 1] <= s)), -1])
  expect_close(predict(f, d[1:5, ], times = tp), refp)
})

test_that("csc() agrees with coxph and riskRegression::CSC", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed")
  f <- csc(Event(time, status) ~ sex + age + hgb + mspike, d)
  r1 <- survival::coxph(survival::Surv(time, status == 1) ~ sex + age + hgb + mspike, d)
  expect_close(coef(f)[[1]], coef(r1))
  skip_if_not_installed("riskRegression"); skip_if_not_installed("prodlim")
  Hist <- prodlim::Hist
  rr <- riskRegression::CSC(Hist(time, status) ~ sex + age + hgb + mspike, data = d)
  expect_close(predict(f, d[1:5, ], times = c(50, 100), cause = 1),
               riskRegression::predictRisk(rr, newdata = d[1:5, ], times = c(50, 100), cause = 1))
})

test_that("regression fits give predictions in [0, 1] and working methods", {
  d <- sim_cr(); nd <- data.frame(x = c(0, 1), z = 0)
  f1 <- cox(Event(time, dead) ~ x + z, d); f2 <- fg(Event(time, status) ~ x + z, d)
  f3 <- csc(Event(time, status) ~ x + z, d)
  for (p in list(predict(f2, nd, times = c(1, 2)), predict(f3, nd, times = c(1, 2)))) expect_true(all(p >= 0 & p <= 1))
  for (f in list(f1, f2, f3)) { expect_output(print(f)); expect_output(print(summary(f))) }
  expect_s3_class(estimates(f1), "data.frame")
  expect_close(estimates(f2)$exp_estimate, exp(coef(f2)), 1e-12)
  expect_type(confint(f3), "list")
})

test_that("character and factor covariates are coded as indicators", {
  d <- sim_cr(); d$grp <- c("a", "b", "c")[1 + (seq_len(nrow(d)) %% 3)]
  f1 <- fg(Event(time, status) ~ grp + z, d)
  expect_length(coef(f1), 3L)
  f2 <- fg(Event(time, status) ~ factor(x) + z, d); d$xf <- factor(d$x); f3 <- fg(Event(time, status) ~ xf + z, d)
  expect_close(unname(coef(f2)), unname(coef(f3)), 1e-12)
  pr <- predict(f1, data.frame(grp = c("c", "a"), z = 0), times = 1)
  expect_true(all(is.finite(pr)) && nrow(pr) == 2)
})
