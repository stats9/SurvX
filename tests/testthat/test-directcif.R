# Independent R code of the joint Gompertz likelihood (Jeong and Fine 2007, eq. 3.9, 4.1-4.2), two causes.
ref_joint_ll <- function(th, time, status, X, alpha) {
  p <- ncol(X); S <- rep(1, length(time)); ll <- 0
  for (k in 1:2) {
    q <- th[(k - 1) * (2 + p) + 1:(2 + p)]; rho <- q[1]; tau <- exp(q[2]); eta <- drop(X %*% q[-(1:2)])
    u <- tau * expm1(rho * time) / rho; du <- tau * exp(rho * time); eu <- exp(eta) * u
    if (alpha == 0) { Fk <- -expm1(-eu); lf <- eta + log(du) - eu }
    else { Fk <- 1 - (1 + alpha * eu)^(-1 / alpha); lf <- eta + log(du) - (1 / alpha + 1) * log1p(alpha * eu) }
    ll <- ll + sum(lf[status == k]); S <- S - Fk
  }
  ll + sum(log(S[status == 0]))
}

test_that("joint Gompertz PH and PO: likelihood, AIC and CIF match the published formulas", {
  d <- sim_cr(400)
  X <- model.matrix(~ x + z, d)[, -1]
  tp <- c(1, 2)
  for (lk in c("PH", "PO")) {
    f <- directcif(Event(time, status) ~ x + z, d, baseline = "gompertz", link = lk)
    expect_true(f@converged)
    expect_close(as.numeric(logLik(f)), ref_joint_ll(coef(f), d$time, d$status, X, if (lk == "PH") 0 else 1), 1e-9)
    expect_close(AIC(f), -2 * as.numeric(logLik(f)) + 2 * length(coef(f)), 1e-12)
    pr <- predict(f, d[1:5, ], times = tp, cause = 1, se = TRUE)
    q <- coef(f)[1:4]
    u <- exp(q[2]) * expm1(q[1] * tp) / q[1]
    eu <- exp(drop(X[1:5, ] %*% q[3:4])) %o% u
    expect_close(pr$cif, if (lk == "PH") -expm1(-eu) else 1 - 1 / (1 + eu), 1e-9)
    expect_true(all(is.finite(pr$se) & pr$se > 0))
  }
})

test_that("joint likelihood agrees with the formulas on mgus2", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed")
  X <- model.matrix(~ sex + age + hgb + mspike, d)[, -1]
  tp <- c(50, 100, Inf)
  for (lk in c("PH", "PO")) {
    f <- directcif(Event(time, status) ~ sex + age + hgb + mspike, d, link = lk)
    expect_close(as.numeric(logLik(f)), ref_joint_ll(coef(f), d$time, d$status, X, if (lk == "PH") 0 else 1), 1e-9)
    # the plateau F(Inf) of eq. (3.10); the Gompertz shape rho is negative in these data
    q <- coef(f)[1:6]
    u <- exp(q[2]) * ifelse(is.finite(tp), expm1(q[1] * tp) / q[1], -1 / q[1])
    eu <- exp(drop(X[1:5, ] %*% q[3:6])) %o% u
    expect_close(predict(f, d[1:5, ], times = tp, cause = 1), if (lk == "PH") -expm1(-eu) else 1 - 1 / (1 + eu), 1e-9)
  }
})

test_that("GOR link, other baselines and the model without covariates run", {
  d <- sim_cr(400)
  fg0 <- directcif(Event(time, status) ~ x + z, d, link = "GOR")
  expect_true(is.finite(AIC(fg0)))
  for (b in c("gompertz3", "weibull", "mix_weibull", "logistic3", "spline")) {
    f <- directcif(Event(time, status) ~ x + z, d, baseline = b)
    pr <- predict(f, d[1:5, ], times = c(1, 2))
    expect_true(all(pr >= 0 & pr <= 1), info = b)
  }
  f0 <- directcif(Event(time, status) ~ 1, d)
  expect_length(coef(f0), 4L)
})

test_that("methods of a direct fit work", {
  d <- sim_cr()
  f <- directcif(Event(time, status) ~ x + z, d)
  expect_length(coef(f), 8L)
  expect_output(print(f)); expect_output(print(summary(f)))
  e <- estimates(f); ci <- confint(f)
  expect_identical(nrow(ci), length(coef(f)))
  expect_close(c(e$lower, e$upper), c(ci[, 1], ci[, 2]), 1e-12)
  expect_close(as.data.frame(summary(f))$estimate, e$estimate, 0)
  fpo <- update(f, link = "PO")
  expect_identical(fpo@link[1], "PO"); expect_identical(nobs(fpo), nobs(f))
  a <- AIC(f, fpo)
  expect_close(a$AIC, sapply(list(f, fpo), function(m) -2 * as.numeric(logLik(m)) + 2 * length(coef(m))), 1e-9)
  expect_identical(ncol(model.matrix(f)), 2L)
  expect_s3_class(formula(f), "formula")
})

test_that("the coherence diagnostic of a joint fit counts patients with summed CIFs above one", {
  d <- sim_cr()
  f <- directcif(Event(time, status) ~ x + z, d)
  n_ex <- f@diagnostics$n_exceed_range
  expect_true(is.numeric(n_ex) && n_ex >= 0 && n_ex <= nrow(d))
})
