# Independent R code of the weighted pseudo-likelihood (Gompertz baseline) and of the constrained
# likelihood of Shi, Cheng and Jeong (2013) with a Gompertz PH model for cause 1 and a Weibull G.
ref_weighted <- function(th, time, status, X, alpha, cause = 1) {
  rho <- th[1]; tau <- exp(th[2]); eta <- drop(X %*% th[-(1:2)])
  eu <- function(t, e) exp(e) * tau * expm1(rho * t) / rho
  lS <- function(t, e) if (alpha == 0) -eu(t, e) else -log1p(alpha * eu(t, e)) / alpha
  lf <- function(t, e) e + log(tau) + rho * t - (if (alpha == 0) eu(t, e) else (1 / alpha + 1) * log1p(alpha * eu(t, e)))
  ct <- sort(unique(time[status == 0])); G <- 1; mass <- numeric(length(ct))
  for (j in seq_along(ct)) {
    dc <- sum(time == ct[j] & status == 0); r <- sum(time > ct[j]) + dc
    Gn <- G * (1 - dc / r); mass[j] <- G - Gn; G <- Gn
  }
  if (G > 1e-14) { ct <- c(ct, max(time)); mass <- c(mass, G) }
  ll <- 0
  for (i in seq_along(time)) {
    if (status[i] == cause) ll <- ll + lf(time[i], eta[i])
    else if (status[i] == 0) ll <- ll + lS(time[i], eta[i])
    else { k <- ct >= time[i]; if (any(k)) ll <- ll + sum(mass[k] / sum(mass[k]) * lS(ct[k], eta[i])) }
  }
  ll
}
ref_constrained <- function(th, time, status, X) {
  p <- ncol(X); rho <- th[1]; tau <- exp(th[2]); eta <- drop(X %*% th[2 + seq_len(p)])
  lam <- exp(th[p + 3]); ka <- exp(th[p + 4])
  Ft <- -expm1(-exp(eta) * tau * expm1(rho * time) / rho); Finf <- -expm1(exp(eta) * tau / rho)
  lf1 <- eta + log(tau) + rho * time - exp(eta) * tau * expm1(rho * time) / rho
  G <- pweibull(time, ka, 1 / lam); g <- dweibull(time, ka, 1 / lam)
  sum(lf1[status == 1]) + sum(log(1 - Finf[status == 2]) + log(g[status == 2])) + sum(log(1 - Ft - (1 - Finf) * G)[status == 0])
}

test_that("weighted estimation matches the pseudo-likelihood formula and is close to Fine-Gray", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed")
  fm <- Event(time, status) ~ sex + age + hgb + mspike
  X <- model.matrix(~ sex + age + hgb + mspike, d)[, -1]
  for (lk in c("PH", "PO")) {
    fw <- directcif(fm, d, link = lk, estimation = "weighted", cause = 1)
    expect_true(fw@converged)
    expect_close(as.numeric(logLik(fw)), ref_weighted(coef(fw), d$time, d$status, X, if (lk == "PH") 0 else 1), 1e-9)
  }
  fw <- directcif(fm, d, estimation = "weighted"); ff <- fg(fm, d)
  expect_close(coef(fw)[-(1:2)], coef(ff), 0.1)
  expect_close(sqrt(diag(vcov(fw)))[-(1:2)], sqrt(diag(vcov(ff))), 0.1)
  pw <- predict(fw, d[1:5, ], times = c(50, 100, Inf), se = TRUE)
  expect_true(all(pw$cif >= 0 & pw$cif <= 1 & is.finite(pw$se)))
})

test_that("constrained model matches Shi et al. and its plateaus add up to one", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed")
  fm <- Event(time, status) ~ sex + age + hgb + mspike
  X <- model.matrix(~ sex + age + hgb + mspike, d)[, -1]
  fc <- directcif(fm, d, estimation = "constrained", cause = 1, cond = "weibull")
  expect_close(as.numeric(logLik(fc)), ref_constrained(coef(fc), d$time, d$status, X), 1e-9)
  p1 <- predict(fc, d[1:5, ], times = c(50, Inf), cause = 1)
  p2 <- predict(fc, d[1:5, ], times = c(50, Inf), cause = 2, se = TRUE)
  expect_close(p1[, 2] + p2$cif[, 2], rep(1, 5), 1e-12)
  expect_true(all(p1[, 1] + p2$cif[, 1] < 1) && all(is.finite(p2$se)))
  for (cd in c("logistic", "loglogistic", "lognormal")) {
    f <- directcif(fm, d, link = "PO", estimation = "constrained", cause = 1, cond = cd)
    expect_true(is.finite(as.numeric(logLik(f))), info = cd)
  }
})

test_that("weighted and constrained fits run on simulated data", {
  d <- sim_cr()
  f6 <- directcif(Event(time, status) ~ x + z, d, estimation = "weighted")
  f7 <- directcif(Event(time, status) ~ x + z, d, estimation = "constrained", cond = "weibull")
  expect_output(print(summary(f6))); expect_output(print(summary(f7)))
})
