# parametric cure models

sim_cure <- function(n = 500, seed = 61, entry = FALSE) {
  set.seed(seed)
  x <- rnorm(n); z <- rbinom(n, 1, 0.5)
  pi <- plogis(-0.5 + 1.0 * z)
  cured <- runif(n) < pi
  t <- ifelse(cured, Inf, rweibull(n, 1.3, 2 * exp(-0.4 * x)))
  cens <- runif(n, 0, 8)
  d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = z)
  if (entry) { d$entry <- runif(n, 0, 0.3) * d$time }
  d
}

# independent R code of the log-likelihood (Weibull PH latency)
ref_cure_ll <- function(th, d, mixture) {
  u <- th[1] + th[2] * d$x; k <- exp(th[3]); c <- th[4] + th[5] * d$z; pi <- plogis(c)
  H <- function(t) exp(u) * t^k
  lf <- u + log(k) + (k - 1) * log(d$time) - H(d$time)
  Su <- exp(-H(d$time)); Se <- if (is.null(d$entry)) 1 else exp(-H(d$entry))
  if (mixture) {
    ll <- ifelse(d$status == 1, log(1 - pi) + lf, log(pi + (1 - pi) * Su)) - log(pi + (1 - pi) * Se)
  } else {
    th_c <- -log(pi)
    ll <- ifelse(d$status == 1, log(th_c) + lf, 0) - th_c * (1 - Su) + th_c * (1 - Se)
  }
  sum(ll)
}

test_that("cure() maximises the mixture and non-mixture likelihoods", {
  for (en in c(FALSE, TRUE)) for (ty in c("mixture", "nonmixture")) {
    d <- sim_cure(entry = en)
    f <- if (en) cure(Event(time, status, entry = entry) ~ x, d, cure = ~ z, type = ty)
         else cure(Event(time, status) ~ x, d, cure = ~ z, type = ty)
    expect_true(f@converged)
    expect_close(as.numeric(logLik(f)), ref_cure_ll(coef(f), d, ty == "mixture"), 1e-10)
    opt <- stats::optim(coef(f) + 0.1, function(th) -ref_cure_ll(th, d, ty == "mixture"), method = "BFGS",
                        control = list(reltol = 1e-14, maxit = 1000))
    expect_lte(as.numeric(logLik(f)) + opt$value, 1e-6)          # no better optimum
    expect_close(coef(f), opt$par, 1e-3)
  }
})

test_that("cure() agrees with flexsurvcure", {
  skip_if_not_installed("flexsurvcure")
  d <- sim_cure()
  for (mix in c(TRUE, FALSE)) {
    f <- cure(Event(time, status) ~ x, d, cure = ~ z, type = if (mix) "mixture" else "nonmixture")
    r <- flexsurvcure::flexsurvcure(survival::Surv(time, status) ~ z, data = d, dist = "weibullPH", mixture = mix,
                                    anc = list(scale = ~ x))
    expect_close(as.numeric(logLik(f)), r$loglik, 1e-5)
    nd <- data.frame(x = c(0, 1), z = c(0, 1))
    tt <- c(0.5, 2, 6)
    for (i in 1:2) {
      ref <- summary(r, newdata = nd[i, , drop = FALSE], t = tt, type = "survival", tidy = TRUE)$est
      expect_close(predict(f, nd[i, , drop = FALSE], times = tt), ref, 1e-3)
    }
  }
})

test_that("cure(): predictions, output and generics", {
  d <- sim_cure()
  f <- cure(Event(time, status) ~ x, d, cure = ~ z)
  pc <- predict(f, data.frame(x = 0, z = 0:1), type = "cured", se = TRUE)
  expect_close(pc$est, plogis(coef(f)[4] + coef(f)[5] * 0:1), 1e-12)
  expect_true(all(pc$lower < pc$est & pc$est < pc$upper))
  s <- predict(f, data.frame(x = 0, z = 1), times = c(1, 50))
  expect_close(s[1, 2], pc$est[2], 1e-6)                          # the survival levels off at the cured fraction
  expect_output(print(f), "Cured fraction")
  expect_output(print(f), "odds ratio of being cured")
  expect_output(print(f), "Farewell \\(1982\\)")
  expect_output(print(summary(f)), "Latency distribution")
  expect_length(residuals(f, "coxsnell"), nrow(d))
  with_null_device(expect_no_error(plot(f)))
  e <- estimates(f)
  expect_setequal(unique(e$role), c("baseline", "beta", "shape", "cure_baseline", "cure"))
  expect_s3_class(update(f, type = "nonmixture"), "SurvX::CureFit")
  # default: the cured fraction uses the latency covariates; ~ 1: a common fraction
  expect_length(coef(cure(Event(time, status) ~ x, d)), 5L)
  expect_length(coef(cure(Event(time, status) ~ x, d, cure = ~ 1)), 4L)
  # other latency distributions
  for (b in c("lognormal", "loglogistic", "exponential"))
    expect_true(cure(Event(time, status) ~ x, d, cure = ~ z, baseline = b)@converged, info = b)
})
