# parametric survival regression, parametric cause-specific hazards, non-parametric hazard

sim_surv <- function(n = 400, seed = 11, entry = FALSE) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.5); z <- rnorm(n)
  t <- rweibull(n, shape = 1.4, scale = 3 * exp(-0.5 * x + 0.2 * z)); cens <- runif(n, 0, 8)
  d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = z)
  if (entry) { d$entry <- runif(n, 0, 1); d <- d[d$time > d$entry, ] }
  d
}

# independent R code of the log-likelihood with delayed entry (Klein & Moeschberger 2003, Sec. 3.5)
ref_loglik <- function(th, d, dist, scale, X) {
  eta <- drop(th[1] + X %*% th[1 + seq_len(ncol(X))]); s <- th[length(th)]
  e <- if (is.null(d$entry)) rep(0, nrow(d)) else d$entry
  if (scale == "PH") {
    H0 <- switch(dist, exponential = function(u) u, weibull = function(u) u^exp(s), gompertz = function(u) expm1(s * u) / s)
    lh0 <- switch(dist, exponential = 0, weibull = s + (exp(s) - 1) * log(d$time), gompertz = s * d$time)
    return(sum(d$status * (eta + lh0) - exp(eta) * (H0(d$time) - H0(e))))
  }
  sig <- if (dist == "exponential") 1 else exp(s)
  w <- (log(d$time) - eta) / sig
  lf <- switch(dist, exponential = , weibull = w - exp(w), lognormal = dnorm(w, log = TRUE), loglogistic = dlogis(w, log = TRUE))
  lS <- function(w) switch(dist, exponential = , weibull = -exp(w), lognormal = pnorm(w, lower.tail = FALSE, log.p = TRUE),
                           loglogistic = plogis(w, lower.tail = FALSE, log.p = TRUE))
  ll <- ifelse(d$status == 1, lf - log(sig) - log(d$time), lS(w))
  we <- ifelse(e > 0, (log(pmax(e, 1e-300)) - eta) / sig, -Inf)
  sum(ll) - sum(ifelse(e > 0, lS(we), 0))
}

test_that("parsurv() AFT models agree with survival::survreg", {
  skip_if_not_installed("survival")
  d <- sim_surv()
  for (dist in c("exponential", "weibull", "lognormal", "loglogistic")) {
    f <- parsurv(Event(time, status) ~ x + z, d, baseline = dist, scale = "AFT")
    r <- survival::survreg(survival::Surv(time, status) ~ x + z, d, dist = dist)
    expect_true(f@converged, info = dist)
    expect_close(coef(f)[1:3], coef(r), 1e-5)
    expect_close(as.numeric(logLik(f)), r$loglik[2], 1e-6)
    if (dist != "exponential") {
      expect_close(exp(coef(f)[4]), r$scale, 1e-5)
      expect_close(sqrt(diag(vcov(f))), sqrt(diag(vcov(r))), 1e-4)
    }
  }
})

test_that("parsurv() Weibull PH is the reparametrised Weibull AFT model", {
  d <- sim_surv()
  ph <- parsurv(Event(time, status) ~ x + z, d, baseline = "weibull", scale = "PH")
  aft <- parsurv(Event(time, status) ~ x + z, d, baseline = "weibull", scale = "AFT")
  expect_close(as.numeric(logLik(ph)), as.numeric(logLik(aft)), 1e-8)
  sig <- exp(coef(aft)[4])
  expect_close(coef(ph)[2:3], -coef(aft)[2:3] / sig, 1e-5)
  expect_close(exp(coef(ph)[4]), 1 / sig, 1e-5)
  expect_close(predict(ph, d[1:5, ], times = c(1, 3)), predict(aft, d[1:5, ], times = c(1, 3)), 1e-6)
})

test_that("parsurv() maximises the published likelihoods, also with delayed entry", {
  for (entry in c(FALSE, TRUE)) {
    d <- sim_surv(entry = entry)
    X <- as.matrix(d[, c("x", "z")])
    fm <- if (entry) Event(time, status, entry = entry) ~ x + z else Event(time, status) ~ x + z
    for (cfg in list(c("exponential", "PH"), c("weibull", "PH"), c("gompertz", "PH"), c("weibull", "AFT"),
                     c("lognormal", "AFT"), c("loglogistic", "AFT"))) {
      f <- parsurv(fm, d, baseline = cfg[1], scale = cfg[2])
      th <- coef(f)
      if (cfg[1] == "exponential") th <- c(th, 0)
      ll <- ref_loglik(th, d, cfg[1], cfg[2], X)
      expect_close(as.numeric(logLik(f)), ll, 1e-8)
      neg <- function(p) -ref_loglik(if (cfg[1] == "exponential") c(p, 0) else p, d, cfg[1], cfg[2], X)
      o <- stats::optim(coef(f), neg, method = "BFGS", control = list(reltol = 1e-14, maxit = 500))
      expect_lte(-o$value - as.numeric(logLik(f)), 1e-6)
    }
  }
})

test_that("parsurv() agrees with flexsurv for the Gompertz model and with delayed entry", {
  skip_if_not_installed("flexsurv")
  d <- sim_surv(entry = TRUE)
  f <- parsurv(Event(time, status, entry = entry) ~ x + z, d, baseline = "gompertz")
  r <- flexsurv::flexsurvreg(survival::Surv(entry, time, status) ~ x + z, data = d, dist = "gompertz")
  expect_close(as.numeric(logLik(f)), as.numeric(logLik(r)), 1e-5)
  expect_close(unname(coef(f)[c(4, 1, 2, 3)]), unname(coef(r)), 1e-4)
})

test_that("predictions match the closed forms", {
  d <- sim_surv()
  f <- parsurv(Event(time, status) ~ x, d, baseline = "weibull", scale = "AFT")
  mu <- coef(f)[1] + coef(f)[2] * c(0, 1); sig <- exp(coef(f)[3])
  tt <- c(0.5, 2, 5)
  expect_close(predict(f, data.frame(x = 0:1), times = tt), t(sapply(mu, function(m) pweibull(tt, 1 / sig, exp(m), lower.tail = FALSE))), 1e-10)
  expect_close(predict(f, data.frame(x = 0:1), type = "quantile", p = c(0.25, 0.5)),
               t(sapply(mu, function(m) qweibull(c(0.25, 0.5), 1 / sig, exp(m)))), 1e-8)
  h <- predict(f, data.frame(x = 0), times = tt, type = "hazard")
  expect_close(h, dweibull(tt, 1 / sig, exp(mu[1])) / pweibull(tt, 1 / sig, exp(mu[1]), lower.tail = FALSE), 1e-8)
  r <- predict(f, data.frame(x = 0), times = 3, type = "rmst")
  expect_close(r, integrate(function(u) pweibull(u, 1 / sig, exp(mu[1]), lower.tail = FALSE), 0, 3, rel.tol = 1e-12)$value, 1e-8)
  ps <- predict(f, data.frame(x = 0:1), times = tt, se = TRUE)
  expect_true(all(ps$se > 0 & ps$lower < ps$est & ps$est < ps$upper))
  avg <- predict(f, times = tt, average = TRUE)
  expect_close(avg, colMeans(predict(f, times = tt)), 1e-12)
})

test_that("parametric cause-specific hazards give coherent CIFs equal to the integral of h_k S", {
  d <- sim_cr(400)
  fp <- csc(Event(time, status) ~ x + z, d, baseline = c("weibull", "lognormal"))
  expect_true(S7::S7_inherits(fp, SurvX:::ParCSCFit))
  nd <- data.frame(x = c(0, 1), z = 0); tt <- c(0.5, 1, 3)
  F1 <- predict(fp, nd, times = tt, cause = 1); F2 <- predict(fp, nd, times = tt, cause = 2)
  S <- predict(fp, nd, times = tt, type = "survival")
  expect_close(F1 + F2 + S, matrix(1, 2, 3), 1e-8)
  # independent integration with the fitted hazards
  m1 <- fp@models[[1]]; m2 <- fp@models[[2]]
  hz <- function(m, u, row) drop(predict(m, nd[row, , drop = FALSE], times = u, type = "hazard"))
  ch <- function(m, u, row) drop(predict(m, nd[row, , drop = FALSE], times = u, type = "cumhaz"))
  ref <- integrate(function(u) hz(m1, u, 2) * exp(-ch(m1, u, 2) - ch(m2, u, 2)), 0, 3, rel.tol = 1e-10)$value
  expect_close(F1[2, 3], ref, 1e-7)
  expect_identical(predict(fp, nd, times = rev(tt), cause = 1)[, 3:1], F1)
  pse <- predict(fp, nd, times = tt, cause = 1, se = TRUE)
  expect_true(all(is.finite(pse$se) & pse$se > 0))
  expect_close(as.numeric(logLik(fp)), as.numeric(logLik(m1)) + as.numeric(logLik(m2)), 1e-12)
  expect_output(print(fp), "coherent by construction")
  expect_output(print(summary(fp)))
  expect_s3_class(estimates(fp), "data.frame")
})

test_that("parametric cause-specific fits work in cif_data(), cif_plot() and cv_compare()", {
  d <- sim_cr(300); fm <- Event(time, status) ~ x + z
  fp <- csc(fm, d, baseline = "weibull")
  cd <- cif_data(fp, times = c(0.5, 1, 2))
  expect_identical(nrow(cd), 6L)
  expect_close(cd$cif[cd$cause == 1], predict(fp, data.frame(x = mean(d$x), z = mean(d$z)), times = c(0.5, 1, 2)), 1e-10)
  ma <- cif_data(fp, type = "marginal", times = c(1, 2), cause = 1, conf.int = FALSE)
  expect_close(ma$cif, colMeans(predict(fp, times = c(1, 2))), 1e-10)
  expect_no_error(with_null_device({
    cif_plot(fp)
    cif_plot(list(weibull = fp, Cox = csc(fm, d)), cause = 1, B = 20)
  }))
  skip_on_cran()
  cc <- cv_compare(fm, d, times = c(1, 2), folds = 3, models = list(
    weibull = function(f, d) csc(f, d, baseline = "weibull"), Cox = function(f, d) csc(f, d)))
  expect_true(all(is.finite(cc@table$brier_cv)))
})

test_that("hazard_np() agrees with survival::survfit, also with delayed entry", {
  skip_if_not_installed("survival")
  for (entry in c(FALSE, TRUE)) {
    d <- sim_surv(entry = entry)
    hz <- if (entry) hazard_np(Event(time, status, entry = entry) ~ x, d, smooth = FALSE)
          else hazard_np(Event(time, status) ~ x, d, smooth = FALSE)
    sf <- if (entry) survival::survfit(survival::Surv(entry, time, status) ~ x, d)
          else survival::survfit(survival::Surv(time, status) ~ x, d)
    ss <- summary(sf, censored = FALSE)
    ours <- do.call(rbind, hz@cumhaz)
    expect_close(ours$cumhaz, ss$cumhaz, 1e-10)
    expect_close(ours$se, ss$std.chaz, 1e-8)
  }
})

test_that("the kernel-smoothed hazard matches the formula of Ramlau-Hansen with boundary kernels", {
  d <- sim_surv()
  hz <- hazard_np(Event(time, status) ~ 1, d, bw = 1)
  h <- hz@hazard[[1]]; inc <- hz@cumhaz[[1]]
  R <- max(inc$time); b <- 1
  K <- function(z) ifelse(abs(z) <= 1, 0.75 * (1 - z^2), 0)
  KB <- function(z, q) ifelse(z >= -1 & z <= q, 12 / (1 + q)^4 * (z + 1) * (z * (1 - 2 * q) + (3 * q^2 - 2 * q + 1) / 2), 0)
  ref <- sapply(h$time, function(x) {
    z <- (x - inc$time) / b
    w <- if (x < b) KB(z, x / b) else if (x > R - b) KB(-z, (R - x) / b) else K(z)
    sum(w * inc$n.event / inc$n.risk) / b
  })
  expect_close(h$hazard, ref, 1e-10)
  auto <- hazard_np(Event(time, status) ~ 1, d)
  expect_true(auto@bandwidth$bandwidth > 0)
  expect_true(max(auto@hazard[[1]]$time) < max(inc$time))       # stops where the risk set is small
  expect_output(print(auto), "bandwidth")
  expect_no_error(with_null_device({ plot(auto); plot(auto, type = "cumhaz") }))
})

test_that("printing, plotting and delayed entry guards", {
  d <- sim_surv(entry = TRUE)
  f <- parsurv(Event(time, status, entry = entry) ~ x + z, d)
  expect_output(print(f), "Reading")
  expect_output(print(f), "delayed entry: yes")
  expect_output(print(summary(f)), "Weibull shape")
  e <- estimates(f)
  expect_identical(sum(e$role == "beta"), 2L)
  expect_s3_class(estimates(cox(Event(time, status, entry = entry) ~ x, d)), "data.frame")   # Cox: (entry, time] since 0.8.0
  expect_error(fg(Event(time, status, entry = entry) ~ x, d), "delayed entry")
  expect_error(parsurv(Event(time, status) ~ x, d, baseline = "lognormal", scale = "PH"), "AFT")
  expect_output(print(Event(c(2, 3), c(1, 0), entry = c(0, 1))), "\\(1")
  fl <- list(Weibull = f, lognormal = parsurv(Event(time, status, entry = entry) ~ x + z, d, baseline = "lognormal"))
  expect_no_error(with_null_device({
    hazard_plot(fl); hazard_plot(fl, type = "cumhaz"); hazard_plot(fl, type = "survival")
    hazard_plot(csc(Event(time, status) ~ x + z, sim_cr(300), baseline = "gompertz"))
  }))
  a <- AIC(f, fl$lognormal)
  expect_identical(nrow(a), 2L)
})

# ---- spline (Royston-Parmar) and piecewise exponential baselines ----------------------------------
rp_basis_R <- function(x, k) {
  kmin <- k[1]; kmax <- k[length(k)]
  B <- cbind(x)
  for (kj in k[-c(1, length(k))]) {
    l <- (kmax - kj) / (kmax - kmin)
    B <- cbind(B, pmax(x - kj, 0)^3 - l * pmax(x - kmin, 0)^3 - (1 - l) * pmax(x - kmax, 0)^3)
  }
  B
}
rp_dbasis_R <- function(x, k) {
  kmin <- k[1]; kmax <- k[length(k)]
  B <- cbind(rep(1, length(x)))
  for (kj in k[-c(1, length(k))]) {
    l <- (kmax - kj) / (kmax - kmin)
    B <- cbind(B, 3 * (pmax(x - kj, 0)^2 - l * pmax(x - kmin, 0)^2 - (1 - l) * pmax(x - kmax, 0)^2))
  }
  B
}

test_that("the Royston-Parmar model maximises its likelihood, also with delayed entry", {
  for (entry in c(FALSE, TRUE)) {
    d <- sim_surv(entry = entry)
    fm <- if (entry) Event(time, status, entry = entry) ~ x + z else Event(time, status) ~ x + z
    f <- parsurv(fm, d, baseline = "rp", df = 3)
    expect_true(f@converged)
    k <- f@knots; X <- as.matrix(d[, c("x", "z")])
    expect_length(k, 4L)
    expect_close(k[c(1, 4)], range(log(d$time[d$status == 1])), 1e-12)
    ll <- function(th) {
      u <- drop(th[1] + X %*% th[2:3]); g <- th[4:6]
      s <- u + drop(rp_basis_R(log(d$time), k) %*% g); sd <- drop(rp_dbasis_R(log(d$time), k) %*% g)
      if (any(sd[d$status == 1] <= 0)) return(-Inf)
      e <- if (entry) d$entry else rep(0, nrow(d))
      se <- ifelse(e > 0, u + drop(rp_basis_R(log(pmax(e, 1e-300)), k) %*% g), -Inf)
      sum(d$status * (s + log(pmax(sd, 1e-300)) - log(d$time)) - exp(s) + ifelse(e > 0, exp(se), 0))
    }
    expect_close(as.numeric(logLik(f)), ll(coef(f)), 1e-8)
    o <- stats::optim(coef(f), function(p) -ll(p), method = "BFGS", control = list(reltol = 1e-14, maxit = 500))
    expect_lte(-o$value - as.numeric(logLik(f)), 1e-6)
  }
  d <- sim_surv()
  w <- parsurv(Event(time, status) ~ x + z, d, baseline = "weibull")
  r1 <- parsurv(Event(time, status) ~ x + z, d, baseline = "rp", df = 1)
  expect_close(as.numeric(logLik(r1)), as.numeric(logLik(w)), 1e-8)       # df = 1 is the Weibull model
  expect_close(coef(r1)[4], exp(coef(w)[4]), 1e-5)
})

test_that("the Royston-Parmar model agrees with flexsurv::flexsurvspline", {
  skip_if_not_installed("flexsurv")
  d <- sim_surv()
  f <- parsurv(Event(time, status) ~ x + z, d, baseline = "rp", df = 3)
  r <- flexsurv::flexsurvspline(survival::Surv(time, status) ~ x + z, data = d, k = 2, scale = "hazard")
  expect_close(as.numeric(logLik(f)), as.numeric(logLik(r)), 1e-5)
  expect_close(unname(coef(f)[c(1, 4, 5, 6, 2, 3)]), unname(coef(r)), 1e-4)
})

test_that("the piecewise exponential model equals a Poisson GLM on split data", {
  d <- sim_surv()
  cuts <- c(1, 2.5, 4)
  f <- parsurv(Event(time, status) ~ x + z, d, baseline = "pwexp", cuts = cuts)
  # split the follow-up at the cut points: one row per subject and interval (lo, hi]
  br <- c(0, cuts, Inf)
  sp <- do.call(rbind, lapply(seq_len(length(br) - 1L), function(j) {
    r <- d[d$time > br[j], , drop = FALSE]
    data.frame(x = r$x, z = r$z, interval = j, expo = pmin(r$time, br[j + 1L]) - br[j],
               status = as.integer(r$status == 1L & r$time <= br[j + 1L]))
  }))
  g <- stats::glm(status ~ factor(interval) + x + z + offset(log(expo)), family = stats::poisson(), data = sp)
  cg <- coef(g)
  expect_close(coef(f)[c("(Intercept)", "x", "z")], cg[c("(Intercept)", "x", "z")], 1e-6)
  expect_close(unname(coef(f)[4:6]), unname(cg[grep("interval", names(cg))]), 1e-6)
  expect_output(print(summary(f)), "rate ratio against the first interval")
})

test_that("spline and piecewise exponential cause-specific models are coherent and comparable", {
  d <- sim_cr(400); fm <- Event(time, status) ~ x + z
  for (b in list("rp", "pwexp", c("rp", "weibull"))) {
    fp <- csc(fm, d, baseline = b)
    tt <- c(0.3, 1, 2.5)
    s <- predict(fp, d[1:4, ], times = tt, cause = 1) + predict(fp, d[1:4, ], times = tt, cause = 2) +
         predict(fp, d[1:4, ], times = tt, type = "survival")
    expect_close(s, matrix(1, 4, 3), 1e-7)
  }
  fq <- parsurv(Event(time, status) ~ x, sim_surv(), baseline = "rp")
  q <- predict(fq, data.frame(x = 0), type = "quantile", p = 0.5)
  expect_close(predict(fq, data.frame(x = 0), times = q), 0.5, 1e-8)
  skip_on_cran()
  cm <- compare_models(fm, d, times = c(1, 2), folds = 3, hazard = list("weibull", c("rp", 2)))
  expect_identical(nrow(cm@table), 2L * 9L)
  expect_true(all(is.finite(cm@table$brier_cv)))
  expect_true(all(is.finite(cm@fit$aic[8:9])))
  expect_output(print(cm), "parametric Weibull cause-specific")
})

test_that("generalized gamma agrees with flexsurv and nests the Weibull and log-normal", {
  skip_if_not_installed("flexsurv")
  d <- sim_surv(entry = TRUE)
  f <- parsurv(Event(time, status, entry = entry) ~ x + z, d, baseline = "gengamma")
  r <- flexsurv::flexsurvreg(survival::Surv(entry, time, status) ~ x + z, data = d, dist = "gengamma")
  expect_true(f@converged)
  expect_equal(f@scale, "AFT")
  expect_close(as.numeric(logLik(f)), r$loglik, 1e-5)
  rc <- r$res.t[, "est"]                                # mu, sigma (log), Q, covariates
  expect_close(coef(f), c(rc["mu"], rc["x"], rc["z"], rc["sigma"], rc["Q"]), 5e-3)
  # the Weibull is Q = 1: the generalized gamma cannot fit worse
  w <- parsurv(Event(time, status, entry = entry) ~ x + z, d, baseline = "weibull", scale = "AFT")
  expect_gte(as.numeric(logLik(f)), as.numeric(logLik(w)) - 1e-6)
  # survival curve and median agree with flexsurv
  nd <- data.frame(x = 1, z = 0.3)
  expect_close(predict(f, nd, times = c(1, 2.5), type = "survival"),
               summary(r, newdata = nd, t = c(1, 2.5), type = "survival", tidy = TRUE)$est, 1e-3)
  expect_output(print(f), "generalized gamma")
})

test_that("generalized gamma works as a cause-specific baseline", {
  d <- sim_cr(400, 3)
  f <- csc(Event(time, status) ~ x + z, d, baseline = "gengamma")
  expect_true(all(vapply(f@models, function(m) m@converged, logical(1))))
  p <- c(predict(f, data.frame(x = 0, z = 0), times = 1:2, cause = 1), predict(f, data.frame(x = 0, z = 0), times = 1:2, cause = 2))
  expect_true(all(is.finite(p) & p > 0 & p < 1))
})
