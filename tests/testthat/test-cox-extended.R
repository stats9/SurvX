# strata, (start, stop] data and the robust variance, against survival::coxph

sim_strata <- function(n = 400, seed = 41, ties = FALSE) {
  set.seed(seed)
  g <- sample(c("A", "B", "C"), n, TRUE); x <- rbinom(n, 1, 0.5); z <- rnorm(n)
  base <- c(A = 0.2, B = 0.5, C = 1)[g]
  t <- rweibull(n, c(A = 0.8, B = 1.2, C = 1.6)[g], 1 / base * exp(-(0.5 * x - 0.3 * z)))
  cens <- runif(n, 0, 6)
  d <- data.frame(time = pmin(t, cens), status = as.integer(t <= cens), x = x, z = z, g = g,
                  id = rep(seq_len(n / 2), each = 2))
  if (ties) d$time <- round(d$time, 1) + 0.1
  d
}

# split every subject at the cut points: (start, stop] rows with the same covariates
split_rows <- function(d, cuts) {
  out <- lapply(seq_len(nrow(d)), function(i) {
    br <- c(0, cuts[cuts < d$time[i]], d$time[i])
    k <- length(br) - 1L
    data.frame(d[rep(i, k), , drop = FALSE], start = br[-(k + 1L)], stop = br[-1L],
               event = c(rep(0L, k - 1L), d$status[i]), row.names = NULL)
  })
  do.call(rbind, out)
}

test_that("stratified Cox agrees with survival::coxph", {
  skip_if_not_installed("survival")
  d <- sim_strata()
  for (ti in c("efron", "breslow")) {
    f <- cox(Event(time, status) ~ x + z + strata(g), d, ties = ti)
    r <- survival::coxph(survival::Surv(time, status) ~ x + z + strata(g), d, ties = ti)
    expect_close(coef(f), coef(r), 1e-6)
    expect_close(vcov(f), vcov(r), 1e-6)
    expect_close(as.numeric(logLik(f)), r$loglik[2], 1e-8)
    expect_close(residuals(f), residuals(r, "martingale"), 1e-6)
  }
  # survival curves per stratum
  f <- cox(Event(time, status) ~ x + z + strata(g), d, ties = "breslow")
  r <- survival::coxph(survival::Surv(time, status) ~ x + z + strata(g), d, ties = "breslow")
  nd <- data.frame(x = c(0, 1, 1), z = c(0, 0.5, -1), g = c("A", "B", "C"))
  tt <- c(0.5, 1, 2)
  ours <- predict(f, nd, times = tt, type = "survival")
  ref <- t(vapply(seq_len(nrow(nd)), function(i) {
    s <- survival::survfit(r, newdata = nd[i, ])
    summary(s, times = tt, extend = TRUE)$surv
  }, numeric(length(tt))))
  expect_close(ours, ref, 1e-6)
  expect_error(predict(f, data.frame(x = 0, z = 0, g = "D"), times = 1, type = "survival"), "not in the data")
  expect_output(print(f), "Strata: strata\\(g\\)")
  expect_output(print(f), "Kalbfleisch & Prentice")
})

test_that("(start, stop] data: split data give the unsplit fit, and agree with coxph", {
  skip_if_not_installed("survival")
  d <- sim_strata(300, 42)
  s <- split_rows(d, c(0.5, 1, 2))
  f0 <- cox(Event(time, status) ~ x + z, d)
  f1 <- cox(Event(stop, event, entry = start) ~ x + z, s)
  expect_close(coef(f1), coef(f0), 1e-8)
  expect_close(vcov(f1), vcov(f0), 1e-8)
  expect_close(as.numeric(logLik(f1)), as.numeric(logLik(f0)), 1e-10)
  # time-varying covariate: x switches on at time 1
  s$xt <- s$x * (s$start >= 1)
  f2 <- cox(Event(stop, event, entry = start) ~ xt + z, s)
  r2 <- survival::coxph(survival::Surv(start, stop, event) ~ xt + z, s)
  expect_close(coef(f2), coef(r2), 1e-6)
  expect_close(vcov(f2), vcov(r2), 1e-6)
  expect_close(residuals(f2), residuals(r2, "martingale"), 1e-6)
  expect_output(print(f2), "Andersen & Gill")
  expect_error(cox(Event(stop, event, entry = stop) ~ x, s), "smaller than time")
})

test_that("delayed entry and strata together agree with coxph", {
  skip_if_not_installed("survival")
  d <- sim_strata(400, 43)
  d$entry <- runif(nrow(d), 0, 0.3) * d$time
  f <- cox(Event(time, status, entry = entry) ~ x + z + strata(g), d)
  r <- survival::coxph(survival::Surv(entry, time, status) ~ x + z + strata(g), d)
  expect_close(coef(f), coef(r), 1e-6)
  expect_close(vcov(f), vcov(r), 1e-6)
  expect_close(residuals(f), residuals(r, "martingale"), 1e-6)
  pt <- ph_test(f, transform = "identity")
  zr <- survival::cox.zph(r, transform = "identity", terms = FALSE)$table
  expect_close(pt@table$chisq, zr[, "chisq"], 1e-5)
})

test_that("robust variance agrees with coxph(cluster = )", {
  skip_if_not_installed("survival")
  d <- sim_strata(400, 44)
  f <- cox(Event(time, status) ~ x + z + cluster(id), d)
  r <- survival::coxph(survival::Surv(time, status) ~ x + z, d, cluster = id)
  expect_close(coef(f), coef(r), 1e-6)
  expect_close(vcov(f), vcov(r), 1e-6)
  expect_close(f@vcov_model, r$naive.var, 1e-6)
  expect_output(print(f), "clustered on 200 clusters")
  expect_output(print(f), "Lin & Wei")
  # robust = TRUE without clusters: each row is its own cluster
  f1 <- cox(Event(time, status) ~ x + z, d, robust = TRUE)
  r1 <- survival::coxph(survival::Surv(time, status) ~ x + z, d, robust = TRUE)
  expect_close(vcov(f1), vcov(r1), 1e-6)
  # tied times, Efron and Breslow
  dt <- sim_strata(400, 45, ties = TRUE)
  for (ti in c("efron", "breslow")) {
    f2 <- cox(Event(time, status) ~ x + z + strata(g) + cluster(id), dt, ties = ti)
    r2 <- survival::coxph(survival::Surv(time, status) ~ x + z + strata(g), dt, cluster = id, ties = ti)
    expect_close(coef(f2), coef(r2), 1e-6)
    expect_close(vcov(f2), vcov(r2), 1e-5)
  }
})

test_that("csc() with strata, entry and clusters", {
  skip_if_not_installed("survival")
  d <- sim_cr(400, 46)
  d$g <- rep(c("u", "v"), 200); d$id <- rep(1:100, each = 4)
  f <- csc(Event(time, status) ~ x + z + strata(g) + cluster(id), d)
  for (k in 1:2) {
    r <- survival::coxph(survival::Surv(time, status == k) ~ x + z + strata(g), d, cluster = id)
    expect_close(f@models[[k]]$coefficients, coef(r), 1e-6)
    expect_close(f@models[[k]]$vcov, vcov(r), 1e-6)
  }
  p <- predict(f, data.frame(x = c(0, 1), z = 0, g = c("u", "v")), times = c(0.5, 1, 2), cause = 1)
  expect_true(all(is.finite(p) & p >= 0 & p <= 1))
  expect_true(all(diff(p[1, ]) >= 0))
  cd <- cif_data(f, cause = 1, B = 20)
  expect_setequal(unique(cd$profile), c("typical covariate values, u", "typical covariate values, v"))
  expect_true(all(cd$lower <= cd$cif + 1e-12 & cd$cif <= cd$upper + 1e-12, na.rm = TRUE))
  cn <- cif_data(f, newdata = data.frame(x = 1, z = 0, g = "v"), cause = 2, B = 0, conf.int = FALSE)
  expect_close(cn$cif, predict(f, data.frame(x = 1, z = 0, g = "v"), times = cn$time, cause = 2), 1e-10)
  with_null_device(expect_no_error(cif_plot(f, B = 10, type = "marginal")))
  # the general CIF code reproduces the C++ code on an unstratified fit
  f0 <- csc(Event(time, status) ~ x + z, d)
  nd <- data.frame(x = c(0, 1), z = c(-1, 1))
  for (pl in c(TRUE, FALSE)) {
    f0@product_limit <- pl
    for (k in 1:2)
      expect_close(SurvX:::.csc_cif_general(f0, SurvX:::.sx_newX(f0, nd), NULL, c(0.3, 1, 2.5), k),
                   predict(f0, nd, times = c(0.3, 1, 2.5), cause = k), 1e-10)
  }
})

test_that("strata() and cluster() are refused where they do not apply", {
  d <- sim_cr(200, 47); d$g <- rep(c("u", "v"), 100)
  expect_error(fg(Event(time, status) ~ x + strata(g), d), "available in cox\\(\\) and csc\\(\\)")
  expect_error(directcif(Event(time, status) ~ x + strata(g), d), "available in cox\\(\\) and csc\\(\\)")
  expect_error(csc(Event(time, status) ~ x + strata(g), d, baseline = "weibull"), "available in cox\\(\\) and csc\\(\\)")
  expect_error(cox(Event(time, dead) ~ x * strata(g), d), "interaction")
})
