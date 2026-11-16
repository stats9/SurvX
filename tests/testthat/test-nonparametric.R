test_that("km() runs, prints and plots", {
  d <- sim_cr()
  k <- km(Event(time, dead) ~ x, d)
  expect_output(print(k))
  expect_output(print(summary(k)))
  s <- k@curves[[1]]$surv
  expect_true(all(diff(s) <= 1e-12) && all(s >= 0 & s <= 1))
  with_null_device(plot(k))
})

test_that("km() agrees with survival::survfit and survdiff", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed")
  k <- km(Event(time, dead) ~ sex, d)
  ss <- summary(survival::survfit(survival::Surv(time, dead) ~ sex, d), censored = FALSE)
  ours <- do.call(rbind, k@curves); ours <- ours[ours$n.event > 0, ]
  expect_close(ours$surv, ss$surv)
  expect_close(ours$se, ss$std.err)
  expect_close(c(ours$lower, ours$upper), c(ss$lower, ss$upper))
  med <- summary(survival::survfit(survival::Surv(time, dead) ~ sex, d))$table
  expect_close(c(k@counts$median, k@counts$median.lower, k@counts$median.upper),
               c(med[, "median"], med[, "0.95LCL"], med[, "0.95UCL"]), 1e-9)
  expect_close(k@test$statistic, survival::survdiff(survival::Surv(time, dead) ~ sex, d)$chisq)
})

test_that("rmst() agrees with survival", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed")
  tb <- summary(survival::survfit(survival::Surv(time, dead) ~ 1, d), rmean = 200)$table
  expect_close(rmst(d$time, d$dead, 200), c(tb[["rmean"]], tb[["se(rmean)"]]))
})

test_that("aj() gives CIFs that add up with the survival to one", {
  d <- sim_cr()
  a <- aj(Event(time, status) ~ 1, d)
  cv <- a@curves[[1]]
  expect_close(cv$surv + cv$cif.1 + cv$cif.2, rep(1, nrow(cv)), 1e-10)
  expect_output(print(a))
  with_null_device(plot(a))
})

test_that("aj() and the Gray test agree with cmprsk", {
  d <- mgus_cr(); skip_if(is.null(d), "survival not installed"); skip_if_not_installed("cmprsk")
  a <- aj(Event(time, status) ~ sex, d)
  ci <- cmprsk::cuminc(d$time, d$status, d$sex)
  tp <- c(50, 100, 200); tt <- cmprsk::timepoints(ci, tp)
  ours <- summary(a, times = tp)@tables
  expect_close(ours[["F"]]$cif.1, tt$est["F 1", ])
  expect_close(ours[["M"]]$cif.2, tt$est["M 2", ])
  expect_close(a@tests$statistic, ci$Tests[, "stat"])
})
