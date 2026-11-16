# Fine-Gray with delayed entry (Geskus 2011)

sim_cr_entry <- function(n = 400, seed = 51) {
  d <- sim_cr(n, seed)
  set.seed(seed + 1)
  d$entry <- runif(n, 0, 0.4) * d$time
  d$id <- seq_len(n)
  d
}

test_that("the weighted-data form reproduces the Fine-Gray fit without entry", {
  d <- sim_cr(400, 52)
  f <- fg(Event(time, status) ~ x + z, d)
  e <- SurvX:::.fg_expand(d$time, d$status, NULL, 1L)
  X <- as.matrix(d[e$id, c("x", "z")])
  r <- SurvX:::.cox_fit(e$stop, e$status, X, e$w, "breslow", 50L, 1e-10, e$start, NULL)
  expect_close(r$coefficients, coef(f), 1e-7)
  expect_close(r$vcov, f@vcov_model, 1e-7)
})

test_that("fg() with delayed entry agrees with survival::finegray + coxph", {
  skip_if_not_installed("survival")
  d <- sim_cr_entry()
  f <- fg(Event(time, status, entry = entry) ~ x + z, d)
  fgd <- survival::finegray(survival::Surv(entry, time, factor(status, 0:2)) ~ x + z + id, data = d, etype = "1", id = id)
  r <- survival::coxph(survival::Surv(fgstart, fgstop, fgstatus) ~ x + z, data = fgd, weights = fgwt,
                       ties = "breslow", cluster = id)
  expect_close(coef(f), coef(r), 1e-6)
  expect_close(vcov(f), vcov(r), 1e-5)
  expect_close(f@vcov_model, r$naive.var, 1e-5)
  p <- predict(f, data.frame(x = 0:1, z = 0), times = c(0.5, 1, 2))
  expect_true(all(is.finite(p) & p >= 0 & p <= 1))
  expect_output(print(f), "Geskus")
  cd <- cif_data(f, B = 20)
  expect_true(all(cd$cause == 1L))
  expect_close(cd$cif, predict(f, SurvX:::.profile(f)$X |> as.data.frame(), times = cd$time)[1, ], 1e-8)
  with_null_device(expect_no_error(cif_plot(f, B = 10)))
})
