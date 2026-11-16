test_that("cif_data() reproduces predict() at typical covariate values and on average", {
  d <- sim_cr(); fm <- Event(time, status) ~ x + z; tp <- c(0.5, 1, 2)
  fj <- directcif(fm, d); ff <- fg(fm, d); fs <- csc(fm, d)
  prof <- data.frame(x = mean(d$x), z = mean(d$z))
  expect_close(cif_data(fj, conf.int = FALSE, times = tp)$cif,
               c(predict(fj, prof, times = tp, cause = 1), predict(fj, prof, times = tp, cause = 2)), 1e-10)
  expect_close(cif_data(ff, cause = 1, conf.int = FALSE, times = tp)$cif, predict(ff, prof, times = tp), 1e-10)
  expect_close(cif_data(fs, conf.int = FALSE, times = tp)$cif,
               c(predict(fs, prof, times = tp, cause = 1), predict(fs, prof, times = tp, cause = 2)), 1e-10)
  expect_close(cif_data(fj, type = "marginal", conf.int = FALSE, cause = 1, times = tp)$cif,
               colMeans(predict(fj, times = tp, cause = 1)), 1e-10)
})

test_that("unadjusted cause-specific Cox equals Aalen-Johansen", {
  d <- sim_cr(); tp <- c(0.5, 1, 2)
  fs <- csc(Event(time, status) ~ x + z, d)
  ajr <- summary(aj(Event(time, status) ~ 1, d), times = tp)@tables[[1]]
  expect_close(cif_data(fs, type = "unadjusted", conf.int = FALSE, times = tp)$cif, c(ajr$cif.1, ajr$cif.2), 1e-10)
})

test_that("bootstrap bands are reproducible and contain the estimate", {
  d <- sim_cr(); ff <- fg(Event(time, status) ~ x + z, d)
  b1 <- cif_data(ff, cause = 1, B = 30, seed = 7, times = c(1, 2))
  b2 <- cif_data(ff, cause = 1, B = 30, seed = 7, times = c(1, 2))
  expect_identical(c(b1$lower, b1$upper), c(b2$lower, b2$upper))
  expect_true(all(b1$lower <= b1$cif & b1$cif <= b1$upper))
})

test_that("plots of every class draw without error", {
  d <- sim_cr(); nd <- data.frame(x = c(0, 1), z = 0); fm <- Event(time, status) ~ x + z
  f1 <- cox(Event(time, dead) ~ x + z, d); f2 <- fg(fm, d); f3 <- csc(fm, d); f4 <- directcif(fm, d)
  expect_no_error(with_null_device({
    plot(f1, newdata = nd); plot(f2, newdata = nd, B = 20); plot(f3, B = 20); plot(f4, newdata = nd)
    cif_plot(list(direct = f4, "Fine-Gray" = f2, CSC = f3), B = 20, layout = "overlay")
    cif_plot(f4, layout = "stacked")
    cif_plot(f4, type = "unadjusted", reference = TRUE)
    cif_plot(f4, main = paste(rep("a long title that must be wrapped", 6), collapse = " "))
  }))
  skip_if_not_installed("ggplot2")
  expect_no_error(with_null_device(print(cif_plot(list(direct = f4, "Fine-Gray" = f2), B = 20, engine = "ggplot2"))))
})
