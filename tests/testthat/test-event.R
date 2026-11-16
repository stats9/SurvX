test_that("Event() builds, subsets and prints", {
  y <- Event(c(1, 2, 3, 4, 5), c(0L, 1L, 2L, 1L, 0L))
  expect_s3_class(y, "Event")
  expect_s3_class(y[2:4], "Event")
  expect_identical(nrow(unclass(y[2:4])), 3L)
  expect_output(print(y))
  expect_false(any(is.na(y)))
})

test_that("Event() in a formula gives the same fit as the data columns", {
  d <- sim_cr()
  f1 <- fg(Event(time, status) ~ x + z, d)
  d2 <- d; names(d2)[1:4] <- c("t", "a", "b", "s")
  f2 <- fg(Event(t, s) ~ a + b, d2)
  expect_close(coef(f1), coef(f2), 1e-12)
})
