# Shared data and helpers for the tests.

# Largest absolute difference scaled by max(1, |reference|): the measure used by
# inst/validation/check_SurvX.R, so that both report the same thing.
expect_close <- function(ours, ref, tol = 1e-6) {
  ours <- as.numeric(ours); ref <- as.numeric(ref)
  testthat::expect_identical(length(ours), length(ref))
  ok <- is.finite(ours) & is.finite(ref)
  testthat::expect_identical(is.finite(ours), is.finite(ref))
  d <- if (any(ok)) max(abs(ours[ok] - ref[ok]) / pmax(1, abs(ref[ok]))) else 0
  testthat::expect_lte(d, tol)
}

# Simulated competing risks data with two causes and independent censoring (no package needed).
sim_cr <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.5); z <- runif(n, -1, 1)
  t1 <- rexp(n, 0.5 * exp(0.5 * x)); t2 <- rexp(n, 0.3 * exp(0.3 * z)); cens <- runif(n, 0, 4)
  d <- data.frame(time = pmin(t1, t2, cens), x = x, z = z)
  d$status <- ifelse(d$time == cens, 0L, ifelse(t1 < t2, 1L, 2L))
  d$dead <- as.integer(d$status > 0)
  d
}

# mgus2 from survival with progression (1) and death without progression (2) as competing events,
# the data of the validation script. NULL when survival is not installed.
mgus_cr <- function() {
  if (!requireNamespace("survival", quietly = TRUE)) return(NULL)
  m <- survival::mgus2
  m$time <- ifelse(m$pstat == 0, m$futime, m$ptime)
  m$status <- ifelse(m$pstat == 0, 2L * m$death, 1L)
  d <- stats::na.omit(data.frame(time = m$time, status = m$status, dead = as.integer(m$status > 0),
                                 sex = m$sex, age = m$age - 70, hgb = m$hgb - 13, mspike = m$mspike))
  d[d$time > 0, ]
}

# Draw on a throw-away device.
with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}
