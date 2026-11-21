# time_effects(): screening, choice of the form, final model, predictions and output

te_data <- function(seed = 4, n = 400) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.5); z <- rnorm(n)
  t <- rweibull(n, ifelse(x == 1, 0.7, 1.4), 3 * exp(-0.3 * z)); cens <- runif(n, 0, 6)
  data.frame(time = round(pmin(t, cens), 3) + 0.001, status = as.integer(t <= cens), x = x, z = z,
             g = sample(c("a", "b"), n, TRUE))
}

test_that("the screening test is the test of ph_test()", {
  d <- te_data()
  f <- cox(Event(time, status) ~ x + z, d)
  te <- time_effects(f, forms = c("constant", "log"))
  pt <- ph_test(f)
  u <- te@units[[1]]
  expect_close(u$screening_columns$chisq, pt@table$chisq[pt@table$term != "GLOBAL"], 1e-6)
  expect_close(u$global[["chisq"]], pt@table$chisq[pt@table$term == "GLOBAL"], 1e-6)
  expect_identical(u$screening$term, c("x", "z"))
  expect_true(u$screening$considered[1])          # x: strongly non-proportional
  expect_false(u$screening$considered[2])
})

test_that("the log form is the model of cox(..., tt = ) and the choice picks a time-varying form for x", {
  d <- te_data()
  f <- cox(Event(time, status) ~ x + z, d)
  te <- time_effects(f, forms = c("constant", "log"), fix = list(x = "log"))
  u <- te@units[[1]]
  tm <- stats::median(sort(unique(d$time[d$status == 1])))
  g <- cox(Event(time, status) ~ x + z + tt(x), d, tt = function(t) log(t / tm))
  expect_close(u$final$coefficients, coef(g), 1e-6)
  expect_close(u$final$vcov, vcov(g), 1e-5)
  expect_close(u$final$loglik, as.numeric(logLik(g)), 1e-7)
  e <- u$effects[u$effects$column == "x", ]
  b <- coef(g)
  expect_close(e$log_hr, b[["x"]] + b[["tt(x)"]] * log(e$time / tm), 1e-6)
  te2 <- time_effects(f)                         # default forms
  expect_true(te2@units[[1]]$selected[["x"]] != "constant")
  expect_identical(te2@units[[1]]$selected[["z"]], "constant")
  cm <- te2@units[[1]]$comparison
  expect_true(all(c("constant", "log", "piecewise", "spline") %in% cm$form[cm$term == "x"]))
  expect_equal(sum(cm$chosen[cm$term == "x"]), 1L)
  expect_true(all(cm$delta[cm$term == "x"] >= 0))
})

test_that("time_effects() agrees with survival::coxph(tt = ) for a spline form", {
  skip_if_not_installed("survival")
  d <- te_data(7)
  f <- cox(Event(time, status) ~ x + z, d)
  te <- time_effects(f, forms = c("constant", "spline"), fix = list(x = "spline"), df = 3)
  et <- sort(unique(d$time[d$status == 1]))
  kn <- stats::quantile(log(et), (1:2) / 3, names = FALSE); bk <- range(log(et))
  r <- survival::coxph(survival::Surv(time, status) ~ x + z + tt(x), data = d,
                       tt = function(x, t, ...) x * splines::ns(log(t), knots = kn, Boundary.knots = bk))
  expect_close(te@units[[1]]$final$coefficients, coef(r), 1e-6)
  expect_close(te@units[[1]]$final$vcov, vcov(r), 1e-5)
})

test_that("with every effect constant the fit and predictions reproduce cox(), with strata", {
  d <- te_data(5)
  f <- cox(Event(time, status) ~ x + z + strata(g), d)
  te <- time_effects(f, fix = list(x = "constant", z = "constant"))
  expect_close(te@units[[1]]$final$coefficients, coef(f), 1e-7)
  nd <- data.frame(x = c(0, 1, 1), z = c(0, 0.5, -1), g = c("a", "b", "a"))
  expect_close(predict(te, nd, times = c(0.5, 1, 2, 4)), predict(f, nd, times = c(0.5, 1, 2, 4), type = "survival"), 1e-7)
})

test_that("Fine-Gray: constant effects reproduce fg(); predictions reproduce the weighted Breslow estimate", {
  d <- sim_cr(400, 52)
  g <- fg(Event(time, status) ~ x + z, d)
  te <- time_effects(g, fix = list(x = "constant", z = "constant"))
  u <- te@units[[1]]
  expect_close(u$final$coefficients, coef(g), 1e-6)
  expect_close(u$final$vcov_model, g@vcov_model, 1e-6)
  expect_true(u$final$robust)
  nd <- data.frame(x = 0:1, z = c(0, 0.5)); tt <- c(0.5, 1, 2)
  e <- SurvX:::.fg_expand(d$time, d$status, NULL, 1L)
  bh <- SurvX:::.cox_baseline_raw(coef(g), e$stop, e$status, as.matrix(d[e$id, c("x", "z")]), e$w, "breslow", e$start, NULL)
  H <- SurvX:::.step_at(bh$time, bh$cumhaz, tt, 0)
  ref <- 1 - exp(-outer(exp(drop(as.matrix(nd) %*% coef(g))), H))
  expect_close(predict(te, nd, times = tt), ref, 1e-6)
  te2 <- time_effects(g)
  expect_identical(te2@kind, "fg")
  expect_output(print(te2), "Fine-Gray")
})

test_that("cause-specific models: one analysis per cause, screening as ph_test(), predictions as predict(csc)", {
  d <- sim_cr(500, 3)
  f <- csc(Event(time, status) ~ x + z, d)
  te <- time_effects(f, fix = list(x = "constant", z = "constant"))
  expect_length(te@units, 2L)
  pt <- ph_test(f)
  for (k in 1:2) {
    ref <- pt@table[pt@table$cause == k & pt@table$term != "GLOBAL", ]
    expect_close(te@units[[k]]$screening_columns$chisq, ref$chisq, 1e-6)
    expect_close(te@units[[k]]$final$coefficients, f@models[[k]]$coefficients, 1e-7)
  }
  nd <- data.frame(x = 0:1, z = c(0, 0.5)); tt <- c(0.5, 1, 2)
  for (k in 1:2) expect_close(predict(te, nd, times = tt, cause = k), predict(f, nd, times = tt, cause = k), 1e-6)
})

test_that("own forms, given breaks, factors and the output methods", {
  d <- te_data(6)
  d$grp <- factor(sample(c("u", "v", "w"), nrow(d), TRUE))
  f <- cox(Event(time, status) ~ x + z + grp, d)
  te <- time_effects(f, forms = list("constant", "log", sqrt = function(t) sqrt(t)), screen = "all")
  cm <- te@units[[1]]$comparison
  expect_true("sqrt" %in% cm$form)
  expect_equal(cm$k[cm$term == "grp" & cm$form == "log"], length(coef(f)) + 2L)   # two levels, one column each
  tp <- time_effects(f, fix = list(x = "piecewise"), breaks = 1)
  e <- tp@units[[1]]$effects
  ex <- e[e$column == "x", ]
  expect_equal(length(unique(round(ex$hr[ex$time <= 1], 10))), 1L)
  expect_equal(length(unique(round(ex$hr[ex$time > 1], 10))), 1L)
  expect_output(print(te), "Test of a constant effect")
  expect_output(summary(te), "Cautions")
  for (w in c("effects", "screening", "forms", "coefficients")) expect_s3_class(as.data.frame(te, what = w), "data.frame")
  with_null_device(expect_no_error(plot(te)))
  with_null_device(expect_no_error(plot(te, term = "z", scale = "log")))
  html <- tempfile(fileext = ".html")
  expect_message(report(te, file = html), "report written")
  s <- paste(readLines(html), collapse = "\n")
  expect_match(s, "Grambsch")
  expect_match(s, "data:image/png;base64")
  expect_true(any(grepl("Grambsch", capture.output(sx_references(te)))))
})

test_that("the bootstrap is reproducible, leaves the RNG alone and gives selection frequencies", {
  d <- te_data(8, 250)
  f <- cox(Event(time, status) ~ x + z, d)
  set.seed(99); before <- runif(1); set.seed(99)
  a <- time_effects(f, forms = c("constant", "log"), bootstrap = 20, seed = 3)
  expect_identical(runif(1), before)
  b <- time_effects(f, forms = c("constant", "log"), bootstrap = 20, seed = 3)
  expect_identical(a@bootstrap, b@bootstrap)
  fq <- as.data.frame(a, what = "stability")
  expect_close(tapply(fq$share, fq$term, sum), c(1, 1), 1e-12)
  expect_output(print(a), "Stability of the choice")
})

test_that("time_effects() refuses what it cannot do", {
  d <- te_data()
  f <- cox(Event(time, status) ~ x + z, d)
  expect_error(time_effects(f, forms = "cubic"), "unknown form")
  expect_error(time_effects(f, forms = list("constant", function(t) t)), "needs a name")
  expect_error(time_effects(f, terms = "w"), "terms not in the model")
  expect_error(time_effects(f, fix = list(x = "spline"), forms = c("constant", "log")), "forms not among")
  expect_error(time_effects(f, bootstrap = 5), "at least 20")
  expect_error(time_effects(cox(Event(time, status) ~ x + tt(x), d, tt = log)), "tt\\(\\) terms")
  expect_error(time_effects(km(Event(time, status) ~ 1, d)), "takes a model")
})
