# Further performance measures: concordance, standard errors and paired differences, O/E, smooth
# calibration, decision curves, likelihood-ratio tests.

# C index over all pairs, written directly from the definition (reference for the O(n log n) core)
ref_cindex <- function(score, time, status, cause = 1L, tau = Inf, uno = FALSE) {
  n <- length(time)
  o <- order(time)
  G <- rep(1, n); g <- 1; atr <- n; i <- 1L
  while (i <= n) {                                  # Kaplan-Meier of censoring, G(t-), events first at ties
    j <- i; while (j < n && time[o[j + 1L]] == time[o[i]]) j <- j + 1L
    blk <- o[i:j]; G[blk] <- g
    ce <- sum(status[blk] == 0L); ev <- sum(status[blk] != 0L)
    if (ce > 0L) g <- g * (1 - ce / (atr - ev))
    atr <- atr - length(blk); i <- j + 1L
  }
  num <- den <- 0; Nu <- Du <- numeric(n)
  for (a in seq_len(n)) {
    if (status[a] != cause || time[a] > tau) next
    w <- if (uno) 1 / G[a]^2 else 1
    for (b in seq_len(n)) {
      if (b == a) next
      ctrl <- time[b] > time[a] || (time[b] == time[a] && status[b] == 0L) || (status[b] != 0L && status[b] != cause)
      if (!ctrl) next
      s <- if (score[a] > score[b]) 1 else if (score[a] == score[b]) 0.5 else 0
      num <- num + w * s; den <- den + w
      Nu[a] <- Nu[a] + w * s; Du[a] <- Du[a] + w; Nu[b] <- Nu[b] + w * s; Du[b] <- Du[b] + w
    }
  }
  C <- num / den
  c(C = C, se = sqrt(sum(((Nu - C * Du) / den)^2)))
}

test_that("cindex() equals the C index over all pairs (Harrell, Uno, Wolbers, ties, tau)", {
  set.seed(11)
  for (r in 1:6) {
    n <- 60
    time <- sample(1:15, n, TRUE); status <- sample(0:2, n, TRUE); score <- round(rnorm(n), 1)
    for (uno in c(FALSE, TRUE)) for (tau in c(Inf, 8)) {
      ours <- cindex(score, time, status, method = if (uno) "uno" else "harrell", tau = tau)
      ref <- ref_cindex(score, time, status, tau = tau, uno = uno)
      expect_close(c(ours$concordance, ours$se), ref, 1e-12)
    }
  }
})

test_that("cindex() agrees with survival::concordance for one event type", {
  skip_if_not_installed("survival")
  d <- sim_cr(300)
  d$ev <- as.integer(d$status == 1L)
  sc <- d$x + d$z
  ours <- cindex(sc, d$time, d$ev)
  ref <- survival::concordance(survival::Surv(time, ev) ~ sc, data = d, reverse = TRUE)
  expect_close(ours$concordance, unname(ref$concordance), 1e-10)
  expect_close(ours$se, sqrt(unname(ref$var)), 1e-8)
  f <- cox(Event(time, ev) ~ x + z, d)
  expect_close(cindex(f)$concordance, cindex(drop(predict(f, type = "lp")), d$time, d$ev)$concordance, 1e-12)
  expect_output(print(ours), "concordance")
  expect_close(c(ours$lower, ours$upper), ours$concordance + c(-1, 1) * stats::qnorm(0.975) * ours$se, 1e-12)
})

test_that("standard errors, paired differences and O/E of compare_models() and cv_compare()", {
  skip_on_cran()
  d <- sim_cr(400); fm <- Event(time, status) ~ x + z
  cm <- compare_models(fm, d, times = c(1, 2), folds = 3)
  tb <- cm@table
  for (nm in c("brier_se", "auc_se", "oe", "ici", "brier_diff_p", "auc_diff_se", "reference_cv")) {
    expect_true(nm %in% names(tb)); expect_true(paste0(sub("_cv$", "", nm), "_cv") %in% names(tb))
  }
  # the reference model at each time has the lowest Brier score and a zero difference
  for (t in c(1, 2)) {
    s <- tb[tb$time == t, ]
    expect_identical(s$reference[1L], s$model[which.min(s$brier)])
    expect_equal(s$brier_diff[s$model == s$reference[1L]], 0)
    expect_true(all(s$brier_diff >= 0, na.rm = TRUE))
  }
  expect_true(all(tb$brier_se > 0 & tb$auc_se > 0, na.rm = TRUE))
  expect_true(all(tb$oe_lower < tb$oe & tb$oe < tb$oe_upper, na.rm = TRUE))
  # the same measures from the risks directly (cv_compare predicts with the R fits)
  cc <- cv_compare(fm, d, times = 2, folds = 3, models = list(
    "Fine-Gray" = function(f, d) fg(f, d), "cause-specific Cox" = function(f, d) csc(f, d)))
  expect_true(all(is.finite(cc@table$brier_se_cv)))
  r <- drop(predict(fg(fm, d), times = 2))
  b <- brier_score(d$time, d$status, r, 2, se = TRUE)
  s <- cc@table[cc@table$model == "Fine-Gray", ]
  expect_close(c(b[["brier"]], b[["se"]]), c(s$brier, s$brier_se), 1e-9)
  a <- auc_t(d$time, d$status, r, 2, se = TRUE)
  expect_close(a[c("auc", "se")], c(s$auc, s$auc_se), 1e-9)
  z <- stats::qnorm(0.975)
  expect_close(a[c("lower", "upper")], c(max(0, a[["auc"]] - z * a[["se"]]), min(1, a[["auc"]] + z * a[["se"]])), 1e-12)
  expect_close(c(s$brier_lower, s$brier_upper), c(max(0, s$brier - z * s$brier_se), s$brier + z * s$brier_se), 1e-12)
  cal <- calibration(d$time, d$status, r, 2)
  expect_close(c(attr(cal, "oe")[1L], attr(cal, "ici")), c(s$oe, s$ici), 1e-9)
  expect_output(print(cm), "difference from best")
  with_null_device({ plot(cm, measure = "oe"); plot(cm, measure = "ici"); plot(cm, measure = "auc") })
})

test_that("calibration(): ICI, E50, E90 from the spline Cox calibration model of Austin et al. (2020)", {
  skip_if_not_installed("survival"); skip_if_not_installed("splines")
  d <- sim_cr(500)
  d$ev <- as.integer(d$status == 1L)
  f <- cox(Event(time, ev) ~ x + z, d)
  r <- 1 - drop(predict(f, d, times = 1, type = "survival"))
  cal <- calibration(d$time, d$ev, r, 1)
  cl <- log(-log(1 - r)); kn <- stats::quantile(cl, c(0.1, 0.5, 0.9))
  B <- splines::ns(cl, knots = kn[2L], Boundary.knots = kn[c(1L, 3L)])
  cf <- survival::coxph(survival::Surv(time, ev) ~ B, data = d, ties = "efron")
  bh <- survival::basehaz(cf, centered = FALSE)
  H1 <- max(c(0, bh$hazard[bh$time <= 1]))
  obs <- 1 - exp(-H1 * exp(drop(B %*% stats::coef(cf))))
  dd <- abs(obs - r)
  expect_close(c(attr(cal, "ici"), attr(cal, "e50"), attr(cal, "e90")),
               c(mean(dd), stats::median(dd), unname(stats::quantile(dd, 0.9))), 1e-5)
  expect_true(is.data.frame(attr(cal, "smooth")))
})

test_that("decision_curve(): net benefit from the Kaplan-Meier risk above the threshold", {
  skip_if_not_installed("survival")
  d <- sim_cr(400); d$ev <- as.integer(d$status == 1L)
  f <- cox(Event(time, ev) ~ x + z, d)
  r <- 1 - drop(predict(f, d, times = 1, type = "survival"))
  th <- c(0.1, 0.2, 0.3)
  dc <- decision_curve(list(cox = r), d$time, d$ev, t = 1, thresholds = th)
  km_risk <- function(i) { s <- summary(survival::survfit(survival::Surv(d$time[i], d$ev[i]) ~ 1), times = 1, extend = TRUE); 1 - s$surv }
  ref <- vapply(th, function(p) { i <- r >= p; pi <- mean(i); if (!any(i)) return(0); F <- km_risk(i); pi * F - pi * (1 - F) * p / (1 - p) }, 0)
  tb <- dc$table
  expect_close(tb$net_benefit[tb$model == "cox"], ref, 1e-10)
  Fa <- km_risk(rep(TRUE, nrow(d)))
  expect_close(tb$net_benefit[tb$model == "treat all"], Fa - (1 - Fa) * th / (1 - th), 1e-10)
  expect_output(print(dc), "net benefit")
  with_null_device(plot(dc))
  dc2 <- decision_curve(list(cox = f), data = d, t = 1, thresholds = th)
  expect_close(dc2$table$net_benefit, tb$net_benefit, 1e-12)
})

test_that("anova(): likelihood-ratio tests of nested models agree with survival::anova.coxph", {
  skip_if_not_installed("survival")
  d <- sim_cr(300); d$ev <- as.integer(d$status == 1L)
  a <- anova(cox(Event(time, ev) ~ x, d), cox(Event(time, ev) ~ x + z, d))
  ref <- anova(survival::coxph(survival::Surv(time, ev) ~ x, d), survival::coxph(survival::Surv(time, ev) ~ x + z, d))
  expect_close(a$Chisq[2L], ref$Chisq[2L], 1e-8)
  expect_close(a[["Pr(>|Chi|)"]][2L], ref[["Pr(>|Chi|)"]][2L], 1e-6)
  expect_error(anova(fg(Event(time, status) ~ x, d), fg(Event(time, status) ~ x + z, d)), "pseudo-likelihood")
  expect_error(anova(cox(Event(time, ev) ~ x + z, d), cox(Event(time, ev) ~ x, d)), "nested")
})

test_that("model_scan() reports the time-averaged AUC and ICI", {
  skip_on_cran()
  d <- sim_cr(300); d$ev <- as.integer(d$status == 1L)
  ms <- model_scan(Event(time, ev) ~ x + z, d, depth = "quick", folds = 3)
  expect_true(all(c("iauc_cv", "ici_cv") %in% names(ms@table)))
  ok <- ms@table$fitted & is.finite(ms@table$ibs_cv)
  expect_true(all(ms@table$iauc_cv[ok] > 0.4 & ms@table$iauc_cv[ok] < 1, na.rm = TRUE))
  expect_true(all(c("brier_se_cv", "auc_diff_p_cv") %in% names(ms@per_time)))
  # the time average of the AUC: trapezoids over the grid points where it is defined
  g <- c(0.5, 1, 2, 3); y <- c(NA, 0.6, 0.8, 0.7)
  expect_close(SurvX:::.scan_mean_curve(g, y), ((2 - 1) * (0.6 + 0.8) / 2 + (3 - 2) * (0.8 + 0.7) / 2) / (3 - 1), 1e-14)
  expect_true(is.na(SurvX:::.scan_mean_curve(g, c(NA, NA, 0.5, NA))))
})
