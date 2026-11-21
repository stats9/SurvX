# Quick check of the new functions before the simulation study. Run from the repository root,
# with the development version installed:
#   Rscript paper2/check_new_functions.R
# Writes paper2/results/check_new.txt, check_new.pdf and check_new_report.html (a few minutes).
library(SurvX)
dir.create(file.path("paper2", "results"), showWarnings = FALSE)
out <- file.path("paper2", "results")
sink(file.path(out, "check_new.txt"), split = TRUE)
pdf(file.path(out, "check_new.pdf"), width = 9, height = 7)
cat("SurvX", as.character(utils::packageVersion("SurvX")), "|", R.version.string, "|", format(Sys.time()), "\n")
step <- function(title, expr) {
  cat("\n", strrep("=", 90), "\n", title, "\n", strrep("=", 90), "\n", sep = "")
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch({ force(expr); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(sprintf("-- %s (%.1f s)\n", r, proc.time()[["elapsed"]] - t0))
}

# ---------------------------------------------------------------- 1. weighted estimator
set.seed(1)
n <- 600
d <- data.frame(x = rbinom(n, 1, 0.5), z = rnorm(n))
E <- exp(0.5 * d$x - 0.3 * d$z); p1 <- 1 - (1 + E * 0.4 / 0.6)^(-1)
is1 <- runif(n) < p1
v <- runif(n) * p1; u <- ((1 - v)^(-1) - 1) / E
t <- ifelse(is1, log1p(-0.6 * u / 0.4) / -0.6, rlnorm(n, 0.2 + 0.4 * d$x, 0.7))
cens <- runif(n, 0, 6)
d$time <- pmax(pmin(t, cens), 1e-4); d$status <- ifelse(t <= cens, ifelse(is1, 1L, 2L), 0L)
cat("events:", table(d$status), "\n")

step("1a. directcif weighted, PO link, default variance (ipcw)", {
  fw <- directcif(Event(time, status) ~ x + z, d, link = "PO", estimation = "weighted", cause = 1)
  print(fw)
})
step("1b. standard errors: ipcw vs censoring treated as known vs model", {
  se <- cbind(ipcw = sqrt(diag(fw@diagnostics$vcov_ipcw)), robust = sqrt(diag(fw@diagnostics$vcov_robust)),
              model = sqrt(diag(fw@vcov_model)))
  print(round(cbind(se, ratio_ipcw_robust = se[, 1] / se[, 2]), 4))
  fr <- directcif(Event(time, status) ~ x + z, d, link = "PO", estimation = "weighted", cause = 1, variance = "robust")
  cat("variance = 'robust' gives the robust matrix:", isTRUE(all.equal(vcov(fr), fw@diagnostics$vcov_robust)), "\n")
  cat("predicted F(inf) for x = 0, 1 with SE (ipcw):\n")
  print(predict(fw, data.frame(x = 0:1, z = 0), times = Inf, cause = 1, se = TRUE))
})
step("1c. weighted with censoring grid, GOR link, and a spline baseline", {
  print(coef(directcif(Event(time, status) ~ x + z, d, link = "GOR", estimation = "weighted", censoring_grid = 50L)))
  fs <- directcif(Event(time, status) ~ x + z, d, baseline = "spline", estimation = "weighted")
  print(round(sqrt(diag(vcov(fs))) / sqrt(diag(fs@diagnostics$vcov_robust)), 3))
})

# ---------------------------------------------------------------- 2. constrained model, K = 3
d3 <- d
d3$status[d3$status == 2L & runif(n) < 0.4 + 0.3 * d3$x] <- 3L
cat("\nevents with three causes:", table(d3$status), "\n")
step("2a. constrained, three causes, Weibull G", {
  fk <- directcif(Event(time, status) ~ x + z, d3, estimation = "constrained", cause = 1, cond = "weibull")
  print(fk)
})
step("2b. summary, estimates, and long-run incidences add up to one", {
  summary(fk)
  print(estimates(fk))
  nd <- data.frame(x = c(0, 1, 0, 1), z = c(0, 0, 2, -2))
  lr <- sapply(1:3, function(k) predict(fk, nd, times = Inf, cause = k))
  print(cbind(nd, round(lr, 4), sum = rowSums(lr)))
  print(predict(fk, nd[1:2, ], times = c(1, 3), cause = 3, se = TRUE))
})
step("2c. one distribution per cause; cif_data and cif_plot", {
  fk2 <- directcif(Event(time, status) ~ x + z, d3, estimation = "constrained", cause = 1, cond = c("weibull", "lognormal"))
  print(coef(fk2))
  print(head(cif_data(fk2, cause = 1:3, times = c(1, 3))))
  cif_plot(fk2, cause = 1:3)
  cat("AIC (Weibull, Weibull+lognormal):", AIC(fk), AIC(fk2), "\n")
})
step("2d. K = 2 unchanged", {
  print(coef(directcif(Event(time, status) ~ x + z, d, estimation = "constrained", cause = 1, cond = "weibull")))
})

# ---------------------------------------------------------------- 3. model_scan
set.seed(2)
ds <- data.frame(x = rbinom(400, 1, 0.5), z = rnorm(400))
ts <- rweibull(400, 1.4, exp(1 - 0.5 * ds$x + 0.3 * ds$z)); cs <- runif(400, 0, 5)
ds$time <- pmin(ts, cs); ds$status <- as.integer(ts <= cs)
step("3a. model_scan, one event type, standard depth", {
  sc1 <- model_scan(Event(time, status) ~ x + z, ds, verbose = TRUE)
  print(sc1)
  plot(sc1)
})
step("3b. model_scan summary and table", {
  summary(sc1)
  print(as.data.frame(sc1)[, c("model", "fitted", "ibs_cv", "se_diff", "rank", "aic", "np_distance", "folds_failed")], digits = 3)
})
step("3c. model_scan, cure data (plateau)", {
  set.seed(3)
  dc <- data.frame(x = rbinom(500, 1, 0.5))
  cured <- runif(500) < plogis(-0.5 + 1 * dc$x)
  tt <- ifelse(cured, Inf, rweibull(500, 1.5, 2)); cc <- runif(500, 0, 10)
  dc$time <- pmin(tt, cc); dc$status <- as.integer(tt <= cc)
  sc3 <- model_scan(Event(time, status) ~ x, dc, depth = "quick", folds = 3)
  print(sc3)
})
step("3d. model_scan, competing risks (fertility), goal = 'effect'", {
  sc2 <- model_scan(Event(time, status) ~ education + age, fertility, cause = 1, goal = "effect", verbose = TRUE)
  summary(sc2)
  plot(sc2)
  report(sc2, file = file.path(out, "check_new_report.html"))
})
step("3e. model_scan, three causes, cause 3", {
  sc4 <- model_scan(Event(time, status) ~ x + z, d3, cause = 3, depth = "quick", folds = 3)
  print(sc4)
})

# ---------------------------------------------------------------- 4. the simulation code, one replication
step("4. one replication of each simulation study", {
  source(file.path("paper2", "sim", "sim_functions.R"))
  set.seed(10); dA <- A_gen(500, 1, 4); print(table(dA$status))
  print(A_fit_all(dA, 1), digits = 3)
  print(A_truth(1), digits = 3)
  set.seed(11); dB <- B_gen(800, "constrained", 12); print(table(dB$status))
  rB <- B_fit_all(dB); print(rB$est, digits = 3); print(rB$coherence, digits = 3)
  print(B_truth("constrained"), digits = 3)
  set.seed(12); dC <- B_gen(800, "csh", 5); print(table(dC$status))
  rC <- B_fit_all(dC); print(rC$coherence, digits = 3)
})

# ---------------------------------------------------------------- 5. flchain, three causes of death
step("5. flchain: constrained model, aged 80 or more vs whole cohort, with the Aalen-Johansen check", {
  data(flchain, package = "survival")
  fl <- flchain
  fl$time <- pmax(fl$futime, 1) / 365.25
  ch <- as.character(fl$chapter)
  fl$status <- ifelse(fl$death == 0, 0L, ifelse(ch %in% "Circulatory", 1L, ifelse(ch %in% "Neoplasms", 2L, 3L)))
  fl$male <- as.integer(fl$sex == "M"); fl$flc_high <- as.integer(fl$flc.grp >= 9)
  for (grp in list(list("aged 80 or more", fl$age >= 80, 85), list("whole cohort", rep(TRUE, nrow(fl)), 65))) {
    d5 <- fl[grp[[2]], ]; d5$age5 <- (d5$age - grp[[3]]) / 5
    cat("\n-- ", grp[[1]], ": n = ", nrow(d5), ", alive at the end = ", round(mean(d5$status == 0), 3), "\n", sep = "")
    a <- aj(Event(time, status) ~ 1, d5); cv <- a@curves[[1L]]
    aj_end <- sapply(1:3, function(k) utils::tail(cv[[paste0("cif.", k)]], 1L))
    for (b in c("gompertz", "mix_weibull")) {
      f <- tryCatch(directcif(Event(time, status) ~ age5 + male + flc_high, d5, baseline = b, estimation = "constrained",
                              cause = 1, cond = "weibull"), error = function(e) { cat(b, "failed:", conditionMessage(e), "\n"); NULL })
      if (is.null(f)) next
      avg <- colMeans(sapply(1:3, function(k) predict(f, d5, times = Inf, cause = k)))
      cat(sprintf("  %-12s converged %s  AIC %.1f\n", b, f@converged, AIC(f)))
      print(round(rbind(longrun_average = avg, aalen_johansen_end = aj_end, difference = avg - aj_end), 3))
    }
  }
})

# ---------------------------------------------------------------- 6. model_scan plots after the layout fix
step("6. model_scan plots (layout)", {
  plot(sc1)
  plot(sc2)
  print(as.data.frame(sc1)[, c("model", "rank", "ibs_cv", "equivalent", "delta_aic")], digits = 3)
  print(sc2)
})
dev.off()
sink()
cat("written paper2/results/check_new.txt, check_new.pdf, check_new_report.html\n")
