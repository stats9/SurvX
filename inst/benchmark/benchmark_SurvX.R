#!/usr/bin/env Rscript
# Accuracy, stability and speed of SurvX against the established R packages.
#
#   Rscript benchmark_SurvX.R                         full run (about 30-60 minutes)
#   Rscript benchmark_SurvX.R quick                   short trial (a few minutes)
#   Rscript benchmark_SurvX.R R=50 sizes=500,2000,10000,50000 max_seconds=120
#
# Compared on the same simulated data sets:
#   Fine-Gray            SurvX::fg      cmprsk::crr, survival::finegray + coxph (coefficients only)
#   cause-specific Cox   SurvX::csc     survival::coxph per cause, riskRegression::CSC (predicted CIF)
#   Cox                  SurvX::cox     survival::coxph
#   Aalen-Johansen       SurvX::aj      cmprsk::cuminc, survival::survfit (multi-state)
#   Brier score, AUC     SurvX          riskRegression::Score
#   direct models        SurvX::directcif (no other package fits them): convergence and speed only
#
# Accuracy: largest absolute difference from the reference, over the coefficients, the standard
# errors and the predicted CIFs. Stability: share of fits that converged, failed or gave a warning,
# in ordinary data, with heavily tied event times, with rare events and with a nearly separating
# covariate. Speed: elapsed time of one fit (median over the replications).
#
# Writes to benchmark_results/: accuracy.csv, timing.csv, stability.csv, summary.txt, timing.pdf.
# Packages: SurvX (required); survival, cmprsk, riskRegression, prodlim (each optional: a comparison
# is skipped when its package is missing).

suppressPackageStartupMessages(library(SurvX))
args <- commandArgs(trailingOnly = TRUE)
opt <- list(R = 20L, sizes = c(500L, 2000L, 10000L, 50000L), max_seconds = 120, p = 5L, seed = 2026L)
if ("quick" %in% args) opt <- modifyList(opt, list(R = 3L, sizes = c(500L, 2000L), max_seconds = 30))
for (a in grep("=", args, value = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  opt[[kv[1]]] <- if (kv[1] == "sizes") as.integer(strsplit(kv[2], ",")[[1]]) else as.numeric(kv[2])
}
has <- function(p) requireNamespace(p, quietly = TRUE)
pk <- c(survival = has("survival"), cmprsk = has("cmprsk"), riskRegression = has("riskRegression") && has("prodlim"))
out_dir <- "benchmark_results"; dir.create(out_dir, showWarnings = FALSE)
cat("SurvX", format(packageVersion("SurvX")), "| R", format(getRversion()), "| replications", opt$R,
    "| sizes", paste(opt$sizes, collapse = ", "), "\nreference packages:",
    paste(names(pk), ifelse(pk, "yes", "missing"), collapse = ", "), "\n\n")

# ---- data ------------------------------------------------------------------------------------------
# Two causes with exponential cause-specific hazards, p covariates (x1 binary, the rest normal),
# uniform censoring. design: "standard"; "ties" (times rounded to 20 distinct values);
# "rare" (about 3% events of cause 1); "separation" (x1 almost determines cause 1).
sim_data <- function(n, design = "standard", p = 5L) {
  X <- cbind(x1 = rbinom(n, 1, 0.4), matrix(rnorm(n * (p - 1)), n, p - 1, dimnames = list(NULL, paste0("x", 2:p))))
  b1 <- c(0.6, seq(0.4, -0.4, length.out = p - 1)); b2 <- c(-0.3, seq(-0.2, 0.3, length.out = p - 1))
  h1 <- if (design == "rare") 0.02 else 0.3
  if (design == "separation") b1[1] <- 4
  t1 <- rexp(n, h1 * exp(drop(X %*% b1))); t2 <- rexp(n, 0.2 * exp(drop(X %*% b2))); cens <- runif(n, 0, 6)
  time <- pmin(t1, t2, cens)
  status <- ifelse(time == cens, 0L, ifelse(t1 < t2, 1L, 2L))
  if (design == "ties") time <- ceiling(time / max(time) * 20)
  data.frame(time = time, status = status, X)
}

# ---- helpers ---------------------------------------------------------------------------------------
timed <- function(expr) {
  warn <- NULL
  t0 <- proc.time()[["elapsed"]]
  val <- withCallingHandlers(tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "fit_error")),
                             warning = function(w) { warn <<- conditionMessage(w); invokeRestart("muffleWarning") })
  list(value = val, seconds = proc.time()[["elapsed"]] - t0, error = inherits(val, "fit_error"),
       warning = !is.null(warn), message = if (inherits(val, "fit_error")) as.character(val) else if (!is.null(warn)) warn else "")
}
maxdiff <- function(a, b) { a <- as.numeric(a); b <- as.numeric(b); if (length(a) != length(b)) return(NA_real_); ok <- is.finite(a) & is.finite(b); if (!any(ok)) NA_real_ else max(abs(a[ok] - b[ok])) }

acc <- list(); tim <- list(); stab <- list()
add_acc <- function(...) acc[[length(acc) + 1L]] <<- data.frame(...)
add_tim <- function(...) tim[[length(tim) + 1L]] <<- data.frame(...)
add_stab <- function(...) stab[[length(stab) + 1L]] <<- data.frame(...)
too_slow <- character(0)        # methods skipped for larger n after exceeding max_seconds

run_method <- function(name, n, design, rep, expr) {
  if (name %in% too_slow) { add_tim(method = name, n = n, design = design, rep = rep, seconds = NA, status = "skipped (too slow)"); return(NULL) }
  r <- timed(expr)
  add_tim(method = name, n = n, design = design, rep = rep, seconds = r$seconds,
          status = if (r$error) "error" else if (r$warning) "warning" else "ok")
  add_stab(method = name, n = n, design = design, rep = rep, error = r$error, warning = r$warning, message = substr(r$message, 1, 120))
  if (r$seconds > opt$max_seconds) too_slow <<- union(too_slow, name)
  if (r$error) NULL else r$value
}

# ---- one data set ----------------------------------------------------------------------------------
one_dataset <- function(n, design, rep) {
  d <- sim_data(n, design, opt$p)
  xn <- grep("^x", names(d), value = TRUE)
  fm <- as.formula(paste("Event(time, status) ~", paste(xn, collapse = " + ")))
  X <- as.matrix(d[, xn])
  tq <- unname(quantile(d$time[d$status > 0], c(0.25, 0.5, 0.75)))
  nd <- d[seq_len(min(50, n)), ]
  # an error while computing a difference is recorded as NA instead of stopping the run
  A <- function(comparison, quantity, max_abs_diff)
    add_acc(n = n, design = design, rep = rep, comparison = comparison, quantity = quantity,
            max_abs_diff = tryCatch(max_abs_diff, error = function(e) NA_real_))

  # Fine-Gray
  ours <- run_method("fg: SurvX", n, design, rep, fg(fm, d))
  if (!is.null(ours)) add_stab(method = "fg: SurvX (converged)", n = n, design = design, rep = rep, error = !isTRUE(ours@converged), warning = FALSE, message = "")
  if (pk[["cmprsk"]]) {
    ref <- run_method("fg: cmprsk::crr", n, design, rep, cmprsk::crr(d$time, d$status, X, failcode = 1, cencode = 0))
    if (!is.null(ours) && !is.null(ref)) {
      A(comparison = "fg vs crr", quantity = "coefficients", max_abs_diff = maxdiff(coef(ours), ref$coef))
      A(comparison = "fg vs crr", quantity = "robust standard errors", max_abs_diff = maxdiff(sqrt(diag(vcov(ours))), sqrt(diag(ref$var))))
      pr <- predict(ref, X[seq_len(nrow(nd)), , drop = FALSE])
      refp <- sapply(tq, function(s) pr[max(which(pr[, 1] <= s)), -1])
      A(comparison = "fg vs crr", quantity = "predicted CIF", max_abs_diff = maxdiff(predict(ours, nd, times = tq), refp))
    }
  }
  if (pk[["survival"]]) {
    ref2 <- run_method("fg: survival::finegray + coxph", n, design, rep, {
      dd <- d; dd$ev <- factor(dd$status, 0:2)
      fd <- survival::finegray(as.formula(paste("survival::Surv(time, ev) ~", paste(xn, collapse = " + "))), data = dd, etype = "1")
      survival::coxph(as.formula(paste("survival::Surv(fgstart, fgstop, fgstatus) ~", paste(xn, collapse = " + "))), data = fd, weights = fgwt)
    })
    if (!is.null(ours) && !is.null(ref2)) A(comparison = "fg vs finegray + coxph", quantity = "coefficients", max_abs_diff = maxdiff(coef(ours), coef(ref2)))
  }

  # cause-specific Cox
  oc <- run_method("csc: SurvX", n, design, rep, csc(fm, d))
  if (pk[["survival"]]) {
    rc <- run_method("csc: survival::coxph x 2", n, design, rep, lapply(1:2, function(k)
      survival::coxph(as.formula(paste0("survival::Surv(time, status == ", k, ") ~ ", paste(xn, collapse = " + "))), data = d)))
    if (!is.null(oc) && !is.null(rc)) {
      A(comparison = "csc vs coxph", quantity = "coefficients", max_abs_diff = maxdiff(unlist(coef(oc)), unlist(lapply(rc, coef))))
      A(comparison = "csc vs coxph", quantity = "standard errors", max_abs_diff = maxdiff(unlist(lapply(seq_along(rc), function(k) sqrt(diag(vcov(oc)[[k]])))), unlist(lapply(rc, function(m) sqrt(diag(vcov(m)))))))
    }
  }
  if (pk[["riskRegression"]]) {
    Hist <- prodlim::Hist
    rr <- run_method("csc: riskRegression::CSC", n, design, rep,
                     riskRegression::CSC(as.formula(paste("Hist(time, status) ~", paste(xn, collapse = " + "))), data = d))
    if (!is.null(oc) && !is.null(rr))
      A(comparison = "csc vs riskRegression::CSC", quantity = "predicted CIF", max_abs_diff =
          maxdiff(predict(oc, nd, times = tq, cause = 1), riskRegression::predictRisk(rr, newdata = nd, times = tq, cause = 1)))
  }

  # Cox for all-cause failure
  d$dead <- as.integer(d$status > 0)
  fmc <- as.formula(paste("Event(time, dead) ~", paste(xn, collapse = " + ")))
  ox <- run_method("cox: SurvX", n, design, rep, cox(fmc, d))
  if (pk[["survival"]]) {
    rx <- run_method("cox: survival::coxph", n, design, rep, survival::coxph(as.formula(paste("survival::Surv(time, dead) ~", paste(xn, collapse = " + "))), data = d))
    if (!is.null(ox) && !is.null(rx)) {
      A(comparison = "cox vs coxph", quantity = "coefficients", max_abs_diff = maxdiff(coef(ox), coef(rx)))
      A(comparison = "cox vs coxph", quantity = "standard errors", max_abs_diff = maxdiff(sqrt(diag(vcov(ox))), sqrt(diag(vcov(rx)))))
    }
  }

  # Aalen-Johansen
  oa <- run_method("aj: SurvX", n, design, rep, aj(Event(time, status) ~ x1, d))
  if (pk[["cmprsk"]]) {
    ra <- run_method("aj: cmprsk::cuminc", n, design, rep, cmprsk::cuminc(d$time, d$status, d$x1))
    if (!is.null(oa) && !is.null(ra)) {
      tt <- cmprsk::timepoints(ra, tq); ours_t <- summary(oa, times = tq)@tables
      A(comparison = "aj vs cuminc", quantity = "CIF of cause 1", max_abs_diff = maxdiff(c(ours_t[["0"]]$cif.1, ours_t[["1"]]$cif.1), c(tt$est["0 1", ], tt$est["1 1", ])))
      A(comparison = "aj vs cuminc", quantity = "Gray test statistic", max_abs_diff = maxdiff(oa@tests$statistic, ra$Tests[, "stat"]))
    }
  }
  if (pk[["survival"]]) run_method("aj: survival::survfit (multi-state)", n, design, rep,
                                   survival::survfit(survival::Surv(time, factor(status, 0:2)) ~ x1, data = d))

  # Brier score and AUC
  if (!is.null(ours)) {
    risk <- drop(predict(ours, times = tq[2]))
    ob <- run_method("brier + auc: SurvX", n, design, rep, c(brier_score(d$time, d$status, risk, tq[2])[["brier"]], auc_t(d$time, d$status, risk, tq[2])))
    if (pk[["riskRegression"]]) {
      Hist <- prodlim::Hist
      sc <- run_method("brier + auc: riskRegression::Score", n, design, rep,
                       riskRegression::Score(list(m = matrix(risk)), formula = Hist(time, status) ~ 1, data = d, times = tq[2], cause = 1,
                                             metrics = c("auc", "brier"), cens.model = "km", conf.int = FALSE, null.model = FALSE))
      if (!is.null(ob) && !is.null(sc)) {
        bs <- as.data.frame(sc$Brier$score); au <- as.data.frame(sc$AUC$score)
        A(comparison = "metrics vs Score", quantity = "Brier score", max_abs_diff = maxdiff(ob[1], bs$Brier[bs$model == "m"]))
        A(comparison = "metrics vs Score", quantity = "AUC", max_abs_diff = maxdiff(ob[2], au$AUC[au$model == "m"]))
      }
    }
  }

  # direct models (no reference package): convergence and speed
  for (lk in c("PH", "PO", "GOR")) {
    f <- run_method(paste("directcif joint", lk), n, design, rep, directcif(fm, d, link = lk))
    if (!is.null(f)) add_stab(method = paste("directcif joint", lk, "(converged)"), n = n, design = design, rep = rep,
                              error = !isTRUE(f@converged) || any(!is.finite(sqrt(diag(vcov(f))))), warning = FALSE, message = "")
  }
  f <- run_method("directcif weighted PH", n, design, rep, directcif(fm, d, estimation = "weighted"))
  f <- run_method("directcif constrained PH", n, design, rep, directcif(fm, d, estimation = "constrained", cond = "weibull"))
}

# ---- run -------------------------------------------------------------------------------------------
set.seed(opt$seed)
designs <- c("standard", "ties", "rare", "separation")
for (n in opt$sizes) for (design in designs) {
  if (design != "standard" && n > 2000) next          # stress designs only at moderate n
  cat(sprintf("n = %6d  %-10s ", n, design))
  for (rep in seq_len(opt$R)) {
    tryCatch(one_dataset(n, design, rep), error = function(e) cat("\n  data set", rep, "stopped:", conditionMessage(e), "\n"))
    cat(".")
  }
  cat("\n")
}

acc <- do.call(rbind, acc); tim <- do.call(rbind, tim); stab <- do.call(rbind, stab)
write.csv(acc, file.path(out_dir, "accuracy.csv"), row.names = FALSE)
write.csv(tim, file.path(out_dir, "timing.csv"), row.names = FALSE)
write.csv(stab, file.path(out_dir, "stability.csv"), row.names = FALSE)

# ---- summaries -------------------------------------------------------------------------------------
fmt <- function(x) formatC(x, format = "g", digits = 3)
acc_s <- aggregate(max_abs_diff ~ comparison + quantity + design, acc, function(x) max(x, na.rm = TRUE))
names(acc_s)[4] <- "largest_difference"
tim_ok <- tim[!is.na(tim$seconds) & tim$status != "error" & tim$design == "standard", ]
tim_s <- aggregate(seconds ~ method + n, tim_ok, median)
tim_w <- reshape(tim_s, idvar = "method", timevar = "n", direction = "wide"); names(tim_w) <- sub("seconds.", "n = ", names(tim_w))
stab_s <- aggregate(cbind(failed_pct = 100 * error, warning_pct = 100 * warning) ~ method + design, stab, mean)

sink(file.path(out_dir, "summary.txt"), split = TRUE)
cat("SurvX", format(packageVersion("SurvX")), "| R", format(getRversion()), "|", opt$R, "replications per cell\n")
cat("\n== Accuracy: largest absolute difference from the reference over all replications ==\n")
print(transform(acc_s, largest_difference = fmt(largest_difference)), row.names = FALSE)
cat("\n== Speed: median elapsed seconds of one fit, standard design ==\n")
print(tim_w, row.names = FALSE, digits = 3)
if (length(too_slow)) cat("\nskipped at larger n after exceeding", opt$max_seconds, "s:", paste(too_slow, collapse = ", "), "\n")
cat("\n== Stability: % of fits that failed (error or no convergence) or warned ==\n")
print(stab_s[stab_s$failed_pct > 0 | stab_s$warning_pct > 0 | grepl("converged", stab_s$method), ], row.names = FALSE, digits = 3)
sink()

grDevices::pdf(file.path(out_dir, "timing.pdf"), width = 8, height = 5.5)
ms <- unique(tim_s$method); cols <- grDevices::hcl.colors(length(ms), "Dark 3")
graphics::par(mar = c(4.5, 4.5, 1, 14), xpd = NA)
graphics::plot(range(tim_s$n), range(pmax(tim_s$seconds, 1e-3)), type = "n", log = "xy", xlab = "sample size n", ylab = "median seconds per fit")
for (i in seq_along(ms)) { s <- tim_s[tim_s$method == ms[i], ]; graphics::lines(s$n, pmax(s$seconds, 1e-3), type = "b", col = cols[i], pch = 16 + (i %% 3), lwd = 2) }
graphics::legend(10^graphics::par("usr")[2], 10^graphics::par("usr")[4], legend = ms, col = cols, pch = 16 + (seq_along(ms) %% 3), lwd = 2, bty = "n", cex = 0.7, xjust = 0)
invisible(grDevices::dev.off())
cat("\nwritten to", normalizePath(out_dir), ": accuracy.csv, timing.csv, stability.csv, summary.txt, timing.pdf\n")
