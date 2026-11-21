# Paper 2: real-data analyses. Run from the repository root:
#   Rscript paper2/realdata/realdata.R
# Writes paper2/results/realdata.txt, CSV tables and PDF figures.
library(SurvX)
out_dir <- file.path("paper2", "results")
dir.create(out_dir, showWarnings = FALSE)
safe <- function(expr) tryCatch(expr, error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
se_tab <- function(fit) {
  e <- estimates(fit)
  data.frame(term = paste0(e$cause, ":", e$term), estimate = e$estimate,
             se_ipcw = sqrt(diag(fit@diagnostics$vcov_ipcw)), se_ignored = sqrt(diag(fit@diagnostics$vcov_robust)),
             ratio = sqrt(diag(fit@diagnostics$vcov_ipcw)) / sqrt(diag(fit@diagnostics$vcov_robust)))
}

sink(file.path(out_dir, "realdata.txt"), split = TRUE)
cat("SurvX", as.character(utils::packageVersion("SurvX")), "\n")

# ------------------------------------------------------------------------------------------
# 1. fertility: two causes (live birth; stillbirth or abortion). Weighted estimator for each
#    cause, links, and the standard errors with and without the censoring term.
#    Education only: age is recorded at the end of follow-up (?fertility), so it is not a
#    baseline covariate.
# ------------------------------------------------------------------------------------------
cat("\n==== fertility ====\n")
data(fertility, package = "SurvX")
print(table(fertility$status))
print(aj(Event(time, status) ~ education, fertility))      # includes Gray's test

for (k in 1:2) {
  cat("\n---- cause", k, "----\n")
  fits <- list()
  for (b in c("gompertz", "mix_weibull", "mix_loglogistic")) for (l in c("PH", "PO", "GOR")) {
    f <- safe(directcif(Event(time, status) ~ education, fertility, baseline = b, link = l,
                        estimation = "weighted", cause = k))
    if (!is.null(f) && f@converged) fits[[paste(b, l)]] <- f
  }
  cat("weighted pseudo-log-likelihood by baseline and link:\n")
  print(data.frame(model = names(fits), loglik = vapply(fits, function(f) f@loglik, 0), npar = vapply(fits, function(f) length(coef(f)), 0)),
        row.names = FALSE)
  best <- fits[[which.max(vapply(fits, function(f) f@loglik - length(coef(f)), 0))]]
  cat("\nselected (largest pseudo-log-likelihood minus number of parameters):\n")
  print(best)
  cat("\nstandard errors with the censoring term (ipcw) and with the censoring distribution treated as known:\n")
  st <- se_tab(best); print(st, digits = 4)
  utils::write.csv(st, file.path(out_dir, paste0("fertility_cause", k, "_se.csv")), row.names = FALSE)
  nd <- data.frame(education = factor(levels(fertility$education), levels = levels(fertility$education)))
  pr <- predict(best, nd, times = Inf, cause = k, se = TRUE)
  cat("\nlong-term incidence by education:\n")
  print(data.frame(education = nd$education, Finf = pr$cif[, 1], se = pr$se[, 1]), digits = 3)
  ff <- safe(fg(Event(time, status) ~ education, fertility, cause = k))
  if (!is.null(ff)) { cat("\nsemi-parametric Fine-Gray for comparison:\n"); print(ff) }
  pdf(file.path(out_dir, paste0("fertility_cause", k, ".pdf")), width = 7, height = 5)
  safe(cif_plot(best, cause = k, type = "marginal", reference = TRUE))
  dev.off()
}

# ------------------------------------------------------------------------------------------
# 2. flchain: three causes of death (circulatory, neoplasms, other).
#    The long-run share of each cause is identified only where most subjects die during
#    follow-up. In the whole cohort (from age 50) 72% are alive at the end of follow-up and the
#    long-run split is extrapolation; the main analysis uses the subjects aged 80 or more at
#    enrolment, most of whom die during the 14 years of follow-up. The whole cohort is shown as a
#    sensitivity analysis. Check: the long-run risk of every cause must be at least its
#    Aalen-Johansen estimate at the end of follow-up.
# ------------------------------------------------------------------------------------------
cat("\n==== flchain ====\n")
data(flchain, package = "survival")
prep <- function(fl, age0) {
  fl$time <- pmax(fl$futime, 1) / 365.25
  ch <- as.character(fl$chapter)
  fl$status <- ifelse(fl$death == 0, 0L, ifelse(ch %in% "Circulatory", 1L, ifelse(ch %in% "Neoplasms", 2L, 3L)))
  fl$age5 <- (fl$age - age0) / 5
  fl$male <- as.integer(fl$sex == "M")
  fl$flc_high <- as.integer(fl$flc.grp >= 9)
  fl
}
causes <- c("circulatory", "neoplasms", "other")
analyse <- function(fl, label, age0) {
  fl$age5 <- (fl$age - age0) / 5                       # centred within the group
  cat("\n---- flchain,", label, "----\n")
  cat("status: 0 alive, 1 circulatory, 2 neoplasms, 3 other\n"); print(table(fl$status))
  cat("share alive at the end of follow-up:", round(mean(fl$status == 0L), 3), "\n")
  a <- aj(Event(time, status) ~ 1, fl)
  tend <- max(fl$time[fl$status > 0L])
  aj_end <- sapply(1:3, function(k) .aj_at(a, k, tend))
  cat("Aalen-Johansen incidences at the last death (", round(tend, 1), " years): ", paste(causes, round(aj_end, 3), collapse = "; "), "\n", sep = "")
  form <- Event(time, status) ~ age5 + male + flc_high
  nd <- expand.grid(age5 = c(0, 1), male = 0:1, flc_high = 0:1)
  cons <- list()
  for (b in c("gompertz", "mix_weibull", "logistic3")) for (g in c("weibull", "loglogistic", "lognormal")) {
    f <- safe(directcif(form, fl, baseline = b, estimation = "constrained", cause = 1, cond = g))
    if (!is.null(f) && f@converged) cons[[paste(b, g)]] <- f
  }
  if (!length(cons)) { cat("no constrained model converged\n"); return(invisible(NULL)) }
  ic <- vapply(cons, function(f) AIC(f), 0)
  avg <- t(sapply(cons, function(f) colMeans(sapply(1:3, function(k) predict(f, fl, times = Inf, cause = k)))))
  colnames(avg) <- paste0("longrun_", causes)
  cat("\nconstrained models: AIC and the average long-run risk of each cause\n")
  print(round(cbind(AIC = ic, dAIC = ic - min(ic), avg), 3))
  best <- cons[[which.min(ic)]]
  cat("\nbest AIC:", names(cons)[which.min(ic)], "\n")
  print(summary(best))
  lr <- sapply(1:3, function(k) predict(best, nd, times = Inf, cause = k)); colnames(lr) <- causes
  cat("\nlong-run risks by profile (age5 = (age - ", age0, ") / 5) and their sum:\n", sep = "")
  print(cbind(nd, round(lr, 3), sum = round(rowSums(lr), 6)))
  cat("\ncheck: average long-run risk minus Aalen-Johansen incidence at the end of follow-up (must be >= 0):\n")
  print(round(avg[which.min(ic), ] - aj_end, 3))
  utils::write.csv(cbind(nd, lr), file.path(out_dir, paste0("flchain_longrun_", gsub("[^a-z0-9]+", "_", label), ".csv")), row.names = FALSE)
  fj <- safe(directcif(form, fl, baseline = "gompertz"))
  if (!is.null(fj)) {
    tot <- Reduce(`+`, lapply(1:3, function(k) predict(fj, fl, times = Inf, cause = k)))
    cat("\njoint model without the constraint: share of subjects whose long-run risks add up to more than 1:",
        round(mean(tot > 1), 3), "  largest sum:", round(max(tot), 3), "\n")
  }
  fcs <- safe(csc(form, fl))
  if (!is.null(fcs)) {
    t10 <- min(10, floor(tend))
    cat("\nincidences at", t10, "years, cause-specific Cox and constrained model:\n")
    print(cbind(nd, csc = round(sapply(1:3, function(k) predict(fcs, nd, times = t10, cause = k)), 3),
                constrained = round(sapply(1:3, function(k) predict(best, nd, times = t10, cause = k)), 3)))
  }
  pdf(file.path(out_dir, paste0("flchain_cif_", gsub("[^a-z0-9]+", "_", label), ".pdf")), width = 10, height = 4)
  safe(cif_plot(best, cause = 1:3, type = "marginal", reference = TRUE))
  dev.off()
  invisible(best)
}
.aj_at <- function(a, k, t) {
  cv <- a@curves[[1L]]; v <- cv[[paste0("cif.", k)]]
  v[max(which(cv$time <= t))]
}
fl_all <- prep(flchain, 65)
analyse(fl_all[fl_all$age >= 80, , drop = FALSE], "aged 80 or more", 85)
analyse(fl_all[fl_all$age >= 70, , drop = FALSE], "aged 70 or more", 75)
analyse(prep(flchain, 65), "whole cohort", 65)
sink()
cat("written", file.path(out_dir, "realdata.txt"), "\n")
